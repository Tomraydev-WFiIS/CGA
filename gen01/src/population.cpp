#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include "../include/population.h"

using namespace std;

// Class constructor
population::population(int n){
    pop_size = n;
    epoch = 0;
    best_phenotype = 0;
    p_mutation = 0.03;

    genotypes = new bool*[pop_size];
    for(int i = 0; i < pop_size; i++){
        genotypes[i] = new bool[genotype_len];
    }
    phenotypes = new double[pop_size];
    fitness = new double[pop_size];
    p_selection = new double[pop_size];
    acc_p_selection = new double[pop_size];
    best_genotype = new bool[genotype_len];
    hist_avg_fit = new double[E];
    hist_fit_best = new double[E];
    hist_sd = new double[E];

    //Generating random genotypes
    for(int i = 0; i < pop_size; i++){ 
        for(int j = 0; j < genotype_len; j++){
            genotypes[i][j] = rand()%2;
        }
    }

    //Translating to phenotypes
    decode_gray_all();
    //Evaluating the fitness function
    eval_all();
    //Saving data to analyze
    record_all();
}
// Class destructor
population::~population(void){
    for(int i = 0; i < pop_size; i++){
        delete[] genotypes[i];
    }
    delete[] genotypes;
    delete[] phenotypes;
    delete[] fitness;
    delete[] p_selection;
    delete[] acc_p_selection;
    delete[] best_genotype;
    delete[] hist_avg_fit;
    delete[] hist_fit_best;
    delete[] hist_sd;
}

double population::decode(bool *genotype, int len){
    int genes = 0;
    for (int i = 0; i < genotype_len; i++){
        genes = genes | (genotype[i] << i);
    }
    double phenotype = (double) (genes/4194304.0)*3 - 1;
    return phenotype;
}

void population::decode_all(void){
    for (int i = 0; i < pop_size; i++){
        phenotypes[i] = decode(genotypes[i]);
    }
}

double population::decode_gray(bool *genotype, int len){
    bool genotype_nbc[22] = {0};
    genotype_nbc[0] = genotype[0];
    // Converting Gray encoded array to nbc encoded array
    for(int i = 1; i < len; i++){
        genotype_nbc[i] = genotype_nbc[i-1] ^ genotype[i];    
    }
    
    unsigned int genes = 0;
    for (int i = 0; i < len; i++){
        genes = genes | (genotype_nbc[i] << i);
    }
    
    double phenotype = (double) (genes/4194303.0)*3 - 1;
    return phenotype;
}

void population::decode_gray_all(void){
    for (int i = 0; i < pop_size; i++){
        phenotypes[i] = decode_gray(genotypes[i]);
    }
}

double population::eval(double x){
    return x * sin(10 * M_PI * x) + 2;
}

void population::eval_all(){
    for(int i = 0; i < pop_size; i++){
        fitness[i] = eval(phenotypes[i]);
    }
}

void population::evolve(int max_epoch){
    while (epoch != (max_epoch-1) ){
        //Allocating memory
        bool **temp_genotypes = new bool*[pop_size];
        for (int i = 0; i < pop_size; i++){
            temp_genotypes[i] = new bool[genotype_len]();
        }

        //Calculating the probabality of crossover
        double total_fitness = 0;
        for (int i = 0; i < pop_size; i++){
            total_fitness += fitness[i];
        }
        for (int i = 0; i < pop_size; i++){
            p_selection[i] = fitness[i] / total_fitness;
        }
        acc_p_selection[0] = p_selection[0];
        for(int i = 1; i < pop_size; i++){
            acc_p_selection[i] = p_selection[i] + acc_p_selection[i-1];
        }
        //Roulette selection and crossover
        for (int i = 0; i < pop_size; i+=2){
            int i_1 = select_i();
            int i_2 = select_i();
            int locus = rand()%22;
            //cout << "locus = " << locus << endl;

            memcpy(&(temp_genotypes[i][0]), &(genotypes[i_1][0]), locus);
            memcpy(&(temp_genotypes[i][locus]), &(genotypes[i_2][locus]), genotype_len-locus);

            memcpy(&(temp_genotypes[i+1][0]), &(genotypes[i_2][0]), locus);
            memcpy(&(temp_genotypes[i+1][locus]), &(genotypes[i_1][locus]), genotype_len-locus);    
        }
        /*
        for(int i = 0; i < pop_size; i++){
            for(int j = 0; j < genotype_len; j++){
                cout << temp_genotypes[i][j];
            }
            cout << "\n";
        }
        */

        //Updating the population
        for (int i = 0; i < pop_size; i++){
            memcpy(&(genotypes[i][0]), &(temp_genotypes[i][0]), genotype_len);
        }
        mutate();
        decode_gray_all();
        eval_all();

        //Deallocating temporary memory
        for(int i = 0; i < pop_size; i++){
            delete[] temp_genotypes[i];
        }
        delete[] temp_genotypes;

        /*
        cout << "p_selection:\n";
        for (int i = 0; i < pop_size; i++){
            cout << p_selection[i] << endl;
        }
        */

        // Updating the best solution
        int i_best = 0;
        for(int i = 1; i < pop_size; i++){
            if (fitness[i] > fitness[i_best]){
                i_best = i;
            }
        }
        if ( fitness[i_best] > eval(best_phenotype) ){
            memcpy( &(best_genotype[0]), &(genotypes[i_best][0]), genotype_len);
            best_phenotype = decode_gray(best_genotype);
        }

        epoch++;
        record_all();
        cout << "x_best = " << get_genotype(genotypes[i_best]) << "    fit_best = " << fitness[i_best] << "    avg_fit = " << hist_avg_fit[epoch] << "   std = " << hist_sd[epoch] << endl;
    }
}

void population::record_all(){
    double total_fitness = 0;
    for (int i = 0; i < pop_size; i++){
        total_fitness += fitness[i];
    }
    
    int i_best = 0;
    for(int i = 1; i < pop_size; i++){
        if (fitness[i] > fitness[i_best]){
            i_best = i;
        }
    }

    double avg_fit = total_fitness/pop_size;
    double variance = 0;
    for(int i = 0; i < pop_size; i++){
        variance += pow(fitness[i] - avg_fit,2);
    }
    variance = variance/pop_size;
    double std = pow(variance,0.5);
    
    hist_avg_fit[epoch] = avg_fit;
    hist_fit_best[epoch] = fitness[i_best];
    hist_sd[epoch] = std;
}
void population::mutate(){
    for (int i = 0; i < pop_size; i++){
        for (int j = 0; j < genotype_len; j++){
            if ( (double)rand()/RAND_MAX < p_mutation ){
                genotypes[i][j] = !(genotypes[i][j]);
            }
        }
    }
}

int population::select_i(){
    double pill = rand()/(double)RAND_MAX;
    for (int i = 0; i < pop_size; i++){
        if (pill < acc_p_selection[i]){
            //cout << "i_win = " << i << endl;
            return i;
        }
    }
    return -1; //if error
}

string population::get_genotype(bool *genotype){
    string return_string = "";
    for (int i = 0; i < genotype_len; i++){
        return_string += ((char)genotype[i] + '0' );
    }
    return return_string;
}

void population::print_genotypes(){
    for(int i = 0; i < pop_size; i++){ 
        for(int j = 0; j < genotype_len; j++){
            cout << genotypes[i][j];
        }
        cout << "\n";
    }
}

void population::print_phenotypes(){
    for(int i = 0; i < pop_size; i++){
        cout << phenotypes[i] << endl;
    }
}

void population::print_fitness(){
    for(int i = 0; i < pop_size; i++){
        cout << fitness[i] << endl;
    }
}