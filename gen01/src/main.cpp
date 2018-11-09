#include <iostream>
#include <fstream>
#include "../include/population.h"

using namespace std;

int main(){
    clock_t t_start = clock();
    srand(time(NULL));

    population pop(20);
    pop.evolve(E);
    cout << "best_phenotype=" << (double)pop.decode_gray(pop.best_genotype) << endl;
    cout << "fitness:" << endl;
    for (int i = 0; i < pop.pop_size; i++){
        cout << pop.fitness[i] << endl;
    }
    
    //********** Writing data to files ***************
    FILE *f1 = fopen("../data/fit_best.csv", "w");
    if(f1 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int i = 0; i < E; i++){
        fprintf(f1, "%lf\n", pop.hist_fit_best[i]);
    }
    fclose(f1);
    //************************************************
    FILE *f2 = fopen("../data/avg_fit.csv", "w");
    if(f2 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int i = 0; i < E; i++){
        fprintf(f2, "%lf\n", pop.hist_avg_fit[i]);
    }
    fclose(f2);
    //************************************************
    FILE *f3 = fopen("../data/sd.csv", "w");
    if(f3 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int i = 0; i < E; i++){
        fprintf(f3, "%lf\n", pop.hist_sd[i]);
    }
    fclose(f3);
    //************************************************

    cout << "Time elapsed = " << (double)(clock()-t_start)/CLOCKS_PER_SEC << " s" << endl;

    return 0;
}