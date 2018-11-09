#include <iostream>
#include <fstream>
#include <cmath>
#include "../include/population.h"

using namespace std;

int main(){
    clock_t t_start = clock();
    srand(time(NULL));

    double global_fit_best[E];
    double global_avg_fit[E];
    double global_sd_fit_best[E];

    population pop[I];
    for(int i = 0; i < I; i++){
        pop[i].evolve(E);
        printf("best_phenotype = %.8f\n", (double)pop[i].decode_gray(pop[i].best_genotype));
    }

    //Combining data for plotting
    for(int e = 0; e < E; e++){
        double A_fit_best = 0.0;
        double A_avg_fit = 0.0;
        for(int i = 0; i < I; i++){
            A_fit_best += pop[i].hist_fit_best[e];
            A_avg_fit += pop[i].hist_avg_fit[e];
        }
        A_fit_best = A_fit_best/I;
        A_avg_fit = A_avg_fit/I;

        double variance = 0.0;
        double sd = 0.0;
        for(int i = 0; i < I; i++){
            variance += pow(pop[i].hist_fit_best[e] - A_fit_best, 2);
        }
        variance = variance/I;
        sd = pow(variance, 0.5);

        global_fit_best[e] = A_fit_best;
        global_avg_fit[e] = A_avg_fit;
        global_sd_fit_best[e] = sd;
    }

    // Printing the best solution
    double global_best_fit = 0.0;
    double global_best_phenotype = 0.0;
    for(int i = 0; i < I; i++){
        if (pop[i].eval(pop[i].best_phenotype) > global_best_fit){
            global_best_fit = pop[i].eval(pop[i].best_phenotype);
            global_best_phenotype = pop[i].best_phenotype;
        }
    }
    cout << "best fitness = " << global_best_fit << " for x = " << global_best_phenotype << endl;

    //************ Writing to files ***********************
    FILE *f1 = fopen("../data/global_fit_best.csv", "w");
    if(f1 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int e = 0; e < E; e++){
        fprintf(f1, "%lf\n", global_fit_best[e]);
    }
    fclose(f1);

    //*****************************************************
    FILE *f2 = fopen("../data/global_avg_fit.csv", "w");
    if(f2 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int e = 0; e < E; e++){
        fprintf(f2, "%lf\n", global_avg_fit[e]);
    }
    fclose(f2);
    //*****************************************************
    FILE *f3 = fopen("../data/global_sd_fit_best.csv", "w");
    if(f3 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int e = 0; e < E; e++){
        fprintf(f3, "%lf\n", global_sd_fit_best[e]);
    }
    fclose(f3);
    //*****************************************************
    FILE *f4 = fopen("../data/all_data.csv", "w");
    if(f4 == NULL){
        cout << "error opening file" << endl; 
        exit(1); 
    }
    for (int e = 0; e < E; e++){
        fprintf(f4, "%d,%lf,%lf\n",e+1, global_fit_best[e], global_sd_fit_best[e]);
    }
    fclose(f4);
    //*****************************************************

    cout << "Time elapsed = " << (double)(clock()-t_start)/CLOCKS_PER_SEC << " s" << endl;
    return 0;
}