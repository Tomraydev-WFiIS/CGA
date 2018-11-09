#ifndef POPULATION_H
#define POPULATION_H

#define E 100

class population {
    public:
        int pop_size;
        int epoch;
        const int genotype_len = 22;
        bool **genotypes;
        double *phenotypes;
        double *fitness;
        double *p_selection;
        double *acc_p_selection;
        double p_mutation;

        //Best solution
        bool *best_genotype;
        double best_phenotype;

        //For plotting
        double *hist_avg_fit;
        double *hist_fit_best;
        double *hist_sd;
        void record_all();

        void print_genotypes();
        void print_phenotypes();
        std::string get_genotype(bool *);
        void print_fitness();
        double decode(bool *, int len=22);
        void decode_all();
        //
        double decode_gray(bool *, int len=22);
        void decode_gray_all();
        //
        double eval(double);
        void eval_all();

        void evolve(int );
        int select_i();
        void mutate();

        population(int );
        ~population(void);   
};



#endif