#ifndef BASE_MUTATION
#define BASE_MUTATION

class BaseMutation
{
public:
    virtual void operator()(double *population_genes, const double * min_value, const double * max_value, int population_size, int genes_number, const int * is_mutation_applicable) = 0;
};

#endif
