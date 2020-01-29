//
// Created by andrey on 1/29/20.
//

#ifndef GA_MALECKAR_H

#include <cmath>
#include <iostream>

#define CONSTANT_ARRAY_SIZE 51
#define ALGEBRAIC_ARRAY_SIZE 70
#define STATE_ARRAY_SIZE 30

struct State {
    double V;
    double Na_c;
    double Na_i;
    double m;
    double h1;
    double h2;
    double Ca_d;
    double d_L;
    double f_L1;
    double f_L2;
    double K_c;
    double K_i;
    double r;
    double s;
    double a_ur;
    double i_ur;
    double n;
    double pa;
    double Ca_c;
    double Ca_i;
    double O_C;
    double O_TC;
    double O_TMgC;
    double O_TMgMg;
    double O;
    double Ca_rel;
    double Ca_up;
    double O_Calse;
    double F1;
    double F2;
};

int action_potential(struct State *initial_state, double *scaling_coefficients, double *AP, float CL, float amp,
                     int current_time, int iso, int baseline_index, int amount_of_baselines, int amount_of_genes);

#define GA_MALECKAR_H

#endif //GA_MALECKAR_H
