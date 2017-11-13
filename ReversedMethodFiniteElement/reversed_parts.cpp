//
//  reversed_parts.cpp
//  ReversedMethodFiniteElement
//
//  Created by Кочнев Андрей Владимирович on 13/11/2017.
//  Copyright © 2017 Кочнев Андрей Владимирович. All rights reserved.
//

#include <stdio.h>
#include "base.h"

// d_delta_epsilons - матрица значений частных производных вектора замера по переменной среды
// размерность: params_count x reciever_counts
matrix_type reversed_matrix_a(matrix_type d_delta_epsilons, vector<double> weight, double alpha){
    matrix_type res;
    for (int i = 0; i < params_count; i++) {
        vector<double> tmp;
        res.push_back(tmp);
        
        for (int j = 0; j < params_count; j++) {
            double val = 0;
            for (int k = 0; k < reciver_count; k++) {
                val += weight[k] * weight[k] * d_delta_epsilons[i][k] * d_delta_epsilons[j][k];
            }
            
            if(i == j) {
                val += 0;
            }
            
            res[i].push_back(val);
        }
    }
    
    return res;
}

vector<double> reversed_vector_f(vector<double> delta_epsilon_with_true, matrix_type d_delta_epsilons, vector<double> params, vector<double> weight, double alpha) {
    vector<double> res;
    
    double mean_param = mean(params);
    
    for (int i = 0; i < params_count; i++) {
        double val = 0;
        
        for (int k = 0; k < reciver_count; k++) {
            val -= weight[k] * weight[k] * delta_epsilon_with_true[k] * d_delta_epsilons[i][k];
        }
        
        val -= 0 * (params[i] - mean_param);
        
        res.push_back(val);
    }
    
    return res;
}

vector<double> get_delta_u(grid_type grid, vector<double> true_epsilon, vector<double> weight, double alpha) {
    vector<double> delta_u;
    vector<double> b = make_global_vector(grid);
    matrix_type A = make_global_matrix(grid);
    boundary_cond_2(grid.second_cond, b);
    boundary_cond_1(grid.first_cond, b, A);
    vector<double> q = solve_equation(A, b);
    vector<double> epsilon_v = epsilon_vec(q, grid);
    vector<double> delta_epsilon_v = delta_epsilon(epsilon_v, true_epsilon);
    vector<double> u = get_params(grid);
    
    matrix_type d_delta_epsilons =  deriative_delta_epsilon(u, grid, delta_epsilon_v, true_epsilon);
    matrix_type r_A = reversed_matrix_a(d_delta_epsilons, weight, alpha);
    vector<double> r_b = reversed_vector_f(delta_epsilon_v, d_delta_epsilons, u, weight, alpha);
    delta_u = solve_equation(r_A, r_b);
    
    
    return delta_u;
}

double calc_alpha(vector<double> u_step, double d_eps_l) {
    double deviation = calc_deviation(u_step);
    return (gamma_val * d_eps_l)/deviation;
    
}
