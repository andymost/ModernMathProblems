#include <iostream>
#include "base.h"
#include <math.h>

int main() {
    cout.unsetf(ios::floatfield );
    cout.precision(10);
    
    // прямая задача и значения в счетчиках
    grid_type grid = read_grid();
    vector<double> b = make_global_vector(grid);
    matrix_type A = make_global_matrix(grid);
    boundary_cond_2(grid.second_cond, b);
    boundary_cond_1(grid.first_cond, b, A);
    vector<double> q = solve_equation(A, b);
    vector<double> eps_true_vec = epsilon_noized_vec(q, grid);
    vector<double> weigth = {1, 1, 1, 1, 1, 1, 1};
    print_vector(eps_true_vec, "Epsilon true");
    
    // процедур расчета
    grid_type grid_c = grid;
    vector<double> u = setStartParams(grid_c), u_prev = u;                                           // вектора параметров среды
    // параметры процедуры
    int iteration_count = 0;
    // расчет начального J(u)
    b = make_global_vector(grid_c);
    A = make_global_matrix(grid_c);
    boundary_cond_2(grid_c.second_cond, b);
    boundary_cond_1(grid_c.first_cond, b, A);
    q = solve_equation(A, b);
    vector<double> eps_vec_prev = epsilon_vec(q, grid_c);
    vector<double> delta_eps_vec_prev = delta_epsilon(eps_vec_prev, eps_true_vec);
    double d_eps_l_prev = vector_length(delta_eps_vec_prev);
    double beta = beta_start;
    double alpha = 0;
    
    while (iteration_count < max_iter_count) {
        auto delta_u = get_delta_u(grid_c, eps_true_vec, weigth, alpha);
        
        // выход по малому изменению параметров
        if(vector_length(delta_u) < delta_u_min_diff){
            return 0;
        }
        
        vector<double> u_step = vector_sum(u_prev, scale_vector(delta_u, beta));
        grid_type step_grid = setParams(grid_c, u_step);
        
        vector<double> b = make_global_vector(step_grid);
        matrix_type A = make_global_matrix(step_grid);
        boundary_cond_2(step_grid.second_cond, b);
        boundary_cond_1(step_grid.first_cond, b, A);
        vector<double> q = solve_equation(A, b);
        vector<double> eps_vec = epsilon_vec(q, step_grid);
        vector<double> delta_eps_vec = delta_epsilon(eps_vec, eps_true_vec);
        double d_eps_l = vector_length(delta_eps_vec);
        
        bool is_lower = d_eps_l < d_eps_l_prev;
        
        if(is_lower) {
            d_eps_l_prev = d_eps_l;
            u_prev = u_step;
            iteration_count++;
            grid_c = step_grid;
            beta = 1;
            
            cout<<endl<<"J(u): "<<d_eps_l<<endl;
            print_vector(u_step, "RESULT VECTOR");
            
            vector<double> diff;
            vector<double> true_u = {2,2,5, 5, 5,3,3};
            diff = vector_sum(u_step, scale_vector(true_u, -1));
            print_vector(diff, "DIFF VECTOR");
            cout<<"MEAN ERR: "<<vector_length(diff)/7<<endl;
        } else {
            if(beta < 1/1024){
                return 1;
            }
            cout<<"Beta split: "<< beta <<endl;
            beta = beta/2;
        }
    }
    
    cout<<"Iter count: "<<iteration_count<<endl;
    return 0;
}
