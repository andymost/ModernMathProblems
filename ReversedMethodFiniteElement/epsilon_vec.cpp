//
//  epsilon_vec.cpp
//  ReversedMethodFiniteElement
//
//  Created by Кочнев Андрей Владимирович on 29/10/2017.
//  Copyright © 2017 Кочнев Андрей Владимирович. All rights reserved.
//

#include "base.h"
#include <math.h>
#include <random>
#include <iostream>


// суть процедуры зашумления
// гененрируем нормально распределенную случайную велечину в районе 0 с отклонением в процентах шума
// домножаем полученную велечину на значение и добавляем к значению
vector<double> epsilon_noized_vec(vector<double> q, grid_type grid) {
    vector<double> result;
    long long seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    
    for (int i = 0; i < reciver_count; ++i) {
        normal_distribution<double> distribution(0, noize_percentage[i]);
        double val = calc_value(q, grid, reciver_x_pos[i], reciver_y_pos[i]);
        double noized_val = val + distribution(generator) * val;
        result.push_back(noized_val);
    }
    return result;
}

vector<double> epsilon_vec(vector<double> q, grid_type grid) {
    vector<double> result;
    
    for (int i = 0; i < reciver_count; ++i) {
        double val = calc_value(q, grid, reciver_x_pos[i], reciver_y_pos[i]);
        result.push_back(val);
    }
    return result;
}

vector<double> delta_epsilon(vector<double> eps_1, vector<double> eps_2) {
    vector<double> res;
    for (int i = 0; i < eps_1.size(); i ++) {
        res.push_back(eps_1[i] - eps_2[i]);
    }
    
    return res;
}

matrix_type deriative_delta_epsilon(
    vector<double> u,                   // значение параметра базовое
    grid_type grid,                     // сетка
    vector<double> delta_epsilon_v,     // delta_epsilon для базового значения параметра
    vector<double> true_epsilon         // эксперементальное значение в приемниках
) {
    // результат - вектор производных векторов
    // нужен для матрицы и правой части обратной задачи
    // по строкам переменные (u)
    // по столбцам значения в приемниках
    
    //          del_eps[0]  del_eps[1]  ...
    //  u[0]        ...         ...
    //  u[1]        ...         ...
    //  ...
    matrix_type res;
    
    
    for (int i = 0; i < params_count; i++) {
        // не засоряем сетку по итерациям
        grid_type tmp_grid = grid;
        
        double increase_val = delta_u_coef * tmp_grid.elements[i].lambda;
        
        tmp_grid.elements[i].lambda = tmp_grid.elements[i].lambda + increase_val;
        
        auto b = make_global_vector(tmp_grid);
        auto A = make_global_matrix(tmp_grid);
        boundary_cond_2(tmp_grid.second_cond, b);
        boundary_cond_1(tmp_grid.first_cond, b, A);
        auto q = solve_equation(A, b);
        auto eps_increased = epsilon_vec(q, tmp_grid);
        
        // delta_epsilon для увеличенного параметра
        auto delta_eps_increased = delta_epsilon(eps_increased, true_epsilon);
        
        // отнимаем и нормируем вектора замеров, складываем в матрицу
        auto negative_delta_epsilon_v = scale_vector(delta_epsilon_v, -1);
        auto diff_v = vector_sum(delta_eps_increased, negative_delta_epsilon_v);
        res.push_back(scale_vector(diff_v, 1/increase_val));
    }
    
    return res;
}
