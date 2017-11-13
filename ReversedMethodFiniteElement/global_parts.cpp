//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//

#include "base.h"

matrix_type init_global_matrix(unsigned long node_count) {
    vector<vector<double>> result;
    vector<double> tmp(node_count, 0);

    for (int i = 0; i < node_count; ++i) {
        result.push_back(tmp);
    }

    return result;
}

matrix_type append_local_to_global_matrix(matrix_type global, matrix_type local, vector<unsigned short> node_map) {
    for (int i = 0; i < local_matrix_size; ++i) {
        for (int j = 0; j < local_matrix_size; ++j) {
            unsigned short i_global = node_map[i];
            unsigned short j_global = node_map[j];
            global[i_global][j_global] += local[i][j];
        }
    }

    return global;
}

matrix_type make_global_matrix(grid_type grid) {
    matrix_type A = init_global_matrix(grid.node_count);
    matrix_type m, g, a;

    for (int i = 0; i < grid.omega_count; ++i) {
        m = mass_matrix(grid, i);
        g = rigidity_matrix(grid, i);
        a = sum_sqrt_matrix(m, g, local_matrix_size);
        A = append_local_to_global_matrix(A, a, grid.elements[i].node_indecies);
    }

    return A;
}

vector<double> append_local_to_global_vector(
        vector<double> global,
        vector<double> local,
        vector<unsigned short> node_map
) {
    for (int i = 0; i < local_matrix_size; ++i) {
        unsigned short g_index = node_map[i];
        global[g_index] += local[i];
    }

    return global;
}

vector<double> make_global_vector(grid_type grid) {
    vector<double> b(grid.node_count, 0), b_local;

    for (int i = 0; i < grid.omega_count; ++i) {
        b_local = b_vector(grid, i);
        b = append_local_to_global_vector(b, b_local, grid.elements[i].node_indecies);
    }

    return b;
}

vector<double> solve_equation(matrix_type global_matrix, vector<double> b_vector){
    unsigned long size = global_matrix.size();
    matrix_type l_matrix = init_global_matrix(size);
    matrix_type u_matrix = init_global_matrix(size);

    lu_decomp(global_matrix, l_matrix, u_matrix);

    vector<double> x_1 = forward_solution(l_matrix, b_vector);
    return backward_solution(u_matrix, x_1);
}
