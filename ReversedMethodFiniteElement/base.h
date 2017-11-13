//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//

#ifndef MPAM_BASE_H
#define MPAM_BASE_H


#include <cstdio>
#include <vector>
using namespace std;

typedef vector<vector<double>> matrix_type;

typedef struct {
    double x;
    double y;
} point_type;

typedef struct {
    unsigned short index;
    point_type node_pos[2];
    unsigned short node_indecies[2];
} first_cond_type;

typedef struct {
    unsigned short index;
    point_type node_pos[2];
    unsigned short node_indecies[2];
    double edge_h;
} second_cond_type;

typedef struct {
    unsigned short omega_index;
    double lambda;
    double gamma;
    double h_x;
    double h_y;
    vector<double> f_values;
    vector<unsigned short> node_indecies;
    vector<point_type> node_pos;
} element_type;

typedef struct {
    unsigned short omega_count;
    unsigned short node_count;
    vector<element_type> elements;
    vector<first_cond_type> first_cond;
    vector<second_cond_type> second_cond;
    double i_amperage;
} grid_type;


double f_value(unsigned short omega_index, double x_pos, double y_pos, double i_amperage);
grid_type read_grid();

matrix_type append_local_to_global_matrix(matrix_type global, matrix_type local, vector<unsigned short> node_map);
matrix_type rigidity_matrix(grid_type grid, unsigned short omega_index);
matrix_type mass_matrix(grid_type grid, unsigned short omega_index);
matrix_type sum_sqrt_matrix(matrix_type a, matrix_type b, unsigned long size);
vector<double> b_vector(grid_type grid, unsigned short omega_index);
matrix_type make_global_matrix(grid_type grid);
vector<double> make_global_vector(grid_type grid);
void boundary_cond_1(vector<first_cond_type> conds, vector<double>  &b, matrix_type  &A);
void boundary_cond_2(vector<second_cond_type> conds, vector<double> &b);
void lu_decomp(matrix_type a, matrix_type &l, matrix_type &u);
vector<double> backward_solution(matrix_type u, vector<double> b);
vector<double> forward_solution(matrix_type l, vector<double> b);
vector<double> solve_equation(matrix_type global_matrix, vector<double> b_vector);
double calc_value(vector<double> q, grid_type grid, double x, double y);
vector<double> epsilon_vec(vector<double> q, grid_type grid);
double vector_length(vector<double> v_1);
vector<double> delta_epsilon(vector<double> eps_1, vector<double> eps_2);
matrix_type deriative_delta_epsilon(vector<double> u, grid_type grid, vector<double> delta_epsilon, vector<double> true_epsilon);
double mean(vector<double> v);
vector<double> reversed_vector_f(vector<double> delta_epsilon_with_true, matrix_type d_delta_epsilons, vector<double> params, vector<double> weight, double alpha);
vector<double> get_params(grid_type grid);
matrix_type reversed_matrix_a(matrix_type d_delta_epsilons, vector<double> weight, double alpha);
vector<double> setStartParams(grid_type &grid);
vector<double> vector_sum(vector<double> v_1, vector<double> v_2);
grid_type setParams(grid_type grid, vector<double> u);
vector<double> get_delta_u(grid_type grid, vector<double> true_epsilon, vector<double> weight, double alpha);
vector<double> epsilon_noized_vec(vector<double> q, grid_type grid);
void print_vector(vector<double> v, string str);
vector<double> get_weights(vector<double> eps_true);
vector<double> scale_vector(vector<double> v, double s);
double calc_alpha(vector<double> u_step, double d_eps_l);
double calc_deviation(vector<double> u_step);

// размер локальной матрицы
unsigned short const local_matrix_size = 4;
// эпсилон для сравнения чисел
double const epsilon = 0.00000001;
// положение источника
double const emitter_x_pos = 0;
double const emitter_y_pos = 6;
// сила тока в источнике
double const amperage = 1;

// количество и положения приемников
unsigned short const reciver_count = 7, N = 7;

// количество восстанавливаемых параметров
unsigned short const params_count = 7, M = 7;

const double reciver_y_pos[7] = {1, 1, 3, 3, 3, 5, 5};
const double reciver_x_pos[7] = {3, 9, 2, 6, 10, 3, 9};

const double sigma_start[7] = {1, 1, 1, 1, 1, 1, 1};

// коэффициен шага численной производной
const double delta_u_coef = 0.005;

// параметр релаксации
const double beta_start = 1;
const int max_iter_count = 10;
// параметр регуляризации
const double gamma_val = 1-0E2;

const double diff_size = 0.01;

const double delta_u_min_diff = 1E-15;
const double noize_percentage[7] = {0, 0.01, 0.01, 0.01, 0, 0, 0.01};

#endif //MPAM_BASE_H
