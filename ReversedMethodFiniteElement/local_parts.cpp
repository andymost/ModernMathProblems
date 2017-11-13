//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//

#include "base.h"

double g_1[4][4] = {
        {2, -2, 1, -1},
        {-2, 2, -1, 1},
        {1, -1, 2, -2},
        {-1, 1, -2, 2}
};

double g_2[4][4] = {
        {2, 1, -2, -1},
        {1, 2, -1, -2},
        {-2, -1, 2, 1},
        {-1, -2, 1, 2}
};

double c[4][4] = {
        {4, 2, 2, 1},
        {2, 4, 1, 2},
        {2, 1, 4, 2},
        {1, 2, 2, 4}
};


// матрица жесткости - G
matrix_type rigidity_matrix(grid_type grid, unsigned short omega_index) {
    matrix_type result;
    vector<double> tmp;

    double lambda_by_6 = grid.elements[omega_index].lambda/6;
    double h_x_by_h_y = grid.elements[omega_index].h_x / grid.elements[omega_index].h_y;
    double h_y_by_h_x = grid.elements[omega_index].h_y / grid.elements[omega_index].h_x;

    for (int i = 0; i < local_matrix_size; ++i) {
        result.push_back(tmp);

        for (int j = 0; j < local_matrix_size; ++j) {
            double item_val = lambda_by_6 * ( h_y_by_h_x * g_1[i][j] + h_x_by_h_y * g_2[i][j]);
            result[i].push_back(item_val);
        }
    }

    return  result;
}


matrix_type mass_matrix(grid_type grid, unsigned short omega_index) {
    matrix_type result;
    vector<double> tmp;
    double koef = grid.elements[omega_index].h_x * grid.elements[omega_index].h_y/36;

    for (int i = 0; i < local_matrix_size; ++i) {
        result.push_back(tmp);

        for (int j = 0; j < local_matrix_size; ++j) {
            double item_val = grid.elements[omega_index].gamma * koef *c[i][j];
            result[i].push_back(item_val);
        }
    }

    return  result;
}

vector<double> b_vector(grid_type grid, unsigned short omega_index) {
    vector<double> result;
    double koef = (grid.elements[omega_index].h_x * grid.elements[omega_index].h_y)/36;
    double sum;

    for (int i = 0; i < local_matrix_size; ++i) {
        sum = 0;
        for (int j = 0; j < local_matrix_size; ++j) {
            sum += c[i][j] * koef * f_value(
                    omega_index,
                    grid.elements[omega_index].node_pos[j].x,
                    grid.elements[omega_index].node_pos[j].y,
                    grid.i_amperage
            );
        }
        result.push_back(sum);
    }

    return result;
}

double psi_1(double x, double x_next, double h_x, double y, double y_next, double h_y) {
    return ((x_next - x)/h_x) * ((y_next - y)/h_y);
}

double psi_2(double x, double x_cur, double h_x, double y, double y_next, double h_y) {
    return ((x - x_cur)/h_x) * ((y_next - y)/h_y);
}

double psi_3(double x, double x_next, double h_x, double y, double y_cur, double h_y) {
    return ((x_next - x)/h_x) * ((y - y_cur)/h_y);
}

double psi_4(double x, double x_cur, double h_x, double y, double y_cur, double h_y) {
    return ((x - x_cur)/h_x) * ((y - y_cur)/h_y);
}

double calc_value(vector<double> q, grid_type grid, double x, double y) {
    unsigned short element_index = 0;
    bool isFinded = false;

    for (int i = 0; i < grid.omega_count && !isFinded; ++i) {
        bool inX = x + epsilon >=  grid.elements[i].node_pos[0].x && x - epsilon <=  grid.elements[i].node_pos[1].x;
        bool inY = y + epsilon >=  grid.elements[i].node_pos[0].y && y - epsilon <=  grid.elements[i].node_pos[2].y;

        if(inX && inY) {
            isFinded = true;
            element_index = i;
        }
    }

    if(!isFinded) {
        throw ;
    }
    
    double
        theta_0 = q[grid.elements[element_index].node_indecies[0]],
        theta_1 = q[grid.elements[element_index].node_indecies[1]],
        theta_2 = q[grid.elements[element_index].node_indecies[2]],
        theta_3 = q[grid.elements[element_index].node_indecies[3]];

    double x_cur = grid.elements[element_index].node_pos[0].x;
    double y_cur = grid.elements[element_index].node_pos[0].y;
    double x_next = x_cur + grid.elements[element_index].h_x;
    double y_next = y_cur + grid.elements[element_index].h_y;
    double h_x = grid.elements[element_index].h_x;
    double h_y = grid.elements[element_index].h_y;

    double res =
            theta_0 * psi_1(x, x_next, h_x, y, y_next, h_y) +
            theta_1 * psi_2(x, x_cur, h_x, y, y_next, h_y) +
            theta_2 * psi_3(x, x_next, h_x, y, y_cur, h_y) +
            theta_3 * psi_4(x, x_cur, h_x, y, y_cur, h_y);

    return res;
}
