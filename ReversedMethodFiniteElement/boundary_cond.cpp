//
// Created by Кочнев Андрей Владимирович on 16/10/2017.
//

#include "base.h"

double b_cond_1 (double x, double y, unsigned short index) {
//    switch (index) {
//        case 0:
//            return 2;
//        case 1:
//            return 1.8 + 0.1 * x;
//        default:{}
//    }
//
//    throw ;

    return  0;
}


double theta_cond_2 (double x, double y, unsigned short index) {
    switch (index) {
        case 0:
        case 1:
            return 1;
        case 2:
            return 0;
        default:{}
    }

    throw ;
}


void boundary_cond_1(vector<first_cond_type> conds, vector<double>  &b, matrix_type  &A) {
    for (unsigned short i = 0; i < conds.size(); i ++) {
        first_cond_type cond = conds[i];
        unsigned short
            first_node_index = cond.node_indecies[0],
            second_node_index = cond.node_indecies[1];

        double
            first_b = b_cond_1(cond.node_pos[0].x, cond.node_pos[0].y, i),
            second_b = b_cond_1(cond.node_pos[1].x, cond.node_pos[1].y, i);

        b[first_node_index] = first_b;
        b[second_node_index] = second_b;

        for (int j = 0; j < A[first_node_index].size(); ++j) {
            A[first_node_index][j] = 0;
        }

        for (int j = 0; j < A[second_node_index].size(); ++j) {
            A[second_node_index][j] = 0;
        }

        A[first_node_index][first_node_index] = 1;
        A[second_node_index][second_node_index] = 1;
    }
}

void boundary_cond_2(vector<second_cond_type> conds, vector<double> &b) {
    for (unsigned short i = 0; i < conds.size(); i ++) {
        second_cond_type cond = conds[i];
        unsigned short
                first_node_index = cond.node_indecies[0],
                second_node_index = cond.node_indecies[1];

        double
            theta_1 = theta_cond_2(cond.node_pos[0].x, cond.node_pos[0].y, i),
            theta_2 = theta_cond_2(cond.node_pos[1].x, cond.node_pos[1].y, i);

        double first_b = (cond.edge_h / 6) * (2 * theta_1 + theta_2);
        double second_b = (cond.edge_h / 6) * (2 * theta_2 + theta_1);


        b[first_node_index] += first_b;
        b[second_node_index] += second_b;
    }
}
