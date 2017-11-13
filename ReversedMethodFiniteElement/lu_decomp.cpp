//
// Created by Кочнев Андрей Владимирович on 16/10/2017.
//


#include "base.h"

void lu_decomp(matrix_type a, matrix_type &l, matrix_type &u) {
    int i = 0, j = 0, k = 0;
    for (i = 0; i < a.size(); i++)
    {
        for (j = 0; j < a.size(); j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < a.size(); j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}

vector<double> backward_solution(matrix_type u, vector<double> b) {
    int size = (int)b.size();
    vector<double> result(b.size(), 0);
    double sum = 0;
    for (int i = size - 1; i >= 0; --i) {
        for (int j = size - 1; j > i; --j) {
            sum += result[j] * u[i][j];
        }
        result[i] = (b[i] - sum)/u[i][i];
        sum = 0;
    }

    return result;
}

vector<double> forward_solution(matrix_type l, vector<double> b) {
    unsigned long size = b.size();
    vector<double> result;
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < i; ++j) {
            sum += result[j] * l[i][j];
        }
        result.push_back((b[i] - sum)/l[i][i]);
        sum = 0;
    }
    return result;
}
