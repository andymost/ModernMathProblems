//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//
#include "base.h"
#include <cmath>
#include <iostream>
#include <string>

matrix_type sum_sqrt_matrix(matrix_type a, matrix_type b, unsigned long size) {
    matrix_type result;
    vector<double> tmp_vector;

    // выделяем память
    for (int i = 0; i < size; ++i) {
        result.push_back(tmp_vector);
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            result[i].push_back(a[i][j] + b[i][j]);
        }
    }

    return result;
}

double mean(vector<double> v) {
    int s = (int)v.size();
    double res = 0;
    
    for (int i = 0; i < s; i++) {
        res += v[i];
    }
    
    return res/s;
}

double vector_length(vector<double> v_1) {
    double acc = 0;
    for (int i = 0; i < v_1.size(); i ++) {
        acc += pow(v_1[i], 2);
    }
    
    return sqrt(acc);
}

vector<double> vector_sum(vector<double> v_1, vector<double> v_2) {
    vector<double> res;
    for (int i = 0; i < v_1.size(); i++) {
        res.push_back(v_1[i] + v_2[i]);
    }
    
    return res;
}

void print_vector(vector<double> v, string str) {
    cout<<"vector: "<<str<<endl;
    for (int i = 0; i < v.size(); i++) {
        cout<<v[i]<<' ';
    }
    cout<<endl;
}

vector<double> scale_vector(vector<double> v, double s) {
    vector<double> res;
    for (int i = 0; i < v.size(); i++) {
        res.push_back(v[i] * s);
    }
    
    return res;
}

double calc_deviation(vector<double> v) {
    double v_mean = mean(v);
    double res = 0;
    
    for (int i = 0; i < v.size(); i++) {
        res += pow(v_mean - v[i], 2);
    }
    
    return res;
}
