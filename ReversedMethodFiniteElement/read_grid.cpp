//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//
#include "base.h"
#include <fstream>
#include <cmath>

/// Считывается сетка в формате
// кол-во элементов
// далее повторяется для каждого элемента
    // индекс элемента
    // лямбда
    // гамма
    // 4 точки (глобальный номер узла, x, y)

grid_type read_grid() {
    grid_type grid;

    fstream grid_file ("./base.txt", ios_base::in);

    grid_file >> grid.omega_count;
    grid_file >> grid.node_count;

    for (short i = 0; i < grid.omega_count; ++i) {
        element_type grid_element {};
        grid.elements.push_back(grid_element);

        grid_file >> grid.elements[i].omega_index;
        grid_file >> grid.elements[i].lambda;
        grid_file >> grid.elements[i].gamma;



        for (int j = 0; j < 4; ++j) {
            unsigned short index;
            double x_pos, y_pos;

            grid_file >> index;
            grid_file >> x_pos;
            grid_file >> y_pos;

            point_type p{
                    .x = x_pos,
                    .y = y_pos
            };

            grid.elements[i].f_values.push_back(f_value(i,x_pos, y_pos, amperage));
            grid.elements[i].node_indecies.push_back(index);
            grid.elements[i].node_pos.push_back(p);
        }

        // Длинна действительна только при проходе против часовой с юго-западного угла
        grid.elements[i].h_x = abs(grid.elements[i].node_pos[0].x - grid.elements[i].node_pos[3].x);
        grid.elements[i].h_y = abs(grid.elements[i].node_pos[0].y - grid.elements[i].node_pos[3].y);
    }

    // краевые условия
    unsigned short cond_count;
    grid_file >> cond_count;

    // первого рода
    for (unsigned short k = 0; k < cond_count; ++k) {
        first_cond_type border_cond {
                .index = k
        };

        grid_file >> border_cond.node_indecies[0];
        grid_file >> border_cond.node_pos[0].x;
        grid_file >> border_cond.node_pos[0].y;

        grid_file >> border_cond.node_indecies[1];
        grid_file >> border_cond.node_pos[1].x;
        grid_file >> border_cond.node_pos[1].y;

        grid.first_cond.push_back(border_cond);
    }

    // второго рода
    grid_file >> cond_count;
    for (unsigned short k = 0; k < cond_count; ++k) {
        second_cond_type border_cond {
                .index = k
        };
        grid_file >> border_cond.node_indecies[0];
        grid_file >> border_cond.node_pos[0].x;
        grid_file >> border_cond.node_pos[0].y;

        grid_file >> border_cond.node_indecies[1];
        grid_file >> border_cond.node_pos[1].x;
        grid_file >> border_cond.node_pos[1].y;
        grid_file >> border_cond.edge_h;

        grid.second_cond.push_back(border_cond);
    }

    grid.i_amperage = amperage;

    return grid;
}

vector<double> get_params(grid_type grid) {
    vector<double> res;
    res.push_back(grid.elements[0].lambda);
    res.push_back(grid.elements[1].lambda);
    res.push_back(grid.elements[2].lambda);
    
    return res;
}

vector<double> setStartParams(grid_type &grid) {
    vector<double> res;
    
    for (int i = 0; i < grid.elements.size(); i++) {
        res.push_back(sigma_start[i]);
        grid.elements[i].lambda = sigma_start[i];
    }

    return res;
}

grid_type setParams(grid_type grid, vector<double> u) {
    grid_type res_grid = grid;
    for (int i = 0; i < res_grid.elements.size(); i++) {
        res_grid.elements[i].lambda = u[i];
    }
    
    return res_grid;
}
