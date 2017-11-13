//
// Created by Кочнев Андрей Владимирович on 15/10/2017.
//
#include "base.h"

double f_value(unsigned short omega_index, double x_pos, double y_pos, double i_amperage) {
    bool inX = x_pos - epsilon < emitter_x_pos && x_pos + epsilon > emitter_x_pos;
    bool inY = y_pos - epsilon < emitter_y_pos && y_pos + epsilon > emitter_y_pos;

    if(inX && inY) {
        return i_amperage;
    }

    return 0;
}
