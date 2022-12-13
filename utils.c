//
// Created by Chris on 12/12/2022.
//

#include <math.h>

double distance(int a_x, int a_y, int b_x, int b_y){
    double x=pow(a_x-b_x,2);
    double y=pow(a_y-b_y,2);
    return sqrt(x+y);
}
