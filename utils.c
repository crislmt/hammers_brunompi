#include <math.h>

double distance(float a_x, float a_y, float b_x, float b_y){
    float x=powf(a_x-b_x,2);
    float y=powf(a_y-b_y,2);
    return (float)sqrtf(x+y);
}
