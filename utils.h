#ifndef NSDS_MPI_UTILS_H
#define NSDS_MPI_UTILS_H

typedef enum infectionStatus{
    SUSCEPTIBLE, INFECTED, IMMUNE, IN_CONTACT
} InfectionStatus;

typedef struct person{
    int id;
    int x,y;
    float vx,vy;
    InfectionStatus status;
    int time;
} Person;
double distance(float, float , float, float);

#endif //NSDS_MPI_UTILS_H
