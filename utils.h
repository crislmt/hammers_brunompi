//
// Created by Chris on 12/12/2022.
//

#ifndef NSDS_MPI_UTILS_H
#define NSDS_MPI_UTILS_H

typedef enum infectionStatus{
    SUSCEPTIBLE, INFECTED, IMMUNE, IN_CONTACT
} InfectionStatus;

typedef struct person{
    int x,y,vx,vy;
    InfectionStatus status;
    int time;
} Person;

double distance(int, int , int, int);

#endif //NSDS_MPI_UTILS_H
