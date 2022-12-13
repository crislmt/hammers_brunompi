#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#include "utils.h"

#define ten_days 864000
#define three_months 77760000
#define ten_minutes 36000

//Global Variables
int N,I,t, W, L, w, l, v, d;
Person p;

//Function headers
bool is_in_range(Person, int, int);
void move ();
void infection();
void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , int speed, int distance, int timestep);

/*
 * Receives as command line arguments:
 *
 * N = argv[1] = number of persons
 * I = argv[2] = number of persons that are initially infected
 * W = argv[3] = width of the rectangular area where persons move (in meters)
 * L = argv[4] = length of the rectangular area where persons move (in meters)
 * w = argv[5] = width of each country (in meters)
 * l = argv[6] = length of each country (in meters)
 * v = argv[7] = moving speed on the x-axis for a person
 * d = argv[8] = maximum spreading distance (in meters): a susceptible person that remains closer than d to at least one infected person becomes infected
 * t = argv[9] = time step (in seconds): the simulation recomputes the position and status (susceptible,infected, immune, in_contact ) of each person with a temporal granularity of t second
 */
int main(int argc, char *argv[]){

    MPI_Init(NULL, NULL);

    initialize(atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
    infection();

    MPI_Finalize();
}


void infection(){
    
    int x_positions[N];
    int y_positions[N];
    int statuses[N];

    while (true)
    {
        move(&p);

        p.time+=t;
        if(p.status==INFECTED && p.time >= ten_days){
            p.status=IMMUNE;
            p.time=0;
        }
        if(p.status==IMMUNE && p.time>=three_months){
            p.status=SUSCEPTIBLE;
            p.time=0;
        }

        MPI_GATHER(&x_positions, 1, MPI_INT, p.x, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_GATHER(&y_positions, 1, MPI_INT, p.y, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_GATHER(&statuses, 1, MPI_INT, p.status, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_BROADCAST(x_positions, N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_BROADCAST(y_positions, N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_BROADCAST(statuses, N, MPI_INT, 0, MPI_COMM_WORLD);

        bool found = false;
        for(int i = 0; i<N; i++){
            if(!found && is_in_range(p, x_positions[i], y_positions[i]) && statuses[i]==INFECTED && i!=rk){
                if(p.status==IN_CONTACT){
                    if(p.time>=ten_minutes){
                        p.status=INFECTED;
                        p.time=0;
                    }
                } else if(p.status==SUSCEPTIBLE){
                    p.status=IN_CONTACT;
                    p.time=0;
                }
                found = true;
            }

        }
        if (!found){
            if(p.status==IN_CONTACT){
                p.status=SUSCEPTIBLE;
                p.time=0;
            }
        }
    }
}
void move(){

    if((p.x + (p.x)*t)<=0){
        p.vx=-p.vx;
        p.x=0;
    }
    else if(p.x + p.x*t>=L){
        p.vx=-p.vx;
        p.x=L;
    }
    else p.x=p.x+p.vx*t;

    if((p.y + p.y*t)<=0){
        p.vy=-p.vy;
        p.y=0;
    }
    else if(p.y + p.y*t>=W){
        p.vy=-p.vy;
        p.y=W;
    }
    else p.y=p.y+p.vy*t;

}

bool is_in_range(Person p, int x, int y){

    if(distance(p.x, p.y, x, y)<d) return true;
    return false;

}

void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , int speed, int distance, int timestep){

    N=number_of_persons;
    I=infected_persons;
    W=rectangular_width;
    I=rectangular_length;
    w=country_width;
    l=country_length;
    v=speed;
    d=distance;
    t=timestep;

    //Setting the initial status and position for this person, the first I persons are marked as infected

    srand(time(NULL));
    int rk = MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    p.x= rand()%L;
    p.y= rand()%W;
    if(rk < I) p.status=INFECTED;
    p.time=0;
    
};





