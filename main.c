#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#include "utils.h"

#define TEN_DAYS 864000
#define THREE_MONTHS 77760000
#define TEN_MINUTES 36000
#define DAY 86400

//Global Variables
int N,I,t, W, L, w, l, v, d, true_N, chunksize, n_country;

//Function headers
bool is_in_range(Person, int, int);
void move (Person* p);
void infection();
void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , int speed, int distance, int timestep);
int compute_country_id(int x, int y);
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

    if (argc != 10) {
        std::cout << "Expected 1 input, got " << argc - 1<< std::endl;
        return 0;
    }

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    initialize(atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
    infection();

    MPI_Finalize();
}


void infection(){

    int x_positions[N];
    int y_positions[N];
    int statuses[N];
    int local_x_positions[chunksize];
    int local_y_positions[chunksize];
    int local_statuses[chunksize];
    int global_timer = DAY;
    int nation_infected[n_country];
    int nation_susceptible[n_country];
    int nation_immune[n_country];
    int nation_infected_total[n_country];
    int nation_susceptible_total[n_country];
    int nation_immune_total[n_country];
    int rk;
    int sz;

    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);



    Person process_people[chunksize];
    srand(time(0));
    for(int i=0; i<chunksize; i++){
        process_people[i].x= rand()%L;
        process_people[i].y= rand()%W;
        process_people[i].id = rk*chunksize+i;
        process_people[i].vx = v;
        process_people[i].vy = v;
        if(process_people[i].id < I) process_people[i].status=INFECTED;
        else process_people[i].status=SUSCEPTIBLE;
        process_people[i].time=0;

        if(process_people[i].id>=true_N){
            process_people[i].x= -1;
            process_people[i].y= -1;
            process_people[i].vx = 0;
            process_people[i].vy = 0;
        }

    }
    while (true)
    {


        global_timer -= t;

        for(int i =0;  i<chunksize; i++){
            if(process_people[i].id<true_N){

                move(&process_people[i]);

                process_people[i].time+=t;
                if(process_people[i].status==INFECTED && process_people[i].time >= TEN_DAYS){
                    process_people[i].status=IMMUNE;
                    process_people[i].time=0;
                }
                if(process_people[i].status==IMMUNE && process_people[i].time>=THREE_MONTHS){
                    process_people[i].status=SUSCEPTIBLE;
                    process_people[i].time=0;
                }
            }

        }

        for(int i=0; i<chunksize; i++){
            local_x_positions[i]=process_people[i].x;
            local_y_positions[i]=process_people[i].y;
            local_statuses[i]=process_people[i].status;
        }

        MPI_Gather(&local_x_positions, chunksize, MPI_INT, x_positions, chunksize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&local_y_positions, chunksize, MPI_INT, y_positions,chunksize, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&local_statuses, chunksize, MPI_INT, statuses, chunksize, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(x_positions, N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(y_positions, N, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(statuses, N, MPI_INT, 0, MPI_COMM_WORLD);


        for(int person=0; person<chunksize; person++){
            if(process_people[person].id<true_N){
                bool found = false;
                for(int i = 0; i<N; i++){
                    if(!found && is_in_range(process_people[person], x_positions[i], y_positions[i]) && statuses[i]==INFECTED && i!=process_people[person].id){
                        if(process_people[person].status==IN_CONTACT){
                            if(process_people[person].time>=TEN_MINUTES){
                                printf("infected\n");
                                process_people[person].status=INFECTED;
                                process_people[person].time=0;
                            }
                        } else if(process_people[person].status==SUSCEPTIBLE){
                            process_people[person].status=IN_CONTACT;
                            process_people[person].time=0;
                        }
                        found = true;
                    }

                }
                if (!found){
                    if(process_people[person].status==IN_CONTACT){
                        process_people[person].status=SUSCEPTIBLE;
                        process_people[person].time=0;
                    }
                }
            }
            if(global_timer==0){
                global_timer = DAY;
                for(int i=0;i<n_country;i++){
                    nation_infected[i]=0;
                    nation_susceptible[i]=0;
                    nation_immune[i]=0;
                }
                for(int person=0;person< chunksize;person++){
                    if(process_people[person].id< true_N){
                        Person p=process_people[person];
                        int country_id= compute_country_id(p.x,p.y);
                        if(p.status==INFECTED){
                            nation_infected[country_id]++;
                        }
                        else if(p.status==SUSCEPTIBLE || p.status==IN_CONTACT){
                            nation_susceptible[country_id]++;
                        }
                        else{
                            nation_immune[country_id]++;
                        }
                    }
                }
                MPI_Reduce(&nation_immune,&nation_immune_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&nation_susceptible,&nation_susceptible_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&nation_infected,&nation_infected_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                if(rk==0){
                    for(int i=0; i<n_country; i++){
                        printf("for nation %d, there are %d infected people, %d susceptible people, %d immune people \n",
                               i, nation_infected_total[i],nation_susceptible_total[i],nation_immune_total[i]);
                    }
                    int cont;
                    printf("pres 0 to stop simulation or 1 to simulate another day\n");
                    //scanf("%d",&cont);
                    //if(cont==0)return;
                }


            }
        }

    }
}
void move(Person* p){

    if((p->x + (p->vx)*t)<=0){
        p->vx=-p->vx;
        p->x=0;
    }
    else if(p->x + p->vx*t>=L){
        p->vx=-p->vx;
        p->x=L;
    }
    else p->x=p->x+p->vx*t;

    if((p->y + p->vy*t)<=0){
        p->vy=-p->vy;
        p->y=0;
    }
    else if(p->y + p->vy*t>=W){
        p->vy=-p->vy;
        p->y=W;
    }
    else p->y=p->y+p->vy*t;

}

bool is_in_range(Person p, int x, int y){

    if(distance(p.x, p.y, x, y)<d) return true;
    return false;

}

void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , int speed, int distance, int timestep){


    N=number_of_persons;
    I=infected_persons;
    W=rectangular_width;
    L=rectangular_length;
    w=country_width;
    l=country_length;
    v=speed;
    d=distance;
    t=timestep;
    true_N=N;
    n_country = W/w * L/l;
    int rk;
    int sz;
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);



    /*
    if(rk==0){
        while (N%sz!=0){
            N++;
        }
        chunksize=N/sz;
    }
    MPI_BROADCAST(chunksize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    */

    while (N%sz!=0){
        N++;
    }
    chunksize=N/sz;

    //Setting the initial status and position for this person, the first I persons are marked as infected

    srand(time(NULL));

    /*
    for(int i=0; i<chunksize; i++){
        p.x= rand()%L;
        p.y= rand()%W;
        p.id = rk*chunksize+i;
        if(p.id < I) p.status=INFECTED;
        p.time=0;
    }
     */
};
int compute_country_id(int x, int y){
    int country_on_x=W/w;
    int country_on_y=L/l;
    int cy, cx;
    if(country_on_x==1)cy=0;
    else cy=y/(country_on_x-1);
    if(country_on_y==1)cx=0;
    else cx=x/(country_on_y-1);
    return (cy*country_on_x)+cx;
}





