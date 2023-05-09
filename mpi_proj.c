#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include "utils.h"
#include <time.h>

#define TEN_DAYS 864000
//#define THREE_MONTHS 77760000
#define THREE_MONTHS 864000
#define TEN_MINUTES 600
#define DAY 86400
//Global Variables
int N,I,t, w, l, true_N, chunksize, n_country;
long  W, L;
float v, d;
bool large_grid;

//Function headers
bool is_in_range(Person p, int x, int y);
void move (Person* p);
void infection();
void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , float speed, float distance, int timestep);
int compute_country_id(float x, float y);
void printWorld(bool* world);
void resetWorld(bool* world);
void putInfectedOnTheWorld(bool* world, int* x_positions, int* y_positions, int num_infected);
void removeInfectedFromTheWorld(bool* world, int* x_positions, int* y_positions, int num_infected);
void computeNewInfected(bool* reduced_world, Person* process_people);
void computeNewInfectedOptLargeDistance(Person* process_people, int* infected_x_positions, int* infected_y_positions, int global_infected);
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
        printf("Wrong number of arguments\n");
        return 0;
    }

    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    initialize(atoi(argv[1]),atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atof(argv[7]), atof(argv[8]), atoi(argv[9]));
    infection();

    MPI_Finalize();
}


void infection(){
    int day = 0;
    int global_timer = DAY;
    int nation_infected[n_country];
    int nation_susceptible[n_country];
    int nation_immune[n_country];
    int nation_pot_infected[n_country];
    int nation_infected_total[n_country];
    int nation_susceptible_total[n_country];
    int nation_immune_total[n_country];
    int nation_pot_infected_total[n_country];
    int rk;
    int sz;
    clock_t initTime, endTime;
    if(W*L > 1000000000)large_grid=true;
    bool* reduced_world;
    if(!large_grid)reduced_world=(bool*)malloc((W*L)*sizeof(bool));
    else reduced_world=NULL;
    resetWorld(reduced_world);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
    int infected_per_process[sz];
    int displ[sz];
    int rcvcount[sz];
    int global_infected=0;
   
    Person process_people[chunksize];
    srand(time(NULL) + rk+1);
    
    for(int i=0; i<chunksize; i++){
        process_people[i].x= rand()%W;
        process_people[i].y= rand()%L;
        process_people[i].id = rk*chunksize+i;
        process_people[i].vx =  (((float)rand()/(float)(RAND_MAX)) * (v*2))-v;
        process_people[i].vy = (((float)rand()/(float)(RAND_MAX)) * (v*2))-v;
        //if(process_people[i].vx>=0&&process_people[i].vx<1)process_people[i].vx=1;
        //if(process_people[i].vx<0&&process_people[i].vx>-1)process_people[i].vx=-1;
        //if(process_people[i].vy>=0&&process_people[i].vy<1)process_people[i].vy=1;
        //if(process_people[i].vy<0&&process_people[i].vy>-1)process_people[i].vy=-1;


        if(process_people[i].id < I)process_people[i].status=INFECTED;
        else process_people[i].status=SUSCEPTIBLE;
        process_people[i].time=0;
        if(process_people[i].id < I) process_people[i].time=rand()%TEN_DAYS;
        
        if(process_people[i].id>=true_N){
            process_people[i].x= -1;
            process_people[i].y= -1;
            process_people[i].vx = 0;
            process_people[i].vy = 0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    while (true)
    {   
        if(global_timer==DAY)initTime=clock();
        global_timer -= t;
        int local_infected = 0;
        global_infected = 0;
        int* local_infected_x_positions = NULL;
        int* local_infected_y_positions= NULL;
        for(int i =0;  i<chunksize; i++){
            if(process_people[i].id<true_N){
                move(&process_people[i]);
                if(global_timer==0){
                  process_people[i].vx =  (((float)rand()/(float)(RAND_MAX)) * (v*2))-v;
                  process_people[i].vy = (((float)rand()/(float)(RAND_MAX)) * (v*2))-v;
                }
                process_people[i].time+=t;
                if(process_people[i].status==INFECTED){
                  if(process_people[i].time >= TEN_DAYS){
                    process_people[i].status=IMMUNE;
                    process_people[i].time=0;
                  }
                  else{
                    local_infected +=1;
                    if(local_infected==1){
                      local_infected_x_positions = (int*)malloc(sizeof(int));
                      local_infected_y_positions = (int*)malloc(sizeof(int));
                    }
                    else{
                      local_infected_x_positions = (int*)realloc(local_infected_x_positions, sizeof(int) * local_infected);
                      local_infected_y_positions = (int*)realloc(local_infected_y_positions, sizeof(int) * local_infected);
                    }
                    local_infected_x_positions[local_infected-1]=(int)process_people[i].x;
                    local_infected_y_positions[local_infected-1]=(int)process_people[i].y;
                  }
                }
                if(process_people[i].status==IMMUNE && process_people[i].time>=THREE_MONTHS){
                    process_people[i].status=SUSCEPTIBLE;
                    process_people[i].time=0;
                }
            }
        }
       
       if(rk==0){
            displ[0]=0;
            displ[1]=local_infected;
            rcvcount[0]=local_infected;
            global_infected=local_infected;
            //getting the number of infected from each process to compute the display and receive count for the gatherv and preparing vector or right size
            for (unsigned remote_process = 1; remote_process < sz; ++remote_process) {
                int recvd_infected;
                //could have used an MPI_reduce to compute global_filtered_lines_number but the filtered lines are need to compute the display and rcvcount too
                MPI_Recv(&recvd_infected, 1, MPI_INT,
                         remote_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                global_infected+=recvd_infected;
                if(remote_process<sz-1)displ[remote_process+1] = displ[remote_process]+recvd_infected;
                rcvcount[remote_process]=recvd_infected;            }
            
        }
        else {
          MPI_Send(&local_infected, 1, MPI_INT,0,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(&global_infected, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int global_infected_x_positions[global_infected];
        int global_infected_y_positions[global_infected];
        MPI_Gatherv(local_infected_x_positions, local_infected, MPI_INT, &global_infected_x_positions, rcvcount ,displ, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(local_infected_y_positions, local_infected, MPI_INT, &global_infected_y_positions, rcvcount ,displ, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_infected_x_positions, global_infected, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_infected_y_positions, global_infected, MPI_INT, 0, MPI_COMM_WORLD);
        if(d*d<global_infected && !large_grid){
           putInfectedOnTheWorld(reduced_world, global_infected_x_positions, global_infected_y_positions, global_infected);
           computeNewInfected(reduced_world, process_people);
           removeInfectedFromTheWorld(reduced_world, global_infected_x_positions, global_infected_y_positions, global_infected); 
        }
        else if(d*d>=global_infected || large_grid){
          computeNewInfectedOptLargeDistance(process_people, global_infected_x_positions, global_infected_y_positions, global_infected);
        } 
        
        if(global_timer==0){
            day++;
            if(rk==0)printf("day %d:\n", day);
            endTime = clock();
            double cpu_time_used = ((double) (endTime - initTime)) / CLOCKS_PER_SEC;
            
            global_timer = DAY;
            for(int i=0;i<n_country;i++){
                  nation_infected[i]=0;
                  nation_infected_total[i]=0;
                  nation_susceptible[i]=0;
                  nation_susceptible_total[i]=0;
                  nation_immune[i]=0;
                  nation_immune_total[i]=0;
                  nation_pot_infected[i]=0;
                  nation_pot_infected_total[i]=0;
            }
            for(int person=0; person< chunksize; person++){
                if(process_people[person].id< true_N){
                    Person p=process_people[person];
                    int country_id=compute_country_id(p.x,p.y);
                    if(p.status==INFECTED){
                        nation_infected[country_id]++;
                        nation_pot_infected[country_id]++;
                    }
                    else if(p.status==IMMUNE){
                        nation_immune[country_id]++;
                    }
                    else {
                        nation_susceptible[country_id]++;
                        if(p.status==IN_CONTACT)nation_pot_infected_total[country_id]++;
                    }
                    
                }
            }
            MPI_Reduce(nation_immune,&nation_immune_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(nation_susceptible,&nation_susceptible_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(nation_infected,&nation_infected_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(nation_pot_infected,&nation_pot_infected_total,n_country,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            int end = 0;
            if(rk==0){
                int infected = 0;
                int immune = 0;
                int susceptible = 0;
                for(int i=0; i<n_country; i++){
                    printf("for nation %d, there are %d infected people, %d susceptible people, %d immune people \n",
                               i, nation_infected_total[i],nation_susceptible_total[i],nation_immune_total[i]);
                    end += nation_pot_infected_total[i]; 
                    infected+=nation_infected_total[i];
                    susceptible+=nation_susceptible_total[i];
                    immune+=nation_immune_total[i];
                }
                printf("in total there are %d infected people, %d susceptible people, %d immune people \n", infected, susceptible, immune);
            }
            if(rk==0)printf("Execution time: %f seconds\n", cpu_time_used);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(end==0)return;
        }
        free(local_infected_x_positions);
        free(local_infected_y_positions);
    }

}

void move(Person* p){

    if((p->x + (p->vx)*t)<0){
        p->vx=-p->vx;
        p->x=0;
    }
    else if(p->x + p->vx*t>=W){
        p->vx=-p->vx;
        p->x=W-1;
    }
    else p->x=p->x+p->vx*t;

    if((p->y + p->vy*t)<0){
        p->vy=-p->vy;
        p->y=0;
    }
    else if(p->y + p->vy*t>=L){
        p->vy=-p->vy;
        p->y=L-1;
    }
    else p->y=p->y+p->vy*t;

}

bool is_in_range(Person p, int x, int y){
    if(distance(p.x, p.y, x, y)<d) return true;
    return false;

}

void initialize(int number_of_persons, int infected_persons, int rectangular_width, int rectangular_length, int country_width, int country_length , float speed, float distance, int timestep){
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

    while (N%sz!=0){
        N++;
    }
    chunksize=N/sz;
    
}

void printWorld(bool* world){
  if(large_grid)return;
  for(int i = 0; i<L; i++){
    for(int j = 0; j<W; j++){
      printf("%d ", world[i*W+j]);
    }
    printf("\n");
  }
}

void resetWorld(bool* world){
  if(large_grid)return;
  for (int i = 0; i < L; i++){
        for(int j = 0; j<W; j++){
          world[i*W+j]=false;
      }
    }
}

void putInfectedOnTheWorld(bool* world, int* x_positions, int* y_positions, int num_infected){
  if(large_grid)return;
  unsigned pos;
  for(int i = 0; i<num_infected; i++){
    pos=(int)y_positions[i]*W+(int)x_positions[i];
    world[pos]=true;
  }
}

void removeInfectedFromTheWorld(bool* world, int* x_positions, int* y_positions, int num_infected){
  if(large_grid)return;
  unsigned pos;
  for(int i = 0; i<num_infected; i++){
    pos=(int)y_positions[i]*W+(int)x_positions[i];
    world[pos]=false;
  }
}

void computeNewInfected(bool *reduced_world, Person *process_people)
{
  for (int person = 0; person < chunksize; person++)
  {
    if (process_people[person].id < true_N)
    {
          int x = (int)process_people[person].x;
          int y = (int)process_people[person].y;
          if (process_people[person].status == INFECTED || process_people[person].status == IMMUNE)
          {
                continue;
          }
          else
          {
                bool found = false;
                for (int xdelta = -d; xdelta < d; xdelta++)
                {
                    for (int ydelta = -d; ydelta < d; ydelta++)
                    {
                      if (!found)
                      {
                        if (x + xdelta >= 0 && x + xdelta < W && y + ydelta >= 0 && y + ydelta < L)
                        {
                          if (reduced_world[(y + ydelta) * W + (x + xdelta)] == true)
                          {
                            if (process_people[person].status == IN_CONTACT)
                            {
                              if (process_people[person].time >= TEN_MINUTES)
                              {
                                process_people[person].status = INFECTED;
                                process_people[person].time = 0;
                              }
                            }
                            else if (process_people[person].status == SUSCEPTIBLE)
                            {
                              process_people[person].status = IN_CONTACT;
                              process_people[person].time = 0;
                            }
                            found = true;
                            break;
                          }
                        }
                      }
                    } 
                    if(found)break;
                }
                if (!found)
                {
                    if (process_people[person].status == IN_CONTACT)
                    {
                        process_people[person].status = SUSCEPTIBLE;
                        process_people[person].time = 0;
                    }
                }
          }
    }
  }
}

void computeNewInfectedOptLargeDistance(Person* process_people, int* infected_x_positions, int* infected_y_positions, int global_infected){
  for(int person=0; person<chunksize; person++){    
      if(process_people[person].id<true_N){
              if(process_people[person].status==INFECTED || process_people[person].status == IMMUNE){
                continue;
              }
              else{
                bool found = false;
                for(int i=0; i<global_infected; i++){
                  if(!found){
                    if(is_in_range(process_people[person], infected_x_positions[i], infected_y_positions[i])){
                      if(process_people[person].status==IN_CONTACT){
                            if(process_people[person].time>=TEN_MINUTES){
                                process_people[person].status=INFECTED;
                                process_people[person].time=0;
                              }
                            } 
                            else if(process_people[person].status==SUSCEPTIBLE){
                              process_people[person].status=IN_CONTACT;
                              process_people[person].time=0;
                            }
                            found = true;
                            break;
                      }
                  }
                }
                if (!found){
                    if(process_people[person].status==IN_CONTACT){
                        process_people[person].status=SUSCEPTIBLE;
                        process_people[person].time=0;
                    }
                }
              }
    }
  }
}

int compute_country_id(float x, float y){
    int country_on_x=W/w;
    int countr_y = (int) ceil(y/l)-1;
    int countr_x = (int) ceil(x/w)-1;
    int country_id = (int) (countr_y * country_on_x + countr_x);
    if (country_id <0) return 0;
    if (country_id > (n_country-1))return(n_country-1);
    else return country_id;
}



