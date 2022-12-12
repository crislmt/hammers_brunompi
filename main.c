#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#define ten_days 864000
#define three_months 77760000
#define ten_minutes 36000

int N,I,t;
double  W, L, w, l, v, d;

typedef struct person{
    double x,y,vx,vy;
    InfectionStatus status;
    int time;
} Person;

typedef enum infectionStatus{
    SUSCEPTIBLE, INFECTED, IMMUNE, IN_CONTACT
} InfectionStatus;

double distance(int a_x, int a_y, int b_x, int b_y){
    double x=pow(a_x-b_x,2);
    double y=pow(a_y-b_y,2);
    return sqrt(x+y);
}


Person move(Person p){
    if((p.x + p.x*t)<=0){
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

    return p;
}
bool is_in_range(int x, int y){
    if(distance(p.x, p.y, x, y)<d)return 1;
    return 0;
}

void infection(){
    Person p;
    srand(time(NULL));
    int rk = MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    p.x= rand()%L;
    p.y= rand()%W;
    if(rk < I)p.status=INFECTED;
    p.time=0;

    int x_positions[N];
    int y_positions[N];
    int statuses[N];

    while (1)
    {
        p=move(p);

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

        int found = 0;
        for(int i = 0; i<N; i++){
            if(!found && is_in_range(x_positions[i], y_positions[i]) && statuses[i]==INFECTED && i!=rk){
                if(p.status==IN_CONTACT){
                    if(p.time>=ten_minutes){
                        p.status=INFECTED;
                        p.time=0;
                    }
                } else if(p.status==SUSCEPTIBLE){
                    p.status=IN_CONTACT;
                    p.time=0;
                }
                found = 1;
            }

        }
        if (!found){
            if(p.status==IN_CONTACT){
                p.status=SUSCEPTIBLE;
                p.time=0;
            }
        }
    }

    return null;
}

void initialize(int N, int I, double W, double L, double vx, double vy, double d, int t){

};

int main(){
    MPI_Init(NULL,NULL);

    initialize(argv[1], argv[2], argv[3], argv[4],argv[5], argv[6],argv[7], argv[8]);

    infection();

    MPI_Finalize();
}