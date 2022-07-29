#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

#define PI 3.1415926
#define N 81
#define delta_t 0.0001
#define tol_vel pow(10, -6)
#define tol_p pow(10, -6)
#define Re 100
double h=1.0/(N-1.0);

double p[N+1][N+1]={0.0}, u[N+1][N]={0.0}, v[N][N+1]={0.0}; //current
double u_fake[N+1][N]={0.0}, v_fake[N][N+1]={0.0};
double p_1[N+1][N+1]={0.0}, u_1[N+1][N]={0.0}, v_1[N][N+1]={0.0}; //previous
double res_vel=1.0, res_p=1.0, dev_p=1.0;
int iteration, timestep=0;


void init_u(), init_v(), init_p();
void cal_fake_u(), cal_fake_v(), cal_p(), cal_u(), cal_v();
bool velocity_not_converge(), pressure_not_converge();
void print(int);
void write(double* a, int x, int y, char name[]);
double div_vel();
char name1[] = "velocity.csv", name2[] = "u.csv", name3[] = "v.csv", name4[] = "p.csv";

int main(){
    //test section
    init_u();
    init_v();

    while(velocity_not_converge()||timestep<100){
        if(timestep%1==0){
            cout<<"Time:    "<< timestep*delta_t<<" residual(vel):   "<<res_vel<<"    div_vel:    "<<div_vel()<<endl;
        }
        timestep++;
        init_p();
        //step 1
        cal_fake_u();
        cal_fake_v();
        //step 2
        cal_p();
        //step 3
        cal_u();
        cal_v();
    }

    double vel[N][N]={0.0};

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            vel[i][j]=pow(u[i][j]*u[i][j]+v[i][j]*v[i][j],0.5);

        }
    }

	write(&vel[0][0], N, N, name1);
	write(&u[0][0], N, N, name2);
	write(&v[0][0], N, N, name3);
	write(&p[0][0], N, N, name4);
    return 0;
}

void init_u(){
    for(int i=0; i<N; i++){
        /*Implement*/
        u[0][i]=2-u[0][i];
        u[N][i]=-u[N-1][i];
        u[i][0]=-u[i][1];
        u[i][N-1]=-u[i][N-2];
    }
}

void init_v(){}

void init_p(){
    /*Implement*/
    for(int i=0; i<N+1; i++){
        p[0][i]=p[1][i];
        p[N][i]=p[N-1][i];
        p[i][0]=p[i][1];
        p[i][N]=p[i][N-1];
    }
}

void cal_fake_u(){
    /*Implement*/
    double C,D;
    for(int i=1; i<N; i++){
        for(int j=1; j<N-1; j++){
            // C=....
            // D=....
            u_fake[i][j]=(D-C)*delta_t+u[i][j];
        }
    }
}

void cal_fake_v(){
    /*Implement*/
    double C,D;
    for(int i=1; i<N; i++){
        for(int j=1; j<N-1; j++){
            // C=....
            // D=....
            v_fake[i][j]=(D-C)*delta_t+v[i][j];
        }
    }
}

void cal_p(){
    int i,j, te=2;
    double term_poisson_left;
    iteration=0;
    do{
        iteration++, dev_p=0;
        /*Implement*/
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                // term_poisson_left=....
                p[i][j]=0.25*(p[i][j-1]+p[i][j+1]+p[i-1][j]+p[i+1][j]-term_poisson_left*h*h);

            }
        }

        res_p

        if(iteration%100==0) //iteration%100
            cout<<"Iteration:   "<<iteration<<" residual(p):   "<<res_p<<"  dev:    "<<dev_p<<endl;

    }while(pressure_not_converge());
    cout<<"Iteration:   "<<iteration<<" residual(p):   "<<res_p<<endl;
}



void cal_u(){
        /*Implement*/
}

void cal_v(){
        /*Implement*/
}

bool pressure_not_converge(){
    if(res_p>tol_p)
        return true;
    else
        return false;
}
bool velocity_not_converge(){
    if(res_vel>tol_vel)
        return true;
    else
        return false;
}
void write(double* a, int x, int y, char name[]) {
	ofstream myfile(name);
	int i, j;
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
			myfile << *(a + i * x + j) << ",";
		}
		myfile << endl;
	}

	myfile.close();
	return;
}
double div_vel(){
    int i,j;
    double div=0.0;
    for(i=1; i<N-1; i++){
        for(j=1; j<N-1; j++){
            div+=fabs((u[i][j+1]-u[i][j-1])+(v[i-1][j]+v[i+1][j]));
        }
    }
    return div;
}