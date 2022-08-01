#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

#define PI 3.1415926
#define N 81
#define delta_t 0.001
#define tol_vel pow(10, -8)
#define tol_p pow(10, -6)
#define Re 100
double h=1.0/(N-1.0);

double p[N+1][N+1]={1.0}, u[N+1][N]={0.0}, v[N][N+1]={0.0}; //current
double u_fake[N+1][N]={0.0}, v_fake[N][N+1]={0.0};
double p_1[N+1][N+1]={1.0}, u_1[N+1][N]={0.0}, v_1[N][N+1]={0.0}; //previous
double p_c[N][N]={1.0}, u_c[N][N]={0.0}, v_c[N][N]={0.0}, vel_c[N][N]={0.0}; //collocated grid
double res_vel=1.0, res_p=1.0, dev_p=1.0;
int iteration, timestep=0;


void init_u(), init_v(), init_p();
void cal_fake_u(), cal_fake_v(), cal_p(), cal_u(), cal_v();
bool velocity_not_converge(), pressure_not_converge();
void output(), collocate(), write(double* a, int x, int y, char name[]);
double div_vel();
char name1[] = "velocity.dat", name2[] = "u.dat", name3[] = "v.dat", name4[] = "p.dat";
char name5[]="time_iterations.dat", name6[]="training_input.dat", name7[]="training_output.dat";
ofstream myfile5(name5), myfile6(name6);

int main(){

    myfile5<< "time, iterations"<<endl;
    //test section
    init_u();
    init_v();
    init_p();

    while(velocity_not_converge()||timestep<100){
        if(timestep%1==0){
            cout<<"Time:    "<< timestep*delta_t<<" residual(vel):   "<<res_vel<<"    div_vel:    "<<div_vel(); 
        }
        timestep++;
        //step 1
        cal_fake_u();
        init_u();
        cal_fake_v();
        init_v();
        //step 2
        cal_p();
        init_p();
        //step 3
        cal_u();
        init_u();
        cal_v();
        init_v();
        if(timestep%10==0)
            output();

    }
    
	write(&u[0][0], N+1, N, name2);
    //output();
	myfile5.close();
    return 0;
}

void init_u(){
    for(int i=0; i<N; i++){
        /*Implement*/
        u[0][i]=2-u[1][i]; u_fake[0][i]=2-u_fake[1][i];
        u[N][i]=-u[N-1][i]; u_fake[N][i]=-u_fake[N-1][i];
        u[i][0]=-u[i][1]; u_fake[i][0]=-u_fake[i][1];
        u[i][N-1]=-u[i][N-2]; u_fake[i][N-1]=-u_fake[i][N-2];
    }
}
void init_v(){
    for(int i=0; i<N; i++){
        /*Implement*/
        v[0][i]=-v[1][i]; v_fake[0][i]=-v_fake[1][i];
        v[N-1][i]=-v[N-2][i]; v_fake[N-1][i]=-v_fake[N-2][i];
        v[i][0]=-v[i][1]; v_fake[i][0]=-v_fake[i][1];
        v[i][N]=-v[i][N-1]; v_fake[i][N]=-v_fake[i][N-1];
    }
    }
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
    double v_p;
    double u_e, u_w, u_n, u_s;
    double C,D;
    for(int i=1; i<N; i++){
        for(int j=1; j<N-1; j++){
            v_p=(v[i-1][j]+v[i][j]+v[i-1][j+1]+v[i][j+1])*0.25;
            u_e=(u[i][j]+u[i][j+1])*0.5; u_w=(u[i][j-1]+u[i][j])*0.5; u_n=(u[i-1][j]+u[i][j])*0.5; u_s=(u[i+1][j]+u[i][j])*0.5;
            // C=u_p*(u_e-u_w)/h+v_p*(u_n-u_s)/h
            C=u[i][j]*(u_e-u_w)/h+v_p*(u_n-u_s)/h;
            // D=....
            D=1.0/Re*(u_e+u_w+u_n+u_s-4*u[i][j])/(h*h);
            u_fake[i][j]=(D-C)*delta_t+u[i][j];
        }
    }
}
void cal_fake_v(){
    /*Implement*/
    double u_p;
    double v_e, v_w, v_n, v_s;
    double C,D;
    for(int i=1; i<N; i++){
        for(int j=1; j<N-1; j++){
            u_p=(u[i][j-1]+u[i][j]+u[i+1][j-1]+u[i+1][j])*0.25;
            v_e=(v[i][j]+v[i][j+1])*0.5; v_w=(v[i][j-1]+v[i][j])*0.5; v_n=(v[i-1][j]+v[i][j])*0.5; v_s=(v[i+1][j]+v[i][j])*0.5;
            // C=....
            C=u_p*(v_e-v_w)/h+v[i][j]*(v_n-v_s)/h;
            // D=....
            D=1.0/Re*(v_e+v_w+v_n+v_s-4*v[i][j])/(h*h);
            v_fake[i][j]=(D-C)*delta_t+v[i][j];
        }
    }
}

void cal_p(){
    int i,j;
    double term_poisson_left;
    double u_fake_e, u_fake_w, v_fake_n, v_fake_s;
    iteration=0;
    do{
        iteration++, dev_p=0; double p_1=0;
        /*Implement*/
        for(i=1; i<N+1; i++){
            for(j=1; j<N+1; j++){
                p_1=p[i][j];
                u_fake_e=u_fake[i][j]; u_fake_w=u_fake[i][j-1]; v_fake_n=v[i-1][j]; v_fake_s=v[i][j];
                // term_poisson_left=....
                term_poisson_left=(u_fake_e-u_fake_w+v_fake_n-v_fake_s)/h/delta_t;
                p[i][j]=0.25*(p[i][j-1]+p[i][j+1]+p[i-1][j]+p[i+1][j]-term_poisson_left*h*h);
                dev_p+=fabs(p[i][j]-p_1);
            }
        }
        //init_p();

        res_p=0;
        for(i=1;i<N+1; i++){
            for(j=1; j<N+1; j++){
                u_fake_e=u_fake[i][j]; u_fake_w=u_fake[i][j-1]; v_fake_n=v[i-1][j]; v_fake_s=v[i][j];
                // term_poisson_left=....
                term_poisson_left=(u_fake_e-u_fake_w+v_fake_n-v_fake_s)/h/delta_t;
                res_p+=fabs(p[i][j]-(0.25*(p[i][j-1]+p[i][j+1]+p[i-1][j]+p[i+1][j]-term_poisson_left*h*h)));
            }
        }

        if(0 ){//iteration%1000
            cout<<"Iteration:   "<<iteration<<" residual(p):   "<<res_p<<"  dev:    "<<dev_p<<endl;
            //init_p();
        } 


    }while(pressure_not_converge());
    cout<<" Iteration:   "<<iteration<<" residual(p):   "<<res_p<<endl;
	myfile5<< timestep*delta_t<< ","<<iteration<<endl;
}
void cal_u(){
    /*Implement*/
    res_vel=0;
    for(int i=1; i<N;i++){
        for(int j=1; j<N-1; j++){
            u_1[i][j]=u[i][j];
            u[i][j]=u_fake[i][j]-delta_t*(p[i][j+1]-p[i][j])/h;
            res_vel=fabs(u[i][j]-u_1[i][j]);
        }
    }
}
void cal_v(){
    /*Implement*/
    for(int i=1; i<N-1;i++){
        for(int j=1; j<N; j++){
            v_1[i][j]=v[i][j];
            v[i][j]=v_fake[i][j]-delta_t*(p[i][j]-p[i+1][j])/h;
            res_vel=fabs(v[i][j]-v_1[i][j]);
        }
    }
    
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
double div_vel(){
    int i,j;
    double div=0.0;
    for(i=1; i<N-1; i++){
        for(j=1; j<N-1; j++){
            div+=fabs((u[i][j]-u[i][j-1])+(v[i-1][j]-v[i][j]));
        }
    }
    return div;
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
void collocate(){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            u_c[i][j]=u[i][j]+u[i+1][j];
            v_c[i][j]=v[i][j]+v[i][j+1];
            vel_c[i][j]=pow(u_c[i][j]*u_c[i][j]+v_c[i][j]*u_c[i][j],0.5);
            p_c[i][j]=0.25*(p[i][j]+p[i][j+1]+p[i+1][j]+p[i+1][j+1]);

        }
    }
    return;
}
void output(){
    collocate();
	write(&vel_c[0][0], N, N, name1);
	write(&u_c[0][0], N, N, name2);
	write(&v_c[0][0], N, N, name3);
	write(&p_c[0][0], N, N, name4);
    return;
}

