#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
//-----//
#include <torch/torch.h>
#include <torch/script.h>
using namespace std;

//User defined area
#define Re 100
#define delta_t 0.001
#define record_per 100 //Record per every 0.1s is fine
#define tol_vel pow(10, -3) //=delta_t
#define tol_p pow(10, -12)
#define N 81

//Data exportatation 
string directory="./raw_data/Re_"+to_string(Re); 
string directory_ML="./ML/data/Re_"+to_string(Re); 
ofstream raw_data(directory+"/flowfield.dat");
ofstream step_time(directory+"/stepstoconvergence_timestep.dat");
ofstream input_v_fake(directory_ML+"/input_v_fake.dat");
ofstream input_u_fake(directory_ML+"/input_u_fake.dat");
ofstream output_p(directory_ML+"/output_p.dat");

//Grid
double p[N+1][N+1]={1.0}, v[N+1][N+1]={0.0}, u[N+1][N+1]={0.0}; //current
double v_fake[N+1][N+1]={0.0}, u_fake[N+1][N+1]={0.0}; //u* and v*
double p_1[N+1][N+1]={1.0}, v_1[N+1][N+1]={0.0}, u_1[N+1][N+1]={0.0}; //previous
double p_c[N][N]={1.0}, v_c[N][N]={0.0}, u_c[N][N]={0.0}, vel_c[N][N]={0.0}; //collocated grid
double pred_p[N+1][N+1]={1.0}, pred_p_c[N][N]={0.0};

#define PI 3.1415926
double h=1.0/(N-1.0);
double res_vel=1.0, res_p=1.0, dev_p=1.0;
int iteration, timestep=0;

void init_v(), init_u(), init_p();
void cal_fake_v(), cal_fake_u(), cal_p(), cal_v(), cal_u();
void guess_p();
bool velocity_not_converge(), pressure_not_converge();
void output(), output_train(), collocate(), write(double* a, int x, int y, char name[]), write_train(), peek();
double div_vel();

int main(){
    //test section
    init_u();
    init_v();
    init_p();

    while(velocity_not_converge()||timestep<100){
        cout<<"Time:    "<< timestep*delta_t<<" residual(vel):   "<<res_vel<<"    div_vel:    "<<div_vel();
        if(timestep%record_per==0){
            output();
            output_train();
            peek();
        }

        //step 1
        cal_fake_u();
        init_u();
        cal_fake_v();
        init_v();
        //step 2
        if(0)//(timestep>100)
            guess_p();
        cal_p();
        init_p();
        //step 3
        cal_u();
        init_u();
        cal_v();
        init_v();

        timestep++;
        //peek();
    }
    
    return 0;
}

void init_u(){
    for(int i=0; i<N; i++){
        /*Implement*/
        u[0][i]=-u[1][i]; u_fake[0][i]=-u_fake[1][i];
        u[N-1][i]=-u[N-2][i]; u_fake[N-1][i]=-u_fake[N-2][i];
        u[i][0]=-u[i][1]; u_fake[i][0]=-u_fake[i][1];
        u[i][N]=2.0-u[i][N-1]; u_fake[i][N]=2.0-u_fake[i][N-1];
    }
}
void init_v(){
    for(int i=0; i<N; i++){
        /*Implement*/
        v[0][i]=-v[1][i]; v_fake[0][i]=-v_fake[1][i];
        v[N][i]=-v[N-1][i]; v_fake[N][i]=-v_fake[N-1][i];
        v[i][0]=-v[i][1]; v_fake[i][0]=-v_fake[i][1];
        v[i][N-1]=-v[i][N-2]; v_fake[i][N-1]=-v_fake[i][N-2];
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
    for(int i=1; i<N-1; i++){
        for(int j=1; j<N; j++){
            v_p=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])*0.25;
            u_e=(u[i][j]+u[i+1][j])*0.5; u_w=(u[i-1][j]+u[i][j])*0.5; u_n=(u[i][j]+u[i][j+1])*0.5; u_s=(u[i][j]+u[i][j-1])*0.5;
            // C=....
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
            u_p=(u[i-1][j]+u[i][j]+u[i-1][j+1]+u[i][j+1])*0.25;
            v_e=(v[i+1][j]+v[i][j])*0.5; v_w=(v[i-1][j]+v[i][j])*0.5; v_n=(v[i][j]+v[i][j+1])*0.5; v_s=(v[i][j]+v[i][j-1])*0.5;
            // C=v_p*(v_e-v_w)/h+u_p*(v_n-v_s)/h
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
                u_fake_e=u_fake[i][j]; u_fake_w=u_fake[i-1][j]; v_fake_n=v[i][j]; v_fake_s=v[i][j-1];
                // term_poisson_left=....
                term_poisson_left=(u_fake_e-u_fake_w+v_fake_n-v_fake_s)/h/delta_t;
                p[i][j]=0.25*(p[i][j-1]+p[i][j+1]+p[i-1][j]+p[i+1][j]-term_poisson_left*h*h);
                dev_p+=fabs(p[i][j]-p_1);
            }
        }

        res_p=0;
        for(i=1;i<N+1; i++){
            for(j=1; j<N+1; j++){
                u_fake_e=u_fake[i][j]; u_fake_w=u_fake[i-1][j]; v_fake_n=v[i][j]; v_fake_s=v[i][j-1];
                // term_poisson_left=....
                term_poisson_left=(u_fake_e-u_fake_w+v_fake_n-v_fake_s)/h/delta_t;
                res_p+=fabs(p[i][j]-(0.25*(p[i][j-1]+p[i][j+1]+p[i-1][j]+p[i+1][j]-term_poisson_left*h*h)));
            }
        }

        if(0){//iteration%1000
            cout<<"Iteration:   "<<iteration<<" residual(p):   "<<res_p<<"  dev:    "<<dev_p<<endl;
            //init_p();
        } 
    }while(pressure_not_converge());//(pressure_not_converge());
    cout<<" Iteration:   "<<iteration<<" residual(p):   "<<res_p<<endl;

}
void cal_u(){
    res_vel=0;
    /*Implement*/
    for(int i=1; i<N-1;i++){
        for(int j=1; j<N; j++){
            u_1[i][j]=u[i][j];
            u[i][j]=u_fake[i][j]-delta_t*(p[i+1][j]-p[i][j])/h;
            res_vel+=fabs(u[i][j]-u_1[i][j]);
        }
    }
    
}
void cal_v(){
    /*Implement*/
    for(int i=1; i<N;i++){
        for(int j=1; j<N-1; j++){
            v_1[i][j]=v[i][j];
            v[i][j]=v_fake[i][j]-delta_t*(p[i][j+1]-p[i][j])/h;
            res_vel+=fabs(v[i][j]-v_1[i][j]);
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
            div+=fabs((u[i][j]-u[i-1][j])+(v[i][j]-v[i][j-1]));
        }
    }
    return div;
}
void collocate(){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            v_c[i][j]=0.5*(v[i][j]+v[i+1][j]);
            u_c[i][j]=0.5*(u[i][j]+u[i][j+1]);
            p_c[i][j]=0.25*(p[i][j]+p[i][j+1]+p[i+1][j]+p[i+1][j+1]);
            pred_p_c[i][j]=0.25*(pred_p[i][j]+pred_p[i][j+1]+pred_p[i+1][j]+pred_p[i+1][j+1]);

        }
    }
    return;
}

void output(){ //Tecplot format
    if(timestep==0){
        raw_data<<"VARIABLES=\"x\", \"y\", \"time\", \"u\", \"v\", \"p\""<<endl;
        raw_data<<"ZONE T=\"1\""<<endl;
        raw_data<<"F=POINT"<<endl;
        raw_data<<"I=81,J=81,K=(define)"<<endl;
        step_time<<"step iteration"<<endl;
    }
    else{
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                raw_data<<1-i*h<<" "<<j*h<<"    "<< timestep*delta_t<<" "<<v_c[i][j]<<"   "<<u_c[i][j]<<"   "<<p_c[i][j]<<endl; 
            }
        }
    }
    step_time<< timestep << "   " << iteration <<endl;

    
    
}
void write(double* a, int x, int y, char* name) {
	ofstream file(name);
	int i, j;
	for (i = 0; i < y; i++) {
		for (j = 0; j < x; j++) {
			file << *(a + i * x + j) << ",";
		}
		file << endl;
	}

	file.close();
	return;
}
void peek(){
    collocate();
	write(&u_c[0][0], N, N, (char*)"./peek/u.dat");
	write(&v_c[0][0], N, N, (char*)"./peek/v.dat");
	write(&p_c[0][0], N, N, (char*)"./peek/p.dat");
	write(&pred_p_c[0][0], N, N, (char*)"./peek/pred_p.dat");
    return;
}

void output_train(){
    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            input_u_fake<< u_fake[i][j]<<"  ";
            input_v_fake<< v_fake[i][j]<<"  ";
            output_p<< p[i][j]<<"  ";
        }
    }
    input_u_fake<< endl;
    input_v_fake<< endl;
    output_p<<endl;
}


void guess_p(){
    /*once*/
    // /home/enchou/git-repo
    // /mnt/c/Users/ENCHOU/Documents/git-repo
    torch::jit::script::Module model=torch::jit::load("/mnt/c/Users/ENCHOU/Documents/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/ML/01_model_jit.pth");
    double test[N+1][N+1]={0.0}; //這奇怪的bug，明明沒用到但不加就會core dump
    double u_st[0][(N+1)*(N+1)];
    auto options = torch::TensorOptions().dtype(torch::kFloat32);
    torch::Tensor x, out;
    vector<torch::jit::IValue> input;
    /*------*/


        /*(input-array) 2D into 1D*/
        for(int i=0; i<N+1; i++){
            for(int j=0; j<N+1; j++){
                u_st[0][i*(N+1)+j]=u_fake[i][j]; //test_u->u_fake
            }
        }

        /*(input) array into tensor into vector*/
        x=torch::from_blob(u_st, {1,(N+1)*(N+1)}, options);
        input.clear();
        input.push_back(x);

        /*input into output(tensor)*/
        out=model.forward(input).toTensor();

        /*(output)tensor into 2D-array*/
        for(int i=0; i<N+1; i++){
            for(int j=0; j<N+1; j++){
                pred_p[i][j]=out[0][i*(N+1)+j].item<float>(); //test_p->p
            }
        }


    
    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            //cout<<p[i][j]<<"   "; //test_p->p
        }
        //cout<<endl;
    }

}