//目前程式只能用cmake compile，沒辦法直接用vs code
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <torch/torch.h>
#include <torch/script.h>
//#include <../libtorch/include/torch/script.h>
using namespace std;
#define N 81
double test_u[N+1][N+1]={0.0}, test_p[N+1][N+1]={0.0};
void guess();

int main(){/*
    torch::jit::script::Module model=torch::jit::load("/home/enchou/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/ML/01_model_jit.pth");
    torch::Tensor x = torch::randn({ 1,6724 });  
	vector<torch::jit::IValue> input;  
	input.push_back(x);  
	auto out = model.forward(input);  

	cout << "out: " << out << endl;  */
    guess();
    return 0;
}

void guess(){
    /*once*/
    torch::jit::script::Module model=torch::jit::load("/home/enchou/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/ML/01_model_jit.pth");
    double test[N+1][N+1]={0.0}; //這奇怪的bug，明明沒用到但不加就會core dump
    double u_st[0][(N+1)*(N+1)];
    auto options = torch::TensorOptions().dtype(torch::kFloat32);
    torch::Tensor x, out;
    vector<torch::jit::IValue> input;
    /*------*/

    for (int test_n=0; test_n<1; test_n++){

        /*(input-array) 2D into 1D*/
        for(int i=0; i<N+1; i++){
            for(int j=0; j<N+1; j++){
                u_st[0][i*(N+1)+j]=test_u[i][j]; //test_u->u_fake
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
                test_p[i][j]=out[0][i*(N+1)+j].item<float>(); //test_p->p
            }
        }

    }
    
    for(int i=0; i<N+1; i++){
        for(int j=0; j<N+1; j++){
            cout<<test_p[i][j]<<"   "; //test_p->p
        }
        cout<<endl;
    }

}
