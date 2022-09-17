#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
//-----//
#include <torch/torch.h>
#include <torch/script.h>
using namespace std;
double u_fake[2][2]={1.0}, pred_p[2][2]={1.0};

int main(){

    torch::jit::script::Module model=torch::jit::load("/mnt/c/Users/ENCHOU/Documents/git-repo/Fractional-Step-FDM-Staggered-Lid-Driven-Cavity-/ML/test_2/model_jit.pth");
    double test[2][2]={0.0}; //這奇怪的bug，明明沒用到但不加就會core dump
    double u_st[1][4]={{1.345, 2.345, 3.876, 4.344}};
    auto options = torch::TensorOptions().dtype(torch::kFloat64);
    torch::Tensor x, out;
    vector<torch::jit::IValue> input;
    /*------*/
    
    /*(input) array into tensor into vector*/
    x=torch::from_blob(u_st, {1,4}, options);

    input.clear();
    input.push_back(x);

    /*input into output(tensor)*/
    out=model.forward(input).toTensor();
    /*(output)tensor into 2D-array*/
    
    cout<< "out[0][j]:   "<<endl;
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            pred_p[i][j]=out[0][i*(2)+j].item<double>(); //test_p->p
        }
    }
    return 0;
    
}