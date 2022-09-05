#include <iostream>
# define debug
# ifdef debug 
# define FLAG cout<<"Successfully print out Hello wo\sfdg
rld"<<endl
#else
# define FLAG 
#endif
using namespace std;

int main(){
    cout<<"Hello world"<<endl;
    FLAG;
    return 0;
}