#include <iostream>
#include <random>

using namespace std;

class ESN {  
  public:      
    float rho; 
    int Nr;  // Number recurrent units
    int Nu; // Number inputs
    int Ny; // Number outputs
    float r_density; // density of recurrent matrix
    float** Win; // input-to-reservoir matrix
    float** W; //recurrent matrix
    float** x; //state


    ESN(int Nr = 100, int Ny = 1, int Nu = 1, float rho = 0.9, float r_d = 0.1) { // Constructor with parameters
      this->Nr = Nr;
      this->Ny = Ny;
      this->Nu = Nu;
      this->rho = rho;
      this->r_density = r_d;
    }
};

int main() {

    ESN n =  ESN() ;
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd());
    uniform_real_distribution<> dis(1.0, 2.0);
    for(int i=0; i<10; i++){
        std::cout<< dis(gen) <<" \n";
    }

    cout << "\nftt!\n";
    return 0;
}