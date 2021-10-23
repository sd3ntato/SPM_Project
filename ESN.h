 
class ESN {  
  public:      
    ESN(int Nr = 100, int Ny = 1, int Nu = 1, float rho = 0.9, float r_d = 0.1); 
    float rho; 
    int Nr;  // Number recurrent units
    int Nu; // Number inputs
    int Ny; // Number outputs
    float r_density; // density of recurrent matrix
    float** Win; // input-to-reservoir matrix
    float** W; //recurrent matrix
    float** x; //state
};
