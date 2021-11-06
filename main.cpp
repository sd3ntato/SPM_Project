#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef io_stream
#define io_stream
#include <iostream>
#endif

#include "ESN.h"
#include <fstream>
#include <string>
#include <sstream>
#include "math.h"

//#include "matplotlib-cpp/matplotlibcpp.h"
//namespace plt = matplotlibcpp;

using namespace std;

Matrix_wrapper read_dataset(string filename, int n_samples);
float sum(int start, int stop, float* v1, float* v2);

int main()
{
  int n_samples = 1000;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("BTCUSDT-1m-data.csv", n_samples);
  Matrix_wrapper dataset_n = normalize(dataset);
  cout << "dataset read" << endl;

  int Nr;
  string input;
  cout << "insert number of recurrent neurons: \n";
  cin >> input;
  Nr = stoi(input);

  int Nu = 4;
  int Ny = 4;
  float lambda = 0.995;
  float nabla = 0.1;
  int i = 0;

  vector<double> errors_norms{};

  ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
  float **W = n.W.m;
  float **Win = n.Win.m;
  float **Wout = zeros(Ny, Nr).m;
  float **Wold = zeros(Ny, Nr).m;
  float **P = (eye(Nr) * (1 / nabla)).m;
  float **Pold = zeros(Ny, Nr).m;

  float *x = zeros(1, Nr).m[0];
  float *x_old = zeros(1, Nr).m[0];
  float *k = zeros(1, Nr).m[0];
  float *zeta = zeros(1, Nr).m[0];
  float *y = zeros(1, 4).m[0];
  float *u;
  float *d;

  while (i < n_samples)
  {

    u = dataset_n.m[i];
    d = dataset.m[i + 1];

    for(int i=0; i<Nr; i++){
      x[i] = tanh( sum( 1, Nr, W[i], x_old ) + sum( 1,Nu, Win[i], u ) );
    }
    /*

    x1 = n.compute_state(u);
    x = vstack(x1, ones(1, 1)); //shape(Nr,1)
    //print_matrix(n.x);

    y = n.compute_output(); //shape(Ny,1)
    //print_matrix(y);

    psi = d - y; // shape(Ny,1)
    //print_matrix(psi); //errore

    zeta = P | x; //shape(Nr,1)
    //print_matrix(zeta);

    xt = x.transpose();
    t1 = xt | zeta;
    float k_den = lambda + t1.to_float();
    //cout << k_den <<endl;

    k = zeta / k_den; //shape(Nr,1)
    //print_matrix(k);

    zeta_t = zeta.transpose();
    p1 = k | zeta_t;
    p2 = P - p1 ;
    free_matrices({P});
    P = p2 * (1 / lambda); //shape(Nr,Nr)
    //print_matrix(P);

    kt = k.transpose();
    delta = psi | kt ; // shape (Ny,1) (1, Nr)
    old_Wout = copy(n.Wout);
    free_matrices({n.Wout});
    n.Wout = old_Wout + delta; // shape(Ny,Nr)
    //print_matrix(n.Wout);

    errors_norms.push_back(norm(psi));
    cout << errors_norms[i] << " " << flush;
    i++;

    free_matrices({u, d, x, psi, zeta, k, x1, xt, p1, p2, kt, delta, zeta_t, u_l, d_l, old_Wout, t1});
  */
  }

  cout << endl;
  //plt::plot(errors_norms);
  //plt::show();

  cout << "\nftt!\n";
  return 0;
}

float sum(int start, int stop, float* v1, float* v2){
  float s = 0.0;
  for(int i=start;i<stop;i++){
    s += v1[i] * v2[i];
  }
  return s;
}


Matrix_wrapper read_dataset(string filename, int n_samples)
{
  string line;
  ifstream myfile(filename);
  Matrix_wrapper dataset = Matrix_wrapper(nullptr, 0, 0); //container for the final dataset
  int lines_read = 0;
  if (myfile.is_open())
  {
    getline(myfile, line); // discard the first line as it contains intestation.
    while (getline(myfile, line) && lines_read < n_samples + 1)
    {                          //get a line from the csv
      istringstream iss(line); //turn it into this thing
      string s;
      string *ss = new string[12]; // each line is divided into 12 tokens (we are interested in token 1 to 5)
      Matrix_wrapper n;
      float *numbers = new float[4]; // keep only ohlc values
      int i = 0;
      while (getline(iss, s, ','))
      {
        ss[i] = s; //store the tokens in the apposite array
        i++;
      }
      for (int i = 0; i < 5; i++)
      {
        numbers[i] = stof(ss[i + 1]);
      }
      n = from_array(numbers, 4);   //dichiarazione va spostata fuori
      dataset = vstack(dataset, n); // put the collected line into the dataset
      delete[] ss;
      free_matrices({n});
      lines_read++;
    }
    myfile.close();
  }
  return dataset;
}
