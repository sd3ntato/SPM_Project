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
#include "matplotlibcpp.h"

//#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

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
  float l = 0.995;
  float nabla = 0.1;
  int i = 0;

  vector<double> errors_norms{};

  ESN n = ESN(Nr = Nr, Nu = Nu, Ny = Ny);
  float **W = n.W.m;
  float **Win = n.Win.m;
  float **Wout = zeros(Ny, Nr).m;
  float **Wold = zeros(Ny, Nr).m;
  float **P = (eye(Nr) * (1 / nabla)).m;
  float **Pold = (eye(Nr) * (1 / nabla)).m;

  float *x = zeros(1, Nr).m[0];
  float *x_old = zeros(1, Nr).m[0];
  float *k = zeros(1, Nr).m[0];
  float *z = zeros(1, Nr).m[0];
  float *y = zeros(1, 4).m[0];
  float *u;
  float *d;

  float k_den;
  float s;

  while (i < n_samples)
  {

    u = dataset_n.m[i];
    d = dataset.m[i + 1];

    for(int i=0; i<Nr; i++){
      x[i] = tanh( sum( 1, Nr, W[i], x_old ) + sum( 1,Nu, Win[i], u ) );
    }

    for(int i=0; i<Nr; i++){
      z[i] = sum( 1, Nr, P[i], x );
    }
    
    for(int i=0; i<Ny; i++){
      y[i] = sum( 1, Nr, Wout[i], x);
    }

    k_den = l + sum( 1, Nr, x , z ) ;

    for(int i=0; i<Nr; i++){
      k[i] = z[i] / k_den ;
    }

    for(int i=0; i<Ny; i++){
      for(int j=0; j<Nr; j++){
        Wout[i][j] = Wold[i][j] + ( d[i] - y[i] ) * k[j] ;
      }
    }

    for(int i=0; i<Nr; i++){
      for(int j=0; j<Nr; j++){
        P[i][j] = ( Pold[i][j] - k[i] * z[j] ) * 1/l;
      }
    }

    for(int i=0; i<Ny; i++){
      for(int j=0; j<Nr; j++){
        Wold[i][j] = Wout[i][j];
      }
    }

    for(int i=0; i<Nr; i++){
      for(int j=0; j<Nr; j++){
        Pold[i][j] = P[i][j] ;
      }
    }

    for(int i=0; i<Nr; i++){
      x_old[i] = x[i];
    }

    s=0;
    for(int i=0; i<Ny; i++){
      s += pow( d[i] - y[i], 2 );
    }
    cout<< sqrt(s) << " " << flush;
    errors_norms.push_back( s );


    i++;
  }

  cout << endl;
  plt::plot(errors_norms);
  plt::show();

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
