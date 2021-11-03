#include <iostream>
#include <random>
#include "ESN.h"
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <assert.h>
#include "matplotlib-cpp/matplotlibcpp.h"
#include "omp.h"
#include <Eigen/Dense>

//BTCUSDT-1m-data.csv
using namespace Eigen;

using namespace std;
namespace plt = matplotlibcpp;

Matrix_wrapper read_dataset(string filename, int n_samples);

int main()
{
  int n_samples = 10;

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
  Matrix_wrapper P = eye(Nr + 1) * (1 / nabla);
  n.Wout = zeros(Ny, Nr + 1);
  Matrix_wrapper y = zeros(4, 1);
  Matrix_wrapper u, d, x, psi, zeta, k;
  while (i < n_samples)
  {

    u = dataset_n.get_line(i).transpose();   // shape(Nu,1)
    d = dataset.get_line(i + 1).transpose(); //shape(Ny,1)

    x = vstack(n.compute_state(u), ones(1, 1)); //shape(Nr,1)
    //print_matrix(n.x);

    y = n.compute_output(); //shape(Ny,1)
    //print_matrix(y);

    psi = d - y; // shape(Ny,1)
    //print_matrix(psi); //errore

    zeta = P | x; //shape(Nr,1)
    //print_matrix(zeta);

    float k_den = lambda + (x.transpose() | zeta).to_float();
    //cout << k_den <<endl;

    k = zeta / k_den; //shape(Nr,1)
    //print_matrix(k);

    P = (P - (k | zeta.transpose())) * (1 / lambda); //shape(Nr,Nr)
    //print_matrix(P);

    n.Wout = n.Wout + (psi | k.transpose()); // shape(Ny,Nr)
    //print_matrix(n.Wout);

    errors_norms.push_back(norm(psi));
    cout << errors_norms[i] << " " << flush;
    i++;

    free_matrices({u, d, x, y, psi, zeta, k});
  }

  cout << endl;
  //plt::plot(errors_norms);
  //plt::show();

  cout << "\nftt!\n";
  return 0;
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
      lines_read++;
    }
    myfile.close();
  }
  return dataset;
}
