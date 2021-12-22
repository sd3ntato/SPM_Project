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
#include "matplotlibcpp.h"

/*
* this was the absolute first, sequential version of the main file
*/

//#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;


using namespace std;

Matrix_wrapper read_dataset(string filename, int n_samples);

int main()
{
  int n_samples = 1000;

  cout << "reading dataset...";
  Matrix_wrapper dataset = read_dataset("./BTCUSDT-1m-data.csv", n_samples);
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
  Matrix_wrapper u, d, x, psi, zeta, k, x1, xt, p1, p2, kt, delta, zeta_t, u_l, d_l, old_Wout, t1;
  while (i < n_samples)
  {

    u_l = dataset_n.get_line(i);
    u = u_l.transpose();   // shape(Nu,1)
    d_l = dataset.get_line(i + 1);
    d = d_l.transpose(); //shape(Ny,1)

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
  }

  cout << endl;
  errors_norms.erase(errors_norms.begin());
  plt::plot(errors_norms);
  plt::show();

  cout << "\nftt!\n";
  return 0;
}

