#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

#ifndef io_stream
#define io_stream
#include <iostream>
#endif

#ifndef prot_queue_h
#define prot_queue_h
#include "prot_queue.h"
#endif

#ifndef pool_h
#define pool_h
#include <pool.h>
#endif

#include "ESN.h"
#include "math.h"
#include "matplotlibcpp.h"

//#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

Matrix_wrapper read_dataset(string filename, int n_samples);

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
  float **Wout = zeros(Ny, Nr + 1).m;
  float **Wold = zeros(Ny, Nr + 1).m;
  float **P = (eye(Nr + 1) * (1 / nabla)).m;
  float **Pold = (eye(Nr + 1) * (1 / nabla)).m;

  float *x = zeros(1, Nr + 1).m[0];
  float *x_rec = zeros(1, Nr + 1).m[0];
  float *x_in = zeros(1, Nr + 1).m[0];
  float *x_old = zeros(1, Nr + 1).m[0];
  x_old[Nr] = 1.0;
  float *k = zeros(1, Nr + 1).m[0];
  float *z = zeros(1, Nr + 1).m[0];
  float *y = zeros(1, 4).m[0];
  float *u;
  float *d;

  float k_den;
  float s;

  Pool p(7);

  while (i < n_samples)
  {

    u = dataset_n.m[i];
    d = dataset.m[i + 1];

    for (int i = 0; i < Nr; i++)
    {
      p.submit({new Dot_task(0, Nr, W[i], x_old, &x_rec[i])});
      // x_rec[i] = dot(0, Nr, W[i], x_old) ;

      p.submit({new Dot_task(0, Nu, Win[i], u, &x_in[i])});
      //x_in[i] = dot(0, Nu, Win[i], u) ;
    }
    p.await_no_tasks_todo();

    // TODO funzione map
    int diff = 100; //TODO mettere Nr/n_workers 
    int start, stop;
    for (int i = 0; i < Nr; i += diff)
    {
      int start = i;
      int stop = min(i + diff, Nr);
      p.submit({new Comp_state_task(start, stop, Nu, x_rec, x_in, Win, x)});
    }
    p.await_no_tasks_todo();

    for (int i = 0; i < Nr + 1; i++)
    {
      p.submit({new Dot_task(0, Nr+1, P[i], x, &z[i])});
      //z[i] = dot(0, Nr + 1, P[i], x);
    }
    p.await_no_tasks_todo();

    for (int i = 0; i < Ny; i++)
    {
      p.submit({new Dot_task(0, Nr+1, Wout[i], x, &y[i])});
      //y[i] = dot(0, Nr + 1, Wout[i], x);
    }
    p.submit({new Dot_task(0,Nr+1,x,z, &k_den)});
    p.await_no_tasks_todo();

    k_den += l;

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      k[i] = z[i] / k_den;
    }

    // map
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;
      }
    }

    // map
    for (int i = 0; i < Ny; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Wold[i][j] = Wout[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      for (int j = 0; j < Nr + 1; j++)
      {
        Pold[i][j] = P[i][j];
      }
    }

    // map
    for (int i = 0; i < Nr + 1; i++)
    {
      x_old[i] = x[i];
    }

    s = 0;
    for (int i = 0; i < Ny; i++)
    {
      s += pow(d[i] - y[i], 2);
    }
    cout << sqrt(s) << " " << flush;
    errors_norms.push_back(s);

    i++;
  }

  cout << endl;
  errors_norms.erase(errors_norms.begin());
  plt::plot(errors_norms);
  plt::show();

  cout << "\nftt!\n";
  return 0;
}

