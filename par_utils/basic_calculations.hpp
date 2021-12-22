/************* STANDARDIZED CALCULATIONS *********************/
// Together with the dot product, these are the low level math operations that constitute the
// heart of the calculations done for this project

// computes i-th entry in the new state vector, copies old entry into apposit vector
#define comp_state_i(x, x_rec, x_in, x_old, Win, Nu, i) \
  {                                                     \
    x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]);       \
    x_old[i] = x[i];                                    \
  }

// divides an entry of a vector by a constant
#define divide_by_const(r, v, k, i) \
  {                                 \
    r[i] = (v[i] / k);              \
  }

// computes a row of the new Wout matrix, copies its values into apposite matrix
#define compute_line_of_wout(start, stop, Wout, Wold, d, y, k, i) \
  {                                                               \
    for (int j = start; j < stop; j++)                            \
    {                                                             \
      Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];             \
      Wold[i][j] = Wout[i][j];                                    \
    }                                                             \
  }

// computes a row of the new P matrix, copies its values into apposite matrix
#define compute_line_of_P(start, stop, P, Pold, k, z, l, i) \
  {                                                         \
    for (int j = start; j < stop; j++)                      \
    {                                                       \
      P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;         \
      Pold[i][j] = P[i][j];                                 \
    }                                                       \
  }
