/************* STANDARDIZED CALCULATIONS *********************/

#define comp_state_i(x, x_rec, x_in, x_old, Win, Nu, i) \
  {                                                     \
    x[i] = tanh(x_rec[i] + x_in[i] + Win[i][Nu]);       \
    x_old[i] = x[i];                                    \
  }

#define divide_by_const(r, v, k, i) \
  {                                 \
    r[i] = (v[i] / k);              \
  }

#define compute_line_of_wout(start, stop, Wout, Wold, d, y, k, i) \
  {                                                               \
    for (int j = start; j < stop; j++)                            \
    {                                                             \
      Wout[i][j] = Wold[i][j] + (d[i] - y[i]) * k[j];             \
      Wold[i][j] = Wout[i][j];                                    \
    }                                                             \
  }

#define compute_line_of_P(start, stop, P, Pold, k, z, l, i) \
  {                                                         \
    for (int j = start; j < stop; j++)                      \
    {                                                       \
      P[i][j] = (Pold[i][j] - k[i] * z[j]) * 1 / l;         \
      Pold[i][j] = P[i][j];                                 \
    }                                                       \
  }

