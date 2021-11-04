#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

class Matrix_wrapper
{
public:
  Matrix_wrapper(float **m, int n1, int n2);
  Matrix_wrapper() = default;
  float **m;
  int n1;
  int n2;

  Matrix_wrapper operator+(const Matrix_wrapper &m2);
  Matrix_wrapper operator|(const Matrix_wrapper &m2);
  Matrix_wrapper operator*(const Matrix_wrapper &m2);
  Matrix_wrapper operator-(const Matrix_wrapper &m2);
  Matrix_wrapper operator*(const float &k);
  Matrix_wrapper operator+(const float &k);
  Matrix_wrapper operator/(const float &k);
  Matrix_wrapper operator-(const float &k);
  Matrix_wrapper get_line(int i);
  Matrix_wrapper transpose();
  float to_float();
};

void print_matrix(Matrix_wrapper m);
Matrix_wrapper zeros(int n1, int n2);
Matrix_wrapper ones(int n1, int n2);
Matrix_wrapper eye(int n);
Matrix_wrapper generate_random_sparse_matrix(int n1, int n2, float d, float interval_start = -1, float interval_stop = 1);
Matrix_wrapper build_sparse_contractive_matrix(int n1, int n2);
Matrix_wrapper vstack(Matrix_wrapper m1, Matrix_wrapper m2);
Matrix_wrapper elementwise_tanh(Matrix_wrapper m1);
Matrix_wrapper copy(Matrix_wrapper m1);
Matrix_wrapper from_array(float *a, int n);
Matrix_wrapper normalize(Matrix_wrapper mat);
double norm(Matrix_wrapper mat);
double mean(Matrix_wrapper mat);
double var(Matrix_wrapper mat);
void free_matrices(std::vector<Matrix_wrapper> matrices);
float spectral_radius(MatrixXd M);