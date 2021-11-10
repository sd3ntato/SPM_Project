#ifndef linear_algebra
#define linear_algebra
#include "linear_algebra.h"
#endif

using namespace std;


Matrix_wrapper::Matrix_wrapper(float **m, int n1, int n2)
{
  this->m = m;
  this->n1 = n1;
  this->n2 = n2;
}

Matrix_wrapper Matrix_wrapper::operator|(const Matrix_wrapper &m2)
{
  if (this->n2 != m2.n1)
  {
    cout << "dot product: invalid dimensions : left matrix " << this->n1 << " " << this->n2 << endl;
    cout << "dot product: invalid dimensions : right matrix " << m2.n1 << " " << m2.n2 << endl;
    assert(false);
  }

  Matrix_wrapper res = zeros(this->n1, m2.n2);
  for (int i = 0; i < this->n1; i++)
  { // scorre le righe della prima matrice
    for (int j = 0; j < m2.n2; j++)
    { // scorren le colonne della seconda matrice
      for (int cnt = 0; cnt < this->n2; cnt++)
      { //contatore per scorrere riga-colonna
        res.m[i][j] += (this->m[i][cnt] * m2.m[cnt][j]);
      }
    }
  }

  return res;
}

Matrix_wrapper Matrix_wrapper::operator+(const Matrix_wrapper &m2)
{
  if (this->n1 != m2.n1 or this->n2 != m2.n2)
  {
    cout << "elementwise sum: invalid dimensions : left matrix " << this->n1 << " " << this->n2 << endl;
    cout << "elementwise sum: invalid dimensions : right matrix " << m2.n1 << " " << m2.n2 << endl;
    assert(false);
  }

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] + m2.m[i][j];
    }
  }

  return res;
}

Matrix_wrapper Matrix_wrapper::operator*(const Matrix_wrapper &m2)
{
  if (this->n1 != m2.n1 or this->n2 != m2.n2)
  {
    cout << "elementwise multiplication: invalid dimensions" << endl;
    assert(false);
  }

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] * m2.m[i][j];
    }
  }

  return res;
}

Matrix_wrapper Matrix_wrapper::operator-(const Matrix_wrapper &m2)
{
  if (this->n1 != m2.n1 or this->n2 != m2.n2)
  {
    cout << "elementwise subtraction: invalid dimensions" << endl;
    assert(false);
  }

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] - m2.m[i][j];
    }
  }

  return res;
}

Matrix_wrapper Matrix_wrapper::operator*(const float &k)
{

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] * k;
    }
  }
  return res;
}

Matrix_wrapper Matrix_wrapper::operator+(const float &k)
{

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] + k;
    }
  }
  return res;
}

Matrix_wrapper Matrix_wrapper::operator/(const float &k)
{

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] / k;
    }
  }
  return res;
}

Matrix_wrapper Matrix_wrapper::operator-(const float &k)
{

  Matrix_wrapper res = zeros(this->n1, this->n2);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[i][j] = this->m[i][j] - k;
    }
  }
  return res;
}

Matrix_wrapper Matrix_wrapper::get_line(int i)
{
  if (i > this->n1)
  {
    cout << "Martix_wrapper::get_line: invalid argument grater than number of lines available" << endl;
    assert(false);
  }
  float **lines = new float *[1];
  lines[0] = this->m[i];
  return Matrix_wrapper(lines, 1, this->n2);
}

Matrix_wrapper Matrix_wrapper::transpose()
{
  Matrix_wrapper res = zeros(this->n2, this->n1);
  for (int i = 0; i < this->n1; i++)
  {
    for (int j = 0; j < this->n2; j++)
    {
      res.m[j][i] = this->m[i][j];
    }
  }
  return res;
}

float Matrix_wrapper::to_float()
{
  if (not(this->n1 == 1 and this->n2 == 1))
  {
    cout << "Matrix_wrapper::operator float: invalid size, not (1,1)" << endl;
    assert(false);
  }
  return this->m[0][0];
}

void print_matrix(Matrix_wrapper mat)
{
  float **m = mat.m;
  int n1 = mat.n1;
  int n2 = mat.n2;
  for (int i = 0; i < n1; i++)
  {
    cout << "[ ";
    for (int ii = 0; ii < n2; ii++)
    {
      cout << m[i][ii] << " , ";
    }
    cout << "]\n";
  }
  cout << "\n";
}

Matrix_wrapper zeros(int n1, int n2)
{
  float **m = new float *[n1];
  for (int i = 0; i < n1; i++)
  {
    m[i] = new float[n2];
    for (int j = 0; j < n2; j++)
    {
      m[i][j] = 0;
    }
  }
  Matrix_wrapper mat = Matrix_wrapper(m, n1, n2);
  return mat;
}

Matrix_wrapper ones(int n1, int n2)
{
  float **m = new float *[n1];
  for (int i = 0; i < n1; i++)
  {
    m[i] = new float[n2];
    for (int j = 0; j < n2; j++)
    {
      m[i][j] = 1;
    }
  }
  Matrix_wrapper mat = Matrix_wrapper(m, n1, n2);
  return mat;
}

bool contains(int *l, int n, int t)
{
  for (int i = 0; i < n; i++)
  {
    if (l[i] == t)
    {
      return true;
    }
  }
  return false;
}

Matrix_wrapper generate_random_sparse_matrix(int n1, int n2, float d, float interval_start, float interval_stop)
{

  Matrix_wrapper mat = zeros(n1, n2);
  float **m = mat.m;

  // compute #elements to be filled out
  int n = ceil(n1 * n2 * d);

  // generate a list of n random values using uniform distribution on floats
  float *numbers = new float[n];
  int *positions = new int[n];

  // put all positions to -1 so that whe rand number generator gets a 0 he put it in only one time
  for (int i = 0; i < n; i++)
  {
    positions[i] = -1;
  }

  random_device rd; // Will be used to obtain a seed for the random number engine
  mt19937 gen(rd());
  uniform_real_distribution<> fdis(interval_start, interval_stop);
  uniform_int_distribution<> idis(0, n1 * n2 - 1);
  for (int i = 0; i < n; i++)
  {
    numbers[i] = fdis(gen); // decido che numero voglio inserire

    /* per la posizione del numero: non posso scegliere la stessa posizione due volte,
        // che significa che la lista positions non puo' avere duplicati.
        // modo molto stupido per inserire senza fare duplicati: prima di inserire controllo se 
        // il numero che sto inserendo l'ho gia messo. Se l'ho gia messo non va bene quindi ne prendo un altro
        */
    int p = idis(gen); // prendo un numero a caso
    while (contains(positions, n1 * n2, p))
    {                //controllo se gia sta nella lista
      p = idis(gen); // se ho verificato che quello che volevo inserire gia sta nella lista ne prendo un altro
    }                // alla fine sono sicuro che not constains(positions,p) quindi posso inserire tranquillamente p nella lista
    positions[i] = p;
  }

  // inserisco i numeri che ho scelto, nelle relative posizioni all'interno della matrice
  for (int i = 0; i < n; i++)
  {
    int p = positions[i];
    int row = (int)floor(p / n2);
    int col = p % n2;
    m[row][col] = numbers[i];
  }
  return mat;
}

Matrix_wrapper vstack(Matrix_wrapper m1, Matrix_wrapper m2)
{
  if (m1.m)
  {
    if (m1.n2 != m2.n2)
    {
      cout << "vstack: invalid dimensions!" << endl;
      assert(false);
    }
    Matrix_wrapper res = zeros(m1.n1 + m2.n1, m1.n2);
    for (int i = 0; i < m1.n1; i++)
    {
      for (int j = 0; j < m1.n2; j++)
      {
        res.m[i][j] = m1.m[i][j];
      }
    }
    for (int i = 0; i < m2.n1; i++)
    {
      for (int j = 0; j < m1.n2; j++)
      {
        res.m[i + m1.n1][j] = m2.m[i][j];
      }
    }
    return res;
  }
  else
  {
    return m2;
  }
  return Matrix_wrapper(nullptr, 0, 0);
}

Matrix_wrapper elementwise_tanh(Matrix_wrapper m1)
{
  Matrix_wrapper res = zeros(m1.n1, m1.n2);
  for (int i = 0; i < m1.n1; i++)
  {
    for (int j = 0; j < m1.n2; j++)
    {
      res.m[i][j] = tanh(m1.m[i][j]);
    }
  }
  return res;
}

Matrix_wrapper copy(Matrix_wrapper m1)
{
  Matrix_wrapper res = zeros(m1.n1, m1.n2);
  for (int i = 0; i < m1.n1; i++)
  {
    for (int j = 0; j < m1.n2; j++)
    {
      res.m[i][j] = m1.m[i][j];
    }
  }
  return res;
}

Matrix_wrapper from_array(float *a, int n)
{
  float **m = new float *[1];
  m[0] = a;
  Matrix_wrapper mat = Matrix_wrapper(m, 1, n);
  return mat;
}

Matrix_wrapper eye(int n)
{
  Matrix_wrapper res = zeros(n, n);
  for (int i = 0; i < n; i++)
  {
    res.m[i][i] = 1.0;
  }
  return res;
}

double norm(Matrix_wrapper mat)
{
  if (not(mat.n2 == 1))
  {
    cout << "norm: this function is for column vectors only! mat.n2 = " << mat.n2 << endl;
  }
  float n = 0.0;
  for (int i = 0; i < mat.n1; i++)
  {
    n = n + pow(mat.m[i][0], 2);
  }
  return sqrt(n);
}

double mean(Matrix_wrapper mat)
{
  float s = 0.0;
  for (int i = 0; i < mat.n1; i++)
  {
    for (int j = 0; j < mat.n2; j++)
    {
      s = s + mat.m[i][j];
    }
  }
  s = s / (mat.n1 * mat.n2);
  return s;
}

Matrix_wrapper square(Matrix_wrapper mat)
{
  Matrix_wrapper res = zeros(mat.n1, mat.n2);
  for (int i = 0; i < mat.n1; i++)
  {
    for (int j = 0; j < mat.n2; j++)
    {
      res.m[i][j] = pow(mat.m[i][j], 2);
    }
  }
  return res;
}

double var(Matrix_wrapper mat)
{
  return mean(square(mat - mean(mat)));
}

Matrix_wrapper normalize(Matrix_wrapper mat)
{
  Matrix_wrapper res = copy(mat);
  auto m = mean(res);
  auto dev = sqrt(var(res));

  res = (res - m) / dev;

  assert(mean(res) < 0.001);
  assert(abs(var(res) - 1) < 0.001);
  return res;
}

void free_matrices(vector<Matrix_wrapper> matrices)
{
  for (int i = 0; i < (int)matrices.size(); i++)
  {
    float **mat = matrices[i].m;
    for (int j = 0; j < matrices[i].n1; j++)
    {
      delete[] mat[j];
    }
    delete[] mat;
  }
}

using namespace Spectra;
using namespace Eigen;

float spectral_radius(MatrixXd M)
{
  // Construct matrix operation object using the wrapper class
  DenseGenMatProd<double> op(M);

  // Construct eigen solver object, requesting the largest
  // (in magnitude, or norm) three eigenvalues
  GenEigsSolver<DenseGenMatProd<double>> eigs(op, 1, 6);

  // Initialize and compute
  eigs.init();
  eigs.compute(SortRule::LargestMagn);

  // Retrieve results
  Eigen::VectorXcd evalues;
  if (eigs.info() == CompInfo::Successful)
  {
    evalues = eigs.eigenvalues();
    //std::cout << "Eigenvalues found:" << evalues(0) << std::endl;
    return abs(evalues(0));
  }
  cout << "eig comp falied";
  return 0;
}

Matrix_wrapper build_sparse_contractive_matrix(int n1, int n2)
{
  MatrixXd mat = (MatrixXd::Random(n1, n2).array() > 0.9).cast<double>() * MatrixXd::Random(n1, n2).array();
  float rho_act = spectral_radius(mat);

  if (rho_act < 0.1)
  {
    cout << "recurrent matrix creation failed bc spectral radius too small" << endl;
    assert(false);
  }

  mat = mat * (0.8 / rho_act);

  cout << "spectral radius of the  recurrent matrix: " << spectral_radius(mat) << endl;

  Matrix_wrapper mm = zeros(n1, n2);
  float **m = mm.m;
  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      m[i][j] = mat(i, j);
    }
  }
  return mm;
}



float dot(int start, int stop, float *v1, float *v2)
{
  float s = 0.0;
  for (int i = start; i < stop; i++)
  {
    s += v1[i] * v2[i];
  }
  return s;
}
