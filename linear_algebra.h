
class Matrix_wrapper{
  public:
   Matrix_wrapper(float** m, int n1, int n2);
   Matrix_wrapper() = default;
   float **m;
   int n1;
   int n2;
};

void print_matrix(Matrix_wrapper m);
Matrix_wrapper zeros(int n1, int n2);
Matrix_wrapper generate_random_sparse_matrix(int n1, int n2, float d, float interval_start=-1, float interval_stop=1);
Matrix_wrapper build_sparse_contractive_matrix(int n1, int n2);