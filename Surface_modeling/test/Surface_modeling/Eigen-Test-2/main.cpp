#include <Eigen/Sparse>
#include <Eigen/SparseLU>

//#define ALSO_TEST_SUPERLU

#ifdef ALSO_TEST_SUPERLU
  #include <Eigen/SuperLUSupport>
#endif

#include <iostream>
#include <fstream>
#include <vector>

template<class SparseMatrix>
void write_matrix_sparse(std::ostream &stream, const SparseMatrix& m)
{
  stream << m.rows() << " " << m.cols() << std::endl;

  for(SparseMatrix::Index r = 0; r < m.rows(); ++r) 
  {
    for(SparseMatrix::Index c = 0; c < m.cols(); ++c) 
    {
      stream << m.coeff(r,c) << " ";
    }
    stream << std::endl;
  }
}

template<class SparseMatrix>
void read_matrix_sparse(std::istream &stream, SparseMatrix& m)
{
  SparseMatrix::Index rows, cols;
  stream >> rows >> cols;

  std::vector<Eigen::Triplet<double> > triplets;
  for(SparseMatrix::Index r = 0; r < rows; ++r)
  for(SparseMatrix::Index c = 0; c < cols; ++c)
  {
    double val;
    stream >> val;
    if(val != 0.0) {
      triplets.push_back(Eigen::Triplet<double>(r, c, val));
    }
  }

  m.resize(rows, cols);
  m.setFromTriplets(triplets.begin(), triplets.end());
  m.makeCompressed();
}

template<class Matrix>
void write_matrix(std::ostream &stream, const Matrix& m)
{
  stream << m.rows() << " " << m.cols() << std::endl;

  for(Matrix::Index r = 0; r < m.rows(); ++r) 
  {
    for(Matrix::Index c = 0; c < m.cols(); ++c) 
    {
      stream << m(r,c) << " ";
    }
    stream << std::endl;
  }
}

template<class Matrix>
void read_matrix(std::istream &stream, Matrix& m)
{
  Matrix::Index rows, cols;
  stream >> rows >> cols;
  m = Matrix(rows, cols);

  for(Matrix::Index r = 0; r < rows; ++r)
  for(Matrix::Index c = 0; c < cols; ++c)
  {
    stream >> m(r,c);    
  }
}

template<class Matrix>
void print_dif(const Matrix& m1, const Matrix& m2)
{
  std::cout << "-----------print dif enter...-----------" << std::endl;
  const double threshold = 1e-5;

  for(Matrix::Index r = 0; r < m1.rows(); ++r)
  for(Matrix::Index c = 0; c < m1.cols(); ++c)
  {
    double abs_dif = std::abs(m1(r, c) - m2(r, c));
    if(abs_dif > threshold)
    {
      std::cout << "Dif greater than: " << threshold << ", dif: " << abs_dif << std::endl; 
    }  
  }
  std::cout << "-----------print dif end!-----------" << std::endl;
}

int main(int argc, char* argv[]) {
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSparseMatrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1>     EigenVector;
  typedef Eigen::SparseLU<EigenSparseMatrix, Eigen::COLAMDOrdering<int> > SparseLUSolver;

  std::ifstream input_M("A.txt");
  EigenSparseMatrix M; // The input matrix M should be in a compressed and column-major form... (from documenetation)
  read_matrix_sparse(input_M, M); // it also compress M
  
  std::ifstream input_vM("Bx.txt");
  EigenVector vM;
  read_matrix(input_vM, vM);

  SparseLUSolver sparseLU;
  sparseLU.compute(M);
  std::cout << ((sparseLU.info() == Eigen::Success) 
                ? "SUCCESS" : "FAIL") << std::endl;
  std::cout << sparseLU.lastErrorMessage() << std::endl;

  EigenVector sparseLU_result = sparseLU.solve(vM);
  EigenVector sparseLU_expected_vM = M * sparseLU_result;

  std::cout << "A*x = b is solved using sparseLu, dif between A*x and b is:" << std::endl;
  print_dif(vM, sparseLU_expected_vM);

#ifdef ALSO_TEST_SUPERLU
  typedef Eigen::SuperLU<EigenSparseMatrix> SuperLUSolver;

  SuperLUSolver superLU;
  superLU.compute(M);

  std::cout << ((superLU.info() == Eigen::Success) 
           ? "SUCCESS" : "FAIL") << std::endl;

  EigenVector superLU_result = superLU.solve(vM);
  EigenVector superLU_expected_vM = M * superLU_result;

  std::cout << "A*x = b is solved using superLU, dif between A*x and b is:" << std::endl;
  print_dif(vM, superLU_expected_vM);

  //print_dif(superLU_result, sparseLU_result);
#endif
  int dum;
  std::cin >> dum;
}