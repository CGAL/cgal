
//license

#ifndef CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
#define CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H

#include <CGAL/basic.h>
#include <Eigen/Dense>
#include <Eigen/QR>


namespace CGAL {


template<typename T>
class Eigen_dense_matrix
{

public:

  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> EigenType;
  typedef Eigen::Vector3d VectorType;

  Eigen_dense_matrix(std::size_t nrows, std::size_t ncols)
    : m_matrix(static_cast<int>(nrows), static_cast<int>(ncols))
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(int nrows, int ncols)
    : m_matrix(nrows, ncols)
  {
    CGAL_assertion(m_matrix.rows() > 0);
    CGAL_assertion(m_matrix.cols() > 0);
  }

  Eigen_dense_matrix(const EigenType& eigen_mat)
    : m_matrix(eigen_mat) {}

  Eigen_dense_matrix() : m_matrix() {}

  ~Eigen_dense_matrix() {}

  T& operator() (int i_, int j_)
  {
    return m_matrix(i_, j_);
  }

  void set_coef(std::size_t i_, std::size_t j_, T val)
  {
    int i = static_cast<int>(i_);
    int j = static_cast<int>(j_);
    CGAL_precondition(i < m_matrix.rows());
    CGAL_precondition(j < m_matrix.cols());

    m_matrix.coeffRef(i,j) = val;
  }

  std::size_t rows() const {return m_matrix.rows();}
  std::size_t cols() const {return m_matrix.cols();}

  void resize(int i_, int j_) { m_matrix.resize(i_, j_);}

  const Eigen_dense_matrix row(int i_) const
  {
    VectorType row = m_matrix.row(i_);
    Eigen_dense_matrix r(row);
    return r;
  }

  const T& coeff(int i_) const
  {
    return m_matrix.coeff(i_);
  }

  const Eigen_dense_matrix transpose() const
  {
    Eigen_dense_matrix dm(m_matrix.transpose());
    return dm;
    //return Eigen_dense_matrix(m_matrix.transpose());
  }

  const double determinant() const
  {
    return m_matrix.determinant();
  }

  void qr_factorization(Eigen_dense_matrix& Q)
  {
    Q.resize(m_matrix.rows(), m_matrix.cols());
    Eigen::HouseholderQR<EigenType> qr(m_matrix);
    Eigen_dense_matrix qr_mat(qr.householderQ());
    Q = qr_mat;
  }

  friend void qr_factorization(std::vector<Eigen_dense_matrix>& simplex)
  {
    for(int i = 0; i < simplex.size(); ++i)
    {
      Eigen::HouseholderQR<EigenType> qr(simplex[i].m_matrix);
      Eigen_dense_matrix Q(qr.householderQ());
      simplex[i] = Q;
    }
  }

  friend const Eigen_dense_matrix operator* (const Eigen_dense_matrix& A,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix product(A.m_matrix * B.m_matrix);
    return product;
  }

  friend const Eigen_dense_matrix operator* (const T& scalar,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix product(scalar * B.m_matrix);
    return product;
  }

  friend const Eigen_dense_matrix operator+ (const Eigen_dense_matrix& A,
                                             const Eigen_dense_matrix& B)
  {
    const Eigen_dense_matrix sum_result(A.m_matrix + B.m_matrix);
    return sum_result;
  }

  friend const Eigen_dense_matrix operator/ (const Eigen_dense_matrix& A,
                                             const double& scalar)
  {
    const Eigen_dense_matrix product(A.m_matrix / scalar);
    return product;
  }

  friend const Eigen_dense_matrix operator/ (const double& scalar ,
                                             const Eigen_dense_matrix& A)
  {
    const Eigen_dense_matrix product(scalar / A.m_matrix);
    return product;
  }



private:
  mutable EigenType m_matrix;
};







} // end namespace



#endif // CGAL_EIGEN_LINEAR_ALGEBRA_TRAITS_H
