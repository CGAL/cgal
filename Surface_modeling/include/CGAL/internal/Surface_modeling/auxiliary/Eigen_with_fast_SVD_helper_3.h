#ifndef CGAL_EIGEN_SVD_HELPER_3_H
#define CGAL_EIGEN_SVD_HELPER_3_H

#include <Eigen/Eigen>

namespace CGAL {

class Eigen_SVD_helper_3{
public:

  typedef Eigen::Matrix3f Matrix;
  typedef Eigen::Vector3f Vector;
  typedef std::pair<Matrix, Matrix> Solver;
  
  /// Equivalent to `result = m1*m2^t`
  void matrix_matrix_transpose_mult(Matrix& result, const Matrix& m1, const Matrix& m2) 
  { 
    result = m1*m2.transpose();
  }

  /// Equivalent to `result += w * (v1*v2^t)`
  void scalar_vector_vector_transpose_mult(Matrix& result, double w, const Vector& v1, const Vector& v2) 
  {
    result += w * (v1*v2.transpose());
  }
  
  /// Equivalent to `result += (w1*m1 + w2*m2) * v`
  void scalar_matrix_scalar_matrix_vector_mult(Vector& result, double w1, const Matrix& m1, double w2, const Matrix& m2, const Vector& v) 
  {
    result += (w1*m1 + w2*m2) * v;
  }
  
  /// Equivalent to `result += w * (m1 + m2 + m3) * v`
  void scalar_mult_with_matrix_sum(Vector& result, double w, const Matrix& m1, const Matrix& m2, const Matrix& m3, const Vector& v) 
  { 
    result += w * (m1 + m2 + m3) * v;
  }
  
  /// Returns the squared norm of `v1 - m*v2`
  double squared_norm_vector_scalar_vector_subs(const Vector& v1, const Matrix& m, const Vector& v2) 
  { 
    return (v1 - m*v2).squaredNorm(); 
  } 

  /// Negates column `i` of matrix `result`
  void negate_column(Matrix& result, int i) 
  {
    result.col(i) *= -1;  
  }

  /// Returns an identity matrix
  Matrix identity_matrix() 
  {
    return Matrix().setIdentity();
  }
  
  /// Returns a zero initialized matrix
  Matrix zero_matrix() 
  {
    return Matrix().setZero();
  }

  /// Returns vector initialized with parameters
  Vector vector(double x, double y, double z)
  {
    return Vector(x, y, z);
  }
  
  /// Returns a coefficient of a vector
  double vector_coeff(const Vector& v, int i)
  {
    return v(i);
  }
  
  /// Returns determinant of m
  double determinant(const Matrix& m) 
  {
    return m.determinant();
  }

  // #define USE_SCALAR_IMPLEMENTATION
  #define USE_SSE_IMPLEMENTATION
  #define COMPUTE_V_AS_MATRIX
  #define COMPUTE_U_AS_MATRIX

  #include "Singular_Value_Decomposition_Preamble.hpp"

  /// Computes the singular value decomposition
  Solver compute_svd(const Matrix& m)
  {
   // return std::make_pair(identity_matrix(), identity_matrix());

	#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"

    ENABLE_SCALAR_IMPLEMENTATION(Sa11.f=m(0,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa21.f=m(1,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa31.f=m(2,0);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa12.f=m(0,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa22.f=m(1,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa32.f=m(2,1);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa13.f=m(0,2);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa23.f=m(1,2);) 
    ENABLE_SCALAR_IMPLEMENTATION(Sa33.f=m(2,2);) 

    ENABLE_SSE_IMPLEMENTATION(Va11=_mm_set1_ps(m(0,0));)
    ENABLE_SSE_IMPLEMENTATION(Va21=_mm_set1_ps(m(1,0));)
    ENABLE_SSE_IMPLEMENTATION(Va31=_mm_set1_ps(m(2,0));)
    ENABLE_SSE_IMPLEMENTATION(Va12=_mm_set1_ps(m(0,1));)
    ENABLE_SSE_IMPLEMENTATION(Va22=_mm_set1_ps(m(1,1));)
    ENABLE_SSE_IMPLEMENTATION(Va32=_mm_set1_ps(m(2,1));)
    ENABLE_SSE_IMPLEMENTATION(Va13=_mm_set1_ps(m(0,2));)
    ENABLE_SSE_IMPLEMENTATION(Va23=_mm_set1_ps(m(1,2));)
    ENABLE_SSE_IMPLEMENTATION(Va33=_mm_set1_ps(m(2,2));)

    #include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"  
    
    Solver solver;

#ifdef USE_SCALAR_IMPLEMENTATION
    solver.first(0,0) = Su11.f;
    solver.first(1,0) = Su21.f;
    solver.first(2,0) = Su31.f;
    solver.first(0,1) = Su12.f;
    solver.first(1,1) = Su22.f;
    solver.first(2,1) = Su32.f;
    solver.first(0,2) = Su13.f;
    solver.first(1,2) = Su23.f;
    solver.first(2,2) = Su33.f;

    solver.second(0,0) = Sv11.f;
    solver.second(1,0) = Sv21.f;
    solver.second(2,0) = Sv31.f;
    solver.second(0,1) = Sv12.f;
    solver.second(1,1) = Sv22.f;
    solver.second(2,1) = Sv32.f;
    solver.second(0,2) = Sv13.f;
    solver.second(1,2) = Sv23.f;
    solver.second(2,2) = Sv33.f;
#endif

#ifdef USE_SSE_IMPLEMENTATION
    float buf[4];
    _mm_storeu_ps(buf,Vu11);solver.first(0,0)=buf[0];
    _mm_storeu_ps(buf,Vu21);solver.first(1,0)=buf[0];
    _mm_storeu_ps(buf,Vu31);solver.first(2,0)=buf[0];
    _mm_storeu_ps(buf,Vu12);solver.first(0,1)=buf[0];
    _mm_storeu_ps(buf,Vu22);solver.first(1,1)=buf[0];
    _mm_storeu_ps(buf,Vu32);solver.first(2,1)=buf[0];
    _mm_storeu_ps(buf,Vu13);solver.first(0,2)=buf[0];
    _mm_storeu_ps(buf,Vu23);solver.first(1,2)=buf[0];
    _mm_storeu_ps(buf,Vu33);solver.first(2,2)=buf[0];

    _mm_storeu_ps(buf,Vv11);solver.second(0,0)=buf[0];
    _mm_storeu_ps(buf,Vv21);solver.second(1,0)=buf[0];
    _mm_storeu_ps(buf,Vv31);solver.second(2,0)=buf[0];
    _mm_storeu_ps(buf,Vv12);solver.second(0,1)=buf[0];
    _mm_storeu_ps(buf,Vv22);solver.second(1,1)=buf[0];
    _mm_storeu_ps(buf,Vv32);solver.second(2,1)=buf[0];
    _mm_storeu_ps(buf,Vv13);solver.second(0,2)=buf[0];
    _mm_storeu_ps(buf,Vv23);solver.second(1,2)=buf[0];
    _mm_storeu_ps(buf,Vv33);solver.second(2,2)=buf[0];
#endif

    return solver;
  }

  int get_smallest_singular_value_index(const Solver& /*solver*/)
  {
    return 2; // singular values are always sorted in decreasing order so use column 2
  }

  /// Gets matrix U from solver
  const Matrix& get_matrixU(const Solver& solver) 
  {
    return solver.first;
  }

  /// Gets matrix V from solver
  const Matrix& get_matrixV(const Solver& solver) 
  {
    return solver.second;
  }
};

}//namespace CGAL
#endif // CGAL_EIGEN_SVD_HELPER_3_H