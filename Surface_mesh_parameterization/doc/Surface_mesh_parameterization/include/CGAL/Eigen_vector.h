
namespace CGAL {

/*!
\ingroup PkgSurfaceParameterizationAlgebra

The class `Eigen_vector` is a C++ wrapper around \ref thirdpartyEigen "Eigen" vector, which is a simple array of numbers. 
The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be available on the system. 

\models `SparseLinearAlgebraTraits_d::Vector`. 

Parameters 
-------------- 

`T`: Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa http://eigen.tuxfamily.org

*/
template< typename T >
class Eigen_vector {
public:

/// \name Types 
/// @{

/*! 
The internal vector type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_vector */
} /* end namespace CGAL */
