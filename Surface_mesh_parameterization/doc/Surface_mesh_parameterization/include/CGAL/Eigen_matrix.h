
namespace CGAL {

/*!
\ingroup PkgSurfaceParameterizationAlgebra

The class `Eigen_sparse_matrix` is a C++ wrapper around \ref thirdpartyEigen "Eigen" matrix type `Eigen::SparseMatrix` 
that represents general matrices, be they symmetric or not. 
The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be available on the system. 

\cgalModels `SparseLinearAlgebraTraits_d::Matrix` 

Parameters 
-------------- 

`T`: Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`
\sa http://eigen.tuxfamily.org

*/
template< typename T >
class Eigen_sparse_matrix {
public:

/// \name Types 
/// @{

/*! 
The internal matrix type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_sparse_matrix */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceParameterizationAlgebra

The class `Eigen_sparse_symmetric_matrix` is a C++ wrapper around \ref thirdpartyEigen "Eigen" matrix type `Eigen::SparseMatrix`. 

As the matrix is symmetric only the lower triangle part is stored. 

\cgalModels `SparseLinearAlgebraTraits_d::Matrix` 

Parameters 
-------------- 

`T`: Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`
\sa http://eigen.tuxfamily.org

*/
template< typename T >
class Eigen_sparse_symmetric_matrix {
public:

/// \name Types 
/// @{

/*! 
The internal matrix type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_sparse_symmetric_matrix */
} /* end namespace CGAL */
