
namespace CGAL {

/*!
\ingroup PkgSurfaceParameterization

The class `Eigen_sparse_matrix` is a C++ wrapper around Eigen's matrix type `Eigen::SparseMatrix` 
that represents general matrices, be they symmetric or not. 
The version 3.1 (or greater) of <span class="textsc">Eigen</span> must be available on the system. 

CONVERRORIsModel: Model of the `SparseLinearAlgebraTraits_d::Matrix` concept. 

Parameters 
-------------- 

`T`: Number type. 

CONVERRORSeeAlso: \ref ::CGAL::Eigen_solver_traits<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_sparse_symmetric_matrix<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_vector<T> 

*/
template< typename T >
class Eigen_sparse_matrix {
public:

/// \name Types 
/// @{

/*! 
The internal matrix type from <span class="textsc">Eigen</span>. 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_sparse_matrix */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSurfaceParameterization

The class `Eigen_sparse_symmetric_matrix` is a C++ wrapper around <span class="textsc">Eigen</span>'s matrix type `Eigen::SparseMatrix`. 

As the matrix is symmetric only the lower triangle part is stored. 

CONVERRORIsModel: Model of the `SparseLinearAlgebraTraits_d::Matrix` concept. 

Parameters 
-------------- 

`T`: Number type. 

CONVERRORSeeAlso: \ref ::CGAL::Eigen_solver_traits<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_sparse_symmetric_matrix<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_vector<T> 

*/
template< typename T >
class Eigen_sparse_symmetric_matrix {
public:

/// \name Types 
/// @{

/*! 
The internal matrix type from <span class="textsc">Eigen</span>. 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_sparse_symmetric_matrix */
} /* end namespace CGAL */
