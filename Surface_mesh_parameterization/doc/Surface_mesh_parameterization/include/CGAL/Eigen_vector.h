
namespace CGAL {

/*!
\ingroup PkgSurfaceParameterization

The class `Eigen_vector` is a C++ wrapper around <span class="textsc">Eigen</span>'s vector, which is a simple array of numbers. 
The version 3.1 (or greater) of <span class="textsc">Eigen</span> must be available on the system. 

CONVERRORIsModel: `SparseLinearAlgebraTraits_d::Vector`. 

Parameters 
-------------- 

`T`: Number type. 

CONVERRORSeeAlso: \ref ::CGAL::Eigen_solver_traits<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_sparse_matrix<T> 
CONVERRORSeeAlso: \ref ::CGAL::Eigen_sparse_symmetric_matrix<T> 

*/
template< typename T >
class Eigen_vector {
public:

/// \name Types 
/// @{

/*! 
The internal vector type from <span class="textsc">Eigen</span>. 
*/ 
typedef Hidden_type EigenType; 

/// @}

}; /* end Eigen_vector */
} /* end namespace CGAL */
