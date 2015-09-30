
namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_vector` is a wrapper around \ref thirdpartyEigen "Eigen" vector
type <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html"> </a>, 
which is a simple array of numbers. 

\cgalModels `SvdTraits::Vector` 
\cgalModels `SparseLinearAlgebraTraits_d::Vector`. 


\tparam T  Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`

*/
template< typename T >
class Eigen_vector {
public:

/// \name Types 
/// @{

/*!
The internal vector type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef unspecified_type EigenType; 

/// @}

}; /* end Eigen_vector */
} /* end namespace CGAL */
