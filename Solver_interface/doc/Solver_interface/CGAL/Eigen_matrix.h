
namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_sparse_matrix` is a wrapper around \ref thirdpartyEigen "Eigen" matrix type <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html">`Eigen::SparseMatrix` </a>
that represents general matrices, be they symmetric or not. 

\cgalModels `SparseLinearAlgebraTraits_d::Matrix` 

\tparam T Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`
*/
template< typename T >
class Eigen_sparse_matrix {
public:

/// \name Types 
/// @{

/*!
The internal matrix type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef unspecified_type EigenType; 

/// @}

}; /* end Eigen_sparse_matrix */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_sparse_symmetric_matrix` is a wrapper around \ref thirdpartyEigen "Eigen" matrix type <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html">`Eigen::SparseMatrix` </a>

As the matrix is symmetric only the lower triangle part is stored. 

\cgalModels `SparseLinearAlgebraTraits_d::Matrix` 

\tparam T Number type. 

\sa `CGAL::Eigen_solver_traits<T>`
\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_vector<T>`

*/
template< typename T >
class Eigen_sparse_symmetric_matrix {
public:

/// \name Types 
/// @{

/*!
The internal matrix type from \ref thirdpartyEigen "Eigen". 
*/ 
typedef unspecified_type EigenType; 

/// @}

}; /* end Eigen_sparse_symmetric_matrix */


/*!
\ingroup PkgSolver

The class `Eigen_matrix` is a wrapper around \ref thirdpartyEigen "Eigen" 
matrix type 
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html">`Eigen::Matrix`</a>. 


\cgalModels `SvdTraits::Matrix` 

\tparam T Number type. 

\sa `CGAL::Eigen_svd`
\sa `CGAL::Eigen_vector<T>`

*/
template< typename T >
class Eigen_matrix {
public:

};


} /* end namespace CGAL */
