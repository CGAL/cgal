
namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_solver_traits` provides an interface to the sparse solvers of \ref thirdpartyEigen "Eigen".
The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be available on the system. 

\cgalModels `SparseLinearAlgebraTraitsWithFactor_d`


\tparam T a sparse solver of \ref thirdpartyEigen "Eigen". The default solver is the iterative bi-congugate gradient stabilized solver 
<a href="http://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html">`Eigen::BiCGSTAB`</a> for `double`. 

\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`


*/
template< typename T >
class Eigen_solver_traits {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef typename T::Scalar NT; 

/*!

*/ 
typedef CGAL::Eigen_vector<NT> Vector; 

/*!
If `T` is <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html">`Eigen::ConjugateGradient<M>`</a> or <a href="http://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialCholesky.html">`Eigen::SimplicialCholesky<M>`</a>, `Matrix` is `CGAL::Eigen_sparse_symmetric_matrix<T>` 
and `CGAL::Eigen_sparse_matrix<T>` otherwise. 

*/ 
typedef unspecified_type Matrix; 

/// @} 

/// \name Operations 
/// @{

/*!

Returns a reference to the internal \ref thirdpartyEigen "Eigen" solver. This function can be used for example to set specific parameters of the solver. 

*/ 
T& solver(); 

/// @}

}; /* end Eigen_solver_traits */
} /* end namespace CGAL */
