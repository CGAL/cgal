
namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_solver_traits` provides an interface to the sparse solvers of \ref thirdpartyEigen "Eigen".
The version 3.1 (or greater) of \ref thirdpartyEigen "Eigen" must be available on the system. 

\cgalModels `SparseLinearAlgebraWithFactorTraits_d` and `NormalEquationSparseLinearAlgebraTraits_d`


\tparam T A sparse solver of \ref thirdpartyEigen "Eigen". The default solver is the iterative bi-congugate gradient stabilized solver  `Eigen::BiCGSTAB` for `double`. 

\sa `CGAL::Eigen_sparse_matrix<T>`
\sa `CGAL::Eigen_sparse_symmetric_matrix<T>`
\sa `CGAL::Eigen_vector<T>`
\sa http://eigen.tuxfamily.org

Example 
-------------- 

The instantiation of this class assumes an \ref thirdpartyEigen "Eigen" sparse solver is provided. Here are few examples: 

\code{.cpp} 

typedef CGAL::Eigen_sparse_matrix<double>::EigenType EigenMatrix; 

//iterative general solver 
typedef CGAL::Eigen_solver_traits< Eigen::BiCGSTAB<EigenMatrix> > Iterative_general_solver; 

//iterative symmetric solver 
typedef CGAL::Eigen_solver_traits< Eigen::ConjugateGradient<EigenMatrix> > Iterative_symmetric_solver; 

//direct symmetric solver 
typedef CGAL::Eigen_solver_traits< Eigen::SimplicialCholesky<EigenMatrix> > Direct_symmetric_solver; 

\endcode 

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
If `T` is `Eigen::ConjugateGradient<M>` or `Eigen::SimplicialCholesky<M>`, `Matrix` is `CGAL::Eigen_sparse_symmetric_matrix<T>` 
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
