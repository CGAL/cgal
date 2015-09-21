
namespace CGAL {

/*!
\ingroup PkgSolver

The class `Eigen_svd` provides an algorithm to solve in the least 
square sense a linear system with a singular value decomposition using 
\ref thirdpartyEigen. 

\cgalModels `SvdTraits`

*/

class Eigen_svd {
public:

  typedef double FT;

  typedef Eigen_vector<FT> Vector;

  typedef Eigen_matrix<FT> Matrix;

}; /* end Eigen_svd */
} /* end namespace CGAL */
