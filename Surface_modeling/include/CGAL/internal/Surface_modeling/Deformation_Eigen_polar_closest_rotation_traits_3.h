#ifndef CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H
#define CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H
/// @cond CGAL_DOCUMENT_INTERNAL

#include <Eigen/Eigen>
#include <Eigen/SVD>

#include <CGAL/Deformation_Eigen_closest_rotation_traits_3.h>
#include <CGAL/Profile_counter.h>
namespace CGAL {
  /// \ingroup PkgSurfaceModeling
  /// A wrapper class to compute closest rotation to a 3x3 Matrix using `Eigen` library.
  /// The internal computation relies on hybrid system of `Eigen::SelfAdjointEigenSolver<>` and `Eigen::JacobiSVD<>` solvers.
  ///
  /// \cgalModels `DeformationClosestRotationTraits_3`
  class Deformation_Eigen_polar_closest_rotation_traits_3 : 
    public Deformation_Eigen_closest_rotation_traits_3{
  public:

    /// \cond SKIP_FROM_MANUAL

    /// Computes closest rotation to `m` and places it into `R`
    /// Warning: it is adapted from previous experimental code, and not checked deeply
    void compute_closest_rotation(const Matrix& m, Matrix& R)
    {
      CGAL_PROFILER(" times closest rotation is computed");
      bool solved = polar_eigen(m, R);

      if(!solved) { 
        CGAL_PROFILER(" times fallback from polar_eigen failed and SVD is called");
        Deformation_Eigen_closest_rotation_traits_3::compute_closest_rotation(m, R); 
      }       
    }

  private:
    // polar decomposition using Eigen, 5 times faster than SVD
    bool polar_eigen(const Matrix& A, Matrix& R)
    {
      if(A.determinant() < 0)
      { return false; }

      typedef Matrix::Scalar Scalar;

      const Scalar th = std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());

      Eigen::SelfAdjointEigenSolver<Matrix> eig;
      feclearexcept(FE_UNDERFLOW);
      eig.computeDirect(A.transpose()*A);
      if(fetestexcept(FE_UNDERFLOW) || eig.eigenvalues()(0)/eig.eigenvalues()(2)<th)
      { return false; }

      Vector S = eig.eigenvalues().cwiseSqrt();
      R = A  * eig.eigenvectors() * S.asDiagonal().inverse()
        * eig.eigenvectors().transpose();

      if(std::abs(R.squaredNorm()-3.) > th || R.determinant() < 0)
      { return false; }

      R.transposeInPlace(); // the optimal rotation matrix should be transpose of decomposition result
      return true;
    }
    /// \endcond

  };
  /// @endcond
}//namespace CGAL
#endif // CGAL_DEFORMATION_EIGEN_POLAR_CLOSEST_ROTATION_TRAITS_3_H