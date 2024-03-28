/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\cgalRefines{DefaultConstructible, CopyConstructible, Assignable}

The concept `IsosurfacingInterpolationScheme_3` describes the set of requirements to be fulfilled
by the interpolation scheme template parameter of the domain classes `CGAL::Isosurfacing::Interpolated_discrete_values_3` and `CGAL::Isosurfacing::Interpolated_discrete_gradients_3`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Trilinear_interpolation}
\cgalHasModelsEnd
*/
class IsosurfacingInterpolationScheme_3
{
public:
  /*!
  * The geometric traits type.
  * Must be a model of `IsosurfacingTraits_3`, equal to that of the grid.
  */
  typedef unspecified_type Geom_traits;

  /*!
  * The scalar type.
  */
  typedef unspecified_type FT;

  /*!
  * The 3D point type.
  */
  typedef unspecified_type Point_3;

  /*!
  * The 3D vector type.
  */
  typedef unspecified_type Vector_3;

  /*!
  * \brief interpolates the value of the value field at the point `p` using values defined by `vr`
  * over the grid `g`.
  *
  * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3`
  * \tparam ValueRange must be a model of `RandomAccessRange` with `FT` as value type
  */
  template <typename Grid, typename ValueRange>
  FT interpolate_values(const Point_3& p, const Grid& g, const ValueRange& vr) const;

  /*!
   * \brief interpolates the gradient of the gradient field at the point `p`
   *        using gradients defined by `gr` over the grid `g`.
   *
   * \tparam Grid must be `CGAL::Isosurfacing::Cartesian_grid_3`
   * \tparam GradientRange must be a model of `RandomAccessRange` with `Vector_3` as value type
   */
  template <typename Grid, typename GradientRange>
  Vector_3 interpolate_gradients(const Point_3& p, const Grid& g, const GradientRange& gr) const;
};
