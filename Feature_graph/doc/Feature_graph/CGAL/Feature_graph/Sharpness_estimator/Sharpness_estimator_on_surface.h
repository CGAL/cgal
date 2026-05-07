namespace CGAL {

namespace Feature_graph {

namespace Sharpness_estimator
{

/*!
* \ingroup PkgFeatureGraphSharpnessEstimator
*
* \brief Functor that assign a sharpness value on elements of a surface.
*
* \cgalModels{SharpnessEstimator}
*/
struct Sharpness_estimator_on_surface
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of of the sharpness value.
  */
  typedef double Sharpness_value_type;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * Constructor that pre-computes the normals on the surface.
  *
  * \tparam Surface a model of `FaceListGraph` that represents a surface mesh.
  * \tparam FT a model of `RealEmbeddable`
  *
  * \param surface the surface where the normals are evaluated.
  * \param selection_threshold a threshold on the sharpness value.
  *     Elements with a sharpness value lower than this threshold are considered flat
  *     and will be given a negative value.
  */
  template <typename Surface, typename FT = Sharpness_value_type>
  Sharpness_estimator_on_surface(const Surface& surface, const FT& selection_threshold);

  /// @}

  /// \name Functor
  /// @{

  /*!
  * returns the sharpness value of the surface element described by a type and an index.
  *
  * \tparam DimensionTag a tag that represent the element type.
  *         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
  * \tparam Index the type of index of the element to evaluate.
  *
  * \param element_index the index of the element to evaluate.
  */
  template <typename DimensionTag, typename Index>
  Sharpness_value_type operator()(const Index& element_index) const;

  /// @}
};

} /* namespace Sharpness_estimator */

} /* namespace Feature_graph */

} /* namespace CGAL */