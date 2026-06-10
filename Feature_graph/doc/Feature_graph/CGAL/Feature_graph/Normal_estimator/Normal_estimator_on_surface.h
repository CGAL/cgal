namespace CGAL {

namespace Feature_graph {

namespace Normal_estimator
{

/*!
* \ingroup PkgFeatureGraphNormalEstimator
*
* \brief Functor that assign a normal on elements of a surface.
*
* \tparam Vector_3 the type of the normal vector model of `Kernel::Vector_3`.
*
* \cgalModels{NormalEstimator}
*/
template <typename Vector_3>
struct Normal_estimator_on_surface
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of the normal vector.
  */
  typedef Vector_3 Normal_type;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * Constructor that pre-computes the normals on the surface.
  *
  * \tparam PolygonMesh a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param pmesh the surface where the normals are evaluated.
  */
  template <typename PolygonMesh>
  Normal_estimator_on_surface(const PolygonMesh& pmesh);

  /// @}

  /// \name Functor
  /// @{

  /*!
  * returns the normal vector of the surface element identified by a dimension and an index.
  *
  * \tparam DimensionTag a tag that represent the element type.
  *         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
  * \tparam Index the type of index used to identify the element to evaluate,
  * which can be a vertex, an edge, or a facet according to the DimesionTag.
  *
  * \param element_index the index of the element to evaluate.
  */
  template <typename DimensionTag, typename Index>
  Normal_type operator()(const Index& element_index) const;

  /// @}
};

} /* namespace Normal_estimator */

} /* namespace Feature_graph */

} /* namespace CGAL */