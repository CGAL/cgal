namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphNormalEstimator
*
* \brief Estimator that assign a normal on elements of a surface.
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

  /// \name Estimator
  /// @{

  /*!
  * returns the normal vector of the surface element identified by a dimension and a descriptor.
  *
  * \tparam DimensionTag a tag that represent the element type.
  *         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
  * \tparam Descriptor the type of descriptor used to identify the element to evaluate,
  * which can be a vertex, an edge, or a facet according to the DimensionTag.
  * If the domain is a model of `FeatureImage_3`, it is a `std::size_t` for element with
  * dimension 0, 1 and 2. If the domain is a model of `FaceListGraph`, it is a
  * `vertex_descriptor` (resp. `halfedge_descriptor`; `face_descriptor `) for element with
  * dimension 0 (resp. 1 ; 2).
  * \tparam Domain the type of the surface where the element is embedded.
  * Can be an Image class model of `FeatureImage_3`, or a model of `FaceListGraph` that represents a surface mesh.
  *
  * \param element_descriptor the descriptor of the element on the surface.
  * \param domain the domain that contains the elements.
  */
  template <typename DimensionTag, typename Descriptor, typename Domain>
  Normal_type operator()(const Descriptor& element_descriptor, const Domain& domain) const;

  /// @}
};

} /* namespace Feature_graph */

} /* namespace CGAL */