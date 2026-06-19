namespace CGAL {

namespace Feature_graph {

/*!
* \ingroup PkgFeatureGraphSharpnessEstimator
*
* \brief Estimator that assign a sharpness value on elements of a surface.
*
* \cgalModels{SharpnessEstimator}
*/
struct Sharpness_estimator_on_surface
{
public:
  /// \name Types
  /// @{

  /*!
  * The type of the sharpness value.
  */
  typedef double Sharpness_value_type;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * Constructor that pre-computes the normals on the surface.
  *
  * \tparam PolygonMesh a model of `FaceListGraph` that represents a surface.
  * \tparam FT a model of `RealEmbeddable`
  *
  * \param pmesh the surface where the normals are evaluated.
  * \param selection_threshold a threshold on the sharpness value.
  *     Elements with a sharpness value lower than this threshold are considered flat
  *     and will be given a negative value.
  */
  template <typename PolygonMesh, typename FT = Sharpness_value_type>
  Sharpness_estimator_on_surface(const PolygonMesh& pmesh, const FT& selection_threshold);

  /// @}

  /// \name Estimator
  /// @{

  /*!
  * returns the sharpness value of the surface element described by a type and a descriptor.
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
  Sharpness_value_type operator()(const Descriptor& element_descriptor, const Domain& domain) const;

  /// @}
};

} /* namespace Feature_graph */

} /* namespace CGAL */