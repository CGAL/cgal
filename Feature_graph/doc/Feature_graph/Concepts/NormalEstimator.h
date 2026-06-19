/*!
*
* \ingroup PkgFeatureGraphConcepts
* \cgalConcept
*
* The concept `NormalEstimator` describes an estimator that
* extracts the normal for a surface element.
*
* \cgalHasModelsBegin
* \cgalHasModels{CGAL::Feature_graph::AmbrosioTortorelli_on_image::Normal_estimator}
* \cgalHasModels{CGAL::Feature_graph::Normal_estimator_on_surface}
* \cgalHasModelsEnd
*
*/

class NormalEstimator {
public:

/// \name Types
/// @{

/*!
* The type of the normal vector model of `Kernel::Vector_3`
*/
typedef unspecified_type Normal_type;

/// @}

/// \name Estimator
/// @{

/*!
* \brief returns the normal vector of the surface element described by a dimension and a descriptor.
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

}; /* end Normal_estimator */