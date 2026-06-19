/*!
* \ingroup PkgFeatureGraphConcepts
* \cgalConcept
*
* The concept `Sharpness_estimator` describes an estimator that
* extracts the sharpness value for a surface element.
*
* \cgalHasModelsBegin
* \cgalHasModels{CGAL::Feature_graph::AmbrosioTortorelli_on_image::Sharpness_estimator}
* \cgalHasModels{CGAL::Feature_graph::Sharpness_estimator_on_surface}
* \cgalHasModelsEnd
*
*/

class SharpnessEstimator {
public:

/// \name Types
/// @{

/*!
* \brief The type of the sharpness value.
* It is a scalar constructible from the integer 0.
* \cgalModels{RealEmbeddable}
*/
typedef unspecified_type Sharpness_value_type;

/// @}

/// \name Estimator
/// @{

/*!
* \brief returns the sharpness value of the surface element described by a type and a descriptor.
* A low sharpness value should represent a flat area,
* while a high value implies a sharp feature.
* For two sharpness values `A` and `B`,
* the value `A` is told to have a higher sharpness value iff `B < A`.
* Negative values represent smooth areas that should be erased in the
* selection step of the feature graph extraction method.
*
* \tparam DimensionTag a tag that represent the element type.
*         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
* \tparam Descriptor the type of descriptor used to identify the element to evaluate,
* which can be a vertex, an edge, or a facet according to the DimensionTag.
* If the domain is of type `CGAL::Image_3`,, it is a `std::size_t` for element with
* dimension 0, 1 and 2. If the domain is a model of `FaceListGraph`, it is a
* `vertex_descriptor` (resp. `halfedge_descriptor`; `face_descriptor `) for element with
* dimension 0 (resp. 1 ; 2).
* \tparam Domain the type of the surface where the element is embedded.
* Can be a `CGAL::Image_3`, or a model of `FaceListGraph` that represents a surface mesh.
*
* \param element_descriptor the descriptor of the element on the surface.
* \param domain the domain that contains the elements.
*/
template <typename DimensionTag, typename Descriptor, typename Domain>
Sharpness_value_type operator()(const Descriptor& element_descriptor, const Domain& domain) const;

/// @}

}; /* end Sharpness_estimator */