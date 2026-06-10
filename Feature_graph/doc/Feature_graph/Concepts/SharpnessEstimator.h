/*!
* \ingroup PkgFeatureGraphConcepts
* \cgalConcept
*
* The concept `Sharpness_estimator` describes a functor that
* extracts the sharpness value for a surface element.
*
* \cgalHasModelsBegin
* \cgalHasModels{CGAL::Feature_graph::AmbrosioTortorelli_on_image::Sharpness_functor}
* \cgalHasModels{CGAL::Feature_graph::Sharpness_estimator::Sharpness_estimator_on_surface}
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

/// \name Functor
/// @{

/*!
* \brief returns the sharpness value of the surface element described by a type and an index.
* A low sharpness value should represent a flat area,
* while a high value implies a sharp feature.
* For two sharpness values `A` and `B`,
* the value `A` is has a higher sharpness value iff `B < A`.
* Negative values represent smooth areas that should be erased in the
* selection step of the feature graph extraction method.
*
* \tparam DimensionTag a tag that represent the element type.
*         Can be `CGAL::Dimension_tag<0>`, `CGAL::Dimension_tag<1>` or `CGAL::Dimension_tag<2>`
* \tparam Index the type of index used to identify the element to evaluate,
* which can be a vertex, an edge, or a facet according to the DimesionTag.
*
* \param element_index the index of the element to evaluate.
*/
template <typename DimensionTag, typename Index>
Sharpness_value_type operator()(const Index& element_index) const;

/// @}

}; /* end Sharpness_estimator */