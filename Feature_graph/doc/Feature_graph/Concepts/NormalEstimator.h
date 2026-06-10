/*!
*
* \ingroup PkgFeatureGraphConcepts
* \cgalConcept
*
* The concept `NormalEstimator` describes a functor that
* extracts the normal for a surface element.
*
* \cgalHasModelsBegin
* \cgalHasModels{CGAL::Feature_graph::AmbrosioTortorelli_on_image::Normal_functor}
* \cgalHasModels{CGAL::Feature_graph::Normal_estimator::Normal_estimator_on_surface}
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

/// \name Functor
/// @{

/*!
* \brief returns the normal vector of the surface element described by a dimension and an index.
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

}; /* end Normal_estimator */