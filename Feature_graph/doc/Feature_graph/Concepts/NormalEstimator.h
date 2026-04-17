/*!
\ingroup PkgFeatureGraphConcepts
\cgalConcept

The concept 'Normal_estimator' describes a functor that
extracts the normal for a surface element.

\cgalHasModelsBegin
\cgalHasModels{CGAL::AmbrosioTortorelli_on_image::Normal_functor}
\cgalHasModels{CGAL::Normal_estimator::Normal_estimator_on_surface}
\cgalHasModelsEnd

*/

class NormalEstimator {
public:

/// \name Types
/// @{

  /*!
  The type of the normal vector model of `Kernel::Vector_3`
  */
  typedef unspecified_type Normal_type;

/// @}

/// \name Functor
/// @{

  /*!
  \brief Returns the normal vector of the surface element described by a type and an index.

  \tparam Element_type_tag a tag that represent the element type.
          Can be `CGAL::Element_type::Point`, `CGAL::Element_type::line` or `CGAL::Element_type::Surface`
  \tparam Index the type of index of the element to evaluate.

  \param element_index the index of the element to evaluate.
  */
  template <typename Element_type_tag, typename Index>
  Normal_type operator()(const Index& element_index) const;

/// @}

}; /* end Normal_estimator */