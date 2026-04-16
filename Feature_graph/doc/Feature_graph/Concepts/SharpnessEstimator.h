/*!
\ingroup PkgFeatureGraphConcepts
\cgalConcept

The concept 'Sharpness_estimator' describes a functor that
extracts the sharpness value for a surface element.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Sharpness_estimator::AmbrosioTortorelli_on_image}
\cgalHasModels{CGAL::Sharpness_estimator::Sharpness_estimator_on_surface}
\cgalHasModelsEnd

*/

class SharpnessEstimator {
public:

/// \name Types
/// @{

  /*!
  The type of the sharpness value.
  It must implement `bool Sharpness_value_type::operator<(const Sharpness_value_type&) const`
  */
  typedef unspecified_type Sharpness_value_type;

/// @}

/// \name Functor
/// @{

  /*!
  \brief Returns the sharpness value of the surface element described by a type and an index.
  A low sharpness value should represent a flat area,
  while a high value implies a sharp feature.
  For two sharpness values `A` and `B`,
  the value `A` is has a higher sharpness value iff `B < A`.

  \tparam Index the type of index of the element to evaluate.

  \param element_type the type of the element (point, line or surface).
  \param element_index the index of the element to evaluate.
  */
  template <typename Index>
  Sharpness_value_type operator()(const CGAL::Surface_element_type& element_type, const Index& element_index) const;

/// @}

}; /* end Sharpness_estimator */