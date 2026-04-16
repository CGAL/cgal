/*!
\ingroup PkgFeatureGraphConcepts
\cgalConcept

The concept 'Normal_estimator' describes a functor that
extracts the normal for a surface element.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Normal_estimator::AmbrosioTortorelli_on_image}
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

  \tparam Index the type of index of the element to evaluate.

  \param element_type the type of the element (point, line or surface).
  \param element_index the index of the element to evaluate.
  */
  template <typename Index>
  Normal_type operator()(const CGAL::Surface_element_type& element_type, const Index& element_index) const;

/// @}

}; /* end Normal_estimator */