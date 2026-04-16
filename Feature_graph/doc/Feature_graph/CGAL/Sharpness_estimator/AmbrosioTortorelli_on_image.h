namespace CGAL {

namespace Sharpness_estimator
{

/*!
\ingroup PkgFeatureGraphSharpnessEstimator

\brief Functor that assign the Ambrosio-Tortorelli sharpness values on elements of an image.

As the Ambrosio-Tortorelli estimate the normals and the sharpness together,
this structure only recovers the sharpness values from a previous computation.

\cgalModels{SharpnessEstimator}
*/
struct AmbrosioTortorelli_on_image
{
public:
  /// \name Types
  /// @{

  /*!
  The type of the normal vector model of `Kernel::Vector_3`
  */
  typedef double Sharpness_value_type;

  ///@}

  /// \name Constructor
  /// @{

  /*!
  Constructor that uses a previous computation of the Ambrosio-Tortorelli energy.

  \param ambrosioTortorelli_computation_result the computation result of the Ambrosio-Tortorelli energy on an image.
  */
  AmbrosioTortorelli_on_image(const unspecified_type& ambrosioTortorelli_computation_result);

  // @}

  /// \name Functor
  /// @{
  /*!
  Returns the sharpness value of the surface element described by a type and an index.

  \tparam Index the type of index of the element to evaluate.

  \param element_type the type of the element (point, line or surface).
  \param element_index the index of the element to evaluate.
  */
  template <typename Index>
  Sharpness_value_type operator()(const CGAL::Surface_element_type& element_type, const Index& element_index) const;

  ///@}
};
} /* namespace Sharpness_estimator */

} /* namespace CGAL */