namespace CGAL {

namespace Normal_estimator
{

/*!
\ingroup PkgFeatureGraphNormalEstimator

\brief Functor that assign the Ambrosio-Tortorelli normals on elements of an image.

As the Ambrosio-Tortorelli estimate the normals and the sharpness together,
this structure only recovers the Normal values from a previous computation.

\tparam Vector_3 the type of the normal vector.

\cgalModels{NormalEstimator}
*/
template <typename Vector_3>
struct AmbrosioTortorelli_on_image
{
public:
  /// \name Types
  /// @{

  /*!
  The type of the normal vector model of `Kernel::Vector_3`
  */
  typedef Vector_3 Normal_type;

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
  Returns the normal of the surface element described by a type and an index.

  \tparam Index the type of index of the element to evaluate.

  \param element_type the type of the element (point, line or surface).
  \param element_index the index of the element to evaluate.
  */
  template <typename Index>
  Normal_type operator()(const CGAL::Surface_element_type& element_type, const Index& element_index) const;

  ///@}
};

} /* namespace Normal_estimator */

} /* namespace CGAL */