namespace CGAL {

namespace Sharpness_estimator
{

/*!
\ingroup PkgFeatureGraphSharpnessEstimator

\brief Functor that assign a sharpness value on elements of a surface.

\cgalModels{SharpnessEstimator}
*/
struct Sharpness_estimator_on_surface
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
  Constructor that pre-computes the normals on the surface.

  \tparam Surface a model of `FaceListGraph` that represents a surface mesh.

  \param surface the surface where the normals are evaluated.
  */
  template <typename Surface>
  Sharpness_estimator_on_surface(const Surface& surface);
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

  /// @}
};

} /* namespace Sharpness_estimator */

} /* namespace CGAL */