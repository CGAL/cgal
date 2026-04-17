namespace CGAL {

namespace Normal_estimator
{

/*!
\ingroup PkgFeatureGraphNormalEstimator

\brief Functor that assign a normal on elements of a surface.

\tparam Vector_3 the type of the normal vector model of `Kernel::Vector_3`.

\cgalModels{NormalEstimator}
*/
template <typename Vector_3>
struct Normal_estimator_on_surface
{
public:
  /// \name Types
  /// @{

  /*!
  The type of the normal vector.
  */
  typedef Vector_3 Normal_type;

  ///@}

  /// \name Constructor
  /// @{

  /*!
  Constructor that pre-computes the normals on the surface.

  \tparam Surface a model of `FaceListGraph` that represents a surface mesh.

  \param surface the surface where the normals are evaluated.
  */
  template <typename Surface>
  Normal_estimator_on_surface(const Surface& surface);
  
  // @}

  /// \name Functor
  /// @{

  /*!
  Returns the normal vector of the surface element described by a type and an index.

  \tparam Element_type_tag a tag that represent the element type.
          Can be `CGAL::Element_type::Point`, `CGAL::Element_type::line` or `CGAL::Element_type::Surface`
  \tparam Index the type of index of the element to evaluate.

  \param element_index the index of the element to evaluate.
  */
  template <typename Element_type_tag, typename Index>
  Normal_type operator()(const Index& element_index) const;

  /// @}
};

} /* namespace Normal_estimator */

} /* namespace CGAL */