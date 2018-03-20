/*!

  \ingroup PkgSetMovableSeparability2Concepts
  \cgalConcept

  This concept generalizes the concept of a 2D Kernel.

  \cgalRefines `DefaultConstructible`
  \cgalRefines `PolygonTraits_2`

  \cgalHasModel Any CGAL kernel, e.g., CGAL::Exact_predicates_exact_constructions_kernel.

*/

class CastingTraits_2 {
public:

  /// \name Types
  /// @{

  //! The direction type. Models the concept `Kernel::Direction_2`.
  typedef unspecified_type Direction_2;

  /// @}

  /// \name Functor Types
  /// @{

  //! Models the concept `Kernel::Counterclockwise_in_between_2`.
  typedef unspecified_type Counterclockwise_in_between_2;

  //! Models the concept `Kernel::Collinear_2`.
  typedef unspecified_type Collinear_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  Counterclockwise_in_between_2 counterclockwise_in_between_2_object() const;
  Collinear_2 collinear_2_object() const;

  /// @}

};
