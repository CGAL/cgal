/*! \ingroup PkgArrangementOnSurface2ConceptsTopologyTraits
 * \cgalConcept
 *
 * The concept `ArrangementBasicTopologyTraits` defines the minimal
 * functionality needed for a model of a topology traits, which can substitutes
 * the `TopolTraits` template parameters when the class template
 * `Arrangement_on_surface_2<GeomTraits, TopolTraits>` is instantiated.  In
 * particular. a model of this concept holds the Dcel data structure used to
 * represent the arrangement cells (i.e., vertices, edges, and facets) and the
 * incident relations between them.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_spherical_topology_traits_2<GeometryTraits_2, Dcel>}
 * \cgalHasModelsEnd
 */

class ArrangementBasicTopologyTraits {
public:
  /// \name Types
  /// @{

  //! models the concept `ArrTraits::Point_2`.
  typedef unspecified_type                              Point_2;

  //! models the concept `ArrTraits::XMonotoneCurve_2`.
  typedef unspecified_type                              X_monotone_curve_2;

  //! models the concept `ArrangementDcel`.
  typedef unspecified_type                              Dcel;

  /// @}

  /// \name Creation
  /// @{

  /// @}

  /// \name Access Functions
  /// @{

  /*! Obtain the DCEL (const version). */
  const Dcel& dcel() const;

  /*! Obtain the DCEL (non-const version). */
  Dcel& dcel();

  /// @}

  /// \name Modifiers
  /// @{
  /// @}

};
