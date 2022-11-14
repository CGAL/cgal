/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `PlacementFilter` describes the requirements for the <I>policy
function object</I> which gets the profile and placement of an edge
and which can filter the placement.  The filtering is only done when
an edge is taken from the priority queue in order to get collapsed,
and neither when the edge is inserted nor when it is updated in the
priority queue.

The placement returned is a `boost::optional` value (i.e., it can
be absent). The value `boost::none` indicates that the edge should not be collapsed.

\cgalRefines `DefaultConstructible`
\cgalRefines `CopyConstructible`

\cgalHasModel `CGAL::Surface_mesh_simplification::Bounded_normal_change_filter<BasePlacementFilter>`
\cgalHasModel `CGAL::Surface_mesh_simplification::Polyhedral_envelope_filter<GeomTraits, BasePlacementFilter>`
*/


class PlacementFilter
{
public:

  /// The class `Edge_profile` regroups useful information about an edge, such as its incident vertices and faces.
  typedef CGAL::Surface_mesh_simplification::Edge_profile Edge_profile;

  /// \name Operations
  /// @{

  /*!
  filters the placement.

  */
  boost::optional<Edge_profile::Point> operator()(const Edge_profile& profile, boost::optional<Edge_profile::Point> placement) const;

  /// @}

}; /* end PlacementFilter */
