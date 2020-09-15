/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `GetFilter` describes the requirements for the <I>policy
function object</I> which gets the profile and placement of an edge.

The placement returned is a `boost::optional` value (i.e., it can
be absent). The value `boost::none` indicates that the edge should not be collapsed.

\cgalRefines `DefaultConstructible`
\cgalRefines `CopyConstructible`

\cgalHasModel `CGAL::Surface_mesh_simplification::Bounded_normal_change_filter<Placement>`
\cgalHasModel `CGAL::Surface_mesh_simplification::FastEnvelopeFilter<GeomTraits,Placement>`
*/


class GetFilter
{
public:

  /// The class Edge_profile regroups useful information about an edge, such as its incident vertices and faces.
  typedef CGAL::Surface_mesh_simplification::Edge_profile Edge_profile;

  /// \name Operations
  /// @{

  /*!
  filters the placement.

  */
  boost::optional<Edge_profile::Point> operator()(const Edge_profile& profile, boost::optional<Edge_profile::Point> placement) const;

  /// @}

}; /* end Filter */
