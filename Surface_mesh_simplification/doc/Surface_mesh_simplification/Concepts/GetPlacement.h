/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `GetPlacement` describes the requirements for the <I>policy
function object</I> which gets the <I>collapse placement</I> of an edge,
that is, the new position of the vertex that remains after a
halfedge-collapse operation.

The placement returned is a `std::optional` value (i.e., it can
be absent). An absent result indicates that the edge should not be collapsed.
This could be the result of a computational limitation (such as an overflow),
or can be intentionally returned to prevent the edge from being collapsed.

\cgalRefines{DefaultConstructible,CopyConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Surface_mesh_simplification::Midpoint_placement<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_simplification::LindstromTurk_placement<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_simplification::GarlandHeckbert_policies<TriangleMesh, GeomTraits>}
\cgalHasModels{CGAL::Surface_mesh_simplification::Bounded_normal_change_placement<Placement>}
\cgalHasModels{CGAL::Surface_mesh_simplification::Constrained_placement<Placement>}
\cgalHasModelsEnd
*/


class GetPlacement
{
public:

  /// The class `Edge_profile` regroups useful information about an edge, such as its incident vertices and faces.
  typedef CGAL::Surface_mesh_simplification::Edge_profile Edge_profile;

  /// \name Operations
  /// @{

  /*!
  Computes and returns the placement, that is, the position of the vertex
  which replaces the collapsing edge (represented by its profile).
  */
  std::optional<Edge_profile::Point>operator()(const Edge_profile& profile) const;

  /// @}

}; /* end GetPlacement */
