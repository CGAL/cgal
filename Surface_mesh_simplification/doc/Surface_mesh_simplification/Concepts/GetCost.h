/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `GetCost` describes the requirements for the <I>policy function object</I>
which gets the <I>collapse cost</I> of an edge.

The cost returned is a `std::optional` value (i.e.\ it can be absent).
An <I>absent</I> cost indicates that the edge should not be collapsed.
This could be the result of a computational limitation (such as an overflow),
or can be intentionally returned to prevent the edge from being collapsed.

\cgalRefines{DefaultConstructible,CopyConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Surface_mesh_simplification::Edge_length_cost<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_simplification::LindstromTurk_cost<TriangleMesh>}
\cgalHasModels{CGAL::Surface_mesh_simplification::GarlandHeckbert_policies<TriangleMesh, GeomTraits>}
\cgalHasModelsEnd

*/
class GetCost
{
public:

  /// The class `Edge_profile` regroups useful information about an edge, such as its incident vertices and faces.
  typedef CGAL::Surface_mesh_simplification::Edge_profile Edge_profile;

  /// \name Operations
  /// @{

  /*!
  Computes and returns the cost of collapsing the edge (represented by its profile),
  using the calculated placement.
  */
  std::optional<typename Edge_profile::FT> operator()(const Edge_profile& edge_profile,
                                                        const std::optional<typename Edge_profile::Point>& placement) const;

/// @}

}; /* end GetCost */
