namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

Simplifies `surface_mesh` in-place by collapsing edges, and returns
the number of edges effectively removed.

@tparam TriangleMesh a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
@tparam StopPolicy a model of `StopPredicate`
@tparam NamedParameters a sequence of \ref sms_namedparameters "Named Parameters"

@param surface_mesh a triangle mesh
@param should_stop the stop-condition policy
@param np optional sequence of \ref sms_namedparameters "Named Parameters" among the ones listed below

\cgalNamedParamsBegin
  \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the mesh.
     If this parameter is omitted, an internal property map for
     `CGAL::vertex_point_t` should be available in `PolygonMesh`
  \cgalParamEnd

  \cgalParamBegin{halfedge_index_map} the property map containing an index for each halfedge,
    initialized 0 to `num_halfedges(graph)`
  \cgalParamEnd

  \cgalParamBegin{get_cost}
    The policy which returns the collapse cost for an edge.
  \cgalParamEnd

  \cgalParamBegin{get_placement}
    The policy which returns the placement (position of the replacemet vertex) for an edge.
  \cgalParamEnd

  \cgalParamBegin{edge_is_constrained_map}
    The property map containing the constrained-or-not status of each edge of `pmesh`
  \cgalParamEnd

  \cgalParamBegin{visitor}
    The visitor that is called by the `edge_collapse` function
    in certain points to allow the user to track the simplification process.
  \cgalParamEnd

  \cgalParamBegin{vertex_index_map}
    is the property map containing the index of each vertex of the input polygon mesh.
    \cgalDebugBegin
    This parameter is only used by debug functions and is usually not needed for users.
    \cgalDebugEnd
  \cgalParamEnd
\cgalNamedParamsEnd

\cgalHeading{Semantics}

The simplification process continues until the `should_stop` policy returns `true` 
or the surface mesh cannot be simplified any further due to topological constraints. 

`get_cost` and `get_placement` are the policies which control 
the <I>cost-strategy</I>, that is, the order in which edges are collapsed 
and the remaining vertex is re-positioned. 

`visitor` is used to keep track of the simplification process. It has several member functions which 
are called at certain points in the simplification code. 
*/
template<class TriangleMesh, class StopPolicy, class NamedParameters>
int edge_collapse(TriangleMesh& surface_mesh,
                  const StopPolicy& should_stop,
                  const NamedParameters& np);

} /* namespace Surface_mesh_simplification */
} /* namespace CGAL */

