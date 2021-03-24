namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplificationRef

Simplifies `tmesh` in-place by collapsing edges, and returns
the number of edges effectively removed.

@tparam TriangleMesh a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.
@tparam StopPolicy a model of `StopPredicate`
@tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

@param tmesh a triangle mesh
@param should_stop the stop-condition policy
@param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

\cgalNamedParamsBegin
  \cgalParamNBegin{vertex_point_map}
    \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
    \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                   as key type and `%Point_3` as value type}
    \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                    should be available for the vertices of `tmesh`.}
  \cgalParamNEnd

  \cgalParamNBegin{geom_traits}
    \cgalParamDescription{an instance of a geometric traits class}
    \cgalParamType{a class model of `Kernel`}
    \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
    \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
  \cgalParamNEnd

   \cgalParamNBegin{halfedge_index_map}
     \cgalParamDescription{a property map associating to each halfedge of `tmesh` a unique index between `0` and `num_halfedges(tmesh) - 1`}
     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%halfedge_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{an automatically indexed internal map}
     \cgalParamExtra{If this parameter is not passed, internal machinery will create and initialize
                     a face index property map, either using the internal property map if it exists
                     or using an external map. The latter might result in  - slightly - worsened performance
                     in case of non-constant complexity for index access.}
   \cgalParamNEnd

  \cgalParamNBegin{get_cost}
    \cgalParamDescription{a policy which returns the collapse cost for an edge}
    \cgalParamType{a model of the concept `GetCost`}
    \cgalParamDefault{`CGAL::Surface_mesh_simplification::LindstromTurk_cost<TriangleMesh>`}
  \cgalParamNEnd

  \cgalParamNBegin{get_placement}
    \cgalParamDescription{a policy which returns the placement (position of the replacemet vertex) for an edge}
    \cgalParamType{a model of the concept `GetPlacement`}
    \cgalParamDefault{`CGAL::Surface_mesh_simplification::LindstromTurk_placement<TriangleMesh>`}
  \cgalParamNEnd

  \cgalParamNBegin{filter}
    \cgalParamDescription{a policy which returns the filter for a placement}
    \cgalParamType{a model of the concept `Filter`}
    \cgalParamDefault{no placement gets filtered}
  \cgalParamNEnd

  \cgalParamNBegin{edge_is_constrained_map}
    \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tmesh`}
    \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
                   as key type and `bool` as value type}
    \cgalParamDefault{a constant property map returning `false` for any edge key}
    \cgalParamExtra{A constrained edge cannot be collapsed.}
  \cgalParamNEnd

  \cgalParamNBegin{visitor}
    \cgalParamDescription{a visitor that is called by the `edge_collapse` function
                          at certain points to allow the user to track the simplification process}
    \cgalParamType{a class model of the concept `EdgeCollapseSimplificationVisitor`}
    \cgalParamDefault{unused}
  \cgalParamNEnd

   \cgalParamNBegin{vertex_index_map}
     \cgalParamDescription{a property map associating to each vertex of `tmesh` a unique index between `0` and `num_vertices(tmesh) - 1`}
     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
                    as key type and `std::size_t` as value type}
     \cgalParamDefault{an automatically indexed internal map}
     \cgalParamExtra{If this parameter is not passed, internal machinery will create and initialize
                     a face index property map, either using the internal property map if it exists
                     or using an external map. The latter might result in  - slightly - worsened performance
                     in case of non-constant complexity for index access.}
     \cgalParamExtra{This parameter is only used by debug functions and is usually not needed for users.}
   \cgalParamNEnd
\cgalNamedParamsEnd

\cgalHeading{Semantics}

The simplification process continues until the `should_stop` policy returns `true`
or until the surface mesh cannot be simplified any further due to topological constraints.

`get_cost` and `get_placement` are the policies which control
the <I>cost-strategy</I>, that is, the order in which edges are collapsed
and the remaining vertex is re-positioned.

`visitor` is used to keep track of the simplification process. It has several member functions which
are called at certain points in the simplification code.
*/
template<class TriangleMesh, class StopPolicy, class NamedParameters>
int edge_collapse(TriangleMesh& tmesh,
                  const StopPolicy& should_stop,
                  const NamedParameters& np);

} // namespace Surface_mesh_simplification
} /* namespace CGAL */
