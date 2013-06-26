namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

Simplifies `surface` in-place by collapsing edges, and returns
the number of edges effectively removed.

The function `Surface_mesh_simplification::edge_collapse` simplifies in-place a triangulated surface mesh by iteratively collapsing edges. 

\cgalHeading{Non-Named Parameters}

`surface` is the surface to simplify. 
It must be a model of the `EdgeCollapsableMesh` concept. 

`should_stop` is the stop-condition policy. 
It must be a model of the `StopPredicate` concept. 

\cgalHeading{Named Parameters}

`named_parameters` holds the list of all the additional parameters 
used by the `edge_collapse` function (including default parameters). 

The named parameters list is a composition of function calls separated by a dot (\f$ .\f$) where 
the name of each function matches the name of an argument and wraps the actual parameter. 

This is an example with 2 arguments: 

`vertex_index_map(the_actual_vertex_index_map).edge_index_map(the_actual_edge_index_map)` 

`the_actual_vertex_index_map` and `the_actual_edge_index_map` are 
the actual parameters, while `vertex_index_map()` and `edge_index_map()` 
are wrapper functions used to designate each formal argument. 

All named parameters have default values so you only need to compose those for which the default 
is inappropriate. Furthermore, since each actual parameter is wrapped in a function whose name 
designates the formal argument, the order of named parameters in the list is totally irrelevant. 

In the following subsections, each named parameter is documented as a helper function. The argument to each helper 
function is the actual parameter to `edge_collapse()`, while the name of the helper 
function designates which formal argument it is. 

\cgalHeading{vertex_index_map(VertexIndexMap vpm)}

Maps each vertex in the surface into an unsigned integer number 
in the range `[0,num_vertices(surface))`. 

`VertexIndexMap` must be a model of
`ReadablePropertyMap` 
whose `key_type` is 
`boost::graph_traits<EdgeCollapsableMesh const>::%vertex_descriptor` 
and whose `value_type` is 
`boost::graph_traits<EdgeCollapsableMesh>::%size_type`, 

<B>%Default</B>: the property map obtained by calling `get(vertex_index,surface)`, 
which requires the surface vertices to have an `id()` member properly initialized to the 
required value. 

If the vertices don't have such an `id()`, you must pass some property map explicitly. 
An external property map can be easily obtained by calling 
`get(vertex_external_index,surface)`. This constructs on the fly, and returns, 
a property map which non-intrusively associates a proper id with each vertex. 

\cgalHeading{edge_index_map(EdgeIndexMap eim)}

Maps each <I>directed</I> edge in the surface into an unsigned integer number 
in the range `[0,num_edges(surface))`. 

`EdgeIndexMap` must be a model of
`ReadablePropertyMap` whose `key_type` is 
`boost::graph_traits<EdgeCollapsableMesh const>::%edge_descriptor` 
and whose `value_type` is 
`boost::graph_traits<EdgeCollapsableMesh>::%size_type` 

<B>%Default</B>: the property map obtained by calling `get(edge_index,surface)`, 
which requires the surface edges to have an `id()` member properly initialized to the 
require value. 

If the edges don't have such an `id()`, you must pass some property map explicitly. 
An external property map can be easily obtained by calling 
`get(edge_external_index,surface)`. This constructs on the fly, and returns, 
a property map which non-intrusively associates a proper id with each edge. 

\cgalHeading{edge_is_border_map(EdgeIsBorderMap ebm)}

Maps each <I>directed</I> edge in the surface into a Boolean value 
which indicates if the edge belongs to the boundary of the surface 
(facing the outside). 
`EdgeIsBorderMap` must be a model
`ReadablePropertyMap` whose `key_type` is 
`boost::graph_traits<EdgeCollapsableMesh const>::%edge_descriptor` 
and whose `value_type` is `bool`. 

<B>%Default</B>: the property map obtained by calling `get(edge_is_border,surface)`. 

\cgalHeading{get_cost(GetCost gc)}

The policy which returns the collapse cost for an edge. 

The type of `gc` must be a model of the `GetCost` concept. 

<B>%Default</B>: 
`CGAL::Surface_mesh_simplification::LindstromTurk_cost<EdgeCollapsableMesh>`. 

\cgalHeading{get_placement(GetPlacement gp)}

The policy which returns the placement (position of the replacemet vertex) 
for an edge. 

The type of `gp` must be a model of the `GetPlacement` concept. 

<B>%Default</B>: 
`CGAL::Surface_mesh_simplification::LindstromTurk_placement<EdgeCollapsableMesh>` 

\cgalHeading{visitor(EdgeCollapseSimplificationVisitor v)}

The visitor that is called by the `edge_collapse` function 
in certain points to allow the user to track the simplification process. 

The type of `v` must be a model of the `EdgeCollapseSimplificationVisitor` concept. 

<B>%Default: an implementation-defined dummy visitor</B>. 

If you wish to provide your own visitor, you can derive from: 
`CGAL::Surface_mesh_simplification::Edge_collapse_visitor_base<EdgeCollapsableMesh>` 
and override only the callbacks you are interested in. 

All these functions naming parameters are defined in 
`namespace CGAL`. Being non-member functions, they could clash 
with equally named functions in some other namespace. If that happens, 
simply qualify the <I>first</I> 
\cgalFootnote{The second and subsequent named parameters shall not be qualified as they are member functions} 
named parameter with `CGAL::`, as shown in the examples in the user manual. 

\cgalHeading{Semantics}


The simplification process continues until the `should_stop` policy returns `true` 
or the surface cannot be simplified any further due to topological constraints. 

`get_cost` and `get_placement` are the policies which control 
the <I>cost-strategy</I>, that is, the order in which edges are collapsed 
and the remaining vertex is re-positioned. 

`visitor` is used to keep track of the simplification process. It has several member functions which 
are called at certain points in the simplification code. 

*/

template<class EdgeCollapsableMesh,class StopPredicate, class P, class T, class R>
int edge_collapse ( EdgeCollapsableMesh& surface
, StopPredicate const& should_stop
, sms_named_params<P,T,R> const& named_parameters
) ;

} /* namespace Surface_mesh_simplification */
} /* namespace CGAL */

