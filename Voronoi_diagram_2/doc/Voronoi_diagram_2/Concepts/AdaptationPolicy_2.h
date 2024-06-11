
/*!
\ingroup PkgVoronoiDiagram2Concepts
\cgalConcept

The concept `AdaptationPolicy_2` defines the requirements on
the predicate functors that determine whether a feature of the
triangulated Delaunay graph should be rejected or not. It also
provides a functor for inserting sites in the Delaunay graph. The last
functor is optional and a tag determines whether it is provided or
not. Note that while the first two functors do not modify the Delaunay
graph they take as an argument, the last ones does.

\cgalRefines{CopyConstructible,Assignable,DefaultConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Identity_policy_2<DG,AT>}
\cgalHasModels{CGAL::Apollonius_graph_degeneracy_removal_policy_2<AG2>}
\cgalHasModels{CGAL::Apollonius_graph_caching_degeneracy_removal_policy_2<AG2>}
\cgalHasModels{CGAL::Delaunay_triangulation_degeneracy_removal_policy_2<DT2>}
\cgalHasModels{CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT2>}
\cgalHasModels{CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT2>}
\cgalHasModels{CGAL::Regular_triangulation_caching_degeneracy_removal_policy_2<RT2>}
\cgalHasModels{CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2>}
\cgalHasModels{CGAL::Segment_Delaunay_graph_caching_degeneracy_removal_policy_2<SDG2>}
\cgalHasModelsEnd

\sa `DelaunayGraph_2`
\sa `AdaptationTraits_2`
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>`

*/

class AdaptationPolicy_2 {
public:

/// \name Types
/// @{

/*!
A type for the sites of the Voronoi diagram.
*/
typedef unspecified_type Site_2;

/*!
A type for the triangulated Delaunay
graph. The type `Delaunay_graph` must be a model of the
`DelaunayGraph_2` concept.
*/
typedef unspecified_type Delaunay_graph;

/*!

*/
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle;

/*!

*/
typedef Delaunay_graph::Face_handle Delaunay_face_handle;

/*!

*/
typedef Delaunay_graph::Edge Delaunay_edge;

/*!

*/
typedef Delaunay_graph::All_edges_iterator All_Delaunay_edges_iterator;

/*!

*/
typedef Delaunay_graph::Finite_edges_iterator Finite_Delaunay_edges_iterator;

/*!

*/
typedef Delaunay_graph::Edge_circulator Delaunay_edge_circulator;

/*!
A type for the predicate functor that is
responsible for rejecting an edge of the Delaunay graph (or
equivalently rejecting its dual edge in the Voronoi diagram). It must be
model of the concepts `DefaultConstructible`,
`CopyConstructible`, `Assignable`, and `AdaptableBinaryFunction`.
It must provide the following operators:

`bool operator()(Delaunay_graph dg, Delaunay_edge e)`

`bool operator()(Delaunay_graph dg, Delaunay_face_handle f, int i)`

`bool operator()(Delaunay_graph dg, Delaunay_edge_circulator ec)`

`bool operator()(Delaunay_graph dg, All_Delaunay_edges_iterator eit)`

`bool operator()(Delaunay_graph dg, Finite_Delaunay_edges_iterator eit)`

The functor returns `true` iff the edge is rejected.
*/
typedef unspecified_type Edge_rejector;

/*!
A type for the predicate functor that is
responsible for rejecting a vertex of the Delaunay graph (or
equivalently its dual face in the Voronoi diagram - hence the name
of the functor). It must be model of the concepts `DefaultConstructible`,
`CopyConstructible`, `Assignable`, `AdaptableBinaryFunction`.
It must provide the following operator:

<CENTER>`bool operator()(Delaunay graph dg, Delaunay_vertex_handle v)`</CENTER>

The functor returns `true` iff the face is rejected.
*/
typedef unspecified_type Face_rejector;

/*!
A tag for determining if the adaptation
policy class provides a functor for inserting sites in the Delaunay
graph. This tag is equal to either `CGAL::Tag_true` (a site
inserter functor is available) or `CGAL::Tag_false` (a site
inserter functor is not available).
*/
typedef unspecified_type Has_inserter;

/*!
A type for a functor that inserts sites
in the Delaunay graph. It must be model of the concepts
`DefaultConstructible`, `CopyConstructible`, `Assignable`,
`AdaptableBinaryFunction`
following operator

<CENTER>`Delaunay_vertex_handle operator()(Delaunay_graph& dg, Site_2 t)`</CENTER>

The vertex handle returned either points to the vertex of the
Delaunay graph corresponding to the site just inserted or is the
default constructed vertex handle. The latter case can happen if the
site inserted is <I>hidden</I>, i.e., it has an empty Voronoi cell.

This type is required only if the `Has_inserter` tag is equal to
`CGAL::Tag_true`.
*/
typedef unspecified_type Site_inserter;

/// @}

/// \name Access to objects
/// @{

/*!

*/
Edge_rejector edge_rejector_object();

/*!

*/
Face_rejector face_rejector_object();

/*!
This method is
required only if `Has_inserter` is equal to `CGAL::Tag_true`.
*/
Site_inserter site_inserter_object();

/// @}

/// \name Miscellaneous
/// The following methods are important when the adaptation policy
/// maintains a state. This can happen if we have a caching adaptation
/// policy, i.e., when we cache the results of the edge and face
/// rejectors.
/// @{

/*!
Clears the state of the adaptation policy.
*/
void clear();

/*!
The adaptation policies
`ap` and `other` are swapped. This method should be
preferred to `ap=other` or `ap(other)` if
`other` is deleted afterwards.
*/
void swap(AdaptationPolicy_2 other);

/*!
Tests the validity of the adaptation policy.
*/
bool is_valid();

/*!
Tests the validity of the
adaptation policy using extra information from the Delaunay graph
`dg`.
*/
bool is_valid(Delaunay_graph dg);

/// @}

}; /* end AdaptationPolicy_2 */

