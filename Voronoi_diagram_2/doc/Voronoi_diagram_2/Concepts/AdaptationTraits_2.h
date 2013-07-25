
/*!
\ingroup PkgVoronoiDiagramAdaptor2Concepts
\cgalConcept

The concept `AdaptationTraits_2` defines the functors required for 
accessing geometric information in the Delaunay graph that is needed by the 
`Voronoi_diagram_2<DG,AT,AP>` class. 
It optionally defines a functor for performing nearest site queries. A 
tag is provided for determining whether this functor is defined or not. 

\cgalRefines `DefaultConstructible,` \cgalRefines `CopyConstructible,` \cgalRefines `Assignable` 

\cgalHasModel `CGAL::Apollonius_graph_adaptation_traits_2<AG2>`
\cgalHasModel `CGAL::Delaunay_triangulation_adaptation_traits_2<DT2>`
\cgalHasModel `CGAL::Regular_triangulation_adaptation_traits_2<RT2>`
\cgalHasModel `CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2>`

\sa `DelaunayGraph_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 

*/

class AdaptationTraits_2 {
public:

/// \name Types 
/// @{

/*!
A type for a point. 
*/ 
typedef unspecified_type Point_2; 

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
The type of the edges of the Delaunay graph 
*/ 
typedef Delaunay_graph::Edge Delaunay_edge; 

/*!
The type of the face handles of the Delaunay graph 
*/ 
typedef Delaunay_graph::Face_handle Delaunay_face_handle; 

/*!
The type of the vertex handles of the Delaunay graph. 
*/ 
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle; 

/*!
A type for a functor that accesses the 
site associated with a vertex. The functor should be a model of the 
concepts `DefaultConstructible`, `CopyConstructible`, 
`Assignable` and `AdaptableFunctor` (with one argument). The 
functor must provide the following operator: 

<CENTER>`result_type operator()(Delaunay_vertex_handle v)`</CENTER> 

where the result type `result_type` must be either `Site_2` or 
`const Site_2&`. 
*/ 
typedef unspecified_type Access_site_2; 

/*!
A type for a functor that 
constructs the dual point of a (triangular) face in the Delaunay 
graph. This point is the Voronoi vertex of the three sites defining 
the face in the Delaunay graph. The functor must be a model of the 
concepts `DefaultConstructible`, `CopyConstructible`, 
`Assignable`, `AdaptableFunctor` (with one argument). It 
must provide the following operator: 

<CENTER>`Point_2 operator()(Delaunay_face_handle f)`</CENTER>. 
The face handle `f` must not correspond to an infinite face. 
*/ 
typedef unspecified_type Construct_Voronoi_point_2; 

/*!
A tag for determining if the adaptation 
traits class provides a functor for performing nearest site queries. 
This tag is equal to either `CGAL::Tag_true` (a nearest site 
query functor is available) or `CGAL::Tag_false` (a nearest site 
query functor is not available). 
*/ 
typedef unspecified_type Has_nearest_site_2; 

/*!
A type for a functor that performs 
nearest site queries. Semantically, the result of the query is 
either a face, edge or vertex of the Delaunay graph. It is a face if 
the query point has at least three closest sites; the returned face 
has closest sites as vertices. It is an edge if the query point is 
equidistant to exactly two vertices of the Delaunay graph, which are 
the source and target vertices of the edge. In all other cases, the 
search result is a vertex, namely, the unique vertex of the 
Delaunay graph closest to the query point. The functor must be a 
model of the concepts `DefaultConstructible`, 
`CopyConstructible`, `Assignable`, `AdaptableFunctor` 
(with two arguments). It must provide the following operator: 

<CENTER>`result_type operator()(Delaunay_graph dg, Point_2 p)`</CENTER> 

where the result type `result_type` is 
`boost::variant<Delaunay_vertex_handle,Delaunay_edge,Delaunay_face_handle>`. 

This type is required only if 
`Has_nearest_site_2` is equal to `CGAL::Tag_true`. 
*/ 
typedef unspecified_type Nearest_site_2; 

/// @} 

/// \name Access to objects 
/// @{

/*!

*/ 
Access_site_2 access_site_2_object(); 

/*!

*/ 
Construct_Voronoi_point_2 construct_Voronoi_point_2_object(); 

/*!
This method is 
required only if `Has_nearest_site_2` is equal to `CGAL::Tag_true`. 
*/ 
Nearest_site_2 nearest_site_2_object(); 

/// @}

}; /* end AdaptationTraits_2 */

