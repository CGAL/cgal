
namespace CGAL {

/*!
\ingroup PkgVoronoiDiagramAdaptor2

The class `Voronoi_diagram_2` provides an adaptor that enables us 
to view a triangulated Delaunay graph as their dual subdivision, the 
Voronoi diagram. The class `Voronoi_diagram_2` is designed to provide an API 
that is similar to that of \cgal's arrangements. 

The first template parameter of the `Voronoi_diagram_2` class corresponds to the 
triangulated Delaunay graph and must be a model of the 
`DelaunayGraph_2` concept. The second template parameter must be a 
model of the `AdaptationTraits_2` concept. The third template 
parameter must be a model of the `AdaptationPolicy_2` concept. The 
third template parameter defaults to `CGAL::Identity_policy_2<DG,AT>`. 

\refines ::DefaultConstructible, \refines ::CopyConstructible, \refines ::Assignable 

Traversal of the Voronoi diagram 
-------------- 

A Voronoi diagram can be seen as a container of faces, vertices and 
halfedges. Therefore the Voronoi diagram provides several iterators 
and circulators that allow to traverse it. 

\sa `DelaunayGraph_2` 
\sa `AdaptationTraits_2` 
\sa `AdaptationPolicy_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Face` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Halfedge` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>::Vertex` 
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>` 
\sa `CGAL::Regular_triangulation_2<Traits,Tds>` 
\sa `CGAL::Triangulation_hierarchy_2<Tr>` provided that `Tr` is a  model of `DelaunayGraph_2` 
\sa `CGAL::Segment_Delaunay_graph_2<Gt,DS>` 
\sa `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>` 
\sa `CGAL::Apollonius_graph_2<Gt,Agds>` 
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>` 
\sa `CGAL::Apollonius_graph_adaptation_traits_2<AG2>` 
\sa `CGAL::Delaunay_triangulation_adaptation_traits_2<DT2>` 
\sa `CGAL::Regular_triangulation_adaptation_traits_2<RT2>` 
\sa `CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2>` 
\sa `CGAL::Identity_policy_2<DG,AT>` 
\sa `CGAL::Apollonius_graph_degeneracy_removal_policy_2<AG2>` 
\sa `CGAL::Apollonius_graph_caching_degeneracy_removal_policy_2<AG2>` 
\sa `CGAL::Delaunay_triangulation_degeneracy_removal_policy_2<DT2>` 
\sa `CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT2>` 
\sa `CGAL::Regular_triangulation_degeneracy_removal_policy_2<RT2>` 
\sa `CGAL::Regular_triangulation_caching_degeneracy_removal_policy_2<RT2>` 
\sa `CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2>` 
\sa `CGAL::Segment_Delaunay_graph_caching_degeneracy_removal_policy_2<SDG2>` 

*/
template< typename DG, typename AT, typename AP >
class Voronoi_diagram_2 {
public:

/// \name Types 
/// @{

/*! 
A type for the dual Delaunay graph. 
*/ 
typedef DG Delaunay_graph; 

/*! 
A type for the adaptation 
traits needed by the Voronoi diagram adaptor. 
*/ 
typedef AT Adaptation_traits; 

/*! 
A type for the adaptation 
policy used. 
*/ 
typedef AP Adaptation_policy; 

/*! 
A type a point. 
*/ 
typedef Adaptation_traits::Point_2 Point_2; 

/*! 
A type 
for the sites of the Voronoi diagram. 
*/ 
typedef Adaptation_traits::Site_2 Site_2; 

/*! 
A type for sizes. 
*/ 
typedef Delaunay_graph::size_type size_type; 

/*! 
A type for the geometric traits of the Delaunay graph. 
*/ 
typedef Delaunay_graph::Geom_traits Delaunay_geom_traits; 

/*! 
A type for the vertex handles of the Delaunay graph. 
*/ 
typedef Delaunay_graph::Vertex_handle Delaunay_vertex_handle; 

/*! 
A type for the face handles of the Delaunay graph. 
*/ 
typedef Delaunay_graph::Face_handle Delaunay_face_handle; 

/*! 
A type for the edges of the Delaunay graph. 
*/ 
typedef Delaunay_graph::Edge Delaunay_edge; 

/*! 
A type for the halfedges of the Voronoi diagram. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
A type for the vertices of the Voronoi diagram. 
*/ 
typedef Hidden_type Vertex; 

/*! 
A type for the faces of the Voronoi diagram. 
*/ 
typedef Hidden_type Face; 

/// @}

/// \name
/// The vertices, edges and faces of the Voronoi diagram are accessed
/// through `handles`, `iterators` and `circulators`. The iterators
/// and circulators are all bidirectional and non-mutable. The
/// circulators and iterators are assignable to the corresponding
/// handle types, and they are also convertible to the corresponding
/// handles.
/// @{

/*! 
Handle for halfedges. 
*/ 
typedef Hidden_type Halfedge_handle; 

/*! 
Handle for vertices. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
Handle for faces. 
*/ 
typedef Hidden_type Face_handle; 

/*! 
A type for an iterator over Voronoi 
edges. Edges are considered non-oriented. Its value type is 
`Halfedge`. 
*/ 
typedef Hidden_type Edge_iterator; 

/*! 
A type for an iterator over Voronoi 
halfedges. Halfedges are oriented and come in pairs. Its value type 
is `Halfedge`. 
*/ 
typedef Hidden_type Halfedge_iterator; 

/*! 
A type for an iterator over Voronoi 
faces. Its value type is `Face`. 
*/ 
typedef Hidden_type Face_iterator; 

/*! 
A type for an iterator over Voronoi 
vertices. Its value type is `Vertex`. 
*/ 
typedef Hidden_type Vertex_iterator; 

/*! 
A type for a 
circulator over the halfedges that have a common vertex as their 
target. Its value type is `Halfedge`. 
*/ 
typedef Hidden_type Halfedge_around_vertex_circulator; 

/*! 
A type for a circulator over 
the halfedges on the boundary of a Voronoi face. Its value type of 
is `Halfedge`. 
*/ 
typedef Hidden_type Ccb_halfedge_circulator; 

/*! 
A type for an iterator over 
the unbounded faces of the Voronoi diagram. Its value type is 
`Face`. 
*/ 
typedef Hidden_type Unbounded_faces_iterator; 

/*! 
A type for an iterator over 
the bounded faces of the Voronoi diagram. Its value type is 
`Face`. 
*/ 
typedef Hidden_type Bounded_faces_iterator; 

/*! 
A type for an iterator over 
the unbounded halfedges of the Voronoi diagram. Its value type is 
`Halfedge`. 
*/ 
typedef Hidden_type Unbounded_halfedges_iterator; 

/*! 
A type for an iterator over 
the bounded halfedges of the Voronoi diagram. Its value type is 
`Halfedge`. 
*/ 
typedef Hidden_type Bounded_halfedges_iterator; 

/*! 
A type for an iterator over the 
sites of the Voronoi diagram. Its value type is `Site_2`. 
*/ 
typedef Hidden_type Site_iterator; 

/*! 
The result type of the point location queries. 
*/ 
typedef boost::variant<Face_handle,Halfedge_handle,Vertex_handle> 
Locate_result; 

/// @} 

/// \name Creation 
/// @{

/*! 
Creates a Voronoi diagram using `at` as adaptation traits and 
`ap` as adaptation policy; the underlying Delaunay graph is 
created using `gt` as geometric traits. 
*/ 
Voronoi_diagram_2(Adaptation_traits 
at = Adaptation_traits(), Adaptation_policy ap = Adaptation_policy(), 
Delaunay_geom_traits gt = Delaunay_geom_traits()); 

/*! 
Creates a Voronoi diagram from the Delaunay graph `dg` and using 
`at` as adaptation traits and `ap` as adaptation policy. The 
Delaunay graph `dg` is fully copied if `swap_dg` is set to 
`false`, or swapped with the one stored internally if 
`swap_dg` is set to `true`. 
*/ 
Voronoi_diagram_2(Delaunay_graph dg, bool swap_dg = false, 
Adaptation_traits at = Adaptation_traits(), Adaptation_policy ap = 
Adaptation_policy()); 

/*! 
Creates a Voronoi diagram using as sites the sites in the iterator 
range `[first, beyond)`, `at` as adaptation traits and 
`ap` as adaptation policy; the underlying Delaunay graph is 
created using `gt` as geometric traits. `Iterator` must be a 
model of the `InputIterator` concept and its value type must be 
`Site_2`. 
*/ 
template<class Iterator> 
Voronoi_diagram_2(Iterator first, Iterator beyond, Adaptation_traits 
at = Adaptation_traits(), Adaptation_policy ap = Adaptation_policy(), 
Delaunay_geom_traits gt = Delaunay_geom_traits()); 

/// @} 

/// \name Access Methods 
/// @{

/*! 
Returns a const reference to the dual graph, i.e., the Delaunay graph. 
*/ 
Delaunay_graph dual(); 

/*! 
Returns a handle to the halfedge in the Voronoi diagram that is dual 
to the edge `e` in the Delaunay graph. 
*/ 
Halfedge_handle dual(Delaunay_edge e); 

/*! 
Returns a handle to the face in the Voronoi diagram that is dual to 
the vertex corresponding to the vertex handle `v` in the 
Delaunay graph. 
*/ 
Face_handle dual(Delaunay_vertex_handle v); 

/*! 
Returns a handle to the vertex in the Voronoi diagram that is dual to 
the face corresponding to the face handle `f` in the Delaunay graph. 
*/ 
Vertex_handle dual(Delaunay_face_handle f); 

/*! 
Returns a reference to the Voronoi traits. 
*/ 
Adaptation_traits adaptation_traits(); 

/*! 
Returns a reference to the adaptation policy. 
*/ 
Adaptation_policy adaptation_policy(); 

/*! 
Returns the number of Voronoi vertices. 
*/ 
size_type number_of_vertices(); 

/*! 
Returns the number of Voronoi faces (bounded and unbounded). 
*/ 
size_type number_of_faces(); 

/*! 
Returns the number of halfedges (bounded and unbounded) in the 
Voronoi diagram. This is always an even number. 
*/ 
size_type number_of_halfedges(); 

/*! 
Returns the number of connected components of the Voronoi skeleton. 
*/ 
size_type number_of_connected_components(); 

/*! 
Returns one of the unbounded 
faces of the Voronoi diagram. If no unbounded faces exist (this can 
happen if the number of sites is zero) the 
default constructed face handle is returned. 
*/ 
Face_handle unbounded_face(); 

/*! 
Returns one of the bounded 
faces of the Voronoi diagram. If no bounded faces exist the default 
constructed face handle is returned. 
*/ 
Face_handle bounded_face(); 

/*! 
Returns one of the unbounded 
halfedges of the Voronoi diagram. If no unbounded halfedges exist the 
default constructed halfedge handle is returned. 
*/ 
Halfedge_handle unbounded_halfedge(); 

/*! 
Returns one of the bounded 
halfedges of the Voronoi diagram. If no bounded halfedges exist the 
default constructed halfedge handle is returned. 
*/ 
Halfedge_handle bounded_halfedge(); 

/// @} 

/// \name Iterators 
/// The following iterators allow respectively to visit the faces (all
/// or only the unbounded/bounded ones), edges, halfedges (all or only
/// the unbounded/bounded ones) and vertices of the Voronoi
/// diagram. These iterators are non-mutable, bidirectional and their
/// value types are respectively `Face`, `Halfedge`, `Halfedge` and
/// `Vertex`. All iterators are convertible to the corresponding
/// handles and are invalidated by any change in the Voronoi
/// diagram. The following iterator provides access to the sites that
/// define the Voronoi diagram. Its value type is `Site_2`. It is
/// invalidated by any change in the Voronoi diagram.
/// @{

/*! 
Starts at an arbitrary Voronoi face. 
*/ 
Face_iterator faces_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Face_iterator faces_end(); 

/*! 
Starts at an arbitrary unbounded Voronoi face. 
*/ 
Unbounded_faces_iterator unbounded_faces_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Unbounded_faces_iterator unbounded_faces_end(); 

/*! 
Starts at an arbitrary bounded Voronoi face. 
*/ 
Bounded_faces_iterator bounded_faces_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Bounded_faces_iterator bounded_faces_end(); 

/*! 
Starts at an arbitrary Voronoi edge. 
*/ 
Edge_iterator edges_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Edge_iterator edges_end(); 

/*! 
Starts at an arbitrary Voronoi halfedge. 
*/ 
Halfedge_iterator halfedges_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Halfedge_iterator halfedges_end(); 

/*! 
Starts at an arbitrary unbounded Voronoi edge. 
*/ 
Unbounded_halfedges_iterator unbounded_halfedges_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Unbounded_halfedges_iterator unbounded_halfedges_end(); 

/*! 
Starts at an arbitrary bounded Voronoi edge. 
*/ 
Bounded_halfedges_iterator bounded_halfedges_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Bounded_halfedges_iterator bounded_halfedges_end(); 

/*! 
Starts at an arbitrary Voronoi vertex. 
*/ 
Vertex_iterator vertices_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Vertex_iterator vertices_end(); 

/*! 
Starts at an arbitrary site. 
*/ 
Site_iterator sites_begin(); 

/*! 
Past-the-end iterator. 
*/ 
Site_iterator sites_end(); 

/// @} 

/// \name Circulators 
/// The Voronoi diagram adaptor also provides circulators that allow
/// to visit all halfedges whose target is a given vertex - this is
/// the `Halfedge_around_vertex_circulator`, as well as all halfedges
/// on the boundary of a Voronoi face - this is the
/// `Ccb_halfedge_circulator`. These circulators are non-mutable and
/// bidirectional. The operator `operator++` moves the former
/// circulator counterclockwise around the vertex while the
/// `operator-` moves clockwise. The latter circulator is moved by the
/// operator `operator++` to the next halfedge on the boundary in the
/// counterclockwise sense, while `operator-` moves clockwise. When
/// the `Ccb_halfedge_circulator` is defined over an infinite Voronoi
/// face `f`, then applying `operator++` to a circulator corresponding
/// to a halfedge whose target is not finite moves to the next
/// infinite (or semi-infinite) halfedge of `f` in the
/// counterclockwise sense. Similarly, applying `operator++` to a
/// circulator corresponding to a halfedge whose source is not finite,
/// moves to the previous infinite (or semi-infinite) halfedge of `f`
/// in the clockwise sense. The `Halfedge_around_vertex_circulator`
/// circulator is invalidated by any modification in the faces
/// adjacent to the vertex over which it is defined. The
/// `Ccb_halfedge_circulator` is invalidated by any modification in
/// the face over which it is defined.
/// @{

/*! 
Returns a circulator over the halfedges on the boundary of `f`. 
The circulator is initialized to an arbitrary halfedge on the 
boundary of the Voronoi face `f`. 
*/ 
Ccb_halfedge_circulator ccb_halfedges(Face_handle f); 

/*! 
Returns a circulator over the halfedges on the boundary of 
`f`. The circulator is initialized with the halfedge `h`. 
\pre The halfedge `h` must lie on the boundary of `f`. 
*/ 
Ccb_halfedge_circulator ccb_halfedges(Face_handle f, 
Halfedge_handle h); 

/*! 
Returns a circulator over the halfedges whose target is the Voronoi 
vertex `v`. The circulator is initialized to an arbitrary 
halfedge incident to `v`. 
*/ 
Halfedge_around_vertex_circulator 
incident_halfedges(Vertex_handle v); 

/*! 
Returns a circulator over the halfedges whose target is the Voronoi 
vertex `v`. The circulator is initialized with the halfedge 
`h`. 
\pre The vertex `v` must be the target vertex of the halfedge `h`. 
*/ 
Halfedge_around_vertex_circulator 
incident_halfedges(Vertex_handle v, Halfedge_handle h); 

/// @} 

/// \name Insertion 
/// @{

/*! 
Inserts the site 
`t` in the Voronoi diagram. A handle to the face corresponding 
to the Voronoi face of `t` in the Voronoi diagram is 
returned. If `t` has an empty Voronoi cell, the default 
constructed face handle is returned. This method is supported only 
if `Voronoi_traits::Has_inserter` is set to 
`CGAL::Tag_true`. 
*/ 
Face_handle insert(Site_2 t); 

/*! 
Inserts, in the 
Voronoi diagram, the sites in the iterator range `[first, beyond)`. The value type of `Iterator` must be 
`Site_2`. The number of sites in the iterator range is 
returned. This method is supported only if 
`Voronoi_traits::Has_inserter` is set to `CGAL::Tag_true`. 
*/ 
template<class Iterator> 
size_type insert(Iterator first, Iterator beyond); 

/// @} 

/// \name Queries 
/// @{

/*! 
Performs point location for 
the query point `p`. In other words, the face, halfedge or 
vertex of the Voronoi diagram is found on which the point `p` 
lies. This method is supported only if 
`Voronoi_traits::Has_nearest_site_2` is set to 
`CGAL::Tag_true`. 
\pre The Voronoi diagram must contain at least one face. 
*/ 
Locate_result locate(Point_2 p); 

/// @} 

/// \name I/O 
/// @{

/*! 
Writes the current state of the Voronoi diagram to the output 
stream `os`. 

The following operator must be defined: 

`std::ostream& operator<<(std::ostream&, Delaunay_graph)` 

*/ 
void file_output(std::ostream& os); 

/*! 
Reads the current state of the Voronoi diagram from the input 
stream `is`. 

The following operator must be defined: 

`std::istream& operator>>(std::istream&, Delaunay_graph)` 

*/ 
void file_input(std::istream& is); 

/*! 
Writes the current state of the Voronoi diagram to the output 
stream `os`. 

The following operator must be defined: 

`std::ostream& operator<<(std::ostream&, Delaunay_graph)` 
*/ 
std::ostream& operator<<(std::ostream& os, Voronoi_diagram_2<DG,AT,AP> vd); 

/*! 
Reads the current state of the Voronoi diagram from the input 
stream `is`. 

The following operator must be defined: 
`std::istream& operator>>(std::istream&, Delaunay_graph)` 
*/ 
std::istream& operator>>(std::istream& is, Voronoi_diagram_2<DG,AT,AP> vd); 

/// @} 

/// \name Validity check 
/// @{

/*! 
Checks the validity of the dual Delaunay 
graph and the Voronoi diagram adaptor. 
*/ 
bool is_valid(); 

/// @} 

/// \name Miscellaneous 
/// @{

/*! 
Clears all contents of the Voronoi diagram. 
*/ 
void clear(); 

/*! 
The Voronoi 
diagrams `other` and `vd` are 
swapped. `vd`.`swap(other)` should be preferred to 
`vd``= other` or to 
`vd``(other)` if `other` is deleted afterwards. 
*/ 
void swap(Voronoi_diagram_2<DG,AT,AP> other); 

/// @}

}; /* end Voronoi_diagram_2 */
} /* end namespace CGAL */
