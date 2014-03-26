
/*!
\ingroup PkgVoronoiDiagramAdaptor2Concepts
\cgalConcept

The concept `DelaunayGraph_2` defines the requirements for the 
first template parameter of the `Voronoi_diagram_2<DG,AT,AP>` 
class. The `DelaunayGraph_2` concept essentially defines the 
requirements that a class representing a Delaunay graph must obey so 
that the Voronoi diagram adaptor can adapt it. 

\cgalRefines `DefaultConstructible,` \cgalRefines `CopyConstructible,` \cgalRefines `Assignable` 

\cgalHeading{Traversal of the Delaunay graph}

A model of the `DelaunayGraph_2` concept must provide several 
iterators and circulators that allow to traverse it (completely or 
partially). All iterators and circulators must be convertible to the 
corresponding handles. 

\cgalHasModel `CGAL::Delaunay_triangulation_2<Traits,Tds>` 
\cgalHasModel `CGAL::Regular_triangulation_2<Traits,Tds>` 
\cgalHasModel `CGAL::Triangulation_hierarchy_2<Tr>` provided that `Tr` is a model of `DelaunayGraph_2` 
\cgalHasModel `CGAL::Segment_Delaunay_graph_2<Gt,DS>` 
\cgalHasModel `CGAL::Segment_Delaunay_graph_hierarchy_2<Gt,STag,DS>` 
\cgalHasModel `CGAL::Apollonius_graph_2<Gt,Agds>` 
\cgalHasModel `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>` 

\sa `AdaptationTraits_2` 
\sa `AdaptationPolicy_2` 
\sa `CGAL::Voronoi_diagram_2<DG,AT,AP>` 

*/

class DelaunayGraph_2 {
public:

/// \name Types 
/// @{

/*!
A type for sizes. 
*/ 
typedef unspecified_type size_type; 

/*!
A type for the geometric traits associated 
with the Delaunay graph. 
*/ 
typedef unspecified_type Geom_traits; 

/*!
A type for the underlying 
triangulation data structure. It must be a model of the concept 
`TriangulationDataStructure_2.` 
*/ 
typedef unspecified_type Triangulation_data_structure; 

/*!
A type for the vertices of the Delaunay graph. 
*/ 
typedef unspecified_type Vertex; 

/*!
A type for the faces of the Delaunay graph. 
*/ 
typedef unspecified_type Face; 

/*!
The type of the 
edges of the Delaunay graph. 
*/ 
typedef std::pair<Face_handle,int> Edge; 

/*!
Handle to the vertices of the Delaunay graph. 
*/ 
typedef unspecified_type Vertex_handle; 

/*!
Handle to the faces of the Delaunay graph. 
*/ 
typedef unspecified_type Face_handle; 

/// @}

/// \name 
/// The following iterators and circulators must be defined. All
/// iterators and circulators must be assignable and convertible to
/// their corresponding handles.
/// @{

/*!
A type for an iterator over all edges of 
the Delaunay graph. Its value type must be `Edge`. 
*/ 
typedef unspecified_type All_edges_iterator; 

/*!
A type for an iterator over the finite 
edges of the Delaunay graph. Its value type must be `Edge`. 
*/ 
typedef unspecified_type Finite_edges_iterator; 

/*!
A type for an iterator over all faces of 
the Delaunay graph. Its value type must be `Face`. 
*/ 
typedef unspecified_type All_faces_iterator; 

/*!
A type for an iterator over the finite 
faces of the Delaunay graph. Its value type must be `Face`. 
*/ 
typedef unspecified_type Finite_faces_iterator; 

/*!
A type for an iterator over all 
vertices of the Delaunay graph. Its value type must be `Vertex`. 
*/ 
typedef unspecified_type All_vertices_iterator; 

/*!
A type for an iterator over 
the finite vertices of the Delaunay graph. Its value type must be 
`Vertex`. 
*/ 
typedef unspecified_type Finite_vertices_iterator; 

/*!
A type for a circulator over the 
adjacent faces of a vertex of the Delaunay graph. Its value type 
must be `Face`. 
*/ 
typedef unspecified_type Face_circulator; 

/*!
A type for a circulator over the 
adjacent vertices of a vertex of the Delaunay graph. Its value type 
must be `Vertex`. 
*/ 
typedef unspecified_type Vertex_circulator; 

/*!
A type for a circulator over the 
adjacent edges of a vertex of the Delaunay graph. Its value type 
must be `Edge`. 
*/ 
typedef unspecified_type Edge_circulator; 

/// @} 

/// \name Creation 
/// In addition to the default and copy constructors, as well as the
/// assignment operator, the following constructors are required.
/// @{

/*!
Constructor that takes an instance of the geometric traits. 
*/ 
DelaunayGraph_2(Geom_traits gt); 

/*!
Constructor that takes an iterator range. The value type of the 
iterator must be the type of the sites of the Delaunay graph. 
*/ 
template<class It> DelaunayGraph_2(It first, It beyond); 

/*!
Constructor that takes an iterator range and an instance of the 
geometric traits. The value type of the iterator must be the type of 
the sites of the Delaunay graph. 
*/ 
template<class It> DelaunayGraph_2(It first, It beyond, 
Geom_traits gt); 

/// @} 

/// \name Access methods 
/// @{

/*!
Returns a reference to the underlying triangulation data structure. 
*/ 
Triangulation_data_structure tds(); 

/*!
Returns a reference to the geometric traits object. 
*/ 
Geom_traits geom_traits(); 

/*!
Returns a handle to the infinite vertex. 
*/ 
Vertex_handle infinite_vertex(); 

/*!
Returns a handle to a finite vertex, provided there exists one. 
*/ 
Vertex_handle finite_vertex(); 

/*!
Returns a handle to a face incident to the infinite vertex. 
*/ 
Face_handle infinite_face(); 

/*!
Returns the dimension of the Delaunay graph. 
*/ 
int dimension(); 

/*!
Returns the number of finite vertices. 
*/ 
size_type number_of_vertices(); 

/*!
Returns the number of faces (both finite and infinite). 
*/ 
size_type number_of_faces(); 

/// @} 

/// \name Finite Face, Edge and Vertex Iterators 
/// The following iterators must allow, respectively, to visit finite
/// faces, finite edges and finite vertices of the Delaunay
/// graph. These iterators must be non-mutable, bidirectional and
/// their value types are respectively `Face`, `Edge` and `Vertex`.
/// @{

/*!
Starts at an arbitrary finite vertex. 
*/ 
Finite_vertices_iterator finite_vertices_begin(); 

/*!
Past-the-end iterator. 
*/ 
Finite_vertices_iterator finite_vertices_end(); 

/*!
Starts at an arbitrary finite edge. 
*/ 
Finite_edges_iterator finite_edges_begin(); 

/*!
Past-the-end iterator. 
*/ 
Finite_edges_iterator finite_edges_end(); 

/*!
Starts at an arbitrary finite face. 
*/ 
Finite_faces_iterator finite_faces_begin(); 

/*!
Past-the-end iterator. 
*/ 
Finite_faces_iterator finite_faces_end() 
const; 

/// @}

/// \name All Face, Edge and Vertex Iterators 
/// The following iterators must allow, respectively, to visit all
/// (both finite and infinite) faces, edges and vertices of the
/// Delaunay graph. These iterators are non-mutable, bidirectional and
/// their value types are respectively `Face`, `Edge` and `Vertex`.
/// @{

/*!
Starts at an arbitrary vertex. 
*/ 
All_vertices_iterator all_vertices_begin(); 

/*!
Past-the-end iterator. 
*/ 
All_vertices_iterator all_vertices_end(); 

/*!
Starts at an arbitrary edge. 
*/ 
All_edges_iterator all_edges_begin(); 

/*!
Past-the-end iterator. 
*/ 
All_edges_iterator all_edges_end(); 

/*!
Starts at an arbitrary face. 
*/ 
All_faces_iterator all_faces_begin(); 

/*!
Past-the-end iterator. 
*/ 
All_faces_iterator all_faces_end(); 

/// @} 

/// \name Face, Edge and Vertex Circulators 
/// A model of the `DelaunayGraph_2` concept must also provide
/// circulators that allow to visit, respectively, all faces or edges
/// incident to a given vertex or all vertices adjacent to a given
/// vertex. These circulators are non-mutable and bidirectional. The
/// operator `operator++` must move the circulator counterclockwise
/// around the vertex while the `operator-` must move the circulator
/// clockwise.
/// @{

/*!
Starts at an arbitrary face incident to `v`. 
*/ 
Face_circulator incident_faces(Vertex_handle v); 

/*!
Starts at face `f`. 
\pre Face `f` must be incident to vertex `v`. 
*/ 
Face_circulator incident_faces(Vertex_handle v, Face_handle f); 

/*!
Starts at an arbitrary edge incident to `v`. 
*/ 
Edge_circulator incident_edges(Vertex_handle v); 

/*!
Starts at the first edge of `f` incident to `v`, in 
counterclockwise order around `v`. 
\pre Face `f` must be incident to vertex `v`. 
*/ 
Edge_circulator incident_edges(Vertex_handle v, Face_handle f); 

/*!
Starts at an arbitrary vertex incident to `v`. 
*/ 
Vertex_circulator incident_vertices(Vertex_handle v); 

/*!
Starts at the first vertex of `f` adjacent to `v` in 
counterclockwise order around `v`. 
\pre Face `f` must be incident to vertex `v`. 
*/ 
Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f); 

/// @} 

/// \name Predicates 
/// A model of the `DelaunayGraph_2` concept must provide methods to
/// test the finite or infinite character of any feature.
/// @{

/*!
`true`, iff `v` is the `infinite_vertex`. 
*/ 
bool is_infinite(Vertex_handle v); 

/*!
`true`, iff face `f` is infinite. 
*/ 
bool is_infinite(Face_handle f); 

/*!
`true`, iff edge `(f,i)` is infinite. 
*/ 
bool is_infinite(Face_handle f, int i); 

/*!
`true`, iff edge `e` is infinite. 
*/ 
bool is_infinite(Edge e); 

/*!
`true`, iff edge `*ec` is infinite. 
*/ 
bool is_infinite(Edge_circulator ec); 

/// @} 

/// \name Validity check 
/// @{

/*!
Checks the validity of the Delaunay graph. If `verbose` is 
`true` a short message is sent to `std::cerr`. 
*/ 
bool is_valid(bool verbose = false); 

/// @} 

/// \name Miscellaneous 
/// @{

/*!
Clears all contents of the Delaunay graph. 
*/ 
void clear(); 

/*!
The Delaunay graphs 
`other` and `dg` are swapped. `dg`.`swap(other)` 
should be preferred to `dg``= other` or to 
`dg``(other)` if `other` is deleted afterwards. 
*/ 
void swap(DelaunayGraph_2 other); 

/// @}

}; /* end DelaunayGraph_2 */

