
namespace CGAL {

/*!
\ingroup PkgApolloniusGraph2Ref

The class `Apollonius_graph_2` represents the Apollonius graph.
It supports insertions and deletions of sites.

\tparam Gt is the geometric traits class and must be a model of `ApolloniusGraphTraits_2`.

\tparam Agds is the Apollonius graph data structure and must be a model of `ApolloniusGraphDataStructure_2`
whose vertex and face must be models of `ApolloniusGraphVertexBase_2` and `TriangulationFaceBase_2`,
respectively.
It defaults to:
\code
  CGAL::Triangulation_data_structure_2<
    CGAL::Apollonius_graph_vertex_base_2<Gt,true>,
    CGAL::Triangulation_face_base_2<Gt> >`
\endcode

\cgalHeading{Traversal of the Apollonius Graph}

An Apollonius graph can be seen as a container of faces and vertices.
Therefore the Apollonius graph provides several iterators and
circulators that allow to traverse it (completely or partially).

\cgalHeading{Traversal of the Convex Hull}

Applied on the `infinite_vertex` the `incident_*` functions allow to
visit the vertices on the convex hull and the infinite edges and
faces. Note that a counterclockwise traversal of the vertices adjacent
to the `infinite_vertex` is a clockwise traversal of the convex hull.

\code{.cpp}
CGAL::Apollonius_graph_2<Gt, Agds> ag;
CGAL::Apollonius_graph_2<Gt, Agds>::Face f;

ag.incident_vertices(ag.infinite_vertex());
ag.incident_vertices(ag.infinite_vertex(), f);

ag.incident_faces(ag.infinite_vertex());
ag.incident_faces(ag.infinite_vertex(), f);

ag.incident_edges(ag.infinite_vertex());
ag.incident_edges(ag.infinite_vertex(), f);
\endcode

\cgalModels `DelaunayGraph_2`

\sa `CGAL::Apollonius_graph_traits_2<K,Method_tag>`
\sa `CGAL::Apollonius_graph_filtered_traits_2<CK,CM,EK,EM,FK,FM>`
\sa `CGAL::Apollonius_graph_hierarchy_2<Gt,Agds>`
*/
template< typename Gt, typename Agds >
class Apollonius_graph_2 {
public:

/// \name Types
/// @{

/*!
A type for the underlying
data structure.
*/
typedef Agds Data_structure;

/*!
Same as the `Data_structure` type. This type has been introduced
in order for the `Apollonius_graph_2` class to be a
model of the `DelaunayGraph_2` concept.
*/
typedef Data_structure Triangulation_data_structure;

/*!
A type for the geometric traits.
*/
typedef Gt Geom_traits;

/*!
A type for the
point defined in the geometric traits.
*/
typedef Gt::Point_2 Point_2;

/*!
A type for the Apollonius site, defined in the geometric traits.
*/
typedef Gt::Site_2 Site_2;

/// @}

/// \name Handles And Iterators
/// The vertices and faces of the Apollonius graph are accessed
/// through `handles`, `iterators`, and `circulators`. The iterators
/// and circulators are all bidirectional and non-mutable. The
/// circulators and iterators are assignable to the corresponding
/// handle types, and they are also convertible to the corresponding
/// handles. The edges of the Apollonius graph can also be visited
/// through iterators and circulators, the edge circulators and
/// iterators are also bidirectional and non-mutable. In the
/// following, we call <I>infinite</I> any face or edge incident to
/// the infinite vertex and the infinite vertex itself. Any other
/// feature (face, edge or vertex) of the Apollonius graph is said to
/// be <I>finite</I>. Some iterators (the `All` iterators ) allow to
/// visit finite or infinite features while the others (the `Finite`
/// iterators) visit only finite features. Circulators visit both
/// infinite and finite features.
/// @{

/*!
the edge type.
The `Edge(f,i)` is the edge common to faces `f` and
`f.neighbor(i)`. It is also the edge joining the vertices
`vertex(cw(i))` and `vertex(ccw(i))` of `f`.
\pre `i` must be `0`, `1` or `2`.
*/
typedef Data_structure::Edge Edge;

/*!
A type for a vertex.
*/
typedef Data_structure::Vertex Vertex;

/*!
A type for a face.
*/
typedef Data_structure::Face Face;

/*!
A type for a handle to a vertex.
*/
typedef Data_structure::Vertex_handle Vertex_handle;

/*!
A type for a handle to a face.
*/
typedef Data_structure::Face_handle Face_handle;

/*!
A type for a circulator over vertices incident to a given vertex.
*/
typedef Data_structure::Vertex_circulator Vertex_circulator;

/*!
A type for a circulator over faces incident to a given vertex.
*/
typedef Data_structure::Face_circulator Face_circulator;

/*!
A type for a circulator over edges incident to a given vertex.
*/
typedef Data_structure::Edge_circulator Edge_circulator;

/*!
A type for an iterator over all vertices.
*/
typedef Data_structure::Vertex_iterator
All_vertices_iterator;

/*!
A type for an iterator over all faces.
*/
typedef Data_structure::Face_iterator
All_faces_iterator;

/*!
A type for an iterator over all edges.
*/
typedef Data_structure::Edge_iterator
All_edges_iterator;

/*!
An unsigned integral type.
*/
typedef Data_structure::size_type size_type;

/*!
A type for an iterator over finite vertices.
*/
typedef unspecified_type Finite_vertices_iterator;

/*!
A type for an iterator over finite faces.
*/
typedef unspecified_type Finite_faces_iterator;

/*!
A type for an iterator over finite edges.
*/
typedef unspecified_type Finite_edges_iterator;

/// @}

/// \name Site Iterators
/// In addition to iterators and circulators for vertices and faces,
/// iterators for sites are provided. In particular there are
/// iterators for the entire set of sites, the hidden sites and the
/// visible sites of the Apollonius graph.
/// @{

/*!
A type for an iterator over all sites.
*/
typedef unspecified_type Sites_iterator;

/*!
A type for an iterator over all visible sites.
*/
typedef unspecified_type Visible_sites_iterator;

/*!
A type for an iterator over all hidden sites.
*/
typedef unspecified_type Hidden_sites_iterator;

/// @}

/// \name Creation
/// @{

/*!
Creates an
Apollonius graph `ag` using `gt` as geometric traits.
*/
Apollonius_graph_2(Gt gt=Gt());

/*!
Creates an Apollonius graph `ag` using `gt` as
geometric traits and inserts all sites in the range
[`first`, `beyond`).
\pre `Input_iterator` must be a model of `InputIterator`. The value type of `Input_iterator` must be `Site_2`.
*/
template< class Input_iterator >
Apollonius_graph_2(Input_iterator first, Input_iterator beyond,
Gt gt=Gt());

/*!
Copy constructor. All faces and vertices are duplicated. After the
construction,
`ag` and `other` refer to two different Apollonius graphs : if
`other` is modified, `ag` is not.
*/
Apollonius_graph_2(const Apollonius_graph_2<Gt,Agds>& other);

/*!
Assignment. If `ag` and `other` are the same object
nothing is done. Otherwise, all the vertices and faces are
duplicated. After the assignment, `ag` and `other` refer to
different Apollonius graphs : if `other` is modified, `ag` is
not.
*/
Apollonius_graph_2<Gt,Agds>
operator=(const Apollonius_graph_2<Gt,Agds>& other);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a reference to the Apollonius graph traits object.
*/
const Geom_traits& geom_traits() const;

/*!
Returns a reference to the
underlying data structure.
*/
const Data_structure& data_structure() const;

/*!
Same as `data_structure()`. This
method has been added in compliance with the `DelaunayGraph_2`
concept.
*/
const Data_structure& tds() const;

/*!
Returns the dimension of the Apollonius graph.
*/
int dimension() const;

/*!
Returns the number of finite vertices.
*/
size_type number_of_vertices() const;

/*!
Returns the number of visible sites.
*/
size_type number_of_visible_sites() const;

/*!
Returns the number of hidden sites.
*/
size_type number_of_hidden_sites() const;

/*!
Returns the number of faces (both finite and infinite) of the
Apollonius graph.
*/
size_type number_of_faces() const;

/*!
Returns a face incident to the `infinite_vertex`.
*/
Face_handle infinite_face() const;

/*!
Returns the `infinite_vertex`.
*/
Vertex_handle infinite_vertex() const;

/*!
Returns a vertex distinct from the `infinite_vertex`.
\pre The number of (visible) vertices in the Apollonius graph must be at least one.
*/
Vertex_handle finite_vertex() const;

/// @}

/// \name Face, Edge and Vertex Iterators
/// The following iterators allow respectively to visit finite faces,
/// finite edges and finite vertices of the Apollonius graph. These
/// iterators are non-mutable, bidirectional and their value types are
/// respectively `Face`, `Edge` and `Vertex`. They are all invalidated
/// by any change in the Apollonius graph. The following iterators
/// allow respectively to visit all (both finite and infinite) faces,
/// edges and vertices of the Apollonius graph. These iterators are
/// non-mutable, bidirectional and their value types are respectively
/// `Face`, `Edge` and `Vertex`. They are all invalidated by any
/// change in the Apollonius graph.
/// @{

/*!
Starts at an arbitrary finite vertex.
*/
Finite_vertices_iterator finite_vertices_begin() const;

/*!
Past-the-end iterator.
*/
Finite_vertices_iterator finite_vertices_end() const;

/*!
Starts at an arbitrary finite edge.
*/
Finite_edges_iterator finite_edges_begin() const;

/*!
Past-the-end iterator.
*/
Finite_edges_iterator finite_edges_end() const;

/*!
Starts at an arbitrary finite face.
*/
Finite_faces_iterator finite_faces_begin() const;

/*!
Past-the-end iterator.
*/
Finite_faces_iterator finite_faces_end() const;

/*!
Starts at an arbitrary vertex.
*/
All_vertices_iterator all_vertices_begin() const;

/*!
Past-the-end iterator.
*/
All_vertices_iterator all_vertices_end() const;

/*!
Starts at an arbitrary edge.
*/
All_edges_iterator all_edges_begin() const;

/*!
Past-the-end iterator.
*/
All_edges_iterator all_edges_end() const;

/*!
Starts at an arbitrary face.
*/
All_faces_iterator all_faces_begin() const;

/*!
Past-the-end iterator.
*/
All_faces_iterator all_faces_end() const;

/// @}

/// \name Site Iterators
/// The following iterators allow respectively to visit all sites, the
/// visible sites and the hidden sites. These iterators are
/// non-mutable, bidirectional and their value type is `Site_2`. They
/// are all invalidated by any change in the Apollonius graph.
/// @{

/*!
Starts at an arbitrary site.
*/
Sites_iterator sites_begin() const;

/*!
Past-the-end iterator.
*/
Sites_iterator sites_end() const;

/*!
Starts at an arbitrary visible site.
*/
Visible_sites_iterator visible_sites_begin() const;

/*!
Past-the-end iterator.
*/
Visible_sites_iterator visible_sites_end() const;

/*!
Starts at an arbitrary hidden site.
*/
Hidden_sites_iterator hidden_sites_begin() const;

/*!
Past-the-end iterator.
*/
Hidden_sites_iterator hidden_sites_end() const;

/// @}

/// \name Face, Edge and Vertex Circulators
/// The Apollonius graph also provides circulators that allow to visit
/// respectively all faces or edges incident to a given vertex or all
/// vertices adjacent to a given vertex. These circulators are
/// non-mutable and bidirectional. The operator `operator++` moves the
/// circulator counterclockwise around the vertex while the
/// `operator--` moves clockwise. A face circulator is invalidated by
/// any modification of the face pointed to. An edge circulator is
/// invalidated by any modification in one of the two faces incident
/// to the edge pointed to. A vertex circulator is invalidated by any
/// modification in any of the faces adjacent to the vertex pointed
/// to.
/// @{

/*!
Starts at an arbitrary face incident
to `v`.
*/
Face_circulator incident_faces(Vertex_handle v) const;

/*!
Starts at face `f`.
\pre Face `f` is incident to vertex `v`.
*/
Face_circulator incident_faces(Vertex_handle v, Face_handle f) const;

/*!
Starts at an arbitrary edge incident
to `v`.
*/
Edge_circulator incident_edges(Vertex_handle v) const;

/*!
Starts at the first edge of `f` incident to
`v`, in counterclockwise order around `v`.
\pre Face `f` is incident to vertex `v`.
*/
Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const;

/*!
Starts at an arbitrary vertex incident
to `v`.
*/
Vertex_circulator incident_vertices(Vertex_handle v) const;

/*!
Starts at the first vertex of `f` adjacent to `v`
in counterclockwise order around `v`.
\pre Face `f` is incident to vertex `v`.
*/
Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f) const;

/// @}

/// \name Predicates
/// The class `Apollonius_graph_2` provides methods to test the finite
/// or infinite character of any feature.
/// @{

/*!
`true`, iff `v` is the `infinite_vertex`.
*/
bool
is_infinite(Vertex_handle v) const;

/*!
`true`, iff face `f` is infinite.
*/
bool
is_infinite(Face_handle f) const;

/*!
`true`, iff edge `(f,i)` is infinite.
*/
bool is_infinite(Face_handle f, int i) const;

/*!
`true`, iff edge `e` is infinite.
*/
bool
is_infinite(const Edge& e) const;

/*!
`true`, iff edge `*ec` is infinite.
*/
bool
is_infinite(Edge_circulator ec) const;

/// @}

/// \name Insertion
/// @{

/*!
Inserts the sites in the range
[`first`,`beyond`). The number of sites in the range
[`first`, `beyond`) is returned.
\pre `Input_iterator` must be a model of `InputIterator` and its value type must be `Site_2`.
*/
template< class Input_iterator >
unsigned int insert(Input_iterator first, Input_iterator beyond);

/*!
Inserts the
site `s` in the Apollonius graph. If `s` is visible then the
vertex handle of `s` is returned, otherwise
`Vertex_handle(nullptr)` is returned.
*/
Vertex_handle insert(const Site_2& s);

/*!
Inserts `s` in the Apollonius graph using the site
associated with `vnear` as an estimate for the nearest neighbor of
the center of `s`. If `s` is visible then the vertex handle of
`s` is returned, otherwise `Vertex_handle(nullptr)` is
returned.
*/
Vertex_handle insert(const Site_2& s, Vertex_handle vnear);

/// @}

/// \name Removal
/// @{

/*!
Removes the site
associated to the vertex handle `v` from the Apollonius
graph.
\pre `v` must correspond to a valid finite vertex of the Apollonius graph.
*/
void remove(Vertex_handle v);

/// @}

/// \name Nearest Neighbor Location
/// @{

/*!
Finds the nearest neighbor of the point `p`. In other words it
finds the site whose Apollonius cell contains `p`. Ties are broken
arbitrarily and one of the nearest neighbors of `p` is
returned. If there are no visible sites in the Apollonius diagram
`Vertex_handle(nullptr)` is returned.
*/
Vertex_handle nearest_neighbor(const Point_2& p) const;

/*!
Finds the nearest neighbor of the point
`p` using the site associated with `vnear` as an
estimate for the nearest neighbor of `p`. Ties are broken
arbitrarily and one of the nearest neighbors of `p` is
returned. If there are no visible sites in the Apollonius diagram
`Vertex_handle(nullptr)` is returned.
*/
Vertex_handle nearest_neighbor(const Point_2& p, Vertex_handle vnear) const;

/// @}

/// \name Access to the Dual
/// The `Apollonius_graph_2` class provides access to the duals of the
/// faces of the graph. The dual of a face of the Apollonius graph is
/// a site. If the originating face is infinite, its dual is a site
/// with center at infinity (or equivalently with infinite weight),
/// which means that it can be represented geometrically as a line. If
/// the originating face is finite, its dual is a site with finite
/// center and weight. In the following three methods the returned
/// object is assignable to either `Site_2` or `Gt::Line_2`, depending
/// on whether the corresponding face of the Apollonius graph is
/// finite or infinite, respectively.
/// @{

/*!
Returns the
dual corresponding to the face handle `f`. The returned object can
be assignable to one of the following: `Site_2`, `Gt::Line_2`.
*/
Gt::Object_2 dual(Face_handle f) const;

/*!
Returns the
dual of the face to which `it` points to. The returned object can
be assignable to one of the following: `Site_2`, `Gt::Line_2`.
*/
Gt::Object_2 dual(All_faces_iterator it) const;

/*!
Returns
the dual of the face to which `it` points to. The returned
object can be assignable to one of the following: `Site_2`,
`Gt::Line_2`.
*/
Gt::Object_2 dual(Finite_faces_iterator it) const;

/// @}

/// \name I/O
/// @{

/*!
Draws the Apollonius graph to
the stream `str`.
\pre The following operators must be defined:
`Stream& operator<<(Stream&, Gt::Segment_2)`,
`Stream& operator<<(Stream&, Gt::Ray_2)`.

*/
template< class Stream >
Stream& draw_primal(Stream& str) const;

/*!
Draws the dual of the
Apollonius graph, i.e., the Apollonius diagram, to the stream
`str`.
\pre The following operators must be defined:
`Stream& operator<<(Stream&, Gt::Segment_2)`,
`Stream& operator<<(Stream&, Gt::Ray_2)`,
`Stream& operator<<(Stream&, Gt::Line_2)`.

*/
template < class Stream >
Stream& draw_dual(Stream& str) const;

/*!
Draws the edge
`e` of the Apollonius graph to the stream `str`.
\pre The following operators must be defined:
`Stream& operator<<(Stream&, Gt::Segment_2)`,
`Stream& operator<<(Stream&, Gt::Ray_2)`.

*/
template< class Stream >
Stream& draw_primal_edge(const Edge& e, Stream& str) const;

/*!
Draws the dual of the
edge `e` to the stream `str`. The dual of `e` is an edge
of the Apollonius diagram.
\pre The following operators must be defined:
`Stream& operator<<(Stream&, Gt::Segment_2)`,
`Stream& operator<<(Stream&, Gt::Ray_2)`,
`Stream& operator<<(Stream&, Gt::Line_2)`.

*/
template< class Stream >
Stream& draw_dual_edge(const Edge& e, Stream& str) const;

/*!
Writes the current
state of the Apollonius graph to an output stream. In particular,
all visible and hidden sites are written as well as the
underlying combinatorial data structure.
*/
void file_output(std::ostream& os) const;

/*!
Reads the state of the
Apollonius graph from an input stream.
*/
void file_input(std::istream& is);

/*!
Writes the current state of the Apollonius graph to an output stream.
*/
std::ostream& operator<<(std::ostream& os, const Apollonius_graph_2<Gt,Agds>& ag) const;

/*!
Reads the state of the Apollonius graph from an input stream.
*/
std::istream& operator>>(std::istream& is, const Apollonius_graph_2<Gt,Agds>& ag);

/// @}

/// \name Validity Check
/// @{

/*!
Checks the validity of the Apollonius graph. If `verbose` is
`true` a short message is sent to `std::cerr`. If `level`
is 0, only the data structure is validated. If `level` is 1, then
both the data structure and the Apollonius graph are
validated. Negative values of `level` always return true, and
values greater than 1 are equivalent to `level` being 1.
*/
bool is_valid(bool verbose = false, int level = 1) const;

/// @}

/// \name Miscellaneous
/// @{

/*!
Clears all contents of the Apollonius graph.
*/
void clear();

/*!
The Apollonius graphs
`other` and `ag` are swapped. `ag.swap(other)` should
be preferred to `ag = other` or to `ag(other)` if
`other` is deleted afterwards.
*/
void swap(Apollonius_graph_2<Gt,Agds>& other);

/// @}

}; /* end Apollonius_graph_2 */
} /* end namespace CGAL */
