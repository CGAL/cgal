
/*!
\ingroup PkgTDS2Concepts
\cgalConcept

The concept `TriangulationDataStructure_2` describes the requirements for
the second template parameter of the basic triangulation class
`Triangulation_2<Traits,Tds>` and of all other 2D triangulation classes.

The concept can be seen as a container for the
faces and vertices of the triangulation.
The concept `TriangulationDataStructure_2` includes two sub-concepts
`TriangulationDataStructure_2::Vertex` and
`TriangulationDataStructure_2::Face`.

The `TriangulationDataStructure_2`
maintains incidence and adjacency relations
among vertices and faces.

Each triangular face gives access to its three incident vertices
and to its three adjacent faces.
Each vertex gives access to one of its incident faces
and through that face to the circular list of its incident faces.

The three vertices of a face are indexed with 0, 1, and 2.
The neighbors of a face are also
indexed with 0, 1, and 2 in such a way that the neighbor indexed by `i`
is opposite to the vertex with the same index.

Each edge has two implicit representations : the edge
of a face `f` which is opposed to the vertex indexed `i`,
can be represented as well as an edge of the `neighbor(i)` of
`f`. See \cgalFigureRef{TDS_2D_Fig_neighbors1}.

The triangulation data structure
is responsible for the combinatorial integrity of the triangulation.
This means that the triangulation data structure
allows to perform some combinatorial operations
on the triangulation and guarantees the maintenance on
proper incidence and adjacency relations among the vertices
and faces. The term combinatorial operations
means that those operations are purely topological
and do not depend on the geometric embedding.
Insertion of a new vertex in a given face, or in a given edge,
suppression of a vertex of degree three, flip of two edges
are examples of combinatorial operations.

\cgalHeading{I/O}

The information output in the `iostream` is:
the dimension, the number of (finite) vertices,
the number of (finite) faces.
Then comes
for each vertex, the non combinatorial information stored in that vertex
if any.
Then comes
for each faces, the indices of its vertices and
the non combinatorial information (if any) stored in this face.
Then comes
for each face again
the indices of the neighboring faces.
The index of an item (vertex of face)
the rank of this item in the output order.
When dimension \f$ <\f$ 2, the same information is output
for faces of maximal dimension instead of faces.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_data_structure_2<Vb,Fb>}
\cgalHasModelsEnd

\sa `TriangulationDataStructure_2::Face`
\sa `TriangulationDataStructure_2::Vertex`
\sa `CGAL::Triangulation_2`

*/

class TriangulationDataStructure_2 {
public:

/// \name Types
/// @{

/*!
Size type (unsigned integral type)
*/
typedef unspecified_type size_type;

/*!
Difference type (signed integral type)
*/
typedef unspecified_type difference_type;

/*!
The vertex type, requirements for this type are described in concept `TriangulationDataStructure_2::Vertex`.
*/
typedef unspecified_type Vertex;

/*!
The face type, requirements for this type are described in concept `TriangulationDataStructure_2::Face`.
*/
typedef unspecified_type Face;

/*!
%Face data type, requirements are described in `TriangulationDataStructure_2::Face_data`.
*/
typedef unspecified_type Face_data;

/*!
Handle to a vertex.
\cgalModels{Handle}
*/
typedef unspecified_type Vertex_handle;

/*!
Handle to a face.
\cgalModels{Handle}
*/
typedef unspecified_type Face_handle;

/*!
The edge type.
The `Edge(f,i)` is edge common to faces `f` and
`f.neighbor(i)`. It is also the edge joining the vertices
`vertex(cw(i))` and `vertex(ccw(i))` of `f`.
*/
typedef std::pair<Face_handle,int> Edge;

/// @}

/// \name Iterators and Circulators
///
/// The iterators allow one to visit all the vertices, edges
/// and faces of a triangulation data structure.
///
/// The circulators
/// allow to visit all the edges or faces incident to a given vertex
/// and all the vertices adjacent to a given vertex.
///
/// The iterators and circulators  are bidirectional and non mutable,  and they are
/// convertible to the corresponding handles, thus they can be passed
/// directly as argument to the functions expecting a handle.
///
/// A face circulator is invalidated by any modification of the face it
/// points to. An edge circulator is invalidated by any modification of
/// any of the two faces incident to the edge pointed to. A vertex
/// circulator that turns around vertex `v` and that has as value a handle
/// to vertex `w`, is invalidated by any modification of anyone of the two
/// faces incident to `v` and `w`.

/// @{


/*!

*/
typedef unspecified_type Face_iterator;

/*!

*/
typedef unspecified_type Edge_iterator;

/*!

*/
typedef unspecified_type Vertex_iterator;

/*!

*/
typedef unspecified_type Face_circulator;

/*!

*/
typedef unspecified_type Edge_circulator;

/*!

*/
typedef unspecified_type Vertex_circulator;

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
TriangulationDataStructure_2();

/*!
Copy constructor, performing a deep copy, that is all vertices and faces are duplicated.
*/
TriangulationDataStructure_2( const
TriangulationDataStructure_2& tds1);

/*!
Assignment. All the vertices and faces of `tds1` are duplicated
in the triangulation data structure. Former faces and vertices of
the triangulation data structure , if any, are
deleted.
*/
TriangulationDataStructure_2& operator=( const
TriangulationDataStructure_2& tds1);

/*!
`tds1` is copied into the triangulation data structure.
If `v != Vertex_handle()`, the vertex of the triangulation data structure
corresponding to `v` is returned, otherwise `Vertex_handle()`
is returned.
\pre The optional argument `v` is a vertex of `tds1`.
*/
Vertex_handle
copy_tds(const TriangulationDataStructure_2 & tds1,
Vertex_handle v = Vertex_handle());

/*!
`tds_src` is copied into `this`. As the vertex and face types might be different
and incompatible, the creation of new faces and vertices is made thanks to the
functors `convert_vertex` and `convert_face`, that convert vertex and face types.
For each vertex `v_src` in `tds_src`, the corresponding vertex `v_tgt` in `this` is a
copy of the vertex returned by `convert_vertex(v_src)`. The same operations are
done for faces with the functor convert_face. If `v != TDS_src::Vertex_handle()`,
a handle to the vertex created in `this` that is the copy of `v` is returned,
otherwise `Vertex_handle()` is returned.

 - A model of `ConvertVertex` must provide two operator()'s that are responsible for converting the source vertex `v_src` into the target vertex:
  - `Vertex operator()(const TDS_src::Vertex& v_src);` This operator is used to create the vertex from `v_src`.
  - `void operator()(const TDS_src::Vertex& v_src, Vertex& v_tgt);` This operator is meant to be used in case heavy data should transferred to `v_tgt`.
 - A model of ConvertFace must provide two operator()'s that are responsible for converting the source face `f_src` into the target face:
  - `Face operator()(const TDS_src::Face& f_src);` This operator is used to create the face from `f_src`.
  - `void operator()(const TDS_src::Face& f_src, Face& f_tgt);` This operator is meant to be used in case heavy data should transferred to `f_tgt`.

\pre The optional argument `v` is a vertex of `tds_src` or is `Vertex_handle()`.
*/
template <class TDS_src, class ConvertVertex, class ConvertFace>
Vertex_handle copy_tds(const TDS_src& tds_src, typename TDS_src::Vertex_handle v, const ConvertVertex& convert_vertex, const ConvertFace& convert_face);

/*!
Swaps the triangulation data structure and `tds1`.
Should be preferred to an assignment or copy constructor
when `tds1` is deleted after that.
*/
void swap( TriangulationDataStructure_2 & tds1);

/*!
Deletes all faces and all finite vertices.
*/
void clear();

/// @}

/// \name Access Functions
/// @{

/*!
  \cgalAdvancedFunction
  \cgalAdvancedBegin
  returns the dimension of the triangulation data structure.
  \cgalAdvancedEnd
*/
int dimension() const;

/*!
returns the number of vertices in the triangulation data structure.
*/
size_type number_of_vertices() const;

/*!
returns the number of two dimensional faces in the triangulation data structure.
*/
size_type number_of_faces() const ;

/*!
returns the number of edges in the triangulation data structure.
*/
size_type number_of_edges() const;

/*!
returns the number of full dimensional faces,
i.e.\ faces of dimension equal to the dimension
of the triangulation data structure. This is the actual
number of faces stored in the triangulation data structure.
*/
size_type number_of_full_dim_faces() const;

/// @}

/// \name Setting
/// @{

/*!
sets the dimension.
*/
void set_dimension (int n);

/// @}

/// \name Queries
/// @{

/*!
returns `true` if
`v` is a vertex of the triangulation data structure.
*/
bool is_vertex(Vertex_handle v) const;

/*!
returns `true` if `(fh,i)` is an edge of  the triangulation data structure.
Returns `false` when `dimension() < 1`.
*/
bool is_edge(Face_handle fh, int i) const;

/*!
returns `true` if
`(va, vb)` is an edge of  the triangulation data structure.
*/
bool is_edge(Vertex_handle va, Vertex_handle vb) const;

/*!
returns `true` if
`(va, vb)` is an edge of  the triangulation data structure.
In addition, if true is returned
`fr` and `i` are set such that the pair `(fr,i)`
is the description
of the ordered edge `(va, vb)`.
*/
bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle &fr,
int &i) const;

/*!
returns `true` if `fh` is a face of  the triangulation data structure.
Returns `false` when `dimension() < 2`.
*/
bool is_face(Face_handle fh) const;

/*!
`true` if there is a face having `v1`, `v2`, and `v3`
as vertices.
*/
bool is_face(Vertex_handle v1, Vertex_handle v2,
Vertex_handle v3)
const;

/*!
`true` if there is a face having `v1`, `v2`, and `v3`
as vertices.
 In addition, if `true` is returned, `fr` is a handle
to the face with `v1`, `v2` and `v3`
as vertices.
*/
bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
Face_handle &fr) const;

/// @}

/// \name Traversing the Triangulation

/// @{

/*!
visits all faces.
*/
Face_iterator faces_begin() const;

/*!

*/
Face_iterator faces_end() const;

/*!
visits all vertices.
*/
Vertex_iterator vertices_begin() const;

/*!

*/
Vertex_iterator vertices_end() const;

/*!
visits all edges.
*/
Edge_iterator edges_begin() const;

/*!

*/
Edge_iterator edges_end() const;

/*!
\pre If the face `f` is given, it has to be incident to be a face incident to `v` and the circulator begins with the vertex `f->vertex(ccw(i))` if `i` is the index of `v` in `f`.
*/
Vertex_circulator
incident_vertices(Vertex_handle v, Face_handle f=Face_handle())
const;

/*!
\pre If the face `f` is given, it has to be a face of incident to `v` and the circulator begins with the edge `(f,cw(i))` of `f` if `i` is the index of `v` in `f`.
*/
Edge_circulator
incident_edges(Vertex_handle v, Face_handle f=Face_handle()) const;

/*!
\pre If the face `f` is given, it has to be a face of incident to `v` and the circulator begins with the face `f`.
*/
Face_circulator
incident_faces(Vertex_handle v, Face_handle f=Face_handle()) const;

/*!
returns vertex of `f->neighbor(i)`.
*/
Vertex_handle mirror_vertex(Face_handle f, int i) const;

/*!
returns the index of `f` as a neighbor of `f->neighbor(i)`.
*/
int mirror_index(Face_handle f, int i) const;

/*!
returns the same edge seen from the other adjacent face.
*/
Edge mirror_edge(Edge e) const;

/// @}

/// \name Modifiers

/// @{

/*!
exchanges the edge incident to
`f` and `f->neighbor(i)` with the other
diagonal of the quadrilateral formed by `f` and `f->neighbor(i)`.

\image html Flip.png "Flip"
\image latex Flip.png "Flip"
*/
void flip(Face_handle f, int i);

/*!
creates the first
vertex and returns a handle to it.
*/
Vertex_handle insert_first();

/*!
creates the second
vertex and returns a handle to it.
*/
Vertex_handle insert_second();

/*!
adds a
vertex `v` splitting
edge `i` of face `f`. Return a handle to `v`.
*/
Vertex_handle insert_in_edge(Face_handle f, int i);

/*!
adds a vertex
`v` splitting face
`f` in three. %Face `f` is modified,
two new faces are created. Return a handle to `v`
*/
Vertex_handle insert_in_face(Face_handle f);

/*!
adds
a vertex `v`, increasing by one the dimension of the triangulation data structure.
%Vertex `v` and the existing vertex `w` are linked to all
the vertices of the triangulation data structure.
The Boolean `orient` decides the final orientation of all
faces. A handle to vertex `v` is returned.

*/
Vertex_handle insert_dim_up(Vertex_handle w, bool
orient=true);

/*!
removes a vertex of degree 3. Two of the incident faces are destroyed,
the third one is modified.
If parameter `f` is specified, it has to be a face incident to `v`
and will be the modified face.
\pre %Vertex `v` is a finite vertex with degree 3 and, if specified, face `f` is incident to `v`.

\image html Three.png "Insertion"
\image latex Three.png "Insertion"
*/
void remove_degree_3(Vertex_handle v,  Face_handle f = Face_handle());

/*!
removes the before last
vertex.
*/
void remove_second(Vertex_handle v);

/*!
removes the last vertex.
*/
void remove_first(Vertex_handle v);

/*!
removes vertex `v` incident to all other vertices
and decreases by one the dimension of the triangulation data structure.
\pre If the dimension is 2, the number of vertices is more than 3, if the dimension is 1, the number of vertices is 2.
*/
void remove_dim_down(Vertex_handle v);

/*!
must be called when the displacement of a vertex decreases the dimension of the triangulation data structure.

The link of a vertex `v` is formed by the edges disjoint from `v`
that are included in the faces incident to `v`.
When the link of `v = f->vertex(i)` contains all the other vertices
of the two-dimensional triangulation data structure (\f$ \mathbb{S}^2\f$), `dim_down()` crushes the two-dimensional
data-structure (\f$ \mathbb{S}^2\f$) onto the one-dimensional data structure (\f$ \mathbb{S}^1\f$) formed by the link of `v`
augmented with the vertex `v` itself; this one is placed on the edge `(f, i)`
(see Figure \ref figtdsdim_down_2).
\pre `dimension()` must be equal to `2`, the degree of `f->vertex(i)` must be equal to the total number of vertices minus 1.

\anchor figtdsdim_down_2
\image html tds-dim_down.png "From a two-dimensional data structure to a one-dimensional data structure."
\image latex tds-dim_down.png "From a two-dimensional data structure to a one-dimensional data structure."

*/
void dim_down(Face_handle f, int i);

/// @}

/// \name Advanced Modifiers
/// The following modifiers are required for convenience of the advanced user.
/// They do not guarantee the combinatorial validity of the resulting triangulation data structure.

/// @{


/*!
creates a new vertex `v` and uses it to star a hole.

Given a set of faces 'F' describing a simply connected hole (i.e., a topological disk),
the function deletes all the faces in `F`, creates a new vertex `v` and for each edge on
the boundary of the hole creates a new face with `v` as a vertex. The input is an iterator
range `[face_begin, face_end[` of `Face_handle`s over the connected faces in `F`. The handle
to the new vertex `v` is returned.

\pre `tds.dimension() = 2` and the set of faces has the topology of a disk.
*/
template< class FaceIt >
Vertex_handle insert_in_hole(FaceIt face_begin, FaceIt face_end);

/*!
same as above, except that `new_v` will be used as the new vertex, which must have been
allocated previously, for example with `create_vertex`.
*/
template< class FaceIt >
void insert_in_hole(Vertex_handle new_v, FaceIt face_begin, FaceIt face_end);


/*!
creates a new vertex `v` and use it to star the hole
whose boundary is described by the sequence of edges `[edge_begin, edge_end)`. Returns a handle to the vertex.
*/
template< class EdgeIt>
Vertex_handle star_hole(EdgeIt edge_begin,EdgeIt edge_end);

/*!
same as above, except that, to build the new faces, the algorithm
first recycles faces in the sequence `[face_begin, face_end)` and create new ones when the sequence is exhausted.
*/
template< class EdgeIt, class FaceIt>
Vertex_handle star_hole(EdgeIt edge_begin,
EdgeIt edge_end,
FaceIt face_begin,
FaceIt face_end);

/*!
uses vertex v to star the hole
whose boundary is described by the sequence of edges `[edge_begin, edge_end)`.

*/
template< class EdgeIt>
void star_hole(Vertex_handle v, EdgeIt edge_begin, EdgeIt edge_end);

/*!
same as above, recycling faces in the sequence `[face_begin, face_end)`.
*/
template< class EdgeIt, class FaceIt>
void star_hole(Vertex_handle v,
EdgeIt edge_begin,
EdgeIt edge_end,
FaceIt face_begin,
FaceIt face_end);

/*!
removes the vertex v, and store in `hole` the list of edges
on the boundary of the hole.
*/
void make_hole(Vertex_handle v, List_edges& hole);

/*!
adds a new vertex.
*/
Vertex_handle create_vertex();

/*!
adds a face which is the neighbor `i1` of `f1`,
`i2` of `f2` and `i3` of `f3`.
*/
Face_handle create_face(Face_handle f1, int i1, Face_handle f2, int i2, Face_handle
f3, int i3);

/*!
adds a face which is the neighbor `i1` of `f1`,
and the neighbor `i2` of `f2`.
*/
Face_handle create_face(Face_handle f1, int i1, Face_handle f2, int i2);

/*!
adds a face which is the neighbor `i1` of `f1`,
and has `v` as vertex.
*/
Face_handle create_face(Face_handle f1, int i1, Vertex_handle v);

/*!
adds a face with vertices `v1`, `v2` and `v3`.
*/
Face_handle create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3);

/*!
adds a face with vertices `v1`, `v2` and `v3`,
and neighbors `f1`, `f2`, `f3`.
*/
Face_handle create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
Face_handle f1, Face_handle f2, Face_handle f3);

/*!
adds a face whose vertices and neighbors are set to `Vertex_handle()` and `Face_handle()`.
*/
Face_handle create_face();

/*!
deletes a face.
*/
void delete_face(Face_handle );

/*!
deletes a vertex.
*/
void delete_vertex(Vertex_handle );

/// @}

/// \name Miscellaneous
/// @{

/*!
returns \f$ i+1\f$ modulo 3, with\f$ 0\leq i \leq2\f$.
*/
int ccw(int i) const;

/*!
returns \f$ i+2\f$ modulo 3, with \f$ 0\leq i \leq2\f$.
*/
int cw(int i) const;

/*!
checks the combinatorial validity of the
triangulation data structure: call the `is_valid()` member function for each vertex and
each face, checks the number of vertices and the Euler relation
between numbers of vertices, faces and edges.
*/
bool is_valid();

/*!
Returns the degree of `v` in the triangulation data structure.
*/
size_type degree(Vertex_handle v) const;

/*!
writes the triangulation data structure into the stream `os`.
If `v` is not `Vertex_handle()`, vertex `v`
is output first or skipped if `skip_first` is true.
*/
void file_output( ostream& os, Vertex_handle v = Vertex_handle(), bool skip_first=false);

/*!
inputs the triangulation data structure from file and returns a handle to the first input vertex.
If `skip_first` is true, it is assumed that the first
vertex has been omitted when output.
*/
Vertex_handle file_input( istream& is, bool skip_first=false);

/*!
reads a combinatorial triangulation data structure from `is` and assigns it to the triangulation data structure.
*/
istream& operator>> (istream& is, TriangulationDataStructure_2 & tds);

/*!
writes `tds` into the stream `os`.
*/
ostream& operator<< (ostream& os, const TriangulationDataStructure_2 & tds);

/// @}

}; /* end TriangulationDataStructure_2 */



/*!
\ingroup PkgTDS2Concepts
\cgalConcept

The concept `TriangulationDataStructure_2::Vertex` describes the type used by a
`TriangulationDataStructure_2` to store the vertices.

Some of the requirements listed below are of geometric nature
and are *optional*
when using the triangulation data structure class alone.
They became required when the triangulation data structure is plugged
into a triangulation.

\cgalHeading{Creation}

In order to obtain new vertices or destruct unused vertices, the user must
call the `create_vertex()` and `delete_vertex()` methods of the
triangulation data structure.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_vertex_base_2<Tds>}
\cgalHasModelsEnd

\sa `TriangulationDSVertexBase_2`
\sa `TriangulationDataStructure_2::Face`
\sa `TriangulationDataStructure_2`

*/
class TriangulationDataStructure_2::Vertex {
public:

/// \name Types
/// The class `TriangulationDataStructure_2::Vertex` defines the same
/// types as the triangulation data structure except the iterators.
/// @{

/*!
<I>Optional for the triangulation data
structure used alone</I>.
*/
typedef unspecified_type Point;

/// @}

/// \name Access Functions
/// @{

/*!
returns the geometric information of the vertex.
*/
Point point() const;

/*!
returns a face of the triangulation having `*this` as a vertex.
*/
Face_handle face() const;

/// @}

/// \name Setting
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
sets the geometric information to `p`.
\cgalAdvancedEnd
*/
void set_point(const Point& p);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
sets the incident face to `f`.
\cgalAdvancedEnd
*/
void set_face(Face_handle f);

/// @}

/// \name Checking
/// @{

/*!
Checks the validity of the vertex. Must check that its incident face
has this vertex. The validity of the base vertex is also checked.

When `verbose` is set to `true`, messages are printed to give
a precise indication on the kind of invalidity encountered.
*/
bool is_valid(bool verbose = false) const;

/*!
Inputs the non-combinatorial information possibly stored in the vertex.
*/
  std::istream& operator>> (std::istream& is, TriangulationDataStructure_2::Vertex & v);

/*!
Outputs the non combinatorial operation possibly stored in the
vertex.
*/
  std::ostream& operator<< (std::ostream& os, const TriangulationDataStructure_2::Vertex & v);

/// @}

}; /* end TriangulationDataStructure_2::Vertex */


/*!
\ingroup PkgTDS2Concepts
\cgalConcept

The concept `TriangulationDataStructure_2::Face` describes the types used to store the faces
face class of a
`TriangulationDataStructure_2`.
A `TriangulationDataStructure_2::Face` stores three handles to its three vertices
and three handles to its three neighbors.
The vertices are indexed 0,1, and 2 in counterclockwise order.
The neighbor indexed `i` lies
opposite to vertex `i`.

In degenerate cases,
when the triangulation data structure stores a
simplicial complex of dimension `0` and `1`,
the type `TriangulationDataStructure_2::Face` is used to store the faces
of maximal dimension of the complex, i.e., a vertex in dimension `0`, an edge in dimension `1`.
Only vertices and neighbors with index `0` are set in the first case,
only vertices and neighbors with index `0` or `1` are set in the second case.

\cgalHeading{Types}

The class `TriangulationDataStructure_2::Face` defines the same types as
the triangulation data structure
except the iterators and the circulators.

\cgalHeading{Creation}

The methods `create_face()` and
`delete_face()`
have to be used to
define new faces and to delete no longer used faces.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_face_base_2<Tds>}
\cgalHasModelsEnd

\sa `TriangulationDSFaceBase_2`
\sa `TriangulationDataStructure_2`
\sa `TriangulationDataStructure_2::Vertex`

*/
class TriangulationDataStructure_2::Face {
public:

/// \name Vertex Access Functions
/// @{

/*!
returns the vertex `i` of the face.
\pre \f$ 0\leq i \leq2\f$.
*/
Vertex_handle vertex(int i) const;

/*!
returns the index of vertex `v` in the face. \pre `v` is a vertex of the face.
*/
int index(Vertex_handle v) const;

/*!
returns `true` if `v` is a vertex of the face.
*/
bool has_vertex(Vertex_handle v) const;

/*!
returns `true` if `v` is a vertex of the face, and
computes the index `i` of `v` in the face.
*/
bool has_vertex(Vertex_handle v, int& i) const;

/// @}

/// \name Neighbor Access Functions
/// The neighbor with index `i` is the neighbor which is opposite to
/// the vertex with index `i`.
/// @{

/*!
returns the neighbor `i` of the face.
\pre \f$ 0\leq i \leq2\f$.

*/
Face_handle neighbor(int i) const;

/*!
returns the index of face `n`.
\pre `n` is a neighbor of the face.
*/
int index(Face_handle n) const;

/*!
returns `true` if `n` is a neighbor of the face.
*/
bool has_neighbor(Face_handle n) const;

/*!
returns `true` if `n` is a neighbor of the face, and
compute the index `i` of `n`.
*/
bool has_neighbor(Face_handle n, int& i) const;

/// @}

/// \name Setting
/// @{

/*!
sets vertex `i` to be `v`.
\pre \f$ 0\leq i \leq2\f$.

*/
void set_vertex(int i, Vertex_handle v);

/*!
sets neighbor `i` to be `n`.
\pre \f$ 0\leq i \leq2\f$.

*/
void set_neighbor(int i, Face_handle n);

/*!
sets the vertex handles to `Vertex_handle()`.
*/
void set_vertices();

/*!
sets the vertex handles.
*/
void set_vertices(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2);

/*!
sets the neighbors handles to `Face_handle()`.
*/
void set_neighbors();

/*!
sets the neighbors handles.
*/
void set_neighbors(Face_handle n0,
Face_handle n1,
Face_handle n2);

/// @}

/// \name Checking
/// @{

/*!
returns `true` if the function
`is_valid()` of the base class
returns `true` and if, for each index \f$ i\f$, \f$ 0 \le i < 3\f$,
the face is a neighbor of its neighboring face `neighbor(i)`
and shares with this neighbor the vertices `cw(i)` and `ccw(i)`
in correct reverse order.
*/
bool is_valid() const;

/// @}

/// \name Miscellaneous
/// @{

/*!
Returns \f$ i+1\f$ modulo 3, with \f$ 0\leq i \leq2\f$.
*/
int ccw(int i) const;

/*!
Returns \f$ i+2\f$ modulo 3, with  \f$ 0\leq i \leq2\f$.
*/
int cw(int i) const;

/*!
Inputs any non combinatorial information possibly stored in the face.
*/
std::istream& operator>> (std::istream& is, TriangulationDataStructure_2::Face & f);

/*!
Outputs any non combinatorial information possibly stored in the face.
*/
std::ostream& operator<< (std::ostream& os, const TriangulationDataStructure_2::Face & f);

/// @}

}; /* end TriangulationDataStructure_2::Face */

/*!
\ingroup PkgTDS2Concepts
\cgalConcept

Various algorithms using a triangulation data structure, such as Delaunay triangulations
or Alpha Shapes, must be able to associate a state to a face elemental.
For efficiency, this information must be stored directly within the face.

This class is only meant to store a state (Boolean). Consequently, the state must be the default
value (i.e. `false`) unless a setting function (`mark_in_conflict()`, etc.) has been called.

The three states are "in conflict", "on boundary", and "processed".
By default, a face is not in conflict, not on boundary, and not processed.

\sa `TriangulationDataStructure_2::Face`

*/

class TriangulationDataStructure_2::Face_data
{
public:
  /// \name Setting
  /// @{

  /// Clears all flags: the face is neither in conflict, nor on the boundary, nor processed.
  void clear();

  /// Sets the "in conflict" state to `true`.
  ///
  /// \post `is_in_conflict()` returns `true`
  void mark_in_conflict();

  /// Sets the "on boundary" state to `true`.
  ///
  /// \post `is_on_boundary()` returns `true`
  void mark_on_boundary();

  /// Sets the "processed" state to `true`.
  ///
  /// \post `processed()` returns `true`
  void mark_processed();

  /// @}

  /// \name Access Functions
  /// @{

  /// Checks whether the face has default state (not in conflict, not on boundary, not processed).
  bool is_clear();

  /// Returns whether the face has been marked as "in conflict".
  bool is_in_conflict();

  /// Returns whether the face has been marked as "on boundary".
  bool is_on_boundary();

  /// Returns whether the face has been marked as "processed".
  bool processed();

  /// @}

}; /* end Face_data */
