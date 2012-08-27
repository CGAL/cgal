
/*!
\ingroup PkgTDS2Concepts
\cgalconcept

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
`f`. See Figure \ref 2D_Triangulation_Fig_neighbors1 

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

I/O 
-------------- 

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

\hasModel `CGAL::Triangulation_data_structure_2<Vb,Fb>` 

\sa `TriangulationDataStructure_2::Face` 
\sa `TriangulationDataStructure_2::Vertex` 
\sa `CGAL::Triangulation_2<Traits,Tds>` 

*/

class TriangulationDataStructure_2 {
public:

/*!
\ingroup PkgTDS2Concepts
\cgalconcept

The concept `TriangulationDataStructure_2::Vertex` describes the type used by a 
`TriangulationDataStructure_2` to store the vertices. 

Some of the requirements listed below are of geometric nature 
and are *optional* 
when using the triangulation data structure class alone. 
They became required when the triangulation data structure is plugged 
into a triangulation. 

Creation 
-------------- 

In order to obtain new vertices or destruct unused vertices, the user must 
call the `create_vertex()` and `delete_vertex()` methods of the 
triangulation data structure. 

\hasModel CGAL::Triangulation_ds_vertex_2<Vb,Fb> 

\sa `TriangulationDataStructure_2` 
\sa `TriangulationDataStructure_2::Face` 

*/
class Vertex {
public:

/// \name Types 
/// The class `TriangulationDataStructure_2::Vertex` defines the same
/// types as the triangulation data structure except the iterators.
/// @{

/*! 
<I>Optional for the triangulation data 
structure used alone</I>.
*/ 
typedef Hidden_type Point; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the geometric information of `v`. 
*/ 
Point point() const; 

/*! 
returns a face of the triangulation having `v` as vertex. 
*/ 
Face_handle face() const; 

/// @} 

/// \name Setting 
/// @{

/*! 
\advanced sets the geometric information to `p`. 
*/ 
void set_point(const Point& p); 

/*! 
\advanced sets the incident face to `f`. 
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
\cgalconcept

The concept `TriangulationDataStructure_2::Face` describes the types used to store the faces 
face class of a 
`TriangulationDataStructure_2`. 
A `TriangulationDataStructure_2::Face` stores three handles to its three vertices 
and three handles to its three neighbors. 
The vertices are indexed 0,1, and 2 in counterclockwise order. 
The neighbor indexed \f$ i\f$ lies 
opposite to vertex i. 

In degenerate cases, 
when the triangulation data structure stores a 
simplicial complex of dimension \f$ 0\f$ and \f$ 1\f$, 
the type `TriangulationDataStructure_2::Face` is used to store the faces 
of maximal dimension of the complex 
: i.e. a vertex in dimension \f$ 0\f$, an edge in dimension \f$ 1\f$. 
Only vertices and neighbors with index \f$ 0\f$ are set in the first case, 
only vertices and neighbors with index \f$ 0\f$ or \f$ 1\f$ are set in the second case. 

Types 
-------------- 

The class `TriangulationDataStructure_2::Face` defines the same types as 
the triangulation data structure 
except the iterators and the circulators. 

Creation 
-------------- 

The methods `create_face` and 
`delete_face()` 
have to be used to 
define new faces and to delete non longer used faces. 

\sa `TriangulationDataStructure_2`
\sa `TriangulationDataStructure_2::Vertex`
\sa `TriangulationFaceBase_2`

*/
class Face {
public:

/// \name Vertex Access Functions 
/// @{

/*! 
returns the vertex `i` of `f`. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
Vertex_handle vertex(int i) const; 

/*! 
returns the index of vertex `v` in `f`. \pre `v` is a vertex of `f`. 
*/ 
int index(Vertex_handle v) const; 

/*! 
returns `true` if `v` is a vertex of `f`. 
*/ 
bool has_vertex(Vertex_handle v) const; 

/*! 
returns `true` if `v` is a vertex of `f`, and 
computes the index `i` of `v` in `f`. 
*/ 
bool has_vertex(Vertex_handle v, int& i) const; 

/// @} 

/// \name Neighbor Access Functions 
/// The neighbor with index `i` is the neighbor which is opposite to
/// the vertex with index `i`.
/// @{

/*! 
returns the neighbor `i` of `f`. 
\pre \f$ 0\leq i \leq2\f$. 

*/ 
Face_handle neighbor(int i) const; 

/*! 
returns the index of face `n`. 
\pre `n` is a neighbor of `f`. 
*/ 
int index(Face_handle n) const; 

/*! 
returns `true` if `n` is a neighbor of `f`. 
*/ 
bool has_neighbor(Face_handle n) const; 

/*! 
returns `true` if `n` is a neighbor of `f`, and 
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
sets the vertex handles to `NULL`. 
*/ 
void set_vertices(); 

/*! 
sets the vertex handles. 
*/ 
void set_vertices(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2); 

/*! 
sets the neighbors handles to `NULL`. 
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
face \f$ f\f$ is a neighbor of its neighboring face `neighbor(i)` 
and shares with this neighbor the vertices `cw(i)` and `ccw(i)` 
in correct reverse order. 
*/ 
bool is_valid() const; 

/// @} 

/// \name Miscellaneous 
/// @{

/*! 
Returns \f$ i+1\f$ modulo 3.\pre \f$ 0\leq i \leq2\f$. 
*/ 
int ccw(int i) const; 

/*! 
Returns \f$ i+2\f$ modulo 3.\pre \f$ 0\leq i \leq2\f$. 
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

/// \name Types 
/// @{

/*! 
Size type (unsigned integral type) 
*/ 
typedef Hidden_type size_type; 

/*! 
Difference type (signed integral type) 
*/ 
typedef Hidden_type difference_type; 

/*! 
Handle to a vertex.
\models ::Handle
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
Handle to a face. 
\models ::Handle
*/ 
typedef Hidden_type Face_handle; 

/*! 
This nested template class allows to get the type of a triangulation 
data structure that only changes the vertex type. It has to define a type 
`Other` which is a <I>rebound</I> triangulation data structure, that is, the 
one whose `TriangulationDSVertexBase_2` will be `Vb2`. 
*/ 
typedef Hidden_type template <typename Vb2> struct Rebind_vertex; 

/*! 
This nested template class allows to get the type of a triangulation 
data structure that only changes the face type. It has to define a type 
`Other` which is a <I>rebound</I> triangulation data structure, that is, the 
one whose `TriangulationDSFaceBase_2` will be `Fb2`. 
*/ 
typedef Hidden_type template <typename Fb2> struct Rebind_face; 

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
/// and faces of a triangulation data structure. They are all
/// bidirectional, non-mutable iterators. 
///
/// The circulators
/// allow to visit all the edges or faces incident to a given vertex
/// and all the vertices adjacent to a given vertex. 
///
/// They are all bidirectional and non mutable. Iterators and circulators are
/// convertible to the corresponding handles, thus they can be passed
/// directly as argument to the functions expecting a handle.
///
/// A face circulator is invalidated by any modification of the face it 
/// points to. An edge circulator is invalidated by any modification of 
/// anyone of the two faces incident to the edge pointed to. A vertex 
/// circulator that turns around vertex `v` and that has as value a handle 
/// to vertex `w`, is invalidated by any modification of anyone of the two 
/// faces incident to `v` and `w`.

/// @{


/*! 

*/ 
typedef Hidden_type Face_iterator; 

/*! 

*/ 
typedef Hidden_type Edge_iterator; 

/*! 

*/ 
typedef Hidden_type Vertex_iterator; 

/*! 

*/ 
typedef Hidden_type Face_circulator; 

/*! 

*/ 
typedef Hidden_type Edge_circulator; 

/*! 

*/ 
typedef Hidden_type Vertex_circulator; 

/// @} 

/// \name Creation 
/// @{

/*! 
Default constructor of `tds`. 
*/ 
TriangulationDataStructure_2(); 

/*! 
Copy constructor. All the vertices and faces are duplicated. 
*/ 
TriangulationDataStructure_2( const 
TriangulationDataStructure_2& tds1); 

/*! 
Assignment. All the vertices and faces of `tds1` are duplicated 
in `tds` . Former faces and vertices of `tds` , if any, are 
deleted.
*/ 
TriangulationDataStructure_2& operator=( const 
TriangulationDataStructure_2& tds1); 

/*! 
`tds1` is copied into `tds`. If \f$ v\, !\!= NULL\f$, the vertex of `tds` 
corresponding to `v` is returned, otherwise `Vertex_handle()` 
is returned. 
\pre The optional argument `v` is a vertex of `tds1`. 
*/ 
Vertex_handle 
copy_tds(const TriangulationDataStructure_2 & tds1, 
Vertex_handle v = Vertex_handle()); 

/*! 
Swaps `tds` and `tds1`. Should be preferred to `tds`=`tds1` or `tds`(`tds1`) 
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
  \advanced returns the dimension of the triangulation data structure. 
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
`v` is a vertex of `tds`. 
*/ 
bool is_vertex(Vertex_handle v) const; 

/*! 
returns `true` if `(fh,i)` is an edge of `tds`. Returns `false` when 
`dimension()` \f$ <1\f$ . 
*/ 
bool is_edge(Face_handle fh, int i) const; 

/*! 
returns `true` if 
`(va, vb)` is an edge of `tds`. 
*/ 
bool is_edge(Vertex_handle va, Vertex_handle vb) const; 

/*! 
as previous. In addition, if true is returned 
`fr` and `i` are set such that the pair `(fr,i)` 
is the description 
of the ordered edge `(va, vb)`. 
*/ 
bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle &fr, 
int &i) const; 

/*! 
returns `true` if `fh` is a face of `tds`. Returns `false` when 
`dimension()` \f$ <2\f$ . 
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
as above. In addition, if `true` is returned, `fr` is a handle 
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
\pre If the face `f` is given, it has to be incident to be a face of `tds` incident to `v` and the circulator begins with the vertex `f->vertex(ccw(i))` if `i` is the index of `v` in `f`. 
*/ 
Vertex_circulator 
incident_vertices(Vertex_handle v, Face_handle f=NULL) 
const; 

/*! 
\pre If the face `f` is given, it has to be a face of `tds` incident to `v` and the circulator begins with the edge `(f,cw(i))` of `f` if `i` is the index of `v` in `f`. 
*/ 
Edge_circulator 
incident_edges(Vertex_handle v, Face_handle f=NULL) const; 

/*! 
\pre If the face `f` is given, it has to be a face of `tds` incident to `v` and the circulator begins with the face `f`. 
*/ 
Face_circulator 
incident_faces(Vertex_handle v, Face_handle f=NULL) const; 

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

\image html Flip.gif "Flip"
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

\image html Three.gif "Insertion"
*/ 
void remove_degree_3(Vertex_handle v, Face *f=NULL); 

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

The link of a vertex \f$ v\f$ is formed by the edges disjoint from \f$ v\f$ 
that are included in the faces incident to \f$ v\f$. 
When the link of `v = f->vertex(i)` contains all the other vertices 
of the two-dimensional triangulation data structure (\f$ \mathbb{S}^2\f$), `dim_down` crushes the two-dimensional 
data-structure (\f$ \mathbb{S}^2\f$) onto the one-dimensional data structure (\f$ \mathbb{S}^1\f$) formed by the link of `v` 
augmented with the vertex `v` itself; this one is placed on the edge `(f, i)`
(see Fig. \ref figtdsdim_down_2). 
\pre `dimension()` must be equal to `2`, the degree of `f->vertex(i)` must be equal to the total number of vertices minus 1. 

\anchor figtdsdim_down_2
\image html tds-dim_down.png "From a two-dimensional data structure to a one-dimensional data structure."

*/
void dim_down(Face_handle f, int i); 

/// @} 

/// \name Advanced Modifiers 
/// The following modifiers are required for convenience of the advanced user. 
/// They do not guarantee the combinatorial validity of the resulting triangulation data structure.

/// @{


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
same as above, recycling faces in the sequence `[face_begin, face_end[ . ` 
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
adds a face whose vertices and neighbors are set to NULL. 
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
returns \f$ i+1\f$ modulo 3.\pre \f$ 0\leq i \leq2\f$. 
*/ 
int ccw(int i) const; 

/*! 
returns \f$ i+2\f$ modulo 3.\pre \f$ 0\leq i \leq2\f$. 
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
writes `tds` into the stream `os`. 
If `v` is not a null handle, vertex `v` 
is output first or skipped if `skip_first` is true. 
*/ 
void file_output( ostream& os, Vertex_handle v = Vertex_handle(), bool skip_first=false); 

/*! 
inputs `tds` from file and returns a handle to the first input vertex. 
If `skip_first` is true, it is assumed that the first 
vertex has been omitted when output. 
*/ 
Vertex_handle file_input( istream& is, bool skip_first=false); 

/*! 
reads a combinatorial triangulation data structure from `is` and assigns it to `tds`.
*/ 
istream& operator>> (istream& is, TriangulationDataStructure_3 & tds); 

/*! 
writes `tds` into the stream `os` 
*/ 
ostream& operator<< (ostream& os, const TriangulationDataStructure_3 & tds); 

/// @}

}; /* end TriangulationDataStructure_2 */

