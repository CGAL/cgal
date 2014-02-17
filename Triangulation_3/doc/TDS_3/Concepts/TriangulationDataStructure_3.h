
/*!
\ingroup PkgTDS3Concepts
\cgalConcept

3D-triangulation data structures are meant to maintain the 
combinatorial information for 3D-geometric triangulations. 

In \cgal, a triangulation data structure is a container of cells 
(\f$ 3\f$-faces) and vertices (\f$ 0\f$-faces). Following the standard 
vocabulary of simplicial complexes, an \f$ i\f$-face \f$ f_i\f$ and a \f$ j\f$-face 
\f$ f_j\f$ \f$ (0 \leq j < i \leq 3)\f$ are said to be <I>incident</I> in the 
triangulation if \f$ f_j\f$ is a (sub)face of \f$ f_i\f$, and two \f$ i\f$-faces \f$ (0 
\leq i \leq 3)\f$ are said to be <I>adjacent</I> if they share a 
common incident (sub)face. 

Each cell gives 
access to its four incident vertices and to its four adjacent 
cells. Each vertex gives direct access to one of its incident cells, which is 
sufficient to retrieve all the incident cells when needed. 

The four vertices of a cell are indexed with 0, 1, 2 and 3. The 
neighbors of a cell are also indexed with 0, 1, 2, 3 
in such a way that the neighbor indexed by \f$ i\f$ is opposite to the vertex 
with the same index (see \cgalFigureRef{TDS3figrepres}). 

Edges (\f$ 1\f$-faces) and facets (\f$ 2\f$-faces) are not explicitly 
represented: a facet is given by a cell and an index (the facet 
`i` of a cell `c` is the facet of `c` that is opposite to 
the vertex of index `i`) and an edge is given by a cell and two 
indices (the edge `(i,j)` of a cell `c` is the edge 
whose endpoints are the vertices of indices `i` and `j` of 
`c`). 

As \cgal explicitly deals with all degenerate cases, a 
3D-triangulation data structure in \cgal can handle the cases when 
the dimension of the triangulation is lower than 3 
(see Section \ref TDS3secintro). 

Thus, a 3D-triangulation data structure can store a triangulation of a 
topological sphere \f$ S^d\f$ of \f$ \mathbb{R}^{d+1}\f$, for any \f$ d \in \{-1,0,1,2,3\}\f$.<BR> 

 

The second template parameter of the basic triangulation class 
(see Chapter \ref chapterTriangulation3 "3D Triangulations") 
`CGAL::Triangulation_3` is a triangulation data structure class. (See 
Chapter \ref chapterTDS3.) 

To ensure all the *flexibility* of the class `CGAL::Triangulation_3`, a 
model of a triangulation data structure must be templated by the base vertex 
and the base cell classes (see \ref TDS3secintro): 
`TriangulationDataStructure_3<TriangulationVertexBase_3,TriangulationCellBase_3>`. 
The optional functionalities related to geometry are compulsory for 
this use as a template parameter of `CGAL::Triangulation_3`.


A class that satisfies the requirements for a triangulation data structure 
class must provide the following types and operations. 

\cgalHeading{I/O}

The information stored in the `iostream` is: 
the dimension, the number of vertices, the number of cells, 
the indices of the vertices of each cell, then the indices of the 
neighbors of each cell, where the index corresponds to the preceding 
list of cells. When dimension < 3, the same information is stored 
for faces of maximal dimension instead of cells. 

\cgalHasModel `CGAL::Triangulation_data_structure_3` 

\sa `TriangulationDataStructure_3::Vertex` 
\sa `TriangulationDataStructure_3::Cell` 

*/

class TriangulationDataStructure_3 {
public:

/// \name Types 
/// @{

/*!
  %Vertex type, requirements are described in `TriangulationDataStructure_3::Vertex`.
*/
typedef unspecified_type Vertex;

/*!
  %Cell type, requirements are described in `TriangulationDataStructure_3::Cell`.
*/
typedef unspecified_type Cell;

/*!
Size type (unsigned integral type) 
*/ 
typedef unspecified_type size_type; 

/*!
Difference type (signed integral type) 
*/ 
typedef unspecified_type difference_type; 

/*!

*/ 
typedef unspecified_type Vertex_handle; 

/*!

*/ 
typedef unspecified_type Cell_handle; 

/*!
Can be `CGAL::Sequential_tag` or `CGAL::Parallel_tag`. If it is 
`CGAL::Parallel_tag`, the following functions can be called concurrently:
`create_vertex`, `create_cell`, `delete_vertex`, `delete_cell`.
*/
typedef unspecified_type Concurrency_tag;

/*!
\cgalAdvancedBegin
This template class allows to get the type of a triangulation 
data structure that only changes the vertex type. It has to define a type 
`Rebind_vertex<Vb2>::%Other` which is a <I>rebound</I> triangulation data structure, that is, the 
one whose `TriangulationDSVertexBase_3` will be `Vb2`.
\note It can be implemented using a nested template class.
\cgalAdvancedEnd
*/ 
template <typename Vb2> 
using Rebind_vertex = unspecified_type;

/*!
\cgalAdvancedBegin
This template class allows to get the type of a triangulation 
data structure that only changes the cell type. It has to define a type 
`Rebind_cell<Cb2>::%Other` which is a <I>rebound</I> triangulation data structure, that is, the 
one whose `TriangulationDSCellBase_3` will be `Cb2`.
\note It can be implemented using a nested template class.
\cgalAdvancedEnd
*/ 
template <typename Cb2> 
using Rebind_cell = unspecified_type;

/*!
`(c,i,j)` is the 
edge of cell `c` whose vertices indices are `i` and 
`j`. (See Section \ref TDS3secintro.) 
*/ 
typedef Triple<Cell_handle, int, int> Edge; 

/*!
`(c,i)` is the facet 
of `c` opposite to the vertex of index `i`. (See 
Section \ref TDS3secintro.) 
*/ 
typedef std::pair<Cell_handle, int> Facet; 

/// @}

/// \name Iterators
/// The following iterators allow one to visit all the vertices,
/// edges, facets and cells of the triangulation data structure. They
/// are all bidirectional, non-mutable iterators. Iterators are convertible to the
/// corresponding handles, thus the user can pass them directly as
/// arguments to the functions.
/// @{

/*!

*/ 
typedef unspecified_type Cell_iterator; 

/*!

*/ 
typedef unspecified_type Facet_iterator; 

/*!

*/ 
typedef unspecified_type Edge_iterator; 

/*!

*/ 
typedef unspecified_type Vertex_iterator; 

/// @}

/// \name Circulators
/// The following circulators allow us to visit all the cells and
/// facets incident to a given edge. They are bidirectional and
/// non-mutable.Circulators are convertible to the
/// corresponding handles, thus the user can pass them directly as
/// arguments to the functions.
/// @{


/*!

*/ 
typedef unspecified_type Facet_circulator; 

/*!

*/ 
typedef unspecified_type Cell_circulator; 

/// @} 

/// \name Creation 
/// @{

/*!
Default constructor. 
*/ 
TriangulationDataStructure_3(); 

/*!
Copy constructor. All vertices and cells are duplicated. 
*/ 
TriangulationDataStructure_3(const TriangulationDataStructure_3 & tds1); 

/*!
Assignment operator. All vertices and cells are duplicated, and the former 
data structure of `tds` is deleted. 
*/ 
TriangulationDataStructure_3& operator= (const TriangulationDataStructure_3 & tds1); 

/*!
`tds1` is copied into `tds`. If `v != Vertex_handle()`, 
the vertex of `tds` corresponding to `v` is returned, 
otherwise `Vertex_handle()` is returned. 
\pre The optional argument `v` is a vertex of `tds1`. 
*/ 
Vertex_handle 
copy_tds(const TriangulationDataStructure_3 & tds1, 
Vertex_handle v = Vertex_handle()); 

/*!
`tds_src` is copied into `this`. As the vertex and cell types might be different
and incompatible, the creation of new cells and vertices is made thanks to the
functors `convert_vertex` and `convert_cell`, that convert vertex and cell types.
For each vertex `v_src` in `tds_src`, the corresponding vertex `v_tgt` in `this` is a
copy of the vertex returned by `convert_vertex(v_src)`. The same operations are
done for cells with the functor convert_cell. If `v != TDS_src::Vertex_handle()`,
a handle to the vertex created in `this` that is the copy of `v` is returned,
otherwise `Vertex_handle()` is returned.

 - A model of `ConvertVertex` must provide two operator()'s that are responsible for converting the source vertex `v_src` into the target vertex:
  - `Vertex operator()(const TDS_src::Vertex& v_src);` This operator is used to create the vertex from `v_src`.
  - `void operator()(const TDS_src::Vertex& v_src, Vertex& v_tgt);` This operator is meant to be used in case heavy data should transferred to `v_tgt`. 
 - A model of ConvertCell must provide two operator()'s that are responsible for converting the source cell `c_src` into the target cell:
  - `Cell operator()(const TDS_src::Cell& c_src);` This operator is used to create the cell from `c_src`.
  - `void operator()(const TDS_src::Cell& c_src, Cell& c_tgt);` This operator is meant to be used in case heavy data should transferred to `c_tgt`.

\pre The optional argument `v` is a vertex of `tds_src` or is `Vertex_handle()`.
*/
template <class TDS_src, class ConvertVertex, class ConvertCell>
Vertex_handle tds.copy_tds(const TDS_src& tds_src, typename TDS_src::Vertex_handle v, const ConvertVertex& convert_vertex, const ConvertCell& convert_cell);

/*!
Swaps `tds` and `tds1`. There is no copy of cells and vertices, 
thus this method runs in constant time. This method should be preferred to 
`tds`=`tds1` or `tds`(`tds1`) when `tds1` is deleted after 
that. 
*/ 
void swap(TriangulationDataStructure_3 & tds1); 

/*!
Deletes all cells and vertices. `tds` is reset as a triangulation 
data structure constructed by the default constructor. 
*/ 
void clear(); 

/// @} 

/// \name Access Functions 
/// @{

/*!
The dimension of the triangulated topological sphere. 
*/ 
int dimension() const; 

/*!
The number of vertices. Note that the triangulation data structure has one 
more vertex than an associated geometric triangulation, if there is 
one, since the infinite vertex is a standard vertex and is thus also 
counted. 
*/ 
size_type number_of_vertices() const; 

/*!
The number of cells. Returns 0 if `tds`.`dimension()`\f$ <3\f$. 
*/ 
size_type number_of_cells() const; 

/// @} 

/// \name Non constant-time access functions 
/// @{

/*!
The number of facets. Returns 0 if `tds`.`dimension()`\f$ <2\f$. 
*/ 
size_type number_of_facets() const; 

/*!
The number of edges. Returns 0 if `tds`.`dimension()`\f$ <1\f$. 
*/ 
size_type number_of_edges() const; 

/// @} 

/// \name Setting 
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Sets the dimension to `n`. 
\cgalAdvancedEnd
*/ 
void set_dimension(int n); 

/// @} 

/// \name Queries 
/// @{

/*!
Tests whether `v` is a vertex of `tds`. 
*/ 
bool is_vertex(Vertex_handle v) const; 

/*!
Tests whether `(c,i,j)` is an edge of `tds`. Answers `false` when 
`dimension()` \f$ <1\f$ . 
\pre \f$ i,j \in\{0,1,2,3\}\f$ 
*/ 
bool is_edge(Cell_handle c, int i, int j) const; 

/*!
Tests whether `(u,v)` is an edge of `tds`. If the edge is found, 
it computes a cell `c` having this edge and the indices `i` 
and `j` of the vertices `u` and `v`, in this order. 
*/ 
bool is_edge(Vertex_handle u, Vertex_handle v, 
Cell_handle & c, int & i, int & j) const; 

/*!
Tests whether `(u,v)` is an edge of `tds`. 
*/ 
bool is_edge(Vertex_handle u, Vertex_handle v) const; 

/*!
Tests whether `(c,i)` is a facet of `tds`. Answers `false` when 
`dimension()` \f$ <2\f$ . 
\pre \f$ i \in\{0,1,2,3\}\f$ 
*/ 
bool is_facet(Cell_handle c, int i) const; 

/*!
Tests whether `(u,v,w)` is a facet of `tds`. If the facet is found, 
it computes a cell `c` having this facet and the indices `i`, 
`j` and `k` of the vertices `u`, `v` and `w`, in 
this order. 
*/ 
bool is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w, 
Cell_handle & c, int & i, int & j, int & k) const; 

/*!
Tests whether `c` is a cell of `tds`. Answers `false` when 
`dimension()` \f$ <3\f$ . 
*/ 
bool is_cell(Cell_handle c) const; 

/*!
Tests whether `(u,v,w,t)` is a cell of `tds`. If the cell 
`c` is found, it computes the indices `i`, `j`, `k` 
and `l` of the vertices `u`, `v`, `w` and `t` in 
`c`, in this order. 
*/ 
bool is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w, Vertex_handle t, 
Cell_handle & c, int & i, int & j, int & k, int & l) const; 

/// \name has_vertex
/// There is a method `has_vertex` in the cell class. The analogous methods for facets are defined here. 
/// @{

/*!
If `v` is a vertex of `f`, then `j` is the index of 
`v` in the cell `f.first`, and the method returns `true`. 
\pre `tds`.dimension()=3 
*/ 
bool has_vertex(const Facet & f, Vertex_handle v, int & j) const; 

/*!
Same for facet `(c,i)`. Computes the index `j` of `v` in 
`c`. 
*/ 
bool has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const; 

/*!

*/ 
bool has_vertex(const Facet & f, Vertex_handle v) const; 

/*!
Same as the first two methods, but these two methods do not return the 
index of the vertex. 
*/ 
bool has_vertex(Cell_handle c, int i, Vertex_handle v) const; 

/// @}

/// \name Equality Tests
/// The following three methods test whether two facets have the same vertices.
/// @{

/*!

*/ 
bool are_equal(const Facet & f, const Facet & g) const; 

/*!

*/ 
bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const; 

/*!
For these three methods: \pre `tds`.dimension()=3. 
*/ 
bool are_equal(const Facet & f, Cell_handle n, int j) const; 

/// @} 

/*!
\name Flips 

Two kinds of flips exist for a three-dimensional triangulation. They
are reciprocal. To be flipped, an edge must be incident to three
tetrahedra. During the flip, these three tetrahedra disappear and two
tetrahedra appear. Figure \ref TDS3figflips (left) shows the edge that
is flipped as bold dashed, and one of its three incident facets is
shaded. On the right, the facet shared by the two new tetrahedra is
shaded.

\anchor TDS3figflips
\image html flips.png "Flips."
\image latex flips.png "Flips."

The following methods guarantee the validity of the resulting 3D
combinatorial triangulation. Moreover the flip operations do not
invalidate the vertex handles, and only invalidate the cell handles of
the affected cells.

<I>Flips for a 2d triangulation are not implemented yet</I>

*/
/// @{

/*!

*/ 
bool flip(Edge e); 

/*!
Before flipping, these methods check that edge `e=(c,i,j)` is 
flippable (which is quite expensive). They return `false` or 
`true` according to this test. 
*/ 
bool flip(Cell_handle c, int i, int j); 

/*!

*/ 
void flip_flippable(Edge e); 

/*!
Should be preferred to the previous methods when the edge is 
known to be flippable. 
\pre The edge is flippable. 
*/ 
void flip_flippable(Cell_handle c, int i, int j); 

/*!

*/ 
bool flip(Facet f); 

/*!
Before flipping, these methods check that facet `f=(c,i)` is 
flippable (which is quite expensive). They return `false` or 
`true` according to this test. 
*/ 
bool flip(Cell_handle c, int i); 

/*!

*/ 
void flip_flippable(Facet f); 

/*!
Should be preferred to the previous methods when the facet is 
known to be flippable. 
\pre The facet is flippable. 
*/ 
void flip_flippable(Cell_handle c, int i); 

/// @} 

/// \name Insertions 
/// The following modifier member functions guarantee the combinatorial validity of the resulting triangulation. 
/// @{

/*!
Creates a new vertex, inserts it in cell `c` and returns its handle. 
The cell `c` is split into four new cells, each of these cells being 
formed by the new vertex and a facet of `c`. 
\pre `tds`.`dimension()` \f$ = 3\f$ and `c` is a cell of `tds`. 
*/ 
Vertex_handle insert_in_cell(Cell_handle c); 

/*!
Creates a new vertex, inserts it in facet `f` and returns its handle. 
In dimension 3, the two incident cells are split into 3 new cells; 
in dimension 2, the facet is split into 3 facets. 
\pre `tds`.`dimension()` \f$ \geq2\f$ and `f` is a facet of `tds`. 
*/ 
Vertex_handle insert_in_facet(const Facet & f); 

/*!
Creates a new vertex, inserts it in facet `i` of `c` and returns its 
handle. 
\pre `tds`.`dimension()` \f$ \geq2\f$, \f$ i \in\{0,1,2,3\}\f$ in dimension 3, \f$ i=3\f$ in dimension 2 and `(c,i)` is a facet of `tds`. 
*/ 
Vertex_handle insert_in_facet(Cell_handle c, int i); 

/*!
Creates a new vertex, inserts it in edge `e` and returns its handle. 
In dimension 3, all the 
incident cells are split into 2 new cells; in dimension 2, the 2 
incident facets are split into 2 new facets; in dimension 1, the edge is 
split into 2 new edges. 
\pre `tds`.`dimension()` \f$ \geq1\f$ and `e` is an edge of `tds`. 
*/ 
Vertex_handle insert_in_edge(Edge e); 

/*!
Creates a new vertex, inserts it in edge \f$ (i,j)\f$ of `c` and returns its 
handle. 
\pre `tds`.`dimension()` \f$ \geq1\f$. \f$ i\neq j\f$, \f$ i,j \in\{0,1,2,3\}\f$ in dimension 3, \f$ i,j \in\{0,1,2\}\f$ in dimension 2, \f$ i,j \in\{0,1\}\f$ in dimension 1 and `(c,i,j)` is an edge of `tds`. 
*/ 
Vertex_handle insert_in_edge(Cell_handle c, int i, int j); 

/*!
Transforms a triangulation of the sphere \f$ S^d\f$ of \f$ \mathbb{R}^{d+1}\f$ into the 
triangulation of the sphere \f$ S^{d+1}\f$ of \f$ \mathbb{R}^{d+2}\f$ by adding a new vertex 
`v`: 
`v` is linked to all the vertices to triangulate one of the two 
half-spheres of dimension \f$ (d+1)\f$. Vertex `star` is used to 
triangulate the second half-sphere (when there is an associated 
geometric triangulation, `star` is in fact the vertex associated with 
its infinite vertex). 
See Figure \ref TDS3figtopoinsert_outside_affine_hull. 

The numbering of the cells is such that, if `f` was a face of 
maximal dimension in the initial triangulation, then `(f,v)` (in 
this order) is the corresponding face in the new triangulation. 
This method can be used to insert the first two vertices in an empty 
triangulation. 

A handle to `v` is returned. 
\pre `tds`.`dimension()` \f$ = d < 3\f$. When `tds`.`number_of_vertices()` \f$ >0\f$, \f$ star \neq\f$ `Vertex_handle()` and `star` is a vertex of `tds`. 


\anchor TDS3figtopoinsert_outside_affine_hull 
\image html topo-insert_outside_affine_hull.png "insert_increase_dimension (1-dimensional case)."
\image latex topo-insert_outside_affine_hull.png "insert_increase_dimension (1-dimensional case)."

*/ 
Vertex_handle 
insert_increase_dimension(Vertex_handle star = Vertex_handle()); 

/*!
Creates a new vertex by starring a hole. It takes an iterator range 
[`cell_begin`; `cell_end`[ of `Cell_handles` which specifies a set 
of connected cells (resp. facets in dimension 2) describing a hole. 
(`begin`, `i`) is a facet (resp. an edge) on the boundary of the hole, 
that is, `begin` belongs to the set of cells (resp. facets) previously 
described, and `begin->neighbor(i)` does not. Then this function deletes 
all the cells (resp. facets) describing the hole, creates a new vertex 
`v`, and for each facet (resp. edge) on the boundary of the hole, creates 
a new cell (resp. facet) with `v` as vertex. `v` is returned. 
\pre `tds`.`dimension()` \f$ \geq2\f$, the set of cells (resp. facets) is connected, and its boundary is connected. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i); 

/*!
Same as above, except that `newv` will be used as the new vertex, which 
must have been allocated previously with e.g. `create_vertex`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i, Vertex_handle newv); 

/// @} 

/// \name Removal 
/// @{

/*!
This operation is the reciprocal of `insert_increase_dimension()`. 
It transforms a triangulation of the sphere \f$ S^d\f$ of \f$ \mathbb{R}^{d+1}\f$ into the 
triangulation of the sphere \f$ S^{d-1}\f$ of \f$ \mathbb{R}^{d}\f$ by removing the vertex 
`v`. Delete the cells incident to `w`, keep the others. 
\pre `tds`.`dimension()` \f$ = d \geq-1\f$. `tds`.`degree(v)` \f$ =\f$ `degree(w)` \f$ =\f$ `tds`.`number_of_vertices()` \f$ -1\f$. 

*/ 
void remove_decrease_dimension(Vertex_handle v, Vertex_handle w = v); 

/*!
Removes `v`. The incident simplices of maximal dimension incident to 
`v` are replaced by a single simplex of the same dimension. This 
operation is exactly the reciprocal to `tds`.`insert_in_cell(v)` in 
dimension 3, `tds`.`insert_in_facet(v)` in dimension 2, and 
`tds`.`insert_in_edge(v)` in dimension 1. 
\pre `tds`.`degree(v)` \f$ =\f$ `tds`.`dimension()+1`. 

*/ 
Cell_handle remove_from_maximal_dimension_simplex(Vertex_handle v); 

/// @} 

/// \name Dimension Manipulation 
/// The following operation, <TT>decrease_dimension</TT>, is necessary when the displacement of a vertex decreases the dimension of the triangulation.
/// @{

/*!
The link of a vertex \f$ v\f$ is formed by the facets 
disjoint from \f$ v\f$ that are included in the cells incident to \f$ v\f$. When the link of `v = c->vertex(i)` contains all the other vertices, `decrease_dimension` crushes the 
triangulation of the sphere \f$ S^d\f$ of \f$ \mathbb{R}^{d+1}\f$ onto the 
triangulation of the sphere \f$ S^{d-1}\f$ of \f$ \mathbb{R}^{d}\f$ formed by the link of `v` 
augmented with the vertex `v` itself, for \f$ d\f$==2,3; this one is placed on the facet `(c, i)` 
(see Fig. \ref TDS3dim_down). 
\pre The dimension must be 2 or 3. The degree of `v` must be equal to the total number of vertices of the triangulation data structure minus 1. 

\anchor TDS3dim_down
\image html tds-dim_down.png
\image latex tds-dim_down.png

<center><b>From an \f$ S^d\f$ data structure to an \f$ S^{d-1}\f$ data
structure (top: \f$ d==2\f$, bottom: \f$ d==3\f$).
</b></center>

*/ 
void decrease_dimension(Cell_handle c, int i); 

/// @} 

/// \name Other modifiers 
/// The following modifiers can affect the validity of the triangulation data structure.
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Changes the orientation of all cells of the triangulation data structure. 
\cgalAdvancedEnd
\pre `tds`.`dimension()` \f$ \geq1\f$. 
*/ 
void reorient(); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Adds a copy of the vertex `v` to the triangulation data structure. 
\cgalAdvancedEnd
*/ 
Vertex_handle create_vertex(const Vertex &v = Vertex()); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Creates a vertex which is a copy of the one pointed to by `v` 
and adds it to the triangulation data structure. 
\cgalAdvancedEnd
*/ 
Vertex_handle create_vertex(Vertex_handle v); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Adds a copy of the cell `c` to the triangulation data structure. 
\cgalAdvancedEnd
*/ 
Cell_handle create_cell(const Cell &c = Cell()); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Creates a cell which is a copy of the one pointed to by `c` 
and adds it to the triangulation data structure. 
\cgalAdvancedEnd
*/ 
Cell_handle create_cell(Cell_handle c); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Creates a cell and adds it into the triangulation data 
structure. Initializes the vertices of the cell, its neighbor handles 
being initialized with the default constructed handle.
\cgalAdvancedEnd
*/
Cell_handle create_cell(Vertex_handle v0, Vertex_handle v1, 
Vertex_handle v2, Vertex_handle v3); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Creates a cell, initializes its vertices and neighbors, and adds it 
into the triangulation data structure. 
\cgalAdvancedEnd
*/ 
Cell_handle create_cell( Vertex_handle v0, Vertex_handle v1, 
Vertex_handle v2, Vertex_handle v3, 
Cell_handle n0, Cell_handle n1, 
Cell_handle n2, Cell_handle n3); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Removes the vertex from the triangulation data structure. 
\cgalAdvancedEnd
\pre The vertex is a vertex of `tds`. 
*/ 
void delete_vertex( Vertex_handle v ); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Removes the cell from the triangulation data structure. 
\cgalAdvancedEnd
\pre The cell is a cell of `tds`.
*/ 
void delete_cell( Cell_handle c ); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Calls `delete_vertex` over an iterator range of value type `Vertex_handle`.
\cgalAdvancedEnd
*/ 
template <class VertexIt> 
void delete_vertices(VertexIt first, VertexIt last); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Calls `delete_cell` over an iterator range of value type `Cell_handle`. 
\cgalAdvancedEnd
*/ 
template <class CellIt> 
void delete_cells(CellIt first, CellIt last); 

/// @} 

/// \name Traversing the triangulation 
/// @{

/*!
Returns `cells_end()` when `tds.dimension()` \f$ <3\f$. 
*/ 
Cell_iterator cells_begin() const; 

/*!

*/ 
Cell_iterator cells_end() const; 

/*!
Low-level access to the cells, does not return `cells_end()` 
when `tds.dimension()` \f$ <3\f$. 
*/ 
Cell_iterator raw_cells_begin() const; 

/*!

*/ 
Cell_iterator raw_cells_end() const; 

/*!
Returns `facets_end()` when `tds.dimension()` \f$ <2\f$. 
*/ 
Facet_iterator facets_begin() const; 

/*!

*/ 
Facet_iterator facets_end() const; 

/*!
Returns `edges_end()` when `tds.dimension()` \f$ <1\f$. 
*/ 
Edge_iterator edges_begin() const; 

/*!

*/ 
Edge_iterator edges_end() const; 

/*!

*/ 
Vertex_iterator vertices_begin() const; 

/*!

*/ 
Vertex_iterator vertices_end() const; 

/// @} 

/// \name Cell and facet circulators 
/// @{

/*!
Starts at an arbitrary cell incident to `e`. 
\pre `tds.dimension()` \f$ =3\f$ 
*/ 
Cell_circulator incident_cells(const Edge & e) const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Cell_circulator incident_cells(Cell_handle c, int i, int j) const; 

/*!
Starts at cell `start`. 
\pre `tds.dimension()` \f$ =3\f$ and `start` is incident to `e`. 
*/ 
Cell_circulator incident_cells(const Edge & e, Cell_handle start) const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Cell_circulator incident_cells(Cell_handle c, int i, int j, Cell_handle start) 
const; 

/*!
Starts at an arbitrary facet incident to `e`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.

\pre `tds.dimension()` \f$ =3\f$ 
*/ 
Facet_circulator incident_facets(Edge e) const; 

/*!
As above for edge `(i,j)` of `c`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j) const; 

/*!
Starts at facet `start`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.

\pre `start` is incident to `e`. 
*/ 
Facet_circulator incident_facets(Edge e, Facet start) const; 

/*!
Starts at facet of index `f` in `start`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.
*/ 
Facet_circulator incident_facets(Edge e, Cell_handle start, int f) const; 

/*!
As above for edge `(i,j)` of `c`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j, 
Facet start) const; 

/*!
As above for edge `(i,j)` of `c` and facet `(start,f)`. 

Only defined in dimension 3, though are defined also in dimension 2:
there are only two facets sahring an edge in dimension 2.
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j, 
Cell_handle start, int f) const; 

/// @} 

/// \name Traversal of the incident cells, facets and edges, and the adjacent vertices of a given vertex
/// @{

/*!
Copies the `Cell_handle`s of all cells incident to `v` to the 
output iterator `cells`. 
Returns the resulting output iterator. 
\pre `tds.dimension()` \f$ =3\f$, `v` \f$ \neq\f$ `Vertex_handle()`, `tds.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_cells(Vertex_handle v, OutputIterator cells) const; 

/*!
Copies the `Facet`s incident to `v` to the output iterator 
`facets`. 
Returns the resulting output iterator. 
\pre `tds.dimension()` \f$ >1\f$, `v` \f$ \neq\f$ `Vertex_handle()`, `tds.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_facets(Vertex_handle v, OutputIterator facets) const; 

/*!
Copies all `Edge`s incident to `v` to the 
output iterator `edges`. Returns the resulting output iterator. 
\pre `tds.dimension()` \f$ >0\f$, `v` \f$ \neq\f$ `Vertex_handle()`, `tds.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_edges(Vertex_handle v, OutputIterator edges) const; 

/*!
Copies the `Vertex_handle`s of all vertices adjacent to `v` to the 
output iterator `vertices`. If `tds.dimension()` \f$ <0\f$, then do 
nothing. Returns the resulting output iterator. 
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `tds.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
adjacent_vertices(Vertex_handle v, OutputIterator vertices) const; 

/*!
Returns the degree of a vertex, that is, the number of incident vertices. 
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `tds.is_vertex(v)`. 
*/ 
size_type degree(Vertex_handle v) const; 

/// @} 

/// \name Traversal between adjacent cells 
/// @{

/*!
Returns the index of `c` in its \f$ i^{th}\f$ neighbor. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
int mirror_index(Cell_handle c, int i) const; 

/*!
Returns the vertex of the \f$ i^{th}\f$ neighbor of `c` that is opposite to 
`c`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
Vertex_handle mirror_vertex(Cell_handle c, int i) const; 

/*!
Returns the same facet seen from the other adjacent cell. 
*/ 
Facet mirror_facet(Facet f) const; 

/// @} 

/// \name Checking 
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the combinatorial validity of the triangulation by checking 
the local validity of all its cells and vertices (see functions below). 
(See Section \ref TDS3secintro.) Moreover, the Euler relation is 
tested. 

When `verbose` is set to `true`, messages are printed to give 
a precise indication on the kind of invalidity encountered. 
\cgalDebugEnd
*/ 
bool is_valid(bool verbose = false) const; 

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the local validity of the adjacency relations of the triangulation. 
It also calls the `is_valid` member function of the vertex. 
When `verbose` is set to `true`, messages are printed to give 
a precise indication on the kind of invalidity encountered. 
\cgalDebugEnd
*/ 
bool is_valid(Vertex_handle v, bool verbose = false) const; 

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the local validity of the adjacency relations of the triangulation. 
It also calls the `is_valid` member function of the cell. 
When `verbose` is set to `true`, messages are printed to give 
a precise indication on the kind of invalidity encountered. 
\cgalDebugEnd
*/ 
bool is_valid(Cell_handle c, bool verbose = false) const; 

/*!
Reads a combinatorial triangulation from `is` and assigns it to `tds` 
*/ 
istream& operator>> (istream& is, TriangulationDataStructure_3 & tds); 

/*!
Writes `tds` into the stream `os` 
*/ 
ostream& operator<< (ostream& os, const TriangulationDataStructure_3 & tds); 

/// @}

}; /* end TriangulationDataStructure_3 */

/*!
\ingroup PkgTDS3Concepts
\cgalConcept

The concept `TriangulationDataStructure_3::Vertex` represents the vertex class of a 3D-triangulation 
data structure. It must define 
the types and operations listed in this section. Some of these 
requirements are of geometric nature, they are <I>optional</I> 
when using the triangulation data structure class alone. They become 
compulsory when the triangulation data structure is used as a layer 
for the geometric triangulation class. (See Section \ref TDS3secdesign.) 

\cgalHeading{Creation}

In order to obtain new vertices or destruct unused vertices, the user must 
call the `create_vertex()` and `delete_vertex()` methods of the 
triangulation data structure. 

\sa `TriangulationDataStructure_3::Cell`

*/
class TriangulationDataStructure_3::Vertex {
public:

/// \name Types 
/// The class `Vertex` defines types that are the same as some of the
/// types defined by the triangulation data structure class
/// `TriangulationDataStructure_3`.
/// @{

/*!
<I>Optional for the triangulation data structure alone</I>.
*/ 
typedef unspecified_type Point; 

/*!

*/ 
typedef TriangulationDataStructure_3 Triangulation_data_structure; 

/*!

*/ 
typedef TriangulationDataStructure_3::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef TriangulationDataStructure_3::Cell_handle Cell_handle; 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns a cell of the triangulation having `v` as vertex. 
*/ 
Cell_handle cell() const; 

/*!
Returns the point stored in the vertex. 
<I>Optional for the triangulation data structure alone.</I> 
*/ 
Point point() const; 

/// @} 

/// \name Setting 
/// @{

/*!
Sets the incident cell to `c`. 
*/ 
void set_cell(Cell_handle c); 

/*!
Sets the point to `p`. <I>Optional for the 
triangulation data structure alone.</I> 
*/ 
void set_point(const Point & p); 

/// @} 

/// \name Checking 
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the validity of the vertex. Must check that its incident cell 
has this vertex. The validity of the base vertex is also checked. 

When `verbose` is set to `true`, messages are printed to give 
a precise indication on the kind of invalidity encountered. 
\cgalDebugEnd
*/ 
bool is_valid(bool verbose = false) const; 

/// @}

}; /* end Vertex */

/*!
\ingroup PkgTDS3Concepts
\cgalConcept

The concept `TriangulationDataStructure_3::Cell` stores 
four `Vertex_handle`s to its four vertices and four `Cell_handle`s 
to its four neighbors. The vertices are indexed 0, 1, 2, and 3 in consistent 
order. The neighbor indexed \f$ i\f$ lies opposite to vertex `i`. 

In degenerate dimensions, cells are used to store faces of maximal 
dimension: in dimension 2, each cell represents only one 
facet of index 3, and 3 edges \f$ (0,1)\f$, \f$ (1,2)\f$ and \f$ (2,0)\f$; in 
dimension 1, each cell represents one edge \f$ (0,1)\f$. (See also 
Section \ref TDS3secintro.) 

\cgalHeading{Creation}

In order to obtain new cells or destruct unused cells, the user must call the 
`create_cell()` and `delete_cell()` methods of the triangulation data 
structure. 

\sa `TriangulationDataStructure_3::Vertex`

*/

class TriangulationDataStructure_3::Cell {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef TriangulationDataStructure_3 Triangulation_data_structure; 

/*!

*/ 
typedef TriangulationDataStructure_3::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef TriangulationDataStructure_3::Cell_handle Cell_handle; 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns the vertex `i` of `c`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
Vertex_handle vertex(int i) const; 

/*!
Returns the index of vertex `v` in `c`. 
\pre `v` is a vertex of `c`. 
*/ 
int index(Vertex_handle v) const; 

/*!
Returns `true` if `v` is a vertex of `c`. 
*/ 
bool has_vertex(Vertex_handle v) const; 

/*!
Returns `true` if `v` is a vertex of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_vertex(Vertex_handle v, int & i) const; 

/*!
Returns the neighbor `i` of `c`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
Cell_handle neighbor(int i) const; 

/*!
Returns the index corresponding to neighboring cell `n`. 
\pre `n` is a neighbor of `c`. 
*/ 
int index(Cell_handle n) const; 

/*!
Returns `true` if `n` is a neighbor of `c`. 
*/ 
bool has_neighbor(Cell_handle n) const; 

/*!
Returns `true` if `n` is a neighbor of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_neighbor(Cell_handle n, int & i) const; 

/// @} 

/// \name Setting 

/// @{

/*!
Sets vertex `i` to `v`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
void set_vertex(int i, Vertex_handle v); 

/*!
Sets the vertex pointers. 
*/ 
void set_vertices(Vertex_handle v0, 
Vertex_handle v1, 
Vertex_handle v2, 
Vertex_handle v3); 

/*!
Sets neighbor `i` to `n`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
void set_neighbor(int i, Cell_handle n); 

/*!
Sets the neighbors pointers. 
*/ 
void set_neighbors(Cell_handle n0, 
Cell_handle n1, 
Cell_handle n2, 
Cell_handle n3); 

/// @} 

/// \name Checking 
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
User defined local validity checking function.
\cgalDebugEnd
*/ 
bool is_valid(bool verbose = false, int level = 0) const; 

/// @}

}; /* end Cell */
