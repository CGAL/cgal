
namespace CGAL {

/*!
\ingroup PkgTriangulation3TriangulationClasses

The class `Triangulation_3` represents a 3-dimensional tetrahedralization 
of points. 

\tparam TriangulationTraits_3 is the geometric traits class.

\tparam TriangulationDataStructure_3 is the triangulation data structure.
It has the default value `Triangulation_data_structure_3< Triangulation_vertex_base_3<TriangulationTraits_3>,Triangulation_cell_base_3<TriangulationTraits_3> >`. 

\tparam SurjectiveLockDataStructure is an optional parameter to specify the type of the spatial lock data structure.
        It is only used if the triangulation data structure used is concurrency-safe (i.e.\ when 
        TriangulationDataStructure_3::Concurrency_tag is Parallel_tag).
        It must be a model of the `SurjectiveLockDataStructure` concept,
        with `Object` being a `Point`.
        The default value is `Spatial_lock_grid_3<Tag_priority_blocking>` if
        the triangulation data structure is concurrency-safe, and `void` otherwise.
        In order to use concurrent operations, the user must provide a
        reference to a `SurjectiveLockDataStructure`
        instance via the constructor or `Triangulation_3::set_lock_data_structure`.

\cgalHeading{Traversal of the Triangulation}

The triangulation class provides several iterators and circulators 
that allow one to traverse it (completely or partially). 

\sa `TriangulationDataStructure_3::Vertex` 
\sa `TriangulationDataStructure_3::Cell` 

*/
template< typename TriangulationTraits_3, typename TriangulationDataStructure_3,
          typename SurjectiveLockDataStructure >
class Triangulation_3 : public Triangulation_utils_3 {
public:

  /*!
  The enum `Locate_type` is defined by `Triangulation_3` to specify 
  which case occurs when locating a point in the triangulation.
  */
  enum Locate_type {VERTEX=0, EDGE, FACET, CELL, OUTSIDE_CONVEX_HULL, OUTSIDE_AFFINE_HULL};

public:

/// \name Types 
/// The class `Triangulation_3` defines the following types:   
/// @{

/*!

*/ 
typedef TriangulationDataStructure_3 Triangulation_data_structure; 

/*!

*/ 
typedef SurjectiveLockDataStructure Lock_data_structure; 

/*!

*/ 
typedef TriangulationTraits_3 Geom_traits; 

/*!

*/ 
typedef TriangulationTraits_3::Point_3 Point; 

/*!

*/ 
typedef TriangulationTraits_3::Segment_3 Segment; 

/*!

*/ 
typedef TriangulationTraits_3::Triangle_3 Triangle; 

/*!

*/ 
typedef TriangulationTraits_3::Tetrahedron_3 Tetrahedron; 

/// @}

/*! \name
Only vertices (0-faces) and cells (3-faces) are stored. Edges (1-faces) and facets (2-faces) are not explicitly represented and thus there are no corresponding classes (see Section \ref Triangulation3secintro).
   */

/// @{

/*!

*/ 
typedef TriangulationDataStructure_3::Vertex Vertex; 

/*!

*/ 
typedef TriangulationDataStructure_3::Cell Cell; 

/*!

*/ 
typedef TriangulationDataStructure_3::Facet Facet; 

/*!

*/ 
typedef TriangulationDataStructure_3::Edge Edge; 


/// @}

/*! \name 

The vertices and faces of the triangulations are accessed through
`handles`, `iterators` and `circulators`. A handle is a model of the
`Handle` concept, and supports the two dereference operators and
`operator->`. A circulator is a model of the concept
`Circulator`. Iterators and circulators are bidirectional and
non-mutable. The edges and facets of the triangulation can also be
visited through iterators and circulators which are bidirectional and
non-mutable. Iterators and circulators are convertible to the
corresponding handles, thus the user can pass them directly as
arguments to the functions.  
*/
/// @{

/*!
handle to a vertex 
*/ 
typedef TriangulationDataStructure_3::Vertex_handle Vertex_handle; 

/*!
handle to a cell 
*/ 
typedef TriangulationDataStructure_3::Cell_handle Cell_handle; 

/*!
Reference to a simplex (vertex, edge, facet or cell) of the triangulation 
*/ 
typedef Triangulation_simplex_3<Self> Simplex; 

/*!
Size type (an unsigned integral type) 
*/ 
typedef TriangulationDataStructure_3::size_type size_type; 

/*!
Difference type (a signed integral type) 
*/ 
typedef TriangulationDataStructure_3::difference_type difference_type; 

/*!
iterator over cells 
*/ 
typedef TriangulationDataStructure_3::Cell_iterator All_cells_iterator; 

/*!
iterator over facets 
*/ 
typedef TriangulationDataStructure_3::Facet_iterator All_facets_iterator; 

/*!
iterator over edges 
*/ 
typedef TriangulationDataStructure_3::Edge_iterator All_edges_iterator; 

/*!
iterator over vertices 
*/ 
typedef TriangulationDataStructure_3::Vertex_iterator All_vertices_iterator; 

/*!
iterator over finite cells 
*/ 
typedef unspecified_type Finite_cells_iterator; 

/*!
iterator over finite facets 
*/ 
typedef unspecified_type Finite_facets_iterator; 

/*!
iterator over finite edges 
*/ 
typedef unspecified_type Finite_edges_iterator; 

/*!
iterator over finite vertices 
*/ 
typedef unspecified_type Finite_vertices_iterator; 

/*!
iterator over the points corresponding to the 
finite vertices of the triangulation. 
*/ 
typedef unspecified_type Point_iterator; 

/*!
circulator over all cells incident to a given edge 
*/ 
typedef TriangulationDataStructure_3::Cell_circulator Cell_circulator; 

/*!
circulator over all facets incident to a given edge 
*/ 
typedef TriangulationDataStructure_3::Facet_circulator Facet_circulator; 

/*! 
Concurrency tag (from the TDS).
*/ 
typedef TriangulationDataStructure_3::Concurrency_tag Concurrency_tag;

/// @} 

/// \name Creation 
/// @{

/*!
Introduces a triangulation `t` having only one vertex which is the 
infinite vertex.
`lock_ds` is an optional pointer to the lock data structure for parallel operations. It
must be provided if concurrency is enabled.
*/ 
Triangulation_3 
(const TriangulationTraits_3 & traits = TriangulationTraits_3(), 
 Lock_data_structure *lock_ds = NULL);

/*! 
Same as the previous one, but with parameters in reverse order.
*/ 
Triangulation_3 
(Lock_data_structure *lock_ds = NULL,
 const TriangulationTraits_3 & traits = TriangulationTraits_3());

/*!
Copy constructor. All vertices and faces are duplicated. 
The pointer to the lock data structure is not copied. Thus, the copy won't be
concurrency-safe as long as the user has not call `Triangulation_3::set_lock_data_structure`.
*/ 
Triangulation_3 (const Triangulation_3 & tr); 

/*!
Equivalent to constructing an empty triangulation with the optional 
traits class argument and calling `insert(first,last)`. 
*/ 
template < class InputIterator> 
Triangulation_3 (InputIterator first, InputIterator last, 
const TriangulationTraits_3 & traits = TriangulationTraits_3(),
Lock_data_structure *lock_ds = NULL); 

/// @} 

/// \name Assignment 
/// @{

/*!
The triangulation `tr` is duplicated, and modifying the copy after the 
duplication does not modify the original. The previous triangulation held 
by `t` is deleted. 
*/ 
Triangulation_3 & operator=(const Triangulation_3 & tr); 

/*!
The triangulations `tr` and `t` are swapped. 
`t.swap(tr)` should be preferred to `t = tr` or to 
`t(tr)` if `tr` is deleted after that. Indeed, there is no 
copy of cells and vertices, thus this method runs in constant time. 
*/ 
void swap(Triangulation_3 & tr); 

/*!
Deletes all finite vertices and all cells of `t`. 
*/ 
void clear(); 

/*!
Equality operator. Returns `true` iff there exist a bijection between the 
vertices of `t1` and those of `t2` and a bijection between the cells of 
`t1` and those of `t2`, which preserve the geometry of the 
triangulation, that is, the points of each corresponding pair of vertices are 
equal, and the tetrahedra corresponding to each pair of cells are equal (up to 
a permutation of their vertices). 
*/ 
template < class GT, class Tds1, class Tds2 > 
bool operator==(const Triangulation_3<GT, Tds1> & t1, const Triangulation_3<GT, Tds2> & t2); 

/*!
The opposite of `operator==`. 
*/ 
template < class GT, class Tds1, class Tds2 > 
bool operator!=(const Triangulation_3<GT, Tds1> & t1, const Triangulation_3<GT, Tds2> & t2); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns a const reference to the geometric traits object. 
*/ 
const TriangulationTraits_3 & geom_traits() const; 

/*!
Returns a const reference to the triangulation data structure. 
*/ 
const TriangulationDataStructure_3 & tds() const; 


/*!
Returns a reference to the triangulation data structure.
\cgalAdvancedBegin
This method is mainly a help for users implementing their own triangulation algorithms.
The responsibility of keeping a valid triangulation belongs to the user when using advanced operations allowing a direct manipulation of the `tds`.
\cgalAdvancedEnd
*/
TriangulationDataStructure_3 & tds(); 

/*!
Returns the dimension of the affine hull. 
*/ 
int dimension() const; 

/*!
Returns the number of finite vertices. 
*/ 
size_type number_of_vertices() const; 

/*!
Returns the number of cells or 0 if `t.dimension() < 3`. 
*/ 
size_type number_of_cells() const; 

/*!
Returns the infinite vertex. 
*/ 
Vertex_handle infinite_vertex(); 

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
This method is meant to be used only if you have done a low-level operation on the underlying tds that invalidated the infinite vertex.
Sets the infinite vertex.
\cgalAdvancedEnd
*/ 
void set_infinite_vertex(Vertex_handle v);

/*!
Returns a cell incident to the infinite vertex. 
*/ 
Cell_handle infinite_cell() const; 

/// @} 

/// \name Non-Constant-Time Access Functions 
/// As previously said, the triangulation is a collection of cells that are either infinite or represent a finite tetrahedra, where an infinite cell is a cell incident to the infinite vertex. Similarly we call an edge (resp. facet) `infinite` if it is incident to the infinite vertex.
/// @{

/*!
The number of facets. Returns 0 if `t.dimension() < 2`. 
*/ 
size_type number_of_facets() const; 

/*!
The number of edges. Returns 0 if `t.dimension() < 1`. 
*/ 
size_type number_of_edges() const; 

/*!
The number of finite cells. Returns 0 if `t.dimension() < 3`. 
*/ 
size_type number_of_finite_cells() const; 

/*!
The number of finite facets. Returns 0 if `t.dimension() < 2`. 
*/ 
size_type number_of_finite_facets() const; 

/*!
The number of finite edges. Returns 0 if `t.dimension() < 1`. 
*/ 
size_type number_of_finite_edges() const; 

/// @} 

/// \name Geometric Access Functions 
/// @{

/*!
Returns the tetrahedron formed by the four vertices of `c`. 
\pre `t.dimension() == 3` and the cell is finite. 
*/ 
Tetrahedron tetrahedron(Cell_handle c) const; 

/*!
Returns the triangle formed by the three vertices of facet 
`(c,i)`. The triangle is oriented so that its normal points to the 
inside of cell `c`. 
\pre `t.dimension()` \f$ \geq2\f$ and \f$ i \in\{0,1,2,3\}\f$ in dimension 3, \f$ i = 3\f$ in dimension 2, and the facet is finite. 
*/ 
Triangle triangle(Cell_handle c, int i) const; 

/*!
Same as the previous method for facet `f`. 
\pre `t.dimension()` \f$ \geq2\f$ and the facet is finite. 
*/ 
Triangle triangle(const Facet & f) const; 

/*!
Returns the line segment formed by the vertices of `e`. 
\pre `t.dimension()` \f$ \geq1\f$ and `e` is finite. 
*/ 
Segment segment(const Edge & e) const; 

/*!
Same as the previous method for edge `(c,i,j)`. 
\pre As above and \f$ i\neq j\f$. Moreover \f$ i,j \in\{0,1,2,3\}\f$ in dimension 3, \f$ i,j \in\{0,1,2\}\f$ in dimension 2, \f$ i,j \in\{0,1\}\f$ in dimension 1, and the edge is finite. 
*/ 
Segment segment(Cell_handle c, int i, int j) const; 

/*!
Returns the point given by vertex `i` of cell `c`. 
\pre `t.dimension()` \f$ \geq0\f$ and \f$ i \in\{0,1,2,3\}\f$ in dimension 3, \f$ i \in\{0,1,2\}\f$ in dimension 2, \f$ i \in\{0,1\}\f$ in dimension 1, \f$ i = 0\f$ in dimension 0, and the vertex is finite. 
*/ 
const Point & point(Cell_handle c, int i) const; 

/*!
Same as the previous method for vertex `v`. 
\pre `t.dimension()` \f$ \geq0\f$ and the vertex is finite. 
*/ 
const Point & point(Vertex_handle v) const; 

/// @} 

/// \name Tests for Finite and Infinite Vertices and Faces 
/// @{

/*!
`true`, iff vertex `v` is the infinite vertex. 
*/ 
bool is_infinite(Vertex_handle v) const; 

/*!
`true`, iff `c` is incident to the infinite vertex. 
\pre `t.dimension() == 3`. 
*/ 
bool is_infinite(Cell_handle c) const; 

/*!
`true`, iff the facet `i` of cell `c` is incident to the 
infinite vertex. 
\pre `t.dimension()` \f$ \geq2\f$ and \f$ i\in\{0,1,2,3\}\f$ in dimension 3, \f$ i=3\f$ in dimension 2. 
*/ 
bool is_infinite(Cell_handle c, int i) const; 

/*!
`true` iff facet `f` is incident to the infinite vertex. 
\pre `t.dimension()` \f$ \geq2\f$. 
*/ 
bool is_infinite(const Facet & f) const; 

/*!
`true`, iff the edge `(i,j)` of cell `c` is incident to 
the infinite vertex. 
\pre `t.dimension()` \f$ \geq1\f$ and \f$ i\neq j\f$. Moreover \f$ i,j \in\{0,1,2,3\}\f$ in dimension 3, \f$ i,j \in\{0,1,2\}\f$ in dimension 2, \f$ i,j \in\{0,1\}\f$ in dimension 1. 
*/ 
bool is_infinite(Cell_handle c, int i, int j) const; 

/*!
`true` iff edge `e` is incident to the infinite vertex. 
\pre `t.dimension()` \f$ \geq1\f$. 
*/ 
bool is_infinite(const Edge & e) const; 

/// @} 

/// \name Queries 

/// @{

/*!
Tests whether `p` is a vertex of `t` by locating `p` in 
the triangulation. If `p` is found, the associated vertex `v` 
is given. 
*/ 
bool is_vertex(const Point & p, Vertex_handle & v) const; 

/*!
Tests whether `v` is a vertex of `t`. 
*/ 
bool is_vertex(Vertex_handle v) const; 

/*!
Tests whether `(u,v)` is an edge of `t`. If the edge is found, 
it gives a cell `c` having this edge and the indices `i` 
and `j` of the vertices `u` and `v` in `c`, in this order. 
\pre `u` and `v` are vertices of `t`. 
*/ 
bool is_edge(Vertex_handle u, Vertex_handle v, 
Cell_handle & c, int & i, int & j) const; 

/*!
Tests whether `(u,v,w)` is a facet of `t`. If the facet is found, 
it computes a cell `c` having this facet and the indices `i`, 
`j` and `k` of the vertices `u`, `v` and `w` in `c`, 
in this order. 
\pre `u`, `v` and `w` are vertices of `t`. 
*/ 
bool is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w, 
Cell_handle & c, int & i, int & j, int & k) const; 

/*!
Tests whether `c` is a cell of `t`. 
*/ 
bool is_cell(Cell_handle c) const; 

/*!
Tests whether `(u,v,w,x)` is a cell of `t`. 
If the cell `c` is found, the method 
computes the indices `i`, `j`, `k` and `l` of the 
vertices `u`, `v`, `w` and `x` in `c`, in this 
order. 
\pre `u`, `v`, `w` and `x` are vertices of `t`. 
*/ 
bool is_cell(Vertex_handle u, Vertex_handle v, 
Vertex_handle w, Vertex_handle x, 
Cell_handle & c, 
int & i, int & j, int & k, int & l) const; 

/*!
Tests whether `(u,v,w,x)` is a cell of `t` and computes 
this cell `c`. 
\pre `u`, `v`, `w` and `x` are vertices of `t`. 
*/ 
bool is_cell(Vertex_handle u, Vertex_handle v, 
Vertex_handle w, Vertex_handle x, 
Cell_handle & c) const; 

/// @}

/*! \name
There is a method `has_vertex()` in the cell class. The analogous methods for facets are defined here. 
*/
/// @{
/*!
If `v` is a vertex of `f`, then `j` is the index of 
`v` in the cell `f.first`, and the method returns `true`. 
\pre `t.dimension() == 3`
*/ 
bool has_vertex(const Facet & f, Vertex_handle v, int & j) const; 

/*!
Same for facet `(c,i)`. Computes the index `j` of `v` in 
`c`. 
*/ 
bool has_vertex(Cell_handle c, int i, 
Vertex_handle v, int & j) const; 

/*!

*/ 
bool has_vertex(const Facet & f, Vertex_handle v) const; 

/*!
Same as the first two methods, but these two methods do not return the 
index of the vertex. 
*/ 
bool has_vertex(Cell_handle c, int i, Vertex_handle v) const; 

/// @}

/*! \name
The following three methods test whether two facets have the same vertices. 
 */

/// @{
/*!

*/ 
bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const; 

/*!

*/ 
bool are_equal(const Facet & f, const Facet & g) const; 

/*!
For these three methods: \pre `t.dimension() == 3`. 
*/ 
bool are_equal(const Facet & f, Cell_handle n, int j) const; 

/// @} 

/// \name Point Location 
/// The class `Triangulation_3` provides two functions to locate a given point with respect to a triangulation. It provides also functions to test if a given point is inside a finite face or not. Note that the class `Delaunay_triangulation_3` also provides a `nearest_vertex()` function.
/// @{

/*!

If the point `query` lies inside the convex hull of the points, the cell 
that contains the query in its interior is returned. If `query` lies on a 
facet, an edge or on a vertex, one of the cells having `query` on 
its boundary is returned. 

If the point `query` lies outside the convex hull of the points, 
an infinite cell with vertices \f$ \{ p, q, r, \infty\}\f$ is returned such that 
the tetrahedron \f$ ( p, q, r, query )\f$ is positively oriented 
(the rest of the triangulation lies on the other side of facet 
\f$ ( p, q, r )\f$). 

Note that locate works even in degenerate dimensions: in dimension 2 
(resp. 1, 0) the `Cell_handle` returned is the one that represents 
the facet (resp. edge, vertex) containing the query point. 

The optional argument `start` is used as a starting place for the search. 

The optional argument `could_lock_zone` is used by the concurrency-safe
version of the triangulation. When the pointer is not null, the locate will
try to lock all the cells along the walk. If it succeeds, `*could_lock_zone`
is `true`, otherwise it is false. In any case, the locked cells are not
unlocked by `locate`, leaving this choice to the user.
*/ 
Cell_handle 
locate(const Point & query, Cell_handle start = Cell_handle(),
       bool *could_lock_zone = NULL) const; 

/*!
Same as above but uses `hint` as the starting place for the search. 
*/ 
Cell_handle 
locate(const Point & query, Vertex_handle hint,
       bool *could_lock_zone = NULL) const; 

/*!
Same as `locate()` but uses inexact predicates. 
This function returns a handle on a cell that is a good approximation of the exact
location of `query`, while being faster. Note that it may return a handle on a cell 
whose interior does not contain `query`.
When the triangulation has dimension smaller than 3, `start` is returned.

Note that this function is available only if the cartesian coordinates of `query`
are accessible with functions `x()`, `y()` and `z()`.
*/
Cell_handle
inexact_locate(const Point & query, Cell_handle start = Cell_handle()) const;

/*!
If `query` lies inside the affine hull of the points, the \f$ k\f$-face 
(finite or infinite) that contains `query` in its interior is 
returned, by means of the cell returned together with `lt`, which 
is set to the locate type of the query (`VERTEX, EDGE, FACET, CELL`, or `OUTSIDE_CONVEX_HULL` if the cell is infinite and `query` 
lies strictly in it) and two indices `li` and `lj` that 
specify the \f$ k\f$-face of the cell containing `query`. 

If the \f$ k\f$-face is a cell, `li` and `lj` have no 
meaning; if it is a facet (resp. vertex), `li` gives the index of 
the facet (resp. vertex) and `lj` has no meaning; if it is and 
edge, `li` and `lj` give the indices of its vertices. 

If the point `query` lies outside the affine hull of the points, 
which can happen in case of degenerate dimensions, `lt` is set to 
`OUTSIDE_AFFINE_HULL`, and the cell returned has no meaning. 
As a particular case, if there is no finite vertex yet in the 
triangulation, `lt` is set to `OUTSIDE_AFFINE_HULL` and 
<I>locate</I> returns the default constructed handle. 

The optional argument `start` is used as a starting place for the search. 

The optional argument `could_lock_zone` is used by the concurrency-safe
version of the triangulation. When the pointer is not null, the locate will
try to lock all the cells along the walk. If it succeeds, `*could_lock_zone`
is `true`, otherwise it is false. In any case, the locked cells are not
unlocked by `locate`, leaving this choice to the user.
*/ 
Cell_handle 
locate(const Point & query, Locate_type & lt, 
int & li, int & lj, Cell_handle start = Cell_handle(),
bool *could_lock_zone = NULL ) const; 

/*!
Same as above but uses `hint` as the starting place for the search. 
*/ 
Cell_handle 
locate(const Point & query, Locate_type & lt, 
int & li, int & lj, Vertex_handle hint,
bool *could_lock_zone = NULL) const; 


/*!
Returns a value indicating on which side of the oriented boundary 
of `c` the point `p` lies. More precisely, it returns: 

- `ON_BOUNDED_SIDE` if `p` is inside the cell. For an infinite 
cell this means that `p` lies strictly in the half space limited by 
its finite facet and not containing any other point of the triangulation. 

- `ON_BOUNDARY` if p on the boundary of the cell. For an infinite 
cell this means that `p` lies on the <I>finite</I> facet. Then 
`lt` together with `li` and `lj` give the precise location 
on the boundary. (See the descriptions of the <I>locate</I> methods.) 

- `ON_UNBOUNDED_SIDE` if `p` lies outside the cell. For an 
infinite cell this means that `p` does not satisfy either of the 
two previous conditions. 
\pre `t.dimension() == 3` 
*/ 
Bounded_side 
side_of_cell(const Point & p, 
Cell_handle c, 
Locate_type & lt, int & li, int & lj) const; 

/*!
Returns a value indicating on which side of the oriented boundary 
of `f` the point `p` lies: 

- `ON_BOUNDED_SIDE` if `p` is inside the facet. For an 
infinite facet this means that `p` lies strictly in the half plane 
limited by its finite edge and not containing any other point of the 
triangulation . 

- `ON_BOUNDARY` if `p` is on the boundary of the facet. 
For an infinite facet this means that `p` lies on the finite 
edge. `lt`, `li` and `lj` give the precise location of 
`p` on the boundary of the facet. `li` and `lj` refer to 
indices in the degenerate cell `c` representing `f`. 

- `ON_UNBOUNDED_SIDE` if `p` lies outside the facet. For 
an infinite facet this means that `p` does not satisfy either of 
the two previous conditions. 

\pre `t.dimension() == 2` and `p` lies in the plane containing the triangulation. `f.second` \f$ =3\f$ (in dimension 2 there is only one facet per cell). 
*/ 
Bounded_side 
side_of_facet(const Point & p, 
const Facet & f, 
Locate_type & lt, int & li, int & lj) const; 

/*!
Same as the previous method for the facet `(c,3)`. 
*/ 
Bounded_side 
side_of_facet(const Point & p, 
Cell_handle c, 
Locate_type & lt, int & li, int & lj) const; 

/*!
Returns a value indicating on which side of the oriented boundary 
of `e` the point `p` lies: 

- `ON_BOUNDED_SIDE` if `p` is inside the edge. For an 
infinite edge this means that `p` lies in the half line defined by 
the vertex and not containing any other point of the triangulation. 

- `ON_BOUNDARY` if `p` equals one of the vertices, 
`li` give the index of the vertex in the cell storing `e` 

- `ON_UNBOUNDED_SIDE` if `p` lies outside the edge. For 
an infinite edge this means that `p` lies on the other half line, 
which contains the other points of the triangulation. 
\pre `t.dimension() == 1` and `p` is collinear with the points of the triangulation. `e.second == 0` and `e.third` \f$ =1\f$ (in dimension 1 there is only one edge per cell). 
*/ 
Bounded_side 
side_of_edge(const Point & p, 
const Edge & e, 
Locate_type & lt, int & li) const; 

/*!
Same as the previous method for edge \f$ (c,0,1)\f$. 
*/ 
Bounded_side 
side_of_edge(const Point & p, 
Cell_handle c, 
Locate_type & lt, int & li) const; 

/// @} 

/*! \name Flips 

Two kinds of flips exist for a three-dimensional triangulation. They
are reciprocal. To be flipped, an edge must be incident to three
tetrahedra. During the flip, these three tetrahedra disappear and two
tetrahedra appear.  \cgalFigureRef{Triangulation3figflips} (left) shows the
edge that is flipped as bold dashed, and one of its three incident
facets is shaded. On the right, the facet shared by the two new
tetrahedra is shaded. Flips are possible only under the following
conditions: - the edge or facet to be flipped is not on the boundary
of the convex hull of the triangulation - the five points involved are
in convex position. 

\cgalFigureBegin{Triangulation3figflips, flips.png}
Flips
\cgalFigureEnd

The following methods guarantee the validity of the resulting 3D
triangulation. Flips for a 2d triangulation are not implemented yet.

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
/// The following operations are guaranteed to lead to a valid triangulation when they are applied on a valid triangulation.

/// @{

/*!
Inserts point `p` in the triangulation and returns the corresponding 
vertex. 

If point `p` coincides with an already existing vertex, this 
vertex is returned and the triangulation remains unchanged. 

If point `p` lies in the convex hull of the points, it is added 
naturally: if it lies inside a cell, the cell is split into four 
cells, if it lies on a facet, the two incident cells are split into 
three cells, if it lies on an edge, all the cells incident to this 
edge are split into two cells. 

If point `p` is strictly outside the convex hull but in the affine 
hull, `p` is linked to all visible points on the convex hull to 
form the new triangulation. See 
Figure \ref Triangulation3figinsert_outside_convex_hull. 

If point `p` is outside the affine hull of the points, `p` is 
linked to all the points, and the dimension of the triangulation is 
incremented. All the points now belong to the boundary of the convex 
hull, so, the infinite vertex is linked to all the points to 
triangulate the new infinite face. See 
Figure \ref Triangulation3figinsert_outside_affine_hull. 
The optional argument `start` is used as a starting place for the search. 
*/ 
Vertex_handle insert(const Point & p, 
Cell_handle start = Cell_handle() ); 

/*!
Same as above but uses `hint` as the starting place for the search. 
*/ 
Vertex_handle insert(const Point & p, Vertex_handle hint); 

/*!
Inserts point `p` in the triangulation and returns the corresponding 
vertex. Similar to the above `insert()` function, but takes as additional 
parameter the return values of a previous location query. See description of 
<I>locate()</I> above. 
*/ 
Vertex_handle insert(const Point & p, Locate_type lt, 
Cell_handle loc, int li, int lj); 

/*!
Inserts the points in the range `[first,last)`. Returns the number of inserted points. 
Note that this function is not guaranteed to insert the points 
following the order of `InputIterator`. 
\tparam InputIterator must be an input iterator with value type `Point`. 
*/ 
template < class InputIterator > 
std::ptrdiff_t 
insert(InputIterator first, InputIterator last); 


/// @} 

/*! \name 
We also provide some other methods that can be used instead of
`Triangulatation_3::insert()` when the place where the new point must
be inserted is already known. They are also guaranteed to lead to a
valid triangulation when they are applied on a valid triangulation.
*/

// @{
/*!
Inserts point `p` in cell `c`. Cell `c` is split into 4 
tetrahedra. 
\pre `t.dimension() == 3` and `p` lies strictly inside cell `c`. 
*/ 
Vertex_handle insert_in_cell(const Point & p, Cell_handle c); 

/*!
Inserts point `p` in facet `f`. In dimension 3, the 2 
neighboring cells are split into 3 tetrahedra; in dimension 2, the facet 
is split into 3 triangles. 
\pre `t.dimension()` \f$ \geq2\f$ and `p` lies strictly inside face `f`. 
*/ 
Vertex_handle insert_in_facet(const Point & p, const Facet & f); 

/*!
As above, insertion in facet `(c,i)`. 
\pre As above and \f$ i \in\{0,1,2,3\}\f$ in dimension 3, \f$ i = 3\f$ in dimension 2. 
*/ 
Vertex_handle insert_in_facet(const Point & p, 
Cell_handle c, int i); 

/*!
Inserts `p` in edge `e`. In dimension 3, 
all the cells having this edge are split into 2 tetrahedra; in 
dimension 2, the 2 neighboring facets are split into 2 triangles; in 
dimension 1, the edge is split into 2 edges. 
\pre `t.dimension()` \f$ \geq1\f$ and `p` lies on edge `e`. 
*/ 
Vertex_handle insert_in_edge(const Point & p, const Edge & e); 

/*!
As above, inserts `p` in edge \f$ (i, j)\f$ of `c`. 
\pre As above and \f$ i\neq j\f$. Moreover \f$ i,j \in\{0,1,2,3\}\f$ in dimension 3, \f$ i,j \in\{0,1,2\}\f$ in dimension 2, \f$ i,j \in\{0,1\}\f$ in dimension 1. 
*/ 
Vertex_handle insert_in_edge(Point p, Cell_handle c, int i, int j); 

/*!
The cell `c` must be an infinite cell containing `p`. 

Links `p` to all points in the triangulation that are visible from 
`p`. Updates consequently the infinite faces. See 
Figure \ref Triangulation3figinsert_outside_convex_hull. 
\pre `t.dimension() > 0`, `c`, and the \f$ k\f$-face represented by `c` is infinite and contains `t`. 

\anchor Triangulation3figinsert_outside_convex_hull
\image html insert_outside_convex_hull.png "insert_outside_convex_hull() (2-dimensional case)"
\image latex insert_outside_convex_hull.png "insert_outside_convex_hull() (2-dimensional case)"
*/ 
Vertex_handle insert_outside_convex_hull(const Point & p, 
Cell_handle c); 

/*!
`p` is linked to all the points, and the infinite vertex is linked 
to all the points (including `p`) to triangulate the new infinite 
face, so that all the points now belong to the boundary of the convex 
hull. See Figure \ref Triangulation3figinsert_outside_affine_hull. 

This method can be used to insert the first point in an empty 
triangulation. 
\pre `t.dimension() < 3` and `p` lies outside the affine hull of the points. 

\anchor Triangulation3figinsert_outside_affine_hull
\image html insert_outside_affine_hull.png "insert_outside_affine_hull() (2-dimensional case)"
\image latex insert_outside_affine_hull.png "insert_outside_affine_hull() (2-dimensional case)"
*/ 
Vertex_handle insert_outside_affine_hull(const Point & p); 

/*!
Creates a new vertex by starring a hole. It takes an iterator range 
`[cell_begin,cell_end)` of `Cell_handle`s which specifies 
a hole: a set of connected cells (resp. facets in dimension 2) which is 
star-shaped wrt `p`. 
(`begin`, `i`) is a facet (resp. an edge) on the boundary of the hole, 
that is, `begin` belongs to the set of cells (resp. facets) previously 
described, and `begin->neighbor(i)` does not. Then this function deletes 
all the cells (resp. facets) describing the hole, creates a new vertex 
`v`, and for each facet (resp. edge) on the boundary of the hole, creates 
a new cell (resp. facet) with `v` as vertex. Then `v->set_point(p)` 
is called and `v` is returned. 

This operation is equivalent to calling 
`tds().insert_in_hole(cell_begin, cell_end, begin, i); v->set_point(p)`. 
\pre `t.dimension()` \f$ \geq2\f$, the set of cells (resp. facets in dimension 2) is connected, its boundary is connected, and `p` lies inside the hole, which is star-shaped wrt `p`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(Point p, CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i); 

/*!
Same as above, except that `newv` will be used as the new vertex, which 
must have been allocated previously with e.g.\ `create_vertex`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(Point p, CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i, Vertex_handle newv); 

/// @} 

/*! \name Cell, Face, Edge and Vertex Iterators 
The following iterators allow the user to visit cells, facets, edges and vertices of the triangulation. These iterators are non-mutable, bidirectional and their value types are respectively `Cell`, `Facet`, `Edge` and `Vertex`. They are all invalidated by any change in the triangulation.
*/
/// @{

/*!
Starts at an arbitrary finite vertex. Then `++` and `--` will 
iterate over finite vertices. Returns `finite_vertices_end()` when 
`t.number_of_vertices() == 0`. 
*/ 
Finite_vertices_iterator finite_vertices_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_vertices_iterator finite_vertices_end() const; 

/*!
Starts at an arbitrary finite edge. Then `++` and `--` will 
iterate over finite edges. Returns `finite_edges_end()` when 
`t.dimension() < 1`. 
*/ 
Finite_edges_iterator finite_edges_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_edges_iterator finite_edges_end() const; 

/*!
Starts at an arbitrary finite facet. Then `++` and `--` will 
iterate over finite facets. Returns `finite_facets_end()` when 
`t.dimension() < 2`. 
*/ 
Finite_facets_iterator finite_facets_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_facets_iterator finite_facets_end() const; 

/*!
Starts at an arbitrary finite cell. Then `++` and `--` will 
iterate over finite cells. Returns `finite_cells_end()` when 
`t.dimension() < 3`. 
*/ 
Finite_cells_iterator finite_cells_begin() const; 

/*!
Past-the-end iterator 
*/ 
Finite_cells_iterator finite_cells_end() const; 

/*!
Starts at an arbitrary vertex. Iterates over all vertices (even the infinite 
one). Returns `vertices_end()` when 
`t.number_of_vertices() == 0`. 
*/ 
All_vertices_iterator all_vertices_begin() const; 

/*!
Past-the-end iterator 
*/ 
All_vertices_iterator all_vertices_end() const; 

/*!
Starts at an arbitrary edge. Iterates over all edges (even infinite 
ones). Returns `edges_end()` when `t.dimension() < 1`. 
*/ 
All_edges_iterator all_edges_begin() const; 

/*!
Past-the-end iterator 
*/ 
All_edges_iterator all_edges_end() const; 

/*!
Starts at an arbitrary facet. Iterates over all facets (even infinite 
ones). Returns `facets_end()` when 
`t.dimension() < 2`. 
*/ 
All_facets_iterator all_facets_begin() const; 

/*!
Past-the-end iterator 
*/ 
All_facets_iterator all_facets_end() const; 

/*!
Starts at an arbitrary cell. Iterates over all cells (even infinite 
ones). Returns `cells_end()` when 
`t.dimension() < 3`. 
*/ 
All_cells_iterator all_cells_begin() const; 

/*!
Past-the-end iterator 
*/ 
All_cells_iterator all_cells_end() const; 

/*!
Iterates over the points of the triangulation. 
*/ 
Point_iterator points_begin() const; 

/*!
Past-the-end iterator 
*/ 
Point_iterator points_end() const; 

/// @} 

/*!\name Cell and Facet Circulators 
The following circulators respectively visit all cells or all facets incident to a given edge. They are non-mutable and bidirectional. They are invalidated by any modification of one of the cells traversed.  
*/
/// @{

/*!
Starts at an arbitrary cell incident to `e`. 
\pre `t.dimension() == 3`. 
*/ 
Cell_circulator incident_cells(Edge e) const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Cell_circulator incident_cells(Cell_handle c, int i, int j) const; 

/*!
Starts at cell `start`. 
\pre `t.dimension() == 3` and `start` is incident to `e`. 
*/ 
Cell_circulator incident_cells(Edge e, Cell_handle start) const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Cell_circulator incident_cells(Cell_handle c, int i, int j, 
Cell_handle start) const; 
  /// @}

/*!
\name
The following circulators on facets are defined only in dimension 3, though facets are defined also in dimension 2: there are only two facets sharing an edge in dimension 2.
*/

/// @{
/*!
Starts at an arbitrary facet incident to `e`. 
\pre `t.dimension() == 3` 
*/ 
Facet_circulator incident_facets(Edge e) const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j) const; 

/*!
Starts at facet `start`. 
\pre `start` is incident to `e`. 
*/ 
Facet_circulator incident_facets(Edge e, Facet start) const; 

/*!
Starts at facet of index `f` in `start`. 
*/ 
Facet_circulator incident_facets(Edge e, Cell_handle start, int f) 
const; 

/*!
As above for edge `(i,j)` of `c`. 
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j, 
Facet start) const; 

/*!
As above for edge `(i,j)` of `c` and facet `(start,f)`. 
*/ 
Facet_circulator incident_facets(Cell_handle c, int i, int j, 
Cell_handle start, int f) const; 

/// @} 

/// \name Traversal of the Incident Cells, Facets and Edges, and the Adjacent Vertices of a Given Vertex
/// @{


/*!
Copies the `Cell_handle`s of all cells incident to `v` to the output 
iterator `cells`. 
Returns the resulting output iterator. 
\pre `t.dimension() == 3`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_cells(Vertex_handle v, OutputIterator cells) const; 

/*! 
Try to lock and copy the `Cell_handle`s of all cells incident to `v` into
`cells`. 
Returns `true` in case of success. Otherwise, `cells` is emptied and the function
returns false. In any case, the locked cells are not unlocked by 
`try_lock_and_get_incident_cells()`, leaving this choice to the user.

\pre `t.dimension() == 3`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/
bool
  try_lock_and_get_incident_cells(Vertex_handle v,
                                  std::vector<Cell_handle>& cells) const;
/*!
Copies the `Cell_handle`s of all finite cells incident to `v` to the output 
iterator `cells`. 
Returns the resulting output iterator. 
\pre `t.dimension() == 3`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
finite_incident_cells(Vertex_handle v, OutputIterator cells) const; 

/*!
Copies all `Facet`s incident to `v` to the output iterator 
`facets`. 
Returns the resulting output iterator. 
\pre `t.dimension() > 1`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_facets(Vertex_handle v, OutputIterator facets) const; 

/*!
Copies all finite `Facet`s incident to `v` to the output iterator 
`facets`. 
Returns the resulting output iterator. 
\pre `t.dimension() > 1`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
finite_incident_facets(Vertex_handle v, OutputIterator facets) const; 

/*!
Copies all `Edge`s incident to `v` to the 
output iterator `edges`. Returns the resulting output iterator. 
\pre `t.dimension() > 0`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
incident_edges(Vertex_handle v, OutputIterator edges) const; 

/*!
Copies all finite `Edge`s incident to `v` to the 
output iterator `edges`. Returns the resulting output iterator. 
\pre `t.dimension() > 0`, `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
finite_incident_edges(Vertex_handle v, OutputIterator edges) const; 

/*!
Copies the `Vertex_handle`s of all vertices adjacent to `v` to the 
output iterator `vertices`. If `t.dimension() < 0`, then do 
nothing. Returns the resulting output iterator. 
\pre `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
adjacent_vertices(Vertex_handle v, OutputIterator vertices) const; 

/*!
Copies the `Vertex_handle`s of all finite vertices adjacent to `v` to the 
output iterator `vertices`. If `t.dimension() < 0`, then do 
nothing. Returns the resulting output iterator. 
\pre `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
template <class OutputIterator> 
OutputIterator 
finite_adjacent_vertices(Vertex_handle v, OutputIterator vertices) const; 

/*!
Returns the degree of a vertex, that is, the number of incident vertices. 
The infinite vertex is counted. 
\pre `v != Vertex_handle()`, `t.is_vertex(v)`. 
*/ 
size_type degree(Vertex_handle v) const; 

/// @} 

/// \name Traversal Between Adjacent Cells 
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
/// The responsibility of keeping a valid triangulation belongs to the
/// user when using advanced operations allowing a direct manipulation
/// of cells and vertices. We provide the user with the following
/// methods to help debugging.
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the combinatorial validity of the triangulation. Checks also the 
validity of its geometric embedding (see 
Section \ref Triangulation3secintro). 
When `verbose` is set to `true`, 
messages describing the first invalidity encountered are printed. 
\cgalDebugEnd
*/ 
bool 
is_valid(bool verbose = false) const; 

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the combinatorial validity of the cell by calling the 
`is_valid` method of the `TriangulationDataStructure_3` cell class. Also checks the 
geometric validity of `c`, if `c` is finite. (See 
Section \ref Triangulation3secintro.) 

When `verbose` is set to `true`, messages are printed to give 
a precise indication of the kind of invalidity encountered. 
\cgalDebugEnd
*/ 
bool 
is_valid(Cell_handle c, bool verbose = false) const; 

/// @}

/*!
\name I/O

\cgal provides an interface to Geomview for a 3D-triangulation, 
see Chapter \ref Chapter_Geomview on `Geomview_stream`. 

\code
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
\endcode

The information in the `iostream` is: the dimension, the number of 
finite vertices, the non-combinatorial information about vertices 
(point, etc; note that the infinite vertex is numbered 0), the number 
of cells, the indices of the vertices of each cell, plus the 
non-combinatorial information about each cell, then the indices of the 
neighbors of each cell, where the index corresponds to the preceding 
list of cells. When dimension \f$ <\f$ 3, the same information is stored 
for faces of maximal dimension instead of cells. 
*/
/// @{

/*!
Reads the underlying combinatorial triangulation from `is` by 
calling the corresponding input operator of the triangulation data 
structure class (note that the infinite vertex is numbered 0), and the 
non-combinatorial information by calling the corresponding input 
operators of the vertex and the cell classes (such as point 
coordinates), which are provided by overloading the stream operators 
of the vertex and cell types. Assigns the resulting triangulation to 
`t`. 
*/ 
istream& operator>> (istream& is, Triangulation_3 &t); 

/*!
Writes the triangulation `t` into `os`. 
*/ 
ostream& operator<< (ostream& os, const Triangulation_3 &t); 

/// @}

/// @} 

/// \name Concurrency 
/// @{

/*! 
Set the pointer to the lock data structure.
*/ 
void set_lock_data_structure(Lock_data_structure *lock_ds) const;

/// @} 

}; /* end Triangulation_3 */
} /* end namespace CGAL */
