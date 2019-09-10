
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3MainClasses

The class `Periodic_3_triangulation_3` represents a 3-dimensional
triangulation of a point set in \f$ \mathbb T_c^3\f$.

\tparam Traits must be a model of the concept `Periodic_3TriangulationTraits_3`.

\tparam TDS must be a model of the concept `TriangulationDataStructure_3`
with some additional functionality in cells and vertices.
Its default value is
`Triangulation_data_structure_3<Triangulation_vertex_base_3<PT,Periodic_3_triangulation_ds_vertex_base_3<>>,Triangulation_cell_base_3<PT,Periodic_3_triangulation_ds_cell_base_3<>>>`.

\sa `Periodic_3_Delaunay_triangulation_3`
\sa `Periodic_3_regular_triangulation_3`
*/
template< typename Traits, typename TDS >
class Periodic_3_triangulation_3 {
public:

/// \name Types

/// @{

/*!

*/
typedef Traits Geom_traits;

/*!

*/
typedef TDS Triangulation_data_structure;

/*!

*/
typedef Geom_traits::Periodic_3_offset_3
Offset;

/*!
A type representing an axis-aligned cuboid. It must be a model of `Traits::Iso_cuboid_3`. Used to represent the original domain.
*/
typedef Geom_traits::Iso_cuboid_3 Iso_cuboid;

/*!
Integer triple to
store the number of sheets in each direction of space.
*/
typedef array<int,3> Covering_sheets;

/*!
The point type of the triangulation.

\note This type is equal to `Geom_traits::Point_3` when considering periodic
Delaunay triangulations and to `Geom_traits::Weighted_point_3` when
considering periodic weighted Delaunay triangulations.
*/
typedef TDS::Vertex::Point Point;

/*!
The geometric basic 3D point type.
*/
typedef Geom_traits::Point_3 Point_3;

/*!

*/
typedef Geom_traits::Segment_3 Segment;

/*!

*/
typedef Geom_traits::Triangle_3 Triangle;

/*!

*/
typedef Geom_traits::Tetrahedron_3 Tetrahedron;

/*!
Represents a point-offset pair. The point in the pair lies in the original domain.
Note that the inner type is `Point`.
*/
typedef std::pair< Point, Offset > Periodic_point;

/*!
Represents a point-offset pair. The point in the pair lies in the original domain.
Note that the inner type is `Point_3`.
*/
typedef std::pair< Point_3, Offset > Periodic_point_3;

/*!
A periodic segment. Note that the inner type is `Periodic_point`.
*/
typedef array< Periodic_point, 2> Periodic_segment;

/*!
A periodic segment. Note that the inner type is `Periodic_point_3`.
*/
typedef array< Periodic_point_3, 2> Periodic_segment_3;

/*!
A periodic triangle. Note that the inner type is `Periodic_point`.
*/
typedef array< Periodic_point, 3> Periodic_triangle;

/*!
A periodic triangle. Note that the inner type is `Periodic_point_3`.
*/
typedef array< Periodic_point_3, 3> Periodic_triangle_3;

/*!
A periodic tetrahedron. Note that the inner type is `Periodic_point`.
*/
typedef array< Periodic_point, 4> Periodic_tetrahedron;

/*!
A periodic tetrahedron. Note that the inner type is `Periodic_point_3`.
*/
typedef array< Periodic_point_3, 4> Periodic_tetrahedron_3;

/// @}

/*! \name
Only vertices (\f$ 0\f$-faces) and cells (\f$ 3\f$-faces) are stored.
Edges (\f$ 1\f$-faces) and facets (\f$ 2\f$-faces) are not explicitly represented
and thus there are no corresponding classes (see Section \ref P3Triangulation3secintro).
*/
/// @{


/*!

*/
typedef Triangulation_data_structure::Vertex Vertex;

/*!

*/
typedef Triangulation_data_structure::Cell Cell;

/*!

*/
typedef Triangulation_data_structure::Edge Edge;

/*!

*/
typedef Triangulation_data_structure::Facet Facet;



/// @}

/*! \name

The vertices and faces of the triangulations are accessed through
`handles`, `iterators` and `circulators`. A handle is a type which
supports the two dereference operators and `operator->`. The Handle
concept is documented in the support library. Iterators and
circulators are bidirectional and non-mutable. The edges and facets of
the triangulation can also be visited through iterators and
circulators which are bidirectional and non-mutable.

Iterators and circulators are convertible to the corresponding
handles, thus the user can pass them directly as arguments to the
functions.

*/
/// @{

/*!
handle to a vertex
*/
typedef Triangulation_data_structure::Vertex_handle Vertex_handle;

/*!
handle to a cell
*/
typedef Triangulation_data_structure::Cell_handle Cell_handle;

/*!
Size type (an unsigned integral type)
*/
typedef Triangulation_data_structure::size_type size_type;

/*!
Difference type (a signed integral type)
*/
typedef Triangulation_data_structure::difference_type difference_type;

/*!
iterator over cells
*/
typedef Triangulation_data_structure::Cell_iterator Cell_iterator;

/*!
iterator over facets
*/
typedef Triangulation_data_structure::Facet_iterator Facet_iterator;

/*!
iterator over edges
*/
typedef Triangulation_data_structure::Edge_iterator Edge_iterator;

/*!
iterator over vertices
*/
typedef Triangulation_data_structure::Vertex_iterator Vertex_iterator;

/*!
iterator over the vertices whose
corresponding points lie in the original domain, i.e.\ for each set
of periodic copies the `Unique_vertex_iterator` iterates over
exactly one representative.
*/
typedef unspecified_type Unique_vertex_iterator;

/*!
circulator over all cells incident to a given edge
*/
typedef Triangulation_data_structure::Cell_circulator Cell_circulator;

/*!
circulator over all facets incident to a given edge
*/
typedef Triangulation_data_structure::Facet_circulator Facet_circulator;

/// @}

/// \name Geometric Iterators:
/// The following iterators have value type `Periodic_segment`, `Periodic_triangle`,
/// and `Periodic_tetrahedron`, which have inner type `Point`.
/// @{

/*!
iterator over the tetrahedra
corresponding to cells of the triangulation.
*/
typedef unspecified_type Periodic_tetrahedron_iterator;

/*!
iterator over the triangles
corresponding to facets of the triangulation.
*/
typedef unspecified_type Periodic_triangle_iterator;

/*!
iterator over the segments
corresponding to edges of the triangulation.
*/
typedef unspecified_type Periodic_segment_iterator;

/*!
iterator over the points
corresponding to vertices of the triangulation.
*/
typedef unspecified_type Periodic_point_iterator;

/// @}

/// \name Enums

/// @{

/*!
The enum `Locate_type` is defined by `Periodic_3_triangulation_3` to
specify which case occurs when locating a point in the
triangulation. If the triangulation does not contain any points
`EMPTY` is returned.

\sa `CGAL::Periodic_3_triangulation_3`

*/
enum Locate_type {VERTEX=0, EDGE, FACET, CELL, EMPTY};

/*!
The enum `Iterator_type` is defined by `Periodic_3_triangulation_3` to
specify the behavior of geometric iterators.

The elements of the enum have the following meaning:


\sa `CGAL::Periodic_3_triangulation_3`
*/
  enum Iterator_type {STORED=0, /*!< Return all geometric primitives as they are stored internally
                                  in `Triangulation_data_structure_3`. */
                      UNIQUE,  /*!< Return only one representative of each geometric
                                 primitive even if the triangulation is computed in a multiply
                                 sheeted covering space. Choose the representative whose maximum
                                 offset is minimal but non-negative in each direction of space.*/
                      STORED_COVER_DOMAIN,  /*!< Same as STORED but return additionally
                                              all primitives whose intersection with the original domain of the
                                              current covering space is non-empty. */
                      UNIQUE_COVER_DOMAIN  /*!< Same as UNIQUE but return additionally
                                             all primitives whose intersection with the original domain is
                                             non-empty. */
  };

/// @}

/// \name Creation
/// @{

/*!
Introduces an empty triangulation `t` with `domain` as
original domain.
\pre `domain` is a cube.
*/
Periodic_3_triangulation_3(const Iso_cuboid & domain = Iso_cuboid(0,0,0,1,1,1),
                           const Geom_traits & traits = Geom_traits());

/*!
Copy constructor. All vertices and faces are duplicated.
*/
Periodic_3_triangulation_3 (const Periodic_3_triangulation_3 & tr);

/// @}

/// \name Assignment
/// @{

/*!
The triangulation `tr` is duplicated, and modifying the copy after the
duplication does not modify the original. The previous triangulation held
by `t` is deleted.
*/
Periodic_3_triangulation_3 & operator=(const Periodic_3_triangulation_3 & tr);

/*!
The triangulations `tr` and `t` are swapped.
`t`.`swap(tr)` should be preferred to `t` = `tr` or to
`t(tr)` if `tr` is deleted after that. Indeed, there is no
copy of cells and vertices, thus this method runs in constant time.
*/
void swap(Periodic_3_triangulation_3 & tr);

/*!
Deletes all vertices and all cells of `t`.
*/
void clear();

/*!
Equality operator. Returns true iff there exist a bijection between the
vertices of `t1` and those of `t2` and a bijection between the cells of
`t1` and those of `t2`, which preserve the geometry of the
triangulation, that is, the points of each corresponding pair of vertices are
equal, and the tetrahedra corresponding to each pair of cells are equal (up to
a permutation of their vertices).
\relates Periodic_3_triangulation_3
*/
template < class Traits, class TDS1, class TDS2 >
bool operator==(const Periodic_3_triangulation_3<Traits, TDS1> & t1, const Periodic_3_triangulation_3<Traits, TDS2> & t2);

/*!
The opposite of `operator==`.
\relates Periodic_3_triangulation_3
*/
template < class Traits, class TDS1, class TDS2 >
bool operator!=(const Periodic_3_triangulation_3<Traits, TDS1> & t1, const Periodic_3_triangulation_3<Traits, TDS2> & t2);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a const reference to the geometric traits object.
*/
const Geom_traits & geom_traits() const;

/*!
Returns a const reference to the triangulation data structure.
*/
const Triangulation_data_structure & tds() const;

/*!
Returns the original domain.
*/
Iso_cuboid domain() const;

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Returns the number of sheets of the covering space the triangulation is
currently computed in.
\cgalAdvancedEnd
*/
Covering_sheets number_of_sheets() const;

/*!
Get the offset between the origins of the internal offset coordinate
systems of two neighboring cells with respect from ch to its i-th neighbor.
*/
Offset neighbor_offset(Cell_handle ch, int i) const;

/// @}

/// \name Non const access
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Changes the domain. Note that this function calls `clear()`,
i.e., it erases the existing triangulation.
\cgalAdvancedEnd
*/
void set_domain(const Iso_cuboid dom);

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Returns a reference to the triangulation data structure.

The responsibility of keeping a valid triangulation belongs
to the user when using advanced operations allowing a direct
manipulation of the `tds`. This method is mainly a help for users
implementing their own triangulation algorithms.
\cgalAdvancedEnd
*/
Triangulation_data_structure & tds();

/// @}

/// \name Non-constant-time queries and conversions
/// @{

/*!
The current triangulation remains a triangulation in the 1-sheeted
covering space even after adding points if this method returns
`true`. This test relies on a heuristic, i.e.\ if it answers
`false` nothing is known. This test runs in constant-time when
not computing in the 1-sheeted covering space. (This test uses the length
of the longest edge in the triangulation as a
criterion \cgalCite{cgal:ct-c3pt-09}.)
*/
bool is_extensible_triangulation_in_1_sheet_h1() const;

/*!
The same as `is_extensible_triangulation_in_1_sheet_h1()` but with
a more precise heuristic, i.e.\ it might answer `true` in cases in which
`is_extensible_triangulation_in_1_sheet_h1()` would not. However, it is
much less time efficient when not computing in the 1-sheeted covering
space. (This test uses the diameter of the largest empty ball in the
input point set as a criterion \cgalCite{cgal:ct-c3pt-09}.)
*/
bool is_extensible_triangulation_in_1_sheet_h2() const;

/*!
Returns `true` if the current triangulation would still be a
triangulation in the 1-sheeted covering space, returns `false` otherwise.
*/
bool is_triangulation_in_1_sheet() const;

/*!
Converts the current triangulation into the same periodic
triangulation in the 1-sheeted covering space.

\attention It is not recommended to interfere with the built-in
covering management. Especially a premature conversion to the
1-sheeted covering space might lead to problems when modifying the
triangulation later.
*/
void convert_to_1_sheeted_covering() const;

/*!
Converts the current triangulation into the same periodic
triangulation in the 27-sheeted covering space.

\attention It is not recommended to interfere with the built-in
covering management. Especially a premature conversion to the
1-sheeted covering space might lead to problems when modifying the
triangulation later.
*/
void convert_to_27_sheeted_covering() const;

/// @}

/// \name Constant-time access functions
/// @{

/*!
Returns the number of vertices. Counts all vertices that are
representatives of the same point in \f$ \mathbb T_c^3\f$ as one vertex.
*/
size_type number_of_vertices() const;

/*!
Returns the number of cells. Counts all cells that are
representatives of the same tetrahedron in \f$ \mathbb T_c^3\f$ as one
cell.
*/
size_type number_of_cells() const;

/*!
Returns the number of vertices in the data structure. This is the
same as the number of sheets times `number_of_vertices()`.
*/
size_type number_of_stored_vertices() const;

/*!
Returns the number of cells in the data structure. This is the same
as the number of sheets times `number_of_cells()`.
*/
size_type number_of_stored_cells() const;

/// @}

/// \name Non-constant-time access functions
/// @{

/*!
Returns the number of edges. Counts all edges that are
representatives of the same segment in \f$ \mathbb T_c^3\f$ as one edge.
*/
size_type number_of_edges() const;

/*!
Returns the number of facets. Counts all facets that are
representatives of the same triangle in \f$ \mathbb T_c^3\f$ as one
facet.
*/
size_type number_of_facets() const;

/*!
Returns the number of edges in the data structure. This is the same
as the number of sheets times `number_of_edges()`.
*/
size_type number_of_stored_edges() const;

/*!
Returns the number of facets in the data structure. This is the same
as the number of sheets times `number_of_facets()`.
*/
size_type number_of_stored_facets() const;

/// @}

/// \name Geometric Access Functions
/// The following functions return object of types `Periodic_segment`,
/// `Periodic_triangle`, and `Periodic_tetrahedron`, which have inner type `Point`.
/// @{

/*!
Returns the periodic point given by vertex `v`. If `t` is
represented in the 1-sheeted covering space, the offset is always
zero. Otherwise `v` can correspond to a periodic copy outside
`domain` of an input point.
*/
Periodic_point periodic_point(const Vertex_handle v) const;

/*!
If `t` is represented in the 1-sheeted covering space, this function
returns the periodic point given by the \f$ i\f$-th vertex of cell
`c`, that is the point in the original domain and the offset of
the vertex in `c`.
If `t` is represented in the 27-sheeted covering space,
this offset is possibly added to another offset determining the periodic copy.
\pre \f$ i \in\{0,1,2,3\}\f$
*/
Periodic_point periodic_point(const Cell_handle c, int i)
const;

/*!
Returns the periodic segment formed by the two point-offset pairs
corresponding to the two vertices of edge `(c,i,j)`.
\pre \f$ i,j \in\{0,1,2,3\}\f$, \f$ i\neq j\f$
*/
Periodic_segment periodic_segment(const Cell_handle c, int
i, int j) const;

/*!
Returns the periodic segment formed by the two point-offset pairs
corresponding to the two vertices of edge `(c,i,j)`.

A translation in accordance with `offset` is applied on the point-offet pairs.
\pre \f$ i,j \in\{0,1,2,3\}\f$, \f$ i\neq j\f$
*/
Periodic_segment periodic_segment(const Cell_handle c, Offset offset,  int
i, int j) const;

/*!
Same as the previous method for edge `e`.
*/
Periodic_segment periodic_segment(const Edge & e) const;

/*!
Returns the periodic triangle formed by the three point-offset pairs
corresponding to the three vertices of facet
`(c,i)`. The triangle is oriented so that its normal points to the
inside of cell `c`.
\pre \f$ i \in\{0,1,2,3\}\f$
*/
Periodic_triangle periodic_triangle(const Cell_handle c, int
i) const;

/*!
Returns the periodic triangle formed by the three point-offset pairs
corresponding to the three vertices of facet
`(c,i)`.

A translation in accordance with `offset` is applied on the point-offet pairs.

The triangle is oriented so that its normal points to the
inside of cell `c`.
\pre \f$ i \in\{0,1,2,3\}\f$
*/
Periodic_triangle periodic_triangle(const Cell_handle c, Offset offset, int i) const;

/*!
Same as the previous method for facet `f`.
*/
Periodic_triangle periodic_triangle(const Facet & f) const;

/*!
Returns the periodic tetrahedron formed by the four point-offset pairs
corresponding to the four vertices of `c`.
*/
Periodic_tetrahedron periodic_tetrahedron(const Cell_handle c) const;

/*!
Returns the periodic tetrahedron formed by the four point-offset pairs
corresponding to the four vertices of `c`.

A translation in accordance with `offset` is applied on the point-offet pairs.
*/
Periodic_tetrahedron periodic_tetrahedron(const Cell_handle c, Offset offset) const;

/// \name
/// \warning The following functions were renamed with %CGAL 4.11 to clarify
/// that they return geometric objects with inner type `Point_3`.
///
/// Note that a traits class providing exact constructions should be
/// used in order to guarantee the following operations to be exact
/// (as opposed to computing the triangulation only, which requires
/// only exact predicates).
///
/// @{

/*!
Converts the periodic point of type `PP` to a `%Point_3`.
The type PP can be either `Periodic_point` or `Periodic_point_3`.
*/
template<typename PP>
Point_3 construct_point(const PP & pp) const;

/*!
Converts the `Point` `p` to a `%Point_3`.
*/
Point_3 construct_point(const Point & p) const;

/*!
Same as above, with offsets.
*/
template<typename P>
Point_3 construct_point(const P& p1, const Offset& o1) const;

/*!
Converts the periodic segment of type `PS` to a `Segment`. The type `PS` can be
either `Periodic_segment` or `Periodic_segment_3`.
*/
template<typename PS>
Segment construct_segment(const PS & s) const;

/*!
Creates a segment from two points. The type `P` can be either `Point` or `Point_3`.
*/
template<typename P>
Segment construct_segment(const P& p1, const P& p2) const;

/*!
Same as above, with offsets.
*/
template<typename P>
Segment construct_segment(const P& p1, const P& p2, const Offset& o1, const Offset& o2) const;

/*!
Converts the periodic triangle of type `PT` to a `Triangle`. The type `PT` can
be either `Periodic_triangle` or `Periodic_triangle_3`.
*/
template<typename PT>
Triangle construct_triangle(const PT & t) const;

/*!
Creates a triangle from three points. The type `P` can be either `Point` or `Point_3`.
*/
template<typename P>
Triangle construct_triangle(const P& p1, const P& p2, const P& p3) const;

/*!
Same as above, with offsets.
*/
template<typename P>
Triangle construct_triangle(const P& p1, const P& p2, const P& p3,
                            const Offset& o1, const Offset& o2, const Offset& o3) const;

/*!
Converts the periodic tetrahedron of type `PT` to a `Tetrahedron`. The type `PT`
can be either `Periodic_tetrahedron` or `Periodic_tetrahedron_3`.
*/
template<typename PT>
Tetrahedron construct_tetrahedron(const PT & t) const;

/*!
Creates a tetrahedron from four points. The type `P` can be either `Point` or `Point_3`.
*/
template<typename P>
Tetrahedron construct_tetrahedron(const P& p1, const P& p2, const P& p3, const P& p4) const;

/*!
Same as above, with offsets.
*/
template<typename P>
Tetrahedron construct_tetrahedron(const P& p1, const P& p2, const P& p3, const P& p4,
                                  const Offset& o1, const Offset& o2, const Offset& o3, const Offset& o4) const;

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
Tests whether `((u,offu),(v,offu))` is an edge of
`t`. If the edge is found, it gives a cell `c` having this
edge and the indices `i` and `j` of the vertices `u` and
`v` in `c`, in this order.
\pre `u` and `v` are vertices of `t`.
*/
bool is_edge(Vertex_handle u, const Offset & offu,
Vertex_handle v, const Offset & offv,
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
Tests whether `((u,offu),(v,offv),(w,offw))`
is a facet of `t`. If the facet is found,
it computes a cell `c` having this facet and the indices `i`,
`j` and `k` of the vertices `u`, `v` and `w` in `c`,
in this order.
\pre `u`, `v` and `w` are vertices of `t`.
*/
bool is_facet(Vertex_handle u, const Offset & offu,
Vertex_handle v, const Offset & offv,
Vertex_handle w, const Offset & offw,
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

/*!
Tests whether
`((u,offu),(v,offv),(w,offv),(x,offx))` is a cell of `t`.
If the cell `c` is found, the method
computes the indices `i`, `j`, `k` and `l` of the
vertices `u`, `v`, `w` and `x` in `c`, in this
order.
\pre `u`, `v`, `w` and `x` are vertices of `t`.
*/
bool is_cell(Vertex_handle u, const Offset & offu,
Vertex_handle v, const Offset & offv,
Vertex_handle w, const Offset & offw,
Vertex_handle x, const Offset & offx,
Cell_handle & c,
int & i, int & j, int & k, int & l) const;

/*!
Tests whether
`((u,offu),(v,offv),(w,offv),(x,offx))` is a
cell of `t` and computes this cell `c`.
\pre `u`, `v`, `w` and `x` are vertices of `t`.
*/
bool is_cell(Vertex_handle u, const Offset & offu,
Vertex_handle v, const Offset & offv,
Vertex_handle w, const Offset & offw,
Vertex_handle x, const Offset & offx,
Cell_handle & c) const;

/// @}

/// \name
/// There is a method `has_vertex` in the cell class. The analogous methods for facets are defined here.
/// @{

/*!
If `v` is a vertex of `f`, then `j` is the index of
`v` in the cell `f.first`, and the method returns `true`.
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

/// \name
/// The following three methods test whether two facets have the same vertices.
/// @{

/*!

*/
bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const;

/*!

*/
bool are_equal(const Facet & f, const Facet & g) const;

/*!

*/
bool are_equal(const Facet & f, Cell_handle n, int j) const;

/// @}

/// \name Point Location
/// The class `Periodic_3_triangulation_3` provides three functions to
/// locate a given point with respect to a triangulation. It provides
/// also functions to test if a given point is inside a face or
/// not. Note that the class `Periodic_3_Delaunay_triangulation_3`
/// also provides a `nearest_vertex()` function.
/// @{

/*!

Returns the cell that contains the query in its interior. If
`query` lies on a facet, an edge or on a vertex, one of the cells
having `query` on its boundary is returned.

The optional argument `start` is used as a starting place for the
search.
\pre `query` lies in the original domain `domain`.

*/
Cell_handle
locate(const Point & query, Cell_handle start = Cell_handle()) const;

/*!
Returns the cell that contains the query in its interior. If
`query` lies on a facet, an edge or on a vertex, one of the cells
having `query` on its boundary is returned.

`locate_offset` is the offset that must be used by the function
periodic_tetrahedron() together with the returned cell, so that the
constructed Periodic_tetrahedron contains the query point.

The optional argument `start` is used as a starting place for the
search.
\pre `query` lies in the original domain `domain`.
*/
Cell_handle
locate(const Point & query, Offset& locate_offset, Cell_handle start = Cell_handle()) const;

/*!
Same as `locate()` but uses inexact predicates.

This function returns a handle on a cell that is a good approximation of the exact
location of `query`, while being faster. Note that it may return a handle on a cell
whose interior does not contain `query`.

Note that this function is available only if the %Cartesian coordinates of `query`
are accessible with functions `x()`, `y()` and `z()`.
*/
Cell_handle
inexact_locate(const Point & query, Cell_handle start = Cell_handle()) const;


/*!
The \f$ k\f$-face that contains `query` in its interior is
returned, by means of the cell returned together with `lt`, which
is set to the locate type of the query (`VERTEX, EDGE, FACET, CELL`) and two indices `li` and `lj` that
specify the \f$ k\f$-face of the cell containing `query`.

If the \f$ k\f$-face is a cell, `li` and `lj` have no
meaning; if it is a facet (resp. vertex), `li` gives the index of
the facet (resp. vertex) and `lj` has no meaning; if it is an
edge, `li` and `lj` give the indices of its vertices.

If there is no vertex in the triangulation yet, `lt` is set to
`EMPTY` and `locate` returns the default constructed handle.

The optional argument `start` is used as a starting place for the
search.
\pre `query` lies in the original domain `domain`.

*/
Cell_handle
locate(const Point & query, Locate_type & lt,
int & li, int & lj, Cell_handle start = Cell_handle() ) const;

/*!
The \f$ k\f$-face that contains `query` in its interior is
returned, by means of the cell returned together with `lt`, which
is set to the locate type of the query (`VERTEX, EDGE, FACET, CELL`) and two indices `li` and `lj` that
specify the \f$ k\f$-face of the cell containing `query`.

If the \f$ k\f$-face is a cell, `li` and `lj` have no
meaning; if it is a facet (resp. vertex), `li` gives the index of
the facet (resp. vertex) and `lj` has no meaning; if it is an
edge, `li` and `lj` give the indices of its vertices.

If there is no vertex in the triangulation yet, `lt` is set to
`EMPTY` and `locate` returns the default constructed handle.

`locate_offset` is the offset that must be used by the function
periodic_tetrahedron() together with the returned cell, so that the
constructed Periodic_tetrahedron contains the query point.

The optional argument `start` is used as a starting place for the
search.
\pre `query` lies in the original domain `domain`.
*/
Cell_handle
locate(const Point & query, Offset& locate_offset, Locate_type & lt,
int & li, int & lj, Cell_handle start = Cell_handle() ) const;

/*!
Returns a value indicating on which side of the oriented boundary
of `c` the point `p` lies. More precisely, it returns:

- `ON_BOUNDED_SIDE` if `p` is inside the cell.

- `ON_BOUNDARY` if `p` on the boundary of the cell. Then
`lt` together with `li` and `lj` give the precise location
on the boundary. (See the descriptions of the `locate` methods.)

- `ON_UNBOUNDED_SIDE` if `p` lies outside the cell.
\pre `query` lies in the original domain `domain`.
*/
Bounded_side
side_of_cell(const Point & p,
Cell_handle c,
Locate_type & lt, int & li, int & lj) const;

/// @}

/*!
\name Traversal of the Triangulation

The periodic triangulation class provides several iterators and circulators
that allow one to traverse it.

\name Cell, Face, Edge and Vertex Iterators

The following iterators allow the user to visit cells, facets, edges
and vertices of the stored triangulation, i.e.\ in case of computing in
a multiply sheeted covering space all stored periodic copies of each
item are returned. These iterators are non-mutable, bidirectional and
their value types are respectively `Cell`, `Facet`, `Edge` and
`Vertex`. They are all invalidated by any change in the triangulation.

*/

/// @{

/*!
Starts at an arbitrary vertex. Iterates over all vertices. Returns
`vertices_end()` if `t`.`number_of_vertices()` \f$ =0\f$.
*/
Vertex_iterator vertices_begin() const;

/*!
Past-the-end iterator
*/
Vertex_iterator vertices_end() const;

/*!
Starts at an arbitrary edge. Iterates over all edges. Returns
`edges_end()` if `t`.`number_of_vertices()` \f$ =0\f$.
*/
Edge_iterator edges_begin() const;

/*!
Past-the-end iterator
*/
Edge_iterator edges_end() const;

/*!
Starts at an arbitrary facet. Iterates over all facets. Returns
`facets_end()` if `t`.`number_of_vertices()` \f$ =0\f$.
*/
Facet_iterator facets_begin() const;

/*!
Past-the-end iterator
*/
Facet_iterator facets_end() const;

/*!
Starts at an arbitrary cell. Iterates over all cells. Returns
`cells_end()` if `t`.`number_of_vertices()` \f$ =0\f$.
*/
Cell_iterator cells_begin() const;

/*!
Past-the-end iterator
*/
Cell_iterator cells_end() const;

/*!
Starts at an arbitrary vertex. Iterates over all vertices whose
corresponding points lie in the original domain, i.e.\ for each set
of periodic copies the `Unique_vertex_iterator` iterates over
exactly one representative. Returns `unique_vertices_end()` if
`t`.`number_of_vertices()` \f$ =0\f$.
*/
Unique_vertex_iterator unique_vertices_begin() const;

/*!
Past-the-end iterator
*/
Unique_vertex_iterator unique_vertices_end() const;

/// @}

/*! \name Geometric Iterators

The following iterators allow the user to obtain geometric primitives
corresponding to cells, facets, edges, and vertices of the
triangulation. These iterators are non-mutable, bidirectional and
their value types are respectively `Periodic_point`,
`Periodic_segment`, `Periodic_triangle`, and
`Periodic_tetrahedron`. They are all invalidated by any change in the
triangulation. If the periodic triangulation is not computed in the
1-sheeted covering space, these iterators can be used to retain only
the geometric primitives in the original domain. This can be
controlled using the enum `Iterator_type`.

\anchor P3Triangulation3figgeom_iterators
<img border=0 src="./it_STORED_small.jpg" align=middle alt="STORED">
<img border=0 src="./it_STORED_COVER_DOMAIN_small.jpg" align=middle alt="STORED_COVER_DOMAIN">
<img border=0 src="./it_UNIQUE_small.jpg" align=middle alt="UNIQUE">
<img border=0 src="./it_UNIQUE_COVER_DOMAIN_small.jpg" align=middle alt="UNIQUE_COVER_DOMAIN">

<CENTER><b> The four different modes of the geometric iterators: STORED, STORED_COVER_DOMAIN, UNIQUE, UNIQUE_COVER_DOMAIN. Note that in case of computing in the 1-sheeted covering space, STORED and UNIQUE give the same result. </b></center>

*/

/// @{

/*!
Iterates over the points of the triangulation. Its behavior is
defined by the `Iterator_type` `it` as described on
\ref CGAL::Periodic_3_triangulation_3::Iterator_type.
*/
Periodic_point_iterator periodic_points_begin(Iterator_type it =
STORED) const;

/*!
Past-the-end iterator. Note that to match another
`Periodic_point_iterator` both must have the same
`Iterator_type` `it`.
*/
Periodic_point_iterator periodic_points_end(Iterator_type it =
STORED) const;

/*!
Iterates over the segments of the triangulation. Its behavior is
defined by the `Iterator_type` `it` as described on
\ref CGAL::Periodic_3_triangulation_3::Iterator_type.
*/
Periodic_segment_iterator periodic_segments_begin(Iterator_type it =
STORED) const;

/*!
Past-the-end iterator. Note that to match another
`Periodic_segment_iterator` both must have the same
`Iterator_type` `it`.
*/
Periodic_segment_iterator periodic_segments_end(Iterator_type it =
STORED) const;

/*!
Iterates over the triangles of the triangulation. Its behavior is
defined by the `Iterator_type` `it` as described on
\ref CGAL::Periodic_3_triangulation_3::Iterator_type.
*/
Periodic_triangle_iterator periodic_triangles_begin(Iterator_type it =
STORED) const;

/*!
Past-the-end iterator. Note that to match another
`Periodic_triangle_iterator` both must have the same
`Iterator_type` `it`.
*/
Periodic_triangle_iterator periodic_triangles_end(Iterator_type it =
STORED) const;

/*!
Iterates over the tetrahedra of the triangulation. Its behavior is
defined by the `Iterator_type` `it` as described on
\ref CGAL::Periodic_3_triangulation_3::Iterator_type.
*/
Periodic_tetrahedron_iterator periodic_tetrahedra_begin(Iterator_type it =
STORED) const;

/*!
Past-the-end iterator. Note that to match another
`Periodic_tetrahedron_iterator` both must have the same
`Iterator_type` `it`.
*/
Periodic_tetrahedron_iterator periodic_tetrahedra_end(Iterator_type it =
STORED) const;

/// @}

/// \name Cell and Facet Circulators
/// The following circulators respectively visit all cells or all
/// facets incident to a given edge. They are non-mutable and
/// bidirectional. They are invalidated by any modification of one of
/// the cells traversed.
/// @{

/*!
Starts at an arbitrary cell incident to `e`.
*/
Cell_circulator incident_cells(Edge e) const;

/*!
As above for edge `(i,j)` of `c`.
*/
Cell_circulator incident_cells(Cell_handle c, int i, int j) const;

/*!
Starts at cell `start`.
\pre `start` is incident to `e`.
*/
Cell_circulator incident_cells(Edge e, Cell_handle start) const;

/*!
As above for edge `(i,j)` of `c`.
*/
Cell_circulator incident_cells(Cell_handle c, int i, int j,
Cell_handle start) const;

/*!
Starts at an arbitrary facet incident to `e`.
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

/// \name Traversal of the Incident Cells and Facets, and the Adjacent Vertices of a Given Vertex
/// @{

/*!
Copies the `Cell_handle`s of all cells incident to `v` to the output
iterator `cells`. Returns the resulting output iterator.
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `t`.`is_vertex(v)`.
*/
template <class OutputIterator>
OutputIterator
incident_cells(Vertex_handle v, OutputIterator cells) const;

/*!
Copies the `Facet`s incident to `v` to the output iterator
`facets`.
Returns the resulting output iterator.
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `t`.`is_vertex(v)`.
*/
template <class OutputIterator>
OutputIterator
incident_facets(Vertex_handle v, OutputIterator facets) const;

/*!
Copies the `Edge`s incident to `v` to the output iterator
`edges`.
Returns the resulting output iterator.
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `t`.`is_vertex(v)`.
*/
template <class OutputIterator>
OutputIterator
incident_edges(Vertex_handle v, OutputIterator edges) const;

/*!
Copies the `Vertex_handle`s of all vertices adjacent to `v` to the
output iterator `vertices`. Returns the resulting output iterator.
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `t`.`is_vertex(v)`.
*/
template <class OutputIterator>
OutputIterator
adjacent_vertices(Vertex_handle v, OutputIterator vertices) const;

/*!
Returns the degree of a vertex, that is, the number of adjacent vertices.
\pre `v` \f$ \neq\f$ `Vertex_handle()`, `t`.`is_vertex(v)`.
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
Returns the same facet viewed from the other adjacent cell.
*/
Facet mirror_facet(Facet f) const;

/// @}

/// \name Checking
/// \cgalAdvancedBegin
/// The responsibility of keeping a valid triangulation belongs to the
/// user when using advanced operations allowing a direct manipulation
/// of cells and vertices. We provide the user with the following
/// methods to help debugging.
/// \cgalAdvancedEnd
/// @{

/*!
Checks the combinatorial validity of the triangulation. Checks also the
validity of its geometric embedding (see
Section \ref P3Triangulation3secintro).
When `verbose`
is set to true, messages describing the first invalidity encountered
are printed.
*/
bool
is_valid(bool verbose = false) const;

/*!
Checks the combinatorial validity of the cell by calling the
`is_valid` method of the `Triangulation_data_structure` cell
class. Also checks the geometric validity of `c`, if `c` is
finite. (See Section \ref P3Triangulation3secintro.)

When `verbose` is set to `true`, messages are printed to give
a precise indication of the kind of invalidity encountered.
*/
bool
is_valid(Cell_handle c, bool verbose = false) const;

/// @}

}; /* end Periodic_3_triangulation_3 */


/*!
Reads a triangulation from `is` and stores it in `t`.
\pre `is` has the below described format.

The information in the `iostream` is:
<UL>
<LI>the original domain
<LI>the number of sheets of the covering space as in
`number_of_sheets()`
<LI>the number of vertices
<LI>the non-combinatorial information of vertices (point
resp. point-offset pairs, etc.)
<LI>the number of cells
<LI>the indices of the vertices of each cell
<LI>the indices of the neighbors of each cell, where the index
corresponds to the preceding list of cells
<LI>the offsets corresponding to the vertices of the cells
<LI>the non-combinatorial information of each cell
</UL>

\relates Periodic_3_triangulation_3
*/
istream& operator>> (istream& is, Periodic_3_triangulation_3 &t);

/*!
Writes the triangulation `t` into `os`.

The information in the `iostream` is:
<UL>
<LI>the original domain
<LI>the number of sheets of the covering space as in
`number_of_sheets()`
<LI>the number of vertices
<LI>the non-combinatorial information of vertices (point
resp. point-offset pairs, etc.)
<LI>the number of cells
<LI>the indices of the vertices of each cell
<LI>the indices of the neighbors of each cell, where the index
corresponds to the preceding list of cells
<LI>the offsets corresponding to the vertices of the cells
<LI>the non-combinatorial information of each cell
</UL>

\relates Periodic_3_triangulation_3
*/
ostream& operator<< (ostream& os, const Periodic_3_triangulation_3 &t);
} /* end namespace CGAL */
