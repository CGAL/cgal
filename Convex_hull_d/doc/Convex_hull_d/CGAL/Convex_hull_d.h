namespace CGAL {

/*!
\ingroup PkgConvexHullDRef

\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

An instance `C` of type `Convex_hull_d<R>` is the convex hull
of a multi-set `S` of points in \f$ d\f$-dimensional space. We call
`S` the underlying point set and \f$ d\f$ or `dim` the dimension of
the underlying space. We use `dcur` to denote the affine dimension
of `S`. The data type supports incremental construction of hulls.

The closure of the hull is maintained as a simplicial complex, i.e.,
as a collection of simplices. The intersection of any two is a face of
both (The empty set if a facet of every simplex). In the
sequel we reserve the word simplex for the simplices of dimension
`dcur`. For each simplex there is a handle of type
`Simplex_handle` and for each vertex there is a handle of type
`Vertex_handle`. Each simplex has `1 + dcur` vertices
indexed from `0` to `dcur`; for a simplex `s` and an index `i`,
`C.vertex(s,i)` returns the `i`-th vertex of `s`. For any simplex
`s` and any index `i` of `s` there may or may not be a simplex `t`
opposite to the `i`-th vertex of `s`. The function
`C.opposite_simplex(s,i)` returns `t` if it exists and returns
`Simplex_handle()` (the undefined handle) otherwise. If `t` exists
then `s` and `t` share `dcur` vertices, namely all but the vertex
with index `i` of `s` and the vertex with index
`C.index_of_vertex_in_opposite_simplex(s,i)` of `t`. Assume that
`t` exists and let `j = C.index_of_vertex_in_opposite_simplex(s,i)`.
Then `s = C.opposite_simplex(t,j)`
and `i = C.index_of_vertex_in_opposite_simplex(t,j)`.

The boundary of the hull is also a simplicial complex. All simplices
in this complex have dimension `dcur - 1`. For each boundary
simplex there is a handle of type `Facet_handle`. Each facet has
`dcur` vertices indexed from `0` to `dcur - 1`. If `dcur> 1`
then for each facet `f` and each index `i`, `0 <= i < dcur`, there
is a facet `g` opposite to the `i`-th vertex of
`f`. The function `C.opposite_facet(f,i)` returns `g`. Two
neighboring facets `f` and `g` share `dcur - 1` vertices, namely
all but the vertex with index `i` of `f` and the vertex with index
`C.index_of_vertex_in_opposite_facet(f,i)` of `g`. Let `j =
C.index_of_vertex_in_opposite_facet(f,i)`. Then
`f = C.opposite_facet(g,j)` and `i =C.index_of_vertex_in_opposite_facet(g,j)`.

\tparam R must be a model of the concept `ConvexHullTraits_d`.

\cgalHeading{Iteration Statements}

<B>forall_ch_vertices</B>(\f$ v,C\f$) \f$ \{\f$ the vertices of \f$ C\f$ are
successively assigned to \f$ v\f$ \f$ \}\f$

<B>forall_ch_simplices</B>(\f$ s,C\f$) \f$ \{\f$ the simplices of \f$ C\f$ are
successively assigned to \f$ s\f$ \f$ \}\f$

<B>forall_ch_facets</B>(\f$ f,C\f$) \f$ \{\f$ the facets of \f$ C\f$ are
successively assigned to \f$ f\f$ \f$ \}\f$

\cgalHeading{Implementation}

The implementation of type `Convex_hull_d` is based on
\cgalCite{cms-frric-93} and \cgalCite{bms-dgc-94}. The details
of the implementation can be found in the implementation document
available at the download site of this package.

The time and space requirements are input dependent. Let \f$C_1\f$, \f$C_2\f$, \f$C_3\f$,
\f$\ldots\f$ be the sequence of hulls constructed and for a point \f$ x\f$
let \f$ k_i\f$ be the number of facets of \f$ C_i\f$ that are visible from \f$ x\f$
and that are not already facets of \f$ C_{i-1}\f$.

Then the time for inserting \f$ x\f$ is \cgalBigO{dim \sum_i k_i} and
the number of new simplices constructed during the insertion of \f$x\f$
is the number of facets of the hull which were not already facets
of the hull before the insertion.

The data type `Convex_hull_d` is derived from
`Regular_complex_d`. The space requirement of regular complexes is
essentially \f$ 12(dim +2 )\f$ bytes times the number of simplices
plus the space for the points. `Convex_hull_d` needs an additional
\f$ 8 + (4 + x)dim\f$ bytes per simplex where \f$ x\f$ is the space
requirement of the underlying number type and an additional \f$ 12\f$ bytes
per point. The total is therefore \f$ (16 + x)dim + 32\f$ bytes times
the number of simplices plus \f$ 28 + x \cdot dim\f$ bytes times the
number of points.

*/
template< typename R >
class Convex_hull_d {
public:

/*! \name Types
Note that each iterator fits the `Handle` concept, i.e.\ iterators can be
used as handles. Note also that all iterator and handle types come
also in a const flavor, e.g., `Vertex_const_iterator` is the
constant version of `Vertex_iterator`. Const correctness requires
to use the const version whenever the convex hull object is
referenced as constant. The `Hull_vertex_iterator` is convertible
to `Vertex_iterator` and `Vertex_handle`.
*/
/// @{

/*!
the representation class.
*/
typedef unspecified_type R;

/*!
the point type.
*/
typedef unspecified_type Point_d;

/*!
the hyperplane type.
*/
typedef unspecified_type Hyperplane_d;

/*!
handle for simplices.
*/
typedef unspecified_type Simplex_handle;

/*!
handle for facets.
*/
typedef unspecified_type Facet_handle;

/*!
handle for vertices.
*/
typedef unspecified_type Vertex_handle;

/*!
iterator for simplices.
*/
typedef unspecified_type Simplex_iterator;

/*!
iterator for facets.
*/
typedef unspecified_type Facet_iterator;

/*!
iterator for vertices.
*/
typedef unspecified_type Vertex_iterator;

/*!
iterator for vertices that are
part of the convex hull.
*/
typedef unspecified_type Hull_vertex_iterator;

/*!
const iterator for all inserted
points.
*/
typedef unspecified_type Point_const_iterator;

/*!
const iterator for all points
that are part of the convex hull.
*/
typedef unspecified_type Hull_point_const_iterator;

/// @}

/// \name Creation
/// The data type `Convex_hull_d` offers neither copy constructor nor assignment operator.
/// @{

/*!
creates an
instance `C` of type `Convex_hull_d`. The dimension of the
underlying space is \f$ d\f$ and `S` is initialized to the empty point
set. The traits class `R` specifies the models of all types and
the implementations of all geometric primitives used by the convex
hull class. The default model is one of the \f$ d\f$-dimensional
representation classes (e.g., `Homogeneous_d`).
*/
Convex_hull_d<R>(int d, R Kernel = R());

/// @}

/// \name Operations
/// All operations below that take a point `x` as argument have the
/// common precondition that `x` is a point of ambient space.
/// @{

/*!
returns the dimension of ambient space.
*/
int dimension() ;

/*!
returns the affine dimension `dcur` of \f$ S\f$.
*/
int current_dimension() ;

/*!
returns the point associated with vertex \f$ v\f$.
*/
Point_d associated_point(Vertex_handle v) ;

/*!
returns the vertex corresponding to the \f$ i\f$-th vertex of \f$ s\f$.

\pre \f$ 0 \leq i \leq dcur\f$.
*/
Vertex_handle vertex_of_simplex(Simplex_handle s, int i)
;

/*!
same as
`C.associated_point(C.vertex_of_simplex(s,i))`.
*/
Point_d point_of_simplex(Simplex_handle s,int i) ;

/*!
returns the simplex opposite to the \f$ i\f$-th vertex of \f$ s\f$
(`Simplex_handle()` if there is no such simplex). \pre \f$ 0 \leq i \leq dcur\f$.
*/
Simplex_handle opposite_simplex(Simplex_handle s,int i)
;

/*!
returns the index of the vertex opposite to the \f$ i\f$-th vertex of
\f$ s\f$. \pre \f$ 0 \leq i \leq dcur\f$ and there is a simplex opposite to the \f$ i\f$-th vertex of \f$ s\f$.
*/
int index_of_vertex_in_opposite_simplex(Simplex_handle s,int
i) ;

/*!
returns a simplex
of which \f$ v\f$ is a node. Note that this simplex is not unique.
*/
Simplex_handle simplex(Vertex_handle v) ;

/*!
returns the index of \f$ v\f$ in
`simplex(v)`.
*/
int index(Vertex_handle v) ;

/*!
returns the vertex corresponding to the \f$ i\f$-th vertex of \f$ f\f$.
\pre \f$ 0 \leq i <  dcur\f$.
*/
Vertex_handle vertex_of_facet(Facet_handle f, int i)
;

/*!
same as
`C.associated_point(C.vertex_of_facet(f,i))`.
*/
Point_d point_of_facet(Facet_handle f, int i) ;

/*!
returns the facet opposite to the \f$ i\f$-th vertex of \f$ f\f$
(`Facet_handle()` if there is no such facet). \pre \f$ 0 \leq i <  dcur\f$ and `dcur > 1`.
*/
Facet_handle opposite_facet(Facet_handle f, int i)
;

/*!
returns the index of the vertex opposite to the \f$ i\f$-th vertex of
\f$ f\f$. \pre \f$ 0 \leq i <  dcur\f$ and `dcur > 1`.
*/
int index_of_vertex_in_opposite_facet(Facet_handle f, int i);

/*!
returns a hyperplane supporting facet `f`. The hyperplane is
oriented such that the interior of `C` is on the negative side of
it. \pre `f` is a facet of `C` and `dcur > 1`.
*/
Hyperplane_d hyperplane_supporting(Facet_handle f);

/*!
adds point `x`
to the underlying set of points. If \f$ x\f$ is equal to (the point
associated with) a vertex of the current hull this vertex is returned
and its associated point is changed to \f$ x\f$. If \f$ x\f$ lies outside the
current hull, a new vertex `v` with associated point \f$ x\f$ is added
to the hull and returned. In all other cases, i.e., if \f$ x\f$ lies in the
interior of the hull or on the boundary but not on a vertex, the
current hull is not changed and `Vertex_handle()` is returned. If
`CGAL_CHECK_EXPENSIVE` is defined then the validity check
`is_valid(true)` is executed as a post condition.
*/
Vertex_handle insert(const Point_d& x);

/*!
adds `S =
set [first,last)` to the underlying set of points. If any point
`S[i]` is equal to (the point associated with) a vertex of the
current hull its associated point is changed to `S[i]`.
*/
template <typename Forward_iterator> void
insert(Forward_iterator first, Forward_iterator last) ;

/*!
returns true if
\f$ x\f$ is not contained in the affine hull of `S`.
*/
bool is_dimension_jump(const Point_d& x) ;

/*!
returns the list of all facets that are visible from `x`.

\pre `x` is contained in the affine hull of `S`.
*/
std::list<Facet_handle> facets_visible_from(const Point_d&
x);

/*!
returns
`ON_BOUNDED_SIDE` (`ON_BOUNDARY`,`ON_UNBOUNDED_SIDE`) if
`x` is contained in the interior (lies on the boundary, is
contained in the exterior) of `C`. \pre `x` is contained in the affine hull of `S`.
*/
Bounded_side bounded_side(const Point_d& x);

/*!
re-initializes `C` to an empty hull
in \f$ d\f$-dimensional space.
*/
void clear(int d) ;

/*!
returns the number of vertices
of `C`.
*/
int number_of_vertices() ;

/*!
returns the number of facets of
`C`.
*/
int number_of_facets() ;

/*!
returns the number of bounded
simplices of `C`.
*/
int number_of_simplices() ;

/*!
gives information about the size
of the current hull and the number of visibility tests performed.
*/
void print_statistics() ;

/*!
checks the
validity of the data structure. If `throw_exceptions == true`
then the program throws the following exceptions to inform about the
problem.
`chull_has_center_on_wrong_side_of_hull_facet` the
hyperplane supporting a facet has the wrong orientation.

`chull_has_local_non_convexity` a ridge is locally non convex.

`chull_has_double_coverage` the hull has a winding number larger
than 1.
*/
bool is_valid(bool throw_exceptions = false) ;

/// @}

/// \name Lists and Iterators
/// @{

/*!
an iterator of `C` to
start the iteration over all vertices of the complex.
*/
Vertex_iterator vertices_begin() ;

/*!
the past the end iterator
for vertices.
*/
Vertex_iterator vertices_end() ;

/*!
an iterator of `C`
to start the iteration over all simplices of the complex.
*/
Simplex_iterator simplices_begin() ;

/*!
the past the end
iterator for simplices.
*/
Simplex_iterator simplices_end() ;

/*!
an iterator of `C` to
start the iteration over all facets of the complex.
*/
Facet_iterator facets_begin() ;

/*!
the past the end iterator for
facets.
*/
Facet_iterator facets_end() ;

/*!
an iterator to
start the iteration over all vertices of `C` that are part of the
convex hull.
*/
Hull_vertex_iterator hull_vertices_begin() ;

/*!
the past the end
iterator for hull vertices.
*/
Hull_vertex_iterator hull_vertices_end() ;

/*!
returns the start
iterator for all points that have been inserted to construct `C`.
*/
Point_const_iterator points_begin() ;

/*!
returns the past the
end iterator for points.
*/
Point_const_iterator points_end() ;

/*!
returns an
iterator to start the iteration over all points in the convex hull
`C`. Included are points in the interior of facets.
*/
Hull_point_const_iterator hull_points_begin() ;

/*!
returns the
past the end iterator for points in the convex hull.
*/
Hull_point_const_iterator hull_points_end() ;

/*!
each facet of `C` is visited by the visitor object
`V`. `V` has to have a function call operator:
`void
operator()(Facet_handle) const`
*/
template <typename Visitor> void visit_all_facets(const
Visitor& V) ;

/*!
returns a list of
all points that have been inserted to construct `C`.
*/
const std::list<Point_d>& all_points() ;

/*!
returns a list of
all vertices of `C` (also interior ones).
*/
std::list<Vertex_handle> all_vertices() ;

/*!
returns a list
of all simplices in `C`.
*/
std::list<Simplex_handle> all_simplices() ;

/*!
returns a list of
all facets of `C`.
*/
std::list<Facet_handle> all_facets() ;

/// @}

}; /* end Convex_hull_d */
} /* end namespace CGAL */
