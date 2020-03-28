namespace CGAL {

/*!
\ingroup PkgConvexHullDRef

\deprecated This package is deprecated since the version 4.6 of \cgal. The package \ref PkgTriangulations should be used instead.

An instance `DT` of type `Delaunay_d< R, Lifted_R >` is the
nearest and furthest site Delaunay triangulation of a set `S` of
points in some `d`-dimensional space. We call `S` the underlying
point set and `d` or `dim` the dimension of the underlying space.
We use `dcur` to denote the affine dimension of `S`. The data
type supports incremental construction of Delaunay triangulations and
various kind of query operations (in particular, nearest and furthest
neighbor queries and range queries with spheres and simplices).

A Delaunay triangulation is a simplicial complex. All simplices in
the Delaunay triangulation have dimension `dcur`. In the nearest
site Delaunay triangulation the circumsphere of any simplex in the
triangulation contains no point of \f$ S\f$ in its interior. In the
furthest site Delaunay triangulation the circumsphere of any simplex
contains no point of \f$ S\f$ in its exterior. If the points in \f$ S\f$ are
co-circular then any triangulation of \f$ S\f$ is a nearest as well as a
furthest site Delaunay triangulation of \f$ S\f$. If the points in \f$ S\f$ are
not co-circular then no simplex can be a simplex of both
triangulations. Accordingly, we view `DT` as either one or two
collection(s) of simplices. If the points in \f$ S\f$ are co-circular there
is just one collection: the set of simplices of some triangulation.
If the points in \f$ S\f$ are not co-circular there are two
collections. One collection consists of the simplices of a nearest
site Delaunay triangulation and the other collection consists of the
simplices of a furthest site Delaunay triangulation.

For each simplex of maximal dimension there is a handle of type
`Simplex_handle` and for each vertex of the triangulation there is
a handle of type `Vertex_handle`. Each simplex has `1 + dcur`
vertices indexed from `0` to `dcur`. For any simplex `s` and any
index `i`, `DT.vertex_of(s,i)` returns the `i`-th vertex of
`s`. There may or may not be a simplex `t` opposite to the vertex of
`s` with index `i`. The function `DT.opposite_simplex(s,i)`
returns `t` if it exists and returns `Simplex_handle()`
otherwise. If `t` exists then `s` and `t` share `dcur` vertices,
namely all but the vertex with index `i` of `s` and the vertex with
index `DT.index_of_vertex_in_opposite_simplex(s,i)` of `t`.
Assume that `t = DT.opposite_simplex(s,i)` exists and let `
j = DT.index_of_vertex_in_opposite_simplex(s,i)`. Then
`s = DT.opposite_simplex(t,j)` and
`i = DT.index_of_vertex_in_opposite_simplex(t,j)`. In general, a vertex
belongs to many simplices.

Any simplex of `DT` belongs either to the nearest or to the
furthest site Delaunay triangulation or both. The test
`DT.simplex_of_nearest(dt_simplex s)` returns true if `s`
belongs to the nearest site triangulation and the test
`DT.simplex_of_furthest(dt_simplex s)` returns true if `s`
belongs to the furthest site triangulation.


\tparam R must be a model of the concept `DelaunayTraits_d`.
\tparam Lifted_R must be a model of the concept `DelaunayLiftedTraits_d`.

\cgalHeading{Implementation}

The data type is derived from `Convex_hull_d` via
the lifting map. For a point `x` in `d`-dimensional space let
`lift(x)` be its lifting to the unit paraboloid of revolution. There
is an intimate relationship between the Delaunay triangulation of a
point set \f$ S\f$ and the convex hull of `lift(S)`: The nearest site
Delaunay triangulation is the projection of the lower hull and the
furthest site Delaunay triangulation is the upper hull. For
implementation details we refer the reader to the implementation
report available from the \cgal server.

The space requirement is the same as for convex hulls. The time
requirement for an insert is the time to insert the lifted point
into the convex hull of the lifted points.

\cgalHeading{Example}

The abstract data type `Delaunay_d` has a default instantiation by
means of the `d`-dimensional geometric kernel.

\code
#include <CGAL/Homogeneous_d.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Delaunay_d.h>

typedef leda_integer RT;
typedef CGAL::Homogeneous_d<RT> Kernel;
typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
typedef Delaunay_d::Point_d Point;
typedef Delaunay_d::Simplex_handle Simplex_handle;
typedef Delaunay_d::Vertex_handle Vertex_handle;

int main()
{
Delaunay_d T(2);
Vertex_handle v1 = T.insert(Point_d(2,11));
...
}
\endcode

\cgalHeading{Traits Requirements}

`Delaunay_d< R, Lifted_R >` requires the following types from the kernel traits `Lifted_R`:

- `RT`
- `Point_d`
- `Vector_d`
- `Ray_d`
- `Hyperplane_d`

and uses the following function objects from the kernel traits:

- `Construct_hyperplane_d`
- `Construct_vector_d`
- `Vector_to_point_d` / `Point_to_vector_d`
- `Orientation_d`
- `Orthogonal_vector_d`
- `Oriented_side_d` / `Has_on_positive_side_d`
- `Affinely_independent_d`
- `Contained_in_simplex_d`
- `Contained_in_affine_hull_d`
- `Intersect_d`
- `Lift_to_paraboloid_d` / `Project_along_d_axis_d`
- `Component_accessor_d`

`Delaunay_d< R, Lifted_R >` requires the following types from the kernel traits `R`:

- `FT`
- `Point_d`
- `Sphere_d`

- `Construct_sphere_d`
- `Squared_distance_d`
- `Point_of_sphere_d`
- `Affinely_independent_d`
- `Contained_in_simplex_d`
*/
template< typename R, typename Lifted_R >
class Delaunay_d : Convex_hull_d<Lifted_R> {
public:

/*! \name Types

*/
/// @{

/*!
handles to the simplices of the complex.

*/
typedef unspecified_type Simplex_handle;

/*!
handles to vertices of the complex.

*/
typedef unspecified_type Vertex_handle;

/*!
the point type

*/
typedef unspecified_type Point_d;

/*!
the sphere type

*/
typedef unspecified_type Sphere_d;

/*!
interface flags

*/
enum Delaunay_voronoi_kind { NEAREST, FURTHEST };

/*!
the iterator for points.

*/
typedef unspecified_type Point_const_iterator;

/*!
the iterator for vertices.

*/
typedef unspecified_type Vertex_iterator;

/*!
the iterator for simplices.

*/
typedef unspecified_type Simplex_iterator;

/// @}

/// \name Creation
/// The data type `Delaunay_d` offers neither copy constructor nor assignment operator.
/// @{

/*!
creates an instance `DT` of type `Delaunay_d`. The
dimension of the underlying space is `d` and `S` is initialized to the
empty point set. The traits class `R` specifies the models of
all types and the implementations of all geometric primitives used by
the Delaunay class. The traits class `Lifted_R` specifies the models of
all types and the implementations of all geometric primitives used by
the base class of `Delaunay_d< R, Lifted_R >`. The second template parameter defaults to
the first: `Delaunay_d<R> = Delaunay_d<R, Lifted_R = R >`.

*/
Delaunay_d< R, Lifted_R >(int d, R k1 = R(), Lifted_R k2 = Lifted_R());

/// @}

/// \name Operations
/// All operations below that take a point `x` as an argument
/// have the common precondition that `x.dimension() == DT.dimension()`.
/// @{

/*!
returns the dimension of ambient space.

*/
int dimension() ;

/*!
returns the affine dimension of the current point set, i.e.,
`-1` is \f$ S\f$ is empty, `0` if \f$ S\f$ consists of a single point,
`1` if all points of \f$ S\f$ lie on a common line, etc.

*/
int current_dimension() ;

/*!
returns true if `s` is a simplex of the nearest site
triangulation.

*/
bool is_simplex_of_nearest(Simplex_handle s) ;

/*!
returns true if `s` is a simplex of the furthest site
triangulation.

*/
bool is_simplex_of_furthest(Simplex_handle s) ;

/*!
returns the vertex associated with the `i`-th node of `s`.
\pre `0 <= i <= dcur`.

*/
Vertex_handle vertex_of_simplex(Simplex_handle s, int i) ;

/*!
returns the point associated with vertex `v`.

*/
Point_d associated_point(Vertex_handle v) ;

/*!
returns the point associated with the `i`-th vertex of `s`.
\pre `0 <= i <= dcur`.

*/
Point_d point_of_simplex(Simplex_handle s,int i) ;

/*!
returns the simplex opposite to the `i`-th vertex of `s`
(`Simplex_handle()` if there is no such simplex).
\pre `0 <= i <= dcur`.

*/
Simplex_handle opposite_simplex(Simplex_handle s, int i) ;

/*!
returns the index of the vertex opposite to the `i`-th vertex
of `s`.
\pre `0 <= i <= dcur`.

*/
int index_of_vertex_in_opposite_simplex(Simplex_handle s,int i) ;

/*!
returns a simplex of the nearest site triangulation incident
to `v`.

*/
Simplex_handle simplex(Vertex_handle v) ;

/*!
returns the index of `v` in `DT.simplex(v)`.

*/
int index(Vertex_handle v) ;

/*!
returns true if `x` is contained in the closure of simplex `s`.

*/
bool contains(Simplex_handle s, const Point_d& x) ;

/*!
decides whether `DT` is empty.

*/
bool empty() ;

/*!
re-initializes `DT` to the empty Delaunay triangulation.

*/
void clear() ;

/*!
inserts point `x` into `DT` and returns the corresponding
`Vertex_handle`. More precisely, if there is already a vertex `v` in
`DT` positioned at `x` (i.e., `associated_point(v)` is equal to
`x`) then `associated_point(v)` is changed to `x` (i.e.,
`associated_point(v)` is made identical to `x`) and if there is no
such vertex then a new vertex `v` with `associated_point(v) = x` is
added to `DT`. In either case, `v` is returned.

*/
Vertex_handle insert(const Point_d& x) ;

/*!
returns a simplex of the nearest site triangulation
containing `x` in its closure (returns `Simplex_handle()` if `x` lies
outside the convex hull of \f$ S\f$).

*/
Simplex_handle locate(const Point_d& x) ;

/*!
if `DT` contains a vertex `v` with `associated_point(v) = x`
the result is `v` otherwise the result is `Vertex_handle()`.

*/
Vertex_handle lookup(const Point_d& x) ;

/*!
computes a vertex `v` of `DT` that is closest to `x`,
i.e.,
`dist(x,associated_point(v)) = min{dist(x, associated_point(u) | u` \f$\in S\f$ `}`.

*/
Vertex_handle nearest_neighbor(const Point_d& x) ;

/*!
returns the list of all vertices contained in the closure of
sphere \f$ C\f$.

*/
std::list<Vertex_handle> range_search(const Sphere_d& C) ;

/*!
returns the list of all vertices contained in the closure of
the simplex whose corners are given by `A`.
\pre `A` must consist of `d+1` affinely independent points in base space.

*/
std::list<Vertex_handle> range_search(const std::vector<Point_d>& A) ;

/*!
returns a list of all simplices of either the nearest or the
furthest site Delaunay triangulation of `S`.

*/
std::list<Simplex_handle> all_simplices(Delaunay_voronoi_kind k = NEAREST) ;

/*!
returns a list of all vertices of either the nearest or the
furthest site Delaunay triangulation of `S`.

*/
std::list<Vertex_handle> all_vertices(Delaunay_voronoi_kind k = NEAREST) ;

/*!
returns \f$ S\f$.

*/
std::list<Point_d> all_points() ;

/*!
returns the start iterator for points in `DT`.

*/
Point_const_iterator points_begin() ;

/*!
returns the past the end iterator for points in `DT`.

*/
Point_const_iterator points_end() ;

/*!
returns the start iterator for simplices of `DT`.

*/
Simplex_iterator simplices_begin(Delaunay_voronoi_kind k = NEAREST) ;

/*!
returns the past the end iterator for simplices of `DT`.

*/
Simplex_iterator simplices_end() ;

/// @}

}; /* end Delaunay_d */

/*!
constructs a LEDA graph representation of the nearest
(`kind = NEAREST` or the furthest (`kind = FURTHEST`) site
Delaunay triangulation.

\pre `dim() == 2`.

\relates Delaunay_d
*/
template <typename R, typename Lifted_R>
template <typename R, typename Lifted_R> void d2_map(const Delaunay_d<R,Lifted_R>& D, GRAPH< typename Delaunay_d<R,Lifted_R>::Point_d, int >& DTG, typename Delaunay_d<R,Lifted_R>::Delaunay_voronoi_kind k = Delaunay_d<R,Lifted_R>::NEAREST) ;
} /* end namespace CGAL */
