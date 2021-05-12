
namespace CGAL {

/*!
\ingroup PkgNef2Ref

An instance of data type `Nef_polyhedron_2<T>` is a subset of the
plane that is the result of forming complements and intersections
starting from a finite set `H` of
halfspaces. `Nef_polyhedron_2` is closed under all binary set
operations `intersection`, `union`, `difference`,
`complement` and under the topological operations `boundary`,
`closure`, and `interior`.

The template parameter `T` is specified via an extended kernel
concept. `T` must be a model of the concept
`ExtendedKernelTraits_2`.

\cgalHeading{Exploration - Point location - Ray shooting}

As Nef polyhedra are the result of forming complements
and intersections starting from a set `H` of halfspaces that are
defined by oriented lines in the plane, they can be represented by
an attributed plane map \f$ M = (V,E,F)\f$. For topological queries
within `M` the following types and operations allow exploration
access to this structure.

\cgalHeading{Implementation}

Nef polyhedra are implemented on top of a halfedge data structure and
use linear space in the number of vertices, edges and facets.
Operations like `empty` take constant time. The operations
`clear`, `complement`, `interior`, `closure`,
`boundary`, `regularization`, input and output take linear
time. All binary set operations and comparison operations take time
\f$ O(n \log n)\f$ where \f$ n\f$ is the size of the output plus the size of the
input.

The point location and ray shooting operations are implemented in two
flavors. The `NAIVE` operations run in linear query time without
any preprocessing, the `DEFAULT` operations (equals `LMWT`)
run in sub-linear query time, but preprocessing is triggered with the
first operation. Preprocessing takes time \f$ O(N^2)\f$, the sub-linear
point location time is either logarithmic when LEDA's persistent
dictionaries are present or if not then the point location time is
worst-case linear, but experiments show often sublinear runtimes. Ray
shooting equals point location plus a walk in the constrained
triangulation overlayed on the plane map representation. The cost of
the walk is proportional to the number of triangles passed in
direction `d` until an obstacle is met. In a minimum weight
triangulation of the obstacles (the plane map representing the
polyhedron) the theory provides a \f$ O(\sqrt{n})\f$ bound for the number
of steps. Our locally minimum weight triangulation approximates the
minimum weight triangulation only heuristically (the calculation of
the minimum weight triangulation is conjectured to be NP hard). Thus
we have no runtime guarantee but a strong experimental motivation for
its approximation.

\cgalHeading{Example}

Nef polyhedra are parameterized by a so-called extended geometric
kernel. There are three kernels, one based on a homogeneous
representation of extended points called
`Extended_homogeneous<RT>` where `RT` is a ring type providing
additionally a `gcd` operation, one based on a %Cartesian
representation of extended points called `Extended_cartesian<NT>`
where `NT` is a field type, and finally
`Filtered_extended_homogeneous<RT>` (an optimized version of the
first). The following example uses the filtered homogeneous kernel to
construct the intersection of two halfspaces.

\cgalExample{Nef_2/nef_2_intersection.cpp}

After line (*) `N3` is the intersection of `N1` and `N2`.
The member types of `Nef_polyhedron_2< Extended_homogeneous<NT> >`
map to corresponding types of the standard \cgal geometry kernel
(type equality in pseudo-code notation):

\code{.cpp}
CGAL::Nef_polyhedron_2< CGAL::Extended_cartesian< FT > >::Point == CGAL::Cartesian<FT>::Point_2

CGAL::Nef_polyhedron_2< CGAL::Extended_homogeneous< RT > >::Point == CGAL::Homogeneous< RT >::Point_2

CGAL::Nef_polyhedron_2< CGAL::Filtered_extended_homogeneous<RT> >::Point == CGAL::Homogeneous<RT>::Point_2
\endcode

The same holds for the types `Line` and `Direction` in the
local scope of `Nef_polyhedron_2<...>`.

*/
template< typename T >
class Nef_polyhedron_2 {
public:

/// \name Types
/// @{

/*!
the oriented lines modeling halfplanes.
*/
typedef unspecified_type Line;

/*!
the affine points of the plane.
*/
typedef unspecified_type Point;

/*!
directions in our plane.
*/
typedef unspecified_type Direction;

/*!
tag for calling polygon constructor.
*/
typedef unspecified_type Polygons_tag;

/*!
tag for calling polyline constructor.
*/
typedef unspecified_type Polylines_tag;

/*!
construction selection.
*/
enum Boundary { EXCLUDED, INCLUDED };

/*!
construction selection
*/
enum Content { EMPTY, COMPLETE };

/// @}

/// \name Creation
/// @{

/*!

creates an instance `N` of type `Nef_polyhedron_2<T>`
and initializes it to the empty set if `plane == EMPTY`
and to the whole plane if `plane == COMPLETE`.

*/
Nef_polyhedron_2<T>(Content plane = EMPTY);

/*!

creates a Nef polyhedron `N` containing the halfplane left of
`l` including `l` if `line==INCLUDED`, excluding `l`
if `line==EXCLUDED`.
*/
Nef_polyhedron_2<T>(Line l, Boundary line = INCLUDED);

/*!

creates a Nef polyhedron `N` from the simple polygon `P`
spanned by the list of points in the iterator range `[it,end)` and
including its boundary if `b = INCLUDED` excluding the boundary
otherwise. `Forward_iterator` has to be an iterator with
value type `Point`. This construction expects that `P` is
simple. The degenerate cases where `P` contains no point, one
point or spans just one segment (two points) are correctly handled. In
all degenerate cases there's only one unbounded face adjacent to the
degenerate polygon. If `b == INCLUDED` then `N` is just the
boundary. If `b == EXCLUDED` then `N` is the whole plane
without the boundary.
*/
template <class Forward_iterator>
Nef_polyhedron_2<T>(Forward_iterator it, Forward_iterator end,
Boundary b = INCLUDED);

/*!
The iterator range [it, end) defines a set point
ranges, each of which defines the boundary of simple polygon.
*/
template <class Forward_iterator>
Nef_polyhedron_2<T>(Forward_iterator it, Forward_iterator end,
Polygons_tag);

/*!
The iterator range [it, end) defines a set point
ranges, each of which defines a polyline.
*/
template <class Forward_iterator>
Nef_polyhedron_2<T>(Forward_iterator it, Forward_iterator end,
Polylines_tag);

/// @}

/// \name Operations
/// @{

/*!

makes `N` the empty set if `plane == EMPTY` and the
full plane if `plane == COMPLETE`.

*/
void clear(Content plane = EMPTY) ;

/*!

returns true if `N` is empty, false otherwise.
*/
bool is_empty() ;

/*!

returns true if `N` is the whole plane, false otherwise.
*/
bool is_plane() ;

/// @}

/// \name Constructive Operations
/// Additionally there are operators `*,+,-,^,!` which implement the binary
/// operations <I>intersection</I>, <I>join</I>, <I>difference</I>,
/// <I>symmetric difference</I>, and the unary operation
/// <I>complement</I>, respectively. There are also the corresponding
/// modification operations There are also comparison operations like
/// `<,<=,>,>=,==,!=` which implement the relations subset, subset or
/// equal, superset, superset or equal, equality, inequality,
/// respectively.
/// @{

/*!

returns the complement of `N` in the plane.
*/
Nef_polyhedron_2<T> complement() ;

/*!

returns the interior of `N`.
*/
Nef_polyhedron_2<T> interior() ;

/*!

returns the closure of `N`.
*/
Nef_polyhedron_2<T> closure() ;

/*!

returns the boundary of `N`.
*/
Nef_polyhedron_2<T> boundary() ;

/*!

returns the regularized polyhedron (closure of interior).
*/
Nef_polyhedron_2<T> regularization() ;

/*!

returns `N` \f$ \cap\f$ `N1`.
*/
Nef_polyhedron_2<T> intersection(const Nef_polyhedron_2<T>& N1) ;

/*!

returns `N` \f$ \cup\f$ `N1`. Note that "union" is a keyword of C++
and cannot be used for this operation.
*/
Nef_polyhedron_2<T> join(const Nef_polyhedron_2<T>& N1) ;

/*!

returns `N` \f$ -\f$ `N1`.
*/
Nef_polyhedron_2<T> difference(const Nef_polyhedron_2<T>& N1) ;

/*!

returns the symmectric difference `N - T` \f$ \cup\f$ `T - N`.
*/
Nef_polyhedron_2<T> symmetric_difference(const Nef_polyhedron_2<T>& N1) ;

/// @}

/// \name Types
/// @{

/*!
\ingroup PkgNef2Ref

An instance `D` of the data type `Topological_explorer` is a
decorator for interfacing the topological structure of a plane map
`P` (read-only).

A plane map `P` consists of a triple \f$ (V, E, F)\f$ of vertices,
edges, and faces. We collectively call them objects. An edge `e`
is a pair of vertices `(v,w)` with incidence operations `v = source(e)`, `w = target(e)`. The list of all edges with source
`v` is called the adjacency list `A(v)`.

Edges are paired into twins. For each edge `e = (v,w)` there's an
edge `twin(e) = (w,v)` and `twin(twin(e)) == e`\cgalFootnote{The existence of the edge pairs makes `P` a bidirected graph, the `twin` links make `P` a map.}.

An edge `e = (v,w)` knows two adjacent edges `en = next(e)`
and `ep = previous(e)` where `source(en) = w`,
`previous(en) = e` and `target(ep) = v` and `next(ep) = e`. By this symmetric `previous-next` relationship all edges are
partitioned into face cycles. Two edges \f$ e\f$ and \f$ e'\f$ are in the same
face cycle if \f$ e = next^*(e')\f$. All edges `e` in the same
face cycle have the same incident face \f$ f = face(e)\f$. The cyclic
order on the adjacency list of a vertex `v = source(e)` is given
by `cyclic_adj_succ(e) = twin(previous(e))` and
`cyclic_adj_pred(e) = next(twin(e))`.

A vertex `v` is embedded via coordinates `point(v)`. By the
embedding of its source and target an edge corresponds to a
segment. `P` has the property that the embedding is always
<I>order-preserving</I>. This means a ray fixed in `point(v)` of
a vertex `v` and swept around counterclockwise meets the
embeddings of `target(e)` (\f$ e \in A(v)\f$) in the cyclic order
defined by the list order of `A`.

The embedded face cycles partition the plane into maximal connected
subsets of points. Each such set corresponds to a face. A face is
bounded by its incident face cycles. For all the edges in the
non-trivial face cycles it holds that the face is left of the edges.
There can also be trivial face cycles in form of isolated vertices in
the interior of a face. Each such vertex `v` knows its surrounding
face `f = face(v)`.

Plane maps are attributed, for each object \f$ u \in V \cup E \cup F\f$ we
attribute an information `mark(u)` of type `Mark`. `Mark`
fits the concepts assignable, default-constructible, and
equal-comparable.

*/

class Topological_explorer {
public:

/// \name Types

/// @{

/*!
The underlying plane map type
*/
typedef unspecified_type Plane_map;

/*!
The point type of vertices.
*/
typedef unspecified_type Point;

/*!

All objects (vertices, edges, faces) are attributed by a `Mark` object.
*/
typedef unspecified_type Mark;

/*!
The size type.
*/
typedef unspecified_type Size_type;

/// @}

/// \name Circulators
/// Local types are handles, iterators and circulators of the
/// following kind: `Vertex_const_handle`, `Vertex_const_iterator`,
/// `Halfedge_const_handle`, `Halfedge_const_iterator`,
/// `Face_const_handle`, `Face_const_iterator`. Additionally the
/// following circulators are defined.
/// @{

/*!

circulating the outgoing halfedges in \f$ A(v)\f$.
*/
typedef unspecified_type Halfedge_around_vertex_const_circulator;

/*!

circulating the halfedges in the face cycle of a face `f`.

*/
typedef unspecified_type Halfedge_around_face_const_circulator;

/*!

iterating all holes of a face `f`. The type is
convertible to `Halfedge_const_handle`.

*/
typedef unspecified_type Hole_const_iterator;

/*!

iterating all isolated vertices of a face `f`.
The type generalizes `Vertex_const_handle`.

*/
typedef unspecified_type Isolated_vertex_const_iterator;

/// @}

/// \name Operations
/// @{

/*!

returns the source of `e`.
*/
Vertex_const_handle source(Halfedge_const_handle e) ;

/*!

returns the target of `e`.
*/
Vertex_const_handle target(Halfedge_const_handle e) ;

/*!

returns the twin of `e`.
*/
Halfedge_const_handle twin(Halfedge_const_handle e) ;

/*!

returns `true` iff \f$ A(v) = \emptyset\f$.
*/
bool is_isolated(Vertex_const_handle v) ;

/*!

returns one halfedge with source `v`. It's the starting point for
the circular iteration over the halfedges with source `v`.
\pre `!is_isolated(v)`.

*/
Halfedge_const_handle first_out_edge(Vertex_const_handle v) ;

/*!

returns the halfedge with source `v` that is the last
in the circular iteration before encountering `first_out_edge(v)`
again. \pre `!is_isolated(v)`.

*/
Halfedge_const_handle last_out_edge(Vertex_const_handle v) ;

/*!

returns the edge after `e` in the cyclic ordered adjacency list of
`source(e)`.

*/
Halfedge_const_handle cyclic_adj_succ(Halfedge_const_handle e) ;

/*!

returns the edge before `e` in the cyclic ordered adjacency list of
`source(e)`.

*/
Halfedge_const_handle cyclic_adj_pred(Halfedge_const_handle e) ;

/*!

returns the next edge in the face cycle containing `e`.
*/
Halfedge_const_handle next(Halfedge_const_handle e) ;

/*!

returns the previous edge in the face cycle containing `e`.
*/
Halfedge_const_handle previous(Halfedge_const_handle e) ;

/*!

returns the face incident to `e`.
*/
Face_const_handle face(Halfedge_const_handle e) ;

/*!

returns the face incident to `v`. \pre `is_isolated(v)`.
*/
Face_const_handle face(Vertex_const_handle v) ;

/*!

returns a halfedge in the bounding face cycle of `f`
(`Halfedge_const_handle()` if there is no bounding face cycle).

*/
Halfedge_const_handle halfedge(Face_const_handle f) ;

/// @}

/// \name Iteration
/// @{

/*!
iterator over vertices of the map.
*/
Vertex_const_iterator vertices_begin();

/*!
past-the-end iterator for vertices.
*/
Vertex_const_iterator vertices_end();

/*!
iterator over halfedges of the map.
*/
Halfedge_const_iterator halfedges_begin();

/*!
past-the-end iterator for halfedges.
*/
Halfedge_const_iterator halfedges_end();

/*!
iterator over faces of the map.
*/
Face_const_iterator faces_begin();

/*!
past-the-end iterator for faces
*/
Face_const_iterator faces_end();

/*!

returns a circulator for the cyclic adjacency list of `v`.
*/
Halfedge_around_vertex_const_circulator
out_edges(Vertex_const_handle v) ;

/*!

returns a circulator for the outer face cycle of `f`.
*/
Halfedge_around_face_const_circulator
face_cycle(Face_const_handle f) ;

/*!

returns an iterator for all holes in the interior of `f`.
A `Hole_iterator` can be assigned to a
`Halfedge_around_face_const_circulator`.
*/
Hole_const_iterator holes_begin(Face_const_handle f) ;

/*!

returns the past-the-end iterator of `f`.
*/
Hole_const_iterator holes_end(Face_const_handle f) ;

/*!

returns an iterator for all isolated vertices in the interior of `f`.

*/
Isolated_vertex_const_iterator isolated_vertices_begin(
Face_const_handle f) ;

/*!

returns the past the end iterator of `f`.
*/
Isolated_vertex_const_iterator isolated_vertices_end(
Face_const_handle f) ;

/// @}

/// \name Associated Information
/// The type `Mark` is the general attribute of an object.
/// @{

/*!

returns the embedding of `v`.
*/
const Point& point(Vertex_const_handle v) ;

/*!

returns the mark of `v`.
*/
const Mark& mark(Vertex_const_handle v) ;

/*!

returns the mark of `e`.
*/
const Mark& mark(Halfedge_const_handle e) ;

/*!

returns the mark of `f`.
*/
const Mark& mark(Face_const_handle f) ;

/// @}

/// \name Statistics and Integrity
/// @{

/*!

returns the number of vertices.
*/
Size_type number_of_vertices() ;

/*!

returns the number of halfedges.
*/
Size_type number_of_halfedges() ;

/*!

returns the number of halfedge pairs.
*/
Size_type number_of_edges() ;

/*!

returns the number of faces.
*/
Size_type number_of_faces() ;

/*!

returns the number of face cycles.
*/
Size_type number_of_face_cycles() ;

/*!

calculates the number of connected components of `P`.
*/
Size_type number_of_connected_components() ;

/*!

print the statistics of `P`: the number of vertices, edges,
and faces.
*/
void print_statistics(std::ostream& os = std::cout) ;

/*!

checks the link structure and the genus of `P`.
*/
void check_integrity_and_topological_planarity(bool faces=true) ;

/// @}

}; /* end Topological_explorer */

/*!
\ingroup PkgNef2Ref

a decorator to examine the underlying plane map.

An instance `E` of the data type `Explorer` is a decorator to
explore the structure of the plane map underlying the Nef
polyhedron. It inherits all topological adjacency exploration
operations from `Topological_explorer`. `Explorer`
additionally allows one to explore the geometric embedding.

The position of each vertex is given by a so-called extended point,
which is either a standard affine point or the tip of a ray touching
an infinimaximal square frame centered at the origin. A vertex `v`
is called a <I>standard</I> vertex if its embedding is a
<I>standard</I> point and <I>non-standard</I> if its embedding is a
<I>non-standard</I> point. By the straightline embedding of their
source and target vertices, edges correspond to either affine
segments, rays or lines or are part of the bounding frame.

\anchor extsegs
\image html extsegs.png
\image latex extsegs.png
<center><b>
Extended geometry: standard vertices are marked by S, non-standard
vertices are marked by N. <B>A</B>: The possible embeddings of edges:
an affine segment s1, an affine ray s2, an affine line s3. <B>B</B>: A
plane map embedded by extended geometry: note that the frame is
arbitrarily large, the 6 vertices on the frame are at infinity, the
two faces represent a geometrically unbounded area, however they are
topologically closed by the frame edges. No standard point can be
placed outside the frame.
</b></center>

\cgalHeading{Creation}

`Explorer` is copy constructable and assignable. An object can be
obtained via the `Nef_polyhedron_2::explorer()` method of
`Nef_polyhedron_2`.

*/
class Explorer : public Topological_explorer {
public:

/// \name Types
/// Iterators, handles, and circulators are inherited from `Topological_explorer`.
/// @{

/*!
the point type of finite vertices.

*/
typedef unspecified_type Point;

/*!
the ray type of vertices on the frame.

*/
typedef unspecified_type Ray;

/// @}

/// \name Operations
/// @{

/*!
returns true iff
`v`'s position is a standard point.
*/
bool is_standard(Vertex_const_handle v) ;

/*!
returns the standard
point that is the embedding of `v`. \pre `E.is_standard(v)`.
*/
Point point(Vertex_const_handle v) ;

/*!
returns the ray defining
the non-standard point on the frame. \pre `!E.is_standard(v)`.
*/
Ray ray(Vertex_const_handle v) ;

/*!
returns true
iff `e` is part of the infinimaximal frame.
*/
bool is_frame_edge(Halfedge_const_handle e) ;

/// @}

}; /* end Explorer */

/*!
a generic handle to an object of the underlying
plane map. The kind of object `(vertex, halfedge, face)` can
be determined and the object can be assigned to a corresponding
handle by the three functions:

`bool assign(Vertex_const_handle& h, Object_handle)`

`bool assign(Halfedge_const_handle& h, Object_handle)`

`bool assign(Face_const_handle& h, Object_handle)`

where each function returns `true` iff the assignment to
`h` was done.

*/
typedef unspecified_type Object_handle;

/*!
selectionflag for the point location
mode. LMWT stands for Locally Minimum Weight Triangulation, a locally
optimized constrained triangulation where the weight corresponds to
the length of the edges of the triangulation.
*/
enum Location_mode { DEFAULT, NAIVE, LMWT};

/// @}

/// \name Operations
/// @{

/*!

returns true iff the object `h` is contained in the set
represented by `N`.
*/
bool contains(Object_handle h) ;

/*!

returns true iff the object `h` is contained in the \f$ 1\f$-skeleton
of `N`.
*/
bool contained_in_boundary(Object_handle h) ;

/*!
returns a generic handle `h` to an object (face,
halfedge, vertex) of the underlying plane map that contains the point
`p` in its relative interior. The point `p` is contained in
the set represented by `N` if `N.contains(h)` is true. The
location mode flag `m` allows one to choose between different
point location strategies.
*/
Object_handle locate(const Point& p, Location_mode m =
DEFAULT) ;

/*!
returns a handle `h` with
`N.contains(h)`, that can be converted to a
`Vertex_/Halfedge_/Face_const_handle` as described above. The
object returned is intersected by the ray starting in `p` with
direction `d` and has minimal distance to `p`. The operation
returns an empty `Object_handle` if the ray shoot along `d` does
not hit any object `h` of `N` with `N.contains(h)`. The
location mode flag `m` allows one to choose between different
point location strategies.
*/
Object_handle ray_shoot(const Point& p, const Direction& d,
Location_mode m = DEFAULT) ;

/*!
returns a handle `h`,
that can be converted to a `Vertex_/Halfedge_const_handle` as
described above. The object returned is part of the \f$ 1\f$-skeleton of
`N`, intersected by the ray starting in `p` with direction
`d` and has minimal distance to `p`. The operation returns
an empty `Object_handle` if the ray shoot along `d` does not hit
any \f$ 1\f$-skeleton object `h` of `N`. The location mode flag
`m` allows one to choose between different point location
strategies.
*/
Object_handle ray_shoot_to_boundary(const Point& p, const
Direction& d, Location_mode m = DEFAULT) ;

/*!

returns a decorator object that allows read-only access of
the underlying plane map. See the manual page `Explorer` for its
usage.
*/
Explorer explorer() ;

/// @}

}; /* end Nef_polyhedron_2 */
} /* end namespace CGAL */
