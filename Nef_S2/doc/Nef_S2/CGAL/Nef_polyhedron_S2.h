namespace CGAL {

/*!
\ingroup PkgNefS2Ref

An instance of data type `Nef_polyhedron_S2<Traits>` is a subset of
the sphere \f$ S_2\f$ that is the result of forming complements and
intersections starting from a finite set `H` of halfspaces bounded by
a plane containing the origin. Halfspaces correspond to hemispheres of
\f$ S_2\f$ and are therefore modeled by oriented great circles of type
`Sphere_circle`. `Nef_polyhedron_S2` is closed under all binary set
operations `intersection`, `union`, `difference`, `complement` and
under the topological operations `boundary`, `closure`, and
`interior`.

\cgalHeading{Parameters}

\code
template< class Nef_polyhedronTraits_S2,
          class Nef_polyhedronItems_S2 = CGAL::SM_items,
          class Nef_polyhedronMarks = bool >
class Nef_polyhedron_S2;
\endcode

The first parameter requires one of the following exact kernels:
`Homogeneous`, `Simple_homogeneous`
parametrized with `Gmpz`, `leda_integer` or any other number type
modeling \f$\mathbb{Z}\f$, or `Cartesian`, `Simple_cartesian` parametrized with
`Gmpq`, `leda_rational`, `Quotient<Gmpz>` or any other number
type modeling \f$\mathbb{Q}\f$.

The second parameter and the third parameter are for future considerations.
Neither `Nef_polyhedronItems_S2` nor `Nef_polyhedronMarks` is
specified, yet. Do not use other than the default types for these two
template parameters.

\cgalHeading{Exploration - Point location - Ray shooting}

As Nef
polyhedra are the result of forming complements and intersections
starting from a set `H` of half-spaces that are defined by
oriented lines in the plane, they can be represented by an attributed
plane map \f$ M = (V,E,F)\f$. For topological queries within `M` the
following types and operations allow exploration access to this
structure.

\cgalHeading{Input and Output}

A Nef polyhedron `N` can be visualized in an open GL window. The
output operator is defined in the file
`CGAL/IO/Nef_polyhedron_2_Window-stream.h`.

\cgalHeading{Implementation}

Nef polyhedra are implemented on top of a halfedge data structure and
use linear space in the number of vertices, edges and facets.
Operations like `empty` take constant time. The operations
`clear`, `complement`, `interior`, `closure`,
`boundary`, `regularization`, input and output take linear
time. All binary set operations and comparison operations take time
\cgalBigO{n \log n} where \f$ n\f$ is the size of the output plus the size of the
input.

The point location and ray shooting operations are implemented in the
naive way. The operations run in linear query time without any
preprocessing.

*/
template< typename Traits >
class Nef_polyhedron_S2 {
public:

/// \name Types
/// @{

/*!
\ingroup PkgNefS2Ref

An object `c` of type `Sphere_circle` is an oriented great
circle on the surface of a unit sphere. Such circles correspond to
the intersection of an oriented plane (that contains the origin) and
the surface of \f$ S_2\f$. The orientation of the great circle is that of a
counterclockwise walk along the circle as seen from the positive
halfspace of the oriented plane.

*/

class Sphere_circle {
public:

/// \name Types
/// @{

/*!
ring type.
*/
typedef unspecified_type RT;

/*!
plane a `Sphere_circle` lies in.
*/
typedef unspecified_type Plane_3;

/// @}

/// \name Creation
/// @{

/*!
creates some great circle.
*/
Sphere_circle();

/*!
If \f$ p\f$ and \f$ q\f$ are
opposite of each other, then we create the unique great circle on \f$ S_2\f$
which contains p and q. This circle is oriented such
that a walk along `c` meets \f$ p\f$ just before the shorter segment
between \f$ p\f$ and \f$ q\f$. If \f$ p\f$ and \f$ q\f$ are opposite of each other then
we create any great circle that contains \f$ p\f$ and \f$ q\f$.
*/
Sphere_circle(const Sphere_point& p,
const Sphere_point& q);

/*!
creates the
circle corresponding to the plane `h`. \pre `h` contains the origin.
*/
Sphere_circle(const Plane_3& h);

/*!
creates the circle orthogonal to the vector \f$ (x,y,z)\f$.
*/
Sphere_circle(const RT& x, const RT& y, const RT& z);

/*!
creates a great circle orthogonal to \f$ c\f$ that contains \f$ p\f$.
\pre \f$ p\f$ is not part of \f$ c\f$.
*/
Sphere_circle(Sphere_circle c, const Sphere_point&
p);

/// @}

/// \name Operations
/// @{

/*!
Returns a sphere circle
in the opposite direction of `c`.
*/
Sphere_circle opposite() ;

/*!
returns true iff
`c` contains `p`.
*/
bool has_on(const Sphere_point& p) ;

/*!
returns the plane supporting `c`.
*/
Plane_3 plane() ;

/*!
returns the point that
is the pole of the hemisphere left of `c`.
*/
Sphere_point orthogonal_pole() ;

/// @}

}; /* end Sphere_circle */
/*!
\ingroup PkgNefS2Ref

An object `p` of type `Sphere_point<R>` is a point on the
surface of a unit sphere. Such points correspond to the nontrivial
directions in space and similarly to the equivalence classes of all
nontrivial vectors under normalization.

\cgalHeading{Operations}

Access to the coordinates is provided by the following operations.
Note that the vector \f$ (x,y,z)\f$ is not normalized.

*/

class Sphere_point {
public:

/// \name Types
/// @{

/*!
ring number type.
*/
typedef unspecified_type RT;

/// @}

/// \name Creation
/// @{

/*!
creates some sphere point.
*/
Sphere_point();

/*!
creates a sphere
point corresponding to the point of intersection of the ray starting
at the origin in direction \f$ (x,y,z)\f$ and the surface of \f$ S_2\f$.
*/
Sphere_point(RT x, RT y, RT z);

/// @}

/// \name Operations
/// @{

/*!
the \f$ x\f$-coordinate.
*/
RT x() ;

/*!
the \f$ y\f$-coordinate.
*/
RT y() ;

/*!
the \f$ z\f$-coordinate.
*/
RT z() ;

/*!
Equality.
*/
bool operator==(const Nef_polyhedron_S2<Traits>::Sphere_point& q) ;

/*!
Inequality.
*/
bool operator!=(const Nef_polyhedron_S2<Traits>::Sphere_point& q) ;

/*!
returns the antipode of `p`.
*/
Sphere_point antipode() ;

/// @}

}; /* end Sphere_point */

/*!
\ingroup PkgNefS2Ref

An object `s` of type `Sphere_segment` is a segment in the
surface of a unit sphere that is part of a great circle through the
origin. Sphere segments are represented by two sphere points \f$ p\f$ and
\f$ q\f$ plus an oriented plane \f$ h\f$ that contains \f$ p\f$ and \f$ q\f$. The plane
determines the sphere segment as follows. Let \f$ c\f$ be the circle in the
intersection of \f$ h\f$ and \f$ S_2\f$. Then \f$ s\f$ is that part of \f$ c\f$ that is
swept, when we rotate \f$ p\f$ into \f$ q\f$ in counterclockwise rotation around
the normal vector of \f$ h\f$ as seen from the positive halfspace.

*/

class Sphere_segment {
public:

/// \name Creation
/// @{

/*!
creates some sphere segment.
*/
Sphere_segment();

/*!
creates a spherical segment spanning
the shorter arc from `p1` to `p2` if `shorter_arc == true`. Otherwise the longer arc is created. \pre `p1 != p2` and `p1 != p2.opposite()`.
*/
Sphere_segment(
const Sphere_point& p1, const Sphere_point& p2, bool shorter_arc=true);

/*!
creates a spherical segment spanning the
arc from `p1` to `p2` as part of the oriented circle `c`
(`p1 == p2` or `p1 == p2.opposite()` are possible.)
\pre `p1` and `p2` are contained in `c`.
*/
Sphere_segment(const Sphere_point& p1,
const Sphere_point& p2, const Sphere_circle& c);

/*!
creates the spherical segment as part of `c1` that is part
of the halfsphere left of the oriented circle `c2`. \pre `c1 != c2` as unoriented circles.
*/
Sphere_segment(const Sphere_circle& c1,
const Sphere_circle& c2);

/// @}

/// \name Operations
/// @{

/*!
the source point of
`s`.
*/
const Sphere_point& source() ;

/*!
the target point of
`s`.
*/
const Sphere_point& target() ;

/*!
the great circle
supporting `s`.
*/
const Sphere_circle& sphere_circle() ;

/*!
returns the spherical
segment oriented from `target()` to `source()` with the same
point set as `s`.
*/
Sphere_segment opposite() ;

/*!
returns the spherical
segment oriented from `target()` to `source()` with the
point set completing `s` to a full circle.
*/
Sphere_segment complement() ;

/*!
a segment is short iff it is shorter
than a half-circle.
*/
bool is_short() ;

/*!
a segment is long iff it is longer than a
half-circle.
*/
bool is_long() ;

/*!
return true iff `s` is degenerate,
i.e.\ source and target are the same.
*/
bool is_degenerate() ;

/*!
return true iff `s` is a perfect half-circle,
i.e.\ `source().antipode == target()`.
*/
bool is_halfcircle() ;

/*!
return true iff
`s` contains `p`.
*/
bool has_on(const Sphere_point& p) ;

/*!

return true iff `s` contains `p` in its relative interior.

*/
bool has_in_relative_interior(const Sphere_point& p) ;

/// @}

}; /* end Sphere_segment */


/*!
\ingroup PkgNefS2Ref

The type `SFace_cycle_iterator` iterates over a list of
`Object_handles`. Each item of that list can either be assigned
to `SVertex_handle`, `SHalfedge_handle` or `SHalfloop_handle`.
To find out which
of these assignment works out, the member functions `is_svertex()`,
`is_shalfedge()` and `is_shalfloop()` are provided.

\sa `CGAL::Nef_polyhedron_S2::SVertex`
\sa `CGAL::Nef_polyhedron_S2::SHalfedge`
\sa `CGAL::Nef_polyhedron_S2::SHalfloop`

*/

class SFace_cycle_iterator {
public:

/// \name Types
/// @{

/*!
const handle to SVertex.
*/
typedef unspecified_type SVertex_handle;

/*!
const handle to SHalfedge.
*/
typedef unspecified_type SHalfedge_handle;

/*!
const handle to SHalfloop.
*/
typedef unspecified_type SHalfloop_handle;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
SFace_cycle_iterator();

/// @}

/// \name Operations
/// @{

/*!
returns true if the iterator represents a `SVertex_handle`.
*/
bool is_svertex() const;

/*!
returns true if the iterator represents a `SHalfedge_handle`.
*/
bool is_shalfedge() const;

/*!
returns true if the iterator represents a `SHalfloop_handle`.
*/
bool is_shalfloop() const;

/*!
casts the iterator to `SVertex_handle`.
*/
operator SVertex_handle() const;

/*!
casts the iterator to `SHalfedge_handle`.
*/
operator SHalfedge_handle() const;

/*!
casts the iterator to `SHalfloop_handle`.
*/
operator SHalfloop_handle() const;

/// @}

}; /* end SFace_cycle_iterator */

/*!
\ingroup PkgNefS2Ref

Figures \ref figureNefS2SVertexIncidences
and \ref figureNefS2SHalfloopIncidences
illustrate the incidences of an sface. An sface is described
by its boundaries. An entry item to each boundary cycle can be accessed
using the iterator range (`sface_cycles_begin()`/`sface_cycles_end()`).
Additionally, `Nef_polyhedron_S2` provides the macro
`CGAL_forall_sface_cylces_of`. The iterators are of type
`SFace_cycle_const_iterator` and represent either a shalfedge, a shalfloop,
or a svertex.

\cgalHeading{Creation}

There is no need for a user to create a `SFace` explicitly. The
class `Nef_polyhedron_S2<Traits>` manages the needed sfaces internally.

\sa `CGAL::Nef_polyhedron_S2`
\sa `CGAL::Nef_polyhedron_S2::SVertex`

*/

class SFace {
public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_S2<Traits>`.
/// @{

/*!
type of mark.
*/
typedef unspecified_type Mark;

/*!
list of Object handles.
*/
typedef unspecified_type Object_list;

/*!
const handle to Vertex.
*/
typedef unspecified_type Vertex_const_handle;

/*!
const handle to Volume.
*/
typedef unspecified_type Volume_const_handle;

/*!
const handle to SFace.
*/
typedef unspecified_type SFace_const_handle;

/*!
const iterator over the entries to all `sface` cycles of a `sface`.
*/
typedef unspecified_type SFace_cycle_const_iterator;

/// @}

/// \name Operations
/// @{

/*!
the mark of the `sface`.
*/
const Mark& mark() const;

/*!
iterator over the entries to all sface cycles of the `sface` .
*/
SFace_cycle_const_iterator sface_cycle_begin() const;

/*!
past-the-end iterator.
*/
SFace_cycle_const_iterator sface_cycle_end() const;

/// @}

}; /* end SFace */

/*!
\ingroup PkgNefS2Ref

A shalfedge is a great arc on a sphere map.
The figure below
depicts the relationship between a shalfedge and its incident
shalfedges, svertices, and sfaces on a sphere map. A shalfedge is
an oriented sedge between two svertices. It is always paired with a
shalfedge pointing in
the opposite direction. The `twin()` member function returns
this shalfedge of opposite orientation.

\anchor figureNefS2SVertexIncidences
\image html shalfedge.png "Incidences of an SHalfedge"
\image latex shalfedge.png "Incidences of an SHalfedge"

The `snext()` member function points
to the successor shalfedge around this sface while the `sprev()` member
function points to the preceding shalfedge. An
successive assignments of the form `se = se->snext()` cycles
counterclockwise around the sface (or hole).

Similarly, the successive
assignments of the form `se = se->snext()->twin()` cycle
clockwise around the svertex and traverse all halfedges incident to
this svertex. The assignment `se = se->cyclic_adj_succ()` can be
used as a shortcut.

A const circulator is provided for each of the two circular orders.
The circulators are bidirectional and assignable to `SHalfedge_const_handle`.

\cgalHeading{Creation}

There is no need for a user to create a `SHalfedge` explicitly. The
class `Nef_polyhedron_S2<Traits>` manages the needed shalfedges internally.

\sa `CGAL::Nef_polyhedron_S2::SVertex`
\sa `CGAL::Nef_polyhedron_S2::SFace`
\sa `CGAL::Nef_polyhedron_S2::Sphere_circle`

*/

class SHalfedge {
public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_S2<Traits>`.
/// @{

/*!
type of mark.
*/
typedef unspecified_type Mark;

/*!
sphere circle type stored in SHalfedge.
*/
typedef unspecified_type Sphere_circle;

/*!
const handle to SVertex.
*/
typedef unspecified_type SVertex_const_handle;

/*!
const handle to SHalfedge.
*/
typedef unspecified_type SHalfedge_const_handle;

/*!
const handle to SFace.
*/
typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
the mark of the `sedge`.
*/
const Mark& mark() const;

/*!
the sphere circle of the `sedge`.
*/
const Sphere_circle& circle() const;

/*!
the twin of the `sedge`.
*/
SHalfedge_const_handle twin() const;

/*!
the source svertex of the `sedge`.
*/
SVertex_const_handle source() const;

/*!
equals `twin()->source()`.
*/
SVertex_const_handle target() const;

/*!
the SHalfedge previous to the `sedge` in a sface cycle.
*/
SHalfedge_const_handle sprev() const;

/*!
the next SHalfedge of the `sedge` in a sface cycle.
*/
SHalfedge_const_handle snext() const;

/*!
the edge before the `sedge` in the cyclic ordered adjacency list of `source()`.
*/
SHalfedge_const_handle cyclic_adj_pred() const;

/*!
the edge after the `sedge` in the cyclic ordered adjacency list of `source()`.
*/
SHalfedge_const_handle cyclic_adj_succ() const;

/*!
the incident `sface` of the `sedge`.
*/
SFace_const_handle incident_sface() const;

/*!
determines whether the `sedge` is
in an outer sface cycle.
*/
bool in_outer_sface_cycle() const;

/*!
determines whether the `sedge` is
in an inner sface cycle.
*/
bool in_inner_sface_cycle() const;

/// @}

}; /* end SHalfedge */

/*!
\ingroup PkgNefS2Ref

A sloop is a great circle on a sphere. A shalfloop is an oriented sloop. It is always paired with a
shalfloop whose supporting `Sphere_circle` is pointing in
the opposite direction. The `twin()` member function returns
this shalfloop of opposite orientation. Each `Nef_polyhedron_S2` can only have one sloop
(resp. two shalfloops).

The figure below
depicts the relationship between a shalfloop and sfaces on a sphere map.

\anchor figureNefS2SHalfloopIncidences
\image html shalfloopB.png "Incidences of an SHalfloop "
\image latex shalfloopB.png "Incidences of an SHalfloop "

\cgalHeading{Creation}

There is no need for a user to create a `SHalfloop` explicitly. The
class `Nef_polyhedron_S2<Traits>` manages the needed shalfloops internally.

\sa `CGAL::Nef_polyhedron_S2::SFace`
\sa `CGAL::Nef_polyhedron_S2::Sphere_circle`

*/

class SHalfloop {
public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_S2<Traits>`.
/// @{

/*!
type of mark.
*/
typedef unspecified_type Mark;

/*!
sphere circle type stored in SHalfloop.
*/
typedef unspecified_type Sphere_circle;

/*!
const handle to SHalfloop.
*/
typedef unspecified_type SHalfloop_const_handle;

/*!
const handle to SFace.
*/
typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
the mark of the halfloop.
*/
const Mark& mark() const;

/*!
the sphere circle of the halfloop.
*/
const Sphere_circle& circle() const;

/*!
the twin of the halfloop.
*/
SHalfloop_const_handle twin() const;

/*!
the incident sface of the halfloop.
*/
SFace_const_handle incident_sface() const;

/// @}

}; /* end SHalfloop */

/*!
\ingroup PkgNefS2Ref

Figure \ref figureNefS2SVertexIncidences illustrates the incidence of a svertex on a sphere map.

The member function
`out_sedge()` returns the first outgoing shalfedge, and `incident_sface()`
returns the incident sface.

\cgalHeading{Creation}

There is no need for a user to create a `SVertex` explicitly. The
class `Nef_polyhedron_S2<Traits>` manages the needed svertices internally.

\sa `CGAL::Nef_polyhedron_S2::SHalfedge`
\sa `CGAL::Nef_polyhedron_S2::SFace`
\sa `CGAL::Nef_polyhedron_S2::Sphere_point`

*/

class SVertex {
public:

/// \name Types
/// The following types are the same as in `Nef_polyhedron_S2<Traits>`.
/// @{

/*!
type of mark.
*/
typedef unspecified_type Mark;

/*!
sphere point type stored in SVertex.
*/
typedef unspecified_type Sphere_point;

/*!
const handle to SVertex.
*/
typedef unspecified_type SVertex_const_handle;

/*!
const handle to SHalfedge.
*/
typedef unspecified_type SHalfedge_const_handle;

/*!
const handle to SFace.
*/
typedef unspecified_type SFace_const_handle;

/// @}

/// \name Operations
/// @{

/*!
the mark of the `svertex`.
*/
const Mark& mark() const;

/*!
the sphere point of the `svertex`.
*/
const Sphere_point& point() const;

/*!
returns |true| if the `svertex` has no adjacent sedges.
*/
bool is_isolated() const;

/*!
the twin of the `svertex`.
*/
SVertex_const_handle twin() const;

/*!
the first out sedge of  the `svertex`.
*/
SHalfedge_const_handle out_sedge() const;

/*!
the incident sface of  the `svertex`.
*/
SFace_const_handle incident_sface() const;

/// @}

}; /* end SVertex */

/*!
non-mutable handle to svertex.
*/
typedef unspecified_type SVertex_const_handle;

/*!
non-mutable handle to shalfedge.
*/
typedef unspecified_type SHalfedge_const_handle;

/*!
non-mutable handle to shalfloop.
*/
typedef unspecified_type SHalfloop_const_handle;

/*!
non-mutable handle to sface.
*/
typedef unspecified_type SFace_const_handle;

/*!
non-mutable iterator over all svertices.
*/
typedef unspecified_type SVertex_const_iterator;

/*!
non-mutable iterator over all shalfedges.
*/
typedef unspecified_type SHalfedge_const_iterator;

/*!
non-mutable iterator over all shalfloops.
*/
typedef unspecified_type SHalfloop_const_iterator;

/*!
non-mutable iterator over all sfaces.
*/
typedef unspecified_type SFace_const_iterator;

/*!
circulating the
adjacency list of an svertex `v`.
*/
typedef unspecified_type SHalfedge_around_svertex_const_circulator;

/*!
circulating the
sface cycle of an sface `f`.
*/
typedef unspecified_type SHalfedge_around_sface_const_circulator;

/*!
iterating all sface cycles of
an sface `f`. The iterator has method `bool is_svertex()`,
`bool is_shalfedge()`, `bool is_shalfloop()`, and can be
converted to the corresponding handles `SVertex_const_handle`,
`SHalfedge_const_handle`, or `SHalfloop_const_handle`.
*/
typedef unspecified_type SFace_cycle_const_iterator;

/*!
attributes of objects (vertices, edges, faces).
*/
typedef unspecified_type Mark;

/*!
size type
*/
typedef unspecified_type size_type;

/*!
construction selection.
*/
enum Boundary { EXCLUDED, INCLUDED };

/*!
construction selection.
*/
enum Content { EMPTY, COMPLETE };

/// @}

/// \name Creation
/// @{

/*!
creates
an instance `N` of type `Nef_polyhedron_S2<K>` and
initializes it to the empty set if `sphere == EMPTY` and to the
whole sphere if `sphere == COMPLETE`.
*/
Nef_polyhedron_S2<K>(Content sphere = EMPTY);

/*!
creates a Nef polyhedron `N` containing the
half-sphere left of `c` including `c` if
`circle==INCLUDED`, excluding `c` if `circle==EXCLUDED`.

*/
Nef_polyhedron_S2<K>(Sphere_circle c, Boundary circle =
INCLUDED);

/*!
creates a Nef polyhedron `N`
from the set of sphere segments in the iterator range
`[first,beyond)`. If the set of sphere segments is a simple
polygon that separates the sphere surface into two regions, then the
polygonal region that is left of the segment `*first` is
selected. The polygonal region includes its boundary if `b = INCLUDED` and excludes the boundary otherwise.
`Forward_iterator` has to be an iterator with value type
`Sphere_segment`.
*/
template <class Forward_iterator>
Nef_polyhedron_S2<K>(Forward_iterator first, Forward_iterator
beyond, Boundary b = INCLUDED);

/// @}

/// \name Operations
/// @{

/*!
makes `N` the
empty set if `plane == EMPTY` and the full plane if `plane == COMPLETE`.
*/
void clear(Content plane = EMPTY) ;

/*!
returns true if `N` is empty, false
otherwise.
*/
bool is_empty() ;

/*!
returns true if `N` is the whole
sphere, false otherwise.
*/
bool is_sphere() ;

/// @}

/// \name Constructive Operations
/// Additionally there are operators `*,+,-,^,!` which implement the
/// binary operations <I>intersection</I>, <I>union</I>,
/// <I>difference</I>, <I>symmetric difference</I>, and the unary
/// operation <I>complement</I> respectively. There are also the
/// corresponding modification operations `<,<=,>,>=,==,!=`.
///
/// There are also comparison operations like `<,<=,>,>=,==,!=` which
/// implement the relations subset, subset or equal, superset,
/// superset or equal, equality, inequality, respectively.

/// @{

/*!
returns the complement
of `N` in the plane.
*/
Nef_polyhedron_S2<K> complement() ;

/*!
returns the interior of
`N`.
*/
Nef_polyhedron_S2<K> interior() ;

/*!
returns the closure of
`N`.
*/
Nef_polyhedron_S2<K> closure() ;

/*!
returns the boundary of
`N`.
*/
Nef_polyhedron_S2<K> boundary() ;

/*!
returns the
regularized polyhedron (closure of interior).
*/
Nef_polyhedron_S2<K> regularization() ;

/*!
returns `N` \f$ \cap\f$ `N1`.
*/
Nef_polyhedron_S2<K> intersection(const
Nef_polyhedron_S2<K>& N1) ;

/*!
returns `N` \f$ \cup\f$ `N1`.
*/
Nef_polyhedron_S2<K> join(const Nef_polyhedron_S2<K>& N1)
;

/*!
returns `N` \f$ -\f$ `N1`.
*/
Nef_polyhedron_S2<K> difference(const Nef_polyhedron_S2<K>&
N1) ;

/*!
returns the symmectric difference
`N - T` \f$ \cup\f$ `T - N`.
*/
Nef_polyhedron_S2<K> symmetric_difference( const
Nef_polyhedron_S2<K>& N1) ;

/// @}

/// \name Statistics and Integrity
/// @{

/*!
returns the number of
svertices.
*/
Size_type number_of_svertices() ;

/*!
returns the number of
shalfedges.
*/
Size_type number_of_shalfedges() ;

/*!
returns the number of sedges.

*/
Size_type number_of_sedges() ;

/*!
returns the number of
shalfloops.
*/
Size_type number_of_shalfloops() ;

/*!
returns the number of sloops.

*/
Size_type number_of_sloops() ;

/*!
returns the number of sfaces.

*/
Size_type number_of_sfaces() ;

/*!
returns the number of
sface cycles.
*/
Size_type number_of_sface_cycles() ;

/*!
calculates the
number of connected components of `P`.
*/
Size_type number_of_connected_components() ;

/*!
print
the statistics of `P`: the number of vertices, edges, and faces.

*/
void print_statistics(std::ostream& os = std::cout) ;

/*!
checks the link structure and the genus of `P`.

*/
void check_integrity_and_topological_planarity(bool
faces=true) ;

/// @}

/// \name Types
/// @{

/*!

a generic handle to an object of the underlying plane map. The kind of
object `(vertex, halfedge, face)` can be determined and the object can
be assigned to a corresponding handle by the three functions:

- `bool assign(Vertex_const_handle& h, Object_handle)`
- `bool assign(Halfedge_const_handle& h, Object_handle)`
- `bool assign(Face_const_handle& h, Object_handle)`

where each function returns `true` iff the assignment to `h` was done.
*/
typedef unspecified_type Object_handle;

/// @}

/// \name Operations
/// @{

/*!
returns true iff the
object `h` is contained in the set represented by `N`.
*/
bool contains(Object_handle h) ;

/*!
returns true
iff the object `h` is contained in the \f$ 1\f$-skeleton of `N`.

*/
bool contained_in_boundary(Object_handle h) ;

/*!
returns a
generic handle `h` to an object (face, halfedge, vertex) of the
underlying plane map that contains the point `p` in its relative
interior. The point `p` is contained in the set represented by
`N` if `N.contains(h)` is true. The location mode flag
`m` allows one to choose between different point location
strategies.
*/
Object_handle locate(const Sphere_point& p) ;

/*!
returns a handle `h` with
`N.contains(h)` that can be converted to a
`Vertex_/Halfedge_/Face_const_handle` as described above. The
object returned is intersected by the ray starting in `p` with
direction `d` and has minimal distance to `p`. The
operation returns an empty `Object_handle` if the ray shoot along
`d` does not hit any object `h` of `N` with
`N.contains(h)`.
*/
Object_handle ray_shoot(const Sphere_point& p, const
Sphere_direction& d) ;

/*!
returns a handle `h` that can be
converted to a `Vertex_/Halfedge_const_handle` as described
above. The object returned is part of the \f$ 1\f$-skeleton of `N`,
intersected by the ray starting in `p` with direction `d`
and has minimal distance to `p`. The operation returns an
empty `Object_handle` if the ray shoot along `d` does not hit any
\f$ 1\f$-skeleton object `h` of `N`. The location mode flag
`m` allows one to choose between different point location
strategies.
*/
Object_handle ray_shoot_to_boundary(const Sphere_point& p,
const Sphere_direction& d) ;

/// @}

/// \name Iteration

/// The list of all objects can be accessed via iterator ranges. For
/// comfortable iteration we also provide iterations macros. The
/// iterator range access operations are of the following kind:
/// - `SVertex_iterator svertices_begin()/svertices_end()`
/// - `SHalfedge_iterator shalfedges_begin()/shalfedges_end()`
/// - `SHalfloop_iterator shalfloops_begin()/shalfloops_end()`
/// - `SFace_iterator sfaces_begin()/sfaces_end()`
///
/// The macros are then
/// - `CGAL_forall_svertices(v,M)`,
/// - `CGAL_forall_shalfedges(e,M)`,
/// - `CGAL_forall_sfaces(f,M)`,
/// - `CGAL_forall_sface_cycles_of(fc,F)`
///
/// where `M` is a sphere map and `F` is a sface.

/// @{

/*!
returns true iff there is
a shalfloop.
*/
bool has_shalfloop() const;

/*!
returns access to the
sloop.
*/
SHalfloop_const_handle shalfloop() const;

/// @}

}; /* end Nef_polyhedron_S2 */

/*!
returns true iff `c1` and `c2` are equal as unoriented
circles.
\relates Nef_polyhedron_S2::Sphere_circle
*/
bool equal_as_sets(const Nef_polyhedron_S2<Traits>::Sphere_circle c1,
                   const Nef_polyhedron_S2<Traits>::Sphere_circle c2);

} /* end namespace CGAL */

