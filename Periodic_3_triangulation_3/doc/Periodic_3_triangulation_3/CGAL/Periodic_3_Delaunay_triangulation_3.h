
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3MainClasses

The class `Periodic_3_Delaunay_triangulation_3` represents a
Delaunay triangulation in three-dimensional periodic space.

\tparam PT must be a model of the concept `Periodic_3DelaunayTriangulationTraits_3`.

\tparam TDS must be a model of the concept `TriangulationDataStructure_3`. Its default value
is `Triangulation_data_structure_3<Triangulation_vertex_base_3<PT,Periodic_3_triangulation_ds_vertex_base_3<>>,Triangulation_cell_base_3<PT,Periodic_3_triangulation_ds_cell_base_3<>>>`.

*/
template< typename PT, typename TDS >
class Periodic_3_Delaunay_triangulation_3 :
  public Periodic_3_triangulation_3<PT, TDS>
{
public:

/// \name Creation
/// @{

/*!
Creates an empty periodic Delaunay triangulation `dt`, with
`domain` as original domain and possibly specifying
a traits class `traits`.
\pre `domain` is a cube.
*/
Periodic_3_Delaunay_triangulation_3(
const Iso_cuboid & domain = Iso_cuboid(0,0,0,1,1,1),
const Geom_traits & traits = Geom_traits());

/*!
Copy constructor.
*/
Periodic_3_Delaunay_triangulation_3 (const
Periodic_3_Delaunay_triangulation_3 & dt1);

/*!
Equivalent to constructing an empty triangulation with the optional
domain and traits class arguments and calling `insert(first,last)`.
\pre The `value_type` of `first` and `last` are `Point`s lying inside `domain` and `domain` is a cube.
*/
template < class InputIterator >
Periodic_3_Delaunay_triangulation_3 (
InputIterator first,
InputIterator last,
const Iso_cuboid & domain = Iso_cuboid(0,0,0,1,1,1),
const Geom_traits & traits = Geom_traits());

/// @}

/*! \name Insertion

The following methods insert points in the triangulation ensuring the
empty sphere property of Delaunay triangulations. The inserted points
need to lie in the original domain (see Section \ref
P3Triangulation3secspace of the user manual). In the degenerate case
when there are co-spherical points, the Delaunay triangulation is
known not to be uniquely defined. In this case, \cgal chooses a
particular Delaunay triangulation using a symbolic perturbation scheme
\cgalCite{cgal:dt-pvr3d-03}. Note that insertion of a new point can cause a
switch from computing in the 27-sheeted covering space to computing in
the 1-sheeted covering space, which invalidates some `Vertex_handle`s
and `Cell_handle`s.

*/

/// @{

/*!
Inserts point `p` in the triangulation and returns the corresponding
vertex. The optional argument `start` is used as a starting place
for the point location.
\pre `p` lies in the original domain `domain`.
*/
Vertex_handle insert(const Point & p,
Cell_handle start = Cell_handle() );

/*!
Inserts point `p` in the triangulation and returns the corresponding
vertex. Similar to the above `insert()` function, but takes as additional
parameter the return values of a previous location query. See description of
`Periodic_3_triangulation_3::locate()`.
\pre `p` lies in the original domain `domain`.
*/
Vertex_handle insert(const Point & p, Locate_type lt,
Cell_handle loc, int li, int lj);

/// @}

/*! \name
The following method allows the insertion of several points and returns
the number of inserted points.
*/
/// @{

/*!
Inserts the points in the iterator range \f$ \left[\right.\f$`first`,
`last`\f$ \left.\right)\f$. Returns the number of inserted points.
This function uses spatial sorting (cf. chapter \ref secspatial_sorting)
and therefore is not guaranteed to insert the points following the
order of `InputIterator`. If the third argument
`is_large_point_set` is set to `true` a heuristic for
optimizing the insertion of large point sets is applied.
\pre The `value_type` of `first` and `last` are `Point`s lying inside `domain`.
*/
template < class InputIterator >
std::ptrdiff_t
insert(InputIterator first, InputIterator last,
bool is_large_point_set = false);

/// @}

/// \name Point moving
/// @{

/*!
Moves the point stored in `v` to `p`, while preserving the Delaunay
property. This performs an action semantically equivalent to `remove(v)`
followed by `insert(p)`, but is supposedly faster when the point has
not moved much. Returns the handle to the new vertex.
\pre `p` lies in the original domain `domain`.
*/
Vertex_handle move_point(Vertex_handle v, const Point & p);

/// @}

/*! \name Removal

The following methods remove points in the triangulation.

When a vertex `v` is removed from a triangulation, all the cells
incident to `v` must be removed, and the polyhedral region consisting
of all the tetrahedra that are incident to `v` must be
re-triangulated. The problem thus reduces to triangulating a polyhedral
region, while preserving its boundary, or to compute a
<I>constrained</I> triangulation. This is known to be sometimes
impossible: the Sch&ouml;nhardt polyhedron cannot be triangulated
\cgalCite{cgal:s-cgehd-98}.
However, when dealing with Delaunay triangulations, the case of such
polyhedra that cannot be re-triangulated cannot happen, so \cgal
proposes a vertex removal.
*/
/// @{

/*!
Removes the vertex `v` from the triangulation. When computing in
the 27-sheeted covering space it removes all 27 copies of `v`.
*/
void remove(Vertex_handle v);

/*!
Removes the vertices specified by the iterator range (`first, beyond`) of value type `Vertex_handle`.
`remove()` is called for each element of the range.
The number of vertices removed is returned; this number does not
account for periodic copies of removed vertices.
\pre The iterator must not iterate over several periodic copies of the same vertex, use e.g. the `Unique_vertex_iterator`.

*/
template < class InputIterator >
std::ptrdiff_t remove(InputIterator first, InputIterator beyond);

/// @}

/// \name Queries

/// @{

/*!
Returns a value indicating on which side of the circumscribed sphere
of `c` the point-offset pair (`p`,`off`) lies. More
precisely, it returns:

- `ON_BOUNDED_SIDE` if (`p`,`off`) is inside the sphere.

- `ON_BOUNDARY` if (`p`,`off`) on the boundary of the sphere.

- `ON_UNBOUNDED_SIDE` if (`p`,`off`) lies outside the
sphere.
\pre `p` lies in the original domain `domain`.
*/
Bounded_side
side_of_sphere(Cell_handle c, const Point & p,
const Offset & off = Offset(0,0,0)) const;

/*!
Returns any nearest vertex to the point `p`, or the default constructed
handle if the triangulation is empty. The optional argument `c` is a hint
specifying where to start the search. It always returns a vertex
corresponding to a point inside `domain` even if computing in a
multiply sheeted covering space.
\pre `c` is a cell of `dt` and `p` lies in the original domain `domain`.

*/
Vertex_handle nearest_vertex(Point p,
Cell_handle c = Cell_handle());

/*!
Returns the vertex of the cell `c` that is nearest to the
point-offset pair (`p`,`off`).
\pre `p` lies in the original domain `domain`.
*/
Vertex_handle nearest_vertex_in_cell(Cell_handle c,
Point p, Offset off = Offset(0,0,0)) const;


/// @}

/*! \name
A point-offset pair (`p`,`off`) is said to be in conflict with a cell `c` iff `dt.side_of_sphere(c, p, off)` returns `ON_BOUNDED_SIDE`. The set of cells that are in conflict with (`p`,`off`) is star-shaped.
*/
/// @{



/*!
Computes the conflict hole induced by `p`. The starting cell
`c` must be in conflict. Then this function returns
respectively in the output iterators:

- `cit`: the cells in conflict.

- `bfit`: the facets on the boundary, that is, the facets
`(t, i)` where the cell `t` is in
conflict, but `t->neighbor(i)` is not.

- `ifit`: the facets inside the hole, that is, delimiting two
cells in conflict.

Returns the pair composed of the resulting output iterators.
\pre `c` is in conflict with `p` and `p` lies in the original domain `domain`.
*/
template <class OutputIteratorBoundaryFacets,
class OutputIteratorCells,
class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets,
OutputIteratorCells,
OutputIteratorInternalFacets>
find_conflicts(Point p, Cell_handle c,
OutputIteratorBoundaryFacets bfit,
OutputIteratorCells cit,
OutputIteratorInternalFacets ifit);

/*!
Similar to `find_conflicts()`, but reports the vertices which are on the
boundary of the conflict hole of `p`, in the output iterator `res`.
Returns the resulting output iterator.
\pre `c` is in conflict with `p` and `p` lies in the original domain `domain`.
*/
template <class OutputIterator>
OutputIterator
vertices_in_conflict(Point p, Cell_handle c,
OutputIterator res);


/// @}

/*! \name

A face (cell, facet or edge) is said to be a Gabriel face iff its
smallest circumscribing sphere do not enclose any vertex of the
triangulation. Any Gabriel face belongs to the Delaunay triangulation,
but the reciprocal is not true. The following member functions test
the Gabriel property of Delaunay faces.

*/
/// @{


/*!

*/
bool is_Gabriel(Cell_handle c, int i);

/*!

*/
bool is_Gabriel(Cell_handle c, int i, int j);

/*!

*/
bool is_Gabriel(const Facet& f);

/*!

*/
bool is_Gabriel(const Edge& e);

/// @}

/// \name Voronoi diagram
/// \cgal offers several functions to display the Voronoi diagram of a
/// set of points in 3D. Note that a traits class providing exact
/// constructions should be used in order to guarantee the computation
/// of the Voronoi diagram (as opposed to computing the triangulation
/// only, which requires only exact predicates).
/// @{

/*!
Returns the representative of the circumcenter of the four vertices
of c that lies in the original domain `domain`.
*/
Point dual(Cell_handle c) const;

/*!
Returns the dual of facet `f`, which is a periodic segment.
*/
Periodic_segment dual(Facet f) const;

/*!
same as the previous method for facet `(c,i)`.
\pre \f$ i\in\{0,1,2,3\}\f$
*/
Periodic_segment dual(Cell_handle c, int i) const;

/*!
Returns in the output iterator the points of the dual polygon of
edge `e` in the same order as the `Facet_circulator` returns
facets incident to the edge `e`. The points form the dual polygon
in \f$ \mathbb R^3\f$, so they do not necessarily all lie inside the
original domain.
*/
template <class OutputIterator>
OutputIterator
dual(Edge e, OutputIterator pts) const;

/*!
same as the previous method for edge `(c,i,j)`.
\pre \f$ i,j\in\{0,1,2,3\}, i\neq j\f$
*/
template <class OutputIterator>
OutputIterator
dual(Cell_handle c, int i, int j, OutputIterator pts) const;

/*!
Returns in the output iterator the points of the dual polyhedron of
vertex `v` in no particular order. The points form the dual
polyhedron in \f$ \mathbb R^3\f$, so they do not necessarily lie all
inside the original domain.
*/
template <class OutputIterator>
OutputIterator
dual(Vertex_handle v, OutputIterator pts) const;

/*!
Sends the set of duals to all the facets of `dt` into `os`.
*/
template <class Stream> Stream & draw_dual(Stream & os);

/*!
Returns the volume of the Voronoi cell dual to `v`.
*/
Geom_traits::FT dual_volume(Vertex_handle v) const;

/*!
Returns the centroid of the Voronoi cell dual to `v`.
*/
Point dual_centroid(Vertex_handle v) const;

/// @}

/// \name Checking
/// These methods are mainly a debugging help for the users of advanced features.
/// @{

/*!
Checks the combinatorial validity of the triangulation and the
validity of its geometric embedding (see
Section \ref P3Triangulation3secintro). Also checks that all the
circumscribing spheres of cells are empty.

When `verbose` is set to true, messages describing the first
invalidity encountered are printed.
*/
bool
is_valid(bool verbose = false) const;

/*!
Checks the combinatorial and geometric validity of the cell (see
Section \ref P3Triangulation3secintro). Also checks that the
circumscribing sphere of cells is empty.

When `verbose` is set to true, messages are printed to give
a precise indication of the kind of invalidity encountered.
*/
bool
is_valid(Cell_handle c, bool verbose = false) const;

/// @}

}; /* end Periodic_3_Delaunay_triangulation_3 */
} /* end namespace CGAL */
