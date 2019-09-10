
namespace CGAL {

/*!
\ingroup PkgPeriodic3Triangulation3MainClasses

The class `Periodic_3_regular_triangulation_3` represents a
weighted Delaunay triangulation in three-dimensional periodic space.

\tparam PT must be a model of the concept `Periodic_3RegularTriangulationTraits_3`.

\tparam TDS must be a model of the concept `TriangulationDataStructure_3`.
Its default value is
`Triangulation_data_structure_3<Regular_triangulation_vertex_base_3<PT,Periodic_3_triangulation_ds_vertex_base_3<> >,
                                Regular_triangulation_cell_base_3<PT,Periodic_3_triangulation_ds_cell_base_3<>>>`.

*/
template< typename PT, typename TDS >
class Periodic_3_regular_triangulation_3 :
  public Periodic_3_triangulation_3<PT, TDS>
{
public:

/// \name Types
/// @{

/*!
*/
typedef TDS::Point_3 Bare_point;

/*!
The type for points
`p` of weighted points \f$ {p}^{(w)}=(p,w_p)\f$
*/
typedef TDS::Weighted_point_3 Weighted_point;

/// @}

/// \name Creation
/// @{

/*!
Creates an empty periodic regular triangulation `rt`, with
`domain` as original domain and possibly specifying
a traits class `traits`.
\pre `domain` is a cube.
*/
Periodic_3_regular_triangulation_3(
const Iso_cuboid & domain = Iso_cuboid(0,0,0,1,1,1),
const Geom_traits & traits = Geom_traits());

/*!
Copy constructor.
*/
Periodic_3_regular_triangulation_3 (const
Periodic_3_regular_triangulation_3 & rt1);

/*!
Equivalent to constructing an empty triangulation with the optional
domain and traits class arguments and calling `insert(first,last)`.
\pre The `value_type` of `first` and `last` are `Weighted_point`s lying inside `domain` and `domain` is a cube.
Their weights are non-negative and smaller than 1/64 times the squared cube edge length.
If the fourth argument
`is_large_point_set` is set to `true` a heuristic for
optimizing the insertion of large point sets is applied.
*/
template < class InputIterator >
Periodic_3_regular_triangulation_3 (
InputIterator first,
InputIterator last,
const Iso_cuboid & domain = Iso_cuboid(0,0,0,1,1,1),
const Geom_traits & traits = Geom_traits(),
bool is_large_point_set = false);

/// @}

/// \name Access Functions
/// @{

/*!
Returns the number of hidden points.
*/
size_type number_of_hidden_points() const;

/// @}

/*! \name Insertion

The following methods insert points in the triangulation ensuring the
property that all power spheres are regular. The inserted weighted points need
to lie in the original domain (see Section \ref
P3Triangulation3secspace of the user manual). Note that insertion of a
new point can cause a switch from computing in the 27-sheeted covering
space to computing in the 1-sheeted covering space, which invalidates
some `Vertex_handle`s and `Cell_handle`s.

*/

/// @{

/*!
Inserts point `p` in the triangulation and returns the corresponding
vertex. The optional argument `start` is used as a starting place
for the point location.

If this insertion creates a vertex, this vertex is returned.

If `p` coincides with an existing vertex and has a greater weight,
then the existing weighted point becomes hidden (see
`Periodic_3RegularTriangulationCellBase_3`) and `p` replaces it as vertex
of the triangulation.

If `p` coincides with an already existing vertex (both point and
weights being equal), then this vertex is returned and the triangulation
remains unchanged.

Otherwise if `p` does not appear as a vertex of the triangulation,
then it is stored as a hidden point and this method returns the default
constructed handle.

\pre `p` lies in the original domain `domain`. Its weight is non-negative and smaller than 1/64 times the squared cube edge length.
*/
Vertex_handle insert(const Weighted_point & p,
Cell_handle start = Cell_handle() );

/*!
Inserts point `p` in the triangulation and returns the corresponding
vertex. Similar to the above `insert()` function, but takes as additional
parameter the return values of a previous location query. See description of
`Periodic_3_triangulation_3::locate()`.
\pre `p` lies in the original domain `domain`. Its weight is non-negative and smaller than 1/64 times the squared cube edge length.
*/
Vertex_handle insert(const Weighted_point & p, Locate_type lt,
Cell_handle loc, int li, int lj);


/// @}

/*! \name
The following method allows one to insert several points.
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
\pre The `value_type` of `first` and `last` are `Weighted_point`s
lying inside `domain`.
Their weights are non-negative and smaller than 1/64 times the
squared cube edge length.
*/
template < class InputIterator >
std::ptrdiff_t
insert(InputIterator first, InputIterator last,
bool is_large_point_set = false);

/// @}

/*! \name Removal

*/
/// @{

/*!
Removes the vertex `v` from the triangulation. When computing in
the 27-sheeted covering space it removes all 27 copies of `v`.
*/
void remove(Vertex_handle v);

/*!
Removes the vertices specified by the iterator range (`first, beyond`) of
value type `Vertex_handle`.
`remove()` is called for each element of the range.
The number of vertices removed is returned; this number does not
account for periodic copies of removed vertices.
\pre The iterator must not iterate over several periodic copies of the
same vertex, use e.g. the `Unique_vertex_iterator`.

*/
template < class InputIterator >
std::ptrdiff_t remove(InputIterator first, InputIterator beyond);

/// @}

/// \name Queries

/// @{

/*!
Returns a value indicating the position of the (weighted point, offset) pair (`p`,`off`)
with respect to the power sphere of `c`. More
precisely, it returns:

- `ON_BOUNDED_SIDE` if the angle between the weighted point (`p`,`off`)
and the power sphere of `c` is larger than \f$ \pi/2\f$ or if the ball of
(`p`,`off`) is included in the power sphere of `c`.

- `ON_BOUNDARY` if the ball (`p`,`off`) is orthogonal to the power sphere of `c`.

- `ON_UNBOUNDED_SIDE` if the angle between the weighted point (`p`,`off`)
and the power sphere of `c` is less than \f$ \pi/2\f$ or if the two balls
do not intersect.

*/
Bounded_side
side_of_power_sphere(Cell_handle c, const Weighted_point & p,
const Offset & off = Offset(0,0,0)) const;

/*!
Returns the vertex of the triangulation which is nearest to \f$ p\f$
with respect to the power distance. This means that the power
of the query point `p` with respect to the weighted point in
the returned vertex is smaller than the power of `p`
with respect to the weighted point
in any other vertex. Ties are broken arbitrarily.
It always returns a vertex
corresponding to a point inside `domain` even if computing in a
multiply sheeted covering space.
The default constructed
handle is returned if the triangulation is empty.
The optional argument `c` is a hint
specifying where to start the search.
\pre `c` is a cell of `rt` and `p` lies in the original domain `domain`.

*/
Vertex_handle nearest_power_vertex(const Bare_point & p,
Cell_handle c = Cell_handle()) const;

/// @}

/*! \name
A point-offset pair (`p`,`off`) is said to be in conflict with a cell `c` iff
`rt.side_of_power_sphere(c, p, off)` returns `ON_BOUNDED_SIDE`.
The set of cells that are in conflict with (`p`,`off`) is star-shaped.
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

Returns the `Triple` composed of the resulting output iterators.
\pre `c` is in conflict with `p` and `p` lies in the original domain `domain`.
Its weight is non-negative and smaller than 1/64 times the squared cube edge length.
*/
template <class OutputIteratorBoundaryFacets,
class OutputIteratorCells,
class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets,
OutputIteratorCells,
OutputIteratorInternalFacets>
find_conflicts(const Weighted_point & p, Cell_handle c,
OutputIteratorBoundaryFacets bfit,
OutputIteratorCells cit,
OutputIteratorInternalFacets ifit) const;

/*!
Similar to `find_conflicts()`, but reports the vertices that are on the
boundary of the conflict hole of `p`, in the output iterator `res`.
Returns the resulting output iterator.
\pre `c` is in conflict with `p` and `p` lies in the original domain `domain`.
Its weight is non-negative and smaller than 1/64 times the squared cube edge length.
*/
template <class OutputIterator>
OutputIterator
vertices_in_conflict(const Weighted_point & p, Cell_handle c,
OutputIterator res) const;


/// @}

/* \name
In the weighted setting, a face (cell, facet, edge or vertex) is said to be a
Gabriel face iff the smallest sphere orthogonal to the weighted points
associated to its vertices, has a positive power product with the weighted
point of any other vertex of the triangulation. Any weighted Gabriel
face belongs to the regular triangulation, but the reciprocal is not true.
The following member functions test the Gabriel property of the faces of the regular triangulation.
*/
/// @{

/*

*/
// bool is_Gabriel(Cell_handle c, int i) const;

/*

*/
// bool is_Gabriel(Cell_handle c, int i, int j) const;

/*

*/
// bool is_Gabriel(const Facet& f) const;

/*

*/
// bool is_Gabriel(const Edge& e) const;

/*

*/
// bool is_Gabriel(Vertex_handle v) const;

/// @}


/// \name Power diagram
/// \cgal offers several functions to display the Voronoi diagram of a
/// set of points in 3D. Note that a traits class providing exact
/// constructions should be used in order to guarantee the computation
/// of the Voronoi diagram (as opposed to computing the triangulation
/// only, which requires only exact predicates).
/// @{

/*!
Returns the representative of the weighted circumcenter of the four vertices
of c that lies in the original domain `domain`.
*/
Bare_point dual(Cell_handle c) const;

/*!
Returns the dual of facet `f`, which is a periodic segment.
*/
Periodic_segment_3 dual(Facet f) const;

/*!
same as the previous method for facet `(c,i)`.
\pre \f$ i\in\{0,1,2,3\}\f$
*/
Periodic_segment_3 dual(Cell_handle c, int i) const;

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
Sends the set of duals to all the facets of `rt` into `os`.
*/
template <class Stream> Stream & draw_dual(Stream & os);

/*!
Returns the volume of the cell of the power diagram dual to `v`.
*/
Geom_traits::FT dual_volume(Vertex_handle v) const;

/*!
Returns the centroid of the cell of the power diagram dual to `v`.
*/
Bare_point dual_centroid(Vertex_handle v) const;

/// @}

/// \name Checking
/// These methods are mainly a debugging help for the users of advanced features.
/// @{

/*!
Checks the combinatorial validity of the triangulation and the
validity of its geometric embedding (see
Section \ref P3Triangulation3secintro). Also checks that the power
spheres of all cells are regular.

When `verbose` is set to true, messages describing the first
invalidity encountered are printed.
*/
bool
is_valid(bool verbose = false) const;

/*!
Checks the combinatorial and geometric validity of the cell (see
Section \ref P3Triangulation3secintro). Also checks that its power
sphere is regular.

When `verbose` is set to true, messages are printed to give
a precise indication of the kind of invalidity encountered.
*/
bool
is_valid(Cell_handle c, bool verbose = false) const;

/// @}

}; /* end Periodic_3_regular_triangulation_3 */
} /* end namespace CGAL */
