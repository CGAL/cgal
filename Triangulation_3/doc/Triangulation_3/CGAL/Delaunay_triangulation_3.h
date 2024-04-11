
namespace CGAL {

/*!
\ingroup PkgTriangulation3TriangulationClasses

The class `Delaunay_triangulation_3` represents a three-dimensional
Delaunay triangulation.

\tparam Traits is the geometric traits class and must be a model of `DelaunayTriangulationTraits_3`.

\tparam TDS is the triangulation data structure and must be a model of `TriangulationDataStructure_3`.
`Default` may be used, with default type `Triangulation_data_structure_3<Triangulation_vertex_base_3<Traits>,
                                                                         Delaunay_triangulation_cell_base_3<Traits> >`.
Any custom type can be used instead of `Triangulation_vertex_base_3`
and `Delaunay_triangulation_cell_base_3`, provided that they are models of the
concepts `TriangulationVertexBase_3` and `DelaunayTriangulationCellBase_3`,
respectively.

\tparam LP is a tag which must be a `Location_policy<Tag>`:
either `CGAL::Fast_location` or `CGAL::Compact_location`.
`CGAL::Fast_location` offers faster (\cgalBigO{\log n} time) point
location, which can be beneficial when performing point locations or random
point insertions (with no good location hint) in large data sets.
It is currently implemented using an additional triangulation
hierarchy data structure \cgalCite{cgal:d-dh-02}.
The default is `CGAL::Compact_location`, which saves memory (3-5%) by avoiding the need for this
separate data structure, and point location is then performed roughly in
\cgalBigO{n^{1/3}} time.
If the triangulation is parallel (see user manual), the default compact
location policy must be used.
Note that this argument can also come in second position, which can be useful when
the default value for the `TDS` parameter is
satisfactory (this is achieved using so-called deduced parameters).
Note that this argument replaces the functionality
provided before \cgal 3.6 by `Triangulation_hierarchy_3`.
An example of use can be found in the user
manual \ref Triangulation3exfastlocation.

\tparam SLDS is an optional parameter to specify the type of the spatial lock data structure.
        It must be a model of the `SurjectiveLockDataStructure` concept,
        with `Object` being a `Point` (as defined below).
        It is only used if the triangulation data structure used is concurrency-safe (i.e.\ when
        `TDS::Concurrency_tag` is `CGAL::Parallel_tag`).
        The default value is `Spatial_lock_grid_3<Tag_priority_blocking>` if
        the triangulation data structure is concurrency-safe, and `void` otherwise.
        In order to use concurrent operations, the user must provide a
        reference to a `SurjectiveLockDataStructure`
        instance via the constructor or `Triangulation_3::set_lock_data_structure()`.

If `TDS::Concurrency_tag` is `CGAL::Parallel_tag`, some operations,
such as insertion/removal of a range of points, are performed in parallel. See
the documentation of the operations for more details.

\sa `CGAL::Triangulation_3`
\sa `CGAL::Regular_triangulation_3`

*/
template< typename Traits, typename TDS, typename LP, typename SLDS >
class Delaunay_triangulation_3 :
    public Triangulation_3<Traits,
                           Delaunay_triangulation_3<
                             Traits, TDS, LP>::Triangulation_data_structure,
                           SLDS>
{
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
typedef LP Location_policy;

/*!

*/
typedef SLDS Lock_data_structure;

/// @}

/*! \name
In addition to those inherited, the following types are defined, for use by the construction of the Voronoi diagram:
*/
/// @{

/*!

*/
typedef Geom_traits::Line_3 Line;

/*!

*/
typedef Geom_traits::Ray_3 Ray;

/*!

*/
typedef Geom_traits::Plane_3 Plane;

/*!

*/
typedef Geom_traits::Object_3 Object;

/// @}

/// \name Creation
/// @{

/*!
Creates an empty Delaunay triangulation, possibly specifying a traits class
`traits`.
`lock_ds` is an optional pointer to the lock data structure for parallel operations. It
must be provided if concurrency is enabled.
*/
Delaunay_triangulation_3(const Geom_traits& traits = Geom_traits(),
                         Lock_data_structure *lock_ds = nullptr);

/*!
Copy constructor.
The pointer to the lock data structure is not copied. Thus, the copy won't be
concurrency-safe as long as the user has not called `Triangulation_3::set_lock_data_structure()`.
*/
Delaunay_triangulation_3 (const Delaunay_triangulation_3 & dt1);

/*!
Equivalent to constructing an empty triangulation with the optional
traits class argument and calling `insert(first,last)`.
If parallelism is enabled, the points will be inserted in parallel.
*/
template < class InputIterator >
Delaunay_triangulation_3 (InputIterator first, InputIterator last,
                          const Geom_traits& traits = Geom_traits(),
                          Lock_data_structure *lock_ds = nullptr);

/*!
Same as before, with last two parameters in reverse order.
*/
template < class InputIterator >
Delaunay_triangulation_3 (InputIterator first, InputIterator last,
                          Lock_data_structure *lock_ds,
                          const Geom_traits& traits = Geom_traits());

/// @}

/// \name Insertion

/// @{

/*!
Inserts the point `p` in the triangulation and returns the corresponding
vertex. Similar to the insertion in a triangulation, but ensures in
addition the empty sphere property of all the created faces.
The optional argument `start` is used as a starting place for the search.

The optional argument `could_lock_zone` is used by the concurrency-safe
version of the triangulation. If the pointer is not null, the insertion will
try to lock all the cells of the conflict zone, i.e.\ all the vertices that are
inside or on the boundary of the conflict zone. If it succeeds, `*could_lock_zone`
is true, otherwise it is false and the return value is Vertex_handle()
(the point is not inserted). In any case, the locked cells are not unlocked by the
function, leaving this choice to the user.
*/
Vertex_handle insert(const Point & p,
Cell_handle start = Cell_handle(), bool *could_lock_zone = nullptr);

/*!
Same as above but uses `hint` as a starting place for the search.
*/
Vertex_handle insert(const Point & p, Vertex_handle hint,
                     bool *could_lock_zone = nullptr);

/*!
Inserts the point `p` in the triangulation and returns the corresponding
vertex. Similar to the above `insert()` function, but takes as additional
parameter the return values of a previous location query. See description of
`Triangulation_3::locate()`.
*/
Vertex_handle insert(const Point & p, Locate_type lt,
Cell_handle loc, int li, int lj,
bool *could_lock_zone = nullptr);

/*!
Inserts the points in the iterator range `[first,last)`. Returns the number of inserted points.
Note that this function is not guaranteed to insert the points
following the order of `PointInputIterator`, as `spatial_sort()`
is used to improve efficiency.
If parallelism is enabled, the points will be inserted in parallel.

\tparam PointInputIterator must be an input iterator with the value type `Point`.

*/
template < class PointInputIterator >
std::ptrdiff_t
insert(PointInputIterator first, PointInputIterator last);

/*!
Inserts the points in the iterator range  `[first,last)`.
Returns the number of inserted points.
Note that this function is not guaranteed to insert the points
following the order of `PointWithInfoInputIterator`, as `spatial_sort()`
is used to improve efficiency.
If parallelism is enabled, the points will be inserted in parallel.
Given a pair `(p,i)`, the vertex `v` storing `p` also stores `i`, that is
`v.point() == p` and `v.info() == i`. If several pairs have the same point,
only one vertex is created, and one of the objects of type `Vertex::Info` will be stored in the vertex.
\pre `Vertex` must be model of the concept `TriangulationVertexBaseWithInfo_3`.

\tparam PointWithInfoInputIterator must be an input iterator with the value type `std::pair<Point,Vertex::Info>`.

*/
template < class PointWithInfoInputIterator >
std::ptrdiff_t
insert(PointWithInfoInputIterator first, PointWithInfoInputIterator last);

/// @}

/// \name Displacement
/// @{

/*!
If there is not already another vertex placed on `p`,
the triangulation is modified such that the new position of vertex `v`
is `p`, and `v` is returned. Otherwise, the triangulation is not
modified and the vertex at point `p` is returned.
\pre Vertex `v` must be finite.
*/
Vertex_handle move_if_no_collision(Vertex_handle v, const Point & p);

/*!
If there is no collision during the move, this function is the same as
`move_if_no_collision` . Otherwise, `v` is removed and the vertex at point `p`
is returned.
\pre Vertex `v` must be finite.
*/
Vertex_handle move(Vertex_handle v, const Point & p);

/// @}

/*! \name Removal

When a vertex `v` is removed from a triangulation, all the cells
incident to `v` must be removed, and the polyhedral region consisting
of all the tetrahedra that are incident to `v` must be
re-triangulated. So, the problem reduces to triangulating a polyhedral
region, while preserving its boundary, or to compute a
*constrained* triangulation. This is known to be sometimes
impossible: the Sch&ouml;nhardt polyhedron cannot be triangulated
\cgalCite{cgal:s-cgehd-98}. However, when dealing with Delaunay
triangulations, the case of such polyhedra that cannot be
re-triangulated cannot happen, so \cgal proposes a vertex removal.

If,due to some point removals, the size of the Delaunay triangulation
decreases drastically, it might be interesting to defragment the
`CGAL::Compact_container` (used by the `Triangulation_data_structure_3`).
*/
/// @{

/*!
Removes the vertex `v` from the triangulation.

\pre `v` is a finite vertex of the triangulation.
*/
void remove(Vertex_handle v);

/*!
Removes the vertex `v` from the triangulation.

This function is concurrency-safe if the triangulation is concurrency-safe.
It will first
try to lock all the cells adjacent to `v`. If it succeeds, `*could_lock_zone`
is true, otherwise it is false (and the point is not removed). In any case,
the locked cells are not unlocked by the function, leaving this choice to the user.

This function will try to remove `v` only if the removal does not
decrease the dimension.

The return value is only meaningful if `*could_lock_zone` is `true`:
  - returns true if the vertex was removed
  - returns false if the vertex wasn't removed since it would decrease
    the dimension.

\pre `v` is a finite vertex of the triangulation.
\pre `dt`.`dimension()` \f$ =3\f$.

*/
bool remove(Vertex_handle v, bool *could_lock_zone);

/*!
Removes the vertices specified by the iterator range `[first, beyond)`.
The number of vertices removed is returned.
If parallelism is enabled, the points will be removed in parallel.
Note that if at some step, the triangulation dimension becomes lower than 3,
the removal of the remaining points will go on sequentially.

\pre (i) all vertices of the range are finite vertices of the triangulation; and (ii) no vertices are repeated in the range.

\tparam InputIterator must be an input iterator with value type `Vertex_handle`.
*/
template < typename InputIterator >
int remove(InputIterator first, InputIterator beyond);

/*!
This function has exactly the same result and the same preconditions as `remove(first, beyond)`.
The difference is in the implementation and efficiency. This version does not re-triangulate the hole after each
point removal but only after removing all vertices. This is more efficient if (and only if) the removed points
are organized in a small number of connected components of the Delaunay triangulation.
Another difference is that there is no parallel version of this function.

\tparam InputIterator must be an input iterator with value type `Vertex_handle`.
*/
template < typename InputIterator >
int remove_cluster(InputIterator first, InputIterator beyond);

/// @}

/// \name Queries

/// @{

/*!
Returns a value indicating on which side of the circumscribed sphere
of `c` the point `p` lies. More precisely, it returns:

- `ON_BOUNDED_SIDE` if `p` is inside the sphere. For an infinite
cell this means that `p` lies strictly either in the half space
limited by its finite facet and not containing any other point of the
triangulation, or in the interior of the disk circumscribing the
<I>finite</I> facet.

- `ON_BOUNDARY` if p on the boundary of the sphere. For an infinite
cell this means that `p` lies on the circle circumscribing
the <I>finite</I> facet.

- `ON_UNBOUNDED_SIDE` if `p` lies outside the sphere. For an
infinite cell this means that `p` does not satisfy either of the
two previous conditions.
\pre `dt`.`dimension()` \f$ =3\f$.
*/
Bounded_side
side_of_sphere(Cell_handle c, const Point & p) const;

/*!
Returns a value indicating on which side of the circumscribed circle
of `f` the point `p` lies. More precisely, it returns:

- in dimension 3:

- For a finite facet, `ON_BOUNDARY` if `p` lies
on the circle, `ON_UNBOUNDED_SIDE` when it lies in the exterior of
the disk, `ON_BOUNDED_SIDE` when it lies in its interior.

- For an infinite facet, it considers the plane defined by the finite
facet of the same cell, and does the same as in dimension 2 in this
plane.

- in dimension 2:

- For a finite facet, `ON_BOUNDARY` if `p` lies
on the circle, `ON_UNBOUNDED_SIDE` when it lies in the exterior of
the disk, `ON_BOUNDED_SIDE` when it lies in its interior.

- For an infinite facet, `ON_BOUNDARY` if the
point lies on the finite edge of `f` (endpoints included),
`ON_BOUNDED_SIDE` for a point in the open half plane defined
by `f` and not containing any other point of the triangulation,
`ON_UNBOUNDED_SIDE` elsewhere.
\pre `dt`.`dimension()` \f$ \geq2\f$ and in dimension 3, `p` is coplanar with `f`.
*/
Bounded_side
side_of_circle(const Facet & f, const Point & p) const;

/*!
Same as the previous method for facet `i` of cell `c`.
*/
Bounded_side
side_of_circle(Cell_handle c, int i, const Point & p);

/*!
Returns any nearest vertex to the point `p`, or the default constructed
handle if the triangulation is empty. The optional argument `c` is a hint
specifying where to start the search.
\pre `c` is a cell of `dt`.

*/
Vertex_handle nearest_vertex(const Point& p,
                             Cell_handle c = Cell_handle());

/*!
Returns the vertex of the cell `c` that is nearest to \f$ p\f$.
*/
Vertex_handle nearest_vertex_in_cell(const Point& p,
                                     Cell_handle c);

/// @}

/*! \name
A point `p` is said to be in conflict with a cell `c` in dimension 3 (resp.\ a facet `f` in dimension 2) iff `dt.side_of_sphere(c, p)` (resp.\ `dt.side_of_circle(f, p)`) returns `ON_BOUNDED_SIDE`. The set of cells (resp.\ facets in dimension 2) which are in conflict with `p` is connected, and it forms a hole.
*/
/// @{



/*!
Computes the conflict hole induced by `p`. The starting cell
(resp. facet) `c` must be in conflict. Then this function returns
respectively in the output iterators:

- `cit`: the cells (resp. facets) in conflict.

- `bfit`: the facets (resp. edges) on the boundary, that is, the facets
(resp. edges) `(t, i)` where the cell (resp. facet) `t` is in
conflict, but `t->neighbor(i)` is not.

- `could_lock_zone`: The optional argument `could_lock_zone` is used by the concurrency-safe
                     version of the triangulation. If the pointer is not null, the algorithm will
                     try to lock all the cells of the conflict zone, i.e.\ all the vertices that are
                     inside or on the boundary of the conflict zone (as a result, the boundary cells become
                     partially locked). If it succeeds, `*could_lock_zone`
                     is true, otherwise it is false (and the returned conflict zone is only partial). In any case,
                     the locked cells are not unlocked by the function, leaving this choice to the user.

This function can be used in conjunction with `insert_in_hole()` in order
to decide the insertion of a point after seeing which elements of the
triangulation are affected.
Returns the pair composed of the resulting output iterators.
\pre `dt`.`dimension()` \f$ \geq2\f$, and `c` is in conflict with `p`.

*/
template <class OutputIteratorBoundaryFacets,
class OutputIteratorCells>
std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
find_conflicts(const Point& p, Cell_handle c,
               OutputIteratorBoundaryFacets bfit,
               OutputIteratorCells cit, bool *could_lock_zone = nullptr);

/*!
Same as the other `find_conflicts()` function, except that it also
computes the internal facets, i.e.\ the facets common to two cells which
are in conflict with `p`.
Then this function returns respectively in the output iterators:

- `cit`: the cells (resp. facets) in conflict.

- `bfit`: the facets (resp. edges) on the boundary, that is, the facets
(resp. edges) `(t, i)` where the cell (resp. facet) `t` is in
conflict, but `t->neighbor(i)` is not.

- `ifit`: the facets (resp. edges) inside the hole, that is, delimiting
two cells (resp facets) in conflict.

- `could_lock_zone`: The optional argument `could_lock_zone` is used by the concurrency-safe
                     version of the triangulation. If the pointer is not null, the algorithm will
                     try to lock all the cells of the conflict zone, i.e.\ all the vertices that are
                     inside or on the boundary of the conflict zone (as a result, the boundary cells become
                     partially locked). If it succeeds, `*could_lock_zone`
                     is true, otherwise it is false (and the returned conflict zone is only partial). In any case,
                     the locked cells are not unlocked by the function, leaving this choice to the user.

Returns the `Triple` composed of the resulting output iterators.
\pre `dt`.`dimension()` \f$ \geq2\f$, and `c` is in conflict with `p`.

*/
template <class OutputIteratorBoundaryFacets,
          class OutputIteratorCells,
          class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets,
       OutputIteratorCells,
       OutputIteratorInternalFacets>
find_conflicts(const Point& p, Cell_handle c,
               OutputIteratorBoundaryFacets bfit,
               OutputIteratorCells cit,
               OutputIteratorInternalFacets ifit,
               bool *could_lock_zone = nullptr);

/*!
\deprecated This function is renamed `vertices_on_conflict_zone_boundary` since CGAL-3.8.
*/
template <class OutputIterator>
OutputIterator
vertices_in_conflict(const Point& p, Cell_handle c, OutputIterator res);

/*!
Similar to `find_conflicts()`, but reports the vertices which are on the
boundary of the conflict hole of `p`, in the output iterator `res`.
Returns the resulting output iterator.
\pre `dt`.`dimension()` \f$ \geq2\f$, and `c` is in conflict with `p`.

*/
template <class OutputIterator>
OutputIterator
vertices_on_conflict_zone_boundary(const Point& p, Cell_handle c, OutputIterator res);

/// @}

/*! \name
A face (cell, facet or edge) is said to be a Gabriel face iff its smallest circumscribing sphere do not enclose any vertex of the triangulation. Any Gabriel face belongs to the Delaunay triangulation, but the reciprocal is not true. The following member functions test the Gabriel property of Delaunay faces.
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

/*! \name Voronoi Diagram
 \cgal offers several functionalities to display the Voronoi diagram of a set of points in 3D. Note that the user should use a kernel with exact constructions in order to guarantee the computation of the Voronoi diagram (as opposed to computing the triangulation only, which requires only exact predicates).
*/
/// @{

/*!
Returns the circumcenter of the four vertices of c.
\pre `dt`.`dimension()`\f$ =3\f$ and `c` is not infinite.
*/
Point dual(Cell_handle c) const;

/*!
Returns the dual of facet `f`, which is

in dimension 3: either a segment, if the two cells incident to `f`
are finite, or a ray, if one of them is infinite;

in dimension 2: a point.
\pre `dt`.`dimension()` \f$ \geq2\f$ and `f` is not infinite.
*/
Object dual(Facet f) const;

/*!
same as the previous method for facet `(c,i)`.
*/
Object dual(Cell_handle c, int i) const;

/*!
returns the line supporting the dual of the facet.
\pre `dt`.`dimension()` \f$ \geq2\f$ and `f` is not infinite.
*/
Line dual_support(Cell_handle c, int i) const;

/*!
Sends the set of duals to all the facets of `dt` into `os`.
*/
template <class Stream> Stream & draw_dual(Stream & os);

/// @}

/// \name Checking
/// These methods are mainly a debugging help for the users of advanced features.
/// @{

/*!
Checks the combinatorial validity of the triangulation and the
validity of its geometric embedding (see
Section \ref Triangulation3secintro). Also checks that all the
circumscribing spheres (resp. circles in dimension 2) of cells
(resp. facets in dimension 2) are empty.
When `verbose` is set to
true, messages describing the first invalidity encountered are
printed.
*/
bool
is_valid(bool verbose = false) const;

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the combinatorial and geometric validity of the cell (see
Section \ref Triangulation3secintro). Also checks that the
circumscribing sphere (resp. circle in dimension 2) of cells
(resp. facet in dimension 2) is empty.

When `verbose` is set to
true, messages are printed to give
a precise indication of the kind of invalidity encountered.
\cgalDebugEnd
*/
bool
is_valid(Cell_handle c, bool verbose = false) const;

/// @}

}; /* end Delaunay_triangulation_3 */
} /* end namespace CGAL */
