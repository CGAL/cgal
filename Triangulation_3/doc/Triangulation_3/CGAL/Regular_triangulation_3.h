
namespace CGAL {
/*!
\ingroup PkgTriangulation3TriangulationClasses

Let \f$ {S}^{(w)}\f$ be a set of weighted points in \f$ \mathbb{R}^3\f$. Let 
\f$ {p}^{(w)}=(p,w_p), p\in\mathbb{R}^3, w_p\in\mathbb{R}\f$ and 
\f$ {z}^{(w)}=(z,w_z), z\in\mathbb{R}^3, w_z\in\mathbb{R}\f$ be two weighted points. 
A weighted point 
\f$ {p}^{(w)}=(p,w_p)\f$ can also be seen as a sphere of center \f$ p\f$ and 
radius \f$ \sqrt{w_p}\f$. 
The <I>power product</I> (or <I>power distance</I> ) 
between \f$ {p}^{(w)}\f$ and \f$ {z}^{(w)}\f$ is 
defined as 
\f[ \Pi({p}^{(w)},{z}^{(w)}) = {\|{p-z}\|^2-w_p-w_z} \f] 
where \f$ \|{p-z}\|\f$ is the Euclidean distance between \f$ p\f$ and \f$ z\f$. 
\f$ {p}^{(w)}\f$ and \f$ {z}^{(w)}\f$ 
are said to be <I>orthogonal</I> if \f$ \Pi{({p}^{(w)}-{z}^{(w)})} 
= 0\f$ (see \cgalFigureRef{Triangulation3figortho}).

Four weighted points have a unique common orthogonal weighted point called 
the <I>power sphere</I>. A sphere \f$ {z}^{(w)}\f$ is said to be 
<I>regular</I> if \f$ \forall {p}^{(w)}\in{S}^{(w)}, 
\Pi{({p}^{(w)}-{z}^{(w)})}\geq 0\f$. 

A triangulation of \f$ {S}^{(w)}\f$ is <I>regular</I> if the power spheres 
of all simplices are regular. 

\tparam Traits is the geometric traits class, and must be a model of `RegularTriangulationTraits_3`

\tparam TDS is the triangulation data structure and must be a model of `TriangulationDataStructure_3`.
TDS has default value `Triangulation_data_structure_3<Regular_triangulation_vertex_base_3<Traits>,
                                                      Regular_triangulation_cell_base_3<Traits> >`.
Any custom type can be used instead of `Regular_triangulation_vertex_base_3`
and `Regular_triangulation_cell_base_3`, provided that they are models of the
concepts `RegularTriangulationVertexBase_3` and `RegularTriangulationCellBase_3`,
respectively.

\tparam SLDS is an optional parameter to specify the type of the spatial lock data structure.
        It must be a model of the `SurjectiveLockDataStructure` concept,
        with `Object` being a `Point`.
        It is only used if the triangulation data structure used is concurrency-safe (i.e.\ when 
        `TDS::Concurrency_tag` is `Parallel_tag`).
        The default value is `Spatial_lock_grid_3<Tag_priority_blocking>` if
        the triangulation data structure is concurrency-safe, and `void` otherwise.
        In order to use concurrent operations, the user must provide a
        reference to a `SLDS`
        instance via the constructor or `Triangulation_3::set_lock_data_structure`.
        
If `TDS::Concurrency_tag` is `Parallel_tag`, some operations,
such as insertion/removal of a range of points, are performed in parallel. See 
the documentation of the operations for more details.

\sa `CGAL::Triangulation_3`
\sa `CGAL::Delaunay_triangulation_3`

*/
template< typename Traits, typename TDS, typename SLDS >
class Regular_triangulation_3
  : public Triangulation_3<Traits, TDS, SLDS> {
public:

/// \name Types 
/// @{

/*!
The type for points 
`p` of weighted points \f$ {p}^{(w)}=(p,w_p)\f$ 
*/ 
typedef Traits::Point_3 Bare_point;

/*!

*/ 
typedef Traits::Weighted_point_3 Weighted_point;

/*!

*/ 
typedef SLDS Lock_data_structure;

/// @} 

/// \name Creation 
/// @{

/*!
Creates an empty regular triangulation, possibly specifying a traits class 
`traits`. 
`lock_ds` is an optional pointer to the lock data structure for parallel operations. It
must be provided if concurrency is enabled.
*/
Regular_triangulation_3(const Traits & traits = Traits(),
                        Lock_data_structure *lock_ds = NULL);

/*!
Copy constructor. 
The pointer to the lock data structure is not copied. Thus, the copy won't be
concurrency-safe as long as the user has not called `Triangulation_3::set_lock_data_structure`.
*/ 
Regular_triangulation_3(const Regular_triangulation_3 & rt1);

/*!
Equivalent to constructing an empty triangulation with the optional 
traits class argument and calling `insert(first,last)`. 
If parallelism is enabled, the points will be inserted in parallel.
\tparam InputIterator must be an input iterator with value type `Weighted_point`. 
*/ 
template < class InputIterator > 
Regular_triangulation_3 (InputIterator first, InputIterator last,
                         const Traits& traits = Traits(),
                         Lock_data_structure *lock_ds = NULL);

/*! 
Same as before, with last two parameters in reverse order.
*/ 
template < class InputIterator > 
Regular_triangulation_3 (InputIterator first, InputIterator last, 
                         Lock_data_structure *lock_ds,
                         const Traits& traits = Traits());
/// @} 

/*!\name Insertion 
The following methods, which already exist in `Triangulation_3`, are overloaded to ensure the property that all power spheres are regular. The following method allows one to insert several points. 
*/
/// @{

/*!
Inserts weighted point `p` in the triangulation. The optional 
argument `start` is used as a starting place for the search. 

If this insertion creates a vertex, this vertex is returned. 

If `p` coincides with an existing vertex and has a greater weight, 
then the existing weighted point becomes hidden (see 
`RegularTriangulationCellBase_3`) and `p` replaces it as vertex 
of the triangulation. 

If `p` coincides with an already existing vertex (both point and 
weights being equal), then this vertex is returned and the triangulation 
remains unchanged. 

Otherwise if `p` does not appear as a vertex of the triangulation, 
then it is stored as a hidden point and this method returns the default 
constructed handle. 

The optional argument `could_lock_zone` is used by the concurrency-safe
version of the triangulation. If the pointer is not null, the insertion will
try to lock all the cells of the conflict zone, i.e.\ all the vertices that are
inside or on the boundary of the conflict zone. If it succeeds, `*could_lock_zone`
is true, otherwise it is false (and the point is not inserted). In any case, 
the locked cells are not unlocked by the function, leaving this choice to the user.
*/ 
Vertex_handle insert(const Weighted_point & p, 
Cell_handle start = Cell_handle(), bool *could_lock_zone = NULL); 

/*!
Same as above but uses `hint` as a starting place for the search. 
*/ 
Vertex_handle insert(const Weighted_point & p, Vertex_handle hint, bool *could_lock_zone = NULL); 

/*!
Inserts weighted point `p` in the triangulation and returns the corresponding 
vertex. Similar to the above `insert()` function, but takes as additional 
parameter the return values of a previous location query. See description of 
`Triangulation_3::locate()`. 
*/ 
Vertex_handle insert(const Weighted_point & p, Locate_type lt, 
Cell_handle loc, int li, int lj, bool *could_lock_zone = NULL); 

/*!
Inserts the weighted points in the range `[first,last)`. 
It returns the difference of the number of vertices between after and 
before the insertions (it may be negative due to hidden points). 
Note that this function is not guaranteed to insert the points 
following the order of `InputIterator`, as `spatial_sort()` 
is used to improve efficiency. 
If parallelism is enabled, the points will be inserted in parallel.

\tparam InputIterator must be an input iterator with value type `Weighted_point`. 
*/ 
template < class InputIterator > 
std::ptrdiff_t 
insert(InputIterator first, InputIterator last); 

/*!

Inserts the weighted points in the iterator range  `[first,last)`.
It returns the difference of the number of vertices between after and 
before the insertions (it may be negative due to hidden points). 
Note that this function is not guaranteed to insert the weighted points 
following the order of `WeightedPointWithInfoInputIterator`, as `spatial_sort()` 
is used to improve efficiency. 
If parallelism is enabled, the points will be inserted in parallel.
Given a pair `(p,i)`, the vertex `v` storing `p` also stores `i`, that is 
`v.point() == p` and `v.info() == i`. If several pairs have the same point, 
only one vertex is created, one of the objects of type `Vertex::Info` will be stored in the vertex. 
\pre `Vertex` must be model of the concept `TriangulationVertexBaseWithInfo_3`.
\tparam (WeightedPointWithInfoInputIterator must be an input iterator with value type  `std::pair<Weighted_point,Vertex::Info>`. 
*/ 
template < class WeightedPointWithInfoInputIterator > 
std::ptrdiff_t 
insert(WeightedPointWithInfoInputIterator first, WeightedPointWithInfoInputIterator last); 



/// @}

/*! \name
The following methods, which already exist in `Triangulation_3`, are overloaded to ensure that hidden points are well created and maintained.
*/
/// @{

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

If the hole contains interior vertices, each of them is hidden by the insertion 
of `p` and is stored in the new cell which contains it. 
\pre `rt`.`dimension()` \f$ \geq2\f$, the set of cells (resp. facets in dimension 2) is connected, not empty, its boundary is connected, and `p` lies inside the hole, which is star-shaped wrt `p`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(Weighted_point p, CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i); 

/*!
Same as above, except that `newv` will be used as the new vertex, which 
must have been allocated previously with, e.g.\ `create_vertex`. 
*/ 
template <class CellIt> 
Vertex_handle insert_in_hole(Weighted_point p, CellIt cell_begin, CellIt cell_end, 
Cell_handle begin, int i, Vertex_handle newv); 

/// @} 

/// \name Removal 
/// @{

/*!
Removes the vertex `v` from the triangulation. 
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
The return value is only meaningful if `*could_lock_zone` is true:
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

/// @} 

/*! \name Queries 
Let us remark that \f$ \Pi({p}^{(w)}-{z}^{(w)}) > 0 \f$ is equivalent to `p` lies outside the sphere with center `z` and radius \f$ \sqrt{w_p^2+w_z^2}\f$. This remark helps provide an intuition about the following predicates.

\anchor Triangulation3figsidedim2
\image html sidedim2.svg side_of_power_circle
\image latex sidedim2.png side_of_power_circle
*/

/// @{

/*!
Returns the position of the weighted point \f$ p\f$ with respect to the 
power sphere of `c`. More precisely, it returns: 

- `ON_BOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})<0\f$ where 
\f$ {z(c)}^{(w)}\f$ is the power sphere of `c`. For an 
infinite cell this means either that `p` lies strictly in the half 
space limited by its finite facet and not containing any other point 
of the triangulation, or that the angle 
between `p` and the power circle of the <I>finite</I> facet of `c` 
is greater than \f$ \pi/2\f$. 

- `ON_BOUNDARY` if p is orthogonal to the power sphere of `c` 
i.e.\ \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})=0\f$. For an infinite cell this means 
that `p` is orthogonal to the power circle of its <I>finite</I> facet. 

- `ON_UNBOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})>0\f$ 
i.e.\ the angle between the weighted point `p` and the power sphere 
of `c` is less than \f$ \pi/2\f$ or if these two spheres do not 
intersect. For an 
infinite cell this means that `p` does not satisfy either of the 
two previous conditions. 
\pre `rt`.`dimension()` \f$ =3\f$. 
*/ 
Bounded_side 
side_of_power_sphere(Cell_handle c, const Weighted_point & p) const; 

/*!
Returns the position of the point `p` with respect to the 
power circle of `f`. More precisely, it returns: 

- in dimension 3: 

- For a finite facet, 

`ON_BOUNDARY` if `p` is orthogonal to the power circle in the 
plane of the facet, 

`ON_UNBOUNDED_SIDE` when their angle is less than \f$ \pi/2\f$, 

`ON_BOUNDED_SIDE` when it is greater than \f$ \pi/2\f$ (see 
Figure \ref Triangulation3figsidedim2). 

- For an infinite facet, it considers the plane defined by the finite 
facet of the cell `f.first`, and does the same as in 
dimension 2 in this plane. 

- in dimension 2: 

- For a finite facet, 

`ON_BOUNDARY` if `p` is orthogonal to the circle, 

`ON_UNBOUNDED_SIDE` when the angle between `p` and the 
power circle of `f` is less than \f$ \pi/2\f$, 
`ON_BOUNDED_SIDE` when it is greater than \f$ \pi/2\f$. 

- For an infinite facet, 

`ON_BOUNDED_SIDE` for a point in the open half plane defined by 
`f` and not containing any other point of the triangulation, 

`ON_UNBOUNDED_SIDE` in the other open half plane. 

If the point `p` is collinear with the finite edge `e` of 
`f`, it returns: 

`ON_BOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(e)}^{(w)})<0\f$, where 
\f$ {z(e)}^{(w)}\f$ is the power segment of `e` in the line supporting 
`e`, 

`ON_BOUNDARY` if \f$ \Pi({p}^{(w)}-{z(e)}^{(w)})=0\f$, 

`ON_UNBOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(e)}^{(w)})>0\f$ . 
\pre `rt`.`dimension()` \f$ \geq2\f$. 
*/ 
Bounded_side 
side_of_power_circle(const Facet & f, 
const Weighted_point & p) const; 

/*!
Same as the previous method for facet `i` of cell `c`. 
*/ 
Bounded_side 
side_of_power_circle(Cell_handle c, int i, 
const Weighted_point & p) const; 

/*!
In dimension 1, returns 

`ON_BOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})<0\f$, where 
\f$ {z(c)}^{(w)}\f$ is the power segment of the edge represented by 
`c`, 

`ON_BOUNDARY` if \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})=0\f$, 

`ON_UNBOUNDED_SIDE` if \f$ \Pi({p}^{(w)}-{z(c)}^{(w)})>0\f$ . 
\pre `rt`.`dimension()` \f$ = 1\f$. 
*/ 
Bounded_side 
side_of_power_segment(Cell_handle c, const Weighted_point & p) 
const; 

/*!
Returns the vertex of the triangulation which is nearest to \f$ p\f$ 
with respect to the power distance. This means that the power 
of the query point `p` with respect to the weighted point in 
the returned vertex is smaller than the power of `p` 
with respect to the weighted point 
in any other vertex. Ties are broken arbitrarily. 
The default constructed 
handle is returned if the triangulation is empty. 
The optional argument `c` is a hint 
specifying where to start the search. 
\pre `c` is a cell of `rt`. 

*/ 
Vertex_handle nearest_power_vertex(Weighted_point p, 
Cell_handle c = Cell_handle()); 

/*!
Returns the vertex of the cell `c` 
that is nearest to \f$ p\f$ 
with respect to the power distance. 

*/ 
Vertex_handle nearest_power_vertex_in_cell(Weighted_point p, 
Cell_handle c); 


/// @}

/*! \name

A weighted point `p` is said to be in conflict with a cell `c` in dimension 3 (resp.\ with a facet `f` in dimension 2) if it has a negative power distance to the power sphere of `c` (resp.\ to the power circle of `f`).\ The set of cells (resp.\ facets in dimension 2) which are in conflict with `p` is connected. 
*/
/// @{

/*!

Compute the conflicts with `p`. 

@param p                  The query point.
@param c                  The starting cell.
@param cit                The cells (resp. facets) in conflict with `p`. 
@param bfit               The facets (resp. edges) on the boundary of the conflict zone, that is, the facets  (resp.\ edges) `(t, i)` where the cell (resp.. facet) `t` is in conflict, but `t->neighbor(i)` is not. 
@param ifit               The facets (resp.\ edges) inside the conflict zone, that facets incident to two cells (resp.\ facets) in conflict. 
@param could_lock_zone    The optional argument `could_lock_zone` is used by the concurrency-safe
                          version of the triangulation. If the pointer is not null, the algorithm will
                          try to lock all the cells of the conflict zone, i.e.\ all the vertices that are
                          inside or on the boundary of the conflict zone (as a result, the boundary cells become
                          partially locked). If it succeeds, `*could_lock_zone`
                          is true, otherwise it is false (and the returned conflict zone is only partial). In any case, 
                          the locked cells are not unlocked by the function, leaving this choice to the user.
@param this_facet_must_be_in_the_cz 
                          If the optional argument `this_facet_must_be_in_the_cz` is not null, the algorithm will check
                          if this facet is in the conflict zone (it may be internal as well as boundary).
@param the_facet_is_in_its_cz 
                          This argument must be not null if the previous `this_facet_must_be_in_the_cz` argument is not null. 
                          The boolean value pointed by this pointer is set to true if *`this_facet_must_be_in_the_cz` is
                          among the internal or boundary facets of the conflict zone, and false otherwise.                         

\pre  The starting cell (resp.\ facet) `c` must be in conflict with `p`. 
\pre `rt`.`dimension()` \f$ \geq2\f$, and `c` is in conflict with `p`. 

\return the `Triple` composed of the resulting output iterators. 


*/ 
template <class OutputIteratorBoundaryFacets, 
class OutputIteratorCells, 
class OutputIteratorInternalFacets> 
Triple<OutputIteratorBoundaryFacets, 
OutputIteratorCells, 
OutputIteratorInternalFacets> 
find_conflicts(const Weighted_point p, Cell_handle c, 
OutputIteratorBoundaryFacets bfit, 
OutputIteratorCells cit, 
OutputIteratorInternalFacets ifit,
bool *could_lock_zone = NULL,
const Facet *this_facet_must_be_in_the_cz = NULL,
bool *the_facet_is_in_its_cz = NULL);

/*!
\deprecated This function is renamed `vertices_on_conflict_zone_boundary` since CGAL-3.8. 
*/ 
template <class OutputIterator> 
OutputIterator 
vertices_in_conflict(Weighted_point p, Cell_handle c, 
OutputIterator res); 

/*!
Similar to `find_conflicts()`, but reports the vertices which are on the 
boundary of the conflict zone of `p`, in the output iterator `res`. 
Returns the resulting output iterator. 
\pre `rt`.`dimension()` \f$ \geq2\f$, and `c` is a cell containing `p`. 

*/ 
template <class OutputIterator> 
OutputIterator 
vertices_on_conflict_zone_boundary(Weighted_point p, Cell_handle c, 
OutputIterator res); 

/*!
Similar to `find_conflicts()`, but reports the vertices which are in 
the interior of the conflict zone of `p`, in the output iterator 
`res`. The vertices that are on the boundary of the conflict zone are 
not reported. 
Returns the resulting output iterator. 
\pre `rt`.`dimension()` \f$ \geq2\f$, and `c` is a cell containing `p`. 

*/ 
template <class OutputIterator> 
OutputIterator 
vertices_inside_conflict_zone(Weighted_point p, Cell_handle c, 
OutputIterator res); 


/// @}

/*! \name
In the weighted setting, a face (cell, facet, edge or vertex) is said to be a Gabriel face iff the smallest sphere orthogonal to the weighted points associated to its vertices, has a positive power product with the weighted point of any other vertex of the triangulation. Any weighted Gabriel face belongs to the regular triangulation, but the reciprocal is not true. The following member functions test the Gabriel property of the faces of the regular triangulation.
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

/*!

*/ 
bool is_Gabriel(Vertex_handle v); 

/// @} 

/*! \name Power Diagram 
\cgal offers several functionalities to display the power diagram of a set of points in 3D. Note that the user should use a kernel with exact constructions in order to guarantee the computation of the Voronoi diagram (as opposed to computing the triangulation only, which requires only exact predicates).
*/
/// @{

/*!
Returns the weighted circumcenter of the four vertices of c. 
\pre `rt`.`dimension()`\f$ =3\f$ and `c` is not infinite. 
*/ 
Bare_point dual(Cell_handle c) const; 

/*!
Returns the dual of facet `f`, which is 

in dimension 3: either a segment, if the two cells incident to `f` 
are finite, or a ray, if one of them is infinite; 

in dimension 2: a point. 
\pre `rt`.`dimension()` \f$ \geq2\f$ and `f` is not infinite. 
*/ 
Object dual(Facet f) const; 

/*!
same as the previous method for facet `(c,i)`. 
*/ 
Object dual(Cell_handle c, int i) const; 

/*!
Sends the set of duals to all the facets of `rt` into `os`. 
*/ 
template <class Stream> Stream & draw_dual(Stream & os); 

/// @} 

/// \name Checking 
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Checks the combinatorial validity of the triangulation and the 
validity of its geometric embedding (see 
Section \ref Triangulation3secintro). Also checks that all the 
power spheres (resp. power circles in dimension 2, power segments in 
dimension 1) of cells (resp. facets in dimension 2, edges in 
dimension 1) are regular. When `verbose` 
is set to true, messages describing the first invalidity encountered 
are printed. 
This method is mainly a debugging help for the users of 
advanced features. 
\cgalDebugEnd
*/ 
bool 
is_valid(bool verbose = false) const; 

/// @}

}; /* end Regular_triangulation_3 */
} /* end namespace CGAL */
