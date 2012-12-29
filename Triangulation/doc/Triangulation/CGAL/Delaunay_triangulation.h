
namespace CGAL {

/*!
\ingroup PkgTriangulations

The class `Delaunay_triangulation` is used to maintain the full cells and vertices of a
Delaunay triangulation in \f$ \real^D\f$. It permits point insertion and
removal. The dimension \f$ D\f$ should be kept reasonably small,
see the performance section in the user manual for what reasonable
means.

Parameters
--------------

`DelaunayTriangulationTraits` is the geometric traits class that provides the geometric types
and predicates needed by Delaunay triangulations. `DelaunayTriangulationTraits` must be a model of
the concept `DelaunayTriangulationTraits`.

`TriangulationDataStructure` is the class used to store the underlying triangulation data
structure. `TriangulationDataStructure` must be a model of the concept
`TriangulationDataStructure`. The class template `Delaunay_triangulation` can
be defined by specifying only the first parameter, or by using the
tag `CGAL::Default` as
the second parameter. In both cases, `TriangulationDataStructure` defaults to
`Triangulation_data_structure<Maximal_dimension<TriangulationTraits::Point_d>::type, Triangulation_vertex<TriangulationTraits>, Triangulation_full_cell<TriangulationTraits>>`.

The class `Delaunay_triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>` inherits all the types
defined in the base class `Triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>`. Additionally, it
defines or overloads the following methods:

\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>`

*/
template< typename DelaunayTriangulationTraits, typename TriangulationDataStructure >
class Delaunay_triangulation : public Triangulation<DelaunayTriangulationTraits, TriangulationDataStructure> {
public:

/// \name Creation
/// @{

/*!
Instantiates a Delaunay triangulation with one vertex (the vertex
at infinity). See the description of the inherited nested type
`Triangulation::Maximal_dimension` for an explanation of
the use of the parameter `dim`. The complex stores a copy of the geometric
traits `gt`.
*/
Delaunay_triangulation(const int dim, const Geom_traits gt = Geom_traits());

/// @}

/// \name Point removal
/// @{

/*!
Remove the vertex `v`
from the Delaunay triangulation. If the current dimension of the triangulation has not
changed after the removal, then the returned full cell `c` geometrically
contains the removed vertex `v` (`c` can be finite or infinite).
Otherwise, the default-constructed `Full_cell_handle` is returned.
\pre `v` is a vertex of the triangulation, different from the
`infinite_vertex()`.
*/
Full_cell_handle remove(Vertex_handle v);

/*!
Remove the points or the vertices (through their
`Vertex_handle`) in the range `[start, end)`.
`*start` must be of type `Vertex_handle`.

*/
template< typename ForwardIterator > void remove(ForwardIterator
start, ForwardIterator end);

/// @}

/// \name Point insertion
/// @{

/*!
Inserts the points found in range `[s,e)` in the Delaunay triangulation
and ensures that the empty-ball property is preserved.
Returns the number of vertices actually inserted. (If more than one vertex share
the same position in space, only one insertion is counted.)
*/
template< typename ForwardIterator >
size_type insert(ForwardIterator s, ForwardIterator e);

/*!
Inserts point `p` in the Delaunay triangulation
and ensures that the empty-ball property is preserved. Returns a
`Vertex_handle` to the vertex of the triangulation with position `p`.
Prior to the actual insertion, `p` is located in the triangulation;
`hint` is used as a starting place for locating `p`.
*/
Vertex_handle insert(const Point & p, Full_cell_handle hint
= Full_cell_handle());

/*!
Same as above but uses a vertex as starting place for the search.
*/
Vertex_handle insert(const Point & p, Vertex_handle hint);

/*!
\cgalAdvanced Inserts the point `p` in the Delaunay triangulation
and ensures that the empty-ball property is preserved.
Returns a handle to the
(possibly newly created) vertex at that position. The behavior depends on the
value of `lt`:

<DL>
<DT><B>`OUTSIDE_AFFINE_HULL`</B><DD> Point
`p` is inserted so as to increase the current dimension of the Delaunay
triangulation. The method `dt`.`insert_outside_affine_hull()` is called.
<DT><B>`ON_VERTEX`</B><DD> The position of the vertex `v` described by `f`
is set to `p`. `v` is returned. <DT><B>Anything else</B><DD> The point `p`
is inserted. the full cell `c` <I>is assumed</I> to be in conflict
with `p`.
(Roughly speaking, the method `dt`.`insert_in_conflicting_cell()`
is called.)
</DL>
The parameters `lt`, `f`, `ft`
and `c` must be consistent with the localization of point `p` in the
Delaunay triangulation e.g. by a call to
`c = locate(p, lt, f, ft)`.
*/
Vertex_handle insert(const Point & p, const Locate_type lt,
const Face & f, const Facet & ft, const Full_cell_handle c);

/*!
\cgalAdvanced Inserts the point `p` in the Delaunay triangulation. Returns a handle to the
(possibly newly created) vertex at that position. \pre The point `p`
must lie outside the affine hull of the Delaunay triangulation. This implies that
`dt`.`current_dimension()` must be less that
`dt`.`maximal_dimension()`.
*/
Vertex_handle insert_outside_affine_hull(const Point & p);

/*!
\cgalAdvanced Inserts the point `p` in the Delaunay triangulation. Returns a handle to the
(possibly newly created) vertex at that position.
\pre The point `p`
must be in conflict with the full cell `c`.
*/
Vertex_handle insert_in_conflicting_cell(const Point & p, const
Full_cell_handle c);

/// @}

/// \name Queries
/// @{

/*!
Returns `true` if and only if the point `p` is in (Delaunay)
conflict with full cell `c` (i.e., the circumscribing ball of
\f$ c\f$ contains \f$ p\f$ in its interior).

*/
bool is_in_conflict(const Point & p, Full_cell_const_handle c)
const;

/*!
\cgalAdvanced Outputs handles to the full cells in confict with
point `p` into the `OutputIterator out`. The full cell `c` is used
as a starting point for gathering the full cells in conflict with
`p`.
A facet `(cc,i)` on the boundary of the conflict zone with
`cc` in conflict is returned.
\pre `c` is in conflict
with `p`.
`dt`.`current_dimension()`\f$ \geq2\f$.

*/
template< typename OutputIterator >
Facet compute_conflict_zone(const Point & p, const Full_cell_handle c,
OutputIterator out) const;

/// @}

}; /* end Delaunay_triangulation */
} /* end namespace CGAL */
