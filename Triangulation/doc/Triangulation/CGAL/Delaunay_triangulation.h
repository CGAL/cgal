
namespace CGAL {

/*!
\ingroup PkgTriangulationsTriangulationClasses

This class is used to maintain the
Delaunay triangulation of a set of points in \f$ \mathbb{R}^D \f$.
It permits point insertion and
removal. The dimension \f$ D\f$ can be specified at compile-time or
run-time. It should be kept reasonably small,
see the performance section in the user manual for what reasonable
means.

In a Delaunay triangulation, each face has the so-called
<I>Delaunay</I> or <I>empty-ball</I> property: there exists a
circumscribing ball whose interior does not contain 
any vertex of the triangulation.
A <I>circumscribing ball</I> of a simplex is a ball
having all vertices of the simplex on its boundary.


\tparam DelaunayTriangulationTraits is the geometric traits class that provides the geometric types
and predicates needed by Delaunay triangulations. `DelaunayTriangulationTraits` must be a model of
the concept `DelaunayTriangulationTraits`.

\tparam TriangulationDataStructure must be a model of the concept
`TriangulationDataStructure`. This model is used to store 
the faces of the triangulation. The parameter `TriangulationDataStructure` defaults to
`Triangulation_data_structure` whose template parameters are instantiated as
follows:
<UL>
<LI>`DelaunayTriangulationTraits::Dimension`</LI>
<LI>`Triangulation_vertex<DelaunayTriangulationTraits>`</LI>
<LI>`Triangulation_full_cell<DelaunayTriangulationTraits>`.</LI>
</UL>

The class template `Delaunay_triangulation` can
be defined by specifying only the first parameter, or by using the
tag `CGAL::Default` as the second parameter. 

The class `Delaunay_triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>` inherits all the types
defined in the base class `Triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>`. Additionally, it
defines or overloads the following methods:

\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>`

*/
template< typename DelaunayTriangulationTraits, typename TriangulationDataStructure >
class Delaunay_triangulation
  : public Triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>
{
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

/// \name Point Removal
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
Remove the vertices pointed by the vertex handles in the range `[start, end)`.
\tparam ForwardIterator must be an input iterator with the value type `Vertex_handle`.
*/
template< typename ForwardIterator > 
void remove(ForwardIterator start, ForwardIterator end);

/// @}

/// \name Point Insertion
/// @{

/*!
Inserts the points found in range `[s,e)` in the Delaunay triangulation
and ensures that the empty-ball property is preserved.
Returns the number of vertices actually inserted. (If more than one vertex share
the same position in space, only one insertion is counted.)
\tparam ForwardIterator must be an input iterator with the value type `Point`. 
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
\cgalAdvancedBegin
Inserts the point `p` in the Delaunay triangulation
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
\cgalAdvancedEnd
*/
Vertex_handle insert(const Point & p, const Locate_type lt,
const Face & f, const Facet & ft, const Full_cell_handle c);

/*!
\cgalAdvancedBegin
Inserts the point `p` in the Delaunay triangulation. Returns a handle to the
(possibly newly created) vertex at that position. 
\pre The point `p`
must lie outside the affine hull of the Delaunay triangulation. This implies that
`dt`.`current_dimension()` must be less that
`dt`.`maximal_dimension()`.
\cgalAdvancedEnd
*/
Vertex_handle insert_outside_affine_hull(const Point & p);

/*!
\cgalAdvancedBegin
Inserts the point `p` in the Delaunay triangulation. Returns a handle to the
(possibly newly created) vertex at that position.
\pre The point `p` must be in conflict with the full cell `c`.
\cgalAdvancedEnd
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
\cgalAdvancedBegin
Outputs handles to the full cells in confict with
point `p` into the `OutputIterator out`. The full cell `c` is used
as a starting point for gathering the full cells in conflict with
`p`.
A facet `(cc,i)` on the boundary of the conflict zone with
`cc` in conflict is returned.
\pre `c` is in conflict
with `p`.
`dt`.`current_dimension()`\f$ \geq2\f$.
\cgalAdvancedEnd
*/
template< typename OutputIterator >
Facet compute_conflict_zone(const Point & p, const Full_cell_handle c,
OutputIterator out) const;

/// @}

}; /* end Delaunay_triangulation */
} /* end namespace CGAL */
