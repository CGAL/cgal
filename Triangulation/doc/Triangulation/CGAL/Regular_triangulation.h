
namespace CGAL {

/*!
\ingroup PkgTriangulationsTriangulationClasses

This class is used to maintain the
regular triangulation -- also known as weighted Delaunay triangulation --
of a set of weighted points in \f$ \mathbb{R}^D \f$.
The maximal dimension \f$ D\f$ can be specified at compile-time or
run-time. It should be kept reasonably small -- see the performance
section in the user manual for what reasonable means.

\warning The removal of points is not supported yet.

\tparam RegularTriangulationTraits_ is the geometric traits class that provides the
geometric types and predicates needed by regular triangulations.
`RegularTriangulationTraits_` must be a model of the concept
`RegularTriangulationTraits`.

\tparam TriangulationDataStructure_ must be a model of the concept
`TriangulationDataStructure`. This model is used to store
the faces of the triangulation. The parameter `TriangulationDataStructure_`
defaults to `Triangulation_data_structure` whose template parameters are
instantiated as follows:
<UL>
<LI>`RegularTriangulationTraits_::Dimension`</LI>
<LI>`Triangulation_vertex<CGAL::Regular_triangulation_traits_adapter<RegularTriangulationTraits_> >`</LI>
<LI>`Triangulation_full_cell<CGAL::Regular_triangulation_traits_adapter<RegularTriangulationTraits_> >`.</LI>
</UL>

`Regular_triangulation` can
be defined by specifying only the first parameter, or by using the
tag `CGAL::Default` as the second parameter.

\sa `Delaunay_triangulation`
\sa `Triangulation_data_structure`
\sa `Regular_triangulation_traits_adapter`

*/
template< typename RegularTriangulationTraits_, typename TriangulationDataStructure_ >
class Regular_triangulation
  : public Triangulation<Regular_triangulation_traits_adapter<RegularTriangulationTraits_>, TriangulationDataStructure_>
{
public:

/// \name Types
/// @{

/*!
A point in Euclidean space with an associated weight.
*/
typedef RegularTriangulationTraits_::Weighted_point_d Weighted_point;

/// @}

/// \name Creation
/// @{

/*!
Instantiates a regular triangulation with one vertex (the vertex
at infinity). See the description of the inherited nested type
`Triangulation::Maximal_dimension` for an explanation of
the use of the parameter `dim`. The triangulation stores a copy of the
geometric traits `gt`.
*/
Regular_triangulation(int dim, const Geom_traits &gt = Geom_traits());

/// @}

/// \name Point Insertion
/// @{

/*!
Inserts weighted point `p` in the triangulation and returns the corresponding
vertex.

If this insertion creates a vertex, this vertex is returned.

If `p` coincides with an existing vertex and has a greater weight,
then the existing weighted point becomes hidden and `p` replaces it as vertex
of the triangulation.

If `p` coincides with an already existing vertex (both point and
weights being equal), then this vertex is returned and the triangulation
remains unchanged.

Otherwise if `p` does not appear as a vertex of the triangulation,
then it is stored as a hidden point and this method returns the default
constructed handle.

Prior to the actual insertion, `p` is located in the triangulation;
`hint` is used as a starting place for locating `p`.
*/
Vertex_handle insert(const Weighted_point & p, Full_cell_handle hint
  = Full_cell_handle());

/*!
Same as above but uses a vertex as starting place for the search.
*/
Vertex_handle insert(const Weighted_point & p, Vertex_handle hint);

/*!
Inserts weighted point `p` in the triangulation.
Similar to the above `insert()` function, but takes as additional
parameter the return values of a previous location query. See description of
`Triangulation::locate()`.
*/
Vertex_handle insert(const Weighted_point & p, Locate_type lt,
  const Face & f, const Facet & ft, Full_cell_handle c);

/*!
Inserts the weighted points found in range `[s,e)` in the regular triangulation.
Returns the difference between the number of vertices after and before
the insertions (it may be negative due to hidden points).
Note that this function is not guaranteed to insert the points
following the order of `ForwardIterator` because `spatial_sort()`
is used to improve efficiency.

\tparam ForwardIterator must be an input iterator with the value
type `Weighted_point`.
*/
template< typename ForwardIterator >
std::ptrdiff_t insert(ForwardIterator s, ForwardIterator e);

/*!
inserts the weighted point `p` in the triangulation on the
condition that the vertex `star_center` appears
in the conflict zone of `p` (that is, `p` would appear
in the star of `star_center` after the insertion).

If the insertion of `p` creates a new vertex, this vertex is returned.
Otherwise, a default constructed handle is returned.

If the dimension of the triangulation is `0` and if `p` coincides
with an existing vertex and has a greater weight,
then `p` replaces it as vertex of the triangulation
and this vertex is returned.

If `p` coincides with an already existing vertex, with both point and
weights being equal, then this vertex is returned and the triangulation
remains unchanged.

By convention, if the insertion of `p` would increase the dimension
of the triangulation, it is inserted.

Prior to the actual insertion, `p` is located in the triangulation;
`start` is used as a starting place for locating `p`.

\sa `Regular_triangulation::compute_conflict_zone()`
*/
Vertex_handle insert_if_in_star(const Weighted_point& p, Vertex_handle star_center,
                                Full_cell_handle start = Full_cell_handle());

/*!
Same as the above `insert_if_in_star()`, but uses a vertex as starting place for the search.
*/
Vertex_handle insert_if_in_star(const Weighted_point& p, Vertex_handle star_center,
                                Vertex_handle hint);

/*!
inserts the weighted point `p` in the triangulation.

This function is similar to the above `insert_if_in_star()` function,
but takes as additional parameters the return values of the location query
of `p`.

\sa `Triangulation::locate()`.
*/
Vertex_handle insert_if_in_star(const Weighted_point& p, Vertex_handle star_center, Locate_type lt,
                                const Face& f, const Facet& ft, Full_cell_handle s);

/// @}

/// \name Queries
/// @{

/*!
Returns `true` if and only if the point `p` is in
conflict with full cell `c` (A weighted point `p` is said to be in conflict
with a cell `c` if it has a negative power distance to the power sphere of `c`.)
*/
bool is_in_conflict(const Weighted_point & p, Full_cell_const_handle c)
const;

/*!
Outputs handles to the full cells in conflict with
point `p` into the `OutputIterator out`. The full cell `c` is used
as a starting point for gathering the full cells in conflict with
`p`.
A facet `(cc,i)` on the boundary of the conflict zone with
`cc` in conflict is returned.
\pre `c` is in conflict with `p` and `rt`.`current_dimension()`\f$ \geq 1\f$.
*/
template< typename OutputIterator >
Facet compute_conflict_zone(const Weighted_point & p, Full_cell_handle c,
OutputIterator out) const;

/// @}

/// \name Access Functions
/// @{

/*!
Returns the number of finite vertices that are not hidden.
*/
size_type number_of_vertices() const;

/*!
Returns the number of hidden vertices.
*/
size_type number_of_hidden_vertices() const;

/// @}


}; /* end regular_triangulation */
} /* end namespace CGAL */
