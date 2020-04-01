
namespace CGAL {

/*!
\ingroup PkgTriangulationsTriangulationClasses

This class implements triangulations of point sets in dimension \f$ d \f$.
The triangulation covers the convex hull of the input points
(the embedded vertices of the triangulation).

To store this triangulation in a triangulation data structure, we turn the set
of its faces into a topological sphere by adding a
fictitious vertex, called the <i>infinite vertex</i>, as well as infinite
simplices incident to boundary faces of the convex hull.
Each infinite \f$ i\f$-simplex is
incident to the infinite vertex and to an \f$ (i-1)\f$-simplex of the
convex hull boundary.


\tparam TriangulationTraits_ is the geometric traits class that provides the geometric types
and predicates needed by triangulations. `TriangulationTraits_` must be a model of the
concept `TriangulationTraits`.

\tparam TriangulationDataStructure_ must be a model of the concept
`TriangulationDataStructure`. This model is used to store
the faces of the triangulation. The parameter `TriangulationDataStructure_` defaults to
`Triangulation_data_structure` whose template parameters are instantiated as
follows:
<UL>
<LI>`TriangulationTraits_::Dimension`</LI>
<LI>`Triangulation_vertex<TriangulationTraits_>`</LI>
<LI>`Triangulation_full_cell<TriangulationTraits_>`.</LI>
</UL>

The triangulation deduces its maximal dimension from the type
`TriangulationTraits_::Dimension`. This dimension has to match
the dimension returned by
`TriangulationDataStructure_::maximal_dimension()`.

\cgalHeading{Input/Output}

The information in the `iostream` is: the current dimension, the number of
finite vertices, the non-combinatorial information about vertices (point,
<I>etc.</I>), the number of full cells, the indices of the vertices of each
full cell, plus the non-combinatorial information about each full cell, then the
indices of the neighbors of each full cell, where the index corresponds to the
preceding list of full cells.

\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex_, TriangulationDSFullCell_>`
\sa `Delaunay_triangulation<DelaunayTriangulationTraits_, TriangulationDataStructure_>`

*/
template< typename TriangulationTraits_, typename TriangulationDataStructure_>
class Triangulation {
public:
/// \name Types
/// @{

/*!
Type for the model of the `TriangulationTraits_` concept.
*/
typedef TriangulationTraits_ Geom_traits;

/*!
A point in Euclidean space. Note that in the context of a
`Regular_triangulation` class (which derives from this class),
`TriangulationTraits_::Point_d` is a weighted point.
*/
typedef TriangulationTraits_::Point_d Point;

/*!
This indicates whether the maximal dimension is static
(i.e.\ if the type of `Maximal_dimension` is `CGAL::Dimension_tag<int dim>`)
or dynamic (i.e.\ if the type of `Maximal_dimension` is
`CGAL::Dynamic_dimension_tag`).
In the latter case, the `dim` parameter passed to the constructor
of the class is used.
*/
typedef TriangulationTraits_::Dimension Maximal_dimension;

/*!
The second template parameter: the triangulation data structure.
*/
typedef TriangulationDataStructure_ Triangulation_ds;

/*!
A model of the concept `TriangulationVertex`.
*/
typedef TriangulationDataStructure_::Vertex Vertex;

/*!
A model of the concept
`TriangulationFullCell`.
*/
typedef TriangulationDataStructure_::Full_cell Full_cell;

/*!
The facet
class
*/
typedef TriangulationDataStructure_::Facet Facet;

/*!
A model of the concept `TriangulationDSFace`.
*/
typedef TriangulationDataStructure_::Face Face;

/// @}

/// \name Handles and Iterators
/// The vertices and full cells of triangulations are accessed through
/// handles and iterators. A handle is a model of the
/// `Handle` concept, and supports the two dereference operators:
/// `operator*` and `operator->`. Iterators are bidirectional and
/// non-mutable. They are convertible to the
/// corresponding handles, thus the user can pass them directly as
/// arguments to the functions.
/// All handles are model of `LessThanComparable` and `Hashable`,
/// that is they can be used as keys in containers such as `std::map`
/// and `boost::unordered_map`.
/// @{

/*!
handle to a a vertex
*/
typedef TriangulationDataStructure_::Vertex_handle
Vertex_handle;

/*!
const handle to a a vertex
*/
typedef TriangulationDataStructure_::Vertex_const_handle
Vertex_const_handle;

/*!
iterator over all vertices (including the infinite one)
*/
typedef TriangulationDataStructure_::Vertex_iterator
Vertex_iterator;

/*!
const iterator over all vertices (including the infinite one)
*/
typedef TriangulationDataStructure_::Vertex_const_iterator
Vertex_const_iterator;

/*!
iterator over finite vertices
*/
typedef unspecified_type Finite_vertex_iterator;

/*!
const iterator over finite vertices
*/
typedef unspecified_type Finite_vertex_const_iterator;

/*!
handle to a full cell
*/
typedef TriangulationDataStructure_::Full_cell_handle
Full_cell_handle;

/*!
const handle to a full cell
*/
typedef TriangulationDataStructure_::Full_cell_const_handle
Full_cell_const_handle;

/*!
iterator over all full cells (including the infinite ones)
*/
typedef
TriangulationDataStructure_::Full_cell_iterator
Full_cell_iterator;

/*!
const iterator over all full cells (including the infinite ones)
*/
typedef
TriangulationDataStructure_::Full_cell_const_iterator
Full_cell_const_iterator;

/*!
iterator over finite full cells
*/
typedef unspecified_type Finite_full_cell_iterator;

/*!
const iterator over finite full cells
*/
typedef unspecified_type Finite_full_cell_const_iterator;

/*!
iterator over all facets (including the infinite ones)
*/
typedef TriangulationDataStructure_::Facet_iterator
Facet_iterator;

/*!
iterator over finite facets
*/
typedef unspecified_type Finite_facet_iterator;

/*!
size type (an unsigned integral type)
*/
typedef TriangulationDataStructure_::size_type size_type;

/*!
difference
type (a signed integral type)
*/
typedef TriangulationDataStructure_::difference_type difference_type;

/*!
\enum Locate_type
\brief Used by `Triangulation` to specify which case occurs when locating a point in the triangulation.
*/
enum Locate_type { ON_VERTEX=0, /*!< when the located point coincides with a vertex of the triangulation */
                   IN_FACE, /*!< when the point is in the interior of a face of dimension equal or less than `maximal_dimension()` - 2 */
                   IN_FACET, /*!< when the point is in the interior of a facet */
                   IN_FULL_CELL, /*!< when the point is in the interior of a full cell */
                   OUTSIDE_CONVEX_HULL, /*!< when the point is outside the convex hull but in the affine hull of the current triangulation */
                   OUTSIDE_AFFINE_HULL /*!< when the point is outside the affine hull of the current triangulation. */
};

/// @}

/// \name Creation
/// @{

/*!
Instantiates a triangulation with one vertex (the vertex at infinity). See the
description of the nested type `Maximal_dimension` above for an
explanation of the use of the parameter `dim`. The triangulation stores a copy
of the geometric traits `gt`.
*/
Triangulation(int dim, const Geom_traits & gt = Geom_traits());

/*!
The copy constructor.
*/
Triangulation(const Triangulation & t2);

/// @}

/// \name Access Functions
/// @{

/*!
Returns a const reference to the underlying triangulation data structure.
*/
const Triangulation_ds & tds() const;

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
Returns a non-const
reference to the underlying triangulation data structure.
\cgalAdvancedEnd
*/
Triangulation_ds & tds();

/*!
Returns a const reference to the geometric traits instance.
*/
const Geom_traits & geom_traits() const;

/*!
Returns the maximal dimension of
the full dimensional cells that can be stored in the triangulation.
*/
int maximal_dimension() const;

/*!
Returns the dimension of the triangulation (as an embedded manifold).
*/
int current_dimension() const;

/*!
Returns `true` if the triangulation has no finite vertex. Returns
`false` otherwise.
*/
bool empty() const;

/*!
Returns the number of finite vertices in the triangulation.
*/
size_type number_of_vertices() const;

/*!
Returns the number of full cells of maximal dimension in the triangulation
(full cells incident to the vertex at infinity are counted).
*/
size_type number_of_full_cells() const;

/*!
Returns a handle to the vertex at infinity.
*/
Vertex_handle infinite_vertex() const;

/*!
Returns a handle to some full cell incident to the vertex at infinity.
*/
Full_cell_handle infinite_full_cell() const;

/// @}

/// \name Non-Constant-Time Access Functions
/// @{

/*!
Returns the number of full cells of maximal dimension that are not
incident to the vertex at infinity.
*/
size_type number_of_finite_full_cells() const;

/// @}

/// \name Tests for Finite and Infinite Elements
/// @{

/*!
Returns `true` if and only if the vertex `v` is the infinite vertex.
*/
bool is_infinite(Vertex_handle v) const;

/*!
Returns `true` if and only if `c` is incident to the infinite vertex.

*/
bool is_infinite(Full_cell_handle c) const;

/*!
Returns `true` if and only if facet `ft` is incident to the infinite
vertex.

*/
bool is_infinite(const Facet & ft) const;

/*!
Returns `true` if and
only if the face `f` is incident to the infinite vertex.

*/
bool is_infinite(const Face & f) const;

/// @}

/// \name Faces and Facets
/// @{

/*!
Returns a full cell containing the facet `f`
*/
Full_cell_handle full_cell(const Facet & f) const;

/*!
Returns the index of the vertex of the full cell
`c=tr.full_cell(f)` which does not belong to `c`.
*/
int index_of_covertex(const Facet & f) const;

/// @}

/// \name Triangulation Traversal
/// @{

/*!
The first vertex of `tr`.
*/
Vertex_iterator vertices_begin();

/*!
The beyond vertex of `tr`.
*/
Vertex_iterator vertices_end();

/*!
The first finite vertex of `tr`.
*/
Finite_vertex_iterator finite_vertices_begin();

/*!
The beyond finite vertex of `tr`.
*/
Finite_vertex_iterator finite_vertices_end();

/*!
The first full cell of `tr`.
*/
Full_cell_iterator full_cells_begin();

/*!
The beyond full cell of `tr`.
*/
Full_cell_iterator full_cells_end();

/*!
The first finite full cell of `tr`.
*/
Finite_full_cell_iterator finite_full_cells_begin();

/*!
The beyond finite full cell of `tr`.
*/
Finite_full_cell_iterator finite_full_cells_end();

/*!
Iterator to the first facet of the triangulation.
*/
Facet_iterator facets_begin();

/*!
Iterator to the beyond facet of the triangulation.
*/
Facet_iterator facets_end();

/*!
Iterator to the first finite facet of the triangulation.
*/
Finite_facet_iterator finite_facets_begin();

/*!
Iterator to the beyond finite facet of the triangulation.
*/
Finite_facet_iterator finite_facets_end();

/// @}

/// \name Point Location
/// The class `Triangulation` provides methods to locate a query point
/// with respect to the triangulation:
/// @{

/*!
The optional argument `hint` is used as a starting place for the search.

If the `query` point lies outside the affine hull of the points (which can
happen when `tr`.`current_dimension() < `
`tr`.`maximal_dimension()`) or if there is no finite vertex yet in the
triangulation, then <I>locate</I> returns a default constructed
`Full_cell_handle()`.

If the point `query` lies in the interior of a bounded (finite) full cell of `tr`,
the latter full cell is returned.

If `query` lies on the boundary of some finite full cells, one of the cells
is returned.

Let \f$ d=\f$`tr`.`current_dimension()`. If the point `query` lies
outside the convex hull of the points, an infinite full cell with vertices \f$ \{
p_1, p_2, \ldots, p_d, \infty\}\f$ is returned such that the full cell \f$ (p_1, p_2,
\ldots, p_d, query)\f$ is positively oriented (the rest of the triangulation lies
on the other side of facet \f$ (p_1, p_2, \ldots, p_d)\f$).
*/
Full_cell_handle locate(const Point & query,
  Full_cell_const_handle hint = Full_cell_handle()) const;

/*!
Same as above but `hint` is a vertex and not a full cell.
*/
Full_cell_handle locate(const Point & query, Vertex_handle hint)
const;

/*!
The optional argument `hint` is used as a starting place for the
search.
If the `query` point lies outside the affine hull of the points
(which can happen when `tr`.`current_dimension() < `
`tr`.`maximal_dimension()`) or if there is no finite vertex yet in the
triangulation, then `loc_type` is set to
`OUTSIDE_AFFINE_HULL`, and <I>locate</I> returns
`Full_cell_handle()`.
If the `query` point lies inside the affine hull
of the points, the function finds the \f$ k\f$-face that
contains `query` in its relative
interior (if the \f$ k\f$-face is finite, it is
unique) and the result is returned as follows:

<DL> <DT><B>\f$ k=0\f$</B><DD> `loc_type` is set to `ON_VERTEX`,
`f` is set to the vertex `v` the `query` lies on and a full cell
having `v` as a vertex is returned.
<DT><B>\f$ 0<k<\f$`c.current_dimension()-1`</B><DD> `loc_type` is set to
`IN_FACE`, `f` is set to the unique finite face containing the
`query` point. A full cell having `f` on its boundary is returned.
<DT><B>\f$ k=\f$`c.current_dimension()-1`</B><DD> `loc_type` is set to
`IN_FACET`, `ft` is set to one of the two representation of
the finite facet containing the
`query` point. The full cell of `ft` is returned.
<DT><B>\f$ k=\f$`c.current_dimension()`</B><DD> If the `query` point lies
<I>outside</I> the convex hull of the points in the triangulation, then
`loc_type` is set to `OUTSIDE_CONVEX_HULL` and a full cell is returned
as in the `locate` method above. If the `query` point lies
<I>inside</I> the convex hull of the points in the triangulation, then
`loc_type` is set to `IN_FULL_CELL` and the unique full cell containing
the `query` point is returned. </DL>
*/
Full_cell_handle locate(const Point & query, Locate_type & loc_type,
Face & f, Facet & ft, Full_cell_handle hint = Full_cell_handle()) const;

/*!
Same as above but `hint`, the starting place for the search, is a vertex.
The parameter `hint` is ignored if it is a default constructed
`Vertex_handle()`.
*/
Full_cell_handle
locate(const Point & query, Locate_type & loc_type,
  Face & f, Vertex_handle hint) const;

/// @}

/// \name Removal
/// @{

/*!
Contracts the `Face f` to a single vertex at
position `p`. Returns a handle to that vertex.

\pre The boundary of the union of full cells incident to `f` must be a triangulation of a
sphere of dimension `tr`.`current_dimension()`).
*/
Vertex_handle collapse_face(const Point & p, const Face & f);

/// @}

/// \name Point Insertion
/// The class `Triangulation` provides functions to insert a given point in the triangulation:
/// @{

/*!
Inserts the points found in range `[s,e)` in the triangulation. Returns
the number of vertices actually inserted. (If several vertices share the
same position in space, only the vertex that was actually inserted is counted.)
\tparam ForwardIterator must be an input iterator with the value type `Point`.
*/
template< typename ForwardIterator >
size_type insert(ForwardIterator s, ForwardIterator e);

/*!
Inserts point `p` in the triangulation. Returns a
`Vertex_handle` to the vertex of the triangulation with position `p`.
Prior to the actual insertion, `p` is located in the triangulation;
`hint` is used as a starting place for locating `p`.
*/
Vertex_handle insert(const Point &p, Full_cell_handle hint =
Full_cell_handle());

/*!
Same as above but uses a vertex `hint` as the starting place for the search.
*/
Vertex_handle insert(const Point &p, Vertex_handle hint);

/*!
Inserts point `p` into the triangulation and returns a handle to the
`Vertex` at that position. The action taken depends on the value of
`loc_type`:

<DL> <DT><B>`ON_VERTEX`</B><DD> The point of the
p`Vertex` described by `f` is set to `p`.
<DT><B>`IN_FACE`</B><DD> The point `p` is inserted in the `Face f`.
<DT><B>`IN_FACET`</B><DD> The point `p` is inserted in the `Facet ft`.
<DT><B>Anything else</B><DD> The point `p` is inserted in the triangulation according to the value
of `loc_type`, using the full cell `c`.
</DL>

This method is used internally by the other `insert()` methods.
*/
Vertex_handle insert(const Point &p, Locate_type loc_type, Face & f, Facet & ft, Full_cell_handle c);

/*!
Removes the full cells in the range \f$ C=\f$`[s, e)`, inserts a vertex
at position `p` and fills the hole by connecting
each face of the boundary to `p`.
A `Vertex_handle` to the new `Vertex` is
returned. The facet `ft` must lie on the boundary of \f$ C\f$ and its
defining full cell, `tr`.`full_cell(ft)` must lie inside \f$ C\f$. Handles
to the newly created full cells are output in the `out` output iterator.
\pre \f$C\f$ must be a topological ball, must contain `p` in its
interior and must not contain any vertex of the triangulation.
*/
template < typename ForwardIterator, typename OutputIterator >
Vertex_handle insert_in_hole(const Point & p, ForwardIterator s,
ForwardIterator e, const Facet & ft, OutputIterator out);

/*!
Same as above, but the newly created full cells are not
retrieved.
*/
template < typename ForwardIterator > Vertex_handle
insert_in_hole(const Point & p, ForwardIterator s, ForwardIterator e, const
Facet & ft);

/*!
Inserts point `p` in the triangulation.
\pre `p` must lie in the relative interior of `f`.
*/
Vertex_handle insert_in_face(const Point & p, const Face & f);

/*!
Inserts point `p` in the triangulation.
\pre `p` must lie in the relative interior of `ft`.
*/
Vertex_handle insert_in_facet(const Point & p, const Facet & ft);

/*!
Inserts point `p` in the triangulation. \pre `p` must lie in the
interior of `c`.
*/
Vertex_handle insert_in_full_cell(const Point & p, Full_cell_handle
c);

/*!
Inserts point `p` in the triangulation.
\pre `p` must lie outside the convex hull of `tr`. The half-space
defined by the infinite full cell `c` must contain `p`.
*/
Vertex_handle insert_outside_convex_hull(const Point &,
Full_cell_handle c);

/*!
Inserts point `p` in the triangulation.
\pre `p` must lie outside the affine hull of `tr`.
*/
Vertex_handle insert_outside_affine_hull(const Point &);

/// @}

/// \name Validity Check
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Partially checks whether `tr` is a triangulation. This function returns
`true` if the combinatorial triangulation data structure's `is_valid()`
test returns `true` and if some geometric tests are passed with success. It
is checked that the orientation of each finite full cell is positive and that
the orientation of each infinite full cell is consistent with their finite
adjacent full cells.
The `verbose` parameter is not used.
\cgalDebugEnd
*/
bool is_valid(bool verbose=false) const;

/*!
\cgalDebugFunction
\cgalDebugBegin
Returns `true` if and only if all
finite full cells incident to `v` have positive orientation.
The `verbose` parameter is not used.
\cgalDebugEnd
*/
bool are_incident_full_cells_valid(Vertex_const_handle v, bool
verbose = false) const;

/// @}

/// \name Input/Output
/// @{

/*!
Reads the underlying combinatorial triangulation from `is` by
calling the corresponding input operator of the triangulation data
structure class (note that the infinite vertex is numbered 0), and the
non-combinatorial information by calling the corresponding input
operators of the vertex and the full cell classes (such as point
coordinates), which are provided by overloading the stream operators
of the vertex and full cell types. Assigns the resulting triangulation to
`t`.
*/
std::istream & operator>> (std::istream & is, Triangulation & t);

/*!
Writes the triangulation `t` into `os`.
*/
std::ostream& operator<< (std::ostream& os, const Triangulation & t);

/// @}

}; /* end Triangulation */
} /* end namespace CGAL */
