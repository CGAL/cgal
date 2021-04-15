
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

The class `Fixed_alpha_shape_3` represents one (fixed)
alpha shape of points in the 3D space for a real
\f$ \alpha\f$. It maintains an underlying triangulation
of the class `Dt` which
represents connectivity and order among its faces. Each
\f$ k\f$-dimensional face of the `Dt` is associated with
a classification that specifies its status in the alpha complex, alpha being fixed.

\tparam Dt must be either `Delaunay_triangulation_3`, `Regular_triangulation_3`,
`Periodic_3_Delaunay_triangulation_3` or `Periodic_3_regular_triangulation_3`.
Note that `Dt::Geom_traits`, `Dt::Vertex`, and `Dt::Face`
must be model the concepts `AlphaShapeTraits_3`,
`AlphaShapeVertex_3` and `AlphaShapeFace_3`, respectively.

Note that this class is used for <I>basic</I>, <I>weighted</I>,
and <I>periodic</I>Alpha Shapes.

The modifying functions `insert` and `remove` will overwrite
the one inherited from the underlying triangulation class `Dt`.
At the moment, only the static version is implemented.

\cgalHeading{I/O}

The I/O operators are defined for `iostream`, and for
the window stream provided by \cgal. The format for the iostream
is an internal format.

*/
template< typename Dt >
class Fixed_alpha_shape_3 : public Dt {
public:

/// \name Types
/// @{

/*!
the alpha shape traits type. It has to derive from a triangulation
traits class. For example `Dt::Point` is a Point class.
*/
typedef unspecified_type Gt;

/*!
the number type of alpha.
*/
typedef Gt::FT FT;

/*!
Enum to classify the simplices of the underlying
triangulation with respect to a given alpha value.

Each k-dimensional simplex of the triangulation
can be classified as `EXTERIOR`, `SINGULAR`, `REGULAR`
or `INTERIOR`.
A \f$ k\f$ simplex is `REGULAR` if it is on the boundary
of the alpha complex and belongs to a \f$ k+1\f$ simplex in this complex
and it is `SINGULAR` if it is a boundary simplex that is not included in a \f$ k+1\f$ simplex of the complex.

*/
enum Classification_type {EXTERIOR, SINGULAR, REGULAR, INTERIOR};

/// @}

/// \name Creation
/// @{

/*!
builds an empty fixed alpha shape and sets the alpha value to `alpha`.
*/
Fixed_alpha_shape_3(FT alpha = 0);

/*!
builds a fixed alpha shape from the triangulation `dt`,
and sets the alpha value to `alpha`.
\attention This operation swaps `*this` and `dt`, that is `dt` is an empty triangulation once the fixed alpha shape is built.
*/
Fixed_alpha_shape_3(Dt& dt,FT alpha = 0);

/*!
builds a fixed alpha shape for the points in the range
`[first,last)` and sets the alpha value to `alpha`.
\tparam InputIterator must be an input iterator with value type `Point` (the type point of the underlying triangulation.)
*/
template < class InputIterator >
Fixed_alpha_shape_3(
InputIterator first,
InputIterator last,
const FT& alpha = 0);

/// @}

/// \name Modifiers
/// @{

/*!

inserts the point `p` in the underlying triangulation and returns the corresponding vertex.
The optional argument `start` is used as a starting place for the search.
The classification types of the new simplices are computed and that of the simplices incident
to the new ones are updated.

*/
Vertex_handle insert (Point p,Cell_handle start = Cell_handle());

/*!

removes the vertex `v` from the underlying triangulation.
The classification types of new simplices and their incident faces are set or reset.

*/
void remove (Vertex_handle v);

/*!
clears the structure.
*/
void
clear();

/// @}

/// \name Query Functions
/// @{

/*!
returns the \f$ \alpha\f$-value.
*/
const FT&
get_alpha(void) const;

/*!
classifies the cell `c` of the underlying triangulation in the alpha complex.
*/
Classification_type
classify(Cell_handle c) const;

/*!
classifies the facet `f` of the underlying triangulation in the alpha complex.
*/
Classification_type classify(Facet f) const;

/*!
classifies the facet of the cell `f` opposite to the vertex with index
`i` of the underlying triangulation in the alpha complex.
*/
Classification_type classify(Cell_handle f, int i) const;

/*!
classifies the edge `e` of the underlying triangulation in the alpha complex.
*/
Classification_type classify(const Edge& e) const;

/*!
classifies the vertex `v` of the underlying triangulation in the alpha complex.
*/
Classification_type classify(Vertex_handle v) const;

/*!
writes the cells which are of type `type` in the alpha complex
to the sequence
pointed to by the output iterator `it`. Returns past the end
of the output sequence.
*/
template<class OutputIterator>
OutputIterator get_alpha_shape_cells(OutputIterator it, Classification_type type);

/*!
writes the facets which are of type `type` in the alpha complex
to the sequence pointed to by the output iterator `it`. Returns past the end
of the output sequence.
*/
template<class OutputIterator>
OutputIterator get_alpha_shape_facets(OutputIterator it, Classification_type type);

/*!
writes the edges which are of type `type` in the alpha complex
to the sequence
pointed to by the output iterator `it`. Returns past the end
of the output sequence.
*/
template<class OutputIterator>
OutputIterator get_alpha_shape_edges(OutputIterator it, Classification_type type);

/*!
writes the vertices which are of type `type` in the alpha complex
to the sequence pointed to by the output iterator `it`. Returns past the end
of the output sequence.
*/
template<class OutputIterator>
OutputIterator get_alpha_shape_vertices(OutputIterator it, Classification_type type);

/// @}

}; /* end Fixed_alpha_shape_3 */

/*!
inserts the fixed alpha shape `A` into the stream `os`.


An overlaoad of `operator<<` must be available for `GT::Point`.
\relates Fixed_alpha_shape_3
*/
  std::ostream& operator<<(std::ostream& os, const Fixed_alpha_shape_3<Dt>& A);


} /* end namespace CGAL */
