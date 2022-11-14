/*! \ingroup PkgPolygon2Concepts
 * \cgalConcept
 *
 * \cgalRefines `DefaultConstructible`
 * \cgalRefines `CopyConstructible`
 * \cgalRefines `Assignable`
 *
 * A model of this concept represents a general polygon with holes. The
 * concept requires the ability to access the general polygon that
 * represents the outer boundary and the general polygons that represent
 * the holes.
 *
 * \cgalHasModel `CGAL::General_polygon_with_holes_2<Polygon>`
 * \cgalHasModel `CGAL::Polygon_with_holes_2<Kernel,Container>`
 */

class GeneralPolygonWithHoles_2 {
public:

/// \name Types
/// @{

/*! the polygon type used to represent the outer boundary and each hole.
 */
typedef unspecified_type Polygon_2;

/*! a bidirectional iterator
 * over the polygonal holes. Its value type is `Polygon_2`.
 */
typedef unspecified_type Hole_const_iterator;


/*!
range type for iterating over holes.
*/
typedef unspecified_type Holes_container;

/// @}

/// \name Creation
/// @{

/*! constructs a general polygon with holes using a given general polygon
 * `outer` as the outer boundary and a given range of holes. If `outer` is an
 * empty general polygon, then an unbounded polygon with holes will be
 * created. The holes must be contained inside the outer boundary, and the
 * polygons representing the holes must be simple and pairwise disjoint, except
 * perhaps at the vertices.
 */
template <typename InputIterator>
GeneralPolygonWithHoles_2(Polygon_2 & outer,
                          InputIterator begin, InputIterator end);

/// @}

/// \name Predicates
/// @{

/*! returns `true` if the outer boundary is empty, and `false` otherwise.
 */
bool is_unbounded();

/// @}

/// \name Access Functions
/// @{

/*! returns the general polygon that represents the outer boundary. Note that
 * this polygon is not necessarily a valid (simple) general polygon because it
 * may be relatively simple.
 */
const Polygon_2& outer_boundary() const;

/*! returns the begin iterator of the holes.
 */
Hole_const_iterator holes_begin() const;

/*! returns the past-the-end iterator of the holes.
 */
Hole_const_iterator holes_end() const;


/*!
returns the range of holes.
*/
const Holes_container& holes() const;

/// @}

}; /* end GeneralPolygonWithHoles_2 */
