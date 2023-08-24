/*! \ingroup PkgPolygonRepairConcepts
 * \cgalConcept
 *
 * \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}
 *
 * A model of this concept represents a multipolygon with holes.
 * The concept requires the ability to store the polygons with holes
 * that the multipolygon is composed of.
 *
 * \cgalHasModel `CGAL::Multipolygon_with_holes_2<Polygon>`
 */

class MultipolygonWithHoles_2 {
public:

/// \name Types
/// @{

//! the polygon type used to represent each polygon with holes of the multipolygon.
typedef unspecified_type Polygon_with_holes_2;

/*! a bidirectional iterator over the polygons.
 * Its value type is `Polygon_with_holes_2`.
 */
typedef unspecified_type Polygon_const_iterator;

//! range type for iterating over holes.
typedef unspecified_type Polygons_container;

/// @}

/// \name Creation
/// @{

/*! constructs a multipolygon using a range of polygons.
 */
template <typename InputIterator>
MultipolygonWithHoles_2(InputIterator begin, InputIterator end);

/// @}

/// \name Predicates
/// @{

/*! returns the number of polygons.
 */
Size number_of_polygons();

/// @}

/// \name Access Functions
/// @{

/*! returns the begin iterator of the polygons.
 */
Polygon_const_iterator polygons_begin() const;

/*! returns the past-the-end iterator of the polygons.
 */
Polygon_const_iterator polygons_end() const;

/*! returns the range of polygons.
 */
const Polygons_container& polygons() const;

/// @}

/// \name Modifiers
/// @{

/*! adds a given polygon to the multipolygon.
 */
void add_polygon(const Polygon_with_holes_2& polygon);

/*! erases the specified polygon.
 */
void erase_polygon(Polygon_iterator pit);

/*! removes all the polygons.
 */
void clear();

/// @}

}; /* end MultipolygonWithHoles_2 */
