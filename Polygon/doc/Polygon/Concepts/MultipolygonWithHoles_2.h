/*! \ingroup PkgPolygon2Concepts
 * \cgalConcept
 *
 * \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}
 *
 * A model of this concept represents a multipolygon with holes.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Multipolygon_with_holes_2<Polygon>}
 * \cgalHasModelsEnd
 */

class MultipolygonWithHoles_2 {
public:

/// \name Types
/// @{

//! the polygon type used to represent each polygon with holes of the multipolygon.
typedef unspecified_type Polygon_with_holes_2;

/*! a bidirectional iterator over the polygons with holes.
 * Its value type is `Polygon_with_holes_2`.
 */
typedef unspecified_type Polygon_with_holes_iterator;

/*! a bidirectional const iterator over the polygons with holes.
 * Its value type is `Polygon_with_holes_2`.
 */
typedef unspecified_type Polygon_with_holes_const_iterator;

//! range type for iterating over polygons with holes.
typedef unspecified_type Polygon_with_holes_container;

//! size type
typedef unsigned int Size;

/// @}

/// \name Creation
/// @{

/*! constructs a multipolygon using a range of polygons with holes.
 */
template <typename InputIterator>
MultipolygonWithHoles_2(InputIterator begin, InputIterator end);

/// @}

/// \name Predicates
/// @{

/*! returns the number of polygons with holes.
 */
Size number_of_polygons_wih_holes();

/// @}

/// \name Access Functions
/// @{

/*! returns the begin iterator of the polygons with holes.
 */
Polygon_with_holes_iterator polygons_with_holes_begin();

/*! returns the past-the-end iterator of the polygons with holes.
 */
  Polygon_with_holes_iterator polygons_with_holes_end();

/*! returns the begin iterator of the polygons with holes.
 */
Polygon_with_holes_const_iterator polygons_with_holes_begin() const;

/*! returns the past-the-end iterator of the polygons with holes.
 */
Polygon_with_holes_const_iterator polygons_with_holes_end() const;

/*! returns the range of polygons with holes.
 */
const Polygon_with_holes_container& polygons_with_holes() const;

/// @}

/// \name Modifiers
/// @{

/*! adds a given polygon with holes to the multipolygon.
 */
void add_polygon_with_holes(const Polygon_with_holes_2& polygon);

/*! erases the specified polygon.
 */
void erase_polygon_with_holes(Polygon_with_holes_const_iterator pit);

/*! removes all the polygons with holes.
 */
void clear();

/// @}

}; /* end MultipolygonWithHoles_2 */
