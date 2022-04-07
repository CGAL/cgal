
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines `GpsTraitsGeneralPolygon_2`

A model of this concept represents a simple general-polygon. The
geometric mapping of the edges of the polygon must be \f$ x\f$-monotone curves.
The concept requires the ability to access these curves.
The general polygon represented must be simple. That is, the
only points of the plane belonging to two curves are the geometric mapping
of the polygon vertices. In addition, the vertices of the represented
polygon must be ordered consistently, and the curved must be directed
accordingly. Only counterclockwise oriented polygons are valid operands
of Boolean set-operations. General polygon that represent holes must be
clockwise oriented.

\cgalHasModel `CGAL::General_polygon_2<ArrTraits>`

*/

class GeneralPolygon_2 {
public:

/// \name Types
/// @{

/*!
represents a planar (weakly) \f$ x\f$-monotone
curve. The type of the geometric mapping of the polygonal edges. It must model the concept ArrTraits::XMonotoneCurve_2.
*/
typedef unspecified_type X_monotone_curve_2;

/*!
an iterator over the geometric mapping of the
polygon edges. Its value type is `X_monotone_curve_2`.
*/
typedef unspecified_type Curve_iterator;

/*!
a const iterator over the geometric
mapping of the polygon edges. Its value type is `X_monotone_curve_2`.
*/
typedef unspecified_type Curve_const_iterator;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
GeneralPolygon_2();

/*!
copy constructor.
*/
GeneralPolygon_2(GeneralPolygon_2 other);

/*!
assignment operator.
*/
GeneralPolygon_2 operator=(other);

/*!
constructs a general polygon from a given range of curves.
*/
template <class InputIterator>
GeneralPolygon_2(InputIterator begin, InputIterator end);

/// @}

/// \name Access Functions
/// @{

/*!
returns the begin iterator of the curves.
*/
Curve_iterator curves_begin();

/*!
returns the past-the-end iterator of the curves.
*/
Curve_iterator curves_end();

/*!
returns the begin const iterator of the curves.
*/
Curve_const_iterator curves_begin();

/*!
returns the past-the-end const iterator of the curves.
*/
Curve_const_iterator curves_end();

/// @}

/// \name Modifiers
/// @{

/*!
initializes the polygon with the polygonal chain given by the range.
The value type of `Iterator` must be `X_monotone_curve_2`.
\pre The curves in the range must define a simple polygon.
*/
template <class Iterator>
void init(Iterator begin, Iterator end);

/// @}

}; /* end GeneralPolygon_2 */

