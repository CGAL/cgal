
namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2Ref

The class `General_polygon_2` models the concept `GeneralPolygon_2`.
It represents a simple general-polygon. It is parameterized with the type
`ArrTraits` that models the concept
`ArrangementDirectionalXMonotoneTraits_2`. The latter is a refinement
of the concept `ArrangementXMonotoneTraits_2`. In addition to the
requirements of the concept `ArrangementXMonotoneTraits_2`, a
model of the concept `ArrangementDirectionalXMonotoneTraits_2` must
support the following functions:
<UL>
<LI>Given an \f$ x\f$-monotone curve, construct its opposite curve.
<LI>Given an \f$ x\f$-monotone curve, compare its two endpoints lexicographically.
</UL>

This class supports a few convenient operations in addition to the
requirements that the concept `GeneralPolygon_2` lists.

\cgalModels `GeneralPolygon_2`

*/
template< typename ArrTraits >
class General_polygon_2 {
public:

/// \name Types
/// @{

/*!
number of edges size type.
*/
typedef unspecified_type Size;

/*!
a general planar curve.
*/
typedef typename ArrTraits::Curve_2 Curve_2;

/// @}

/// \name Construction
/// @{

/*! constructs a `General_polygon_2` from a \f$X\f$-monotone curves.
 */
template <class CurveIterator>
General_polygon_2(CurveIterator begin, CurveIterator end);

/// @}

/// \name Operations
/// @{

/*!
returns the number of edges of the general polygon.
*/
Size size();

/// @}

/// \name Modifiers
/// @{

/*!
clears the polygon.
*/
void clear();

/*!
reverses the orientation of the polygon.
\pre `is_simple()`.
*/
void reverse_orientation();

/// @}

/// \name Predicates
/// @{

/*!
returns `true` if the polygon is empty, and `false` otherwise.
*/
bool is_empty();

/*!
returns the orientation of the polygon.
\pre `is_simple()`.
*/
Orientation orientation();

/// @}

}; /* end General_polygon_2 */

/*!
This operator imports a general polygon from the input stream `in`.

Both \ascii and binary formats are supported, and the format is automatically detected.

The format consists of the number of points of the outer boundary followed
by the points themselves in counterclockwise order, followed by the number of holes,
and for each hole, the number of points of the outer boundary is followed
by the points themselves in clockwise order.

\relates General_polygon_2
*/
template <class ArrTraits>
std::istream& operator>>(std::istream& in, CGAL::General_polygon_2<ArrTraits>& P);


/*!
This operator exports a general polygon to the output stream `out`.

An \ascii and a binary format exist. The format can be selected with
the \cgal modifiers for streams, `set_ascii_mode` and `set_binary_mode`
respectively. The modifier `set_pretty_mode` can be used to allow for (a
few) structuring comments in the output. Otherwise, the output would
be free of comments. The default for writing is \ascii without
comments.

The number of curves of the outer boundary is exported followed by the
curves themselves in counterclockwise order.

\relates General_polygon_2
*/
template <class ArrTraits>
std::ostream& operator<<(std::ostream& out, CGAL::General_polygon_2<ArrTraits>& P);



} /* end namespace CGAL */
