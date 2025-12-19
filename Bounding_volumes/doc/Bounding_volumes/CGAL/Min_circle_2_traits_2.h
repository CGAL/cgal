
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class `Min_circle_2_traits_2` is a traits class for `Min_circle_2<Traits>`
using the two-dimensional \cgal kernel.

\tparam K must be a model for `Kernel`.

\cgalModels{MinCircle2Traits}

\sa `CGAL::Min_circle_2<Traits>`
\sa `MinCircle2Traits`

*/
template< typename K >
class Min_circle_2_traits_2 {
public:

/// \name Types
/// @{

/*!
typedef to `K::Point_2`.
*/
typedef unspecified_type Point ;

/*!
internal type.
*/
typedef unspecified_type Circle;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Min_circle_2_traits_2( );

/*!
copy constructor.
*/
Min_circle_2_traits_2(
const Min_circle_2_traits_2<K>&);

/// @}

/// \name Operations
/// @{

/*!

returns `CGAL::orientation( p, q, r)`.
*/
CGAL::Orientation
orientation( const Point& p,
const Point& q,
const Point& r) const;

/// @}

}; /* end Min_circle_2_traits_2 */
} /* end namespace CGAL */
