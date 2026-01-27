
namespace CGAL {

/*!
\ingroup PkgBoundingVolumesRef

The class `Min_ellipse_2_traits_2` is a traits class for `CGAL::Min_ellipse_2<Traits>`
using the two-dimensional \cgal kernel.

The template parameter `K` must be a model for `Kernel`.

\cgalModels{MinEllipse2Traits}

\sa `CGAL::Min_ellipse_2<Traits>`
\sa `MinEllipse2Traits`

*/
template< typename K >
class Min_ellipse_2_traits_2 {
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
typedef unspecified_type Ellipse;

/// @}

/// \name Access Functions
/// The Ellipse type provides the following access methods not
/// required by the concept `MinEllipse2Traits`.
/// @{

/*!
tests whether the ellipse is a circle.
*/
bool is_circle();

/*!
gives a double approximation of the
ellipse's conic equation. If `K` is a %Cartesian kernel, the ellipse
is the set of all points \f$ (x,y)\f$ satisfying \f$ rx^2+sy^2+txy+ux+vy+w=0\f$. In the
Homogeneous case, the ellipse is the set of points \f$ (hx,hy,hw)\f$ satisfying
\f$ r(hx)^2+s(hy)^2+t(hx)(hy)+u(hx)(hw)+v(hy)(hw)+w(hw)^2=0\f$.
*/
void double_coefficients (double &r, double &s, double &t, double &u, double &v, double &w);

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Min_ellipse_2_traits_2( );

/*!
copy constructor.
*/
Min_ellipse_2_traits_2(
const Min_ellipse_2_traits_2<K>&);

/// @}

}; /* end Min_ellipse_2_traits_2 */
} /* end namespace CGAL */
