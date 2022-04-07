
namespace CGAL {

/*!
\ingroup kernel_classes2

An object `b` of the class `Bbox_2` is a bounding
box in the two-dimensional Euclidean plane \f$ \E^2\f$. This class is not templated.

\cgalModels `Hashable`

\sa `CGAL::Bbox_3`

*/

class Bbox_2 {
public:

/// \name Creation
/// @{

/*!
introduces an \em empty bounding box with lower left
corner point at \f$ (\infty, \infty) \f$
and with upper right corner point at
\f$ (-\infty, -\infty) \f$, \f$ \infty \f$ being
`std::numeric_limits<double>::%infinity()`.
*/
  Bbox_2();
/*!
introduces a bounding box `b` with lower left corner at
`(xmin, ymin)` and with upper right corner at
`(xmax, ymax)`.
*/
Bbox_2(double x_min, double y_min,
double x_max, double y_max);

/// @}

/// \name Operations
/// @{

/*!
Test for equality.
*/
bool operator==(const Bbox_2 &c) const;

/*!
Test for inequality.
*/
bool operator!=(const Bbox_2 &q) const;

/*!
Returns 2.
*/
int dimension() const;

/*!

*/
double xmin() const;

/*!

*/
double ymin() const;

/*!

*/
double xmax() const;

/*!

*/
double ymax() const;

/*!
Returns `xmin()` if `i==0` or `ymin()` if `i==1`.
\pre i==0 or i==1
*/
double min(int i) const;

/*!
Returns `xmax()` if `i==0` or `ymax()` if `i==1`.
\pre i==0 or i==1
*/
double max(int i) const;

/*!
returns a bounding box of `b` and `c`.
*/
Bbox_2 operator+(const Bbox_2 &c) const;

/*!
updates `b` to be the bounding box of `b` and `c` and returns itself.
*/
Bbox_2& operator+=(const Bbox_2 &c);

/*!
dilates the bounding box by a specified number of ULP.
*/
void dilate(int dist);

/// @}

}; /* end Bbox_2 */

/// \ingroup do_overlap_grp
/// @{

/*!
returns `true` iff `bb1` and `bb2` overlap, i.e., iff their
intersection is non-empty.

\relates Bbox_2
*/
bool do_overlap(const Bbox_2 &bb1, const Bbox_2 &bb2);

/// @}

/*!
returns the bounding box of the objects in the range `[first,past_end[`.
Each object in the range must have a member function `BBox_2 bbox()`
returning its bounding box.

\relates Bbox_2
*/
template<class InputIterator>
Bbox_2 bbox_2(InputIterator begin, InputIterator past_end);

/*!
returns the bounding box of the objects in the range `[first,past_end[`.
`Traits` must provide a functor `Traits::Construct_bbox_2` having an
operator returning the bounding box of each object in the range.
`Traits` must also have a member function
`Traits::Construct_bbox_2 construct_bbox_2_object() const`.

\relates Bbox_2
*/
template<class InputIterator, class Traits>
Bbox_2 bbox_2(InputIterator begin, InputIterator past_end, const Traits& traits);


} /* end namespace CGAL */
