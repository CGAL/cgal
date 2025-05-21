
namespace CGAL {

/*!
\ingroup kernel_classes3

An object `b` of the class `Bbox_3` is a bounding
box in the three-dimensional Euclidean space \f$ \E^3\f$.

\cgalModels{Hashable}

\sa `CGAL::Bbox_2`

*/

class Bbox_3 {
public:

/// \name Creation
/// @{

/*!
introduces an \em empty bounding box with lower left
corner point at \f$ (\infty, \infty, \infty) \f$
and with upper right corner point at
\f$ (-\infty, -\infty, -\infty) \f$, \f$ \infty \f$ being
`std::numeric_limits<double>::%infinity()`.
*/
  Bbox_3();

/*!
introduces a bounding box `b` with lexicographically
smallest corner point at `(xmin, ymin, zmin)`
and lexicographically largest corner point at
`(xmax, ymax, zmax)`.
*/
Bbox_3(double x_min, double y_min, double z_min,
double x_max, double y_max, double z_max);

/// @}

/// \name Operations
/// @{

/*!
Test for equality.
*/
bool operator==(const Bbox_3 &c) const;

/*!
Test for inequality.
*/
bool operator!=(const Bbox_3 &q) const;

/*!
Returns 3.
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
double zmin() const;

/*!

*/
double xmax() const;

/*!

*/
double ymax() const;

/*!

*/
double zmax() const;

/*!
Returns `xmin()` if `i==0` or `ymin()` if `i==1`
or `zmin()` if `i==2`.
\pre `i>=0` and `i<=2`
*/
double min(int i) const;

/*!
Returns `xmax()` if `i==0` or `ymax()` if `i==1`
or `zmax()` if `i==2`.
\pre `i>=0` and `i<=2`
*/
double max(int i) const;

/*!
returns a bounding box of `b` and `c`.
*/
Bbox_3 operator+(const Bbox_3 &c) const;

/*!
updates `b` to be the bounding box of `b` and `c` and returns itself.
*/
Bbox_3& operator+=(const Bbox_3 &c);

/*!
dilates the bounding box by a specified number of ULP.
*/
void dilate(int dist);

/*!
scales the bounding box by `factor`, while keeping its center fixed.
\pre `factor > 0`
*/
void scale(double factor);

/// @}

}; /* end Bbox_3 */

/// \ingroup do_overlap_grp
/// @{

/*!
returns `true` iff `bb1` and `bb2` overlap, i.e., iff their
intersection is non-empty.

\relates Bbox_3
*/
bool do_overlap(const Bbox_3 &bb1, const Bbox_3 &bb2);

/// @}

/*!
returns the bounding box of the objects in the range `[first,past_end[`.
Each object in the range must have a member function `BBox_3 bbox()`
returning its bounding box.

\relates Bbox_3
*/
template<class InputIterator>
Bbox_3 bbox_3(InputIterator begin, InputIterator past_end);

/*!
returns the bounding box of the objects in the range `[first,past_end[`.
`Traits` must provide a functor `Traits::Construct_bbox_3` having an
operator returning the bounding box of each object in the range.
`Traits` must also have a member function
`Traits::Construct_bbox_3 construct_bbox_3_object() const`.

\relates Bbox_3
*/
template<class InputIterator, class Traits>
Bbox_3 bbox_3(InputIterator begin, InputIterator past_end, const Traits& traits);


} /* end namespace CGAL */
