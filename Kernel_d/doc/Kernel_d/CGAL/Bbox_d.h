
namespace CGAL {

/*!
\ingroup PkgKernelDKernelObjs

An object `b` of the class `Bbox_d` is a bounding
box in the d-dimensional Euclidean plane \f$ \E^d\f$. This class is templated with a dimension tag.

\cgalModels{Hashable}

\sa `CGAL::Bbox_2`
\sa `CGAL::Bbox_3`

*/

template <typename DimensionTag>
class Bbox_d {
public:

/// \name Creation
/// @{

/*!
introduces an \em empty bounding box with lower left
corner coordinates at \f$ \infty \f$
and with upper right corner coordinates at
\f$ -\infty \f$, \f$ \infty \f$ being
`std::numeric_limits<double>::%infinity()`.
*/
  Bbox_d();

/*!
introduces a d-dimensional bounding box from a 2d bounding box.
\pre the dimension must be 2D
*/
Bbox_d(const Bbox_2& b);

/*!
introduces a d-dimensional bounding box from a range.
\pre the range must have the size of the dimension.
\tparam I an iterator model of `InputIterator` with value type double
*/
template <typename I>
Bbox_d(int d, I b, I e);

/*!
introduces a d-dimensional bounding box from a 3d bounding box.
\pre the dimension must be 3D
*/
Bbox_d(const Bbox_3& b);

/// @}

/// \name Operations
/// @{

/*!
ests for equality.
*/
bool operator==(const Bbox_d &c) const;

/*!
tests for inequality.
*/
bool operator!=(const Bbox_d &q) const;

/*!
returns the dimension.
*/
int dimension() const;


/*!
returns an iterator for the %Cartesian coordinates of the lower left and the upper right corner.
*/
Cartesian_const_iterator cartesian_begin() const;


/*!
returns the past-the-end iterator for the %Cartesian coordinates of the lower left and the upper right corner.
*/
Cartesian_const_iterator cartesian_begin() const;

/*!
returns a bounding box of `b` and `c`.
*/
Bbox_d operator+(const Bbox_d &c) const;

/*!
updates `b` to be the bounding box of `b` and `c` and returns itself.
*/
Bbox_d& operator+=(const Bbox_d &c);

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

}; /* end Bbox_d */

/// \ingroup do_overlap_grp
/// @{

/*!
returns `true`, iff `bb1` and `bb2` overlap, i.e., iff their
intersection is non-empty.

\relates Bbox_d
*/
template <typename DimensionTag>
bool do_overlap(const Bbox_d &bb1, const Bbox_d &bb2);

/// @}


} /* end namespace CGAL */
