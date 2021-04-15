
namespace CGAL {

/*!
  \ingroup PkgArrangementOnSurface2TraitsClasses

  The traits class `Arr_linear_traits_2` is a model of the
  `ArrangementTraits_2` concept, which enables the construction and maintenance of arrangements of
  linear objects. The linear objects may be bounded (line segments) or
  unbounded (rays and lines). Thus, it is also a model of the concept
  `ArrangementOpenBoundaryTraits_2`. The traits class is parameterized
  parameterized with a \cgal-kernel model; see the reference page of
  `Arr_segment_traits_2<Kernel>` for further explanations and
  recommendations on choosing a kernel.

  `Arr_linear_traits_2` defines `Kernel::Point_2` as its point type. The nested
  `X_monotone_curve_2` and `Curve_2` types defined by the traits class
  (as is the case with the various segment-traits classes, both types
  refer to the same class, as <I>every</I> linear object is (weakly) \f$ x\f$-monotone),
  are constructible from a point, a line segment, a ray and from a line
  (objects of types `Kernel::Point_2`, `Kernel::Segment_2`,
  `Kernel::Ray_2` and `Kernel::Line_2`, respectively). On the other
  hand, when we are given a curve we can find out its actual type and convert
  it to the respective kernel object (say, to a `Kernel::Ray_2`).

  \cgalModels `ArrangementTraits_2`
  \cgalModels `ArrangementLandmarkTraits_2`
  \cgalModels `ArrangementOpenBoundaryTraits_2`
*/
template< typename Kernel >
class Arr_linear_traits_2 {
public:

  /*!
    The `Curve_2` (and the `X_monotone_curve_2`) class nested within
    the linear-traits can represent all types of linear objects. The class
    supports the copy and default constructor and the assignment operator.
    In addition, the \link PkgArrangementOnSurface2op_left_shift `operator<<` \endlink and \link PkgArrangementOnSurface2op_right_shift `operator>>` \endlink for linear objects
    are defined for standard output and input streams.
  */
  class Curve_2 {
  public:

    /// \name Types
    /// @{

    /*!

     */
    typedef typename Kernel::Point_2 Point_2;

    /*!

     */
    typedef typename Kernel::Segment_2 Segment_2;

    /*!

     */
    typedef typename Kernel::Ray_2 Ray_2;

    /*!

     */
    typedef typename Kernel::Line_2 Line_2;

    /// @}

    /// \name Creation
    /// @{

    /*!
      constructs an line segment.
    */
    Curve_2 (const Segment_2& seg);

    /*!
      constructs a ray.
    */
    Curve_2 (const Ray_2& ray);

    /*!
      constructs a line.
    */
    Curve_2 (const Line_2& ln);

    /// @}

    /// \name Access Functions
    /// @{

    /*!
      returns whether `c` is a segment (has two bounded endpoints).
    */
    bool is_segment() const;

    /*!
      returns a segment object equivalent to `c`.
      \pre `c` is a segment.
    */
    Segment_2 segment() const;

    /*!
      returns whether `c` is a ray (has a bounded source point and its
      other end is unbounded).
    */
    bool is_ray() const;

    /*!
      returns a ray object equivalent to `c`.
      \pre `c` is a ray.
    */
    Ray_2 ray() const;

    /*!
      returns whether `c` is a line (has only unbounded ends).
    */
    bool is_line() const;

    /*!
      returns a line object equivalent to `c`.
      \pre `c` is a line.
    */
    Line_2 line() const;

    /*!
      returns the supporting line of `c`.
      \pre `c` is not a point (a degenerate curve).
    */
    Line_2 supporting_line() const;

    /*!
      returns `c`'s source point.
      \pre `c` is a point, a segment or a ray.
    */
    Point_2 source() const;

    /*!
      returns `c`'s target point.
      \pre `c` is a point or a segment.
    */
    Point_2 target() const;

    /// @}

  }; /* end Arr_linear_traits_2::Curve_2 */

class Trim_2{
public:
/// \name Creation
/// @{

/*!
Trims the given x-monotone curve to an from src to tgt.
\ pre `src` and `tgt` lies on the curve
*/

X_monotone_curve_2(const X_monotone_curve_2& xcv,
                                const Point_2& src,
                                const Point_2& tgt)const
/// @}
}/* end Arr_linear_traits_2::Trim_2 */

}; /* end Arr_linear_traits_2 */
} /* end namespace CGAL */
