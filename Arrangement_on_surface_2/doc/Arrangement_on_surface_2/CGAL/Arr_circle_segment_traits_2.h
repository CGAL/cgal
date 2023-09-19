namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * The class `Arr_circle_segment_traits_2` is a model of the
 * `ArrangementTraits_2` concept and can be used to construct and maintain
 * arrangements of circular arcs and line segments.
 *
 * The traits class must be instantiated with a geometric kernel, such that the
 * supporting circles of the circular arcs are of type `Kernel::Circle_2` and
 * the supporting lines of the line segments are of type `Kernel::Line_2`.
 * Thus, the coordinates of the center of supporting circles, and its squared
 * radius are of type `Kernel::FT`, which should be an exact rational
 * number-type; similarly, the coefficients of each supporting line \f$ ax + by
 * + c = 0\f$ are also of type `Kernel::FT`. Note however that the intersection
 * point between two such arcs do not have rational coordinates in general. For
 * this reason, we do not require the endpoints of the input arcs and segments
 * to have rational coordinates.
 *
 * The nested `Point_2` type defined by the traits class is therefore
 * <I>different</I> than the `Kernel::Point_2` type. Its coordinates are of type
 * `CoordNT`, which an instantiation of `Sqrt_extension<NT,ROOT>` where `NT =
 * ROOT = Kernel::FT`.  Moreover, the third and fourth (hidden) template
 * parameters of `Sqrt_extension<NT,ROOT>` are set to `CGAL::Tag_true`, which
 * enables efficient comparison among different extensions.
 *
 * For more details see the documentation of `Sqrt_extension<NT,ROOT>`.
 *
 * While `Arr_circle_segment_traits_2` models the concept
 * `ArrangementDirectionalXMonotoneTraits_2`, the implementation of the
 * `Are_mergeable_2` operation does not enforce the input curves to have the
 * same direction as a precondition. Moreover, `Arr_circle_segment_traits_2`
 * supports the merging of curves of opposite directions.
 *
 * \cgalModels{ArrangementTraits_2,ArrangementDirectionalXMonotoneTraits_2}
 *
 */
template< typename Kernel >
class Arr_circle_segment_traits_2 {
public:

  /*! The `Curve_2` class nested within the traits class can represent
   * arbitrary circular arcs, full circles and line segments and support their
   * construction in various ways.  The copy and default constructor as well as
   * the assignment operator are provided. In addition, an `operator<<` for the
   * curves is defined for standard output streams.
   */
  class Curve_2 {
  public:

    /// \name Creation
    /// @{

    /*! constructs an curve corresponding to the line segment `seg`.
     */
    Curve_2(const typename Kernel::Segment_2& seg);

    /*! constructs an curve corresponding to the line segment directed
     * from `source` to `target`, both having rational coordinates.
     */
    Curve_2(const typename Kernel::Point_2& source,
            const typename Kernel::Point_2& target);

    /*! constructs an curve corresponding to the line segment supported by
     * the given line, directed from `source` to `target`. Note that the two
     * endpoints may have one-root coordinates.
     *
     * \pre Both endpoints must lie on the given supporting line.
     */
    Curve_2(const typename Kernel::Line_2& line,
            const Point_2& source,
            const Point_2& target);

    /*! constructs an curve corresponding to the given circle. `circ`
     * has a center point with rational coordinates and its <I>squared</I>
     * radius is rational.
     */
    Curve_2(const typename Kernel::Circle_2& circ);

    /*! constructs an curve corresponding to a circle centered at the rational
     * point `c` whose radius `r` is rational.
     */
    Curve_2(const typename Kernel::Point_2& c,
            const typename Kernel::FT& r,
            Orientation orient = COUNTERCLOCKWISE);

    /*! constructs a circular arc supported by `circ`, which has a
     * center point with rational coordinates and whose <I>squared</I> radius is
     * rational, with the given endpoints. The orientation of the arc is the
     * same as the orientation of `circ`.

     * \pre Both endpoints must lie on the given supporting circle.
     */
    Curve_2(const typename Kernel::Circle_2& circ,
            const Point_2& source, const Point_2& target);

    /*! constructs a circular arc supported by a circle centered at the rational
     * point `c` whose radius `r` is rational, directed from `source` to
     * `target` with the given orientation.
     *
     * \pre Both endpoints must lie on the supporting circle.
     */
    Curve_2(const typename Kernel::Point_2& c,
            const typename Kernel::FT& r,
            Orientation orient,
            const Point_2& source, const Point_2& target);

    /*! constructs an circular arc whose endpoints are `source` and
     * `target` that passes through `mid`. All three points have
     * rational coordinates.
     *
     * \pre The three points must not be collinear.
     */
    Curve_2 (const typename Kernel::Point_2& source,
             const typename Kernel::Point_2& mid,
             const typename Kernel::Point_2& target);

    /// @}

    /// \name Access Functions
    /// @{

    /*! indicates whether the curve represents a full circle.
     */
    bool is_full() const;

    /*! returns the source point.
     *
     * \pre `cv` is not a full circle.
     */
    const Point_2& source() const;

    /*! returns the target point.
     *
     * \pre `cv` is not a full circle.
    */
    const Point_2& target() const;

    /*! returns the orientation of the curve (`COLLINEAR` in case of line
     * segments).
     */
    Orientation orientation() const;

    /*! determines whether `cv` is a line segment.
     */
    bool is_linear() const;

    /*! determines whether `cv` is a circular arc.
     */
    bool is_circular() const;

    /*! returns the supporting line of `cv`.
     *
     * \pre `cv` is a line segment.
     */
    typename Kernel::Line_2 supporting_line() const;

    /*! returns the supporting circle of `cv`.
     *
     * \pre `cv` is a circular arc.
     */
    typename Kernel::Circle_2 supporting_circle() const;

    /// @}

  }; /* end Arr_circle_segment_traits_2::Curve_2 */


  /*! The `Point_2` number-type nested within the traits class represents
   * a %Cartesian point whose coordinates are algebraic numbers of type
   * `CoordNT`.
   */
  class Point_2 {
  public:

    /// \name Types
    /// @{

    /*! the `Kernel::FT` type.
     */
    typedef unspecified_type Rational;

    /*! the algebraic number-type.
     */
    typedef unspecified_type CoordNT;

    /// @}

    /// \name Creation
    /// @{

    /*! default constructor.
     */
    Point_2();

    /*! creates the point \f$ (x,y)\f$.
     */
    Point_2(const Rational& x, const Rational& y);

    /*! creates the point \f$ (x,y)\f$.
     */
    Point_2(const CoordNT& x, const CoordNT& y);

    /// @}

    /// \name Access Functions
    /// @{

    /*! returns the \f$ x\f$-coordinate.
     */
    CoordNT x() const;

    /*! returns the \f$ y\f$-coordinate.
     */
    CoordNT y() const;

    /// @}

  }; /* end Arr_circle_segment_traits_2::Point_2 */


  /*! The `X_monotone_curve_2` class nested within the traits class can
   * represent \f$ x\f$-monotone and line segments (which are always weakly
   * \f$x\f$-monotone).  The copy and default constructor as well as the
   * assignment operator are provided. In addition, an `operator<<` for the
   * curves is defined for standard output streams.
   */
  class X_monotone_curve_2 {
  public:

    /// \name Creation
    /// @{

    /*! constructs an curve corresponding to the line segment directed
     * from `source` to `target`, both having rational coordinates.
     */
    X_monotone_curve_2 (const typename Kernel::Point_2& source,
                        const typename Kernel::Point_2& target);

    /*! constructs an curve corresponding to the line segment supported by
     * the given line, directed from `source` to `target`.  Note that the two
     * endpoints may have one-root coordinates.
     *
     * \pre Both endpoints must lie on the given supporting line.
     */
    X_monotone_curve_2(const typename Kernel::Line_2& line,
                       const Point_2& source,
                       const Point_2& target);

    /*! constructs a circular arc supported by `circ`, which has a
     * center point with rational coordinates and whose <I>squared</I>
     * radius is rational, with the given endpoints. The orientation of the
     * arc is determined by `orient`.
     *
     * \pre Both endpoints must lie on the given supporting circle.
     *
     * \pre The circular arc is \f$ x\f$-monotone.
     */
    X_monotone_curve_2(const typename Kernel::Circle_2& circ,
                       const Point_2& source, const Point_2& target,
                       Orientation orient);

    /// @}

    /// \name Access Functions
    /// @{

    /*! returns the source point of `xcv`.
     */
    const Point_2& source() const;

    /*! returns the target point of `xcv`.
     */
    const Point_2& target() const;

    /*! returns true if `xcv` is directed right, false otherwise.
     */
    bool is_directed_right () const;

    /*! returns the left (lexicographically smaller) endpoint of `xcv`.
     */
    const Point_2& left() const;

    /*! returns the right (lexicographically larger) endpoint of `xcv`.
     */
    const Point_2& right() const;

    /*! returns the orientation of the curve (`COLLINEAR` in case of line
     * segments).
     */
    Orientation orientation() const;

    /*! determines whether `xcv` is a line segment.
     */
    bool is_linear () const;

    /*! determines whether `xcv` is a circular arc.
     */
    bool is_circular () const;

    /*! returns the supporting line of `xcv`.
     *
     * \pre `xcv` is a line segment.
     */
    typename Kernel::Line_2 supporting_line() const;

    /*! returns the supporting circle of `xcv`.
     *
     * \pre `xcv` is a circular arc.
     */
    typename Kernel::Circle_2 supporting_circle() const;

    /*! returns a bounding box of the arc `xcv`.
     */
    Bbox_2 bbox() const;

    /// @}

  }; /* end Arr_circle_segment_traits_2::X_monotone_curve_2 */

  class Trim_2 {
  public:
    /// \name Creation
    /// @{

    /*! Trims the given x-monotone curve to an from src to tgt.
     * \ pre `src` and `tgt` lies on the curve
     */

    X_monotone_curve_2(const X_monotone_curve_2& xcv,
                       const Point_2& src,
                       const Point_2& tgt)const
    /// @}
      } /* end Arr_circle_segment_traits_2::Trim_2 */

}; /* end Arr_circle_segment_traits_2 */

} /* end namespace CGAL */
