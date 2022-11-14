namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * The class `Arr_conic_traits_2` is a model of the `ArrangementTraits_2`
 * concept and can be used to construct and maintain arrangements of bounded
 * segments of algebraic curves of degree \f$ 2\f$ at most, also known as
 * <I>conic curves</I>.
 *
 * A general conic curve \f$ C\f$ is the locus of all points \f$ (x,y)\f$
 * satisfying the equation: \f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$,
 * where:
 *
 * <UL>

 * <LI>If \f$ 4 r s - t^2 > 0\f$, \f$ C\f$ is an ellipse.  A special case occurs
 * when \f$ r = s\f$ and \f$ t = 0\f$, when \f$ C\f$ becomes a circle.
 *
 * <LI>If \f$ 4 r s - t^2 < 0\f$, \f$ C\f$ is a hyperbola.
 *
 * <LI>If \f$ 4 r s - t^2 = 0\f$, \f$ C\f$ is a parabola.  A degenerate case
 * occurs when \f$ r = s = t = 0\f$, when \f$ C\f$ is a line.
 *
 * </UL>
 *
 * A <I>bounded conic arc</I> is defined as either of the following:
 *
 * <UL>
 *
 * <LI>A full ellipse (or a circle) \f$ C\f$.
 *
 * <LI>The tuple \f$ \langle C, p_s, p_t, o \rangle\f$, where \f$ C\f$ is the
 * supporting conic curve, with the arc endpoints being \f$ p_s\f$ and \f$
 * p_t\f$ (the source and target points, respectively). The orientation \f$ o\f$
 * indicates whether we proceed from \f$ p_s\f$ to \f$ p_t\f$ in a clockwise or
 * in a counterclockwise direction. Note that \f$ C\f$ may also correspond to a
 * line or to pair of lines - in this case \f$ o\f$ may specify a `COLLINEAR`
 * orientation.
 *
 * </UL>
 *
 * A very useful subset of the set of conic arcs are line segments and circular
 * arcs, as arrangements of circular arcs and line segments have some
 * interesting applications (e.g. offsetting polygons, motion planning for a
 * disc robot, etc.). Circular arcs and line segment are simpler objects and can
 * be dealt with more efficiently than arbitrary arcs. For these reasons, it is
 * possible to construct conic arcs from segments and from circles. Using these
 * constructors is highly recommended: It is more straightforward and also
 * speeds up the arrangement construction. However, in case the set of input
 * curves contain only circular arcs and line segments, it is recommended to use
 * the `Arr_circle_segment_2` class to achieve faster running times.
 *
 * In our representation, all conic coefficients (namely \f$ r, s, t, u, v,
 * w\f$) must be rational numbers. This guarantees that the coordinates of all
 * arrangement vertices (in particular, those representing intersection points)
 * are algebraic numbers of degree \f$ 4\f$ (a real number \f$ \alpha\f$ is an
 * algebraic number of degree \f$ d\f$ if there exist a polynomial \f$ p\f$ with
 * <I>integer</I> coefficient of degree \f$ d\f$ such that \f$ p(\alpha) =
 * 0\f$).  We therefore require separate representations of the curve
 * coefficients and the point coordinates. The `NtTraits` should be instantiated
 * with a class that defines nested `Integer`, `Rational` and `Algebraic` number
 * types and supports various operations on them, yielding certified computation
 * results (for example, it can convert rational numbers to algebraic numbers
 * and can compute roots of polynomials with integer coefficients).  The other
 * template parameters, `RatKernel` and `AlgKernel` should be geometric kernels
 * templated with the `NtTraits::Rational` and `NtTraits::Algebraic` number
 * types, respectively.  It is recommended to instantiate the
 * `CORE_algebraic_number_traits` class as the `NtTraits` parameter, with
 * `Cartesian<NtTraits::Rational>` and `Cartesian<NtTraits::Algebraic>`
 * instantiating the two kernel types, respectively.  The number types in this
 * case are provided by the \core library, with its ability to exactly represent
 * simple algebraic numbers.
 *
 * The traits class inherits its point type from `AlgKernel::Point_2`,
 * and defines a curve and \f$ x\f$-monotone curve types, as detailed below.
 *
 * While the `Arr_conic_traits_2` models the concept
 * `ArrangementDirectionalXMonotoneTraits_2`, the implementation of
 * the `Are_mergeable_2` operation does not enforce the input curves
 * to have the same direction as a precondition. Moreover, `Arr_conic_traits_2`
 * supports the merging of curves of opposite directions.
 *
 * \cgalModels `ArrangementTraits_2`
 * \cgalModels `ArrangementLandmarkTraits_2`
 * \cgalModels `ArrangementDirectionalXMonotoneTraits_2`
 *
 * \cgalHeading{Types}
 */
template <typename RatKernel, typename AlgKernel, typename NtTraits>
class Arr_conic_traits_2 {
public:

  /// \name Types
  /// @{

  /*! the `NtTraits::Rational` type (and also the `RatKernel::FT` type).
   */
  typedef unspecified_type      Rational;

  /*! the `NtTraits::Algebraic` type (and also the `AlgKernel::FT` type).
   */
  typedef unspecified_type      Algebraic;

  /// @}

  /*! The `Curve_2` class nested within the conic-arc traits can represent
   * arbitrary conic arcs and support their construction in various ways.  The
   * copy and default constructor as well as the assignment operator are
   * provided for conic arcs. In addition, an `operator<<` for the curves is
   * defined for standard output streams.
   */
  class Curve_2 {
  public:

    /// \name Creation
    /// @{

    /*! constructs an arc corresponding to the line segment `seg`.
     */
    Curve_2(const typename RatKernel::Segment_2& seg);

    /*! constructs an arc corresponding to the full circle `circ`
     * (note that this circle has a center point with rational coordinates
     * and rational squared radius).
     */
    Curve_2(const typename RatKernel::Circle_2& circ);

    /*! constructs a circular arc supported by the circle `circ`, going
     * in the given orientation `o` from the source point `ps` to its target
     * point `pt`.
     *
     * \pre `ps` and `pt` both lie on the circle `circ`.
     *
     * \pre `o` is not `COLLINEAR`.
     */
    Curve_2(const typename RatKernel::Circle_2& circ, Orientation o,
            const Point_2& ps, const Point_2& pt);

    /*! constructs a circular arc going from `p1` (its source point)
     * through `p2` to `p3` (its target point). Note that all three points have
     * rational coordinates. The orientation of the arc is determined
     * automatically.
     *
     * \pre The three points are not collinear.
     */
    Curve_2(const typename RatKernel::Point_2& p1,
            const typename RatKernel::Point_2& p2,
            const typename RatKernel::Point_2& p3);

    /*! constructs a conic arc that corresponds to the full conic curve
     * \f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$.
     *
     * \pre As a conic arc must be bounded, the given curve must be an ellipse,
     * that is \f$ 4 r s - t^2 > 0\f$.
     */
    Curve_2(const Rational& r, const Rational& s,
            const Rational& t, const Rational& u,
            const Rational& v, const Rational& w);

    /*! constructs a conic arc supported by the conic curve
     * \f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$, going in the given
     * orientation `o` from the source point `ps` to its target point `pt`.  * *
     *
     * \pre `ps` and `pt` both satisfy the equation of the supporting conic
     * curve and define a bounded segment of this curve (e.g. in case of a
     * hyperbolic arc, both point should be located on the same branch of the
     * hyperbola).
     *
     * \pre `o` is not `COLLINEAR` if the supporting conic is curves, and must
     * be `COLLINEAR` if it is not curved (a line or a line-pair).
     */
    Curve_2(const Rational& r, const Rational& s,
            const Rational& t, const Rational& u,
            const Rational& v, const Rational& w,
            Orientation o,
            const Point_2& ps, const Point_2& pt);

    /*! constructs a conic arc going from `p1` (its source point)
     * through `p2`, `p3` and `p4` (in this order) to `p5` (its target
     * point). Note that all five points have rational coordinates.  The
     * orientation of the arc is determined automatically.
     *
     * \pre No three points of the five are not collinear.
     *
     * \pre The five points define a valid arc, in their given order.
     */
    Curve_2(const typename RatKernel::Point_2& p1,
            const typename RatKernel::Point_2& p2,
            const typename RatKernel::Point_2& p3,
            const typename RatKernel::Point_2& p4,
            const typename RatKernel::Point_2& p5);

    /*! constructs a conic arc supported by the conic curve
     * \f$ r x^2 + s y^2 + t x y + u x + v y + w = 0\f$, going in the given
     * orientation `o` from its source point to its target Point. In this case
     * only some approximations of the endpoints (`app_ps` and `app_pt`,
     * respectively) is available, and their exact locations are given
     * implicitly, specified by the intersections of the supporting conic curve
     * with \f$ r_1 x^2 + s_1 y^2 + t_1 x y + u_1 x + v_1 y + w_1 = 0\f$ and \f$
     * r_2 x^2 + s_2 y^2 + t_2 x y + u_2 x + v_2 y + w_2 = 0\f$, respectively.
     *
     * \pre The two auxiliary curves specifying the endpoints really intersect
     * with the supporting conic curve, such that the arc endpoints define a
     * bounded segment of the supporting curve (e.g. in case of a hyperbolic
     * arc, both point should be located on the same branch of the hyperbola).
     *
     * \pre `o` is not `COLLINEAR` if the supporting conic is curves, and must
     * be `COLLINEAR` if it is not curved (a line or a line-pair).
     */
    Curve_2(const Rational& r, const Rational& s,
            const Rational& t, const Rational& u,
            const Rational& v, const Rational& w,
            Orientation o,
            const Point_2& app_ps,
            const Rational& r1, const Rational& s1,
            const Rational& t1, const Rational& u1,
            const Rational& v1, const Rational& w1,
            const Point_2& app_pt,
            const Rational& r2, const Rational& s2,
            const Rational& t2, const Rational& u2,
            const Rational& v2, const Rational& w2);

    /// @}

    /// \name Access Functions
    /// @{

    /*! indicates whether `a` is a valid conic arc. As the precondition to
     * some of the constructor listed above are quite complicated, their
     * violation does not cause the program to abort. Instead, the constructed
     * arc is invalid (a defaultly constructed arc is also invalid).  It is
     * however recommended to check that a constructed arc is valid before
     * inserting it to an arrangement, as this operation <I>will</I> cause the
     * program to abort.
     */
    bool is_valid() const;

    /*! determines whether the arc is \f$ x\f$-monotone, namely each vertical
     * line intersects it at most once. A vertical line segment is also
     * considered (weakly) \f$ x\f$-monotone.
     */
    bool is_x_monotone() const;

    /*! determines whether the arc is \f$ y\f$-monotone, namely each horizontal
     * line intersects it at most once. A horizontal line segment is also
     * considered (weakly) \f$ x\f$-monotone.
     */
    bool is_y_monotone() const;

    /*! indicates whether the arc represents a full conic curve (en ellipse or
     * a circle).
     */
    bool is_full_conic() const;

    /// @}

    /*! \name
     * The six following methods return the coefficients of the supported conic,
     * after their conversion to integer number (represented by the `Integer`
     * type of the `NtTraits` class):
     */
    /// @{

    /*! returns the coefficient of \f$ x^2\f$.
     */
    const typename NtTraits::Integer& r() const;

    /*! returns the coefficient of \f$ t^2\f$.
     */
    const typename NtTraits::Integer& s() const;

    /*! returns the coefficient of \f$ x y\f$.
     */
    const typename NtTraits::Integer& t() const;

    /*! returns the coefficient of \f$ x\f$.
     */
    const typename NtTraits::Integer& u() const;

    /*! returns the coefficient of \f$ y\f$.
     */
    const typename NtTraits::Integer& v() const;

    /*! returns the free coefficient.
     */
    const typename NtTraits::Integer& w() const;

    /*! returns the source point of the arc.
     * \pre `a` is not a full conic curve.
     */
    const Point_2& source() const;

    /*! returns the target point of the arc.
     * \pre `a` is not a full conic curve.
     */
    const Point_2& target() const;

    /*! returns the orientation of the arc.
     */
    Orientation orientation() const;

    /*! return a bounding box of the arc `a`.
     */
    Bbox_2 bbox() const;

    /// @}

    /// \name Operations
    /// @{

    /*! sets a new source point for the conic arc.
     * \pre `ps` lies on the supporting conic of `a`.
     */
    void set_source(const Point_2 & ps);

    /*! sets a new target point for the conic arc.
     * \pre `pt` lies on the supporting conic of `a`.
     */
    void set_target(const Point_2 & pt);

    /// @}

  }; /* end Arr_conic_traits_2::Curve_2 */

  /*! The `X_monotone_curve_2` class nested within the conic-arc traits is
   * used to represent \f$x\f$-monotone conic arcs. It inherits from the
   * `Curve_2` type, therefore supports the access methods and the operations
   * listed above.
   *
   * For efficiency reasons, we recommend users not to construct \f$
   * x\f$-monotone conic arc directly, but rather use the `Make_x_monotone_2`
   * functor supplied by the conic-arc traits class to convert conic curves to
   * \f$ x\f$-monotone curves.
   */
  class X_monotone_curve_2 {
  public:

    /// \name Creation
    /// @{

    /*! converts the given arc to an \f$ x\f$-monotone arc.
     * \pre `arc` is \f$ x\f$-monotone.
     */
    X_monotone_curve_2(const Curve_2& arc);

    /// @}

    /// \name Access Functions
    /// @{

    /*! returns the left (lexicographically smaller) endpoint of `xa`.
     */
    const Point_2& left() const;

    /*! returns the right (lexicographically larger) endpoint of `xa`.
     */
    const Point_2& right() const;

    /// @}

  }; /* end Arr_conic_traits_2::X_monotone_curve_2 */

  /*! The `Point_2` class nested within the conic-arc traits is
   * used to represent points. It inherits from the algebraic kernel point.
   */
  class Point_2 : public Alg_kernel::Point_2 {
  public:
    /// \name Creation
    /// @{

    /*! construct a default point.
     */
    Point_2();

    /*! construct a point from an algebraic point.
     */
    Point_2(const typename Alg_kernel::Point_2& p);

    /*! constructs from homegeneous coordinates.
     */
    Point_2(const Algebraic& hx, const Algebraic& hy, const Algebraic& hz);

    /*! constructs from Cartesian coordinates.
     */
    Point_2(const Algebraic& x, const Algebraic& y);:

    /// @}
  };

  class Trim_2 {
  public:
    /// \name Creation
    /// @{

    /*! Trims the given x-monotone curve to an from src to tgt.
     * \ pre `src` and `tgt` lies on the curve
     */
    X_monotone_curve_2(const X_monotone_curve_2& xcv,
                       const Point_2& src,
                       const Point_2& tgt) const;

    /// @}
  }/* end Arr_conic_traits_2::Trim_2 */


}; /* end Arr_conic_traits_2 */
} /* end namespace CGAL */
