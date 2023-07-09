namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2TraitsClasses
 *
 * The traits class `Arr_segment_traits_2` is a model of the
 * `ArrangementTraits_2` concept, which allows the construction and maintenance
 * of arrangements of line segments. It should be parameterized with a
 * \cgal-kernel model that is templated in turn with a number type. To avoid
 * numerical errors and robustness problems, the number type should support
 * exact rational arithmetic - that is, the number type should support the
 * arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and \f$ \div\f$
 * carried out without loss of precision.
 *
 * For example, instantiating the traits template with kernels that support
 * exact predicates and exact constructions, such as
 * `Exact_predicates_exact_constructions_kernel`, ensures the exact and robust
 * operation of the application. While selecting an inexact kernel when
 * developing a program usually leads to shorter running times, it causes
 * robustness problems in most cases, and thus renders the program useless. The
 * `Exact_predicates_exact_constructions_kernel` type is a shortcut for a kernel
 * that uses an exact number type, that is, `Gmpq` (More precisely,
 * `Lazy_exact_nt<Gmpq>`), when `Gmpq` is available. Indeed, this configuration
 * achieves the shortest running times in many cases.  If no type is provided
 * for the kernel template parameter, then
 * `Exact_predicates_exact_constructions_kernel` will be used by default.
 *
 * `Arr_segment_traits_2` defines `Kernel::Point_2` as its point type. However,
 * it does <I>not</I> define `Kernel::Segment_2` as its curve type, as one may
 * expect. The reason is that the kernel segment is represented by its two
 * endpoints only, while the traits class needs to store extra data with its
 * segments, in order to efficiently operate on them. Nevertheless, the nested
 * `X_monotone_curve_2` and `Curve_2` types (in this case both types refer to
 * the same class, as <I>every</I> line segment is (weakly) \f$ x\f$-monotone)
 * can however be converted to the type `Kernel::Segment_2`.
 *
 * `Arr_segment_traits_2` achieves faster running times than the
 * `Arr_non_caching_segment_traits_2<Kernel>` traits-class, when arrangements
 * with relatively many intersection points are constructed.  It also allows for
 * working with less accurate, yet computationally efficient number types, such
 * as `Quotient<MP_Float>`, which represents floating-point numbers with an
 * unbounded mantissa, but with a bounded exponent. Using this traits class is
 * therefore highly recommended for almost all applications that rely on
 * arrangements of line segments.  On the other hand, `Arr_segment_traits_2`
 * uses more space and stores extra data with each segment, so constructing
 * arrangements of huge sets of non-intersecting segments (or segments that
 * intersect very sparsely) could be more efficient with the
 * `Arr_non_caching_segment_traits_2` traits-class.
 *
 * While `Arr_segment_traits_2` models the concept
 * `ArrangementDirectionalXMonotoneTraits_2`, the implementation of the
 * `Are_mergeable_2` operation does not enforce the input curves to have the
 * same direction as a precondition. Moreover, `Arr_segment_traits_2` supports
 * the merging of curves of opposite directions.
 *
 * \cgalModels{ArrangementTraits_2,ArrangementLandmarkTraits_2,ArrangementDirectionalXMonotoneTraits_2}
 */
template <typename Kernel>
class Arr_segment_traits_2 : public Kernel {
public:

  //! \name Types
  //! @{

  //! the segment type.
  typedef typename Kernel::Segment_2            Segment_2;

  //! the line type.
  typedef typename Kernel::Line_2               Line_2;

  //! the point type.
  typedef typename Kernel::Point_2              Point_2;

  //! @}

  /*! The `X_monotone_curve_2` class nested within the traits class is
   * used to represent segments.
   */
  class X_monotone_curve_2 {
  public:
    //! \name Creation
    //! @{

    /*! constructs default. */
    X_monotone_curve_2();

    //! @}

    //! \name Access Functions
    //! @{

    //! obtains the (lexicographically) left endpoint.
    const Point_2& left() const;

    //! obtains the (lexicographically) right endpoint.
    const Point_2& right() const;

    //! obtains the supporting line.
    const Line_2& line() const;

    //! determines whether the curve is vertical.
    bool is_vertical() const;

    //! determines whether the curve is directed lexicographic from left to right
    bool is_directed_right() const;

    //! @}
  };

  //! The curve type.
  typedef X_monotone_curve_2                    Curve_2;

  //! A functor that trims curves.
  class Trim_2 {
  public:
    //! \name Creation
    //! @{

    /*! trims the given x-monotone curve to an from src to tgt.
     * \ pre `src` and `tgt` lies on the curve
     */
    X_monotone_curve_2(const X_monotone_curve_2& xcv,
                       const Point_2& src, const Point_2& tgt) const;

  } /* end Arr_segment_traits_2::Trim_2 */

}; /* end Arr_segment_traits_2 */
} /* end namespace CGAL */
