/// \ingroup PkgArrangementOnSurface2Macros
/// @{
/*!
 * If the macro is set to one, then \f$x\f$-monotone curves are always
 * directed from left-to-right.
 */
#define CGAL_ALWAYS_LEFT_TO_RIGHT
/// @}

namespace CGAL {
  /*! \ingroup PkgArrangementOnSurface2TraitsClasses
   *
   * The traits class `Arr_polyline_traits_2` handles piecewise linear
   * curves, commonly referred to as polylines. Each polyline is a
   * chain of segments, where each two neighboring segments in the
   * chain share a common endpoint; that is, the polyline is
   * continuous. Furthermore, the target of the \f$i\f$th segement of
   * a polyline has to coincide with the source of the \f$i+1\f$st
   * segment; that is, the polyline has to be \a well-oriented. Note
   * that it is possible to construct general polylines that are
   * neither continuous nor well-oriented, as it is impossible to
   * enforce this precondition (using the set of predicates required by
   * the relevant concepts, see below). However, such polylines cannot
   * be used for the actual computation of arrangements. The traits
   * class template exploits the functionality of the `SegmentTraits_2`
   * template-parameter to handle the segments that compose the
   * polyline curves.
   *
   * The type substituting the template parameter `SegmentTraits_2` when
   * the template Arr_polyline_traits_2 is instantiated must be a model
   * of the concepts
   *   - `ArrangementTraits_2`,
   *   - `ArrangementDirectionalXMonotoneTraits_2`,
   *   - `ArrangementConstructXMonotoneCurveTraits_2`.
   *
   * If, in addition, the GeometryTraits_2 models the concept
   * `ArrangementApproximateTraits_2` then `Arr_polycurve_traits_2` models
   * this concept as well. The same holds for the concept
   * `ArrangementOpenBoundaryTraits_2`. If no type is provided, then
   * `Arr_segment_traits_2` (instantiated with
   * `Exact_predicates_exact_constructions_kernel` as the kernel) is used.
   * Otherwise,
   * `Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>`,
   * `Arr_circle_segment_traits_2<Kernel>`,
   * `Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>`
   * `Arr_linear_traits_2<Kernel>`
   * `Arr_non_caching_segment_traits_2<Kernel>`,
   * `Arr_segment_traits_2<Kernel>`,
   * or any other model of the concepts above can be used.
   *
   * The number type used by the injected segment traits should support
   * exact rational arithmetic (that is, the number type should support
   * the arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and
   * \f$ \div\f$ carried out without loss of precision), in order to
   * avoid robustness problems, although other inexact number types
   * could be used at the user's own risk.
   *
   * A polyline that comprises \f$n > 0\f$ segments has \f$ n+1 \f$
   * points, and they are represented as objects of type
   * `SegmentTraits_2::Point_2`. Since the notion of a \a vertex is
   * reserved to 0-dimensional elements of an arrangement, we use, in
   * this context, the notion of \a points in order to refer to the
   * vertices of a polyline. For example, an arrangement induced by a
   * single non-self intersecting polyline has exactly two vertices
   * regardless of the number of points. Finally, the types `Segment_2`
   * and `X_monotone_segment_2` nested in `Arr_polyline_traits_2` are
   * nothing but `SegmentTraits_2::Curve_2` and
   * `SegmentTraits_2::X_monotone_curve_2`, respectively.
   *
   * \cgalHeading{A note on Backwards compatibility} In \cgal version
   * 4.2 (and earlier) any object of the `X_monotone_curve_2` type
   * nested in `Arr_polyline_traits_2` maintained a direction
   * invariant; namely, its vertices were ordered in an \a ascending
   * lexicographical \f$(xy)\f$-order.  This restriction is no longer
   * imposed and `X_monotone_curve_2` can be now directed either from
   * right-to-left \a or left-to-right. If you wish to maintain a
   * left-to-right orientations of the \f$x\f$-monotone polylines, set
   * the macro `CGAL_ALWAYS_LEFT_TO_RIGHT` to 1 before any \cgal header
   * is included.
   *
   * \cgalModels `ArrangementTraits_2`
   * \cgalModels `ArrangementDirectionalXMonotoneTraits_2`
   * \cgalModels `ArrangementConstructXMonotoneCurveTraits_2`
   * \cgalModels `ArrangementApproximateTraits_2` (if the type that substitutes
   *   the template parameter `SegmentTraits_2` models the concept as well)
   *
   * \sa `Arr_polycurve_traits_2<SubcurveTraits_2>`
   * \sa `Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>`
   * \sa `Arr_circle_segment_traits_2<Kernel>`
   * \sa `Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>`
   * \sa `Arr_linear_traits_2<Kernel>`
   * \sa `Arr_non_caching_segment_traits_2<Kernel>`
   * \sa `Arr_segment_traits_2<Kernel>`
   * \sa `CGAL_ALWAYS_LEFT_TO_RIGHT`
   */

  template <typename SegmentTraits_2>
  class Arr_polyline_traits_2 : public Arr_polycurve_traits_2<SegmentTraits_2>{
  public:

    /// \name Types
    /// @{
    /*!
     */
    typedef SegmentTraits_2                             Segment_traits_2;
    // TODO: Have to turn these into links, so whenever I mention Point_2 it
    //       will point here and *not* to Kernel::Point_2 for instance.
    typedef SegmentTraits_2::Point_2                    Point_2;

    /*!
     */
    typedef SegmentTraits_2::Curve_2                    Segment_2;
    typedef SegmentTraits_2::X_monotone_curve_2         X_monotone_segment_2;
    /// @}

    /*! Construction functor of a general (not necessarily \f$x\f$-monotone)
     * polyline.
     *
     * This functor constructs general polylines. Its `operator()` is
     * oveloaded to support various input types.
     *
     * Note that the composing segments, depending on the `SegmentTraits_2`,
     * might not be \f$x\f$-monotone.
     */
    class Construct_curve_2 {
    public:
      /// \name Operations
      /// @{

      /*! Obtain a polyline connecting the two given endpoints.
       * \param p The first point.
       * \param q The second point.
       * \pre `p` and `q` are distinct.
       * \return A segment connecting `p` and `q`.
       */
      Curve_2 operator()(const Point_2& p, const Point_2& q) const;

      /*! Obtain a polyline that comprises of one given segment.
       * \param seg input segment
       * \pre `seg` is not degenerated (not tested)
       * \return A polyline with one segment, namely `seg`.
       */
      Curve_2 operator()(const Segment_2& seg) const;

      /*! Construct a well-oriented polyline from a range of either
       * `SegmentTraits_2::Point_2` or `SegmentTraits_2::Segment_2`.
       *
       * \param begin iterator pointing to the first element in the range.
       * \param end iterator pointing to the past-the-end element in the range.
       * \pre The given range form a continuous and well-oriented polyline
       *      (not tested).
       * \pre Contains no degenerated segments (not tested)
       * \return A polyline using the corresponding construction implementation.
       */
      template <typename ForwardIterator>
      Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const;

      /// @} /* end of operations */
    }; /* end of Arr_polyline_traits_2::Construct_curve_2 */

    /*! Construction functor of \f$x\f$-monotone polyline.
     *
     * Similar to `Construct_curve_2`, only returns \f$x\f$-monotone
     * polylines.  Thus, have the same overloads of the
     * `operator()`. Note that when constructing `X_monotone_curve_2`
     * all preconditions are tested.

     * If `CGAL_ALWAYS_LEFT_TO_RIGHT` is defined, then the resulting
     * \f$x\f$-monotone polyline will be oriented from left-to-right.
     */
    class Construct_x_monotone_curve_2 {};

    /*! Functor to augment a polyline by either adding a vertex or a segment
     * at the back.
     */
    class Push_back_2 {
    public:
      /// \name Operations
      /// @{

      /*! Append a point `p` to an existing polyline `cv` at the back.
       * \param cv a polyline. Note, `cv` is not (necessarily)
       *        \f$ x\f$-monotone.
       * \param p a point to be appended to `cv` at the back.
       * \pre `cv` contains at least one segment.
       */
      void operator()(Curve_2& cv, const Point_2& p) const;

      /*! Append a segment `seg` to an existing polyline `cv` at the back.
       * If `cv` is empty, `seg` will be its first segment.
       * \param cv a polyline. Note, `cv` is (not necessarily) \f$x\f$-monotone.
       * \param seg a segment (not necessarily \f$x\f$-monotone) to be appended
       *        to `cv`
       */
      void operator()(Curve_2& cv, const Segment_2& seg) const;

      /*! Append a point `p` to an existing \f$x\f$-monotone polyline `xcv` at
       * the back.
       * \param xcv the existing \f$x\f$-monotone polyline
       * \param p the point to be pushed back.
       * \pre `xcv` contains at least one segment
       * \pre `p` is either to the right of `xcv` if it is oriented
       *      left-to-right or it is to its left if `xcv` is oriented
       *      right-to-left.
       */
      void operator()(const X_monotone_curve_2& xcv, Point_2& p) const;

      /*! Append a segment `seg` to an existing \f$x\f$-monotone polyline `xcv`
       * at the back. If `xcv` is empty, `seg` will be its first segment.
       * \param xcv existing \f$x\f$-monotone polyline
       * \param seg the segment to be added
       * \pre If `xcv` is not empty then `seg` extends `xcv` to the right if
       *      `xcv` is oriented right-to-left. Otherwise, `seg` extends `xcv` to
       *      the left.
       * \pre `seg` is not degenerated.
       * \pre `xcv` and `seg` should have the same orientation
       */
      void operator()(const X_monotone_curve_2& xcv, Segment_2& seg) const;

      /// @} /* end of operations */
    }; /* end of Arr_polyline_traits_2::Push_back_2 */

    /*! Functor to augment a polyline by either adding a vertex or a segment
     * at the front.
     */
    class Push_front_2 {
    public:
      /// \name Operations
      /// @{

      /*! Append a point `p` to an existing polyline `cv` at the front.
       * \param cv a polyline. Note, `cv` is not (necessarily)
       *        \f$ x\f$-monotone.
       * \param p a point to be appended to `cv` at the back.
       * \pre `cv` contains at least one segment.
       */
      void operator()(Curve_2& cv, const Point_2& p) const;

      /*! Append a segment `seg` to an existing polyline `cv` at the front.
       * If `cv` is empty, `seg` will be its first segment.
       * \param cv a polyline. Note, `cv` is (not necessarily) \f$x\f$-monotone.
       * \param seg a segment (not necessarily \f$x\f$-monotone) to be appended
       * to `cv`
       */
      void operator()(Curve_2& cv, const Segment_2& seg) const;

      /*! Append a point `p` to an existing \f$x\f$-monotone polyline `xcv` at
       * the front.
       * \param xcv the existing \f$x\f$-monotone polyline
       * \param p the point to be pushed back.
       * \pre `xcv` contains at least one segment
       * \pre `p` is either to the left of `xcv` if it is oriented
       *      left-to-right or it is to its right if `xcv` is oriented
       *      right-to-left.
       */
      void operator()(const X_monotone_curve_2& xcv, Point_2& p) const;

      /*! Append a segment `seg` to an existing \f$x\f$-monotone polyline `xcv`
       * at the front. If `xcv` is empty, `seg` will be its first segment.
       * \param xcv existing \f$x\f$-monotone polyline
       * \param seg the segment to be added
       * \pre If `xcv` is not empty then `seg` extends `xcv` to the left if
       *      `xcv` is oriented right-to-left. Otherwise, `seg` extends `xcv` to
       *      the right.
       * \pre `seg` is not degenerated.
       * \pre `xcv` and `seg` should have the same orientation
       */
      void operator()(const X_monotone_curve_2& xcv, Segment_2& seg) const;

      /// @} /* end of operations */
    }; /* end of Arr_polyline_traits_2::Push_front_2 */

    /// \name Accessing Functor Objects
    /// @{

    /*!
     */
    Construct_curve_2 construct_curve_2_object() const;

    /*!
     */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

    /*!
     */
    Push_back_2 push_back_2_object() const;

    /*!
     */
    Push_front_2 push_front_2_object() const;

    /// @} /* End Accessing Functor Objects */

  }; /* end Arr_polyline_traits_2 */
} /* end namespace CGAL */
