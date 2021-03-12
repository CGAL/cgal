namespace CGAL {
  /*! \ingroup PkgArrangementOnSurface2TraitsClasses
   *
   * Note: The `SubcurveTraits_2` can comprise of Line_segments, Conic_arcs,
   *       Circular_arc, Bezier_curves, or Linear_curves. A portion or a part
   *       of any of the above mentioned geometric traits is called a subcurve.
   *
   * The traits class `Arr_polycurve_traits_2` handles piecewise curves that are
   * not necessarily linear, such as conic arcs, circular arcs, Bezier curves,
   * or line segments. We call such a compound curve a polycurve.  A polycurve
   * is a chain of subcurves, where each two neighboring subcurves in the chain
   * share a common endpoint; that is, the polycurve is continuous. Furthermore,
   * the target of the \f$i\f$th segement of a polycurve has to coincide with
   * the source of the \f$i+1\f$st segment; that is, the polycurve has to be
   * \a well-oriented. Note that it is possible to construct general polycurves
   * that are neither continuous nor well-oriented, as it is impossible to
   * enforce this precondition (using the set of predicates required by the
   * relevant concepts, see below). However, such polycurves cannot be used for
   * the actual computation of arrangements. The traits class template exploits
   * the functionality of the `SubcurveTraits_2` template-parameter to handle
   * the subcurves that compose the polycurve.
   *
   * The type substituting the template parameter `SubcurveTraits_2` when
   * the template Arr_polycurve_traits_2 is instantiated must be a model
   * of the concepts
   *   - `ArrangementTraits_2` and
   *   - `ArrangementDirectionalXMonotoneTraits_2`.
   *
   * If, in addition, the SubcurveTraits_2 models the concept
   * `ArrangementApproximateTraits_2` then `Arr_polycurve_traits_2` models this
   * concept as well. The same holds for the concept
   * `ArrangementOpenBoundaryTraits_2`. If no type is provided, then
   * `Arr_segment_traits_2` (instantiated with
   * `Exact_predicates_exact_constructions_kernel` as the kernel) is used.
   * Otherwise,
   * `Arr_algebraic_segment_traits_2<Coefficient>`,
   * `Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>`,
   * `Arr_circle_segment_traits_2<Kernel>`,
   * `Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>`,
   * `Arr_linear_traits_2<Kernel>`,
   * `Arr_non_caching_segment_traits_2<Kernel>`,
   * `Arr_segment_traits_2<Kernel>`,
   * `Arr_rational_function_traits_2<AlgebraicKernel_d_1>`,
   * or any other model of the concepts above can be used.
   *
   * The number type used by the injected subcurve traits should support exact
   * rational arithmetic (that is, the number type should support the arithmetic
   * operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and \f$ \div\f$ carried out
   * without loss of precision), in order to avoid robustness problems, although
   * other inexact number types could be used at the user's own risk.
   *
   * A polycurve that comprises \f$n > 0\f$ subcurves has \f$ n+1 \f$ subcurve
   * end-points, and they are represented as objects of type
   * `SubcurveTraits_2::Point_2`. Since the notion of a \a vertex is reserved to
   * 0-dimensional elements of an arrangement, we use, in this context, the
   * notion of \a points in order to refer to the vertices of a polycurve. For
   * example, an arrangement induced by a single non-self intersecting polycurve
   * has exactly two vertices regardless of the number of subcurve
   * end-points. Finally, the types `Subcurve_2` and `X_monotone_subcurve_2`
   * nested in `Arr_polycurve_traits_2` are nothing but
   * `SubcurveTraits_2::Curve_2` and `SubcurveTraits_2::X_monotone_curve_2`,
   * respectively.
   *
   * \cgalHeading{A note on Backwards compatibility} In \cgal version 4.2 (and
   * earlier) any object of the `X_monotone_curve_2` type nested in
   * `Arr_polycurve_traits_2` which in that version was called
   * `Arr_polyline_tratis_2` maintained a direction invariant; namely, its
   * vertices were ordered in an \a ascending lexicographical \f$(xy)\f$-order.
   * This restriction is no longer imposed and `X_monotone_curve_2` can be now
   * directed either from right-to-left \a or left-to-right. If you wish to
   * maintain a left-to-right orientations of the \f$x\f$-monotone polycurve,
   * set the macro `CGAL_ALWAYS_LEFT_TO_RIGHT` to 1 before any \cgal header is
   * included.
   *
   * \cgalModels `ArrangementTraits_2`
   * \cgalModels `ArrangementDirectionalXMonotoneTraits_2`
   * \cgalModels `ArrangementApproximateTraits_2` (if the type that substitutes
   *   the template parameter `SubcurveTraits_2` models the concept as well)
   *
   * \sa `Arr_algebraic_segment_traits_2<Coefficient>`
   * \sa `Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>`
   * \sa `Arr_circle_segment_traits_2<Kernel>`
   * \sa `Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>`
   * \sa `Arr_linear_traits_2<Kernel>`
   * \sa `Arr_non_caching_segment_traits_2<Kernel>`
   * \sa `Arr_segment_traits_2<Kernel>`
   * \sa `Arr_rational_function_traits_2<AlgebraicKernel_d_1>`
   * \sa `CGAL_ALWAYS_LEFT_TO_RIGHT`
   */

  template <typename SubcurveTraits_2>
  class Arr_polycurve_traits_2 {
  public:

    /// \name Types
    /// @{
    /*!
     */
    // TODO: Have to turn these into links, so whenever I mention Point_2 it
    //       will point here and *not* to Kernel::Point_2 for instance.
    typedef SubcurveTraits_2::Point_2                   Point_2;

    /*!
     */
    typedef SubcurveTraits_2::Curve_2                   Subcurve_2;
    typedef SubcurveTraits_2::X_monotone_curve_2        X_monotone_subcurve_2;
    /// @}

    /*! Construction functor of a general (not necessarily \f$x\f$-monotone)
     * polycurve.
     *
     * This functor constructs general polycurve. Its `operator()` is
     * oveloaded to support various input types.
     *
     * Note that the composing subcurves, depending on the `SubcurveTraits_2`,
     * might not be \f$x\f$-monotone.
     */
    class Construct_curve_2 {
    public:
      /// \name Operations
      /// @{

      /*! Obtain a polycurve that comprises of one given subcurve.
       * \param subcurve input subcurve.
       * \pre `subcurve` is not degenerated (not tested).
       * \return A polycurve with one subcurve, namely `subcurve`.
       */
      Curve_2 operator()(const Subcurve_2& subcurve) const;

      /*! Construct a well-oriented polycurve from a range of either
       * `SubcurveTraits_2::Point_2` or `SubcurveTraits_2::Curve_2`.
       *
       * \param begin iterator pointing to the first element in the
       *        range.
       * \param end iterator pointing to the past-the-end
       *        element in the range.
       * \pre The given range form a continuous and well-oriented polycurve
       *      (not tested).
       * \pre Contains no degenerated subcurves (not tested)
       * \return A polycurve using the corresponding construction implementation.
       */
      template <typename ForwardIterator>
      Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const;

      /// @} /* end of operations */
    }; /* end of Arr_polycurve_traits_2::Construct_curve_2 */

    /*! Construction functor of \f$x\f$-monotone polycurve.
     *
     * Similar to `Construct_curve_2`, only returns \f$x\f$-monotone
     * polycurve.  Thus, have the same overloads of the
     * `operator()`. Note that when constructing `X_monotone_curve_2`
     * all preconditions are tested.
     *
     * If `CGAL_ALWAYS_LEFT_TO_RIGHT` is defined, then the resulting
     * \f$x\f$-monotone polycurve will be oriented from left-to-right.
     */
    class Construct_x_monotone_curve_2 {};

    /*! Function object which returns the number of subcurve end-points of a
     * polycurve.
     */
    class Number_of_points_2 {};

    /*! Functor to augment a polycurve by either adding a vertex or a subcurve
     * at the back.
     */
    class Push_back_2 {
    public:
      /// \name Operations
      /// @{

      /*! Append a subcurve `subcurve` to an existing polycurve `cv` at the back.
       * If `cv` is empty, `subcurve` will be its first subcurve.
       * \param cv a polycurve. Note, `cv` is (not necessarily) \f$x\f$-monotone.
       * \param subcurve a subcurve (not necessarily \f$x\f$-monotone) to be
       *        appended to `cv`
       */
      void operator()(Curve_2& cv, const Subcurve_2& subcurve) const;

      /*! Append a subcurve `subcurve` to an existing \f$x\f$-monotone polycurve
       * `xcv` at the back. If `xcv` is empty, `subcurve` will be its first
       * subcurve.
       * \param xcv existing \f$x\f$-monotone polycurve
       * \param subcurve the subcurve to be added
       * \pre If `xcv` is not empty then `subcurve` extends `xcv` to the right
       *      if `xcv` is oriented right-to-left. Otherwise, `subcurve` extends
       *      `xcv` to the left.
       * \pre `subcurve` is not degenerated.
       * \pre `xcv` and `subcurve` should have the same orientation
       */
      void operator()(X_monotone_curve_2& xcv,
                      const X_monotone_subcurve_2& subcurve) const;
      /// @} /* end of operations */
    }; /* end of Arr_polycurve_traits_2::Push_back_2 */

    /*! Functor to augment a polycurve by either adding a vertex or a subcurve
     * at the front.
     */
    class Push_front_2 {
    public:
      /// \name Operations
      /// @{

      /*! Append a subcurve `subcurve` to an existing polycurve `cv` at the
       * front. If `cv` is empty, `subcurve` will be its first subcurve.
       * \param cv a polycurve. Note, `cv` is (not necessarily) \f$x\f$-monotone.
       * \param subcurve a subcurve (not necessarily \f$x\f$-monotone) to be
       *        appended to `cv`
       */
      void operator()(Curve_2& cv, const Subcurve_2& subcurve) const;

      /*! Append a subcurve `subcurve` to an existing \f$x\f$-monotone polycurve
       * `xcv` at the front. If `xcv` is empty, `subcurve` will be its first
       * subcurve.
       * \param xcv existing \f$x\f$-monotone polycurve
       * \param subcurve the subcurve to be added
       * \pre If `xcv` is not empty then `subcurve` extends `xcv` to the left if
       *      `xcv` is oriented right-to-left. Otherwise, `subcurve` extends
       *      `xcv` to the right.
       * \pre `subcurve` is not degenerated.
       * \pre `xcv` and `subcurve` should have the same orientation
       */
      void operator()(X_monotone_curve_2& xcv,
                      const X_monotone_subcurve_2& subcurve) const;

      /// @} /* end of operations */
    }; /* end of Arr_polycurve_traits_2::Push_front_2 */

    class Trim_2 {
    public:
      /*! Obtain a trimmed version of the polycurve with src and tgt as end
        * vertices.
        * Src and tgt will be swaped if they do not conform to the direction of
        * the polycurve.
        */
      X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv,
                                    const Point_2& src,
                                    const Point_2& tgt) const;
    };

    /*! Subdivide a given subcurve into x-monotone subcurves and insert them
     * into a given output iterator.
     */
    class Make_x_monotone_2 {
    public:
      /*!
      * \pre if `cv` is not empty then it must be continuous and well-oriented.
      * \param cv the subcurve.
      * \param oi an output iterator for the result. Its value type is a variant
      *           that wraps Point_2 or an X_monotone_curve_2 objects.
      * \return the past-the-end iterator.
      */
      template <typename OutputIterator>
      OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const;
    };

    /*! The `Curve_2` type nested in the `Arr_polycurve_traits_2` represents
     * general continuous piecewise-linear subcurves (a polycurve can be
     * self-intersecting) and support their construction from range of
     * subcurves. Construction of polycurves in various ways is supported using
     * the construction functors. <em>It is strongly recommended to avoid
     * construction of `Curve_2` objects directly and prefer the usage of the
     * construction functors.</em> The type `Curve_2` has two template
     * parameters, namely `SubcurveType_2` and `PointType_2`, which are
     * `SubcurveTraits_2::Curve_2` and `SubcurveTraits_2::Point_2` types. Thus,
     * in general, the subcurves that a `Curve_2` instance comprises could be
     * either \f$x\f$-monotone or not!
     *
     * The copy and default constructor as well as the assignment operator are
     * provided for polycurve subcurves. In addition, an \link
     * PkgArrangementOnSurface2op_left_shift `operator<<` \endlink for the subcurves is
     * defined for standard output streams, and an \link
     * PkgArrangementOnSurface2op_right_shift `operator>>` \endlink for the subcurves is
     * defined for standard input streams.
     */
    template <typename SubcurveType_2, typename PointType_2>
    class Curve_2 {
    public:

    protected:
      /// \name Types
      /// @{

      /*! The container of the subcurves that comprises the polycurve.
       */
      typedef typename std::vector<Subcurve_2>        Subcurves_container;

    public:
      /*! The size of the container that comprises the polycurve.
       */
      typedef typename Subcurves_container::size_type Size;
      typedef typename Subcurves_container::size_type size_type;

      /*! \deprecated
       * A bidirectional iterator that allows traversing the points
       * that comprise a polycurve's subcurves.
       */
      typedef unspecified_type const_iterator;

      /*! \deprecated
       * A bidirectional iterator that allows traversing the points
       * that comprise a polycurve's subcurves.
       */
      typedef unspecified_type const_reverse_iterator;

      /*! A bidirectional constant iterator that allows traversing
       * the subcurves that comprise the polycurve.
       */
      typedef unspecified_type Subcurve_const_iterator;

      /*! A bidirectional constant iterator that allows traversing
       * the subcurves that comprise the polycurve.
       */
      typedef unspecified_type Subcurve_const_reverse_iterator;

      /// @} /* End of Types */

      /// \name Creation
      /// @{

      /*! Default constructor that constructs an empty polycurve.
       */
      Curve_2();

      /*! Construct a polycurve from one subcurve.
       */
      Curve_2(const Subcurve_2 subcurve);

      /*! Construct a polycurve defined by the given range of subcurves
       * `[first, last)` (the value-type of `InputIterator` must be
       * `SubcurveTraits_2::Curve_2`. In general, the subcurves might not
       * be \f$x\f$-monotone, furthermore, they might not form a
       * continuous polycurve.
       *
       * \pre The subcurves in the range should form a continuous and
       *       well-oriented polycurve.
       *
       * \deprecated For backwards compatibility, it is
       * possible to call this constructor with a range whose
       * value-type is `SubcurveTraits_2::Point_2`. In this case, the
       * constructed polycurve will concatenate the \f$n\f$th point
       * with the \f$(n+1)\f$-st point in the range (using a
       * `SubcurveTraits_2::Subcurve_2`'s). This functionality is \a deprecated.
       * Instead use the `Construct_curve_2` functors.
       */
      template <typename InputIterator>
      Curve_2(Iterator first, Iterator last);

      /// @} /* End of Creation */

      /// \name Access Functions
      /// @{

      /*! \deprecated
       * Obtain the number of subcurve end-points that comprise the polycurve.
       * Note that for a bounded polycurve, if there are \f$ n\f$ points in the
       * polycurve, it is comprised of \f$ (n - 1)\f$ subcurves.
       * Currently, only bounded polycurves are supported.
       */
      unsigned_int points() const;

      /*! \deprecated
       * Obtain an iterator pointing at the source point of the polycurve.
       */
      const_iterator begin() const;

      /*! Obtain an iterator pointing at the first subcurve of the polycurve.
       */
      Subcurve_const_iterator begin_subcurves() const;

      /*! \deprecated
       * Obtain an iterator pointing after the end of the polycurve.
       */
      const_iterator end() const;

      /*! Get an iterator pointing at the past-the-end subcurve of the polycurve.
       */
      Subcurve_const_iterator end_subcurves() const;

      /*! \deprecated
       * Obtain an iterator pointing at the target point of the polycurve.
       */
      const_iterator rbegin() const;

      /*! Obtain an iterator pointing at the last subcurve of the polycurve.
       */
      Subcurve_const_reverse_iterator rbegin_subcurves() const;

      /*! \deprecated
       * Obtain an iterator pointing before the beginning of the polycurve.
       */
      const_iterator rend() const;

      /*! Obtain an iterator pointing at the past-the-end subcurve of
       * the polycurve in reverse order.
       */
      Subcurve_const_reverse_iterator rend_subcurves() const;

      /*! \deprecated
       * Obtain the number of subcurves composing the polycurve
       * (equivalent to `pi.points() - 1`). Was replaced by number_of_subcurves()
       */
      size_type size() const;

      /*! Obtain the number of subcurves that comprise the polycurve.
       */
      size_type number_of_subcurves() const;

      /*! Obtain the \f$ k\f$th subcurve of the polycurve.
       * \pre \f$k\f$ is not greater than or equal to \f$n-1\f$, where
       *      \f$n\f$ is the number of subcurves.
       */
      typename SubcurveTraits_2::X_monotone_curve_2
      operator[](size_t k) const;

      /*! Obtain the bounding box of the polycurve.
       */
      Bbox_2 bbox() const;

      /// @} /* End of Access functions */

      /// \name Operations
      /// @{

      /*! Append a subcurve to the polycurve at the back.
       * \a Warning: This function does not preform the precondition test
       *             that the `Push_back_2` functor does. Thus, it is
       *             recommended to use the latter.
       * \param subcurve The new subcurve to be appended to the polycurve.
       * \pre If the polycurve is not empty, the source of `subcurve` must
       *      coincide with the target point of the last subcurve in the
       *      polycurve.
       */
      inline void push_back(const Subcurve_2& subcurve);

      /*! Append a subcurve to the polycurve at the front.
       * \a Warning: This is a risky function! Don't use it! Prefer the
       *             corresponding functor which is provided in the traits
       *             class.
       * \param subcurve The new subcurve to be appended to the polycurve.
       * \pre If the polycurve is not empty, the target of `subcurve` must
       *      coincide with the source point of the first subcurve in the
       *      polycurve.
       */
      inline void push_front(const Subcurve_2& subcurve);

      /*! \deprecated
       * Add a new point to the polycurvs, which becomes the new target point
       * of `pi`.
       * \pre SubcurveTraits_2 is a model of
       * ArrangementConstructXMonotoneCurveTraits_2.
       */
      void push_back(const Point_2 & p);

      /*! Reset the polycurve.
       */
      void clear();

      /// @} /* End of Operations */

    }; /* end Arr_polycurve_traits_2::Curve_2 */


    /*! The `X_monotone_curve_2` class nested within the polycurve
     * traits is used to represent \f$ x\f$-monotone piecewise linear subcurves.
     *
     * It inherits from the `Curve_2` type. `X_monotone_curve_2` can be
     * constructed just like `Curve_2`. However, there is precondition
     * (which is not tested) that the input defines an \f$
     * x\f$-monotone polycurve. Furthermore, in contrast to the general
     * `Curve_2` type, in this case, the subcurves that an
     * `X_monotone_curve_2` comprises have to be instances of the type
     * `SubcurveTraits_2::X_monotone_curve_2`. Note that the \f$
     * x\f$-monotonicity ensures that an \f$ x\f$-monotone polycurve
     * is not self-intersecting. (A self-intersecting polycurve is
     * subdivided into several interior-disjoint \f$x\f$-monotone subcurves).
     *
     * The defined \f$ x\f$-monotone polycurve can be directed either from
     * right-to-left (and in turn its vertices are stored in an ascending
     * lexicographical \f$ xy\f$-order) or left-to-right (and in this case the
     * vertices are stored in a descending lexicographical \f$ xy\f$-order).
     */
    template <typename SubcurveType_2, typename PointType_2>
    class X_monotone_curve_2 {

    }; /* end Arr_polycurve_traits_2::X_monotone_curve_2 */

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
    Number_of_points_2 number_of_points_2_object() const;

    /*!
     */
    Push_back_2 push_back_2_object() const;

    /*!
     */
    Push_front_2 push_front_2_object() const;

    /*!
     */
    Make_x_monotone_2 make_x_monotone_2_object() const;

    /// @} /* End Accessing Functor Objects */

  }; /* end Arr_polycurve_traits_2 */
} /* end namespace CGAL */
