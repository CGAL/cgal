namespace CGAL {

  /*!
    \ingroup PkgArrangement2TraitsClasses

    The traits class `Arr_polyline_traits_2` is a model of the
    `ArrangementTraits_2` concept. It handles piecewise linear curves,
    commonly referred to as polylines. Each polyline is a chain of
    segments, where each two neighboring segments in the chain share a
    common endpoint. The traits class exploits the functionality of
    the `SegmentTraits` template-parameter to handle the segments that
    comprise the polyline curves.

    The type substituting the template parameter `SegmentTraits` when
    the template Arr_polyline_traits_2 is instantiated must be a model
    of the `ArrangementTraits_2` and
    `ArrangementDirectionalXMonotoneTraits_2` concepts that handles
    line segments. If no type is provided then
    `Arr_segment_traits_2` (with `Exact_predicates_exact_constructions_kernel`
    as the kernel)
    will be used. Otherwise, `Arr_segment_traits_2<Kernel>` or
    `Arr_non_caching_segment_traits_2<Kernel>` can be used, where the first
    alternative is recommended.

    The number type used by the injected segment traits should support
    exact rational arithmetic (that is, the number type should support
    the arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$ and
    \f$ \div\f$ that should be carried out without loss of precision),
    in order to avoid robustness problems, although other inexact
    number types could be used at the user's own risk.

    `Arr_polyline_traits_2` uses `SegmentTraits::Point_2` as its point
    type. Its nested type `Segment_2` is nothing but `SegmentTraits::Curve_2`.

    \cgalHeading{A note on Backwards compatibility}
    In \cgal version 4.2 (and earlier) the `X_monotone_curve_2` nested type of
    `Arr_polyline_traits_2` maintained a direction invariant, namely, its
    vertices were ordered in an ascending lexicographical \f$(xy)\f$-order.
    This restriction is no longer valid and `X_monotone_curve_2` can be now
    directed either from right-to-left \a or left-to-right. If you wish to
    maintain a left-to-right orientations of the \f$x\f$-monotone polylines set
    the macro `CGAL_ALWAYS_LEFT_TO_RIGHT` to 1 before any \cgal header is
    included.

    \cgalModels `ArrangementTraits_2`
    \cgalModels `ArrangementLandmarkTraits_2`
    \cgalModels `ArrangementDirectionalXMonotoneTraits_2`

    \sa `Arr_segment_traits_2<Kernel>`
    \sa `Arr_non_caching_segment_traits_2<Kernel>`
    \sa `CGAL_ALWAYS_LEFT_TO_RIGHT`
  */

  template< typename SegmentTraits >
  class Arr_polyline_traits_2 {
  public:

    /// \name Types
    /// @{
    /*! */
    // TODO: Have to turn these into links, so whenever I mention Point_2 it
    //       will point here and *not* to Kernel::Point_2 for instance.
    typedef SegmentTraits::Point_2            Point_2;

    /*! */
    typedef SegmentTraits::Curve_2            Segment_2;
    /// @}

    /*!
      Construction functor of a general (not necessarily \f$x\f$-monotone)
      polyline.

      This functor constructs general polylines given various possibilities of
      input, as can be seen by the various overloads of the `operator()`.
     */
    class Construct_curve_2 {
    public:
      /// \name Operations
      /// @{

      /*! Returns an polyline connecting the two given endpoints.
       * \param p The first point.
       * \param q The second point.
       * \pre `p` and `q` are distinct.
       * \return A segment connecting `p` and `q`.
       */
      Curve_2 operator()(const Point_2& p, const Point_2& q) const;

      /*! Returns a polyline consists of one given segment.
       * \param seg input segment
       * \pre `seg` is not degenerated
       * \return A polyline with one segment, namely `seg`.
       */
      Curve_2 operator()(const Segment_2& seg) const;

      /*! Construct a well-oriented polyline from a range of either
       *  the type `SegmentTraits::Point_2` or the type
       *  `SegmentTraits::Segment_2` and \a only one of them.
       *  \param begin iterator pointing to the first element in the range.
       *  \param end iterator pointing to the past-the-end element in the range.
       *  \pre The given range form a continuous and well-oriented polyline.
       *  \return A polyline using the corresponding construction
       *          implementation.
       */
      template <typename ForwardIterator>
      Curve_2 operator()(ForwardIterator begin, ForwardIterator end) const;

      /// @} /* end of operations */
    }; /* end of Arr_polyline_traits_2::Construct_curve_2 */

    /*!
      Construction functor of \f$x\f$-monotone polyline.

      Similar to `Construct_curve_2`, only returns \f$x\f$-monotone polylines.
      Thus, have the same overloads of the `operator()`.
     */
    class Construct_x_monotone_curve_2 {};

    /*!
      Function object which returns the number of vertices of a polyline.
    */
    class Number_of_points_2{};


    /*!
      Functor to augment a polyline by either adding a vertex or a segment.
    */
    class Push_back_2 {
    public:
      /// \name Operations
      /// @{

      /*!
       * Append a point `p` to an existing polyline `cv`.
       * \param cv a polyline. Note, `cv` is not (necessarily)
       *        \f$ x\f$-monotone.
       * \param p a point to be appended to `cv`.
       * \pre `cv` contains at least one segment.
       */
      void operator()(Curve_2& cv, const Point_2& p) const;

      /*!
        Append a segment `seg` to an existing polyline `cv`. If `cv` is
        empty, `seg` will be its first segment.
        \param cv a polyline. Note, `cv` is not (necessarily)
        \f$x\f$-monotone.
        \param seg a segment to be appended to `cv`
        \pre If `cv` is not empty then target of `cv` should coincide
        with the source of `seg` or(!) the source of `cv` should coincide
        with the target of `seg`.
        \post The resulting `cv` is well oriented.
      */
      void operator()(Curve_2& cv, const Segment_2& seg) const;

      /*!
       * Append a point `p` to an existing x-monotone polyline `xcv`.
       * \param xcv the existing \f$x\f$-monotone polyline
       * \param p the point to be pushed back.
       * \pre `xcv` contains at least one segment
       * \pre `p` is either to the right of `xcv` if it is oriented
       *      left-to-right or it is to its left if `xcv` is oriented
       *      right-to-left.
       */
      void operator()(const X_monotone_curve_2& xcv, Point_2& p) const;

      /*!
       * Append a segment `seg` to an existing x-monotone polyline `xcv`. If
       * `xcv` is empty, `seg` will be its first segment.
       * \param xcv existing \f$x\f$-monotone polyline
       * \param seg the segment to be added
       * \pre If `xcv` is not empty then `seg` extends `xcv` to the right if
       *      `xcv` is oriented right-to-left. Otherwise, `seg` extends `xcv` to
       *      the left.
       * \pre `xcv` and `seg` should have the same orientation
       */
      void operator()(const X_monotone_curve_2& xcv, Segment_2& seg) const;

      /// @} /* end of operations */
    }; /* end of Arr_polyline_traits_2::Push_back_2 */

    /*!  The `Curve_2` class nested within the polyline traits is used
      to represent general continuous piecewise-linear curves (a
      polyline can be self-intersecting) and support their
      construction from any range of segments.  Construction of
      polylines from other inputs is supported using the construction
      functors.

      The copy and default constructor as well as the assignment
      operator are provided for polyline curves. In addition, an \link
      PkgArrangement2op_left_shift `operator<<` \endlink for the
      curves is defined for standard output streams, and an \link
      PkgArrangement2op_right_shift `operator>>` \endlink for the
      curves is defined for standard input streams.

    */
    class Curve_2 {
    public:

      /// \name Types
      /// @{

      /*!
        The container of the segments that comprises the polyline.
       */
      typedef typename std::vector<Segment_2>        Segments_container;

      /*!
        The size of the container that comprises the polylines.
       */
      typedef typename Segments_container::size_type Segments_container_size;

      /*!  \deprecated
        A bidirectional iterator that allows traversing the points
        that comprise a polyline curve.
      */
      typedef unspecified_type const_iterator;

      /*! \deprecated
        A bidirectional iterator that allows traversing the points
        that comprise a polyline curve.
      */
      typedef unspecified_type const_reverse_iterator;

      /*!  A bidirectional constant iterator that allows traversing
        the segments the comprise the polyline.
       */
      typedef unspecified_type Segment_const_iterator;

      /*!  A bidirectional constant iterator that allows traversing
        the segments the comprise the polyline.
      */
      typedef unspecified_type Segment_const_reverse_iterator;

      /// @}

      /// \name Creation
      /// @{

      /*!  Default constructor that constructs an empty polyline.
      */
      Curve_2 ();

      /*! Construct a polyline from one segment.
       */
      Curve_2 (const Segment_2 seg);

      /*!  constructs a polyline defined by the given range of points
        `[first, last)` (the value-type of `InputIterator` must be
        `SegmentTraits::Point_2`.  If the range contains \f$ (n +
        1)\f$ points labeled \f$ (p_{0},p_{1},\ldots,p_{n})\f$, the
        generated polyline consists of \f$ n\f$ segments, where the
        \f$ k\f$th segment is defined by the endpoints \f$
        [p_{k-1},p_{k}]\f$. The first point in the range is considered
        as the source point of the polyline while the last point is
        considered as its target.

        \pre There are at least two points in the range.
      */
      template <class InputIterator>
      Curve_2 (Iterator first, Iterator last);

      /// @}

      /// \name Access Functions
      /// @{

      /*!
        \deprecated Returns the number of points that comprise the polyline.
        Note that for a bounded polyline, if there are \f$ n\f$ points in the
        polyline, it is comprised of \f$ (n - 1)\f$ segments.
      */
      size_t points() const;

      /*!
        \deprecated Returns an iterator pointing at the source point of the
        polyline.
      */
      const_iterator begin() const;

      /*!
        \deprecated Returns an iterator pointing after the end of the polyline.
      */
      const_iterator end() const;

      /*!
        \deprecated Returns an iterator pointing at the target point of the
        polyline.
      */
      const_iterator rbegin() const;

      /*!
        \deprecated Returns an iterator pointing before the beginning of the
        polyline.
      */
      const_iterator rend() const;

      /*!
        \deprecated Returns the number of line segments comprising the polyline
        (equivalent to `pi.points() - 1`). Was replaced by number_of_segments()
      */
      size_t size() const;

      /*!
        Returns the number of segments that comprise the polyline.
       */
      Segments_container_size number_of_segments() const;

      /*!
        Returns the \f$ k\f$th segment of the polyline.
        \pre `k` is not greater or equal to `pi.size() - 1`.
      */
      typename SegmentTraits::X_monotone_curve_2
      operator[] (size_t k) const;

      /*!
        Return a bounding box of the polyline.
      */
      Bbox_2 bbox() const;

      /// @}

      /// \name Operations
      /// @{

      /*!
       * Append a segment to the (x-monotone) polyline.
       * Warning: This is a risky function! Don't use it! Prefer the
       *          provided implementation in the traits class.
       * \param seg The new segment to be appended to the polyline.
       * \pre If the polyline is not empty, the source of `seg` must be the
       *      same as the target point of the last segment in the polyline.
       */
      inline void push_back (const Segment_2& seg);

      /*!
        \deprecated adds a new point to the polyline, which becomes the new
        target point of `pi`.
      */
      void push_back (const Point_2 & p);

      /*!
        Resets the polyline.
      */
      void clear();

      /// @}

    }; /* end Arr_polyline_traits_2::Curve_2 */


    /*!  The `X_monotone_curve_2` class nested within the polyline
      traits is used to represent \f$ x\f$-monotone piecewise linear
      curves.

      It inherits from the `Curve_2` type. `X_monotone_curve_2` can be
      constructed just like `Curve_2`. However, there is precondition
      that the input defines an \f$ x\f$-monotone polyline. Note that
      the \f$ x\f$-monotonicity ensures that an \f$ x\f$-monotone
      polyline is never self-intersecting (thus, a self-intersecting
      polyline will be subdivided to several interior-disjoint \f$
      x\f$-monotone subcurves).

      The defined \f$ x\f$-monotone polyline can be directed either from
      right-to-left (and in turn its vertices are stored in an ascending
      lexicographical \f$ xy\f$-order) or left-to-right (and in this case the
      vertices are stored in a descending lexicographical \f$ xy\f$-order).
    */
    class X_monotone_curve_2 {

    }; /* end Arr_polyline_traits_2::X_monotone_curve_2 */

    /// \name Accessing Functor Objects
    /// @{

    /*! */
    Construct_curve_2 construct_curve_2_object() const;

    /*! */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

    /*! */
    Number_of_points_2 number_of_points_2_object() const;

    /*! */
    Push_back_2 push_back_2_object() const;

    /// @}

  }; /* end Arr_polyline_traits_2 */
} /* end namespace CGAL */
