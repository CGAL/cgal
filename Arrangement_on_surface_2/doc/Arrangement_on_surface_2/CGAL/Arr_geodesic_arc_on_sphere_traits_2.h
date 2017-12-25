namespace CGAL {

  /*! \ingroup PkgArrangement2TraitsClasses
   *
   * The traits class `Arr_geodesic_arc_on_sphere_traits_2` is a model of the
   * `ArrangementTraits_2` concept. It enables the construction and
   * maintenance of arrangements of arcs of great circles (also known as
   * geodesic arcs) embedded on the sphere (centered at the origin). Almost
   * all operations on arrangements require a kernel that supports exact
   * predicates. Most operations also require the kernel to support exact
   * constructions. However, all operations on such arrangements can be
   * computed efficiently, since all calculations are performed with
   * rational arithmetic.
   *
   * There is an analogy between this class of arrangements and the class of
   * planar arrangements induced by linear curves (i.e., segments, rays, and
   * lines), as properties of linear curves in the plane often, but not always,
   * hold for geodesic arcs on the sphere. For example, given any two
   * non-antipodal points on the sphere there exists a unique great circle
   * connecting the two points.
   *
   * We use the following parameterization of the unit sphere \f$S =
   * \phi_S(\Phi)\f$: \f$\Phi = [\alpha, 2\pi + \alpha] \times [-\frac{\pi}{2},
   * \frac{\pi}{2}]\f$, \f$\phi_S(x, y) = (\cos y \cos x, \sin y \cos x, \sin
   * x)\f$, where \f$\alpha = atan2(X, Y)\f$. By deafult, \f$X = -1, Y = 0\f$,
   * which implies \f$\alpha = \pi\f$, and a default parameterization \f$\Phi =
   * [-\pi, \pi] \times [-\frac{\pi}{2}, \frac{\pi}{2}]\f$. The equator curve,
   * for example, is given by \f$\gamma(t) = (\pi(2t - 1) + \alpha, 0)\f$, for
   * \f$t \in [0,1]\f$.  This parameterization induces two contraction points
   * \f$p_s = (0, 0, -1) = \phi_S(y,-\frac{\pi}{2})\f$ and \f$p_n = (0, 0, 1) =
   * \phi_S(y,\frac{\pi}{2})\f$, referred to as the south and north poles,
   * respectively, and an identification curve \f$\{\phi_S(\pi +
   * \alpha,x)\,|\,-\frac{\pi}{2} \leq v \leq \frac{\pi}{2}\}\f$, as
   * \f$\phi_S(-\pi + \alpha,v) = \phi_S(+\pi + \alpha,v)\f$ for all \f$x\f$
   * (which coincides with the opposite Prime (Greenwich) Meridian when
   * \f$\alpha = \pi\f$).

   * \cgalModels `ArrangementTraits_2`
   * \cgalModels `ArrangementLandmarkTraits_2`
   * \cgalModels `ArrangementSphericalBoundaryTraits_2`
   */

  template <typename Kernel, typename X, typename Y>
  class Arr_geodesic_arc_on_sphere_traits_2 {
  public:
    /*! The `Point_2` class nested within the traits is used to represent a
     * point on a sphere centered at the origin. The point is in fact a
     * not-necessarily normalized 3D direction extended with information that
     * specifies the location of the point pre-image in the parameter space.
     *
     * \cgalModels `Assignable`
     * \cgalModels `DefaultConstructible`
     * \cgalModels `CopyConstructible`
     */
    class Point_2 {
    public:
      /// \name Enumeration types
      /// @{

      /*! The location type indicates a location in the parameter space.
       */
      enum Location_type {
        /// Internal to the parameter space.
        NO_BOUNDARY_LOC = 0,

        /// The bottom side boundary of the parameter space (the south pole).
        MIN_BOUNDARY_LOC,

        /// The identified left and right side boundaries of the parameter space.
        MID_BOUNDARY_LOC,

        /// The top side boundary of the parameter space (the north pole).
        MAX_BOUNDARY_LOC
      };
      /// @}

      /// \name Types
      /// @{
      typedef Kernel::Direction_3 Direction_3;
      /// @}

      /// \name Creation
      /// @{

      /*! Constructs a point from a direction and a location.
       * \param[in] dir the direction.
       * \param[in] location indicates the location of the point pre-image
       *            in the parameter space.
       */
      Point_2(const Direction_3& dir, Location_type location);

      /// @}

      /// \name Operations
      /// @{

      /*! Set the location of the point pre-image in the parameter space.
       * \param[in] location the updated location of the point pre-image in
       *            the parameter space.
       */
      void set_location(Location_type location);

      /*! Obtain the location of the point.
       * \return the location of the point pre-image in the parameter space.
       */
      Location_type location() const;

      /// @}
    };

    /*! The `X_monotone_curve_2` class nested within the traits is used to
     * represent an \f$x\f$-monotone geodesic arc on the a sphere centered at
     * the origin. The pre-image of an \f$x\f$-monotone geodesic arc does not
     * intersect the identified left and right sides of the boundary of the
     * parameter space.
     *
     * \cgalModels `Assignable`
     * \cgalModels `DefaultConstructible`
     * \cgalModels `CopyConstructible`
     */
    class X_monotone_curve_2 {
    public:
      /// \name Types
      /// @{
      typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>::Point_2 Point_2;
      /// @}

      /// \name Creation
      /// @{

      /*! Constructs an \f$x\f$-monotone geodesic arc.
       * \param[in] source the source point of the arc.
       * \param[in] target the target point of the arc.
       * \param[in] normal the normal of the plane that contains the arc.
       * \param[in] is_vertical is the arc vertical ?
       * \param[in] is_directed_right is the arc directed from left to right?
       * \param[in] is_full is the arc a full great circle?
       * \param[in] is_degenerate is the arc degenerate (single point)?
       * \param[in] is_empty is the arc empty?
       * \pre Both endpoint lie on the given plane.
       */
      X_monotone_curve_2(const Point_2& source,
                         const Point_2& target,
                         const Direction_3& normal,
                         bool is_vertical,
                         bool is_directed_right,
                         bool is_full = false,
                         bool is_degenerate = false,
                         bool is_empty = false);

      /*! Construct an \f$x\f$-monotone geodesic arc.
       * \param[in] normal the normal of the plane containing the arc.
       * \param[in] source the source-point direction.
       * \param[in] target the target-point direction.
       * \pre Both endpoint lie on the given plane.
       */
      X_monotone_curve_2(const Point_2& source,
                         const Point_2& target,
                         const Direction_3& normal);

      /*! Construct a full great-circle.
       * \param[in] point the endpoint of the full great-circle.
       * \param[in] normal the normal of the plane containing the arc.
       * \pre the point lies on the given plane.
       * \pre the point pre-image lies on the identified left and right sides
       *      of the boundary of the parameter space.
       */
      X_monotone_curve_2(const Point_2& point,
                         const Direction_3& normal);

      /// @}

      /// \name Operations
      /// @{

      /*! Sets the source endpoint.
       * \param[in] source the updated source endpoint.
       */
      void set_source(const Point_2& source);

      /*! Sets the target endpoint.
       * \param[in] target the updated target endpoint.
       */
      void set_target(const Point_2& target);

      /*! Sets the normal of the underlying plane.
       * \param[in] normal the updated normal of the underlying plane.
       */
      void set_normal(const Direction_3& normal);

      /*! Sets the flag that indicates whether the arc is vertical.
       * \param[in] flag indicates whether the arc pre-image in the parameter
       *            space is vertical.
       */
      void set_is_vertical(bool flag);

      /*! Sets the flag that indicates whether the direction of the arc
       * pre-image in the parameter space is from left to right.
       * \param flag indicates whether the arc pre-image in the parameter
       *             space is from left to right.
       */
      void set_is_directed_right(bool flag);

      /*! Sets the flag that indicates whether the arc is a full great circle.
       * \param[in] flag indicates whether the arc is a full great circle.
       */
      void set_is_full(bool flag);

      /*! Sets the flag that indicates whether the arc degenerates to a point.
       * \param[in] flag indicates whether the arc degenerates to a point.
       */
      void set_is_degenerate(bool flag);

      /*! Sets the flag that indicates whether the arc is empty.
       * \param[in] flag indicates whether the arc is empty.
       */
      void set_is_empty(bool flag);

      /*! Obtains the source point.
       */
      const Point_2& source() const;

      /*! Obtains the target point.
       */
      const Point_2& target() const;

      /*! Obtains the normal to the containing plane.
       */
      const Direction_3& normal() const;

      /*! Obtains the (lexicographically) left endpoint direction.
       */
      const Point_2& left() const;

      /*! Obtains the (lexicographically) right endpoint.
       */
      const Point_2& right() const;

      /*! Determines whether the arc is vertical.
       */
      bool is_vertical() const;

      /*! Determines whether the arc is directed lexicographically from left to
       * right.
       */
      bool is_directed_right() const;

      /*! Determines whether the arc is a great circle.
       */
      bool is_full() const;

      /*! Determines whether the arc is degenerate.
       */
      bool is_degenerate() const;

      /*! Determines whether the arc is empty. */
      bool is_empty() const;

      /*! Determines whether the arc is a meridian.
       */
      bool is_meridian() const;

      /// @}
    };

    /*!
     */
    class Curve_2 {
    public:
      /// \name Types
      /// @{
      typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>::Point_2 Point_2;
      /// @}

      /// \name Creation
      /// @{
      /// @}

      /// \name Operations
      /// @{
      /// @}
    };

    /*!
     */
    class Construct_point_2 {
    public:
      /// \name Types
      /// @{
      /// @}

      /// \name Operations
      /// @{
      /// @}
    };

    /*! Construction functor of \f$x\f$-monotone geodesic arcs.
     *
     * \cgalModels `Assignable`
     * \cgalModels `CopyConstructible`
     * \cgalModels `AdaptableUnaryFunction`
     * \cgalModels `AdaptableBinaryFunction`
     * \cgalModels `AdaptableTernaryFunction`
     */
    class Construct_x_monotone_curve_2 {
    public:
      /// \name Types
      /// @{
      typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>::Point_2 Point_2;
      typedef Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>::X_monotone_curve_2 result_type;
      typedef Kernel::Direction_3 Direction_3;
      typedef Direction_3 argument_type;
      typedef Point_2 first_argument_type;
      typedef Point_2 second_argument_type;
      typedef Direction_3 third_argument_type;
      /// @}

      /// \name Operations
      /// @{

      /*! Construct the minor geodesic arc from two endpoints. The minor arc
       * is the one with the smaller angle among the two geodesic arcs with
       * the given endpoints.
       * 1. Find out whether the arc is x-monotone.
       * 2. If it is x-monotone,
       *    2.1 Find out whether it is vertical, and
       *    2.2 whether the target is larger than the source (directed right).
       * The arc is vertical, iff
       * 1. one of its endpoint direction pierces a pole, or
       * 2. the projections onto the xy-plane coincide.
       * \param[in] p the first endpoint.
       * \param[in] q the second endpoint.
       * \pre p and q must not coincide.
       * \pre p and q cannot be antipodal.
       */
      X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q);

      /*! Construct a full great circle from a normal to a plane.
       * \param normal the normal to the plane containing the great circle.
       * \pre the plane is not vertical.
       */
      X_monotone_curve_2 operator()(const Direction_3& normal);

      /*! Construct a geodesic arc from two endpoints contained
       * in a plane.
       * \param[in] p the first endpoint.
       * \param[in] q the second endpoint.
       * \param[in] normal the normal to the plane containing the arc.
       * \pre Both endpoint lie on the given plane.
       * \pre Both endpoint lie on the given plane.
       */
      X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q,
                                    const Direction_3& normal)
      /// @} /* end of operations */
    };

    class Construct_Curve_2 {
    public:
      /// \name Types
      /// @{
      /// @}

      /// \name Operations
      /// @{
      /// @}
    };

    /*! Returns an instance of `Construct_point_2`.
     */
    Construct_point_2 construct_point_2_object() const;

    /*! Returns an instance of `Construct_x_monotone_curve_2`.
     */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

    /*! Returns an instance of `Construct_curve_2`.
     */
    Construct_curve_2 construct_curve_2_object() const;
  };

} /* end namespace CGAL */
