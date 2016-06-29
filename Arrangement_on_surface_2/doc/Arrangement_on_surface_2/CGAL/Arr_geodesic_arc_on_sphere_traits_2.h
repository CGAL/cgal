namespace CGAL {

  /*!
   * \ingroup PkgArrangement2TraitsClasses
   *
   * The traits class `Arr_geodesic_arc_on_sphere_traits_2` is a model of the
   * `ArrangementTraits_2` concept that enables the construction and
   * maintenance of arrangements of arcs of great circles (also known as
   * geodesic arcs) embedded on the sphere (centered at the origin).
   * Such arrangements can be computed efficiently, since all calculations
   * are performed with (exact) rational arithmetic.
   *
   * There is an analogy between this class of arrangements and the class of
   * planar arrangements induced by linear curves (i.e., segments, rays, and
   * lines), as properties of linear curves in the plane often, but not always,
   * hold for geodesic arcs on the sphere. For example, given any two
   * non-antipodal points on the sphere there exists a unique great circle
   * connecting the two points.
   *
   * We use the following parameterization of the unit sphere: \f$\Phi = [-\pi +
   * \alpha, \pi + \alpha] \times [-\frac{\pi}{2}, \frac{\pi}{2}]\f$,
   * \f$\phi_S(u, v) = (\cos u \cos v, \sin u \cos v, \sin v)\f$, where
   * \f$\alpha\f$ must be substituted with an angle, the arctangent of which is
   * rational, and defaults to 0 when the class is instantiated (at compile
   * time).\cgalFootnote{The actual template parameters of the class are two
   * integers specifying the arctangent of \f$\alpha\f$.}  The equator curve,
   * for example, is given by \f$\gamma(t) = (\pi(2t - 1) + \alpha, 0)\f$, for
   * \f$t \in [0,1]\f$.  This parameterization induces two contraction points
   * \f$p_s = (0, 0, -1) = \phi_S(u,-\frac{\pi}{2})\f$ and \f$p_n = (0, 0, 1)
   * = \phi_S(u,\frac{\pi}{2})\f$, referred to as the south and north poles,
   * respectively, and an identification curve \f$\{\phi_S(\pi +
   * \alpha,v)\,|\,-\frac{\pi}{2} \leq v \leq \frac{\pi}{2}\}\f$, as
   * \f$\phi_S(-\pi + \alpha,v) = \phi_S(+\pi + \alpha,v)\f$ for all \f$v\f$
   * (which coincides with the opposite Prime (Greenwich) Meridian when
   * \f$\alpha = 0\f$).

   * \cgalModels `ArrangementTraits_2`
   * \cgalModels `ArrangementLandmarkTraits_2`
   * \cgalModels `ArrangementSphericalBoundaryTraits_2`
   */

  template <typename Kernel, typename X, typename Y>
  class Arr_geodesic_arc_on_sphere_traits_2 {
  public:
    /*! Construction functor of \f$x\f$-monotone geodesic arcs .
     */
    class Construct_x_monotone_curve_2 {
      /// \name Operations
      /// @{

      /*! Obtain an x-monotone curve connecting the two given endpoints.
       * \param p the first point.
       * \param q the second point.
       * \pre p and q must not coincide.
       * \return a geodesic arc connecting p and q.
       */
      X_monotone_curve_2 operator()(const Point_2& p, const Point_2& q) const;
      /// @} /* end of operations */
    };

    /*!
     */
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;
  };

} /* end namespace CGAL */
