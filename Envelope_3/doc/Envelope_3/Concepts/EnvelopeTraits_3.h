/*! \ingroup PkgEnvelope3Concepts
 * \cgalConcept
 *
 * This concept defines the minimal set of geometric predicates and operations
 * needed to compute the envelope of a set of arbitrary surfaces in \f$
 * \mathbb{R}^3\f$. It refines the `ArrangementXMonotoneTraits_2` concept. In
 * addition to the `Point_2` and `X_monotone_curve_2` types and the
 * `Has_boundary_category` category tag listed in the base concept, it also
 * lists the `Surface_3` and `Xy_monotone_surface_3` types, which represent
 * arbitrary surfaces and \f$ xy\f$-monotone surfaces, respectively, and some
 * constructions and predicates on these types.  Note however, that these
 * operations usually involve the projection of 3D objects onto the \f$
 * xy\f$-plane.
 *
 * \cgalRefines{ArrangementXMonotoneTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Env_triangle_traits_3<Kernel, ArrLinearTraits>}
 * \cgalHasModels{CGAL::Env_sphere_traits_3<ConicTraits>}
 * \cgalHasModels{CGAL::Env_plane_traits_3<Kernel, ArrLinearTraits>}
 * \cgalHasModels{CGAL::Env_surface_data_traits_3<Traits,XyData,SData,Cnv>}
 * \cgalHasModelsEnd
 */

class EnvelopeTraits_3 {
public:

  /// \name Types
  /// @{

  /*! represents an arbitrary surface in \f$ \mathbb{R}^3\f$.
   */
  typedef unspecified_type Surface_3;

  /*! represents a weakly \f$ xy\f$-monotone surface in \f$ \mathbb{R}^3\f$.
   */
  typedef unspecified_type Xy_monotone_surface_3;

  /// @}

  /// \name Functor Types
  /// @{

  /*! provides the operator (templated by the `OutputIterator` type) :
   * <UL>
   * <LI>`OutputIterator operator() (Surface_3 S, bool is_lower, OutputIterator oi)`
   * <BR>
   * which subdivides the given surface `S` into \f$ xy\f$-monotone parts
   * and inserts them into the output iterator. The value of
   * `is_lower` indicates whether we compute the lower or the upper
   * envelope, so that \f$ xy\f$-monotone surfaces that are irrelevant to the
   * lower-envelope (resp. upper-envelope) computation may be discarded.
   * The value-type of `OutputIterator` is `Xy_monotone_surface_3`.
   * The operator returns a past-the-end iterator for the output sequence.
   * </UL>
   */
  typedef unspecified_type Make_xy_monotone_3;

  /*! provides the operator (templated by the `OutputIterator` type) :
   * <UL>
   * <LI>`OutputIterator operator() (Xy_monotone_surface_3 s, OutputIterator oi)`
   * <BR>
   * which computes all planar \f$ x\f$-monotone curves and possibly isolated
   * planar points that form the projection of the boundary of the given
   * \f$ xy\f$-monotone surface \f$ s\f$ onto the \f$ xy\f$-plane, and inserts
   * them into the output iterator.
   * The value-type of `OutputIterator` is `Object`, where `Object` wraps either
   * a `Point_2`, or a `pair<X_monotone_curve_2, Oriented_side>`. In the former
   * case, the object represents an isolated point of the projected boundary. In
   * the latter, more general, case the object represents an \f$ x\f$-monotone
   * boundary curve along with an enumeration value which is either
   * `ON_NEGATIVE_SIDE` or `ON_POSITIVE_SIDE`, indicating whether whether the
   * projection of the surface onto the \f$ xy\f$-plane lies below or above this
   * \f$ x\f$-monotone curve, respectively. In degenerate case, namely when the
   * surface itself is vertical, and its projection onto the plane is
   * \f$ 1\f$-dimensional, the `Oriented_side` value is `ON_ORIENTED_BOUNDARY`.
   * The operator returns a past-the-end iterator for the output sequence.
   * </UL>
   */
  typedef unspecified_type Construct_projected_boundary_2;

  /*! provides the operator (templated by the `OutputIterator` type) :
   * <UL>
   * <LI>`OutputIterator operator() (Xy_monotone_surface_3 s1, Xy_monotone_surface_3 s2, OutputIterator oi)`
   * <BR>
   * which computes the projection of the intersections of the
   * \f$ xy\f$-monotone surfaces `s1` and `s2` onto the \f$ xy\f$-plane,
   * and inserts them into the output iterator.
   * The value-type of `OutputIterator` is `Object`, where
   * each `Object` either wraps a `pair<X_monotone_curve_2,Multiplicity>`
   * instance, which represents a projected intersection curve with its
   * multiplicity (in case the multiplicity is undefined or not known, it
   * should be set to \f$ 0\f$) or an `Point_2` instance, representing the
   * projected image of a degenerate intersection (the projection of an
   * isolated intersection point, or of a vertical intersection curve).
   * The operator returns a past-the-end iterator for the output sequence.
   * </UL>
   */
  typedef unspecified_type Construct_projected_intersections_2;

  /*! provides the operators :
   * <UL>
   * <LI> `Comparison_result operator() (Point_2 p, Xy_monotone_surface_3 s1, Xy_monotone_surface_3 s2)`
   * <BR>
   * which determines the relative \f$ z\f$-order of the two given
   * \f$ xy\f$-monotone surfaces at the \f$ xy\f$-coordinates of the point `p`,
   * with the precondition that both surfaces are defined over `p`. Namely, it
   * returns the comparison result of \f$ s_1(p)\f$ and \f$ s_2(p)\f$.
   * <LI>`Comparison_result operator() (X_monotone_curve_2 c, Xy_monotone_surface_3 s1, Xy_monotone_surface_3 s2)`
   * <BR>
   * which determines the relative \f$ z\f$-order of the two given
   * \f$ xy\f$-monotone surfaces over the interior of a given \f$ x\f$-monotone
   * curve \f$ c\f$, with the precondition that \f$ c\f$ is fully contained in
   * the \f$ xy\f$-definition range of both \f$ s_1\f$ and \f$ s_2\f$, and that
   * the surfaces do not intersect over \f$ c\f$. The functor should therefore
   * return the comparison result of \f$ s_1(p')\f$ and \f$ s_2(p')\f$ for some
   * point \f$ p'\f$ in the interior of \f$ c\f$.
   * <LI>`Comparison_result operator() (Xy_monotone_surface_3 s1, Xy_monotone_surface_3 s2)`
   * <BR>
   * which determines the relative \f$ z\f$-order of the two given unbounded
   * \f$ xy\f$-monotone surfaces, which are defined over the entire
   * \f$ xy\f$-plane and have no boundary, with the precondition that the
   * surfaces do not intersect at all. The functor should therefore return the
   * comparison result of \f$ s_1(p)\f$ and \f$ s_2(p)\f$ for some planar point
   * \f$ p \in\mathbb{R}^2\f$. This operator is required iff the category tag
   * `Has_boundary_category` is defined as `Tag_true`.
   * </UL>
   */
  typedef unspecified_type Compare_z_at_xy_3;

  /*! provides the operator :
   * <UL>
   * <LI>`Comparison_result operator() (X_monotone_curve_2 c, Xy_monotone_surface_3 s1,  Xy_monotone_surface_3 s2)`
   * <BR>
   * which determines the relative \f$ z\f$-order of the two given
   * \f$ xy\f$-monotone surfaces immediately above their projected intersection
   * curve \f$ c\f$ (a planar point \f$ p\f$ is <I>above</I> an
   * \f$ x\f$-monotone curve \f$ c\f$ if it is in the \f$ x\f$-range of
   * \f$ c\f$, and lies to its left when the curve is traversed from its
   * \f$ xy\f$-lexicographically smaller endpoint to its larger endpoint). We
   * have the precondition that both surfaces are defined "above" \f$ c\f$, and
   * their relative \f$ z\f$-order is the same for some small enough
   * neighborhood of points above \f$ c\f$.
   * </UL>
   */
  typedef unspecified_type Compare_z_at_xy_above_3;

  /*! provides the operator :
   * <UL>
   * <LI>`Comparison_result operator() (X_monotone_curve_2 c,  Xy_monotone_surface_3 s1, Xy_monotone_surface_3 s2)`
   * <BR>
   * which determines the relative \f$ z\f$-order of the two given
   * \f$ xy\f$-monotone surfaces immediately below their projected intersection
   * curve \f$ c\f$ (a planar point \f$ p\f$ is <I>below</I> an
   * \f$ x\f$-monotone curve \f$ c\f$ if it is in the \f$ x\f$-range of
   * \f$ c\f$, and lies to its right when the curve is traversed from its
   * \f$ xy\f$-lexicographically smaller endpoint to its larger endpoint). We
   * have the precondition that both surfaces are defined "below" \f$ c\f$, and
   * their relative \f$ z\f$-order is the same for some small enough
   * neighborhood of points below \f$ c\f$.
   * </UL>
   */
  typedef unspecified_type Compare_z_at_xy_below_3;

  /// @}

  /// \name Creation
  /// @{

  /*! default constructor.
   */
  EnvelopeTraits_3();

  /*! copy constructor.
   */
  EnvelopeTraits_3(EnvelopeTraits_3 other);

  /*! assignment operator.
   */
  EnvelopeTraits_3 operator=(other);

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  /*!
   */
  Make_xy_monotone_3 make_xy_monotone_3_object();

  /*!
   */
  Construct_projected_boundary_2 construct_projected_boundary_2_object();

  /*!
   */
  Construct_projected_intersections_2 construct_projected_intersections_2_object();

  /*!
   */
  Compare_z_at_xy_3 compare_z_at_xy_3_object();

  /*!
   */
  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object();

  /*!
   */
  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object();

  /// @}

}; /* end EnvelopeTraits_3 */
