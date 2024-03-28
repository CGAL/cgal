// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein  <wein@post.tau.ac.il>
//            Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_CONIC_X_MONOTONE_ARC_2_H
#define CGAL_CONIC_X_MONOTONE_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Header file for the Conic_x_monotone_arc_2<Conic_arc_2> class.
 */

#include <ostream>

namespace CGAL {

/*! Representation of an x-monotone conic arc.
 * The class is templated by a representation of a general bounded conic arc.
 */
template <typename ConicArc>
class Conic_x_monotone_arc_2 : public ConicArc {
public:
  typedef ConicArc                              Conic_arc_2;
  typedef Conic_x_monotone_arc_2<Conic_arc_2>   Self;

  typedef typename Conic_arc_2::Alg_kernel      Alg_kernel;
  typedef typename Conic_arc_2::Algebraic       Algebraic;

  typedef typename Conic_arc_2::Alg_point_2     Alg_point_2;
  typedef typename Conic_arc_2::Point_2         Point_2;

  // Type definition for the intersection points mapping.
  typedef typename Point_2::Conic_id            Conic_id;

  using Conic_arc_2::sign_of_extra_data;
  using Conic_arc_2::_is_between_endpoints;

  using Conic_arc_2::IS_VALID;
  using Conic_arc_2::IS_FULL_CONIC;

  using Conic_arc_2::flag_mask;
  using Conic_arc_2::reset_flags;
  using Conic_arc_2::set_flag;
  using Conic_arc_2::reset_flag;
  using Conic_arc_2::flip_flag;
  using Conic_arc_2::test_flag;

private:
  typedef Conic_arc_2                             Base;

  typedef typename Conic_arc_2::Integer           Integer;

  Algebraic m_alg_r;    // The coefficients of the supporting conic curve:
  Algebraic m_alg_s;    //
  Algebraic m_alg_t;    //   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0 ,
  Algebraic m_alg_u;    //
  Algebraic m_alg_v;    // converted to algebraic numbers.
  Algebraic m_alg_w;    //

  Conic_id m_id;        // The ID number of the supporting conic curve.

public:
  // Bit masks for the m_info field (the two least significant bits are already
  // used by the base class).
  enum {
    IS_VERTICAL_SEGMENT = Conic_arc_2::LAST_INFO,
    IS_DIRECTED_RIGHT,
    DEGREE_1,
    DEGREE_2,
    PLUS_SQRT_DISC_ROOT,
    FACING_UP,
    FACING_DOWN,
    IS_SPECIAL_SEGMENT,

    DEGREE_MASK = (0x1 << DEGREE_1) | (0x1 << DEGREE_2),
    FACING_MASK = (0x1 << FACING_UP) | (0x1 << FACING_DOWN)
  };

  /// \name Deprecated Constructions.
  //@{

  // /*! Construct an x-monotone arc from a conic arc.
  //  * \param arc The given (base) arc.
  //  * \param id The ID of the base arc.
  //  */
  // CGAL_DEPRECATED Conic_x_monotone_arc_2(const Base& arc, const Conic_id& id) :
  //   Base(arc),
  //   m_id(id)
  // {
  //   CGAL_precondition(arc.is_valid() && id.is_valid());
  //   _set();
  // }

  /*! Construct an x-monotone sub-arc from a conic arc.
   * \param arc The given (base) arc.
   * \param source The source point.
   * \param target The target point.
   * \param id The ID of the base arc.
   */
  CGAL_DEPRECATED Conic_x_monotone_arc_2(const Base& arc,
                                         const Point_2& source,
                                         const Point_2& target,
                                         const Conic_id& id) :
    Base(arc),
    m_id(id)
  {
    CGAL_precondition(arc.is_valid() && id.is_valid());
    set_endpoints(source, target);
    _set();
  }

  /*! Construct a special segment connecting to given endpoints (for the usage
   * of the landmarks point-location strategy).
   * \param source The source point.
   * \param target The target point.
   */
  CGAL_DEPRECATED
  Conic_x_monotone_arc_2(const Point_2& source, const Point_2& target) :
    Base(source, target)
  {
    set_flag(DEGREE_1);

    Alg_kernel ker;
    auto cmp_xy = ker.compare_xy_2_object();
    Comparison_result dir_res = cmp_xy(source, target);
    if (dir_res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Check if the segment is vertical.
    if (CGAL::sign(this->extra_data->b) == ZERO) set_flag(IS_VERTICAL_SEGMENT);

    // Mark that this is a special segment.
    set_flag(IS_SPECIAL_SEGMENT);
  }

  /*! Construct a special segment of a given line connecting to given
   * endpoints.
   * \param a, b, c The coefficients of the supporting line (ax + by + c = 0).
   * \param source The source point.
   * \param target The target point.
   */
  CGAL_DEPRECATED
  Conic_x_monotone_arc_2(const Algebraic& a, const Algebraic& b,
                         const Algebraic& c,
                         const Point_2& source, const Point_2& target) :
    Base()
  {
    // Make sure the two endpoints lie on the supporting line.
    CGAL_precondition(CGAL::sign(a * source.x() + b * source.y() + c) ==
                      CGAL::ZERO);

    CGAL_precondition(CGAL::sign(a * target.x() + b * target.y() + c) ==
                      CGAL::ZERO);

    // Set the basic properties and clear the _info bits.
    this->set_endpoints(source, target);
    this->set_orientation(COLLINEAR);
    this->reset_flags();                // inavlid arc

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result res = ker.compare_x_2_object()(source, target);

    set_flag(Base::IS_VALID);
    set_flag(DEGREE_1);
    if (res == EQUAL) {
      // Mark that the segment is vertical.
      set_flag(IS_VERTICAL_SEGMENT);

      // Compare the endpoints lexicographically.
      res = ker.compare_y_2_object()(source, target);
      CGAL_precondition(res != EQUAL);
      if (res == EQUAL) {
        reset_flags();  // inavlid arc
        return;
      }
    }

    if (res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Store the coefficients of the line.
    this->m_extra_data = new typename Base::Extra_data;
    this->m_extra_data->a = a;
    this->m_extra_data->b = b;
    this->m_extra_data->c = c;
    this->m_extra_data->side = ZERO;

    // Mark that this is a special segment.
    set_flag(IS_SPECIAL_SEGMENT);
  }
  //@}

private:
  /*! Set the properties of the x-monotone conic arc (for the usage of the
   * constructors).
   */
  CGAL_DEPRECATED void _set() {
    // Convert the coefficients of the supporting conic to algebraic numbers.
    typename Base::Nt_traits nt_traits;

    m_alg_r = nt_traits.convert(this->r());
    m_alg_s = nt_traits.convert(this->s());
    m_alg_t = nt_traits.convert(this->t());
    m_alg_u = nt_traits.convert(this->u());
    m_alg_v = nt_traits.convert(this->v());
    m_alg_w = nt_traits.convert(this->w());

    // Set the generating conic ID for the source and target points.
    this->m_source.set_generating_conic(m_id);
    this->m_target.set_generating_conic(m_id);

    // Clear the _info bits.
    this->set_flag(Base::IS_VALID);
    this->reset_flag(Base::IS_FULL_CONIC);

    // Check if the arc is directed right (the target is lexicographically
    // greater than the source point), or to the left.
    Alg_kernel ker;
    Comparison_result dir_res =
      ker.compare_xy_2_object()(this->source(), this->target());

    CGAL_assertion(dir_res != EQUAL);

    if (dir_res == SMALLER) set_flag(IS_DIRECTED_RIGHT);

    // Compute the degree of the underlying conic.
    if (CGAL::sign(this->r()) != ZERO ||
        CGAL::sign(this->s()) != ZERO ||
        CGAL::sign(this->t()) != ZERO)
    {
      set_flag(DEGREE_2);

      if (this->m_orient == COLLINEAR) {
        set_flag(IS_SPECIAL_SEGMENT);

        // The arc is a vertical segment:
        if (ker.compare_x_2_object()(this->source(), this->target()) == EQUAL)
          set_flag(IS_VERTICAL_SEGMENT);
        return;
      }
    }
    else {
      CGAL_assertion(CGAL::sign(this->m_u) != ZERO ||
                     CGAL::sign(this->m_v) != ZERO);
      // The supporting curve is of the form: _u*x + _w = 0
      if (CGAL::sign(this->m_v) == ZERO) this->set_flag(IS_VERTICAL_SEGMENT);
      this->set_flag(DEGREE_1);
      return;
    }

    if (this->m_orient == COLLINEAR) return;

    // Compute a midpoint between the source and the target and get the y-value
    // of the arc at its x-coordiante.
    Point_2 p_mid =
      ker.construct_midpoint_2_object()(this->source(), this->target());
    Algebraic ys[2];
    CGAL_assertion_code(int n_ys = )
      this->_conic_get_y_coordinates(p_mid.x(), ys);

    CGAL_assertion(n_ys != 0);

    // Check which solution lies on the x-monotone arc.
    Point_2 p_arc_mid(p_mid.x(), ys[0]);

    if (this->_is_strictly_between_endpoints(p_arc_mid)) {
      // Mark that we should use the -sqrt(disc) root for points on this
      // x-monotone arc.
      this->reset_flag(PLUS_SQRT_DISC_ROOT);
    }
    else {
      CGAL_assertion(n_ys == 2);
      p_arc_mid = Point_2(p_mid.x(), ys[1]);

      CGAL_assertion(this->_is_strictly_between_endpoints(p_arc_mid));

      // Mark that we should use the +sqrt(disc) root for points on this
      // x-monotone arc.
      this->set_flag(PLUS_SQRT_DISC_ROOT);
    }

    // Check whether the conic is facing up or facing down:
    // Check whether the arc (which is x-monotone of degree 2) lies above or
    // below the segement that contects its two end-points (x1,y1) and (x2,y2).
    // To do that, we find the y coordinate of a point on the arc whose x
    // coordinate is (x1+x2)/2 and compare it to (y1+y2)/2.
    Comparison_result res = ker.compare_y_2_object()(p_arc_mid, p_mid);
    // The arc is above the connecting segment, so it is facing upwards.
    if (res == LARGER) this->set_flag(FACING_UP);
    // The arc is below the connecting segment, so it is facing downwards.
    else if (res == SMALLER) this->set_flag(FACING_DOWN);
  }

public:
  /// \name Public constrcutors, assignment operators, and destructors.
  //@{

  /*! Default constructor.
   */
  Conic_x_monotone_arc_2() : Base(), m_id() {}

  /*! Copy constructor.
   * \param arc The copied arc.
   */
  Conic_x_monotone_arc_2(const Self& arc) :
    Base(arc),
    m_alg_r(arc.m_alg_r),
    m_alg_s(arc.m_alg_s),
    m_alg_t(arc.m_alg_t),
    m_alg_u(arc.m_alg_u),
    m_alg_v(arc.m_alg_v),
    m_alg_w(arc.m_alg_w),
    m_id(arc.m_id)
  {}

  /*! Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc) {
    CGAL_precondition (arc.is_valid());

    if (this == &arc) return (*this);

    // Copy the base arc.
    Base::operator= (arc);

    // Set the rest of the properties.
    m_alg_r = arc.m_alg_r;
    m_alg_s = arc.m_alg_s;
    m_alg_t = arc.m_alg_t;
    m_alg_u = arc.m_alg_u;
    m_alg_v = arc.m_alg_v;
    m_alg_w = arc.m_alg_w;

    m_id = arc.m_id;

    return (*this);
  }
  //@}

private:
  template <typename, typename, typename> friend class Arr_conic_traits_2;

  /// \name private constrcutors to be used only by the traits class template.
  //@{

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \pre The given arc is x-monotone.
   */
  Conic_x_monotone_arc_2(const Base& arc) : Base(arc), m_id() {}

  /*! Construct an x-monotone arc from a conic arc.
   * \param arc The given (base) arc.
   * \param id The ID of the base arc.
   */
  Conic_x_monotone_arc_2(const Base& arc, const Conic_id& id) :
    Base(arc),
    m_id(id)
  {}
  //@}

public:
  /// \name Accessing the arc properties.
  //@{

  /*! Obtain the facing mask.
   */
  size_t facing_mask() const { return this->m_info & FACING_MASK; }

  /*! Obtain the degree mask.
   */
  size_t degree_mask() const { return this->m_info & DEGREE_MASK; }

  /*! Obtain the coefficients of the underlying conic.
   */
  const Integer& r() const { return ( this->m_r); }
  const Integer& s() const { return ( this->m_s); }
  const Integer& t() const { return ( this->m_t); }
  const Integer& u() const { return ( this->m_u); }
  const Integer& v() const { return ( this->m_v); }
  const Integer& w() const { return ( this->m_w); }

  /*! Obtain the left endpoint of the arc.
   */
  const Point_2& left() const {
    if (this->test_flag(IS_DIRECTED_RIGHT)) return this->source();
    else return this->target();
  }

  /*! Obtain the right endpoint of the arc.
   */
  const Point_2& right() const {
    if (this->test_flag(IS_DIRECTED_RIGHT)) return this->target();
    else return this->source();
  }

  /*! Determine whether the conic arc is directed iexicographically right.
   */
  bool is_directed_right() const { return this->test_flag(IS_DIRECTED_RIGHT); }

  /*! Determine whether the conic arc is a vertical segment.
   */
  bool is_vertical() const { return this->test_flag(IS_VERTICAL_SEGMENT); }

  /*! Determine whether the conic arc is a facing up.
   */
  bool is_upper() const { return this->test_flag(FACING_UP); }

  /*! Determine whether the conic arc is a facing down.
   */
  bool is_lower() const { return this->test_flag(FACING_DOWN); }

  /*! Check whether the arc is a special segment connecting two algebraic
   * endpoints (and has no undelying integer conic coefficients).
   */
  bool is_special_segment() const { return this->test_flag(IS_SPECIAL_SEGMENT); }

  /*! Obtain the mask of the DEGREE_1 flag.
   */
  static constexpr size_t degree_1_mask() { return flag_mask(DEGREE_1); }

  /*! Obtain the mask of the DEGREE_1 flag.
   */
  static constexpr size_t degree_2_mask() { return flag_mask(DEGREE_2); }

  /*! Obtain the algebraic coefficients.
   */
  Algebraic alg_r() const { return m_alg_r; }
  Algebraic alg_s() const { return m_alg_s; }
  Algebraic alg_t() const { return m_alg_t; }
  Algebraic alg_u() const { return m_alg_u; }
  Algebraic alg_v() const { return m_alg_v; }
  Algebraic alg_w() const { return m_alg_w; }

  /*! Obtain the conic id.
   */
  Conic_id id() const { return m_id; }

  /*! Check whether the given point lies on the arc.
   * \param p The qury point.
   * \param (true) if p lies on the arc; (false) otherwise.
   */
  CGAL_DEPRECATED
  bool contains_point(const Point_2& p) const {
    const auto& xcv = *this;
    // First check if p lies on the supporting conic. We first check whether
    // it is one of p's generating conic curves.
    bool p_on_conic(false);
    if (p.is_generating_conic(xcv.id())) p_on_conic = true;
    else {
      // Check whether p satisfies the supporting conic equation.
      p_on_conic = xcv.is_on_supporting_conic(p.x(), p.y());
      if (p_on_conic) {
        // As p lies on the supporting conic of our arc, add its ID to
        // the list of generating conics for p.
        Point_2& p_non_const = const_cast<Point_2&>(p);
        p_non_const.set_generating_conic(xcv.id());
      }
    }

    if (! p_on_conic) return false;

    // Check if p is between the endpoints of the arc.
    return _is_between_endpoints(p);
  }

  /*! Obtain a bounding box for the conic arc.
   * \return The bounding box.
   */
  CGAL_DEPRECATED Bbox_2 bbox() const { return Base::bbox(); }
  //@}

  // Setters
  //@{

  /*! Set the algebraic coefficients.
   */
  void set_alg_coefficients(const Algebraic& alg_r, const Algebraic& alg_s,
                            const Algebraic& alg_t, const Algebraic& alg_u,
                            const Algebraic& alg_v, const Algebraic& alg_w)
  {
    m_alg_r = alg_r;
    m_alg_s = alg_s;
    m_alg_t = alg_t;
    m_alg_u = alg_u;
    m_alg_v = alg_v;
    m_alg_w = alg_w;
  }

  /*! Add a generating conic ID.
   */
  void set_generating_conic(const Conic_id& id) {
    this->m_source.set_generating_conic(id);
    this->m_target.set_generating_conic(id);
  }
  //@}

  /*! Compute a point on an arc with the same \f$x\f$-coordiante as the given
   * point.
   * \param p The given point.
   * \pre The arc is not vertical and `p` is in the \f$x\f$-range of the arc.
   * \return A point on the arc with the same \f$x\f$-coordiante as `p`.
   */
  CGAL_DEPRECATED
  Point_2 point_at_x(const Point_2& p) const {
    const auto& xcv = *this;
    Alg_kernel alg_kernel;

    // Make sure that p is in the x-range of the arc.
    CGAL_precondition(! xcv.is_vertical());

    CGAL_precondition_code(auto cmp_x = alg_kernel.compare_x_2_object());
    CGAL_precondition((cmp_x(p, xcv.left()) != SMALLER) &&
                      (cmp_x(p, xcv.right()) != LARGER));

    if (xcv.is_special_segment()) {
      // In case of a special segment, the equation of the supported line
      // (a*x + b*y + c) = 0 is stored with the extra data field, and we
      // simply have:
      const auto& extra_data = xcv.extra_data();
      Algebraic y = -(extra_data->a*p.x() + extra_data->c) / extra_data->b;

      // Return the computed point.
      return Point_2(p.x(), y);
    }

    // Compute the y-coordinate according to the degree of the supporting
    // conic curve.
    typename Base::Nt_traits nt_traits;
    Algebraic y;

    if (xcv.degree_mask() == Self::degree_1_mask()) {
      // In case of a linear curve, the y-coordinate is a simple linear
      // expression of x(p) (note that v is not 0 as the arc is not vertical):
      //   y = -(u*x(p) + w) / v
      y = -(xcv.alg_u()*p.x() + xcv.alg_w()) / xcv.alg_v();
    }
    else if (xcv.orientation() == COLLINEAR) {
      const auto& extra_data = xcv.extra_data();
      CGAL_assertion(extra_data != nullptr);

      // In this case the equation of the supporting line is given by the
      // extra data structure.
      y = -(extra_data->a * p.x() + extra_data->c) / extra_data->b;
    }
    else {
      CGAL_assertion(xcv.degree_mask() == Self::degree_2_mask());

      // In this case the y-coordinate is one of solutions to the quadratic
      // equation:
      //  s*y^2 + (t*x(p) + v)*y + (r*x(p)^2 + u*x(p) + w) = 0
      Algebraic A = xcv.alg_s();
      Algebraic B = xcv.alg_t()*p.x() + xcv.alg_v();
      Algebraic C = (xcv.alg_r()*p.x() + xcv.alg_u())*p.x() + xcv.alg_w();

      if (CGAL::sign(xcv.s()) == ZERO) {
        // In this case A is 0 and we have a linear equation.
        CGAL_assertion(CGAL::sign(B) != ZERO);

        y = -C / B;
      }
      else {
        // Solve the quadratic equation.
        Algebraic disc = B*B - 4*A*C;

        CGAL_assertion(CGAL::sign(disc) != NEGATIVE);

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        y = (xcv.test_flag(Self::PLUS_SQRT_DISC_ROOT)) ?
          (nt_traits.sqrt(disc) - B) / (2*A) :
          -(B + nt_traits.sqrt(disc)) / (2*A);
      }
    }

    // Return the computed point.
    return Point_2(p.x(), y);
  }

  /// \name Constructing points on the arc.
  //@{

  /*! Obtain a polyline approximating the conic arc.
   * \param n The maximal number of sample points.
   * \param oi An output iterator, whose value-type is pair<double,double>
   *           (representing an approximated point).
   *           In case the arc is a line segment, there are 2 output points,
   *           otherwise the arc is approximated by the polyline defined by
   *           (p_0, p_1, ..., p_n), where p_0 and p_n are the left and right
   *           endpoints of the arc, respectively.
   */
  template <typename OutputIterator>
  OutputIterator polyline_approximation(size_t n, OutputIterator oi) const {
    CGAL_precondition (n != 0);

    const double x_left = CGAL::to_double (left().x());
    const double y_left = CGAL::to_double (left().y());
    const double x_right = CGAL::to_double (right().x());
    const double y_right = CGAL::to_double (right().y());

    if (this->m_orient == COLLINEAR) {
      // In case of a line segment, return the two endpoints.
      *oi++ = std::pair<double, double>(x_left, y_left);
      *oi++ = std::pair<double, double>(x_right, y_right);
      return oi;
    }

    // Otherwise, sample (n - 1) equally-spaced points in between.
    const double app_r = CGAL::to_double(this->m_r);
    const double app_s = CGAL::to_double(this->m_s);
    const double app_t = CGAL::to_double(this->m_t);
    const double app_u = CGAL::to_double(this->m_u);
    const double app_v = CGAL::to_double(this->m_v);
    const double app_w = CGAL::to_double(this->m_w);
    const double x_jump = (x_right - x_left) / n;
    double x, y;
    const bool A_is_zero = (CGAL::sign(this->m_s) == ZERO);
    double A = app_s, B, C;
    double disc;
    size_t i;

    *oi++ = std::pair<double, double>(x_left, y_left);      // the left point
    for (i = 1; i < n; ++i) {
      x = x_left + x_jump*i;

      // Solve the quadratic equation: A*x^2 + B*x + C = 0:
      B = app_t*x + app_v;
      C = (app_r*x + app_u)*x + app_w;

      if (A_is_zero) y = -C / B;
      else {
        disc = B*B - 4*A*C;
        if (disc < 0) disc = 0;

        // We take either the root involving -sqrt(disc) or +sqrt(disc)
        // based on the information flags.
        y = (test_flag(PLUS_SQRT_DISC_ROOT)) ?
          (std::sqrt(disc) - B) / (2*A) : -(B + std::sqrt (disc)) / (2*A);
      }
      *oi++ = std::pair<double, double>(x, y);
    }
    *oi++ = std::pair<double, double>(x_right, y_right);    // the right point

    return oi;
  }

  //@}

  /// \name Constructing x-monotone arcs.
  //@{

  /*! Flip the arc.
   * \return An arc with swapped source and target and a reverse orienation.
   */
  Self flip() const {
    // Make a copy of the current arc.
    Self arc = *this;

    // Reverse the orientation.
    if (this->m_orient == CLOCKWISE) arc.m_orient = COUNTERCLOCKWISE;
    else if (this->m_orient == COUNTERCLOCKWISE) arc.m_orient = CLOCKWISE;

    // Swap the source and the target.
    arc.m_source = this->m_target;
    arc.m_target = this->m_source;

    // Change the direction bit among the information flags.
    arc.flip_flag(IS_DIRECTED_RIGHT);

    return arc;
  }

  /*! Trim the arc given its new endpoints.
   * \param ps The new source point.
   * \param pt The new target point.
   * \return The new trimmed arc.
   * \pre Both ps and pt lies on the arc and must conform with the current
   *      direction of the arc.
   */
  CGAL_DEPRECATED
  Self trim(const Point_2& ps, const Point_2& pt) const {
    auto& xcv = *this;
    // Make sure that both ps and pt lie on the arc.
    CGAL_precondition(contains_point(ps) && contains_point(pt));

    // Make sure that the endpoints conform with the direction of the arc.
    Alg_kernel alg_kernel;
    Self res_xcv = xcv;          // make a copy of the current arc
    auto eq = alg_kernel.equal_2_object();
      auto set_source = [&](const Point_2 ps)->void {
                          if (! eq(ps, xcv.source())) {
                            res_xcv.set_source(ps);
                            if (! ps.is_generating_conic(xcv.id()))
                              res_xcv.source().set_generating_conic(xcv.id());
                          }
                        };
      auto set_target = [&](const Point_2 pt)->void {
                          if (! eq(pt, xcv.target())) {
                            res_xcv.set_target(pt);
                            if (! pt.is_generating_conic(xcv.id()))
                              res_xcv.target().set_generating_conic(xcv.id());
                          }
                        };

    auto cmp_xy = alg_kernel.compare_xy_2_object();
    auto res = cmp_xy(ps, pt);
    CGAL_assertion(res != EQUAL);
    if ((xcv.test_flag(Self::IS_DIRECTED_RIGHT) && (res == LARGER)) ||
        (! xcv.test_flag(Self::IS_DIRECTED_RIGHT) && (res == SMALLER))) {
      set_source(pt);
      set_target(ps);
    }
    else {
      set_source(ps);
      set_target(pt);
    }

    return res_xcv;
  }

  //@}

  /*! Compare two arcs immediately to the leftt of their intersection point.
   * \param xcv1 The first compared arc.
   * \param xcv2 The second compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the left of `p`.
   * \pre Both arcs we compare are not vertical segments.
   */
  CGAL_DEPRECATED
  Comparison_result compare_to_left(const Self& xcv2, const Point_2& p) const {
    const auto& xcv1 = *this;
    CGAL_precondition(! xcv1.is_vertical() && ! xcv2.is_vertical());

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic(xcv2)) {
      if (xcv1.test_flag(Self::FACING_UP) &&
          xcv2.test_flag(Self::FACING_DOWN))
        return LARGER;
      else if (xcv1.test_flag(Self::FACING_DOWN) &&
               xcv2.test_flag(Self::FACING_UP))
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion(xcv1.facing_mask() == xcv2.facing_mask());

      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    xcv1.derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    xcv2.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign (slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign (slope2_denom) == ZERO);

    if (! is_vertical_slope1 && ! is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                                  slope1_numer*slope2_denom);

      if (slope_res != EQUAL) return slope_res;

      // Use the second-order derivative.
      xcv1.derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res = CGAL::compare(slope1_numer*slope2_denom,
                                slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      xcv1.derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res = CGAL::compare(slope2_numer*slope1_denom,
                                slope1_numer*slope2_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      return slope_res;
    }
    else if (! is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv1.facing_mask() != 0);

      return (xcv1.test_flag(Self::FACING_UP)) ?
        LARGER : SMALLER;
    }
    else if (! is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv2.facing_mask() != 0);

      return (xcv2.test_flag(Self::FACING_UP)) ?
        SMALLER : LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if (xcv1.test_flag(Self::FACING_UP) &&
        xcv2.test_flag(Self::FACING_DOWN))
      return LARGER;
    else if (xcv1.test_flag(Self::FACING_DOWN) &&
             xcv2.test_flag(Self::FACING_UP))
      return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    xcv1.derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    xcv2.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      xcv1.derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    // Check whether both are facing up.
    if (xcv1.test_flag(Self::FACING_UP) &&
        xcv2.test_flag(Self::FACING_UP))
      return ((slope_res == LARGER) ? SMALLER : LARGER);

    // Both are facing down.
    return slope_res;
  }

  /*! Compare two arcs immediately to the right of their intersection point.
   * \param xcv1 The first compared arc.
   * \param xcv2 The second compared arc.
   * \param p The reference intersection point.
   * \return The relative position of the arcs to the right of `p`.
   * \pre Both arcs we compare are not vertical segments.
   */
  CGAL_DEPRECATED
  Comparison_result compare_to_right(const Self& xcv2, const Point_2& p) const {
    const auto& xcv1 = *this;
    CGAL_precondition(! xcv1.is_vertical() && ! xcv2.is_vertical());

    // In case one arc is facing upwards and another facing downwards, it is
    // clear that the one facing upward is above the one facing downwards.
    if (_has_same_supporting_conic(xcv2)) {
      if (xcv1.test_flag(Self::FACING_UP) &&
          xcv2.test_flag(Self::FACING_DOWN))
        return LARGER;
      else if (xcv1.test_flag(Self::FACING_DOWN) &&
               xcv2.test_flag(Self::FACING_UP))
        return SMALLER;

      // In this case the two arcs overlap.
      CGAL_assertion(xcv1.facing_mask() == xcv2.facing_mask());
      return EQUAL;
    }

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    xcv1.derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    xcv2.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    // Check if any of the slopes is vertical.
    const bool is_vertical_slope1 = (CGAL::sign(slope1_denom) == ZERO);
    const bool is_vertical_slope2 = (CGAL::sign(slope2_denom) == ZERO);

    if (! is_vertical_slope1 && ! is_vertical_slope2) {
      // The two derivatives at p are well-defined: use them to determine
      // which arc is above the other (the one with a larger slope is below).
      Comparison_result slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the second-order derivative.
      xcv1.derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 2, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      if (slope_res != EQUAL) return (slope_res);

      // Use the third-order derivative.
      xcv1.derive_by_x_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_x_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

      // \todo Handle higher-order derivatives:
      CGAL_assertion(slope_res != EQUAL);

      return slope_res;
    }
    else if (! is_vertical_slope2) {
      // The first arc has a vertical slope at p: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv1.facing_mask() != 0);

      return (xcv1.test_flag(Self::FACING_UP)) ? LARGER : SMALLER;
    }
    else if (! is_vertical_slope1) {
      // The second arc has a vertical slope at p_int: check whether it is
      // facing upwards or downwards and decide accordingly.
      CGAL_assertion(xcv2.facing_mask() != 0);

      return (xcv2.test_flag(Self::FACING_UP)) ? SMALLER : LARGER;
    }

    // The two arcs have vertical slopes at p_int:
    // First check whether one is facing up and one down. In this case the
    // comparison result is trivial.
    if (xcv1.test_flag(Self::FACING_UP) &&
        xcv2.test_flag(Self::FACING_DOWN)) return LARGER;
    else if (xcv1.test_flag(Self::FACING_DOWN) &&
             xcv2.test_flag(Self::FACING_UP)) return SMALLER;

    // Compute the second-order derivative by y and act according to it.
    xcv1.derive_by_y_at(p, 2, slope1_numer, slope1_denom);
    xcv2.derive_by_y_at(p, 2, slope2_numer, slope2_denom);

    Comparison_result slope_res =
      CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom);

    // If necessary, use the third-order derivative by y.
    if (slope_res == EQUAL) {
      // \todo Check this!
      xcv1.derive_by_y_at(p, 3, slope1_numer, slope1_denom);
      xcv2.derive_by_y_at(p, 3, slope2_numer, slope2_denom);

      slope_res =
        CGAL::compare(slope2_numer*slope1_denom, slope1_numer*slope2_denom);
    }

    // \todo Handle higher-order derivatives:
    CGAL_assertion(slope_res != EQUAL);

    if (xcv1.test_flag(Self::FACING_UP) &&
        xcv2.test_flag(Self::FACING_UP))
      return (slope_res == LARGER) ? SMALLER : LARGER;  // both are facing up
    return slope_res;                                   // both are facing down
  }

  /*! Check whether two arcs are equal (have the same graph).
   * \param xcv1 The first compared arc.
   * \param xcv2 The second compared arc.
   * \return `true` if the two arcs have the same graph; `false` otherwise.
   */
  CGAL_DEPRECATED
  bool equals(const Self& xcv2) const {
    const auto& xcv1 = *this;
    Alg_kernel alg_kernel;

    // The two arc must have the same supporting conic curves.
    if (! _has_same_supporting_conic(xcv2)) return false;

    auto eq = alg_kernel.equal_2_object();

    // Check that the arc endpoints are the same.
    if (xcv1.orientation() == COLLINEAR) {
      CGAL_assertion(xcv2.orientation() == COLLINEAR);
      return((eq(xcv1.source(), xcv2.source()) &&
              eq(xcv1.target(), xcv2.target())) ||
             (eq(xcv1.source(), xcv2.target()) &&
              eq(xcv1.target(), xcv2.source())));
    }

    if (xcv1.orientation() == xcv2.m_orient) {
      // Same orientation - the source and target points must be the same.
      return (eq(xcv1.source(), xcv2.source()) &&
              eq(xcv1.target(), xcv2.target()));
    }

    // Reverse orientation - the source and target points must be swapped.
    return (eq(xcv1.source(), xcv2.target()) &&
            eq(xcv1.target(), xcv2.source()));
  }

  /*! Check whether it is possible to merge the arc with the given arc.
   * \param xcv1 The first arc.
   * \param xcv2 The second arc.
   * \return `true` if it is possible to merge the two arcs;
   *         `false` otherwise.
   */
  CGAL_DEPRECATED
  bool can_merge_with(const Self& xcv2) const {
    const auto& xcv1 = *this;
    Alg_kernel alg_kernel;

    // In order to merge the two arcs, they should have the same supporting
    // conic.
    if (! _has_same_supporting_conic(xcv2)) return false;

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    auto eq = alg_kernel.equal_2_object();
    return (eq(xcv1.right(), xcv2.left()) || eq(xcv1.left(), xcv2.right()));
  }

  /*! Merge the current arc with the given arc.
   * \param xcv1 The first arc to merge with.
   * \param xcv2 The second arc to merge with.
   * \pre The two arcs are mergeable.
   */
  CGAL_DEPRECATED
  void merge(const Self& xcv2) const {
    const auto& xcv1 = *this;
    Alg_kernel alg_kernel;

    // Check whether we should extend the arc to the left or to the right.
    auto eq = alg_kernel.equal_2_object();
    if (eq(xcv1.right(), xcv2.left())) {
      // Extend the arc to the right.
      if (xcv1.test_flag(Self::IS_DIRECTED_RIGHT))
        xcv1.set_target(xcv2.right());
      else xcv1.set_source(xcv2.right());
    }
    else {
      CGAL_precondition(eq(xcv1.left(), xcv2.right()));

      // Extend the arc to the left.
      if (xcv1.test_flag(Self::IS_DIRECTED_RIGHT))
        xcv1.set_source(xcv2.left());
      else xcv1.set_target(xcv2.left());
    }
  }

private:
  /// \name Auxiliary (private) functions.
  //@{

  /*! Check whether the given point lies on the supporting conic of the arc.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return (true) if p lies on the supporting conic; (false) otherwise.
   */
  bool is_on_supporting_conic(const Algebraic& px, const Algebraic& py) const {
    CGAL::Sign my_sign = (! is_special_segment()) ?
      // Check whether p satisfies the conic equation.
      // The point must satisfy: r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0.
      CGAL::sign((m_alg_r*px + m_alg_t*py + m_alg_u) * px +
                 (m_alg_s*py + m_alg_v) * py + m_alg_w) :
      // Check whether p satisfies the equation of the line stored with the
      // extra data.
      sign_of_extra_data(px, py);
    return (my_sign == ZERO);
  }

public:
  /*! Obtain the i-th order derivative by x of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void derive_by_x_at(const Alg_point_2& p, const unsigned int& i,
                      Algebraic& slope_numer, Algebraic& slope_denom) const {
    if (is_special_segment()) {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -a/b. The higher-order
      // derivatives are all 0.
      if (i == 1) {
        if (CGAL::sign(this->m_extra_data->b) != NEGATIVE) {
          slope_numer = - this->m_extra_data->a;
          slope_denom = this->m_extra_data->b;
        }
        else {
          slope_numer = this->m_extra_data->a;
          slope_denom = - this->m_extra_data->b;
        }
      }
      else {
        slope_numer = 0;
        slope_denom = 1;
      }

      return;
    }

    // The derivative by x of the conic
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2r*x + t*y + u       alpha
    //   y' = - ---------------- = - -------
    //           2s*y + t*x + v       beta
    //
    const Algebraic two = 2;
    const Algebraic sl_numer = two*m_alg_r*p.x() + m_alg_t*p.y() + m_alg_u;
    const Algebraic sl_denom = two*m_alg_s*p.y() + m_alg_t*p.x() + m_alg_v;

    if (i == 1) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second-order derivative is given by:
    //
    //             s*alpha^2 - t*alpha*beta + r*beta^2     gamma
    //   y'' = -2 ------------------------------------- = -------
    //                           beta^3                    delta
    //
    const Algebraic sl2_numer = m_alg_s * sl_numer*sl_numer -
      m_alg_t * sl_numer*sl_denom + m_alg_r * sl_denom*sl_denom;
    const Algebraic sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else {
        slope_numer = two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // The third-order derivative is given by:
    //
    //              (2s*alpha - t*beta) * gamma
    //   y''' = -6 ------------------------------
    //                    beta^2 * delta
    //
    const Algebraic sl3_numer =
      (two * m_alg_s * sl_numer - m_alg_t * sl_denom) * sl2_numer;
    const Algebraic sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
  }

  /*! Obtain the i-th order derivative by y of the conic at the point p=(x,y).
   * \param p The point where we derive.
   * \param i The order of the derivatives (either 1, 2 or 3).
   * \param slope_numer The numerator of the slope.
   * \param slope_denom The denominator of the slope.
   * \todo Allow higher order derivatives.
   */
  void derive_by_y_at(const Alg_point_2& p, const int& i,
                      Algebraic& slope_numer, Algebraic& slope_denom) const {
    if (is_special_segment()) {
      // Special treatment for special segments, given by (a*x + b*y + c = 0),
      // so their first-order derivative by x is simply -b/a. The higher-order
      // derivatives are all 0.
      if (i == 1) {
        if (CGAL::sign(this->m_extra_data->a) != NEGATIVE) {
          slope_numer = - this->m_extra_data->b;
          slope_denom = this->m_extra_data->a;
        }
        else {
          slope_numer = this->m_extra_data->b;
          slope_denom = - this->m_extra_data->a;
        }
      }
      else {
        slope_numer = 0;
        slope_denom = 1;
      }

      return;
    }

    // The derivative by y of the conic
    //   C: {r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0}
    // at the point p=(x,y) is given by:
    //
    //           2s*y + t*x + v     alpha
    //   x' = - ---------------- = -------
    //           2r*x + t*y + u      beta
    //
    const Algebraic two = 2;
    const Algebraic sl_numer = two*m_alg_s*p.y() + m_alg_t*p.x() + m_alg_v;
    const Algebraic sl_denom = two*m_alg_r*p.x() + m_alg_t*p.y() + m_alg_u;

    if (i == 1) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -sl_numer;
        slope_denom = sl_denom;
      }
      else {
        slope_numer = sl_numer;
        slope_denom = -sl_denom;
      }

      return;
    }

    // The second-order derivative is given by:
    //
    //             r*alpha^2 - t*alpha*beta + s*beta^2
    //   x'' = -2 -------------------------------------
    //                           beta^3
    //
    const Algebraic sl2_numer = m_alg_r * sl_numer*sl_numer -
      m_alg_t * sl_numer*sl_denom + m_alg_s * sl_denom*sl_denom;
    const Algebraic sl2_denom = sl_denom*sl_denom*sl_denom;

    if (i == 2) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign (sl_denom) != NEGATIVE) {
        slope_numer = -two *sl2_numer;
        slope_denom = sl2_denom;
      }
      else {
        slope_numer = two *sl2_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // The third-order derivative is given by:
    //
    //              (2t*alpha - t*beta) * gamma
    //   y''' = -6 ------------------------------
    //                    beta^2 * delta
    //
    const Algebraic sl3_numer =
      (two * m_alg_r * sl_numer - m_alg_t * sl_denom) * sl2_numer;
    const Algebraic sl3_denom = sl_denom*sl_denom * sl2_denom;

    if (i == 3) {
      // Make sure that the denominator is always positive.
      if (CGAL::sign(sl_denom) != NEGATIVE) {
        slope_numer = -6 * sl3_numer;
        slope_denom = sl3_denom;
      }
      else {
        slope_numer = 6 * sl3_numer;
        slope_denom = -sl2_denom;
      }

      return;
    }

    // \todo Handle higher-order derivatives as well.
    CGAL_error();
  }

private:
  /*! Compute the multiplicity of an intersection point.
   * \param arc The arc to intersect with.
   * \param p The intersection point.
   * \return The multiplicity of the intersection point.
   */
  size_t multiplicity_of_intersection_point(const Self& xcv,
                                            const Alg_point_2& p) const {
    CGAL_assertion(! is_special_segment() || ! xcv.is_special_segment());

    // Compare the slopes of the two arcs at p, using their first-order
    // partial derivatives.
    Algebraic slope1_numer, slope1_denom;
    Algebraic slope2_numer, slope2_denom;

    derive_by_x_at(p, 1, slope1_numer, slope1_denom);
    xcv.derive_by_x_at(p, 1, slope2_numer, slope2_denom);

    if (CGAL::compare(slope1_numer*slope2_denom, slope2_numer*slope1_denom) !=
        EQUAL) {
      // Different slopes at p - the mutiplicity of p is 1:
      return 1;
    }

    if (CGAL::sign(slope1_denom) != ZERO &&
        CGAL::sign(slope2_denom) != ZERO) {
      // The curves do not have a vertical slope at p.
      // Compare their second-order derivative by x:
      derive_by_x_at(p, 2, slope1_numer, slope1_denom);
      xcv.derive_by_x_at(p, 2, slope2_numer, slope2_denom);
    }
    else {
      // Both curves have a vertical slope at p.
      // Compare their second-order derivative by y:
      derive_by_y_at(p, 2, slope1_numer, slope1_denom);
      xcv.derive_by_y_at(p, 2, slope2_numer, slope2_denom);
    }

    if (CGAL::compare(slope1_numer*slope2_denom,
                      slope2_numer*slope1_denom) != EQUAL)
    {
      // Different curvatures at p - the mutiplicity of p is 2:
      return 2;
    }

    // If we reached here, the multiplicity of the intersection point is 3:
    return 3;
  }
  //@}

};

/*! Exporter for x-monotone conic arcs.
 */
template <typename Conic_arc_2>
std::ostream& operator<<(std::ostream& os,
                         const Conic_x_monotone_arc_2<Conic_arc_2>& xcv)
{
  // Output the supporting conic curve.
  os << "{" << CGAL::to_double(xcv.r()) << "*x^2 + "
     << CGAL::to_double(xcv.s()) << "*y^2 + "
     << CGAL::to_double(xcv.t()) << "*xy + "
     << CGAL::to_double(xcv.u()) << "*x + "
     << CGAL::to_double(xcv.v()) << "*y + "
     << CGAL::to_double(xcv.w()) << "}";

  // Output the endpoints.
  os << " : (" << CGAL::to_double(xcv.source().x()) << ","
     << CGAL::to_double(xcv.source().y()) << ") ";

  if (xcv.orientation() == CLOCKWISE) os << "--cw-->";
  else if (xcv.orientation() == COUNTERCLOCKWISE) os << "--ccw-->";
  else os << "--l-->";

  os << " (" << CGAL::to_double(xcv.target().x()) << ","
     << CGAL::to_double(xcv.target().y()) << ")";

  return os;
}

} //namespace CGAL

#endif
