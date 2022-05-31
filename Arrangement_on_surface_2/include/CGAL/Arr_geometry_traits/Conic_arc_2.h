// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_CONIC_ARC_2_H
#define CGAL_CONIC_ARC_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Header file for the _Conic_arc_2<Int_kernel, Alg_kernel, Nt_traits> class.
 */

#include <ostream>

#include <CGAL/Arr_geometry_traits/Conic_point_2.h>
#include <CGAL/Arr_geometry_traits/Conic_intersections_2.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {

/*! Representation of a conic arc -- a bounded segment that lies on a conic
 * curve, the loci of all points satisfying the equation:
 *   r*x^2 + s*y^2 + t*xy + u*x + v*y +w = 0
 *
 * The class is templated with three parameters:
 * Rat_kernel A kernel that provides the input objects or coefficients.
 *            Rat_kernel::FT should be an integral or a rational type.
 * Alg_kernel A geometric kernel, where Alg_kernel::FT is the number type
 *            for the coordinates of arrangement vertices, which are algebraic
 *            numbers of degree up to 4 (preferably it is CORE::Expr).
 * Nt_traits A traits class for performing various operations on the integer,
 *           rational and algebraic types.
 */

template <typename RatKernel, typename AlgKernel, typename NtTraits>
class _Conic_arc_2 {
public:
  typedef RatKernel                                        Rat_kernel;
  typedef AlgKernel                                        Alg_kernel;
  typedef NtTraits                                         Nt_traits;

  typedef _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>  Self;

  typedef typename Rat_kernel::FT                          Rational;
  typedef typename Rat_kernel::Point_2                     Rat_point_2;
  typedef typename Rat_kernel::Segment_2                   Rat_segment_2;
  typedef typename Rat_kernel::Circle_2                    Rat_circle_2;

  typedef typename Nt_traits::Integer                      Integer;

  typedef typename Alg_kernel::FT                          Algebraic;
  typedef typename Alg_kernel::Point_2                     Point_2;
  typedef _Conic_point_2<Alg_kernel>                       Conic_point_2;

  /*! \struct
   * For arcs whose base is a hyperbola we store the axis (a*x + b*y + c = 0)
   * which separates the two bracnes of the hyperbola. We also store the side
   * (NEGATIVE or POSITIVE) that the arc occupies.
   * In case of line segments connecting two algebraic endpoints, we use this
   * structure two store the coefficients of the line supporting this segment.
   * In this case we set the side field to be ZERO.
   */
  struct Extra_data {
    Algebraic a;
    Algebraic b;
    Algebraic c;
    Sign side;
  };

  // Bit masks for the m_info field.
  enum {
    IS_VALID = 0,
    IS_FULL_CONIC,
    LAST_INFO,
  };

protected:
  Integer m_r;          //
  Integer m_s;          // The coefficients of the supporting conic curve:
  Integer m_t;          //
  Integer m_u;          //
  Integer m_v;          //   r*x^2 + s*y^2 + t*xy + u*x + v*y + w = 0
  Integer m_w;          //

  Orientation m_orient;         // The orientation of the conic.
  int m_info;                   // does the arc represent a full conic curve.
  Conic_point_2 m_source;       // the source of the arc (if not a full curve).
  Conic_point_2 m_target;       // the target of the arc (if not a full curve).
  Extra_data* m_extra_data;     // The extra data stored with the arc
                                // (may be nullptr).

public:
  /// \name Construction and destruction functions.
  //@{

  /*! Default constructor.
   */
  _Conic_arc_2() :
    m_r(0), m_s(0), m_t(0), m_u(0), m_v(0), m_w(0),
    m_orient(COLLINEAR),
    m_info(0),
    m_extra_data(nullptr)
  {}

  /*! Copy constructor.
   * \param arc The copied arc.
   */
  _Conic_arc_2(const Self& arc) :
    m_r(arc.m_r), m_s(arc.m_s), m_t(arc.m_t),
    m_u(arc.m_u), m_v(arc.m_v), m_w(arc.m_w),
    m_orient(arc.m_orient),
    m_info(arc.m_info),
    m_source(arc.m_source),
    m_target(arc.m_target)
  {
    m_extra_data = (arc.m_extra_data != nullptr) ?
      new Extra_data(*(arc.m_extra_data)) : nullptr;
  }

  /*! Destructor.
   */
  virtual ~_Conic_arc_2() { if (m_extra_data != nullptr) delete m_extra_data; }

  /*! Assignment operator.
   * \param arc The copied arc.
   */
  const Self& operator=(const Self& arc) {
    if (this == &arc) return (*this);

    // Free any existing data.
    if (m_extra_data != nullptr) delete m_extra_data;

    // Copy the arc's attributes.
    m_r = arc.m_r;
    m_s = arc.m_s;
    m_t = arc.m_t;
    m_u = arc.m_u;
    m_v = arc.m_v;
    m_w = arc.m_w;

    m_orient = arc.m_orient;
    m_info = arc.m_info;
    m_source = arc.m_source;
    m_target = arc.m_target;

    // Duplicate the extra data, if necessary.
    m_extra_data = (arc.m_extra_data != nullptr) ?
      new Extra_data(*(arc.m_extra_data)) : nullptr;

    return (*this);
  }
  //@}

  /// \name Get the arc properties.
  //@{

  /*! Determine wheather the arc is valid.
   */
  bool is_valid() const { return test_flag(IS_VALID); }

  /*! Determine whether the arc represents a full conic curve.
   */
  bool is_full_conic() const { return test_flag(IS_FULL_CONIC); }

  /*! Obtain the coefficients of the underlying conic.
   */
  const Integer& r() const { return (m_r); }
  const Integer& s() const { return (m_s); }
  const Integer& t() const { return (m_t); }
  const Integer& u() const { return (m_u); }
  const Integer& v() const { return (m_v); }
  const Integer& w() const { return (m_w); }

  /*! Obtain the arc's source.
   * \return The source point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& source() const {
    CGAL_precondition(! is_full_conic());
    return m_source;
  }

  /*! Obtain the arc's target.
   * \return The target point.
   * \pre The arc does not represent a full conic curve.
   */
  const Point_2& target() const {
    CGAL_precondition(! is_full_conic());
    return m_target;
  }

  /*! Obtain the orientation of the arc.
   * \return The orientation.
   */
  Orientation orientation() const { return m_orient; }

  /*! Obtain the extra data.
   */
  const Extra_data* extra_data() const { return m_extra_data; }
  //@}

protected:
  template <typename, typename, typename> friend class Arr_conic_traits_2;

  /// \name Flag manipulation functions.
  //@{
  template <typename T>
  static constexpr size_t flag_mask(const T flag) { return 0x1 << flag; }

  void reset_flags() { m_info = 0; }

  template <typename T>
  void set_flag(const T flag) { m_info |= flag_mask(flag); }

  template <typename T>
  void reset_flag(const T flag) { m_info &= ~flag_mask(flag); }

  template <typename T>
  void flip_flag(const T flag) { m_info ^= flag_mask(flag);}

  template <typename T>
  bool test_flag(const T flag) const
  { return (m_info & flag_mask(flag)); }
  //@}

public:
  /// \name Modifying functions (setters);
  // only friends have the priviledge to use.
  //@{

  /*! Set the source point of the conic arc.
   * \param ps The new source point.
   */
  void set_source(const Point_2& ps) { m_source = ps; }

  /*! Set the target point of the conic arc.
   * \param pt The new source point.
   */
  void set_target(const Point_2& pt) { m_target = pt; }

  /*! Set the coefficients.
   */
  void set_coefficients(Integer r, Integer s, Integer t,
                        Integer u, Integer v, Integer w) {

    m_r = r;
    m_s = s;
    m_t = t;
    m_u = u;
    m_v = v;
    m_w = w;
  }

  /*! Set the orientation.
   */
  void set_orientation(Orientation orient) { m_orient = orient; }

  /*! Set the endpoints.
   */
  void set_endpoints(const Conic_point_2& source, const Conic_point_2& target) {
    m_source = source;
    m_target = target;
  }

  /*! Set the extra data field.
   */
  void set_extra_data(Extra_data* extra_data) { m_extra_data = extra_data; }
  //@}

public:
  /// \name Auxiliary functions.
  //@{

  /*! Evaluate the sign of (a*x + b*y + c) stored with the extra data field
   * at a given point.
   * \param px The x-coordinate of query point.
   * \param py The y-coordinate of query point.
   * \return The sign of (a*x + b*y + c).
   */
  Sign sign_of_extra_data(const Algebraic& px, const Algebraic& py) const {
    CGAL_assertion(m_extra_data != nullptr);
    if (m_extra_data == nullptr) return ZERO;
    Algebraic val = m_extra_data->a*px + m_extra_data->b*py + m_extra_data->c;
    return CGAL::sign(val);
  }

  //@}
};

/*! Exporter for conic arcs.
 */
template <typename Rat_kernel, typename Alg_kernel, typename Nt_traits>
std::ostream&
operator<< (std::ostream& os,
            const _Conic_arc_2<Rat_kernel, Alg_kernel, Nt_traits>& arc)
{
  os << "{" << CGAL::to_double(arc.r()) << "*x^2 + "
     << CGAL::to_double(arc.s()) << "*y^2 + "
     << CGAL::to_double(arc.t()) << "*xy + "
     << CGAL::to_double(arc.u()) << "*x + "
     << CGAL::to_double(arc.v()) << "*y + "
     << CGAL::to_double(arc.w()) << "}";

  if (arc.is_full_conic()) os << " - Full curve";
  else {
    os << " : (" << CGAL::to_double(arc.source().x()) << ","
       << CGAL::to_double(arc.source().y()) << ") ";

    if (arc.orientation() == CLOCKWISE)
      os << "--cw-->";
    else if (arc.orientation() == COUNTERCLOCKWISE)
      os << "--ccw-->";
    else
      os << "--l-->";

    os << " (" << CGAL::to_double(arc.target().x()) << ","
       << CGAL::to_double(arc.target().y()) << ")";
  }

  return (os);
}

} //namespace CGAL

#endif
