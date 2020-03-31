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
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Baruch Zukerman   <baruchzu@post.tau.ac.il>
//            Waqar Khan        <wkhan@mpi-inf.mpg.de>
//            Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_CIRCLE_SEGMENT_2_H
#define CGAL_CIRCLE_SEGMENT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Header file for the _Circle_segment_2<Kernel, Filter> class.
 */
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Handle_for.h>
#include <list>
#include <map>
#include <ostream>

namespace CGAL {

// Forward declaration:
template <typename NumberType_, bool Filter_> class _One_root_point_2;

/*! \class
 * Representation of a point whose coordinates are one-root numbers.
 */
template <typename NumberType_, bool Filter_>
class _One_root_point_2_rep {
  friend class _One_root_point_2<NumberType_, Filter_>;

public:
  typedef NumberType_                               NT;
  typedef _One_root_point_2_rep<NT, Filter_>        Self;
  typedef Sqrt_extension<NT, NT, Tag_true,Boolean_tag<Filter_> >    CoordNT;

private:
  CoordNT _x;            // The coordinates.
  CoordNT _y;

public:
  /*! Default constructor. */
  _One_root_point_2_rep() :
    _x(0),
    _y(0)
  {}

  /*! Constructor of a point with one-root coefficients.
     This constructor of a point can also be used with rational coefficients
     thanks to convertor of CoordNT. */
  _One_root_point_2_rep(const CoordNT& x, const CoordNT& y) :
    _x(x),
    _y(y)
  {}
};

/*! \class
 * A handle for a point whose coordinates are one-root numbers.
 */
template <typename NumberType_, bool Filter_>
class _One_root_point_2 :
  public Handle_for<_One_root_point_2_rep<NumberType_, Filter_> >
{
public:
  typedef NumberType_                           NT;
  typedef _One_root_point_2<NT, Filter_>        Self;

private:
  typedef _One_root_point_2_rep<NT, Filter_>    Point_rep;
  typedef Handle_for<Point_rep>                 Point_handle;

public:
  typedef typename Point_rep::CoordNT           CoordNT;

  /*! Default constructor. */
  _One_root_point_2() : Point_handle(Point_rep()) {}

  /*! Copy constructor. */
  _One_root_point_2(const Self& p) : Point_handle(p) {}

  _One_root_point_2& operator=(const _One_root_point_2&)=default;

  /*! Constructor of a point with one-root coefficients.
     This constructor of a point can also be used with rational coefficients
     thanks to convertor of CoordNT. */
  _One_root_point_2(const CoordNT& x, const CoordNT& y) :
    Point_handle(Point_rep(x, y))
  {}

  /*! Get the x-coordinate. */
  const CoordNT& x() const { return (this->ptr()->_x); }

  /*! Get the y-coordinate. */
  const CoordNT& y() const { return (this->ptr()->_y); }

  /*! Check for equality. */
  bool equals(const Self& p) const
  {
    if (this->identical(p)) return (true);

    return (CGAL::compare(this->ptr()->_x, p.ptr()->_x) == EQUAL &&
            CGAL::compare(this->ptr()->_y, p.ptr()->_y) == EQUAL);
  }

  bool operator != (const Self& p) const { return ! equals(p); }

  bool operator == (const Self& p) const { return equals(p); }

  /*! Set the point coordinates. */
  void set(const NT& x, const NT& y)
  {
    this->copy_on_write();
    this->ptr()->_x = CoordNT(x);
    this->ptr()->_y = CoordNT(y);
    return;
  }

  /*! Set the point coordinates. */
  void set(const CoordNT& x, const CoordNT& y)
  {
    this->copy_on_write();
    this->ptr()->_x = x;
    this->ptr()->_y = y;
    return;
  }
};

/*!
 * Exporter for conic arcs.
 */
template <typename NT, bool Filter>
std::ostream& operator<<(std::ostream& os,
                         const _One_root_point_2<NT, Filter>& p)
{
  os << CGAL::to_double(p.x()) << ' ' << CGAL::to_double(p.y());
  return (os);
}

/*
template <class NT, bool Filter>
std::istream & operator >> (std::istream & is,
                            _One_root_point_2<NT, Filter>& p)
{
  typename _One_root_point_2<NT, Filter>::CoordNT ort1,ort2;
  is >> ort1 >> ort2;
  p=_One_root_point_2<NT, Filter>(ort1,ort2);
  return is;
}
*/

/*! \class
 * Representation of a circle, a circular arc or a line segment.
 */
template <typename Kernel_, bool Filter_>
class _Circle_segment_2 {
public:
  typedef Kernel_                                          Kernel;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT, Filter_>                   Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Segment_2                       Segment_2;
  typedef typename Kernel::Line_2                          Line_2;

protected:
  typedef typename Point_2::CoordNT                        CoordNT;

  // Data members:
  Line_2 m_line;              // The supporting line (for line segments).
  Circle_2 m_circ;            // The supporting circle (for circular arcs).
  bool m_is_full;             // Whether we have a full circle.
  bool m_has_radius;          // Is the radius (not just the squared radius)
                              // explicitly specified).
  NT m_radius;                // The radius, in case it is specified.
  Point_2 m_source;           // The source point.
  Point_2 m_target;           // The target point.
  Orientation m_orient;       // The orientation (COLLINEAR for line segments).

public:
  /*! Default constructor. */
  _Circle_segment_2() :
    m_is_full(false),
    m_has_radius(false),
    m_orient(COLLINEAR)
  {}

  /*! Constructor from a line segment.
   * \param seg The segment.
   */
  _Circle_segment_2(const Segment_2& seg) :
    m_line(seg),
    m_is_full(false),
    m_has_radius(false),
    m_source(seg.source().x(), seg.source().y()),
    m_target(seg.target().x(), seg.target().y()),
    m_orient(COLLINEAR)
  {}

  /*! Constructor from of a line segment.
   * \param ps The source point.
   * \param pt The target point.
   */
  _Circle_segment_2(const typename Kernel::Point_2& ps,
                    const typename Kernel::Point_2& pt) :
    m_line(ps, pt),
    m_is_full(false),
    m_has_radius(false),
    m_source(ps.x(), ps.y()),
    m_target(pt.x(), pt.y()),
    m_orient(COLLINEAR)
  {}

  /*! Constructor of a segment, given a supporting line and two endpoints,
   * which need not necessarily have rational coordinates.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting line.
   */
  _Circle_segment_2(const Line_2& line,
                    const Point_2& source, const Point_2& target) :
    m_line(line),
    m_is_full(false),
    m_has_radius(false),
    m_source(source),
    m_target(target),
    m_orient(COLLINEAR)
  {
    CGAL_precondition(CGAL::compare(source.x() * line.a() + line.c(),
                                    -source.y() * line.b()) == EQUAL);

    CGAL_precondition(CGAL::compare(target.x() * line.a() + line.c(),
                                    -target.y() * line.b()) == EQUAL);
  }

  /*! Constructor from a circle.
   * \param circ The circle.
   */
  _Circle_segment_2(const Circle_2& circ) :
    m_circ(circ),
    m_is_full(true),
    m_has_radius(false),
    m_orient(circ.orientation())
  { CGAL_assertion(m_orient != COLLINEAR); }

  /*! Constructor from a circle.
   * \param c The circle center.
   * \param r The radius.
   * \param orient The orientation of the circle.
   */
  _Circle_segment_2(const typename Kernel::Point_2& c, const NT& r,
                    Orientation orient = COUNTERCLOCKWISE) :
    m_circ(c, r*r, orient),
    m_is_full(true),
    m_has_radius(true),
    m_radius(r),
    m_orient(orient)
  { CGAL_assertion (orient != COLLINEAR); }

  /*! Constructor of a circular arc, given a supporting circle and two
   * endpoints, which need not necessarily have rational coordinates.
   * The orientation of the circle determines the orientation of the arc.
   * \param circ The supporting circle.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting circle.
   */
  _Circle_segment_2(const Circle_2& circ,
                    const Point_2& source, const Point_2& target) :
    m_circ(circ),
    m_is_full(false),
    m_has_radius(false),
    m_source(source),
    m_target(target),
    m_orient(circ.orientation())
  {
    CGAL_assertion(m_orient != COLLINEAR);

    CGAL_precondition
      (CGAL::compare(CGAL::square(source.x() - circ.center().x()),
                     circ.squared_radius() -
                     CGAL::square(source.y() - circ.center().y())) == EQUAL);

    CGAL_precondition
      (CGAL::compare(CGAL::square(target.x() - circ.center().x()),
                     circ.squared_radius() -
                     CGAL::square(target.y() - circ.center().y())) == EQUAL);
  }

  /*! Constructor of a circular arc, given a supporting circle and two
   * endpoints, which need not necessarily have rational coordinates.
   * \param c The circle center.
   * \param r The radius.
   * \param orient The orientation of the circle.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting circle.
   */
  _Circle_segment_2(const typename Kernel::Point_2& c,
                    const NT& r, Orientation orient,
                    const Point_2& source, const Point_2& target) :
    m_circ(c, r*r, orient),
    m_is_full(false),
    m_has_radius(true),
    m_radius(r),
    m_source(source),
    m_target(target),
    m_orient(orient)
  {
    CGAL_assertion(orient != COLLINEAR);

    CGAL_precondition
      (CGAL::compare(CGAL::square(source.x() - c.x()),
                     CGAL::square(r) -
                     CGAL::square(source.y() - c.y())) == EQUAL);

    CGAL_precondition
      (CGAL::compare(CGAL::square(target.x() - c.x()),
                     CGAL::square(r) -
                     CGAL::square(target.y() - c.y())) == EQUAL);
  }

  /*! Constructor of a circular arc, from the given three points, in case of
   * three collinear points, a segment will be constructed.
   * \param p1 The arc source.
   * \param p2 A point in the interior of the arc.
   * \param p3 The arc target.
   * \pre p1 and p3 are not equal.
   */
   _Circle_segment_2(const typename Kernel::Point_2& p1,
                     const typename Kernel::Point_2& p2,
                     const typename Kernel::Point_2& p3) :
     m_is_full(false),
     m_has_radius(false),
     m_source(p1.x(), p1.y()),
     m_target(p3.x(), p3.y())
  {
    // Set the source and target.
    NT x1 = p1.x();
    NT y1 = p1.y();
    NT x2 = p2.x();
    NT y2 = p2.y();
    NT x3 = p3.x();
    NT y3 = p3.y();

    // Make sure that the source and the target are not the same.
    CGAL_precondition(Kernel().compare_xy_2_object()(p1, p3) != EQUAL);

    // Compute the lines: A1*x + B1*y + C1 = 0,
    //               and: A2*x + B2*y + C2 = 0,
    // where:
    const NT _two  = 2;

    const NT A1 = _two*(x1 - x2);
    const NT B1 = _two*(y1 - y2);
    const NT C1 =
      CGAL::square(y2) - CGAL::square(y1) + CGAL::square(x2) - CGAL::square(x1);

    const NT A2 = _two*(x2 - x3);
    const NT B2 = _two*(y2 - y3);
    const NT C2 =
      CGAL::square(y3) - CGAL::square(y2) + CGAL::square(x3) - CGAL::square(x2);

    // Compute the coordinates of the intersection point between the
    // two lines, given by (Nx / D, Ny / D), where:
    const NT Nx = B1*C2 - B2*C1;
    const NT Ny = A2*C1 - A1*C2;
    const NT D = A1*B2 - A2*B1;

    // Make sure the three points are not collinear.
    const bool points_collinear = (CGAL::sign (D) == ZERO);

    if (points_collinear) {
      m_line  = Line_2(p1, p3);
      m_orient = COLLINEAR;
      return;
    }

    // The equation of the underlying circle is given by:

    NT x_center = Nx / D;
    NT y_center = Ny / D;

    typename Kernel::Point_2 circ_center(x_center, y_center);

    NT sqr_rad =
      (CGAL::square(D*x2 - Nx) + CGAL::square(D*y2 - Ny)) / CGAL::square(D);

    // Determine the orientation: If the mid-point forms a left-turn with
    // the source and the target points, the orientation is positive (going
    // counterclockwise).
    // Otherwise, it is negative (going clockwise).
    Kernel ker;
    typename Kernel::Orientation_2 orient_f = ker.orientation_2_object();

    if (orient_f(p1, p2, p3) == LEFT_TURN) m_orient = COUNTERCLOCKWISE;
    else m_orient = CLOCKWISE;
     m_circ = Circle_2(circ_center, sqr_rad, m_orient);
  }

  /*! Get the orientation of the curve.
   * \return COLLINEAR in case of a line segment,
   *         CLOCKWISE or COUNTERCLOCKWISE for circular curves.
   */
  inline Orientation orientation() const { return (m_orient); }

  /*! Check if the arc is linear. */
  inline bool is_linear() const { return (m_orient == COLLINEAR); }

  /*! Check if the arc is circular. */
  inline bool is_circular() const { return (m_orient != COLLINEAR); }

  /*! Obtain the supporting line.
   * \pre The curve orientation is COLLINEAR.
   */
  const Line_2& supporting_line() const
  {
    CGAL_precondition(m_orient == COLLINEAR);
    return m_line;
  }

  /*! Obtain the supporting circle.
   * \pre The curve orientation is not COLLINEAR.
   */
  const Circle_2& supporting_circle() const
  {
    CGAL_precondition(m_orient != COLLINEAR);
    return m_circ;
  }

  /*! Check if the curve is a full circle. */
  bool is_full() const { return (m_is_full); }

  /*!
   * Get the source point.
   * \pre The curve is not a full circle.
   */
  const Point_2& source() const
  {
    CGAL_precondition(! m_is_full);
    return (m_source);
  }

  /*!
   * Get the target point.
   * \pre The curve is not a full circle.
   */
  const Point_2& target() const
  {
    CGAL_precondition(! m_is_full);
    return (m_target);
  }

  /*! Get the vertical tangency points the arc contains.
   * \param vpts Output: The vertical tangency points.
   * \pre The curve is circular.
   * \return The number of points (0, 1, or 2).
   */
  unsigned int vertical_tangency_points(Point_2* vpts) const
  {
    CGAL_precondition(m_orient != COLLINEAR);
    unsigned int n_vpts = 0;

    if (m_is_full) {
      // In case of a full circle, create both vertical tangency points:
      const NT& x0 = m_circ.center().x();
      const NT& y0 = m_circ.center().y();
      CoordNT xv_left;
      CoordNT xv_right;

      if (m_has_radius) {
        // In case the radius is explicitly given:
        xv_left = CoordNT(x0 - m_radius);
        xv_right = CoordNT(x0 + m_radius);
      }
      else {
        // In case only the squared root is given:
        xv_left = CoordNT(x0, NT(-1), m_circ.squared_radius());
        xv_right = CoordNT(x0, NT(1), m_circ.squared_radius());
      }

      vpts[0] = Point_2(xv_left, y0);
      vpts[1] = Point_2(xv_right, y0);
      return (2);
    }

    if (m_orient == COUNTERCLOCKWISE) {
      // Compute the vertical tangency points for the arc:
      n_vpts = _ccw_vertical_tangency_points(m_source, m_target, vpts);
    }
    else {
      // Compute the vertical tangency points for the opposite arc:
      n_vpts = _ccw_vertical_tangency_points(m_target, m_source, vpts);

      // Swap their order, if necessary.
      if (n_vpts == 2) {
        Point_2 temp = vpts[0];
        vpts[0] = vpts[1];
        vpts[1] = temp;
      }
    }

    return (n_vpts);
  }

private:
  /*! Obtain the vertical tangency points the arc contains, assuming it is
   * counterclockwise oriented.
   * \param vpts Output: The vertical tangency points.
   * \return The number of points (0, 1, or 2).
   */
  unsigned int _ccw_vertical_tangency_points(const Point_2& src,
                                             const Point_2& trg,
                                             Point_2* vpts) const
  {
    unsigned int n_vpts = 0;
    const NT& x0 = m_circ.center().x();
    const NT& y0 = m_circ.center().y();
    int qs = _quart_index(src);
    int qt = _quart_index(trg);

    if (qs == qt) {
      if ((qs == 0 || qs == 1) && CGAL::compare(src.x(), trg.x()) == LARGER)
        // We have an x-monotone arc lying on the upper half of the circle:
        return (0);

      if ((qs == 2 || qs == 3) && CGAL::compare(src.x(), trg.x()) == SMALLER)
        // We have an x-monotone arc lying on the lower half of the circle:
        return (0);
    }

    // Make sure the target quarter is larger than the source quarter, by
    // adding 4 to its index, if necessary.
    if (qt <= qs) qt += 4;

    // Start traversing the quarter-planes and collect the vertical tangency
    // points we encounter.
    while (qs < qt) {
      if ((qs % 4) == 1) {
        // We collect the left tangency point when going from Q[1] to Q[2]:
        if (CGAL::compare(x0, trg.x()) != LARGER ||
            CGAL::compare(y0, trg.y()) != EQUAL)
        {
          if (m_has_radius)
            vpts[n_vpts] = Point_2(CoordNT(x0 - m_radius), y0);
          else
            vpts[n_vpts] =
              Point_2(CoordNT(x0, NT(-1), m_circ.squared_radius()), y0);

          n_vpts++;
        }
      }
      else if ((qs % 4) == 3) {
        // We collect the right tangency point when going from Q[3] to Q[0]:
        if (CGAL::compare(x0, trg.x()) != SMALLER ||
            CGAL::compare(y0, trg.y()) != EQUAL)
        {
          if (m_has_radius)
            vpts[n_vpts] = Point_2(CoordNT(x0 + m_radius), y0);
          else
            vpts[n_vpts] =
              Point_2(CoordNT(x0, NT(1), m_circ.squared_radius()), y0);
          n_vpts++;
        }
      }

      qs++;
    }

    return (n_vpts);
  }

  /*! Obtain the index of the quarter-plane containing the given point,
   * where the circle center is considered to be the origin.
   */
  int _quart_index(const Point_2& p) const
  {
    // The plane looks like:
    //
    //      Q[1] :  |   Q[0]:
    //      x <= 0  |   x >  0
    //      y >  0  |   y >= 0
    //    ----------+-----------
    //      Q[2] :  |   Q[3]:
    //      x <  0  |   x >= 0
    //      y <= 0  |   y <  0
    //
    const CGAL::Sign sign_x = CGAL::sign(p.x() - m_circ.center().x());
    const CGAL::Sign sign_y = CGAL::sign(p.y() - m_circ.center().y());

    if (sign_x == POSITIVE) return ((sign_y == NEGATIVE) ? 3 : 0);
    else if (sign_x == NEGATIVE) return ((sign_y == POSITIVE) ? 1 : 2);

    CGAL_assertion (sign_y != ZERO);
    return ((sign_y == POSITIVE) ? 1 : 3);
  }
};

/*! Exporter for line segments and circular arcs.
 */
template <typename Kernel, bool Filter>
std::ostream&
operator<<(std::ostream& os, const _Circle_segment_2<Kernel, Filter>& c)
{
  if (c.orientation() == COLLINEAR) {
    os<< "segment: " << c.source() << " -> " << c.target();
  }
  else {
    if (! c.is_full()) {
      os << "circular arc: " << c.supporting_circle() << ' '
         << c.source() << " -> " << c.target();
    }
    else {
      os << "circular arc: " << c.supporting_circle();
    }
  }

  return (os);
}

/*! \class
 * Representation of an x-monotone circular arc.
 */
template <typename Kernel_, bool Filter_>
class _X_monotone_circle_segment_2 {
public:
  typedef Kernel_                                          Kernel;
  typedef _X_monotone_circle_segment_2<Kernel, Filter_>    Self;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT, Filter_>                   Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Line_2                          Line_2;
  typedef typename Point_2::CoordNT                        CoordNT;

  // Type definition for the intersection points mapping.
  typedef std::pair<unsigned int, unsigned int>         Curve_id_pair;
  typedef unsigned int                                  Multiplicity;
  typedef std::pair<Point_2, Multiplicity>              Intersection_point;
  typedef std::list<Intersection_point>                 Intersection_list;

  /*!
   * \struct Less functor for Curve_id_pair.
   */
  struct Less_id_pair {
    bool operator()(const Curve_id_pair& ip1, const Curve_id_pair& ip2) const
    {
      // Compare the pairs of IDs lexicographically.
      return (ip1.first < ip2.first ||
              (ip1.first == ip2.first && ip1.second < ip2.second));
    }
  };

  typedef std::map<Curve_id_pair, Intersection_list, Less_id_pair>
                                                  Intersection_map;
  typedef typename Intersection_map::value_type   Intersection_map_entry;
  typedef typename Intersection_map::iterator     Intersection_map_iterator;

protected:
  NT m_first;           // The x-coordinate of the circle center.
                        // Or: the coefficient of x in the line equation.
  NT m_second;          // The y-coordinate of the circle center.
                        // Or: the coefficient of y in the line equation.
  NT m_third;           // The squared radius of the supporting circle.
                        // Or: the free coefficient in the line equation.
  Point_2 m_source;     // The source point.
  Point_2 m_target;     // The target point.

  enum {
    IS_DIRECTED_RIGHT_MASK = 1,
    IS_VERTICAL_SEGMENT_MASK = 2,
    COUNTERCLOCKWISE_CODE = 4,
    CLOCKWISE_CODE = 8,
    ORIENTATION_MASK = 4 + 8,
    INDEX_SHIFT_BITS = 4
  };

  unsigned int m_info;       // A bit vector, where:
                             // Bit 0 (the LSB): marks if the arc is directed
                             //                  from left to right.
                             // Bit 1: marks if the arc is a vertical segment.
                             // Bits 2-3: mark the orientation.
                             // The rest of the bits represent the curve index.

public:
  /*! Default constructor.
   */
  _X_monotone_circle_segment_2() :
    m_first(),
    m_second(),
    m_third(),
    m_source(),
    m_target(),
    m_info(0)
  {}

  /*! Construct an arc from a line segment.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   */
  _X_monotone_circle_segment_2(const Line_2& line,
                               const Point_2& source, const Point_2& target,
                               unsigned int index = 0) :
    m_first(line.a()),
    m_second(line.b()),
    m_third(line.c()),
    m_source(source),
    m_target(target),
    m_info(index << INDEX_SHIFT_BITS)
  {
    // Check if the segment is directed left or right:
    Comparison_result res = CGAL::compare(source.x(), target.x());

    if (res == EQUAL) {
      CGAL_precondition(CGAL::sign(m_second) == ZERO);

      // We have a vertical segment - compare the points by their
      // y-coordinates:
      m_info = (m_info | IS_VERTICAL_SEGMENT_MASK);
      res = CGAL::compare(source.y(), target.y());
    }

    CGAL_precondition(res != EQUAL);
    if (res == SMALLER) m_info = (m_info | IS_DIRECTED_RIGHT_MASK);
  }

  /*! Construct a segment arc from two kernel points
   * \param source the source point.
   * \ param target the target point.
   * \pre source and target are not equal.
   */
  _X_monotone_circle_segment_2(const typename Kernel::Point_2& source,
                               const typename Kernel::Point_2& target) :
    m_source(source.x(), source.y()),
    m_target(target.x(), target.y()),
    m_info(0)
  {
    Line_2 line(source, target);
    m_first  = line.a();
    m_second = line.b();
    m_third  = line.c();

    // Check if the segment is directed left or right:
    Comparison_result res = CGAL::compare(source.x(), target.x());

    if (res == EQUAL) {
      CGAL_precondition(CGAL::sign(m_second) == ZERO);

      // We have a vertical segment - compare the points by their
      // y-coordinates:
      m_info = (m_info | IS_VERTICAL_SEGMENT_MASK);
      res = CGAL::compare(source.y(), target.y());
    }

    CGAL_precondition(res != EQUAL);
    if (res == SMALLER) m_info = (m_info | IS_DIRECTED_RIGHT_MASK);
  }

  /*! Construct a circular arc.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \param orient The orientation of the arc.
   */
  _X_monotone_circle_segment_2(const Circle_2& circ,
                               const Point_2& source, const Point_2& target,
                               Orientation orient,
                               unsigned int index = 0) :
    m_first(circ.center().x()),
    m_second(circ.center().y()),
    m_third(circ.squared_radius()),
    m_source(source),
    m_target(target),
    m_info(index << INDEX_SHIFT_BITS)
  {
    // Check if the segment is directed left or right:
    Comparison_result res = CGAL::compare (source.x(), target.x());

    CGAL_precondition(res != EQUAL);
    if (res == SMALLER) m_info = (m_info | IS_DIRECTED_RIGHT_MASK);

    // Set the orientation.
    CGAL_precondition (orient != COLLINEAR);
    if (orient == COUNTERCLOCKWISE) m_info = (m_info | COUNTERCLOCKWISE_CODE);
    else m_info = (m_info | CLOCKWISE_CODE);
  }

  /*! Check if the arc is linear. */
  inline bool is_linear () const { return ((m_info & ORIENTATION_MASK) == 0); }

  /*! Check if the arc is circular. */
  inline bool is_circular () const
  { return ((m_info & ORIENTATION_MASK) != 0); }

  /*! Obtain the supporting line.
   * \pre The arc is linear (a line segment).
   */
  Line_2 supporting_line() const
  {
    CGAL_precondition (is_linear());
    return (Line_2 (a(), b(), c()));
  }

  /*! Obtain the supporting circle.
   * \pre The arc is circular.
   */
  Circle_2 supporting_circle() const
  {
    CGAL_precondition (is_circular());

    typename Kernel::Point_2 center(x0(), y0());
    return (Circle_2(center , sqr_r(), orientation()));
  }

  /*! Get the source point. */
  inline const Point_2& source() const { return (m_source); }

  /*! Get the target point. */
  inline const Point_2& target() const { return (m_target); }

  /*! True if the arc is directed right, false otherwise. */
  bool is_directed_right() const
  { return ((m_info & IS_DIRECTED_RIGHT_MASK) != 0); }

  bool has_left() const { return true; }

  bool has_right() const { return true; }

  /*! Get the left endpoint of the arc. */
  inline const Point_2& left() const
  { return (((m_info & IS_DIRECTED_RIGHT_MASK) != 0) ? m_source : m_target); }

  /*! Get the right endpoint of the arc. */
  inline const Point_2& right() const
  { return (((m_info & IS_DIRECTED_RIGHT_MASK) != 0) ? m_target : m_source); }

  /*! Check whether the given point is in the x-range of the arc.
   */
  bool is_in_x_range(const Point_2& p) const
  {
    Comparison_result res = CGAL::compare (p.x(), left().x());

    if (res == SMALLER) return false;
    else if (res == EQUAL) return true;

    return (CGAL::compare (p.x(), right().x()) != LARGER);
  }

  /*! Check if the arc is a vertical segment. */
  inline bool is_vertical() const
  { return ((m_info & IS_VERTICAL_SEGMENT_MASK) != 0); }

  /*! Get the orientation of the arc. */
  inline Orientation orientation() const
  {
    unsigned int or_ = (m_info & ORIENTATION_MASK);
    if (or_ == COUNTERCLOCKWISE_CODE) return (CGAL::COUNTERCLOCKWISE);
    else if (or_ == CLOCKWISE_CODE) return (CGAL::CLOCKWISE);

    CGAL_assertion(or_ == 0);
    return (CGAL::COLLINEAR);
  }

  /*! Check the position of a given point with respect to the arc.
   */
  Comparison_result point_position(const Point_2& p) const
  {
    if (is_linear()) return (_line_point_position(p));
    else return (_circ_point_position (p));
  }

  /*! Compare the two arcs to the right of their intersection point.
   */
  Comparison_result compare_to_right(const Self& cv, const Point_2& p) const
  {
    if (is_linear()) {
      if (cv.is_linear()) return (_lines_compare_to_right (cv, p));
      Comparison_result res = cv._circ_line_compare_to_right (*this, p);
      if (res != EQUAL) res = (res == SMALLER) ? LARGER : SMALLER;
      return (res);
    }
    else {
      if (cv.is_linear()) return (_circ_line_compare_to_right (cv, p));
      return (_circs_compare_to_right (cv, p));
    }
  }

  /*! Compare the two arcs to the left of their intersection point.
   */
  Comparison_result compare_to_left(const Self& cv, const Point_2& p) const
  {
    if (is_linear()) {
      if (cv.is_linear()) return (_lines_compare_to_left (cv, p));
      Comparison_result res = cv._circ_line_compare_to_left(*this, p);
      if (res != EQUAL) res = (res == SMALLER) ? LARGER : SMALLER;
      return (res);
    }
    else {
      if (cv.is_linear()) return (_circ_line_compare_to_left(cv, p));
      return (_circs_compare_to_left (cv, p));
    }
  }

  /*! Check whether the two arcs have the same supporting curve.
   */
  bool has_same_supporting_curve(const Self& cv) const
  {
    // Check if the curve indices are the same.
    if (_index() != 0 && _index() == cv._index()) return true;

    // Make sure that the supporting curves are of the same type.
    if (is_linear() && ! cv.is_linear()) return false;

    if (! is_linear() && cv.is_linear()) return false;

    // Compare the curve coefficients.
    if (! is_linear()) {
      // The two circles must have the same center and the same radius.
      return (CGAL::compare(x0(), cv.x0()) == EQUAL &&
              CGAL::compare(y0(), cv.y0()) == EQUAL &&
              CGAL::compare(sqr_r(), cv.sqr_r()) == EQUAL);
    }

    // Compare the line equations: Note that these may be scaled.
    NT fact1;
    NT fact2;

    if (is_vertical()) {
      if (! cv.is_vertical()) return false;

      fact1 = a();
      fact2 = cv.a();
    }
    else {
      fact1 = b();
      fact2 = cv.b();
    }

    return (CGAL::compare(fact2*a(), fact1*cv.a()) == EQUAL &&
            CGAL::compare(fact2*b(), fact1*cv.b()) == EQUAL &&
            CGAL::compare(fact2*c(), fact1*cv.c()) == EQUAL);
  }

  /*! Check whether the two curves are equal.
   */
  bool equals(const Self& cv) const
  {
    if (! this->has_same_supporting_curve(cv)) return false;

    if (is_linear()) {
      // In case of line segments we can swap the source and target:
      return ((m_source.equals(cv.m_source) && m_target.equals(cv.m_target)) ||
              (m_source.equals(cv.m_target) && m_target.equals(cv.m_source)));
    }

    // Once again, opposite circular arcs are considered to be equal:
    return ((orientation() == cv.orientation() &&
             m_source.equals(cv.m_source) && m_target.equals(cv.m_target)) ||
            (orientation() != cv.orientation() &&
             m_source.equals(cv.m_target) && m_target.equals(cv.m_source)));
  }

  /*! Split the curve at a given point into two sub-arcs.
   */
  void split(const Point_2& p, Self& c1, Self& c2) const
  {
    // Copy the properties of this arc to the sub-arcs.
    c1 = *this;
    c2 = *this;

    // Change the endpoint, such that c1 lies to the right of c2:
    if (is_directed_right()) {
      c1.m_target = p;
      c2.m_source = p;
    }
    else {
      c1.m_source = p;
      c2.m_target = p;
    }
  }

  /*! Compute the intersections between the two arcs or segments.
   */
  template <typename OutputIterator>
  OutputIterator intersect(const Self& cv, OutputIterator oi,
                           Intersection_map* inter_map = nullptr) const
  {
    typedef std::pair<Point_2, Multiplicity>            Intersection_point;
    typedef boost::variant<Intersection_point, Self>    Intersection_result;

    // First check whether the two arcs have the same supporting curve.
    if (has_same_supporting_curve(cv)) {
      // Check for overlaps between the two arcs.
      Self overlap;

      if (_compute_overlap(cv, overlap)) {
        // There can be just a single overlap between two x-monotone arcs:
        *oi++ = Intersection_result(overlap);
        return oi;
      }

      // In case there is not overlap and the supporting curves are the same,
      // there cannot be any intersection points, unless the two arcs share
      // a common end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      Multiplicity mult = 0;
      if (left().equals(cv.left()) || left().equals(cv.right())) {
        *oi++ = Intersection_result(std::make_pair(left(), mult));
      }

      if (right().equals(cv.right()) || right().equals(cv.left())) {
        *oi++ = Intersection_result(std::make_pair(right(), mult));
      }

      return oi;
    }

    // Before computing the intersection points between the two supporting
    // curves, check if their intersection has already been computed and
    // cached.
    Curve_id_pair id_pair;
    Intersection_map_iterator map_iter;
    Intersection_list inter_list;
    bool invalid_ids = false;

    if (inter_map != nullptr && _index() != 0 && cv._index() != 0) {
      if (_index() < cv._index())
        id_pair = Curve_id_pair (_index(), cv._index());
      else id_pair = Curve_id_pair (cv._index(), _index());

      map_iter = inter_map->find(id_pair);
    }
    else {
      // In case one of the IDs is invalid, we do not look in the map neither
      // we cache the results.
      if (inter_map != nullptr) map_iter = inter_map->end();
      invalid_ids = true;
    }

    if ((inter_map == nullptr) || (map_iter == inter_map->end())) {
      // Compute the intersections points between the two supporting curves.
      if (is_linear()) {
        if (cv.is_linear()) _lines_intersect(cv, inter_list);
        else cv._circ_line_intersect(*this, inter_list);
      }
      else {
        if (cv.is_linear()) _circ_line_intersect(cv, inter_list);
        else _circs_intersect(cv, inter_list);
      }

      // Cache the result.
      if (! invalid_ids) (*inter_map)[id_pair] = inter_list;
    }
    else {
      // Obtain the precomputed intersection points from the map.
      inter_list = (*map_iter).second;
    }

    // Report only the intersection points that lie on both arcs.
    for (auto iter = inter_list.begin(); iter != inter_list.end(); ++iter) {
      if (this->_is_between_endpoints (iter->first) &&
          cv._is_between_endpoints (iter->first))
      {
        *oi++ = Intersection_result(*iter);
      }
    }

    return oi;
  }

  /*! Check whether it is possible to merge our arc with the given arc.
   */
  bool can_merge_with(const Self& cv) const
  {
    // In order to merge the two arcs, they should have the same supporting
    // curve.
    if (! this->has_same_supporting_curve(cv)) return false;

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    return (right().equals(cv.left()) || left().equals(cv.right()));
  }

  /*! Merge our arc with the given arc.
   * \pre The two arcs are mergeable.
   */
  void merge(const Self& cv)
  {
    CGAL_precondition(this->can_merge_with (cv));

    // Check if we should extend the arc to the left or to the right.
    if (right().equals(cv.left())) {
      // Extend the arc to the right.
      if (is_directed_right()) this->m_target = cv.right();
      else this->m_source = cv.right();
    }
    else {
      CGAL_precondition(left().equals(cv.right()));

      // Extend the arc to the left.
      if (is_directed_right()) this->m_source = cv.left();
      else this->m_target = cv.left();
    }
  }

  /*! construct an opposite arc. */
  Self construct_opposite() const
  {
    Self opp_cv;
    opp_cv.m_first = this->m_first;
    opp_cv.m_second = this->m_second;
    opp_cv.m_third = this->m_third;
    opp_cv.m_source = this->m_target;
    opp_cv.m_target = this->m_source;

    // Take care of the information bits: We flip the orientation bits and
    // the bits that marks the direction.
    if (is_linear()) opp_cv.m_info = (this->m_info ^ IS_DIRECTED_RIGHT_MASK);
    else
      opp_cv.m_info = (this->m_info ^ IS_DIRECTED_RIGHT_MASK ^ ORIENTATION_MASK);

    return (opp_cv);
  }

  Bbox_2 bbox() const
  {
    double x_min = to_double(left().x());
    double x_max = to_double(right().x());
    double y_min = to_double(left().y());
    double y_max = to_double(right().y());
    if (y_min > y_max) std::swap(y_min, y_max);
    if (is_circular()) {
      const Circle_2& circ = this->supporting_circle();
      if (_is_upper()) {
        y_max = to_double(circ.center().y())+
          std::sqrt(to_double(circ.squared_radius()));
      }
      else {
        y_min = to_double(circ.center().y()) -
          std::sqrt(to_double(circ.squared_radius()));
      }
    }

    return Bbox_2(x_min, y_min, x_max, y_max);
  }

protected:
  /*! Get the curve index. */
  inline unsigned int _index() const { return (m_info >> INDEX_SHIFT_BITS); }

  /// \name Accessors for circular arcs.
  //@{

  /*! Get the x-coordinate of the center of the supporting circle. */
  inline const NT& x0() const { return (m_first); }

  /*! Get the y-coordinate of the center of the supporting circle. */
  inline const NT& y0() const { return (m_second); }

  /*! Get the squared radius of the supporting circle. */
  inline const NT& sqr_r() const { return (m_third); }

  /*!
   * Check if the circular arc lies on the upper half of the supporting circle.
   */
  inline bool _is_upper() const
  {
    Orientation orient = orientation();
    bool dir_right = ((m_info & IS_DIRECTED_RIGHT_MASK) != 0);

    CGAL_precondition(orient != COLLINEAR);

    return ((orient == COUNTERCLOCKWISE && !dir_right) ||
            (orient == CLOCKWISE && dir_right));
  }
  //@}

  /// \name Accessors for line segments.
  //@{

  /*! Get the coefficient of x in the equation of the supporting line. */
  inline const NT& a() const { return (m_first); }

  /*! Get the coefficient of y in the equation of the supporting line. */
  inline const NT& b() const { return (m_second); }

  /*! Get the free coefficient in the equation of the supporting line. */
  inline const NT& c() const { return (m_third); }
  //@}

  /// \name Auxiliary functions for the point_position predicate.
  //@{

  /*!
   * Check the position of a given point with respect to a line segment.
   */
  Comparison_result _line_point_position(const Point_2& p) const
  {
    // Check if we have a vertical segment.

    CGAL_precondition(is_in_x_range(p));

    Comparison_result res;

    if (is_vertical()) {
      // left() is the lower endpoint:
      res = CGAL::compare(p.y(), left().y());

      if (res != LARGER) return (res);

      // left() is the upper endpoint:
      res = CGAL::compare (p.y(), right().y());

      if (res != SMALLER) return (res);

      // p lies in the interior of the vertical segment:
      return (EQUAL);
    }

    // Compute the y-coordinate of the vertical projection of p onto the
    // supporting line.
    const CoordNT        y_proj = (a()*p.x() + c()) / (-b());

    return (CGAL::compare (p.y(), y_proj));
  }

  /*!
   * Check the position of a given point with respect to a circular arc.
   */
  Comparison_result _circ_point_position(const Point_2& p) const
  {

    Comparison_result c_res = CGAL::compare (p.y(), y0());

    if (_is_upper()) {
      // Check if p lies below the "equator" (while the arc lies above it):
      if (c_res == SMALLER)
        return (SMALLER);
    }
    else {
      // Check if p lies above the "equator" (while the arc lies below it):
      if (c_res == LARGER)
        return (LARGER);
    }

    // Check if p lies inside the supporting circle, namely we have to check
    // whether (p.x() - x0)^2 + (p.y() - y0)^2 < r^2:
    Comparison_result res =
      CGAL::compare(CGAL::square (p.x() - x0()),
                    sqr_r() - CGAL::square (p.y() - y0()));

    // p lies on the circle:
    if (res == EQUAL) return (EQUAL);

    if (_is_upper()) {
      // If p is inside the circle, it lies below the upper arc:
      return (res);
    }
    else {
      // If p is inside the circle, it lies above the lower arc:
      return (res == SMALLER ? LARGER : SMALLER);
    }
  }
  //@}

  /// \name Auxiliary functions for the compare_to_right predicate.
  //@{

  /*!
   * Compare two line segments to the right of their intersection point.
   */
  Comparison_result _lines_compare_to_right(const Self& cv,
                                            const Point_2& /* p */) const
  {
    if (_index() != 0 && _index() == cv._index()) return (EQUAL);

    // Special treatment for vertical segments: a vertical segment is larger
    // than any other non-vertical segment.
    if (is_vertical()) {
      if (cv.is_vertical()) return (EQUAL);
      return (LARGER);
    }
    else if (cv.is_vertical()) return (SMALLER);

    // Compare the slopes: -A1/B1 and -A2/B2. We actually negate the slopes
    // and swap the result.
    return (CGAL::compare (cv.a()/cv.b(), a()/b()));
  }

  /*!
   * Compare a circular arcs (this) and a line segment (cv) to the right of
   * their intersection point.
   */
  Comparison_result _circ_line_compare_to_right(const Self& cv,
                                                const Point_2& p) const
  {
    // A vertical segment lies above any other circle to the right of p:
    if (cv.is_vertical()) return (SMALLER);

    // We have to compare the slopes of the supporting circles and the
    // supporting line at p:
    //
    //    p.x() - x0(1)           A(2)
    //   ---------------  and  - ------
    //    y0(1) - p.y()           B(2)
    //
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());

    // Check the case of a vertical tangent.
    if (sign_denom1 == ZERO) {
      // The arc lies above any line segment if it is an upper arc, or below
      // any segment if it is a lower arc.
      return (_is_upper() ? LARGER : SMALLER);
    }

    // Compare (p.x() - x0(1)) and (A(2)/B(2)*(p.y() - y0(1)).
    // Note that if the denominator is negative, we have to swap the result.
    const bool swap_res = (sign_denom1 == NEGATIVE);
    Comparison_result slope_res =
      CGAL::compare (p.x() - x0(), (p.y() - y0())*cv.a()/cv.b());

    if (slope_res != EQUAL) {
      // Swap the comparison result, if necessary:
      if (swap_res) slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      return (slope_res);
    }

    // In this case we have a tangency point at p. If the circular arc is an
    // upper arc, it must lie below the tangent line, and if it is a lower arc
    // it must lie above the tangent line.
    return (_is_upper() ? SMALLER : LARGER);
  }

  /*! Compare two circular arcs to the right of their intersection point.
   */
  Comparison_result _circs_compare_to_right(const Self& cv,
                                            const Point_2& p) const
  {
    if (_index() != 0 && _index() == cv._index()) {
      // Check the case of comparing two circular arcs that originate from the
      // same supporting circle. Their comparison result is not EQUAL only if
      // one is an upper arc and the other is a lower arc.
      if (_is_upper() && ! cv._is_upper()) return (LARGER);
      else if (! _is_upper() && cv._is_upper()) return (SMALLER);
      else return (EQUAL);
    }

    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    const CGAL::Sign  sign_numer1 = CGAL::sign(p.x() - x0());
    const CGAL::Sign  sign_denom1 = CGAL::sign(y0() - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign(p.x() - cv.x0());
    const CGAL::Sign  sign_denom2 = CGAL::sign(cv.y0() - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO) {
      if (sign_denom2 == ZERO) {
        if (_is_upper()) {
          if (cv._is_upper()) {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (sqr_r(), cv.sqr_r()));
          }
          else {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else {
          if (cv._is_upper()) {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv.sqr_r(), sqr_r()));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper() ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO) {
      return (cv._is_upper() ? SMALLER : LARGER);
    }

    // Try to act according to the slope signs.
    CGAL::Sign sign_slope1;
    CGAL::Sign sign_slope2;

    if (sign_numer1 == sign_denom1) sign_slope1 = POSITIVE;
    else if (sign_numer1 == ZERO) sign_slope1 = ZERO;
    else sign_slope1 = NEGATIVE;

    if (sign_numer2 == sign_denom2) sign_slope2 = POSITIVE;
    else if (sign_numer2 == ZERO) sign_slope2 = ZERO;
    else sign_slope2 = NEGATIVE;

    if ((sign_slope1 == POSITIVE && sign_slope2 != POSITIVE) ||
        (sign_slope1 == ZERO && sign_slope2 == NEGATIVE))
      return (LARGER);

    if ((sign_slope2 == POSITIVE && sign_slope1 != POSITIVE) ||
        (sign_slope2 == ZERO && sign_slope1 == NEGATIVE))
      return (SMALLER);

    // Compare the slopes of the two tangents to the circles.
    Comparison_result  slope_res;

    if (sign_slope1 == ZERO && sign_slope2 == ZERO)
    {
      // Special case were both circles have a horizontal tangent:
      slope_res = EQUAL;
    }
    else
    {
      // Actually compare the slopes.
      const bool swap_res = (sign_denom1 != sign_denom2);
      const CoordNT A = (cv.y0() - y0())*p.x() + (y0()*cv.x0() - cv.y0()*x0());
      const CoordNT B = (cv.x0() - x0())*p.y();

      slope_res = CGAL::compare (A, B);

      if (slope_res != EQUAL && swap_res) {
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      }
    }

    // In case the two circles have different tangent slopes at p:
    if (slope_res != EQUAL) return (slope_res);

    // In this case we have a tangency point at p.
    if (_is_upper()) {
      if (cv._is_upper()) {
        // The circle with a larger radius is above the other.
        return (CGAL::compare(sqr_r(), cv.sqr_r()));
      }
      else {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else {
      if (cv._is_upper()) {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare (cv.sqr_r(), sqr_r()));
      }
    }
  }
  //@}

  /// \name Auxiliary functions for the compare_to_left predicate.
  //@{

  /*!
   * Compare two line segments to the left of their intersection point.
   */
  Comparison_result _lines_compare_to_left(const Self& cv,
                                           const Point_2& ) const
  {
    if (_index() != 0 && _index() == cv._index()) return (EQUAL);

    // Special treatment for vertical segments: a vertical segment is smaller
    // than any other non-vertical segment.
    if (is_vertical()) {
      if (cv.is_vertical()) return (EQUAL);
      return (SMALLER);
    }
    else if (cv.is_vertical()) {
      return (LARGER);
    }

    // Compare the slopes: -A1/B1 and -A2/B2 and swap the result.
    //  We actually negate the slopes and compare them.
    return (CGAL::compare (a()/b(), cv.a()/cv.b()));
  }

  /*! Compare a circular arcs (this) and a line segment (cv) to the left of
   * their intersection point.
   */
  Comparison_result _circ_line_compare_to_left(const Self& cv,
                                               const Point_2& p) const
  {
    // A vertical segment lies below any other circle to the left of p:
    if (cv.is_vertical()) return (LARGER);

    // We have to compare the slopes of the supporting circles and the
    // supporting line at p, and return the swapped result:
    //
    //    p.x() - x0(1)           A(2)
    //   ---------------  and  - ------
    //    y0(1) - p.y()           B(2)
    //
    const CGAL::Sign  sign_denom1 = CGAL::sign (y0() - p.y());

    // Check the case of a vertical tangent.
    if (sign_denom1 == ZERO) {
      // The arc lies above any line segment if it is an upper arc, or below
      // any segment if it is a lower arc.
      return (_is_upper() ? LARGER : SMALLER);
    }

    // Compare (p.x() - x0(1)) and (A(2)/B(2)*(p.y() - y0(1)).
    // Note that if the denominator is negative, we have to swap the result.
    const bool swap_res = (sign_denom1 == NEGATIVE);
    Comparison_result slope_res =
      CGAL::compare(p.x() - x0(), (p.y() - y0()) * cv.a() / cv.b());

    if (slope_res != EQUAL) {
      if (swap_res)
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;

      // Swap at any case to get the position to the left:
      return ((slope_res == SMALLER) ? LARGER : SMALLER);
    }

    // In this case we have a tangency point at p. If the circular arc is an
    // upper arc, it must lie below the tangent line, and if it is a lower arc
    // it must lie above the tangent line.
    return (_is_upper() ? SMALLER : LARGER);
  }

  /*!
   * Compare the two arcs to the left of their intersection point.
   */
  Comparison_result _circs_compare_to_left(const Self& cv,
                                           const Point_2& p) const
  {
    if (_index() != 0 && _index() == cv._index()) {
      // Check the case of comparing two circular arcs that originate from the
      // same supporting circle. Their comparison result is not EQUAL only if
      // one is an upper arc and the other is a lower arc.
      if (_is_upper() && ! cv._is_upper()) return (LARGER);
      else if (! _is_upper() && cv._is_upper()) return (SMALLER);
      else return (EQUAL);
    }

    // We have to compare the slopes of the two supporting circles at p:
    //
    //    p.x() - x0(1)         p.x() - x0(2)
    //   ---------------  and  ---------------
    //    y0(1) - p.y()         y0(2) - p.y()
    //
    // Eventually, we should take the opposite result.
    const CGAL::Sign  sign_numer1 = CGAL::sign(p.x() - x0());
    const CGAL::Sign  sign_denom1 = CGAL::sign(y0() - p.y());
    const CGAL::Sign  sign_numer2 = CGAL::sign(p.x() - cv.x0());
    const CGAL::Sign  sign_denom2 = CGAL::sign(cv.y0() - p.y());

    // Check the case of vertical tangents.
    if (sign_denom1 == ZERO) {
      if (sign_denom2 == ZERO) {
        if (_is_upper()) {
          if (cv._is_upper()) {
            // The two circles have a vertical tangent:
            // The one with a larger radius is above the other.
            return (CGAL::compare (sqr_r(), cv.sqr_r()));
          }
          else {
            // The other curve is directed downwards:
            return (LARGER);
          }
        }
        else {
          if (cv._is_upper()) {
            // The other curve is directed upwards:
            return (SMALLER);
          }
          else {
            // The two circles have a vertical tangent:
            // The one with a smaller radius is above the other.
            return (CGAL::compare (cv.sqr_r(), sqr_r()));
          }
        }
      }

      // The other arc does not have a vertical tangent.
      return (_is_upper() ? LARGER : SMALLER);
    }
    else if (sign_denom2 == ZERO) {
      return (cv._is_upper() ? SMALLER : LARGER);
    }

    // Try to act according to the slope signs.
    CGAL::Sign sign_slope1;
    CGAL::Sign sign_slope2;

    if (sign_numer1 == sign_denom1) sign_slope1 = POSITIVE;
    else if (sign_numer1 == ZERO) sign_slope1 = ZERO;
    else sign_slope1 = NEGATIVE;

    if (sign_numer2 == sign_denom2) sign_slope2 = POSITIVE;
    else if (sign_numer2 == ZERO) sign_slope2 = ZERO;
    else sign_slope2 = NEGATIVE;

    if ((sign_slope1 == POSITIVE && sign_slope2 != POSITIVE) ||
        (sign_slope1 == ZERO && sign_slope2 == NEGATIVE))
      return (SMALLER);

    if ((sign_slope2 == POSITIVE && sign_slope1 != POSITIVE) ||
        (sign_slope2 == ZERO && sign_slope1 == NEGATIVE))
      return (LARGER);

    // Compare the slopes of the two tangents to the circles.
    Comparison_result  slope_res;

    if (sign_slope1 == ZERO && sign_slope2 == ZERO)
    {
      // Special case were both circles have a horizontal tangent:
      slope_res = EQUAL;
    }
    else {
      // Actually compare the slopes.
      const bool swap_res = (sign_denom1 != sign_denom2);
      const CoordNT A = (cv.y0() - y0())*p.x() + (y0()*cv.x0() - cv.y0()*x0());
      const CoordNT B = (cv.x0() - x0())*p.y();

      slope_res = CGAL::compare(A, B);

      if (slope_res != EQUAL && swap_res) {
        // Swap the comparison result, if necessary:
        slope_res = (slope_res == SMALLER) ? LARGER : SMALLER;
      }
    }

    // In case the two circles have different tangent slopes at p, return
    // the opposite of the slope result (since the slope result is the
    // comparison result to the right of the intersection point):
    if (slope_res != EQUAL) return ((slope_res == SMALLER) ? LARGER : SMALLER);

    // In this case we have a tangency point at p.
    if (_is_upper()) {
      if (cv._is_upper()) {
        // The circle with a larger radius is above the other.
        return (CGAL::compare(sqr_r(), cv.sqr_r()));
      }
      else {
        // The other curve is above our curve:
        return (SMALLER);
      }
    }
    else {
      if (cv._is_upper()) {
        // Out curve is above the other curve:
        return (LARGER);
      }
      else {
        // The circle with a smaller radius is above the other.
        return (CGAL::compare(cv.sqr_r(), sqr_r()));
      }
    }
  }
  //@}

  /// \name Auxiliary functions for computing intersections.
  //@{

  /*! Compute the intersections between two line segments.
   */
  void _lines_intersect(const Self& cv,
                        Intersection_list& inter_list) const
  {
    // The intersection of the lines:
    //   a1*x + b1*y + c1 = 0   and   a2*x + b2*y + c2 = 0 ,
    // is given by:
    //
    //      b1*c2 - c1*b2     c1*a2 - a1*c2
    //   ( --------------- , --------------- )
    //      a1*b2 - b1*a2     a1*b2 - b1*a2
    //
    unsigned int mult = 1;
    const NT denom = a()*cv.b() - b()*cv.a();

    // Make sure the supporting lines are not parallel.
    if (CGAL::sign(denom) == ZERO) return;

    const NT x = (b()*cv.c() - c()*cv.b()) / denom;
    const NT y = (c()*cv.a() - a()*cv.c()) / denom;
    Point_2 p (x, y);

    inter_list.push_back(Intersection_point(p, mult));
  }

  /*! Compute the intersections between the supporting circle of (*this) and
   * the supporting line of the segement cv.
   */
  void _circ_line_intersect(const Self& cv,
                            Intersection_list& inter_list) const
  {
    Point_2 p;
    unsigned int mult;

    // First check the special cases of vertical and horizontal lines.
    if (cv.is_vertical()) {
      // The equation of the vertical line is x = -c / a.
      // The y-coordinates of the intersection points are:
      //   y =  y0 +/- sqrt(r^2 - (x - x0)^2)
      //
      const NT vx = -cv.c() / cv.a();
      const NT vdisc = sqr_r() - CGAL::square (vx - x0());
      CGAL::Sign sign_vdisc = CGAL::sign (vdisc);

      if (sign_vdisc == NEGATIVE) {
        // The circle and the vertical line do not intersect.
        return;
      }
      else if (sign_vdisc == ZERO) {
        // A single tangency point, given by:
        mult = 2;
        p = Point_2 (vx, y0());
        inter_list.push_back(Intersection_point(p, mult));

        return;
      }

      // Compute the two intersection points:
      mult = 1;

      p = Point_2(CoordNT (vx), CoordNT (y0(), NT(-1), vdisc));
      inter_list.push_back (Intersection_point(p, mult));

      p = Point_2(CoordNT (vx), CoordNT (y0(), NT(1), vdisc));
      inter_list.push_back(Intersection_point(p, mult));

      return;
    }
    else if (CGAL::sign (cv.a()) == ZERO) {
      // The equation of the horizontal line is y = -c / b.
      // The y-coordinates of the intersection points are:
      //   x =  x0 +/- sqrt(r^2 - (y - y0)^2)
      //
      const NT hy = -cv.c() / cv.b();
      const NT hdisc = sqr_r() - CGAL::square (hy - y0());
      CGAL::Sign sign_hdisc = CGAL::sign (hdisc);

      // The circle and the vertical line do not intersect.
      if (sign_hdisc == NEGATIVE) return;
      else if (sign_hdisc == ZERO) {
        // A single tangency point, given by:
        mult = 2;
        p = Point_2(x0(), hy);
        inter_list.push_back(Intersection_point (p, mult));

        return;
      }

      // Compute the two intersection points:
      mult = 1;

      p = Point_2(CoordNT(x0(), NT(-1), hdisc), CoordNT (hy));
      inter_list.push_back(Intersection_point (p, mult));

      p = Point_2(CoordNT(x0(), NT(1), hdisc), CoordNT(hy));
      inter_list.push_back(Intersection_point (p, mult));

      return;
    }

    // Compute the squared distance between the line and the circle center,
    // inducing the discriminant of the quadratic equations we have to solve.
    const NT line_factor = CGAL::square(cv.a()) + CGAL::square(cv.b());
    const NT disc =
      line_factor*sqr_r() - CGAL::square(cv.a()*x0() + cv.b()*y0() + cv.c());
    CGAL::Sign sign_disc = CGAL::sign (disc);

    // The circle and the line do not intersect:
    if (sign_disc == NEGATIVE) return;

    // Compare the square-free part of the solution:
    const NT aux = cv.b()*x0() - cv.a()*y0();
    const NT x_base = (aux*cv.b() - cv.a()*cv.c()) / line_factor;
    const NT y_base = (-aux*cv.a() - cv.b()*cv.c()) / line_factor;

    if (sign_disc == ZERO) {
      // A single tangency point, given by:
      mult = 2;
      p = Point_2(x_base, y_base);
      inter_list.push_back(Intersection_point(p, mult));

      return;
    }

    // We have two intersection points, whose coordinates are one-root numbers.
    bool minus_root_first = (CGAL::sign(cv.b()) == POSITIVE);
    const NT x_root_coeff = cv.b() / line_factor;
    const NT y_root_coeff = cv.a() / line_factor;

    mult = 1;
    if (minus_root_first) {
      p = Point_2(CoordNT(x_base, -x_root_coeff, disc),
                  CoordNT(y_base, y_root_coeff, disc));
      inter_list.push_back(Intersection_point (p, mult));

      p = Point_2(CoordNT(x_base, x_root_coeff, disc),
                  CoordNT(y_base, -y_root_coeff, disc));
      inter_list.push_back(Intersection_point(p, mult));
    }
    else {
      p = Point_2(CoordNT(x_base, x_root_coeff, disc),
                  CoordNT(y_base, -y_root_coeff, disc));
      inter_list.push_back(Intersection_point(p, mult));

      p = Point_2(CoordNT(x_base, -x_root_coeff, disc),
                  CoordNT(y_base, y_root_coeff, disc));
      inter_list.push_back(Intersection_point(p, mult));
    }
  }

  /*! Compute the intersections between two circles.
   */
  void _circs_intersect(const Self& cv, Intersection_list& inter_list) const
  {
    Point_2 p;
    unsigned int mult;

    // Compute the squared distance between the circle centers, inducing the
    // discriminant of the quadratic equations we have to solve.
    const NT diff_x = cv.x0() - x0();
    const NT diff_y = cv.y0() - y0();
    const NT sqr_dist = CGAL::square(diff_x) + CGAL::square(diff_y);
    const NT diff_sqr_rad = sqr_r() - cv.sqr_r();
    const NT disc = 2 * sqr_dist * (sqr_r() + cv.sqr_r()) -
      (CGAL::square(diff_sqr_rad) + CGAL::square(sqr_dist));
    CGAL::Sign sign_disc = CGAL::sign (disc);

    // The two circles do not intersect.
    if (sign_disc == NEGATIVE) return;

    // Compare the square-free part of the solution:
    const NT x_base = ((x0() + cv.x0()) + diff_x * diff_sqr_rad / sqr_dist) / 2;
    const NT y_base = ((y0() + cv.y0()) + diff_y * diff_sqr_rad / sqr_dist) / 2;

    if (sign_disc == ZERO) {
      // A single tangency point, given by:
      mult = 2;
      p = Point_2(x_base, y_base);
      inter_list.push_back(Intersection_point(p, mult));
      return;
    }

    // We have two intersection points, whose coordinates are one-root numbers.
    CGAL::Sign sign_diff_y = CGAL::sign (diff_y);
    bool minus_root_first;

    if (sign_diff_y == ZERO)
      minus_root_first = (CGAL::sign (diff_x) == NEGATIVE);
    else
      minus_root_first = (sign_diff_y == POSITIVE);

    const NT x_root_coeff = diff_y / (2 * sqr_dist);
    const NT y_root_coeff = diff_x / (2 * sqr_dist);

    mult = 1;
    if (minus_root_first) {
      p = Point_2(CoordNT(x_base, -x_root_coeff, disc),
                  CoordNT(y_base, y_root_coeff, disc));
      inter_list.push_back(Intersection_point (p, mult));

      p = Point_2(CoordNT(x_base, x_root_coeff, disc),
                  CoordNT(y_base, -y_root_coeff, disc));
      inter_list.push_back(Intersection_point (p, mult));
    }
    else {
      p = Point_2(CoordNT(x_base, x_root_coeff, disc),
                  CoordNT(y_base, -y_root_coeff, disc));
      inter_list.push_back(Intersection_point (p, mult));

      p = Point_2(CoordNT(x_base, -x_root_coeff, disc),
                  CoordNT(y_base, y_root_coeff, disc));
      inter_list.push_back(Intersection_point(p, mult));
    }
  }

  /*! Check if the given point lies on the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_between_endpoints(const Point_2& p) const
  {
    if (is_linear()) {
      if (is_vertical()) {
        // Check if the point is in the y-range of the arc.
        // Note that left() is the lower endpoint and right() is the upper
        // endpoint of the segment in this case.
        Comparison_result res = CGAL::compare(p.y(), left().y());

        if (res == SMALLER) return false;
        else if (res == EQUAL) return true;

        return (CGAL::compare(p.y(), right().y()) != LARGER);
      }

      // For non-vertical segments, it is sufficient to check if the point
      // is in the x-range of the arc.
      return (this->is_in_x_range(p));
    }

    // The supporting curve is a circle:
    // Check whether p lies on the upper or on the lower part of the circle.
    Comparison_result c_res = CGAL::compare(p.y(), y0());

    if ((_is_upper() && c_res == SMALLER) || (! _is_upper() && c_res == LARGER))
    {
      // The point lies on the other half of the circle:
      return false;
    }

    // Check if the point is in the x-range of the arc.
    return (this->is_in_x_range(p));
  }

  /*! Check whether the given point lies in the interior of the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_strictly_between_endpoints(const Point_2& p) const
  {
    if (p.equals (m_source) || p.equals (m_target)) return false;
    return (_is_between_endpoints(p));
  }

  /*! Compute the overlap with a given arc having the same supporting curve.
   * \param cv The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \return Whether we found an overlap.
   */
  bool _compute_overlap(const Self& cv, Self& overlap) const
  {
    // Check if the two arcs are identical.
    if (is_linear()) {
      // In case of line segments we can swap the source and target:
      if (((m_source.equals(cv.m_source) && m_target.equals(cv.m_target)) ||
           (m_source.equals(cv.m_target) && m_target.equals(cv.m_source))))
      {
        overlap = cv;
        return true;
      }
    }
    else {
      if ((orientation() == cv.orientation() &&
           m_source.equals(cv.m_source) && m_target.equals(cv.m_target)) ||
          (orientation() != cv.orientation() &&
           m_source.equals(cv.m_target) && m_target.equals(cv.m_source)))
      {
        overlap = cv;
        return true;
      }
    }

    // Check for other overlaps:
    if (_is_strictly_between_endpoints(cv.left())) {
      if (_is_strictly_between_endpoints(cv.right())) {
        // Case 1 - *this:     +----------->
        //             cv:       +=====>
        overlap = cv;
        return true;
      }
      else {
        // Case 2 - *this:     +----------->
        //             cv:               +=====>
        overlap = *this;

        if (overlap.is_directed_right()) overlap.m_source = cv.left();
        else overlap.m_target = cv.left();

        return true;
      }
    }
    else if (_is_strictly_between_endpoints(cv.right())) {
      // Case 3 - *this:     +----------->
      //             cv:   +=====>
      overlap = *this;

      if (overlap.is_directed_right()) overlap.m_target = cv.right();
      else overlap.m_source = cv.right();

      return true;
    }
    else if (cv._is_between_endpoints(m_source) &&
             cv._is_between_endpoints(m_target) &&
             (cv._is_strictly_between_endpoints(m_source) ||
              cv._is_strictly_between_endpoints(m_target)))
    {
      // Case 4 - *this:     +----------->
      //             cv:   +================>
      overlap = *this;
      return true;
    }

    // If we reached here, there are no overlaps:
    return false;
  }

  public:
  template <class OutputIterator>
  void approximate(OutputIterator oi, unsigned int n) const
  {
    const double x_left = CGAL::to_double(this->source().x());
    const double y_left = CGAL::to_double(this->source().y());

    const double x_right = CGAL::to_double(this->target().x());
    const double y_right = CGAL::to_double(this->target().y());
    if (this->is_linear()) {
      *oi++ = std::make_pair(x_left, y_left);
      *oi++ = std::make_pair(x_right, y_right);
      return;
    }

    // Otherwise, sample (n - 1) equally-spaced points in between.
    const double app_xcenter = CGAL::to_double(this->m_first);
    const double app_ycenter = CGAL::to_double(this->m_second);
    const double app_sqr_rad = CGAL::to_double(this->m_third);

    const double x_jump = (x_right - x_left) / n;
    double x, y;
    double disc;
    unsigned int i;

    const bool is_up = this->_is_upper();
    *oi++ = std::make_pair (x_left, y_left);   // The left point.
    for (i = 1; i < n; ++i) {
      x = x_left + x_jump*i;
      disc = app_sqr_rad - CGAL::square(x - app_xcenter);
      if (disc < 0) disc = 0;
      if(is_up) y = app_ycenter + std::sqrt(disc);
      else y = app_ycenter - std::sqrt(disc);

      *oi++ = std::make_pair(x, y);
    }
    *oi++ = std::make_pair(x_right, y_right);   // The right point.
  }

  /*! Trim the arc given its new endpoints.
   * \param ps The new source point.
   * \param pt The new target point.
   * \return The new trimmed arc.
   * \pre Both ps and pt lies on the arc and must conform with the current
   *      direction of the arc.
   */
  Self trim(const Point_2& ps, const Point_2& pt) const
  {
    Self arc = *this;

    arc.m_source = ps;
    arc.m_target = pt;

    return arc;
  }

  //@}
};

/*! Exporter for circular arcs (or line segments).
 */
template <class Kernel, bool Filter>
std::ostream&
operator<<(std::ostream& os,
           const _X_monotone_circle_segment_2<Kernel, Filter> & arc)
{
  if (! arc.is_linear())
    os << "(" << arc.supporting_circle() << ") ";

  os << "[" << arc.source() << " --> " << arc.target() << "]";
  return (os);
}

} //namespace CGAL

#endif
