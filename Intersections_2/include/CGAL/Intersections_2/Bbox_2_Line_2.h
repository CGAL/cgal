// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_BBOX_2_LINE_2_H
#define CGAL_INTERSECTIONS_2_BBOX_2_LINE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <typename K>
class Bbox_2_Line_2_pair;

template <typename K>
class Bbox_2_Line_2_pair_impl
{
  typedef typename K::Line_2                      Line_2;

public:
  Bbox_2_Line_2_pair_impl() {}
  Bbox_2_Line_2_pair_impl(const Bbox_2& bb, const Line_2& line)
    : _bbox(bb), _line(line), _known(false)
  {}

  Bbox_2 _bbox;
  Line_2 _line;

  mutable bool _known;
  mutable typename Bbox_2_Line_2_pair<K>::Intersection_results _result;
  mutable double _min, _max;
};

template <typename K>
class Bbox_2_Line_2_pair
{
  typedef Bbox_2_Line_2_pair<K>                     Self;

  typedef typename K::Point_2                       Point_2;
  typedef typename K::Vector_2                      Vector_2;
  typedef typename K::Line_2                        Line_2;

public:
  enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};

  Bbox_2_Line_2_pair()
  {
    pimpl = new Bbox_2_Line_2_pair_impl<K>;
    pimpl->_known = false;
  }
  Bbox_2_Line_2_pair(const Self& o) { pimpl = new Bbox_2_Line_2_pair_impl<K>(*o.pimpl); }
  Bbox_2_Line_2_pair(const Bbox_2& bbox, double line_a, double line_b, double line_c) {
    pimpl = new Bbox_2_Line_2_pair_impl<K>(bbox, Line_2(line_a, line_b, line_c));
  }

  ~Bbox_2_Line_2_pair() { delete pimpl; }

  Self& operator=(const Self& o)
  {
    *pimpl = *o.pimpl;
    return *this;
  }

  Intersection_results intersection_type() const
  {
    if (pimpl->_known)
      return pimpl->_result;

    // The non const this pointer is used to cast away const.
    pimpl->_known = true;
    const Point_2 &ref_point = pimpl->_line.point();
    const Vector_2 &dir = pimpl->_line.direction().to_vector();
    bool to_infinity = true;

    // first on x value
    if (dir.x() == 0.0)
    {
      if (ref_point.x() < pimpl->_bbox.xmin())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      if (ref_point.x() > pimpl->_bbox.xmax())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }
    }
    else
    {
      double newmin, newmax;
      if (dir.x() > 0.0)
      {
        newmin = (pimpl->_bbox.xmin()-ref_point.x())/dir.x();
        newmax = (pimpl->_bbox.xmax()-ref_point.x())/dir.x();
      }
      else
      {
        newmin = (pimpl->_bbox.xmax()-ref_point.x())/dir.x();
        newmax = (pimpl->_bbox.xmin()-ref_point.x())/dir.x();
      }

      if (to_infinity)
      {
        pimpl->_min = newmin;
        pimpl->_max = newmax;
      }
      else
      {
        if (newmin > pimpl->_min)
          pimpl->_min = newmin;
        if (newmax < pimpl->_max)
          pimpl->_max = newmax;
        if (pimpl->_max < pimpl->_min)
        {
          pimpl->_result = NO_INTERSECTION;
          return pimpl->_result;
        }
      }

      to_infinity = false;
    }

    // now on y value
    if (dir.y() == 0.0)
    {
      if (ref_point.y() < pimpl->_bbox.ymin())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      if (ref_point.y() > pimpl->_bbox.ymax())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }
    }
    else
    {
      double newmin, newmax;
      if (dir.y() > 0.0)
      {
        newmin = (pimpl->_bbox.ymin()-ref_point.y())/dir.y();
        newmax = (pimpl->_bbox.ymax()-ref_point.y())/dir.y();
      }
      else
      {
        newmin = (pimpl->_bbox.ymax()-ref_point.y())/dir.y();
        newmax = (pimpl->_bbox.ymin()-ref_point.y())/dir.y();
      }

      if (to_infinity)
      {
        pimpl->_min = newmin;
        pimpl->_max = newmax;
      }
      else
      {
        if (newmin > pimpl->_min)
          pimpl->_min = newmin;
        if (newmax < pimpl->_max)
          pimpl->_max = newmax;
        if (pimpl->_max < pimpl->_min)
        {
          pimpl->_result = NO_INTERSECTION;
          return pimpl->_result;
        }
      }

      to_infinity = false;
    }

    CGAL_kernel_assertion(!to_infinity);
    if (pimpl->_max == pimpl->_min)
    {
      pimpl->_result = POINT;
      return pimpl->_result;
    }

    pimpl->_result = SEGMENT;
    return pimpl->_result;
  }

  bool intersection(double& x, double& y) const
  {
    if (!pimpl->_known)
      intersection_type();
    if (pimpl->_result != POINT)
      return false;

    Point_2 pt(pimpl->_line.point() + pimpl->_min*pimpl->_line.direction().to_vector());
    x = pt.x();
    y = pt.y();

    return true;
  }

  bool intersection(double& x1, double& y1, double& x2, double& y2) const
  {
    if (!pimpl->_known)
      intersection_type();
    if (pimpl->_result != SEGMENT)
      return false;

    Point_2 p1(pimpl->_line.point() + pimpl->_min*pimpl->_line.direction().to_vector());
    Point_2 p2(pimpl->_line.point() + pimpl->_max*pimpl->_line.direction().to_vector());
    x1 = p1.x();
    y1 = p1.y();
    x2 = p2.x();
    y2 = p2.y();

    return true;
  }

protected:
  Bbox_2_Line_2_pair_impl<K> *pimpl;
};

template <typename K>
inline bool do_intersect_line_2(const Bbox_2 &box, double line_a, double line_b, double line_c, const K& k = K())
{
  Bbox_2_Line_2_pair<K> pair(box, line_a, line_b, line_c);

  return pair.intersection_type() != Bbox_2_Line_2_pair<K>::NO_INTERSECTION;
}

template <typename K>
bool do_intersect_line_2(const Bbox_2& bbox, const Line_2<K>& line, const K& k = K())
{
  return do_intersect_line_2(bbox, to_double(line.a()), to_double(line.b()), to_double(line.c()), k);
}

template <typename K>
bool do_intersect_line_2(const Line_2<K>& line, const Bbox_2& bbox, const K& k = K())
{
  return do_intersect_line_2(bbox, to_double(line.a()), to_double(line.b()), to_double(line.c()), k);
}

template <typename K>
inline bool do_intersect(const typename K::Line_2& line, const Bbox_2& bbox, const K& k = K())
{
  return do_intersect_line_2(bbox, line, k);
}

template <typename K>
inline bool do_intersect(const Bbox_2& bbox, const typename K::Line_2& line, const K& k = K())
{
  return do_intersect_line_2(bbox, line, k);
}

} // namespace internal
} // namespace Intersections

template <typename R>
inline bool do_intersect(const Line_2<R>& line, const Bbox_2& box)
{
  return Intersections::internal::do_intersect(box, line);
}

template <typename R>
inline bool do_intersect(const Bbox_2& box, const Line_2<R>& line)
{
  return Intersections::internal::do_intersect(box, line);
}

} //namespace CGAL

#endif
