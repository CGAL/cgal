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

#ifndef CGAL_INTERSECTIONS_BBOX_2_RAY_2_H
#define CGAL_INTERSECTIONS_BBOX_2_RAY_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <typename K>
class Bbox_2_Ray_2_pair;

template <typename K>
class Bbox_2_Ray_2_pair_impl
{
  typedef typename K::Point_2                                           Point_2;
  typedef typename K::Vector_2                                          Vector_2;
  typedef typename K::Ray_2                                             Ray_2;

public:
    Bbox_2_Ray_2_pair_impl():_known(false) {}
    Bbox_2_Ray_2_pair_impl(const Bbox_2& bbox,
                           const Point_2& pt,
                           const Vector_2& dir)
        :_box(bbox), _known(false), _ref_point(pt), _dir(dir), _min(0.0)
    {}

    Ray_2 _ray;
    Bbox_2 _box;
    bool _known;

    typename Bbox_2_Ray_2_pair<K>::Intersection_results _result;
    Point_2 _ref_point;
    Vector_2 _dir;

    double _min, _max;
};

template <typename K>
class Bbox_2_Ray_2_pair
{
  typedef Bbox_2_Ray_2_pair<K>                                          Self;
  typedef typename K::Point_2                                           Point_2;
  typedef typename K::Vector_2                                          Vector_2;

public:
  enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};

  ~Bbox_2_Ray_2_pair() { delete pimpl; }
  Bbox_2_Ray_2_pair() { pimpl = new Bbox_2_Ray_2_pair_impl<K>; }
  Bbox_2_Ray_2_pair(const Self& o) { pimpl = new Bbox_2_Ray_2_pair_impl<K>(*o.pimpl); }
  Bbox_2_Ray_2_pair(Bbox_2 const &bbox, double x, double y, double dx, double dy) {
    pimpl = new Bbox_2_Ray_2_pair_impl<K>(bbox, Point_2(x,y), Vector_2(dx,dy));
  }

  Self& operator=(const Self& o)
  {
    *pimpl = *o.pimpl;
    return *this;
  }

  Intersection_results intersection_type() const
  {
    if(pimpl->_known)
      return pimpl->_result;

    pimpl->_known = true;
    bool to_infinity = true;

    // first on x value
    if(pimpl->_dir.x() == 0.0)
    {
      if(pimpl->_ref_point.x() < pimpl->_box.xmin())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      if(pimpl->_ref_point.x() > pimpl->_box.xmax())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }
    }
    else
    {
      double newmin, newmax;
      if(pimpl->_dir.x() > 0.0)
      {
        newmin =(pimpl->_box.xmin()-pimpl->_ref_point.x())/pimpl->_dir.x();
        newmax =(pimpl->_box.xmax()-pimpl->_ref_point.x())/pimpl->_dir.x();
      }
      else
      {
        newmin =(pimpl->_box.xmax()-pimpl->_ref_point.x())/pimpl->_dir.x();
        newmax =(pimpl->_box.xmin()-pimpl->_ref_point.x())/pimpl->_dir.x();
      }

      if(newmin > pimpl->_min)
        pimpl->_min = newmin;

      if(to_infinity)
      {
        pimpl->_max = newmax;
      }
      else
      {
        if(newmax < pimpl->_max)
          pimpl->_max = newmax;
      }

      if(pimpl->_max < pimpl->_min)
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      to_infinity = false;
    }

    // now on y value
    if(pimpl->_dir.y() == 0.0)
    {
      if(pimpl->_ref_point.y() < pimpl->_box.ymin())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      if(pimpl->_ref_point.y() > pimpl->_box.ymax())
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }
    }
    else
    {
      double newmin, newmax;
      if(pimpl->_dir.y() > 0.0)
      {
        newmin =(pimpl->_box.ymin()-pimpl->_ref_point.y())/pimpl->_dir.y();
        newmax =(pimpl->_box.ymax()-pimpl->_ref_point.y())/pimpl->_dir.y();
      }
      else
      {
        newmin =(pimpl->_box.ymax()-pimpl->_ref_point.y())/pimpl->_dir.y();
        newmax =(pimpl->_box.ymin()-pimpl->_ref_point.y())/pimpl->_dir.y();
      }

      if(newmin > pimpl->_min)
        pimpl->_min = newmin;

      if(to_infinity)
      {
        pimpl->_max = newmax;
      }
      else
      {
        if(newmax < pimpl->_max)
          pimpl->_max = newmax;
      }

      if(pimpl->_max < pimpl->_min)
      {
        pimpl->_result = NO_INTERSECTION;
        return pimpl->_result;
      }

      to_infinity = false;
    }

    CGAL_kernel_assertion(!to_infinity);

    if(pimpl->_max == pimpl->_min)
    {
      pimpl->_result = POINT;
      return pimpl->_result;
    }

    pimpl->_result = SEGMENT;
    return pimpl->_result;
  }

  bool intersection(double& x, double& y) const
  {
    if(!pimpl->_known)
        intersection_type();

    if(pimpl->_result != POINT)
        return false;

    Point_2 pt = pimpl->_ref_point + pimpl->_min*pimpl->_dir;
    x = pt.x();
    y = pt.y();

    return true;
  }

  bool intersection(double& x1, double& y1, double& x2, double& y2) const
  {
    if(!pimpl->_known)
      intersection_type();

    if(pimpl->_result != SEGMENT)
      return false;

    Point_2 p1(pimpl->_ref_point + pimpl->_min*pimpl->_dir);
    Point_2 p2(pimpl->_ref_point + pimpl->_max*pimpl->_dir);
    x1 = p1.x();
    y1 = p1.y();
    x2 = p2.x();
    y2 = p2.y();

    return true;
  }

protected:
  Bbox_2_Ray_2_pair_impl<K>* pimpl;
};

template <typename K>
bool do_intersect_ray_2(const Bbox_2& bbox, double x, double y, double dx, double dy, const K& k)
{
  Bbox_2_Ray_2_pair<K> pair(bbox, x, y, dx, dy);
  return pair.intersection_type() != Bbox_2_Ray_2_pair<K>::NO_INTERSECTION;
}

template <typename K>
bool do_intersect_ray_2(const Bbox_2 &box, const Ray_2<K> &ray, const K& k = K())
{
  double startx = to_double(ray.start().x());
  double starty = to_double(ray.start().y());
  double dx = to_double(ray.direction().to_vector().x());
  double dy = to_double(ray.direction().to_vector().y());

  return do_intersect_ray_2(box, startx, starty, dx, dy, k);
}

template <typename K>
inline bool do_intersect_ray_2(const Ray_2<K> &ray, const Bbox_2 &box, const K& k = K())
{
  return do_intersect_ray_2(box, ray, k);
}

template <typename K>
inline bool do_intersect(const typename K::Ray_2& line, const Bbox_2& bbox, const K& k = K())
{
  return do_intersect_ray_2(bbox, line, k);
}

template <typename K>
inline bool do_intersect(const Bbox_2& bbox, const typename K::Ray_2& line, const K& k = K())
{
  return do_intersect_ray_2(bbox, line, k);
}

} // namespace internal
} // namespace Intersections

template <typename R>
inline bool do_intersect(const Ray_2<R>& line, const Bbox_2& box)
{
  return Intersections::internal::do_intersect(box, line);
}

template <typename R>
inline bool do_intersect(const Bbox_2& box, const Ray_2<R>& line)
{
  return Intersections::internal::do_intersect(box, line);
}

} // namespace CGAL

#endif
