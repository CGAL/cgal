// Copyright (c) 1997-2010
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael.Hemmer@sophia.inria.fr

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/NT_converter.h>
#include <CGAL/number_utils.h>

#include <optional>
#include <variant>

namespace CGAL {
namespace Intersections {
namespace internal {

// Intersects a bbox with a segment, ray or line given by a point (lp) and a
// direction vector (ld). The flags select the kind of object:
//   segment: min_infinite = max_infinite = false
//   ray:     min_infinite = false, max_infinite = true
//   line:    min_infinite = max_infinite = true
//
// The bbox coordinates are doubles while the query coordinates are K::FT. To
// avoid losing precision on either side -- the old version converted the query
// to double via to_double() (issue #7124), and converting the bbox to K::FT
// loses information when K::FT is lower precision than double (e.g. float) --
// the clipping is computed in Coercion_traits<K::FT, double>::Type, the common
// type of both. This matches what do_intersect(Bbox_3, Segment_3) already does.
// The resulting point coordinates are converted back to K::FT.
template <class K>
typename std::optional< std::variant<typename K::Segment_3, typename K::Point_3 > >
intersection_bl(const Bbox_3& box,
                const typename K::FT& lpx, const typename K::FT& lpy, const typename K::FT& lpz,
                const typename K::FT& ldx, const typename K::FT& ldy, const typename K::FT& ldz,
                bool min_infinite, bool max_infinite)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename std::optional<std::variant< Segment_3, Point_3 > > result_type;

  typedef typename Coercion_traits<typename K::FT, double>::Type CFT;
  typename Coercion_traits<typename K::FT, double>::Cast to_CFT;

  const CFT clpx = to_CFT(lpx), clpy = to_CFT(lpy), clpz = to_CFT(lpz);
  const CFT cldx = to_CFT(ldx), cldy = to_CFT(ldy), cldz = to_CFT(ldz);

  CFT seg_min = 0, seg_max = 1;

  // first on x value
  if(cldx == 0) {
    if(clpx < to_CFT(box.xmin()))
      return result_type();
    if(clpx > to_CFT(box.xmax()))
      return result_type();
  } else {
    CFT newmin, newmax;
    if(cldx > 0) {
      newmin = (to_CFT(box.xmin())-clpx)/cldx;
      newmax = (to_CFT(box.xmax())-clpx)/cldx;
    } else {
      newmin = (to_CFT(box.xmax())-clpx)/cldx;
      newmax = (to_CFT(box.xmin())-clpx)/cldx;
    }

    if(min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if(newmin > seg_min)
        seg_min = newmin;
    }

    if(max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if(newmax < seg_max)
        seg_max = newmax;
    }

    if(seg_max < seg_min)
      return result_type();
  }

  // now on y value
  if(cldy == 0) {
    if(clpy < to_CFT(box.ymin()))
      return result_type();
    if(clpy > to_CFT(box.ymax()))
      return result_type();
  } else {
    CFT newmin, newmax;
    if(cldy > 0) {
      newmin = (to_CFT(box.ymin())-clpy)/cldy;
      newmax = (to_CFT(box.ymax())-clpy)/cldy;
    } else {
      newmin = (to_CFT(box.ymax())-clpy)/cldy;
      newmax = (to_CFT(box.ymin())-clpy)/cldy;
    }

    if(min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if(newmin > seg_min)
        seg_min = newmin;
    }

    if(max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if(newmax < seg_max)
        seg_max = newmax;
    }

    if(seg_max < seg_min)
      return result_type();
  }

  // now on z value
  if(cldz == 0) {
    if(clpz < to_CFT(box.zmin()))
      return result_type();
    if(clpz > to_CFT(box.zmax()))
      return result_type();
  } else {
    CFT newmin, newmax;
    if(cldz > 0) {
      newmin = (to_CFT(box.zmin())-clpz)/cldz;
      newmax = (to_CFT(box.zmax())-clpz)/cldz;
    } else {
      newmin = (to_CFT(box.zmax())-clpz)/cldz;
      newmax = (to_CFT(box.zmin())-clpz)/cldz;
    }

    if(min_infinite) {
      min_infinite = false;
      seg_min = newmin;
    } else {
      if(newmin > seg_min)
        seg_min = newmin;
    }

    if(max_infinite) {
      max_infinite = false;
      seg_max = newmax;
    } else {
      if(newmax < seg_max)
        seg_max = newmax;
    }

    if(seg_max < seg_min)
      return result_type();
  }

  if(min_infinite || max_infinite)
    seg_max = seg_min;

  NT_converter<CFT, typename K::FT> to_FT;

  if(seg_max == seg_min)
    return result_type(Point_3(to_FT(clpx + cldx*seg_min),
                               to_FT(clpy + cldy*seg_min),
                               to_FT(clpz + cldz*seg_min)));

  return result_type(Segment_3(Point_3(to_FT(clpx + cldx*seg_min),
                                       to_FT(clpy + cldy*seg_min),
                                       to_FT(clpz + cldz*seg_min)),
                               Point_3(to_FT(clpx + cldx*seg_max),
                                       to_FT(clpy + cldy*seg_max),
                                       to_FT(clpz + cldz*seg_max))));
}

template <class K>
typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type
intersection(const typename K::Segment_3& seg,
             const Bbox_3& box,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  const Point_3& linepoint = seg.source();
  const Vector_3 diffvec = seg.target() - linepoint;

  return intersection_bl<K>(box,
                            linepoint.x(), linepoint.y(), linepoint.z(),
                            diffvec.x(), diffvec.y(), diffvec.z(),
                            false, false);
}

template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Segment_3>::result_type
intersection(const Bbox_3& box,
             const typename K::Segment_3& seg,
             const K& k)
{
  return intersection(seg, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H
