// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Coercion_traits.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/Same_uncertainty.h>
#include <CGAL/kernel_assertions.h>

#include <type_traits>

// inspired from https://people.csail.mit.edu/amy/papers/box-jgt.pdf

// This algorithm intersects the line with the x-, y-, and z-slabs of the
// bounding box, and computes the interval [t1, t2], in the
// parameterization of the line given by the segment (for t=0, that is the
// source of the segment, and for t=1 that is its target), where the line
// intersects the three slabs of the bounding box.

// For a segment, the intersection is non-empty iff
//    [t1, t2] intersects [0, 1].

namespace CGAL {
namespace Intersections {
namespace internal {

template <typename FT, bool bounded_0, bool use_static_filters = false>
struct Do_intersect_bbox_segment_aux_is_greater
{
  typedef typename Same_uncertainty<bool, FT>::type result_type;

  void register_new_input_values(const FT&, const FT&) {}
  void compute_new_error_bound() { }
  bool bound_overflow() { return false; }
  bool value_might_underflow() { return false; }

  static result_type uncertain()
  {
    return true;
  }

  result_type operator()(const FT& a, const FT& b) const
  {
    return a > b;
  }
};

template <typename FT, bool bounded_0>
class Do_intersect_bbox_segment_aux_is_greater<FT, bounded_0, true>
{
  double error;
  double tmax;
  double dmax;

public:
  static_assert(std::is_same<FT, double>::value);

  Do_intersect_bbox_segment_aux_is_greater() : error(0.), tmax(0.), dmax(0.) {}

  void register_new_input_values(const double t, const double d)
  {
    if(bounded_0)
    {
      if(t > tmax)
        tmax = t;
      if(d > dmax)
        dmax = d;
    }
    else
    {
      const double at = CGAL::abs(t);
      const double ad = CGAL::abs(d);
      if(at > tmax)
        tmax = at;
      if(ad > dmax)
        dmax = ad;
    }
  }

  void compute_new_error_bound()
  {
    const double EPS = 8.8872057372592798e-16;
    error = tmax * dmax * EPS;
  }

  bool bound_overflow()
  {
    const double OVERF = 1e153;
    return dmax > OVERF || tmax > OVERF;
  }

  bool value_might_underflow()
  {
    const double UNDERF = 1e-146;
    return dmax < UNDERF || tmax < UNDERF;
  }

  typedef Uncertain<bool> result_type;

  static result_type uncertain()
  {
    return result_type::indeterminate();
  }

  result_type operator()(const FT& a, const FT& b) const
  {
    const FT x = a - b;
    if(x > error)
      return true;
    else if(x < -error)
      return false;
    else
      return uncertain();
  }
}; // end specialization Do_intersect_bbox_segment_aux_is_greater<FT, true>

template <typename FT,
          typename BFT,
          bool bounded_0,
          bool bounded_1,
          bool use_static_filters>
inline
typename Do_intersect_bbox_segment_aux_is_greater<FT, bounded_0, use_static_filters>::result_type
do_intersect_bbox_segment_aux(const FT& px, const FT& py, const FT& pz,
                              const FT& qx, const FT& qy, const FT& qz,
                              const BFT& bxmin, const BFT& bymin, const BFT& bzmin,
                              const BFT& bxmax, const BFT& bymax, const BFT& bzmax)
{
  if(((px >= bxmin) && (px <= bxmax) &&
      (py >= bymin) && (py <= bymax) &&
      (pz >= bzmin) && (pz <= bzmax)) ||
     ((qx >= bxmin) && (qx <= bxmax) &&
      (qy >= bymin) && (qy <= bymax) &&
      (qz >= bzmin) && (qz <= bzmax)))
  {
    return true;
  }

  // The following code encode t1 and t2 by:
  //    t1 = tmin/dmin
  //    t2 = tmax/dmax
  // For the first lines, dmax==dmin and is not explicitly defined.

  // -----------------------------------
  // treat x coord
  // -----------------------------------
  typedef typename Coercion_traits<FT,BFT>::Type CFT;
  typename Coercion_traits<FT,BFT>::Cast to_CFT;
  CFT dmin, tmin, tmax, dmax;
  if(qx >= px)
  {
    if(bounded_0 && compare(px, bxmax) == LARGER)
      return false; // segment on the right of bbox
    if(bounded_1 && compare(qx, bxmin) == SMALLER)
      return false; // segment on the left of bbox

    if(bounded_1 && bxmax > qx)
    {
      tmax = 1;
      dmax = 1;
    }
    else
    {
      tmax = to_CFT(bxmax) - px;
      dmax = qx - px;
    }

    tmin = CFT(bxmin) - px;
    dmin = qx - px;
  }
  else
  {
    if(bounded_1 && compare(qx, bxmax) == LARGER)
      return false; // segment on the right of bbox
    if(bounded_0 && compare(px, bxmin) == SMALLER)
      return false; // segment on the left of bbox

    if(bounded_1 && compare(bxmin, qx) == SMALLER)
    {
      tmax = 1;
      dmax = 1;
    }
    else
    {
      tmax = px - to_CFT(bxmin);
      dmax = px - qx;
    }

    tmin = px - to_CFT(bxmax);
    dmin = px - qx;
  }

  if(bounded_0)
    tmin = (CGAL::max)(CFT(0), tmin);

  // If the query is vertical for x, then check its x-coordinate is in
  // the x-slab.
  if((px == qx) &&  // <=> (dmin == 0)
     (!(bounded_0 && bounded_1))) // do not check for a segment
  {
    if(compare(px, bxmax) == LARGER || compare(px, bxmin) == SMALLER)
      return false;

    // Note: for a segment the condition has already been tested by the two
    // previous tests tmax<0 || tmin>dmin (with dmin==0).
  }

  // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
  // is a NaN. But the case with NaNs is treated as if the interval
  // [t1, t2] was ]-inf, +inf[.

  CGAL_assertion(! is_negative(dmin));
  CGAL_assertion(! is_negative(dmax));
  CGAL_assertion_code(if(bounded_0) {)
    CGAL_assertion(! is_negative(tmin));
    CGAL_assertion(! is_negative(tmax));
  CGAL_assertion_code(})

  // -----------------------------------
  // treat y coord
  // -----------------------------------
  CFT dymin, tymin, tymax, dymax;
  if(qy >= py)
  {
    if(bounded_0 && compare(py, bymax) == LARGER)
      return false; // segment on the right of bbox
    if(bounded_1 && compare(qy, bymin) == SMALLER)
      return false; // segment on the left of bbox

    if(bounded_1 && compare(bymax, qy) == LARGER)
    {
      tymax = 1;
      dymax = 1;
    }
    else
    {
      tymax = to_CFT(bymax) - py;
      dymax = qy - py;
    }

    tymin = to_CFT(bymin) - py;
    dymin = qy - py;
  }
  else
  {
    if(bounded_1 && compare(qy, bymax) == LARGER)
      return false; // segment on the right of bbox
    if(bounded_0 && compare(py, bymin) == SMALLER)
      return false; // segment on the left of bbox

    if(bounded_1 && compare(bymin, qy) == SMALLER)
    {
      tymax = 1;
      dymax = 1;
    }
    else
    {
      tymax = py - to_CFT(bymin);
      dymax = py - qy;
    }

    tymin = py - CFT(bymax);
    dymin = py - qy;
  }

  if(bounded_0)
    tymin = (CGAL::max)(CFT(0), tymin);

  // If the query is vertical for y, then check its y-coordinate is in
  // the y-slab.
  if((py == qy) &&  // <=> (dmin == 0)
     (! (bounded_0 && bounded_1))) // do not check for a segment
  {
    if(py > to_CFT(bymax) || compare(py, bymin) == SMALLER)
      return false;
  }

  // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
  // is a NaN. But the case with NaNs is treated as if the interval
  // [t1, t2] was ]-inf, +inf[.

  CGAL_assertion(! is_negative(dymin));
  CGAL_assertion(! is_negative(dymax));
  if(bounded_0)
  {
    CGAL_assertion(! is_negative(tymin));
    CGAL_assertion(! is_negative(tymax));
  }

  // -----------------------------------
  // treat z coord
  // -----------------------------------
  CFT dzmin, tzmin, tzmax, dzmax;
  if(qz >= pz)
  {
    if(bounded_0 && compare(pz, bzmax)== LARGER)
      return false; // segment on the right of bbox
    if(bounded_1 && compare(qz, bzmin) == SMALLER)
      return false; // segment on the left of bbox

    if(bounded_1 && compare(bzmax, qz) == LARGER)
    {
      tzmax = 1;
      dzmax = 1;
    }
    else
    {
      tzmax = to_CFT(bzmax) - pz;
      dzmax = qz - pz;
    }

    tzmin = to_CFT(bzmin) - pz;
    dzmin = qz - pz;
  }
  else
  {
    if(bounded_1 && compare(qz, bzmax) == LARGER)
      return false; // segment on the right of bbox
    if(bounded_0 && compare(pz, bzmin) == SMALLER)
    return false; // segment on the left of bbox

    if(bounded_1 && compare(bzmin, qz) == SMALLER)
    {
      tzmax = 1;
      dzmax = 1;
    }
    else
    {
      tzmax = pz - to_CFT(bzmin);
      dzmax = pz - qz;
    }

    tzmin = pz - to_CFT(bzmax);
    dzmin = pz - qz;
  }

  if(bounded_0)
    tzmin = (CGAL::max)(CFT(0), tzmin);

  // If the query is vertical for z, then check its z-coordinate is in
  // the z-slab.
  if((pz == qz) &&  // <=> (dmin == 0)
     (! (bounded_0 && bounded_1))) // do not check for a segment
  {
    if(compare(pz, bzmax) == LARGER || compare(pz, bzmin) == SMALLER)
      return false;
  }

  // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
  // is a NaN. But the case with NaNs is treated as if the interval
  // [t1, t2] was ]-inf, +inf[.

  CGAL_assertion(! is_negative(dzmin));
  CGAL_assertion(! is_negative(dzmax));
  CGAL_assertion_code(if(bounded_0) {)
    CGAL_assertion(! is_negative(tzmin));
    CGAL_assertion(! is_negative(tzmax));
  CGAL_assertion_code(})

  typedef Do_intersect_bbox_segment_aux_is_greater<CFT, bounded_0, use_static_filters> Is_greater;
  typedef typename Is_greater::result_type Is_greater_value;
  Is_greater is_greater;

  is_greater.register_new_input_values(tmin, dmin);
  is_greater.register_new_input_values(tymin, dymin);
  is_greater.register_new_input_values(tmax, dmax);
  is_greater.register_new_input_values(tymax, dymax);

  is_greater.compute_new_error_bound();
  if(is_greater.bound_overflow() || is_greater.value_might_underflow())
    return Is_greater::uncertain();

  // If t1 > tymax/dymax || tymin/dymin > t2, return false.
  if( py != qy && px != qx) // dmin > 0, dymax >0, dmax > 0, dymin > 0
  {
    const Is_greater_value b1 = is_greater(dymax* tmin,  dmin*tymax);
    if(possibly(b1))
      return !b1; // if(is_greater) return false; // or uncertain

    const Is_greater_value b2 = is_greater( dmax*tymin, dymin* tmax);
    if(possibly(b2))
      return !b2;
  }

  Is_greater_value b = Is_greater_value();

  // If tymin/dymin > t1, set t1 = tymin/dymin.
  if((px == qx) || // <=> (dmin == 0)
     ((py != qy) && // <=> (dymin > 0)
      certainly(b = is_greater(dmin*tymin, dymin*tmin))))
  {
    tmin = tymin;
    dmin = dymin;
  }

  if(is_indeterminate(b))
    return b; // Note that the default-constructed Is_greater_value cannot be indeterminate.

  // If tymax/dymax < t2, set t2 = tymax/dymax.
  if((px == qx) || // <=> (dmax > 0)
     ((py != qy) && // <=> dymax > 0
      certainly(b = is_greater(dymax*tmax, dmax*tymax))))
  {
    tmax = tymax;
    dmax = dymax;
  }

  if(is_indeterminate(b))
    return b;

  CGAL_assertion(! is_negative(dmin));
  CGAL_assertion(! is_negative(dmax));

  // If t1 > tzmax || tzmin > t2, return false.
  if((px != qx || py != qy) &&
     (pz != qz)) // dmin > 0, dmax > 0, dzmax > 0, dzmin > 0
  {
    is_greater.register_new_input_values(tzmin, dzmin);
    is_greater.register_new_input_values(tzmax, dzmax);

    is_greater.compute_new_error_bound();
    if(is_greater.bound_overflow() || is_greater.value_might_underflow())
      return Is_greater::uncertain();

    const Is_greater_value b1 = is_greater(dzmax*tmin, dmin*tzmax);
    if(possibly(b1))
      return !b1; // if(is_greater) return false; // or uncertain

    const Is_greater_value b2 = is_greater( dmax*tzmin, dzmin* tmax);
    if(possibly(b2))
      return !b2; // if(is_greater) return false; // or uncertain
  }

  return true;
}

template <typename FT,
          bool bounded_0,
          bool bounded_1,
          bool use_static_filters>
inline
typename Do_intersect_bbox_segment_aux_is_greater<FT, bounded_0, use_static_filters>::result_type
do_intersect_bbox_segment_aux(const FT& px, const FT& py, const FT& pz,
                              const FT& qx, const FT& qy, const FT& qz,
                              const Bbox_3& bb)

{
  return do_intersect_bbox_segment_aux<FT,double,bounded_0,bounded_1,use_static_filters>(px, py, pz,
                                                                                         qx, qy, qz,
                                                                                         bb.xmin(), bb.ymin(), bb.zmin(),
                                                                                         bb.xmax(), bb.ymax(), bb.zmax());
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Segment_3& segment,
             const CGAL::Bbox_3& bbox,
             const K&)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  const Point_3& source = segment.source();
  const Point_3& target = segment.target();

  return do_intersect_bbox_segment_aux<FT, true, true, false>(source.x(), source.y(), source.z(),
                                                              target.x(), target.y(), target.z(),
                                                              bbox);
}

template <class K>
typename K::Boolean
do_intersect(const CGAL::Bbox_3& bbox,
             const typename K::Segment_3& segment,
             const K& k)
{
  return do_intersect(segment, bbox, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
