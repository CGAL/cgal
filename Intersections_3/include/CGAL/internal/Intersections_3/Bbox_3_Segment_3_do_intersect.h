// Copyright (c) 2012  GeometryFactory Sarl (France)
// All rights reserved.
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
//
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Segment_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/Same_uncertainty.h>
#include <CGAL/assertions.h>
#include <CGAL/Coercion_traits.h>
#include <boost/type_traits/is_same.hpp>

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

// This algorithm intersects the line with the x-, y-, and z-slabs of the
// bounding box, and computes the interval [t1, t2], in the
// parameterization of the line given by the segment (for t=0, that is the
// source of the segment, and for t=1 that is its target), where the line
// intersects the three slabs of the bounding box.

// For a segment, the intersection is non-empty iff 
//    [t1, t2] intersects [0, 1].

namespace CGAL {

namespace internal {

  template <typename FT, bool bounded_0, bool use_static_filters = false>
  struct Do_intersect_bbox_segment_aux_is_greater 
  {
    typedef typename Same_uncertainty<bool, FT>::type result_type;

    void register_new_input_values(const FT&, const FT&) {}
    void compute_new_error_bound() {}
    bool bound_overflow() { return false; }
    bool value_might_underflow() { return false; }

    static result_type uncertain() {
      return true;
    }

    result_type operator()(const FT& a, const FT& b) const {
      return a > b;
    }
  }; // end struct template Do_intersect_bbox_segment_aux_is_greater

  template <typename FT, bool bounded_0>
  class Do_intersect_bbox_segment_aux_is_greater<FT, bounded_0, true>
  {
    double error;
    double tmax;
    double dmax;

  public:
    CGAL_static_assertion((boost::is_same<FT, double>::value));

    Do_intersect_bbox_segment_aux_is_greater() : error(0.), tmax(0.), dmax(0.) {}

    void register_new_input_values(const double& t, const double& d) {
      if(bounded_0) {
        if(t > tmax) tmax = t;
        if(d > dmax) dmax = d;
      } else {
        const double at = CGAL::abs(t);
        const double ad = CGAL::abs(d);
        if(at > tmax) tmax = at;
        if(ad > dmax) dmax = ad;
      }
    }

    void compute_new_error_bound() {
      const double EPS = 8.8872057372592798e-16;
      error = tmax * dmax * EPS;
    }

    bool bound_overflow() {
      const double OVERF = 1e153;
      return dmax > OVERF || tmax > OVERF;
    }

    bool value_might_underflow() {
      const double UNDERF = 1e-146;
      return dmax < UNDERF || tmax < UNDERF;
    }

    typedef Uncertain<bool> result_type;

    static result_type uncertain() {
      return result_type::indeterminate();
    }

    result_type operator()(const FT& a, const FT& b) const {
      const FT  x = a - b;
      if(x > error) return true;
      else if(x < -error) return false;
      else return uncertain();
    }

  }; // end specialization Do_intersect_bbox_segment_aux_is_greater<FT, true>

  template <typename FT,
            bool bounded_0,
            bool bounded_1,
            bool use_static_filters>
  inline
  typename Do_intersect_bbox_segment_aux_is_greater
  <
    FT,
    bounded_0,
    use_static_filters
  >::result_type
  do_intersect_bbox_segment_aux(
                          const FT& px, const FT& py, const FT& pz,
                          const FT& qx, const FT& qy, const FT& qz,
                          const Bbox_3& bbox)
  {
    const double& bxmin = bbox.xmin();
    const double& bymin = bbox.ymin();
    const double& bzmin = bbox.zmin();
    const double& bxmax = bbox.xmax();
    const double& bymax = bbox.ymax();
    const double& bzmax = bbox.zmax();


    if( ( (px >= bxmin) && (px <= bxmax) &&
          (py >= bymin) && (py <= bymax) &&
          (pz >= bzmin) && (pz <= bzmax) ) ||
        ( (qx >= bxmin) && (qx <= bxmax) &&
          (qy >= bymin) && (qy <= bymax) &&
          (qz >= bzmin) && (qz <= bzmax) ) ) {
      return true;
    }

    // The following code encode t1 and t2 by:
    //    t1 = tmin/dmin
    //    t2 = tmax/dmax
    // For the first lines, dmax==dmin and is not explicitly defined.

    // -----------------------------------
    // treat x coord
    // -----------------------------------
    typedef typename Coercion_traits<double,FT>::Type CFT;
    CFT dmin, tmin, tmax, dmax;
    if ( qx >= px )
    {
      if(bounded_0 && px > bxmax) return false; // segment on the right of bbox
      if(bounded_1 && qx < bxmin) return false; // segment on the left of bbox

      if(bounded_1 && bxmax > qx) {
        tmax = 1;
        dmax = 1;
      } else {
        tmax = bxmax - px;
        dmax = qx - px;
      }

      tmin = bxmin - px;
      dmin = qx - px;
    }
    else
    {
      if(bounded_1 && qx > bxmax) return false; // segment on the right of bbox
      if(bounded_0 && px < bxmin) return false; // segment on the left of bbox

      if(bounded_1 && bxmin < qx) {
        tmax = 1;
        dmax = 1;
      } else {
        tmax = px - bxmin;
        dmax = px - qx;
      }

      tmin = px - bxmax;
      dmin = px - qx;
    }

    if(bounded_0) tmin = (CGAL::max)(CFT(0), tmin);

    // If the query is vertical for x, then check its x-coordinate is in
    // the x-slab.
    if( (px == qx) &&  // <=> (dmin == 0)
        (! (bounded_0 && bounded_1) ) ) // do not check for a segment
    {
      if(px > bxmax || px < bxmin) return false;
      // Note: for a segment the condition has already been tested by the two
      // previous tests tmax<0 || tmin>dmin (with dmin==0).
    }

    // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
    // is a NaN. But the case with NaNs is treated as if the interval
    // [t1, t2] was ]-inf, +inf[.

    CGAL_assertion(dmin >= 0);
    CGAL_assertion(dmax >= 0);
    if(bounded_0) {
      CGAL_assertion(tmin >= 0);
      CGAL_assertion(tmax >= 0);
    }


    // -----------------------------------
    // treat y coord
    // -----------------------------------
    CFT dymin, tymin, tymax, dymax;
    if ( qy >= py )
    {
      if(bounded_0 && py > bymax) return false; // segment on the right of bbox
      if(bounded_1 && qy < bymin) return false; // segment on the left of bbox

      if(bounded_1 && bymax > qy) {
        tymax = 1;
        dymax = 1;
      } else {
        tymax = bymax - py;
        dymax = qy - py;
      }

      tymin = bymin - py;
      dymin = qy - py;
    }
    else
    {
      if(bounded_1 && qy > bymax) return false; // segment on the right of bbox
      if(bounded_0 && py < bymin) return false; // segment on the left of bbox

      if(bounded_1 && bymin < qy) {
        tymax = 1;
        dymax = 1;
      } else {
        tymax = py - bymin;
        dymax = py - qy;
      }

      tymin = py - bymax;
      dymin = py - qy;
    }

    if(bounded_0) tymin = (CGAL::max)(CFT(0), tymin);

    // If the query is vertical for y, then check its y-coordinate is in
    // the y-slab.
    if( (py == qy) &&  // <=> (dmin == 0)
        (! (bounded_0 && bounded_1) ) ) // do not check for a segment
    {
      if(py > bymax || py < bymin) return false;
    }

    // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
    // is a NaN. But the case with NaNs is treated as if the interval
    // [t1, t2] was ]-inf, +inf[.

    CGAL_assertion(dymin >= 0);
    CGAL_assertion(dymax >= 0);
    if(bounded_0) {
      CGAL_assertion(tymin >= 0);
      CGAL_assertion(tymax >= 0);
    }


    // -----------------------------------
    // treat z coord
    // -----------------------------------
    CFT dzmin, tzmin, tzmax, dzmax;
    if ( qz >= pz )
    {
      if(bounded_0 && pz > bzmax) return false; // segment on the right of bbox
      if(bounded_1 && qz < bzmin) return false; // segment on the left of bbox

      if(bounded_1 && bzmax > qz) {
        tzmax = 1;
        dzmax = 1;
      } else {
        tzmax = bzmax - pz;
        dzmax = qz - pz;
      }

      tzmin = bzmin - pz;
      dzmin = qz - pz;
    }
    else
    {
      if(bounded_1 && qz > bzmax) return false; // segment on the right of bbox
      if(bounded_0 && pz < bzmin) return false; // segment on the left of bbox

      if(bounded_1 && bzmin < qz) {
        tzmax = 1;
        dzmax = 1;
      } else {
        tzmax = pz - bzmin;
        dzmax = pz - qz;
      }

      tzmin = pz - bzmax;
      dzmin = pz - qz;
    }

    if(bounded_0) tzmin = (CGAL::max)(CFT(0), tzmin);

    // If the query is vertical for z, then check its z-coordinate is in
    // the z-slab.
    if( (pz == qz) &&  // <=> (dmin == 0)
        (! (bounded_0 && bounded_1) ) ) // do not check for a segment
    {
      if(pz > bzmax || pz < bzmin) return false;
    }

    // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
    // is a NaN. But the case with NaNs is treated as if the interval
    // [t1, t2] was ]-inf, +inf[.

    CGAL_assertion(dzmin >= 0);
    CGAL_assertion(dzmax >= 0);
    if(bounded_0) {
      CGAL_assertion(tzmin >= 0);
      CGAL_assertion(tzmax >= 0);
    }


    typedef Do_intersect_bbox_segment_aux_is_greater
      <CFT,
       bounded_0,
       use_static_filters
       > Is_greater;
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
    if( py != qy && px != qx ) { // dmin > 0, dymax >0, dmax > 0, dymin > 0
      const Is_greater_value b1 = is_greater(dymax* tmin,  dmin*tymax);
      if(possibly(b1)) return !b1; // if(is_greater) return false; // or uncertain
      const Is_greater_value b2 = is_greater( dmax*tymin, dymin* tmax);
      if(possibly(b2)) return !b2;
    }

    Is_greater_value b = Is_greater_value();
    // If tymin/dymin > t1, set t1 = tymin/dymin.
    if( (px == qx) || // <=> (dmin == 0)
        ( (py != qy) && // <=> (dymin > 0)
          certainly(b = is_greater( dmin*tymin, dymin* tmin)) ) )
    {
      tmin = tymin;
      dmin = dymin;
    }
    if(is_indeterminate(b)) return b; // Note that the default-constructed
                                      // Is_greater_value cannot be
                                      // indeterminate.

    // If tymax/dymax < t2, set t2 = tymax/dymax.
    if( (px == qx) || // <=> (dmax > 0)
        ( (py != qy) && // <=> dymax > 0
          certainly(b = is_greater(dymax* tmax,  dmax*tymax)) ) )
    {
      tmax = tymax;
      dmax = dymax;
    }
    if(is_indeterminate(b)) return b;

    CGAL_assertion(dmin >= 0);
    CGAL_assertion(dmax >= 0);

    // If t1 > tzmax || tzmin > t2, return false.
    if( (px != qx ||
         py != qy ) &&
        (pz != qz) ) // dmin > 0, dmax > 0, dzmax > 0, dzmin > 0
    {
      is_greater.register_new_input_values(tzmin, dzmin);
      is_greater.register_new_input_values(tzmax, dzmax);

      is_greater.compute_new_error_bound();
      if(is_greater.bound_overflow() || is_greater.value_might_underflow())
        return Is_greater::uncertain();

      const Is_greater_value b1 = is_greater(dzmax* tmin,  dmin*tzmax);
      if(possibly(b1)) return !b1; // if(is_greater) return false; // or uncertain
      const Is_greater_value b2 = is_greater( dmax*tzmin, dzmin* tmax);
      if(possibly(b2)) return !b2; // if(is_greater) return false; // or uncertain
    }
    return true;
  }

  template <class K>
  bool do_intersect(const typename K::Segment_3& segment,
    const CGAL::Bbox_3& bbox,
    const K&)
  {
    typedef typename K::FT FT;
    typedef typename K::Point_3 Point_3;

    const Point_3& source = segment.source();
    const Point_3& target = segment.target();

    return do_intersect_bbox_segment_aux<FT, true, true, false>(
                          source.x(), source.y(), source.z(),
                          target.x(), target.y(), target.z(),
                          bbox);
  }

  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bbox,
                    const typename K::Segment_3& segment,
                    const K& k)
  {
    return do_intersect(segment, bbox, k);
  }

} // namespace internal

template<typename K>
bool do_intersect(const CGAL::Bbox_3 a,
                  const Segment_3<K>& b) {
  return K().do_intersect_3_object()(a, b);
}

template<typename K>
bool do_intersect(const Segment_3<K>& a,
                  const CGAL::Bbox_3& b) {
  return K().do_intersect_3_object()(a, b);
}


} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
