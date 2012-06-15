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

  template <typename FT, bool use_static_filters = false>
  struct Do_intersect_bbox_segment_aux_is_greater 
  {
    typedef typename Same_uncertainty<bool, FT>::type result_type;

    void register_new_input_values(const FT& t, const FT& d) {}
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

  template <typename FT>
  class Do_intersect_bbox_segment_aux_is_greater<FT, true>
  {
    double error;
    double tmax;
    double dmax;

    static const double EPS = 8.8872057372592798e-16;
    static const double OVERF = 1e153;
    static const double UNDERF = 1e-146;

  public:
    CGAL_static_assertion((boost::is_same<FT, double>::value));

    Do_intersect_bbox_segment_aux_is_greater() : error(0.), tmax(0.), dmax(0.) {}

    void register_new_input_values(const double& t, const double& d) {
      const double at = CGAL::abs(t);
      const double ad = CGAL::abs(d);
      if(at > tmax) tmax = at;
      if(ad > dmax) dmax = ad;
    }

    void compute_new_error_bound() {
      error = tmax * dmax * EPS;
    }

    bool bound_overflow() {
      return dmax > OVERF && tmax > OVERF;
    }

    bool value_might_underflow() {
      return dmax < UNDERF || dmax < UNDERF;
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
// __attribute__ ((noinline))
  inline
  typename Do_intersect_bbox_segment_aux_is_greater<FT, use_static_filters>::result_type
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

    // The following code encode t1 and t2 by:
    //    t1 = tmin/dmin
    //    t2 = tmax/dmax
    // For the first lines, dmax==dmin and is not explicitly defined.

    // -----------------------------------
    // treat x coord
    // -----------------------------------
    FT dmin, tmin, tmax, dmax, tmin_minus_dmin, tmax_minus_dmax;
    if ( qx >= px )
    {
      tmin = bxmin - px;
      tmax = bxmax - px;
      dmax = dmin = qx - px;
      tmin_minus_dmin = bxmin - qx;
      tmax_minus_dmax = bxmax - qx;
    }
    else
    {
      tmin = px - bxmax;
      tmax = px - bxmin;
      dmax = dmin = px - qx;
      tmin_minus_dmin = qx - bxmax;
      tmax_minus_dmax = qx - bxmin;
    }

    if( bounded_0 && tmax < 0. ) return false; // t2 < 0 
    if( bounded_1 && tmin_minus_dmin > 0. ) return false; // t1 > 1

    if(bounded_0) tmin = CGAL::max(tmin, FT(0)); // t1 = max(t1,0)
    if(bounded_1 && tmax_minus_dmax > 0)     {   // t2 = min(t2,1);
      tmax = FT(1);
      dmax = FT(1);
    }

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
    // CGAL_assertion(!bounded_0 || ( (dmin == 0) == (px == qx && (px == bxmax || px == bxmin)) ) );
    // CGAL_assertion(!bounded_1 || ( (dmax == 0) == (px == qx && (px == bxmax || px == bxmin)) ) );

    // -----------------------------------
    // treat y coord
    // -----------------------------------
    FT dymin, tymin, tymax, dymax, tymin_minus_dmin, tymax_minus_dmax;
    if ( qy >= py )
    {
      tymin = bymin - py;
      tymax = bymax - py;
      dymax = dymin = qy - py;
      tymin_minus_dmin = bymin - qy;
      tymax_minus_dmax = bymax - qy;
    }
    else
    {
      tymin = py - bymax;
      tymax = py - bymin;
      dymax = dymin = py - qy;
      tymin_minus_dmin = qy - bymax;
      tymax_minus_dmax = qy - bymin;
    }

    if( bounded_0 && tymax < 0. ) return false; // t2 < 0 
    if( bounded_1 && tymin_minus_dmin > 0. ) return false; // t1 > 1

    if(bounded_0) tymin = CGAL::max(tymin, FT(0)); // t1 = max(t1,0)
    if(bounded_1 && tymax_minus_dmax > 0)     {    // t2 = min(t2,1);
      tymax = FT(1);
      dymax = FT(1);
    }

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
    // CGAL_assertion(!bounded_0 || (dmin == 0) == (!bounded_0 && px == qx && py == qy));
    // CGAL_assertion(!bounded_1 || (dmax == 0) == (!bounded_1 && px == qx && py == qy));


    // -----------------------------------
    // treat z coord
    // -----------------------------------
    FT dzmin, tzmin, tzmax, dzmax, tzmax_minus_dmax, tzmin_minus_dmin;

    if ( qz >= pz )
    {
      tzmin = bzmin - pz;
      tzmax = bzmax - pz;
      dzmax = dzmin = qz - pz;
      tzmin_minus_dmin = bzmin - qz;
      tzmax_minus_dmax = bzmax - qz;
    }
    else
    {
      tzmin = pz - bzmax;
      tzmax = pz - bzmin;
      dzmax = dzmin = pz - qz;
      tzmin_minus_dmin = qz - bzmax;
      tzmax_minus_dmax = qz - bzmin;
    }

    if( bounded_0 && tzmax < 0. ) return false; // t2 < 0 
    if( bounded_1 && tzmin_minus_dmin > 0. ) return false; // t1 > 1

    if(bounded_0) tzmin = CGAL::max(tzmin, FT(0)); // t1 = max(t1,0)
    if(bounded_1 && tzmax_minus_dmax > 0)     {    // t2 = min(t2,1);
      tzmax = FT(1);
      dzmax = FT(1);
    }

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

    typedef Do_intersect_bbox_segment_aux_is_greater<FT, 
                                                     use_static_filters> Is_greater;
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

    Is_greater_value b;
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
    // CGAL_assertion((dmin == 0) == (px == qx && py == qy));
    // CGAL_assertion((dmax == 0) == (px == qx && py == qy));

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

  template <typename FT,
            bool bounded_0,
            bool bounded_1>
  inline
  bool
  do_intersect_bbox_segment_aux(const CGAL::cpp0x::array<FT, 6> seg,
                                const CGAL::cpp0x::array<double, 6> box)

  {
    const FT& px = seg[0];
    const FT& py = seg[1];
    const FT& pz = seg[2];
    const FT& qx = seg[3];
    const FT& qy = seg[4];
    const FT& qz = seg[5];
    const double& bxmin = box[0];
    const double& bymin = box[1];
    const double& bzmin = box[2];
    const double& bxmax = box[3];
    const double& bymax = box[4];
    const double& bzmax = box[5];
    // for(int i = 0; i < 3; ++i) {
    //   const int sign = seg[3+i] > seg[i]; // (qx > px)?
    //   if(bounded_0 && seg[3*(1-sign) + i] > box[3+i]) return false; // segment on the right of bbox
    //   if(bounded_1 && seg[3*sign + i] < box[i]) return false; // segment on the left of bbox
    // }

    return do_intersect_bbox_segment_aux<FT, bounded_0, bounded_1>
      (px, py, pz,
       qy, qy, qz,
       bxmin, bymin, bymax,
       bxmax, bymax, bzmax);
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
                          // bbox.xmin(), bbox.ymin(), bbox.zmin(),
                          // bbox.xmax(), bbox.ymax(), bbox.zmax() );

    // const CGAL::cpp0x::array<FT, 6> seg  = {source.x(), source.y(), source.z(),
    //                                         target.x(), target.y(), target.z() };
    // return do_intersect_bbox_segment_aux<FT, true, true>
    //   ( seg,
    //     *reinterpret_cast<const CGAL::cpp0x::array<double, 6>*>(&*bbox.cartesian_begin()) );
  }

  template <class K>
  bool do_intersect(const CGAL::Bbox_3& bbox,
                    const typename K::Segment_3& segment,
                    const K& k)
  {
    return do_intersect(segment, bbox, k);
  }

} // namespace internal

template <class K>
bool do_intersect(const CGAL::Segment_3<K>& segment,
		  const CGAL::Bbox_3& bbox)
{
  return typename K::Do_intersect_3()(segment, bbox);
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox,
		  const CGAL::Segment_3<K>& segment)
{
  return typename K::Do_intersect_3()(segment, bbox);
}

} //namespace CGAL

#endif  // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
