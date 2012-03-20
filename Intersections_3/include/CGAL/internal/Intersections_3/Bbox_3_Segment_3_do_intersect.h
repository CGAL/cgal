// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
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
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez, Stephane Tayeb


#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/Segment_3.h>
#include <CGAL/Bbox_3.h>

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

  template <typename FT,
            bool bounded_0,
            bool bounded_1>
// __attribute__ ((noinline))
  inline
  bool
  do_intersect_bbox_segment_aux(
                          const FT& px, const FT& py, const FT& pz,
                          const FT& qx, const FT& qy, const FT& qz,
                          const double& bxmin, const double& bymin, const double& bzmin,
                          const double& bxmax, const double& bymax, const double& bzmax)
  {
    // The following code encode t1 and t2 by:
    //    t1 = tmin/dmin
    //    t2 = tmax/dmax
    // For the first lines, dmax==dmin and is not explicitly defined.

    // -----------------------------------
    // treat x coord
    // -----------------------------------
    FT dmin, tmin, tmax, dmax;
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

      if(bounded_0 && bxmin < px) // tmin < 0 means px is in the x-range of bbox
      {
        tmin = 0;
        dmin = 1;
      } else {
        tmin = bxmin - px;
        dmin = qx - px;
      }
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

      if(bounded_0 && px < bxmax) // tmin < 0 means px is in the x-range of bbox
      {
        tmin = 0;
        dmin = 1;
      } else {
        tmin = px - bxmax;
        dmin = px - qx;
      }
    }

    // If the query is vertical for x, then check its x-coordinate is in
    // the x-slab.
    if( (px == qx) &&  // <=> (dmin == 0)
        (! (bounded_0 && bounded_1) ) &&
        ( CGAL::sign(tmin) * CGAL::sign(tmax) ) > 0 )
      return false;
    // Note: for a segment the condition sign(tmin)*sign(tmax) > 0 has
    // already been tested by the two previous tests tmax<0 || tmin>dmin
    // (with dmin==0).

    // If dmin == 0, at this point, [t1, t2] == ]-inf, +inf[, or t1 or t2
    // is a NaN. But the case with NaNs is treated as if the interval
    // [t1, t2] was ]-inf, +inf[.

    CGAL_assertion(dmin >= 0);
    CGAL_assertion(dmax >= 0);

    // -----------------------------------
    // treat y coord
    // -----------------------------------
    FT d_, tmin_, tmax_;

    // Say:
    //   tymin = tmin_ / d_
    //   tymax = tmax_ / d_
    if ( qy >= py )
    {
      tmin_ = bymin - py;
      tmax_ = bymax - py;
      d_ = qy - py;
    }
    else
    {
      tmin_ = py - bymax;
      tmax_ = py - bymin;
      d_ = py - qy;
    }

    CGAL_assertion(d_ >= 0);

    // If the segment is vertical for y, then check its y coordinate is in
    // the y-slab.
    if( (py == qy) && // <=> d_ == 0
        ( sign(tmin_) * sign(tmax_) ) > 0 ) return false;

    // If t1 > tymax || tymin > t2, return false.
    if ( dmin > 0 && d_ >  0) {
      if( (dmin*tmax_) < (d_*tmin) ) return false;
      if( (dmax*tmin_) > (d_*tmax) ) return false;
    }

    // If tymin > t1, set t1 = tymin.
    if( dmin == 0 || 
        ( d_ > 0 && (dmin*tmin_) > (d_*tmin) ) )
    {
      tmin = tmin_;
      dmin = d_;
      if(bounded_1 && tmin > dmin) return false; // if t1 > 1, for a segment
      if(bounded_0 && tmin < FT(0)) {
        // set t1=max(t1, 0), for a ray or a segment
        tmin = FT(0);
        dmin = FT(1);
      }
    }

    // If tymax < t2, set t2 = tymax.
    if( dmax == 0 || 
        ( d_ > 0 && (dmax*tmax_) < (d_*tmax) ) )
    {
      tmax = tmax_;
      dmax = d_;
      if( bounded_0 && tmax < FT(0)) return false; // if t2 < 0, for a segment or a ray
      if( bounded_1 && tmax > dmax )
      {
        // set t2=min(t2, 1), for a segment
        tmax = FT(1);
        dmax = FT(1);
      }
    }

    CGAL_assertion(dmin >= 0);
    CGAL_assertion(dmax >= 0);

    // -----------------------------------
    // treat z coord
    // -----------------------------------

    // Say:
    //   tzmin = tmin_ / d_
    //   tzmax = tmax_ / d_
    if ( qz >= pz )
    {
      tmin_ = bzmin - pz;
      tmax_ = bzmax - pz;
      d_ = qz - pz;
    }
    else
    {
      tmin_ = pz - bzmax;
      tmax_ = pz - bzmin;
      d_ = pz - qz;
    }

    CGAL_assertion(d_ >= 0);

    // If the query is vertical for z, then check its z-coordinate is in
    // the z-slab.
    if( (pz == qz) && // <=> (d_ == 0)
        ( CGAL::sign(tmin_) * CGAL::sign(tmax_) ) > 0 ) return false;

    // If t1 > tzmax || tzmin > t2, return false.
    if ( dmin > 0 && d_ >  0) {
      if( (dmin*tmax_) < (d_*tmin) ) return false;
      if( (dmax*tmin_) > (d_*tmax) ) return false;
    }

    return ( d_ == 0 ||
             (!bounded_1 || tmin_ <= d_) &&
             (!bounded_0 || tmax_ >= FT(0)) ); // t1 <= 1 && t2 >= 0
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

    return do_intersect_bbox_segment_aux<FT, true, true>(
                          source.x(), source.y(), source.z(),
                          target.x(), target.y(), target.z(),
                          bbox.xmin(), bbox.ymin(), bbox.zmin(),
                          bbox.xmax(), bbox.ymax(), bbox.zmax() );

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
