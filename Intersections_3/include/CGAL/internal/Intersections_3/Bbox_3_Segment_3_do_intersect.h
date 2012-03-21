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
    FT dymin, tymin, tymax, dymax;
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

      if(bounded_0 && bymin < py) // tmin < 0 means py is in the y-range of bbox
      {
        tymin = 0;
        dymin = 1;
      } else {
        tymin = bymin - py;
        dymin = qy - py;
      }
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

      if(bounded_0 && py < bymax) // tmin < 0 means py is in the y-range of bbox
      {
        tymin = 0;
        dymin = 1;
      } else {
        tymin = py - bymax;
        dymin = py - qy;
      }
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
    FT dzmin, tzmin, tzmax, dzmax;
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

      if(bounded_0 && bzmin < pz) // tmin < 0 means pz is in the z-range of bbox
      {
        tzmin = 0;
        dzmin = 1;
      } else {
        tzmin = bzmin - pz;
        dzmin = qz - pz;
      }
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

      if(bounded_0 && pz < bzmax) // tmin < 0 means pz is in the z-range of bbox
      {
        tzmin = 0;
        dzmin = 1;
      } else {
        tzmin = pz - bzmax;
        dzmin = pz - qz;
      }
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

    // If t1 > tymax/dymax || tymin/dymin > t2, return false.
    if( py != qy && px != qx) { // dmin > 0, dymax >0, dmax > 0, dymin > 0
      if( (dymax* tmin) > ( dmin*tymax) ) return false; // TEST TO FILTER
      if( ( dmax*tymin) > (dymin* tmax) ) return false; // TEST TO FILTER
    }

    // If tymin/dymin > t1, set t1 = tymin/dymin.
    if( (px == qx) || // <=> (dmin == 0)
        ( (py != qy) && // <=> (dymin > 0)
          ( dmin*tymin) > (dymin* tmin) ) )  // TEST TO FILTER
    {
      tmin = tymin;
      dmin = dymin;
    }

    // If tymax/dymax < t2, set t2 = tymax/dymax.
    if( (px == qx) || // <=> (dmax > 0)
        ( (py != qy) && // <=> dymax > 0
          (dymax* tmax) > ( dmax*tymax) ) ) // TEST TO FILTER
    {
      tmax = tymax;
      dmax = dymax;
    }

    CGAL_assertion(dmin >= 0);
    CGAL_assertion(dmax >= 0);
    // CGAL_assertion((dmin == 0) == (px == qx && py == qy));
    // CGAL_assertion((dmax == 0) == (px == qx && py == qy));

    // If t1 > tzmax || tzmin > t2, return false.
    if( (px != qx ||
         py != qy ) &&
        (pz != qz) ) // dmin > 0, dmax > 0, dzmax > 0, dzmin > 0
    {
      if( (dzmax* tmin) > ( dmin*tzmax) ) return false; // TEST TO FILTER
      if( ( dmax*tzmin) > (dzmin* tmax) ) return false; // TEST TO FILTER
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
