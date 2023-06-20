// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2011 INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_DO_INTERSECT_H

// inspired from http://cag.csail.mit.edu/~amy/papers/box-jgt.pdf

#include <CGAL/Bbox_3.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <typename LFT, typename BFT>
inline
bool
bbox_line_do_intersect_aux(const LFT px, const LFT py, const LFT pz,
                           const LFT vx, const LFT vy, const LFT vz,
                           const BFT bxmin, const BFT bymin, const BFT bzmin,
                           const BFT bxmax, const BFT bymax, const BFT bzmax)
{
  // use CGAL::compare to access the Coercion_traits between LFT and BFT

  if(compare(px, bxmin) != SMALLER && compare(px, bxmax) != LARGER && // (px >= bxmin) && (px <= bxmax)
     compare(py, bymin) != SMALLER && compare(py, bymax) != LARGER && // (py >= bymin) && (py <= bymax)
     compare(pz, bzmin) != SMALLER && compare(pz, bzmax) != LARGER)   // (pz >= bzmin) && (pz <= bzmax)
  {
    return true;
  }

  typedef typename Coercion_traits<LFT, BFT>::Type FT;
  typename Coercion_traits<LFT, BFT>::Cast to_FT;

  // -----------------------------------
  // treat x coord
  // -----------------------------------
  FT dmin, tmin, tmax;
  if(is_negative(vx))
  {
    tmin = to_FT(px) - to_FT(bxmax);
    tmax = to_FT(px) - to_FT(bxmin);
    dmin = - to_FT(vx);
  }
  else // vx >= LFT(0)
  {
    tmin = to_FT(bxmin) - to_FT(px);
    tmax = to_FT(bxmax) - to_FT(px);
    dmin = to_FT(vx);
  }

  //if px is not in the x-slab
  if(is_zero(dmin) && (is_positive(tmin) || is_negative(tmax)))
    return false;

  FT dmax = dmin;

  // -----------------------------------
  // treat y coord
  // -----------------------------------
  FT d_, tmin_, tmax_;
  if(is_negative(vy))
  {
    tmin_ = to_FT(py) - to_FT(bymax);
    tmax_ = to_FT(py) - to_FT(bymin);
    d_ = - to_FT(vy);
  }
  else // vy >= LFT(0)
  {
    tmin_ = to_FT(bymin) - to_FT(py);
    tmax_ = to_FT(bymax) - to_FT(py);
    d_ = to_FT(vy);
  }

  if(is_zero(d_))
  {
    //if py is not in the y-slab
    if((is_positive(tmin_) || is_negative(tmax_)))
      return false;
  }
  else
  {
    if((dmin*tmax_) < (d_*tmin) || (dmax*tmin_) > (d_*tmax))
      return false;
  }

  if((dmin*tmin_) > (d_*tmin))
  {
    tmin = tmin_;
    dmin = d_;
  }

  if((dmax*tmax_) < (d_*tmax))
  {
    tmax = tmax_;
    dmax = d_;
  }

  // -----------------------------------
  // treat z coord
  // -----------------------------------
  if(is_negative(vz))
  {
    tmin_ = to_FT(pz) - to_FT(bzmax);
    tmax_ = to_FT(pz) - to_FT(bzmin);
    d_ = - to_FT(vz);
  }
  else // vz >= LFT(0)
  {
    tmin_ = to_FT(bzmin) - to_FT(pz);
    tmax_ = to_FT(bzmax) - to_FT(pz);
    d_ = to_FT(vz);
  }

  // if pz is not in the z-slab
  // if(d_ == FT(0) && (tmin_ > FT(0) || tmax_ < FT(0))) return false;
  // The previous line is not needed as either dmin or d_ are not 0
  // (otherwise the direction of the line would be null).
  // The following is equivalent to the in z-slab test if d_=0.

  return ((dmin*tmax_) >= (d_*tmin) && (dmax*tmin_) <= (d_*tmax));
}

template <class K>
bool do_intersect(const typename K::Line_3& line,
                  const CGAL::Bbox_3& bbox,
                  const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  const Point_3& point = line.point();
  const Vector_3& v = line.to_vector();

  return bbox_line_do_intersect_aux(point.x(), point.y(), point.z(),
                                    v.x(), v.y(), v.z(),
                                    bbox.xmin(), bbox.ymin(), bbox.zmin(),
                                    bbox.xmax(), bbox.ymax(), bbox.zmax());
}

template <class K>
bool do_intersect(const CGAL::Bbox_3& bbox,
                  const typename K::Line_3& line,
                  const K& k)
{
  return do_intersect(line, bbox, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_DO_INTERSECT_H
