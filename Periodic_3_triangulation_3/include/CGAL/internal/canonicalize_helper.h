// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

// The purpose of this file is to implement a function 'canonicalize_point()'
// which transforms a point into its canonical representative (the representative
// that belongs to the domain). To improve robustness, the function is allowed
// to modify slightly the position of the point, snapping it to the domain if
// it is epsilon-close (and outside).
//
// Although 'canoncalize_point()' is used by Periodic_3_mesh_3, this file is here
// (in Periodic_3_triangulation_3) because of 'construct_periodic_point()',
// which is a function used in P3T3.h and also needed by 'canonicalize_point()'.
// However, P3M3 needs 'canoncalize_point()' without having access to a triangulation
// and to avoid duplicating it, the function is here.

// Geom_traits must be a model of the concept 'P3T3Traits' for 'construct_periodic_point()'.
// Geom_traits must be a model of the concept 'P3T3RegularTraits' for everything else.

#ifndef CGAL_PERIODIC_3_TRIANGULATION_3_CANONICALIZE_HELPER_H
#define CGAL_PERIODIC_3_TRIANGULATION_3_CANONICALIZE_HELPER_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Periodic_3_triangulation_traits_3.h>

#include <CGAL/Exact_kernel_selector.h>

#include <utility>

namespace CGAL {
namespace P3T3 { // can't name it Periodic_3_triangulation_3 because it's already a class...
namespace internal {

template <typename Gt_>
std::pair<typename Gt_::Point_3, typename Gt_::Periodic_3_offset_3>
construct_periodic_point_exact(const typename Gt_::Point_3& p,
                               const Gt_& gt)
{
  typedef Gt_                                                  Geom_traits;
  typedef typename Geom_traits::Periodic_3_offset_3            Offset;
  typedef typename Geom_traits::Iso_cuboid_3                   Iso_cuboid;

  const Iso_cuboid& domain = gt.get_domain();

  typedef typename Geom_traits::Kernel                         K;
  typedef typename Exact_kernel_selector<K>::Exact_kernel      EK;
  typedef typename Exact_kernel_selector<K>::C2E               C2E;

  C2E to_exact;

  typedef Periodic_3_triangulation_traits_3<EK> Exact_traits;
  Exact_traits etraits(to_exact(domain));

  Offset transl(0, 0, 0);
  typename EK::Point_3 ep = to_exact(p);
  typename EK::Point_3 dp;

  const typename EK::Iso_cuboid_3& exact_domain = etraits.get_domain();

  while(true) /* while not in */
  {
    dp = etraits.construct_point_3_object()(ep, transl);

    if(dp.x() < exact_domain.xmin())
      transl.x() += 1;
    else if(dp.y() < exact_domain.ymin())
      transl.y() += 1;
    else if(dp.z() < exact_domain.zmin())
      transl.z() += 1;
    else if(!(dp.x() < exact_domain.xmax()))
      transl.x() -= 1;
    else if(!(dp.y() < exact_domain.ymax()))
      transl.y() -= 1;
    else if(!(dp.z() < exact_domain.zmax()))
      transl.z() -= 1;
    else
      break;
  }

  return std::make_pair(p, transl);
}

// Given a point `p` in space, compute its offset `o` with respect
// to the canonical instance and returns (p, o)
template <typename Gt_>
std::pair<typename Gt_::Point_3, typename Gt_::Periodic_3_offset_3>
construct_periodic_point(const typename Gt_::Point_3& p, bool& encountered_issue, const Gt_& gt)
{
  typedef Gt_                                                  Geom_traits;
  typedef typename Geom_traits::Point_3                        Point;
  typedef typename Geom_traits::Periodic_3_offset_3            Offset;
  typedef typename Geom_traits::Iso_cuboid_3                   Iso_cuboid;

  const Iso_cuboid& domain = gt.get_domain();

  // Check if p lies within the domain. If not, translate.
  if(!(p.x() < domain.xmin()) && p.x() < domain.xmax() &&
     !(p.y() < domain.ymin()) && p.y() < domain.ymax() &&
     !(p.z() < domain.zmin()) && p.z() < domain.zmax())
    return std::make_pair(p, Offset());

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();

  // Numerical approximations might create inconsistencies between the constructions
  // and the comparisons. For example in a cubic domain of size 2:
  // 1. initial point: P(2+1e-17, 0, 0)
  // 2. the function computes an offset(1, 0, 0),
  // 3. P + (-1, 0, 0) * domain_size constructs Q(-1e-17, 0, 0) // numerical approximation
  // 4. the function computes an offset of (-1, 0, 0)
  // 5. Q + (1, 0, 0) * domain_size constructs (2+1e-17, 0, 0) (that is P)
  // And the function is looping...
  //
  // If this is happening the 'Last_change' enum will break this infinite
  // loop and return the wrong point and the 'encountered_issue' bool will be
  // set to 'true'. An exact version of this function should then be called.

  enum Last_change {
    NO_LAST_CHANGE,
    INCREASED_X, DECREASED_X, INCREASED_Y, DECREASED_Y, INCREASED_Z, DECREASED_Z
  };

  Last_change lc = NO_LAST_CHANGE;
  bool in = false;

  Offset transl(0, 0, 0);
  Point dp;

  while(!in)
  {
    dp = cp(p, transl);

    if(dp.x() < domain.xmin())
    {
      if(lc == DECREASED_X) // stuck in a loop
        break;

      lc = INCREASED_X;
      transl.x() += 1;
    }
    else if(dp.y() < domain.ymin())
    {
      if(lc == DECREASED_Y) // stuck in a loop
        break;

      lc = INCREASED_Y;
      transl.y() += 1;
    }
    else if(dp.z() < domain.zmin())
    {
      if(lc == DECREASED_Z) // stuck in a loop
        break;

      lc = INCREASED_Z;
      transl.z() += 1;
    }
    else if(!(dp.x() < domain.xmax()))
    {
      if(lc == INCREASED_X) // stuck in a loop
        break;

      lc = DECREASED_X;
      transl.x() -= 1;
    }
    else if(!(dp.y() < domain.ymax()))
    {
      if(lc == INCREASED_Y) // stuck in a loop
        break;

      lc = DECREASED_Y;
      transl.y() -= 1;
    }
    else if(!(dp.z() < domain.zmax()))
    {
      if(lc == INCREASED_Z) // stuck in a loop
        break;

      lc = DECREASED_Z;
      transl.z() -= 1;
    }
    else
    {
      in = true;
    }
  }

  std::pair<Point, Offset> pp(p, transl);

  if(dp.x() < domain.xmin() || !(dp.x() < domain.xmax()) ||
     dp.y() < domain.ymin() || !(dp.y() < domain.ymax()) ||
     dp.z() < domain.zmin() || !(dp.z() < domain.zmax()))
  {
    encountered_issue = true;
    pp = construct_periodic_point_exact(p, gt);
  }

  return pp;
}

template <typename Gt_>
bool
is_point_too_close_to_border(const std::pair<typename Gt_::Point_3,
                                             typename Gt_::Periodic_3_offset_3>& pbp,
                             const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::FT                    FT;
  typedef typename Geom_traits::Point_3               Bare_point;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();

  const Bare_point p = cp(pbp.first /*point*/, pbp.second /*offset*/);
  const FT px = p.x();
  const FT py = p.y();
  const FT pz = p.z();

  const Iso_cuboid& domain = gt.get_domain();
  const FT dxm = domain.xmin();
  const FT dym = domain.ymin();
  const FT dzm = domain.zmin();
  const FT dxM = domain.xmax();
  const FT dyM = domain.ymax();
  const FT dzM = domain.zmax();

  // simply comparing to FT::epsilon() is probably not completely satisfactory
  const FT eps = std::numeric_limits<FT>::epsilon();

  FT diff = CGAL::abs(px - dxm);
  if(diff < eps && diff > 0) return true;
  diff = CGAL::abs(px - dxM);
  if(diff < eps && diff > 0) return true;
  diff = CGAL::abs(py - dym);
  if(diff < eps && diff > 0) return true;
  diff = CGAL::abs(py - dyM);
  if(diff < eps && diff > 0) return true;
  diff = CGAL::abs(pz - dzm);
  if(diff < eps && diff > 0) return true;
  diff = CGAL::abs(pz - dzM);
  if(diff < eps && diff > 0) return true;
  return false;
}

template <typename Gt_>
typename Gt_::Point_3
snap_to_domain_border(const typename Gt_::Point_3& p, const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::FT                    FT;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  const FT px = p.x();
  const FT py = p.y();
  const FT pz = p.z();
  FT sx = px, sy = py, sz = pz;

  const Iso_cuboid& domain = gt.get_domain();
  const FT dxm = domain.xmin();
  const FT dym = domain.ymin();
  const FT dzm = domain.zmin();
  const FT dxM = domain.xmax();
  const FT dyM = domain.ymax();
  const FT dzM = domain.zmax();

  // simply comparing to FT::epsilon() is probably not completely satisfactory
  const FT eps = std::numeric_limits<FT>::epsilon();

  if(CGAL::abs(px - dxm) < eps) sx = dxm;
  if(CGAL::abs(px - dxM) < eps) sx = dxM;
  if(CGAL::abs(py - dym) < eps) sy = dym;
  if(CGAL::abs(py - dyM) < eps) sy = dyM;
  if(CGAL::abs(pz - dzm) < eps) sz = dzm;
  if(CGAL::abs(pz - dzM) < eps) sz = dzM;

  return gt.construct_point_3_object()(sx, sy, sz);
}

template <typename Gt_>
typename Gt_::Weighted_point_3
snap_to_domain_border(const typename Gt_::Weighted_point_3& p,
                      const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::Point_3               Bare_point;

  typename Geom_traits::Compute_weight_3 cw = gt.compute_weight_3_object();

  const Bare_point snapped_p = snap_to_domain_border(gt.construct_point_3_object()(p), gt);

  return gt.construct_weighted_point_3_object()(snapped_p, cw(p));
}

/// transform a bare point (living anywhere in space) into the canonical
/// instance of the same bare point that lives inside the base domain
template <typename Gt_>
typename Gt_::Point_3
robust_canonicalize_point(const typename Gt_::Point_3& p, const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::Point_3               Bare_point;
  typedef typename Geom_traits::Periodic_3_offset_3   Offset;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  const Iso_cuboid& domain = gt.get_domain();
  if(p.x() >= domain.xmin() && p.x() < domain.xmax() &&
     p.y() >= domain.ymin() && p.y() < domain.ymax() &&
     p.z() >= domain.zmin() && p.z() < domain.zmax())
    return p;

  bool should_snap = false;
  std::pair<Bare_point, Offset> pbp = construct_periodic_point(p, should_snap, gt);

  if(!should_snap)
  {
    // Even if there is no issue while constructing the canonical point,
    // snap the point if it's too close to a border of the domain
    should_snap = is_point_too_close_to_border(pbp, gt);
  }

  if(should_snap)
  {
    Bare_point sp = snap_to_domain_border(p, gt);

    // might have snapped to a 'max' of the domain, which is not in the domain
    // note: we could snap to 'min' all the time in 'snap_to_domain_border'
    // but this is clearer like that (and costs very little since we should
    // not have to use exact computations too often)
    return robust_canonicalize_point(sp, gt);
  }

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();

  Bare_point canonical_p = cp(pbp.first /*point*/, pbp.second /*offset*/);
  CGAL_postcondition( !(canonical_p.x() < domain.xmin()) &&
                       (canonical_p.x() < domain.xmax()) );
  CGAL_postcondition( !(canonical_p.y() < domain.ymin()) &&
                       (canonical_p.y() < domain.ymax()) );
  CGAL_postcondition( !(canonical_p.z() < domain.zmin()) &&
                       (canonical_p.z() < domain.zmax()) );

  return canonical_p;
}

/// transform a weighted point (living anywhere in space) into the canonical
/// instance of the same weighted point that lives inside the base domain
template <typename Gt_>
typename Gt_::Weighted_point_3
robust_canonicalize_point(const typename Gt_::Weighted_point_3& wp, const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::Point_3               Bare_point;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  const Iso_cuboid& domain = gt.get_domain();
  if(wp.x() >= domain.xmin() && wp.x() < domain.xmax() &&
     wp.y() >= domain.ymin() && wp.y() < domain.ymax() &&
     wp.z() >= domain.zmin() && wp.z() < domain.zmax())
    return wp;

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();
  typename Geom_traits::Compute_weight_3 cw = gt.compute_weight_3_object();
  typename Geom_traits::Construct_weighted_point_3 cwp = gt.construct_weighted_point_3_object();

  const Bare_point& bp = cp(wp);
  Bare_point canonical_point = robust_canonicalize_point(bp, gt);

  return cwp(canonical_point, cw(wp));
}

} // end namespace internal
} // end namespace Periodic_3_triangulation_3
} // end namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_3_CANONICALIZE_HELPER_H
