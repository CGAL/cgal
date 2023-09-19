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

// Given a point `p` in space, compute its offset `o` with respect to the canonical
// domain (i.e., p + o * d is in the canonical domain) and returns `(p, o)`
template <typename Gt_>
std::pair<typename Gt_::Point_3, typename Gt_::Periodic_3_offset_3>
construct_periodic_point(const typename Gt_::Point_3& p,
                         bool& encountered_issue,
                         const Gt_& gt)
{
  typedef Gt_                                                  Geom_traits;
  typedef typename Geom_traits::Point_3                        Point;
  typedef typename Geom_traits::Periodic_3_offset_3            Offset;
  typedef typename Geom_traits::Iso_cuboid_3                   Iso_cuboid;

  const Iso_cuboid& domain = gt.get_domain();

  // Use these rather than Construct_point_3 to avoid construction inaccuracies
  typename Geom_traits::Compare_x_3 cmp_x3 = gt.compare_x_3_object();
  typename Geom_traits::Compare_y_3 cmp_y3 = gt.compare_y_3_object();
  typename Geom_traits::Compare_z_3 cmp_z3 = gt.compare_z_3_object();

  // Check if p lies within the domain. If not, translate.
  if(!(p.x() < domain.xmin()) && p.x() < domain.xmax() &&
     !(p.y() < domain.ymin()) && p.y() < domain.ymax() &&
     !(p.z() < domain.zmin()) && p.z() < domain.zmax())
  {
    return std::make_pair(p, Offset());
  }

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
  // set to 'true'. An exact version of this function is then be called.

  enum Last_change
  {
    NO_LAST_CHANGE,
    INCREASED_X, DECREASED_X, INCREASED_Y, DECREASED_Y, INCREASED_Z, DECREASED_Z
  };

  Last_change lc = NO_LAST_CHANGE;
  bool in = false;

  Offset transl(0, 0, 0);
  const Offset null_off(0, 0, 0);

  Point domain_m(domain.xmin(), domain.ymin(), domain.zmin());
  Point domain_M(domain.xmax(), domain.ymax(), domain.zmax());

  while(!in)
  {
    if(cmp_x3(p, domain_m, transl, null_off) == SMALLER)
    {
      if(lc == DECREASED_X) // stuck in a loop
        break;

      lc = INCREASED_X;
      transl.x() += 1;
    }
    else if(cmp_y3(p, domain_m, transl, null_off) == SMALLER)
    {
      if(lc == DECREASED_Y) // stuck in a loop
        break;

      lc = INCREASED_Y;
      transl.y() += 1;
    }
    else if(cmp_z3(p, domain_m, transl, null_off) == SMALLER)
    {
      if(lc == DECREASED_Z) // stuck in a loop
        break;

      lc = INCREASED_Z;
      transl.z() += 1;
    }
    else if(!(cmp_x3(p, domain_M, transl, null_off) == SMALLER))
    {
      if(lc == INCREASED_X) // stuck in a loop
        break;

      lc = DECREASED_X;
      transl.x() -= 1;
    }
    else if(!(cmp_y3(p, domain_M, transl, null_off) == SMALLER))
    {
      if(lc == INCREASED_Y) // stuck in a loop
        break;

      lc = DECREASED_Y;
      transl.y() -= 1;
    }
    else if(!(cmp_z3(p, domain_M, transl, null_off) == SMALLER))
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

  if(cmp_x3(p, domain_m, transl, null_off) == SMALLER || // < min
     cmp_y3(p, domain_m, transl, null_off) == SMALLER ||
     cmp_z3(p, domain_m, transl, null_off) == SMALLER ||
     !(cmp_x3(p, domain_M, transl, null_off) == SMALLER) || // >= max
     !(cmp_y3(p, domain_M, transl, null_off) == SMALLER) ||
     !(cmp_z3(p, domain_M, transl, null_off) == SMALLER))
  {
    encountered_issue = true;
  }

  return pp;
}

template <typename Gt_>
typename Gt_::Point_3
constrain_to_canonical_domain(const typename Gt_::Point_3& p,
                              const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::FT                    FT;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();

  const Iso_cuboid& domain = gt.get_domain();
  FT x = p.x();
  FT y = p.y();
  FT z = p.z();

  if(p.x() < domain.xmin() || p.x() >= domain.xmax())
    x = domain.xmin();
  if(p.y() < domain.ymin() || p.y() >= domain.ymax())
    y = domain.ymin();
  if(p.z() < domain.zmin() || p.z() >= domain.zmax())
    z = domain.zmin();

  return cp(x, y, z);
}

/// transform a bare point (living anywhere in space) into the canonical
/// instance of the same bare point that lives inside the base domain
template <typename Gt_>
typename Gt_::Point_3
robust_canonicalize_point(const typename Gt_::Point_3& p,
                          const Gt_& gt)
{
  typedef Gt_                                         Geom_traits;
  typedef typename Geom_traits::Point_3               Bare_point;
  typedef typename Geom_traits::Periodic_3_offset_3   Offset;
  typedef typename Geom_traits::Iso_cuboid_3          Iso_cuboid;

  typename Geom_traits::Construct_point_3 cp = gt.construct_point_3_object();

  const Iso_cuboid& domain = gt.get_domain();
  if(p.x() >= domain.xmin() && p.x() < domain.xmax() &&
     p.y() >= domain.ymin() && p.y() < domain.ymax() &&
     p.z() >= domain.zmin() && p.z() < domain.zmax())
    return p;

  bool encountered_issue = false;
  std::pair<Bare_point, Offset> pbp = construct_periodic_point(p, encountered_issue, gt);
  Bare_point canonical_p = cp(pbp.first /*point*/, pbp.second /*offset*/);

  if(encountered_issue)
  {
    // If we encountered an issue, there's no guarantee that the double construction gives a point
    // in the domain (even if we computed it exactly beforehand). So, forcefully put it into the domain.
    canonical_p = constrain_to_canonical_domain(canonical_p, gt);
  }

  CGAL_postcondition( !(canonical_p.x() < domain.xmin()) && (canonical_p.x() < domain.xmax()));
  CGAL_postcondition( !(canonical_p.y() < domain.ymin()) && (canonical_p.y() < domain.ymax()));
  CGAL_postcondition( !(canonical_p.z() < domain.zmin()) && (canonical_p.z() < domain.zmax()));

  return canonical_p;
}

/// transform a weighted point (living anywhere in space) into the canonical
/// instance of the same weighted point that lives inside the base domain
template <typename Gt_>
typename Gt_::Weighted_point_3
robust_canonicalize_point(const typename Gt_::Weighted_point_3& wp,
                          const Gt_& gt)
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
