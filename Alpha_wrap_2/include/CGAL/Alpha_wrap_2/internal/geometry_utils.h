// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_GEOMETRY_UTILS_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_GEOMETRY_UTILS_H

#include <CGAL/license/Alpha_wrap_2.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename Tr>
bool
less_squared_radius_of_min_empty_circle(typename Tr::Geom_traits::FT sq_alpha,
                                        const typename Tr::Edge& e,
                                        const Tr& tr)
{
  using Face_handle = typename Tr::Face_handle;

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
  std::cout << "Checking for traversability of edge" << std::endl;
#endif

  auto angle = tr.geom_traits().angle_2_object();

  const Face_handle f = e.first;
  const int ic = e.second;
  const Face_handle n = f->neighbor(ic);

  decltype(auto) p1 = tr.point(f, Tr::ccw(ic));
  decltype(auto) p2 = tr.point(f, Tr::cw(ic));

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
  std::cout << "diametrical squared radius: " << squared_radius(p1, p2) << std::endl;
#endif

  // This is not actually possible in the context of alpha wrapping, but keeping it for genericity
  // and because it does not cost anything.
  if(tr.is_infinite(n))
  {
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cerr << "Warning: computing less_squared_radius_of_min_empty_circle() with an infinite neighbor?" << std::endl;
#endif
    CGAL_assertion(!tr.is_infinite(f));

    // In 2D, we can use the angle at the opposite corner to determine
    // if the circumcenter is within the circumcircle.
    const Angle ang = angle(p2, tr.point(f, ic), p1);

    if(ang != OBTUSE)
    {
      // the neighbor is infinite, the CC is within the facet, thus the dual edge crosses the gate
      // and the min empty circle is the diametral circle
      const Comparison_result cr = compare_squared_distance(p1, p2, 4 * sq_alpha);
      return (cr == LARGER);
    }
    else
    {
      const Comparison_result cr = compare_squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2),
                                                          sq_alpha);
      return (cr == LARGER);
    }
  }

  if(tr.is_infinite(f))
  {
    const Angle ang = angle(p1, tr.point(n, tr.mirror_index(f, ic)), p2);

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "Face 'f' is infinite; Angle: " << ang << std::endl;
#endif

    if(ang != OBTUSE)
    {
      const Comparison_result cr = compare_squared_distance(p1, p2, 4 * sq_alpha);
      return (cr == LARGER);
    }
    else
    {
      const Comparison_result cr = compare_squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2),
                                                          sq_alpha);
      return (cr == LARGER);
    }
  }

  // both f and n are finite
  const Angle ang_f = angle(p2, tr.point(f, ic), p1);
  const Angle ang_n = angle(p1, tr.point(n, tr.mirror_index(f, ic)), p2);
  CGAL_assertion(ang_f != OBTUSE || ang_n != OBTUSE);

  if(ang_f != OBTUSE && ang_n != OBTUSE)
  {
    Comparison_result cr = compare_squared_distance(p1, p2, 4 * sq_alpha);
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual crosses the face; CR: "
              << typename Tr::Geom_traits().compute_squared_distance_2_object()(p1, p2)
              << " sq alpha " << sq_alpha << std::endl;
#endif
    return cr == LARGER;
  }
  else
  {
    Comparison_result cr = compare_squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2),
                                                  sq_alpha);
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(f): "
              << typename Tr::Geom_traits().compute_squared_radius_2_object()(tr.point(f, 0), tr.point(f, 1),
                                                                              tr.point(f, 2))
              << " sq alpha " << sq_alpha << std::endl;
    std::cout << "cr = " << cr << std::endl; // @tmp
    std::cout << "check = " << CGAL::compare(typename Tr::Geom_traits().compute_squared_radius_2_object()(
                                 tr.point(f, 0), tr.point(f, 1), tr.point(f, 2)), sq_alpha) << std::endl;
#endif

    if(cr != LARGER)
      return false;

    cr = compare_squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2),
                                sq_alpha);
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(n): "
              << typename Tr::Geom_traits().compute_squared_radius_2_object()(
                   tr.point(n, 0), tr.point(n, 1), tr.point(n, 2))
              << " sq alpha " << sq_alpha << std::endl;
    std::cout << "cr = " << cr << std::endl; // @tmp
    std::cout << "check = " << CGAL::compare(typename Tr::Geom_traits().compute_squared_radius_2_object()(
                                 tr.point(n, 0), tr.point(n, 1), tr.point(n, 2)), sq_alpha) << std::endl;
#endif

    return cr == LARGER;
  }
}

template <typename Tr>
typename Tr::Geom_traits::FT
smallest_squared_radius_2(const typename Tr::Edge& e,
                          const Tr& tr)
{
  using Face_handle = typename Tr::Face_handle;
  using FT = typename Tr::Geom_traits::FT;

  auto angle = tr.geom_traits().angle_2_object();
  auto squared_radius = tr.geom_traits().compute_squared_radius_2_object();

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
  std::cout << "Computing circumradius of edge" << std::endl;
#endif

  CGAL_precondition(!tr.is_infinite(e));

  const Face_handle f = e.first;
  const int ic = e.second;
  const Face_handle n = f->neighbor(ic);

  decltype(auto) p1 = tr.point(f, Tr::ccw(ic));
  decltype(auto) p2 = tr.point(f, Tr::cw(ic));

  // This is not actually possible in the context of alpha wrapping, but keeping it for genericity
  // and because it does not cost anything.
  if(tr.is_infinite(n))
  {
    CGAL_assertion(!tr.is_infinite(f));

    const Angle ang = angle(p2, tr.point(f, ic), p1);

    if(ang != OBTUSE)
      return squared_radius(p1, p2);
    else
      return squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2));
  }

  if(tr.is_infinite(f))
  {
    const Angle ang = angle(p1, tr.point(n, tr.mirror_index(f, ic)), p2);

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "Face 'f' is infinite; Angle: " << ang << std::endl;
#endif

    if(ang != OBTUSE)
      return squared_radius(p1, p2);
    else
      return squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2));
  }

  // both f and n are finite
  const Angle ang_f = angle(p2, tr.point(f, ic), p1);
  const Angle ang_n = angle(p1, tr.point(n, tr.mirror_index(f, ic)), p2);
  CGAL_assertion(ang_f != OBTUSE || ang_n != OBTUSE);

  if(ang_f != OBTUSE && ang_n != OBTUSE)
  {
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual crosses the face; CR: " << squared_radius(p1, p2) << std::endl;
#endif

    return squared_radius(p1, p2);
  }
  else
  {
    const FT cr = squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2));
    const FT cnr = squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2));

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(f): " << cr << " CR(n): " << cnr << std::endl;
#endif

    return (CGAL::min)(cr, cnr);
  }
}


} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_GEOMETRY_UTILS_H
