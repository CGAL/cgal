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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

template <typename K>
struct Orientation_of_circumcenter
{
  typedef typename K::Point_2 Point_2;

  typedef Orientation result_type;

  Orientation operator()(const Point_2& p, const Point_2& q,
                         const Point_2& ccp, const Point_2& ccq, const Point_2& ccr) const
  {
    Point_2 cc = circumcenter(ccp, ccq, ccr);
    return orientation(p, q, cc);
  }
};

template <typename Tr>
bool
less_squared_radius_of_min_empty_circle(typename Tr::Geom_traits::FT sq_alpha,
                                        const typename Tr::Edge& e,
                                        const Tr& tr)
{
  using Face_handle = typename Tr::Face_handle;
  using Point = typename Tr::Point;

  using CK = typename Tr::Geom_traits;
  using Exact_kernel = typename Exact_kernel_selector<CK>::Exact_kernel;
  using Approximate_kernel = Simple_cartesian<Interval_nt_advanced>;
  using C2A = Cartesian_converter<CK, Approximate_kernel>;
  using C2E = typename Exact_kernel_selector<CK>::C2E;

  using Orientation_of_circumcenter = Filtered_predicate<Orientation_of_circumcenter<Exact_kernel>,
                                                         Orientation_of_circumcenter<Approximate_kernel>,
                                                         C2E, C2A>;

  Orientation_of_circumcenter orientation_of_circumcenter;

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
  std::cout << "Checking for traversability of edge" << std::endl;
#endif

  const Face_handle f = e.first;
  const int ic = e.second;
  const Face_handle n = f->neighbor(ic);

  const Point& p1 = tr.point(f, Tr::ccw(ic));
  const Point& p2 = tr.point(f, Tr::cw(ic));

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

    Orientation ori = orientation_of_circumcenter(p1, p2,
                                                  tr.point(f, 0), tr.point(f, 1), tr.point(f, 2));

    if(ori == POSITIVE)
    {
      Comparison_result cr = compare_squared_distance(p1, p2, 4 * sq_alpha);
      return cr == LARGER;
    }
    else
    {
      Comparison_result cr = compare_squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2),
                                                    sq_alpha);
      return cr == LARGER;
    }
  }

  if(tr.is_infinite(f))
  {
    Orientation ori = orientation_of_circumcenter(p1, p2,
                                                  tr.point(n, 0), tr.point(n, 1), tr.point(n, 2));

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "Face 'f' is infinite; Orientation: " << ori << std::endl;
#endif

    if(ori == NEGATIVE)
    {
      Comparison_result cr = compare_squared_distance(p1, p2, 4 * sq_alpha);
      return cr == LARGER;
    }
    else
    {
      Comparison_result cr = compare_squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2),
                                                    sq_alpha);
      return cr == LARGER;
    }
  }

  // both f and n are finite
  if(orientation_of_circumcenter(p1, p2, tr.point(f, 0), tr.point(f, 1), tr.point(f, 2)) !=
     orientation_of_circumcenter(p1, p2, tr.point(n, 0), tr.point(n, 1), tr.point(n, 2)))
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
    std::cout << "check = " << CGAL::compare(typename Tr::Geom_traits().compute_squared_radius_2_object()(tr.point(f, 0), tr.point(f, 1),
                        tr.point(f, 2)), sq_alpha) << std::endl;
#endif

    if(cr != LARGER)
      return false;

    cr = compare_squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2),
                                sq_alpha);
#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(n): "
              << typename Tr::Geom_traits().compute_squared_radius_2_object()(tr.point(n, 0), tr.point(n, 1),
                                                                              tr.point(n, 2))
              << " sq alpha " << sq_alpha << std::endl;
    std::cout << "cr = " << cr << std::endl; // @tmp
    std::cout << "check = " << CGAL::compare(typename Tr::Geom_traits().compute_squared_radius_2_object()(tr.point(n, 0), tr.point(n, 1),
    tr.point(n, 2)), sq_alpha) << std::endl;
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
  using Point = typename Tr::Point;
  using FT = typename Tr::Geom_traits::FT;

  using CK = typename Tr::Geom_traits;
  using Exact_kernel = typename Exact_kernel_selector<CK>::Exact_kernel;
  using Approximate_kernel = Simple_cartesian<Interval_nt_advanced>;
  using C2A = Cartesian_converter<CK, Approximate_kernel>;
  using C2E = typename Exact_kernel_selector<CK>::C2E;

  using Orientation_of_circumcenter = Filtered_predicate<Orientation_of_circumcenter<Exact_kernel>,
                                                         Orientation_of_circumcenter<Approximate_kernel>,
                                                         C2E, C2A>;

  Orientation_of_circumcenter orientation_of_circumcenter;

  auto squared_radius = tr.geom_traits().compute_squared_radius_2_object();

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
  std::cout << "Computing circumradius of edge" << std::endl;
#endif

  CGAL_precondition(!tr.is_infinite(e));

  const Face_handle f = e.first;
  const int ic = e.second;
  const Face_handle n = f->neighbor(ic);

  const Point& p1 = tr.point(f, Tr::ccw(ic));
  const Point& p2 = tr.point(f, Tr::cw(ic));

  // This is not actually possible in the context of alpha wrapping, but keeping it for genericity
  // and because it does not cost anything.
  if(tr.is_infinite(n))
  {
    CGAL_assertion(!tr.is_infinite(f));

    Orientation ori = orientation_of_circumcenter(p1, p2,
                                                  tr.point(f, 0), tr.point(f, 1), tr.point(f, 2));
    if(ori == POSITIVE)
      return squared_radius(p1, p2);
    else
      return squared_radius(tr.point(f, 0), tr.point(f, 1), tr.point(f, 2));
  }

  if(tr.is_infinite(f))
  {
    Orientation ori = orientation_of_circumcenter(p1, p2,
                                                  tr.point(n, 0), tr.point(n, 1), tr.point(n, 2));

#ifdef CGAL_AW2_DEBUG_TRAVERSABILITY
    std::cout << "Face 'f' is infinite; Orientation: " << ori << std::endl;
#endif

    if(ori == NEGATIVE)
      return squared_radius(p1, p2);
    else
      return squared_radius(tr.point(n, 0), tr.point(n, 1), tr.point(n, 2));
  }

  // both f and n are finite
  if(orientation_of_circumcenter(p1, p2, tr.point(f, 0), tr.point(f, 1), tr.point(f, 2)) !=
     orientation_of_circumcenter(p1, p2, tr.point(n, 0), tr.point(n, 1), tr.point(n, 2)))
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
