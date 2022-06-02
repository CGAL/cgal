// Copyright (c) 2019-2022 Google LLC (USA).
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
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_GEOMETRY_UTILS_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_GEOMETRY_UTILS_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

template <typename K>
struct Orientation_of_circumcenter
{
  typedef typename K::Point_3 Point_3;

  typedef Orientation result_type;

  Orientation operator()(const Point_3& p, const Point_3& q, const Point_3& r,
                         const Point_3& ccp, const Point_3& ccq, const Point_3& ccr, const Point_3& ccs) const
  {
    Point_3 cc = circumcenter(ccp, ccq, ccr, ccs);
    return orientation(p, q, r, cc);
  }
};

template <typename Dt>
bool
less_squared_radius_of_min_empty_sphere(typename Dt::Geom_traits::FT sq_alpha,
                                        const typename Dt::Facet& fh,
                                        const Dt& dt)
{
  using Cell_handle = typename Dt::Cell_handle;
  using Point = typename Dt::Point;

  using CK = typename Dt::Geom_traits;
  using Exact_kernel = typename Exact_kernel_selector<CK>::Exact_kernel;
  using Approximate_kernel = Simple_cartesian<Interval_nt_advanced>;
  using C2A = Cartesian_converter<CK, Approximate_kernel>;
  using C2E = typename Exact_kernel_selector<CK>::C2E;

  using Orientation_of_circumcenter = Filtered_predicate<Orientation_of_circumcenter<Exact_kernel>,
                                                         Orientation_of_circumcenter<Approximate_kernel>,
                                                         C2E, C2A>;

  Orientation_of_circumcenter orientation_of_circumcenter;

  const Cell_handle c = fh.first;
  const int ic = fh.second;
  const Cell_handle n = c->neighbor(ic);

  const Point& p1 = dt.point(c, Dt::vertex_triple_index(ic,0));
  const Point& p2 = dt.point(c, Dt::vertex_triple_index(ic,1));
  const Point& p3 = dt.point(c, Dt::vertex_triple_index(ic,2));

  // This is not actually possible in the context of alpha wrapping, but keeping it for genericity
  // and because it does not cost anything.
  if(dt.is_infinite(n))
  {
    Orientation ori = orientation_of_circumcenter(p1, p2, p3,
                                                  dt.point(c, 0), dt.point(c, 1),
                                                  dt.point(c, 2), dt.point(c, 3));

    if(ori == POSITIVE)
    {
      Comparison_result cr = compare_squared_radius(p1, p2, p3, sq_alpha);
      return cr == LARGER;
    }
    else
    {
      Comparison_result cr = compare_squared_radius(dt.point(c, 0), dt.point(c, 1),
                                                    dt.point(c, 2), dt.point(c, 3),
                                                    sq_alpha);
      return cr == LARGER;
    }
  }

  if(dt.is_infinite(c))
  {
    Orientation ori = orientation_of_circumcenter(p1, p2, p3,
                                                  dt.point(n, 0), dt.point(n, 1),
                                                  dt.point(n, 2), dt.point(n, 3));

    if(ori == NEGATIVE)
    {
      Comparison_result cr = compare_squared_radius(p1, p2, p3, sq_alpha);
      return cr == LARGER;
    }
    else
    {
      Comparison_result cr = compare_squared_radius(dt.point(n, 0), dt.point(n, 1),
                                                    dt.point(n, 2), dt.point(n, 3),
                                                    sq_alpha);
      return cr == LARGER;
    }
  }

  // both c and n are finite
  if(orientation_of_circumcenter(p1, p2, p3,
                                 dt.point(c, 0), dt.point(c, 1), dt.point(c, 2), dt.point(c, 3)) !=
     orientation_of_circumcenter(p1, p2, p3,
                                 dt.point(n, 0), dt.point(n, 1), dt.point(n, 2), dt.point(n, 3)))
  {
    Comparison_result cr = compare_squared_radius(p1, p2, p3, sq_alpha);
#ifdef CGAL_AW3_DEBUG_TRAVERSABILITY
    std::cout << "dual crosses the face; CR: "
              << typename Dt::Geom_traits().compute_squared_radius_3_object()(p1, p2, p3)
              << " sq alpha " << sq_alpha << std::endl;
#endif
    return cr == LARGER;
  }
  else
  {
    Comparison_result cr = compare_squared_radius(dt.point(c, 0), dt.point(c, 1),
                                                  dt.point(c, 2), dt.point(c, 3),
                                                  sq_alpha);
#ifdef CGAL_AW3_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(c): "
              << typename Dt::Geom_traits().compute_squared_radius_3_object()(dt.point(c, 0), dt.point(c, 1),
                                                                              dt.point(c, 2), dt.point(c, 3))
              << " sq alpha " << sq_alpha << std::endl;
#endif

    if(cr != LARGER)
      return false;

    cr = compare_squared_radius(dt.point(n, 0), dt.point(n, 1),
                                dt.point(n, 2), dt.point(n, 3),
                                sq_alpha);
#ifdef CGAL_AW3_DEBUG_TRAVERSABILITY
    std::cout << "dual does not cross the face; CR(n): "
              << typename Dt::Geom_traits().compute_squared_radius_3_object()(dt.point(n, 0), dt.point(n, 1),
                                                                              dt.point(n, 2), dt.point(n, 3))
              << " sq alpha " << sq_alpha << std::endl;
#endif

    return cr == LARGER;
  }
}

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_GEOMETRY_UTILS_H
