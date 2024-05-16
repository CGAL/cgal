// Copyright (c) 2024 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_TETRAHEDRAL_REMESHING_AABB_TREES_H
#define CGAL_TETRAHEDRAL_REMESHING_AABB_TREES_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>

#include <vector>
#include <optional>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  template<typename Triangulation>
  class AABB_triangle_tree
  {
    using Tr = Triangulation;
    using Gt = typename Tr::Geom_traits;
    using FT = typename Gt::FT;
    using Point_3 = typename Gt::Point_3;
    using Ray_3 = typename Gt::Ray_3;

    using Triangle_vec = std::vector<typename Tr::Triangle>;
    using Triangle_iter = typename Triangle_vec::iterator;
    using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
    using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
    using AABB_tree = CGAL::AABB_tree<AABB_triangle_traits>;

    Triangle_vec m_aabb_triangles;
    AABB_tree m_triangles_aabb_tree;
    FT m_aabb_epsilon;

  public:
    template<typename C3t3>
    void build_from_c3t3(const C3t3& c3t3)
    {
      for (const auto& f : c3t3.facets_in_complex())
      {
        m_aabb_triangles.push_back(c3t3.triangulation().triangle(f));
      }
      m_triangles_aabb_tree.rebuild(m_aabb_triangles.begin(), m_aabb_triangles.end());
      m_triangles_aabb_tree.accelerate_distance_queries();

      // compute epsilon for AABB tree of facets
      const CGAL::Bbox_3& bb = m_triangles_aabb_tree.bbox();
      m_aabb_epsilon = CGAL::square(1e-2 * (std::min)(bb.xmax() - bb.xmin(),
                                           (std::min)(bb.ymax() - bb.ymin(),
                                                      bb.zmax() - bb.zmin())));
    }

    Point_3 closest_point(const Point_3& p) const
    {
      return m_triangles_aabb_tree.closest_point(p);
    }

    Point_3 project(const Point_3& p, const Ray_3& ray) const
    {
      if (m_triangles_aabb_tree.squared_distance(p) < m_aabb_epsilon)
        return m_triangles_aabb_tree.closest_point(p);

      using Projection = std::optional<
        typename AABB_tree::template Intersection_and_primitive_id<Ray_3>::Type>;

      auto get_intersection_point =
        [](const Projection& proj) -> std::optional<Point_3>
        {
          const auto intersection = proj.value().first;
          if (const Point_3* pt = std::get_if<Point_3>(&intersection))
            return *pt;
          else
            return std::nullopt;
        };

      // this lambda is called only when we are sure that proj is a Segment
      auto get_intersection_midpoint =
        [](const Projection& proj) -> std::optional<Point_3>
        {
          const auto intersection = proj.value().first;
          using Segment = typename Tr::Geom_traits::Segment_3;
          if (const Segment* s = std::get_if<Segment>(&intersection))
            return CGAL::midpoint(s->source(), s->target());
          else
          {
            CGAL_assertion(false);
            return std::nullopt;
          }
        };

      const Projection proj = m_triangles_aabb_tree.first_intersection(ray);
      const Projection proj_opp = m_triangles_aabb_tree.first_intersection(
        Gt().construct_opposite_ray_3_object()(ray));

      if (proj == std::nullopt && proj_opp == std::nullopt)
      {
        return p;
      }
      else if (proj_opp == std::nullopt)
      {
        const auto pt = get_intersection_point(proj);
        if (pt != std::nullopt)
          return pt.value();
        else
          return get_intersection_midpoint(proj).value();
      }
      else if (proj == std::nullopt)
      {
        const auto pt = get_intersection_point(proj_opp);
        if (pt != std::nullopt)
          return pt.value();
        else
          return get_intersection_midpoint(proj_opp).value();
      }
      else //both intersections are valid
      {
        const auto op1 = get_intersection_point(proj);
        const auto op2 = get_intersection_point(proj_opp);

        const FT sqd1 = (op1 == std::nullopt) ? 0.
                        : CGAL::squared_distance(p, op1.value());
        const FT sqd2 = (op2 == std::nullopt) ? 0.
                        : CGAL::squared_distance(p, op2.value());

        if (sqd1 != 0. && sqd1 < sqd2)
          return op1.value();
        else if (sqd2 != 0)
          return op2.value();
        else
          return p;
      }
    }
  };

  template<typename Triangulation>
  class AABB_segment_tree
  {
    using Tr = Triangulation;
    using Gt = typename Tr::Geom_traits;
    using Point_3 = typename Gt::Point_3;

    using Segment_vec = std::vector<typename Gt::Segment_3>;
    using Segment_iter = typename Segment_vec::iterator;
    using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
    using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
    using AABB_tree = CGAL::AABB_tree<AABB_segment_traits>;

    Segment_vec m_aabb_segments;
    AABB_tree m_segments_aabb_tree;

  public:
    // build AABB tree of edges in complex
    template<typename C3t3>
    void build_from_c3t3(const C3t3& c3t3)
    {
      for (const auto& e : c3t3.edges_in_complex())
      {
        m_aabb_segments.push_back(c3t3.triangulation().segment(e));
      }
      m_segments_aabb_tree.rebuild(m_aabb_segments.begin(), m_aabb_segments.end());
      m_segments_aabb_tree.accelerate_distance_queries();
    }

    Point_3 closest_point(const Point_3& p) const
    {
      return m_segments_aabb_tree.closest_point(p);
    }
  };

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_TETRAHEDRAL_REMESHING_AABB_TREES_H
