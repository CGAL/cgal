// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   algo/3d/SelfIntersection.h
 * author Gernot Walzl
 * date   2012-07-18
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_self_intersection.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/vertex_splitters.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/vertex_splitters/Combinatorial_vertex_splitter.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/enum.h>
#include <CGAL/Random.h>

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
# include "ortools/sat/cp_model.h"
# include "ortools/sat/sat_parameters.pb.h"
#endif

#ifdef CGAL_SS3_DUMP_FILES
# include <fstream>
#endif
#include <limits>
#include <list>
#include <random>
#include <unordered_map>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename GeomTraits>
class Polyhedron_perturbation
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Segment_3 = typename GeomTraits::Segment_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Ray_3 = typename GeomTraits::Ray_3;
  using Line_3 = typename GeomTraits::Line_3;
  using Triangle_3 = typename GeomTraits::Triangle_3;
  using Plane_3 = typename GeomTraits::Plane_3;
  using Iso_cuboid_3 = typename GeomTraits::Iso_cuboid_3;

private:
  using Polyhedron = HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;

private:
  using Kernel_wrapper = kernel::Kernel_wrapper<GeomTraits>;
  using Hds_utils = algorithm::Hds_utils<GeomTraits>;
  using Transformation = algorithm::Polyhedron_transformation<GeomTraits>;
  using Self_intersection = algorithm::Self_intersection<GeomTraits>;

private:
  using Mesh = CGAL::Surface_mesh<Point_3>;

private:
  struct Size_shenanigans
  {
    static std::size_t length(const FT& n)
    {
      std::stringstream ss;
      ss << CGAL::exact(n);
      std::string str = ss.str();
      return str.size();
    }

    static std::size_t length(const Point_3& p)
    {
      // returns the maximum number of digits between the x, y and z coordinates
      return (std::max)({length(p.x()), length(p.y()), length(p.z())});
    }

    static std::size_t length(const Plane_3& plane)
    {
      // returns the maximum number of digits between the a, b, c and d coefficients
      return (std::max)({length(plane.a()), length(plane.b()), length(plane.c()), length(plane.d())});
    }
  };

private:
  static std::array<double, 3> rand_vec(double min, double max)
  {
    static std::random_device rd;
    unsigned int s = 0; // rd()
    static std::mt19937 gen(s);
    std::uniform_real_distribution<> rdist(min, max);

    return { rdist(gen), rdist(gen), rdist(gen) };
  }

public:
  /**
    * checks for degeneracies: all pairs of planes should intersect in a non-degenerate line.
    */
  static bool do_all_plane_pairs_intersect(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Check all 2-combinations of planes");
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    typename std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
      FacetSPtr facet1 = *it_f1++;
      typename std::list<FacetSPtr>::iterator it_f2 = it_f1;
      while (it_f2 != polyhedron->facets().end()) {
        FacetSPtr facet2 = *it_f2++;

        // Do not use CGAL::do_intersect: here we want to check that the result is a point
        if (!Kernel_wrapper::intersection(facet1->get_plane(), facet2->get_plane())) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Degenerate facet pair:");
          CGAL_SS3_TRANSF_TRACE_V(1, "  " << facet1->to_string());
          CGAL_SS3_TRANSF_TRACE_V(1, "  " << facet2->to_string());
          result = false;
          break;
        }
      }
      if (!result) {
        break;
      }
    }
    return result;
  }

  /**
    * checks for degeneracies: all triplets of planes should intersect in a point.
    */
  static bool do_all_plane_triplets_intersect(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Check all 3-combinations of planes");

    CGAL_SS3_DEBUG_SPTR(polyhedron);

    bool result = true;
    typename std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
      FacetSPtr facet1 = *it_f1++;
      typename std::list<FacetSPtr>::iterator it_f2 = it_f1;
      while (it_f2 != polyhedron->facets().end()) {
        FacetSPtr facet2 = *it_f2++;
        typename std::list<FacetSPtr>::iterator it_f3 = it_f2;
        while (it_f3 != polyhedron->facets().end()) {
          FacetSPtr facet3 = *it_f3++;

          // Do not use CGAL::do_intersect: here we want to check that the result is a point
          if (!Kernel_wrapper::intersection(facet1->get_plane(), facet2->get_plane(), facet3->get_plane())) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Degenerate facet triplet:");
            CGAL_SS3_TRANSF_TRACE_V(1, "  " << facet1->to_string());
            CGAL_SS3_TRANSF_TRACE_V(1, "  " << facet2->to_string());
            CGAL_SS3_TRANSF_TRACE_V(1, "  " << facet3->to_string());
            result = false;
            break;
          }
        }
        if (!result) {
          break;
        }
      }
      if (!result) {
        break;
      }
    }
    return result;
  }

  /**
   * Checks that the positions of the vertices in two polyhedra are close.
   * Preconditions:
   *   - Both polyhedra are non-null and have the same number of vertices.
   *   - Vertices are assumed to be in the same order.
   * Returns true if all corresponding vertices are within a given epsilon.
   */
  static bool check_perturbed_positions_proximity(const PolyhedronSPtr& poly1,
                                                  const PolyhedronSPtr& poly2,
                                                  double epsilon = 1e-4) // @fixme hardcoded...
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Check vertex promixity");

    CGAL_precondition(poly1 && poly2);
    CGAL_precondition(poly1->vertices().size() == poly2->vertices().size());

    auto it1 = poly1->vertices().begin();
    auto it2 = poly2->vertices().begin();
    for (; it1 != poly1->vertices().end() && it2 != poly2->vertices().end(); ++it1, ++it2) {
      const Point_3& p1 = (*it1)->point();
      const Point_3& p2 = (*it2)->point();
      double dx = CGAL::to_double(p1.x()) - CGAL::to_double(p2.x());
      double dy = CGAL::to_double(p1.y()) - CGAL::to_double(p2.y());
      double dz = CGAL::to_double(p1.z()) - CGAL::to_double(p2.z());
      double dist2 = dx*dx + dy*dy + dz*dz;
      if (dist2 > epsilon * epsilon) {
        CGAL_SS3_TRANSF_TRACE_V(2, "Vertex positions too far: " << (*it1)->id() << " d2=" << dist2 << " epsilon=" << CGAL::square(epsilon));
        return false;
      }
    }
    return true;
  }

  static bool are_all_vertices_degree_3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex->degree() != 3) {
        CGAL_SS3_TRANSF_TRACE_V(32, "High-degree vertex: " << vertex->to_string());
        result = false;
        break;
      }
    }
    return result;
  }

  static bool is_triangle_polyhedron(const PolyhedronSPtr& polyhedron)
  {
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->vertices().size() != 3) {
        return false;
      }
    }
    return true;
  }

  static bool is_stable_vertex(const VertexSPtr& vertex,
                               const PolyhedronSPtr& polyhedron,
                               const FT& sq_tolerance)
  {
    CGAL_SS3_TRANSF_TRACE_V(16, "Check stability of " << vertex->to_string());

    CGAL_SS3_DEBUG_SPTR(polyhedron);

    FT max_sq_displacement = 0;

    for (auto it_wf1 = vertex->facets().begin() ; it_wf1 != vertex->facets().end(); ++it_wf1) {
      if (FacetSPtr f1 = it_wf1->lock()) {
        for (auto it_wf2 = std::next(it_wf1); it_wf2 != vertex->facets().end(); ++it_wf2) {
          if (FacetSPtr f2 = it_wf2->lock()) {
            for (auto it_wf3 = std::next(it_wf2); it_wf3 != vertex->facets().end(); ++it_wf3) {
              if (FacetSPtr f3 = it_wf3->lock()) {
                std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                if (!p_new.has_value()) {
                  CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
                  CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                  continue;
                }

                const FT sqd = CGAL::squared_distance(vertex->point(), p_new.value());
                max_sq_displacement = (std::max)(max_sq_displacement, sqd);
                if (sqd > sq_tolerance) {
                  CGAL_SS3_TRANSF_TRACE_V(1, "  vertex is unstable and would move to " << p_new.value());
                  CGAL_SS3_TRANSF_TRACE_V(1, "    (" << CGAL::approximate_sqrt(CGAL::squared_distance(vertex->point(), p_new.value())) << " away)");
                  CGAL_SS3_TRANSF_TRACE_V(1, "    when recomputed with faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                  return false;
                }
              }
            }
          }
        }
      }
    }

    CGAL_SS3_TRANSF_TRACE_V(16, "Vertex is stable (max displacement = " << CGAL::approximate_sqrt(max_sq_displacement) << ")");
    return true;
  }

  // @fixme this whole approach does not yield a solid perturbed state as the planes
  // go through the points, but the **normalized** planes do not...
  static void rand_move_points(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Moving points randomly...");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // If we are applying a random point perturbation, the mesh must be a triangle mesh.
    // Otherwise, points will no longer be on the supporting planes of their incident facets.
    CGAL_precondition(is_triangle_polyhedron(polyhedron));

    ConfigurationSPtr config = Configuration::get_instance();
    double range = config->get_double("Preprocessing", "perturbation_epsilon");
    CGAL_SS3_TRANSF_TRACE("Points will be moved randomly...");
    CGAL_SS3_TRANSF_TRACE("  perturbation_epsilon = " << range);

    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      const Point_3& p = vertex->point();
      const std::array<double, 3> v_r = rand_vec(-range, range);
      // since it's random, move to doubles to get static filters and avoid DAGs
      FT rx = p.x() + FT(v_r[0]);
      FT ry = p.y() + FT(v_r[1]);
      FT rz = p.z() + FT(v_r[2]);
      Point_3 new_pos{rx, ry, rz};
      CGAL_SS3_TRANSF_TRACE_V(16, "Point from " << vertex->point() << " to " << new_pos);
      vertex->set_point(new_pos);
    }

    // recompute normalized planes to ensure points are on the supporting planes
    polyhedron->init_planes();
    Transformation::normalize_facet_planes(polyhedron);
    CGAL_postcondition(polyhedron && polyhedron->is_consistent());

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/rand_moved.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif
  }

  /**
  * nudges the plane coefficients by a random value in the range [low, high].
  */
  static void perturbPlaneCoefficientsNudge(const FacetSPtr& facet,
                                            const double range)
  {
    CGAL_precondition(Transformation::has_normalized_plane(facet));

    CGAL_SS3_TRANSF_TRACE_V(32, "Nudging (Nudge) Facet " << facet->id());
    CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                        << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    auto nudge = [&](const FT& v)
    {
      static std::random_device rd;
      unsigned int s = 0; // rd()
      // CGAL_SS3_TRANSF_TRACE("seed = " << s);
      static std::mt19937 gen(s);
      static std::uniform_real_distribution<> rdist(-range, range);

      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // the point is also that once 'nv' is a double, its interval will be a singleton,
      // and we will have access to static filters
      double step = rdist(gen);
      double nv = CGAL::to_double(v) + step;
      return nv;
    };

    double na = nudge(facet->get_plane().a());
    double nb = nudge(facet->get_plane().b());
    double nc = nudge(facet->get_plane().c());
    double nd = nudge(facet->get_plane().d()); // @todo do not nudge 'd'? (mind the 'to_double()')

    double n = CGAL::approximate_sqrt(square(na) + square(nb) + square(nc));
    CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

    // below doesn't seem to matter? Probably need specific static filters...
#if 0
    facet->set_plane(Plane_3{na/n, nb/n, nc/n, nd/n});
#else
    // cast to_double() *after* the normalization to have double coordinates in the planes
    // the downside is that we won't have a^2 + b^2 + c^2 == 1,
    // but then again, who does...
    const double a = CGAL::to_double(na/n);
    const double b = CGAL::to_double(nb/n);
    const double c = CGAL::to_double(nc/n);
    const double d = CGAL::to_double(nd/n);
    facet->set_plane(Plane_3{a, b, c, d});
#endif

    CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                      << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    CGAL_postcondition(Transformation::has_normalized_plane(facet));
  }

  // Note that this is meant to be used as V4 perturbation and thus we do NOT reset points
  static void perturbPlanesCoefficientsControlled(const PolyhedronSPtr& polyhedron,
                                                  const double initial_range,
                                                  const double max_displacement,
                                                  const int max_iterations = 50,
                                                  const double reduction_factor = 0.5)
  {
    CGAL::unordered_flat_map<FacetSPtr, Plane_3> original_planes;
    // 'amplitudes' stores the current scale of perturbation for each facet.
    // It starts at 1.0 and shrinks on failure.
    CGAL::unordered_flat_map<FacetSPtr, FT> amplitudes;

    for (const FacetSPtr& f : polyhedron->facets()) {
      original_planes[f] = f->get_plane();
      amplitudes[f] = FT(1);
    }

    CGAL::unordered_flat_map<VertexSPtr, Point_3> original_points;
    CGAL::unordered_flat_map<VertexSPtr, FT> sq_max_displacements;

    // Precompute max displacement per vertex
    const FT global_sq_max = CGAL::square(FT(max_displacement));
    for (const VertexSPtr& v : polyhedron->vertices()) {
      original_points[v] = v->point();

      // Initialize with global max
      FT local_sq_max = global_sq_max;

      // Check neighbors: min of (10% dist) and current max
      for (const EdgeWPtr& we : v->edges()) {
        if (const EdgeSPtr e = we.lock()) {
          VertexSPtr ov = e->other(v);
          const FT sq_dist = CGAL::squared_distance(v->point(), ov->point());
          const FT limit_sq = 0.01 * sq_dist;
          if (limit_sq < local_sq_max) {
            local_sq_max = limit_sq;
          }
        }
      }
      sq_max_displacements[v] = local_sq_max;
    }

    // Random number generation
    static std::random_device rd;
    unsigned int s = 0; // rd();
    std::mt19937 gen(s);
    std::uniform_real_distribution<double> rdist(-1.0, 1.0);

    // Set of facets that need perturbation in the current iteration.
    // Initially, all facets need perturbation.
    std::unordered_set<FacetSPtr> active_facets;
    for (const FacetSPtr& f : polyhedron->facets())
      active_facets.insert(f);

    for (int iter = 0; iter < max_iterations; ++iter)
    {
      CGAL_SS3_TRANSF_TRACE_V(16, "controlled perturbation - iteration #" << iter);

      if (active_facets.empty()) {
         CGAL_SS3_TRANSF_TRACE_V(8, "perturbControlled: converged (no active facets) at iteration #" << iter);
         return;
      }

      // Apply perturbations
      // We only perturb facets that are in the 'active_facets' set.
      // We generate a NEW random nudge every iteration, scaled by the current amplitude.

      for (const FacetSPtr& f : active_facets)
      {
        const Plane_3& orig_pl = original_planes[f];
        const FT amplitude = amplitudes[f];
        const FT scale = FT(initial_range) * amplitude;

        // Generate random nudge in arbitrary precision (FT)
        // We generate in double and convert, which is standard for randoms.
        FT nudge_a = FT(rdist(gen)) * scale;
        FT nudge_b = FT(rdist(gen)) * scale;
        FT nudge_c = FT(rdist(gen)) * scale;
        FT nudge_d = FT(rdist(gen)) * scale;

        FT na = orig_pl.a() + nudge_a;
        FT nb = orig_pl.b() + nudge_b;
        FT nc = orig_pl.c() + nudge_c;
        FT nd = orig_pl.d() + nudge_d;

        FT sqlen = CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc);
        FT len = CGAL::approximate_sqrt(sqlen);
        CGAL_assertion(!is_zero(len));

        CGAL_SS3_TRANSF_TRACE_V(32, "Nudging F" << f->id());
        CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << orig_pl.a() << " " << orig_pl.b() << " "
                                                            << orig_pl.c() << " " << orig_pl.d() << "]");

        f->set_plane(Plane_3(na / len, nb / len, nc / len, nd / len));

        CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                          << f->get_plane().c() << " " << f->get_plane().d() << "]");
      }

      // Check displacement of every vertex ----
      // We check ALL vertices, because a vertex can be affected by a neighbor's perturbation
      // even if its own incident facets were 'inactive'.

      std::unordered_set<FacetSPtr> next_active_facets;

      // Check displacement of a single vertex
      auto is_too_large_displacement = [&](const VertexSPtr& v) {
        const Point_3& p_old = original_points.at(v);

        for (auto it_wf1 = v->facets().begin() ; it_wf1 != v->facets().end(); ++it_wf1) {
          if (FacetSPtr f1 = it_wf1->lock()) {
            for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
              if (FacetSPtr f2 = it_wf2->lock()) {
                for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                  if (FacetSPtr f3 = it_wf3->lock()) {
                    std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                    if (!p_new.has_value()) {
                      CGAL_SS3_TRANSF_TRACE_V(1, "Error: triplet of planes does not define a point!");
                      CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                      std::abort();
                      return true;
                    }

                    if (CGAL::squared_distance(p_old, p_new.value()) > sq_max_displacements[v]) {
                      CGAL_SS3_TRANSF_TRACE_V(16, "  Vertex " << v->id() << " moved too far");
                      CGAL_SS3_TRANSF_TRACE_V(16, "  from " << v->point() << " to " << p_new.value());
                      CGAL_SS3_TRANSF_TRACE_V(16, "  F1 [" << f1->id() << "] " << f1->get_plane().a() << " " << f1->get_plane().b() << " "
                                                                               << f1->get_plane().c() << " " << f1->get_plane().d());
                      CGAL_SS3_TRANSF_TRACE_V(16, "  F2 [" << f2->id() << "] " << f2->get_plane().a() << " " << f2->get_plane().b() << " "
                                                                               << f2->get_plane().c() << " " << f2->get_plane().d());
                      CGAL_SS3_TRANSF_TRACE_V(16, "  F3 [" << f3->id() << "] " << f3->get_plane().a() << " " << f3->get_plane().b() << " "
                                                                               << f3->get_plane().c() << " " << f3->get_plane().d());

                      std::optional<Point_3> p_beg = Kernel_wrapper::intersection(original_planes[f1], original_planes[f2], original_planes[f3]);
                      if (!p_beg.has_value()) {
                        CGAL_SS3_TRANSF_TRACE_V(1, "triplet of ORIGINAL planes does not define a point!");
                      } else {
                        CGAL_SS3_TRANSF_TRACE_V(1, "original intersection is: " << p_beg.value());
                      }

                      next_active_facets.insert(f1);
                      next_active_facets.insert(f2);
                      next_active_facets.insert(f3);
                      return true;
                    }
                  }
                }
              }
            }
          }
        }

        return false;
      };

      bool all_ok = true;

      std::unordered_set<FacetSPtr> overly_perturbed_facets;

      for (const VertexSPtr& v : polyhedron->vertices())
      {
        if (is_too_large_displacement(v)) {
          all_ok = false;
          for (const FacetWPtr& wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
              overly_perturbed_facets.insert(f);
            }
          }
        }
      }

      if (all_ok) {
        CGAL_SS3_TRANSF_TRACE_V(8, "perturbControlled: converged at iteration #" << iter);
        return;
      }

      // Prepare for next iteration: only perturb the facets that caused issues.
      active_facets = std::move(next_active_facets);

      for (FacetSPtr f : active_facets)
        amplitudes[f] *= reduction_factor;
    }

    CGAL_SS3_TRANSF_TRACE_V(1, "Warning: perturbControlled did not converge in " << max_iterations << " iterations");
  }

  /**
  * nudges the plane coefficients but ensure that the perturbed plane goes through 0, 1, or 2 fixed points.
  * If 0 points: nudge all coefficients independently.
  * If 1 point: nudge (a, b, c), recompute d so the plane passes through the point.
  * If 2 points: nudge (a, b, c) with the constraint that the new plane passes through both points.
  */
  static void perturbPlaneCoefficientsFixedPoints(const FacetSPtr& facet,
                                                  const double range,
                                                  const std::vector<const Point_3*>& fixed_points)
  {
    CGAL_precondition(Transformation::has_normalized_plane(facet));
    CGAL_precondition(fixed_points.size() <= 2);

    CGAL_SS3_TRANSF_TRACE_V(16, "Nudging Facet " << facet->id());
    CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                        << facet->get_plane().c() << " " << facet->get_plane().d() << "]");
    CGAL_SS3_TRANSF_TRACE_V(16, "  with " << fixed_points.size() << " fixed points");
    CGAL_SS3_TRANSF_TRACE_CODE(for (const Point_3* fp : fixed_points))
    CGAL_SS3_TRANSF_TRACE_V(16, "    " << *fp);

    static std::random_device rd;
    unsigned int s = 0; // rd()
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(-range, range);

    auto nudge = [&](const FT& v)
    {
      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // the point is also that once 'nv' is a double, its interval will be a singleton,
      // and we will have access to static filters
      double step = rdist(gen);
      double nv = CGAL::to_double(v) + step;
      return nv;
    };

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
    auto nudge_to_simplest_rational_in_interval = [&](const FT& v)
    {
      double d1 = nudge(v);
      double d2 = nudge(v);
      if (d2 < d1) {
        std::swap(d1, d2);
      }
      FT nv = CGAL::simplest_rational_in_interval<CGAL::K::Exact_kernel::FT>(d1, d2);
      return nv;
    };
#endif

    // 0 fixed points: nudge all coefficients independently
    if (fixed_points.size() == 0) {
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      FT na = nudge_to_simplest_rational_in_interval(facet->get_plane().a());
      FT nb = nudge_to_simplest_rational_in_interval(facet->get_plane().b());
      FT nc = nudge_to_simplest_rational_in_interval(facet->get_plane().c());
      FT nd = nudge_to_simplest_rational_in_interval(facet->get_plane().d());
#else
      double na = nudge(facet->get_plane().a());
      double nb = nudge(facet->get_plane().b());
      double nc = nudge(facet->get_plane().c());
      double nd = nudge(facet->get_plane().d());
#endif
      facet->set_plane(Plane_3{na, nb, nc, nd});
    } else if (fixed_points.size() == 1) {
      // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
      const Point_3& p0 = *(fixed_points[0]);

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      FT na = nudge_to_simplest_rational_in_interval(facet->get_plane().a());
      FT nb = nudge_to_simplest_rational_in_interval(facet->get_plane().b());
      FT nc = nudge_to_simplest_rational_in_interval(facet->get_plane().c());
#else
      double na = nudge(facet->get_plane().a());
      double nb = nudge(facet->get_plane().b());
      double nc = nudge(facet->get_plane().c());
#endif
      const FT& x0 = p0.x();
      const FT& y0 = p0.y();
      const FT& z0 = p0.z();
      FT d = - (na * x0 + nb * y0 + nc * z0);
      facet->set_plane(Plane_3{na, nb, nc, d});
      CGAL_postcondition(facet->get_plane().has_on(p0));
    } else if (fixed_points.size() == 2) {
      // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
      const Point_3& p0 = *(fixed_points[0]);
      const Point_3& p1 = *(fixed_points[1]);
      CGAL_assertion(p0 != p1);

      const FT& p0x = p0.x();
      const FT& p0y = p0.y();
      const FT& p0z = p0.z();
      const FT& p1x = p1.x();
      const FT& p1y = p1.y();
      const FT& p1z = p1.z();

      // Direction vector between points
      FT ux = p1x - p0x;
      FT uy = p1y - p0y;
      FT uz = p1z - p0z;
      FT uu = ux*ux + uy*uy + uz*uz;

      // Original normal
      const FT& a0 = facet->get_plane().a();
      const FT& b0 = facet->get_plane().b();
      const FT& c0 = facet->get_plane().c();

      // Project original normal onto plane orthogonal to u
      FT dot = a0*ux + b0*uy + c0*uz;
      FT ab = a0 - dot * ux / uu;
      FT bb = b0 - dot * uy / uu;
      FT cb = c0 - dot * uz / uu;

      // Find a direction to nudge (cross product)
      FT vx = uy * cb - uz * bb;
      FT vy = uz * ab - ux * cb;
      FT vz = ux * bb - uy * ab;

      // Nudge the normal
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      FT epsilon = nudge_to_simplest_rational_in_interval(rdist(gen));
#else
      double epsilon = rdist(gen);
#endif

      FT a1 = ab + epsilon * vx;
      FT b1 = bb + epsilon * vy;
      FT c1 = cb + epsilon * vz;

      // Compute 'd' such that the plane passes through 'p0'
      FT d1 = - (a1 * p0x + b1 * p0y + c1 * p0z);
      facet->set_plane(Plane_3{a1, b1, c1, d1});

      CGAL_postcondition(facet->get_plane().has_on(p0));
      CGAL_postcondition(facet->get_plane().has_on(p1));
    } else {
      CGAL_SS3_TRANSF_TRACE("Error: called fixed point facet perturbation with > 2 fixed points");
    }

    CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                      << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    CGAL_postcondition(Transformation::has_normalized_plane(facet));
  }

  static void perturbPlaneCoefficientsHighDegrees(const FacetSPtr& facet,
                                                  const double range)
  {
    std::vector<const Point_3*> high_degree_points;
    for (const VertexSPtr& v : facet->vertices()) {
      if (v->degree() > 3) {
        high_degree_points.push_back(&(v->point()));
      }
    }
    return perturbPlaneCoefficientsFixedPoints(facet, range, high_degree_points);
  }

  /**
    * checks that all faces have at most two high-degree vertices: a facet with fewer than two high-degree
    * vertices can be perturbed by nudging the high-degree vertices, and pivoting the facet randomly
    * around these fixed points.
    */
  static bool can_trivially_tilt_facets(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->num_high_degree_vertices() > 2) {
        CGAL_SS3_TRANSF_TRACE_V(4, "facet " << facet->id() << " has too many high-degree vertices "
                                             << "(" << facet->num_high_degree_vertices() << ")");
        result = false;
        break;
      }
    }
    return result;
  }

  static void apply_rand_plane_tilts(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    ConfigurationSPtr config = Configuration::get_instance();
    double range = config->get_double("Preprocessing", "perturbation_epsilon");

    // If we only nudged planes with fixed point constraints, we might not ensure generic position,
    // for example if two pairs of constraints are along the same line.
    //
    // @todo could restrict to only high-degree vertices in facets that have 2 high-degree vertices
    for (VertexSPtr vertex : polyhedron->vertices()) {
      const Point_3& p = vertex->point();
      const std::array<double, 3> v_r = rand_vec(-range, range);

      double px = CGAL::to_double(p.x()) + v_r[0];
      double py = CGAL::to_double(p.y()) + v_r[1];
      double pz = CGAL::to_double(p.z()) + v_r[2];

      vertex->set_point(Point_3{px, py, pz});
    }

    for (FacetSPtr facet : polyhedron->facets()) {
      perturbPlaneCoefficientsHighDegrees(facet, range);
    }
  }

  static void apply_rand_plane_tilts_V3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Random Plane Tilt (v3)");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    ConfigurationSPtr config = Configuration::get_instance();
    double range = config->get_double("Preprocessing", "perturbation_epsilon");
    CGAL_SS3_TRANSF_TRACE_V(4, "  perturbation_epsilon = " << range);

    if (can_trivially_tilt_facets(polyhedron)) {
      CGAL_SS3_TRANSF_TRACE_V(4, "Polyhedron can simply be tilted immediately");
      apply_rand_plane_tilts(polyhedron);
      CGAL_assertion_code(bool success =)
      Transformation::reset_points(polyhedron);
      CGAL_assertion(success);
      return;
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/tilt-v3_input.obj", polyhedron, parameters::do_not_triangulate_faces(true).stream_precision(17));
#endif

    CGAL_SS3_TRANSF_TRACE_CODE(unsigned int had_to_triangulate_n = 0;)

    // high-degree vertex --> first 3 incident facets determining the vertex
    std::map<VertexSPtr, std::set<FacetSPtr> > determining_facets;

    // facet --> first 2 determined high-degree vertices
    //
    // The facet becomes fixed at 2 vertices and not 3 vertices despite the vertices being perturbed
    // because if we do 3 random perturbations of vertices, the normal can vary wildly.
    //
    // Ideally, it could be fixed with 3 high-degree vertices and a smarter perturbation (which needs
    // to take into account all incident facets of these 3 fixing vertices...)
    std::map<FacetSPtr, std::set<VertexSPtr> > fixing_vertices;

#ifdef CGAL_SS3_DUMP_FILES
    auto dump_facet = [](const std::string& filename, const FacetSPtr& f)
    {
      using CDT2_Tag = CGAL::No_constraint_intersection_tag; // CGAL::Exact_intersections_tag;
      auto pcdt = Transformation::template construct_facet_triangulation<CDT2_Tag>(f);

      using PCDT = decltype(pcdt);
      using PCDT_VH = typename PCDT::Vertex_handle;
      using PCDT_FH = typename PCDT::Face_handle;

      std::unordered_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(pcdt, in_domain);

      std::map<PCDT_VH, std::size_t> point_to_id;
      std::vector<Point_3> points;
      std::vector<std::vector<std::size_t> > triangles;
      for (PCDT_VH vh : pcdt.finite_vertex_handles()) {
        point_to_id[vh] = points.size();
        points.push_back(vh->point());
      }

      for (PCDT_FH fh : pcdt.finite_face_handles()) {
        if (!get(in_domain, fh)) {
          continue;
        }

        triangles.push_back({point_to_id[fh->vertex(0)],
                             point_to_id[fh->vertex(1)],
                             point_to_id[fh->vertex(2)]});
      }

      CGAL::IO::write_OFF(filename, points, triangles, CGAL::parameters::stream_precision(17));
    };
#endif

    auto is_high_degree = [&](const VertexSPtr& v) -> bool
    {
      return (v->degree() > 3);
    };

    auto has_high_degree_vertices = [](const FacetSPtr& f) -> bool
    {
      for (const VertexSPtr& v : f->vertices()) {
        if (v->degree() > 3) {
          return true;
        }
      }
      return false;
    };

    auto is_vertex_determined = [&](const VertexSPtr& v) -> bool
    {
      CGAL_SS3_TRANSF_TRACE_V(16, "Checking if V" << v->id() << " (deg: " << v->degree() << ") is fixed");
      auto it = determining_facets.find(v);
      return (it != determining_facets.end() && it->second.size() == 3);
    };

    auto is_facet_fixed = [&](const FacetSPtr& f) -> bool
    {
      CGAL_SS3_TRANSF_TRACE_V(16, "Checking if F" << f->id() << " (" << f->vertices().size() << " nv) is fixed");
      CGAL_SS3_TRANSF_TRACE_V(16, "  fixing_vertices size: " << fixing_vertices[f].size());
      CGAL_assertion(fixing_vertices[f].size() <= 3);
      return (f->is_triangle() && fixing_vertices[f].size() == 3) ||
              (!f->is_triangle() && fixing_vertices[f].size() == 2);
    };

#ifdef CGAL_SS3_DUMP_FILES
    unsigned int visited_face_id = 0;
    unsigned int nudged_face_id = 0;
#endif

    // Sort by number of high-degree vertices as to avoid triangulating as much as possible
    auto facet_sorter = [&](const FacetSPtr& a, const FacetSPtr& b)
    {
      auto hdv_count = [&](const FacetSPtr& f) -> unsigned int {
        unsigned int hdv_n = 0;
        for (const VertexSPtr& v : f->vertices()) {
          if (is_high_degree(v)) {
            ++hdv_n;
          }
        }
        return hdv_n;
      };

      // Give priority to facets with no determined vertices.
      // If both or neither have constrained vertices, give priority to the largest hdv count.
      //
      // The point is to avoid cascading exact number types, even if we have to triangulate a little more
      auto get_determined_count = [&](const FacetSPtr& f) -> unsigned int
      {
        unsigned int res = 0;
        for (const VertexSPtr& v : f->vertices()) {
          if (is_vertex_determined(v)) {
            ++res;
          }
        }
        return res;
      };

      unsigned int a_determined_n = get_determined_count(a);
      unsigned int b_determined_n = get_determined_count(b);

      CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->id() << " has " << a_determined_n << " determined vertices");
      CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->id() << " has " << b_determined_n << " determined vertices");

      if (a_determined_n != b_determined_n) {
        // Give priority to the one with the least determined vertices
        return a_determined_n < b_determined_n;
      }

      // same number of determined vertices, give priority to the facet with the most high-degree vertices
      unsigned int a_hdv_n = hdv_count(a);
      unsigned int b_hdv_n = hdv_count(b);

      CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->id() << " has " << a_hdv_n << " high-degree vertices");
      CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->id() << " has " << b_hdv_n << " high-degree vertices");

      if (a_hdv_n != b_hdv_n) {
        // Give priority to the one with the most high-degree vertices
        return a_hdv_n > b_hdv_n;
      }

      // same number of determined vertices and high-degree vertices, give priority to the largest facet
      return a->vertices().size() > b->vertices().size();
    };

    // If the facet has no high-degree vertices, we can just tilt it randomly and it will
    // be fine because by definition all of its vertices are degree 3 and will be stable
    // because an unstable configuration results from almost coplanar facets, which have been
    // merged ahead of randomization.
    // UNLESS we have to triangulate a facet incident to one vertex of this facet without
    // high-degree vertices and then the facet now has a high-degree vertex. If that happens,
    // we want the high-degree vertex to constrain to the (up to 2) facets with no high-degree
    // vertices which we are constraining here ahead of the flooding process.
    // Hence, we mark this as fixed with dummy vertices and add 'v' (a non high-degree vertex)
    // to the 'determining_facets' map.
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (!has_high_degree_vertices(facet)) {
        CGAL_SS3_TRANSF_TRACE_V(32, "Nudge and fix F" << facet->id());
        perturbPlaneCoefficientsNudge(facet, range);

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_low_degree.OFF", facet);
#endif

        // A low degree facet is a constraining place when nudging a vertex incident to it.
        // Use dummy vertices to get that effect.
        for (const VertexSPtr& v : facet->vertices()) {
          fixing_vertices[facet].insert(v);
          if (is_facet_fixed(facet)) {
            break;
          }
        }

        // the point of this is that if the vertex becomes high degree after triangulation,
        // one (or two) facet with low degree vertices will appear in the determining facets
        for (const VertexSPtr& v : facet->vertices()) {
          determining_facets[v].insert(facet);
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is determined by F" << facet->id() << " (a)");
        }
      }
    }

    // This is the main list of facets that we will process
    std::list<FacetSPtr> facets_to_process;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->is_triangle() || !has_high_degree_vertices(facet)) {
        continue;
      }

      facets_to_process.push_back(facet);
    }

    // Some preprocessing: if two faces share more than 2 high-degree vertices, we have to triangulate
    // one of them to ensure generic positioning.
    //
    // We triangulate rather than seemingly smarter method of splitting the facet because
    // perturbing splitted facets is very difficult because they are coplanar and perturbations
    // are thus unstable unless performed around the splitting edge, which adds a ton of constraints
    // and complexity.
    for (;;) // reset every time we triangulate something to avoid needless subdivisions
    {
      bool did_something = false;

      std::list<FacetSPtr> facets_to_exclude;
      for (const FacetSPtr& f : facets_to_process) {
        for (const EdgeSPtr& e : f->edges()) {
          VertexSPtr sv = e->src(f);
          VertexSPtr tv = e->tgt(f);
          if (sv->degree() == 3 || tv->degree() == 3) {
            continue;
          }

          // Find the facets { f' } which appear in both sets of incident facets for the vertices
          std::set<FacetSPtr> sv_facets, tv_facets, common_facets;
          for (FacetWPtr wf : sv->facets()) {
            if (FacetSPtr fptr = wf.lock()) {
              sv_facets.insert(fptr);
            }
          }

          for (FacetWPtr wf : tv->facets()) {
            if (FacetSPtr fptr = wf.lock()) {
              tv_facets.insert(fptr);
            }
          }

          // Find intersection (common facets)
          std::set_intersection(sv_facets.begin(), sv_facets.end(),
                                tv_facets.begin(), tv_facets.end(),
                                std::inserter(common_facets, common_facets.begin())
          );

          for (const FacetSPtr& fprime : common_facets) {
            if (fprime == f) {
              continue;
            }

            bool has_edge = (tv->next(fprime) == sv);
            if (!has_edge) {
              // Mark for triangulation
              facets_to_exclude.push_back(fprime);

              CGAL_SS3_TRANSF_TRACE_V(32, "Facet F" << fprime->id() << " needs triangulating due to missing high-degree edge between V" << sv->id() << " and V" << tv->id());

              CGAL_SS3_TRANSF_TRACE_V(32, "Triangulate F" << fprime->id());
              CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

              Transformation::triangulate_facet(fprime, polyhedron);

              did_something = true;
              break;
            }
          }
          if (did_something) {
            break;
          }
        }
        if (did_something) {
          break; // restart
        }
      }

      for (const FacetSPtr& f : facets_to_exclude) {
        facets_to_process.remove(f);
      }

      if (!did_something) {
        break;
      }
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/tilt-v3_preprocessed.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

    // Forward declarations for mutually recursive lambdas
    std::function<bool(FacetSPtr, VertexSPtr)> add_fixing_vertex;
    std::function<void(VertexSPtr)> determine_vertex;

    auto nudge_constrained_vertex = [&](const VertexSPtr& v)
    {
      CGAL_SS3_TRANSF_TRACE_V(32, "  Nudging V" << v->id() << " from " << v->point());

      std::vector<const Plane_3*> constraining_planes;
      for (const FacetSPtr& df : determining_facets[v]) {
        if (is_facet_fixed(df)) {
          constraining_planes.push_back(&(df->get_plane()));
          CGAL_SS3_TRANSF_TRACE_V(32, "    F" << df->id() << " constrains the nudge");
        }
      }

      CGAL_assertion(constraining_planes.size() <= 3);

      const size_t n_fixed = constraining_planes.size();
      if (n_fixed == 3) {
        Transformation::reset_point(v, { constraining_planes[0],
                                         constraining_planes[1],
                                         constraining_planes[2] });
        CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " reset to " << v->point());
        return;
      }

      const Point_3& p = v->point();

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      This does not look quite ready yet: sometimes, the projections are way off
      and with this, sometimes we produce polyhedra in degenerate positions, probably
      because the smallest rational is the same despite the random intervals.

      // @todo clean of all that stuff + duplicate code in facet.cpp
      static std::random_device rd;
      unsigned int s = 0; // rd()
      // CGAL_SS3_TRANSF_TRACE("seed = " << s);
      static std::mt19937 gen(s);
      static std::uniform_real_distribution<> rdist(-range, range);

      auto nudge = [&](const FT& v)
      {
        // Since we are perturbing, we might as well collapse the DAG of 'v'.
        // the point is also that once 'nv' is a double, its interval will be a singleton,
        // and we will have access to static filters
        double step = rdist(gen);
        double nv = CGAL::to_double(v) + step;
        return nv;
      };

      auto nudge_to_simplest_rational_in_interval = [&](const FT& v)
      {
        double d1 = nudge(v);
        double d2 = nudge(v);
        if (d2 < d1) {
          std::swap(d1, d2);
        }
        FT nv = CGAL::simplest_rational_in_interval<typename GeomTraits::Exact_kernel::FT>(d1, d2);
        return nv;
      };

      FT x = nudge_to_simplest_rational_in_interval(p.x());
      FT y = nudge_to_simplest_rational_in_interval(p.y());
      FT z = nudge_to_simplest_rational_in_interval(p.z());
#else
      std::array<double, 3> v_r = rand_vec(-range/2.0, range/2.0);
      double x = CGAL::to_double(p.x()) + v_r[0];
      double y = CGAL::to_double(p.y()) + v_r[1];
      double z = CGAL::to_double(p.z()) + v_r[2];
#endif
      Point_3 p_nudged { x, y, z };
      CGAL_SS3_TRANSF_TRACE_V(32, "base nudge: " << x << " " << y << " " << z);

      Point_3 p_new;

      if (n_fixed == 0) {
        p_new = p_nudged;
      } else if (n_fixed == 1) {
        const Plane_3& plane = *(constraining_planes[0]);
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        // something similar but a little more subtle:
        // 1. project the point onto the plane
        // 2. express the point as a linear combination of the plane's origin and basis: pp = o + l1 * b1 + l2 * b2
        // 3. nudge l1 and l2 to l1' and l2' with a random interval around l1 and l2, and
        //    simplest_rational_in_interval
        // 4. recompute the point as pp = o + l1' * b1 + l2' * b2
        Point_3 pp = plane.projection(p_nudged);
        const Point_3& o = plane->point();
        const Vector_3& b1 = plane->base1();
        const Vector_3& b2 = plane->base2();
        FT l1 = CGAL::scalar_product(*pp - o, b1);
        FT l2 = CGAL::scalar_product(*pp - o, b2);
        FT nl1 = nudge_to_simplest_rational_in_interval(l1);
        FT nl2 = nudge_to_simplest_rational_in_interval(l2);
        p_new = Point_3 { o.x() + nl1 * b1.x() + nl2 * b2.x(),
                          o.y() + nl1 * b1.y() + nl2 * b2.y(),
                          o.z() + nl1 * b1.z() + nl2 * b2.z() };
#else
        p_new = plane.projection(p_nudged);
#endif
      } else if (n_fixed == 2) {
        const Plane_3& plane1 = *(constraining_planes[0]);
        const Plane_3& plane2 = *(constraining_planes[1]);
        std::optional<Line_3> line = Kernel_wrapper::intersection(plane1, plane2);
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        // something similar but a little more subtle:
        // 1. project the point onto the line
        // 2. express the point as a linear combination of the line's origin and basis: pp = o + l * v
        // 3. nudge l to l' with a random interval around l, and simplest_rational_in_interval
        // 4. recompute the point as pp = o + l' * v
        Point_3 pp = line->projection(p_nudged);
        const Point_3& o = line->point();
        const Vector_3& d = line->to_vector();
        FT l = CGAL::scalar_product(*pp - o, d);
        FT nl = nudge_to_simplest_rational_in_interval(l);
        p_new = Point_3{o.x() + nl * d.x(),
                        o.y() + nl * d.y(),
                        o.z() + nl * d.z()};
#else
        // std::cout << "  Constraint Line" << std::endl;
        // std::cout << "    " << line->point(0) << std::endl;
        // std::cout << "    " << line->point(1) << std::endl;
        p_new = line->projection(p_nudged);
#endif
      }

      CGAL_SS3_TRANSF_TRACE_V(32, "  Nudged V" << v->id() << " to " << p_new);

      v->set_point(p_new);
    };

    // sometimes we could fix as a polygon, but we need triangualte for other reasons
    auto should_triangulate_facet = [&](const FacetSPtr& f) -> bool
    {
      if (f->is_triangle() || is_facet_fixed(f)) {
        return false;
      }

      // force triangulation if the exact stack is getting too deep
      for (const VertexSPtr& v : f->vertices()) {
        // consider only determined or almost-determined vertices
        auto it = determining_facets.find(v);
        if (it == determining_facets.end()) {
          continue;
        }

        std::size_t max_length = 100;

        // if the vertex is determined, it has been recomputed so we can check its length
        if (it->second.size() == 3) {
          std::size_t l = Size_shenanigans::length(v->point());
          if (l > max_length) {
            CGAL_SS3_TRANSF_TRACE_V(32, "Vertex V" << v->id() << " is too long");
            CGAL_SS3_TRANSF_TRACE_V(32, CGAL::exact(v->point()) << " (l=" << l << ")");
            CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " should be triangulated");
            return true;
          }
        }

        // if the vertex will be determined by the fixation of this facet, check the facets length
        if (it->second.size() == 2) {
          for (const FacetSPtr& of : determining_facets[v]) {
            std::size_t l = Size_shenanigans::length(of->get_plane());
            if (l > max_length) {
              CGAL_SS3_TRANSF_TRACE_V(32, "Facet F" << of->id() << " is too long");
              CGAL_SS3_TRANSF_TRACE_V(32, CGAL::exact(of->get_plane()) << " (l=" << l << ")");
              CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " should be triangulated");
              return true;
            }
          }
        }
      }

      return false;
    };

    auto is_facet_overconstrained = [&](const FacetSPtr& f) -> bool
    {
      if (f->is_triangle() || is_facet_fixed(f)) {
        return false;
      }

      // we cannot fix that facet if adding the facet to high-degree vertices
      // would create too many determined vertices (> 2) in any unfixed facet incident to
      // the determined high-degree vertices of this facet
      std::map<FacetSPtr, unsigned int> facets_to_test; // facets --> number of appearances
      for (const VertexSPtr& hdv : f->vertices()) {
        if (is_high_degree(hdv)) {
          for (FacetWPtr inc_f : hdv->facets()) {
            if (FacetSPtr f = inc_f.lock()) {
              if (!is_facet_fixed(f)) {
                ++facets_to_test[f];
              }
            }
          }
        }
      }

      for (const auto& [ft, count] : facets_to_test) {
        // Count the number of high-degree vertices with either:
        // - 3 determining facets
        // - 2 determining facets and incident to 'facet'
        // These are vertices that are determined, or would be determined once we 'add'
        // the facet to its high-degree vertices.
        unsigned int constrain_n = 0;
        for (const VertexSPtr& v : f->vertices()) {
          if (is_high_degree(v)) {
            if (is_vertex_determined(v)) {
              ++constrain_n;
            } else if (determining_facets[v].size() == 2 && ft->has_vertex(v)) {
              ++constrain_n;
            }
          }

          if (constrain_n > 2) {
            CGAL_SS3_TRANSF_TRACE_V(32, "F" << ft->id() << " would be over constrained by fixing of F" << f->id());
            return true;
          }
        }
      }

      return false;
    };

    auto triangulate_facet = [&](const FacetSPtr& facet_tt)
    {
      CGAL_SS3_TRANSF_TRACE_V(32, "Triangulate F" << facet_tt->id());

      CGAL_assertion(!is_facet_fixed(facet_tt));

      // the facet is not yet fixed, so no vertex can have it as determining facet
      CGAL_assertion_code(for (const VertexSPtr& v : facet_tt->vertices()) {)
      CGAL_assertion(determining_facets[v].size() <= 3);
      CGAL_assertion(determining_facets[v].count(facet_tt) == 0);
      CGAL_assertion_code(})

      CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

      auto [local_vertices, new_facets] = Transformation::triangulate_facet(facet_tt, polyhedron);

      for (const VertexSPtr& v : local_vertices) {
        CGAL_SS3_TRANSF_TRACE_V(64, "local vertex " << v->id() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)");

        if (is_vertex_determined(v)) {
          CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is already determined, skipping");
          continue;
        }

        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr fptr = wf.lock()) {
            if (is_facet_fixed(fptr)) {
              CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is determined by F" << fptr->id() << " (c)");
              determining_facets[v].insert(fptr);

              if (is_vertex_determined(v)) {
                determine_vertex(v);
                break;
              }
            }
          }
        }
      }

      // already-determined vertices are fixed points for the new facets
      for (const FacetSPtr& nf : new_facets) {
        CGAL_SS3_TRANSF_TRACE_V(32, "spawned F" << nf->id());

        for (const VertexSPtr& iv : nf->vertices()) {
          if (is_vertex_determined(iv)) {
            CGAL_SS3_TRANSF_TRACE_V(64, "newborn F" << nf->id() << " is constrained by V" << iv->id());
            fixing_vertices[nf].insert(iv);
          }
        }
      }
    };

    add_fixing_vertex = [&](const FacetSPtr& f, const VertexSPtr& v) -> bool
    {
      CGAL_precondition(fixing_vertices[f].size() <= 3);

      CGAL_SS3_TRANSF_TRACE_V(64, "  Fix F" << f->id() << " with V" << v->id());

      if (is_facet_fixed(f)) {
        CGAL_SS3_TRANSF_TRACE_V(64, "  F" << f->id() << " is already fixed");
        return true;
      }

      fixing_vertices[f].insert(v);

      if (!is_facet_fixed(f)) {
        // nothing to do yet, there are still degrees of freedom in the facet
        return true;
      }

      // Here, the facet has just had enough determined vertices to now be fixed.
      // So, fix it: compute its random perturbation, and add the facet ID to its vertices.

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "F" << f->id() << " is now fixed by");
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : fixing_vertices[f]))
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->id() << " [measure=" << Size_shenanigans::length(fv->point()) << "]");
      CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

      if (f->is_triangle()) {
        CGAL_assertion(fixing_vertices[f].size() == 3); // just to be clear

        // for triangles, all vertices are determined, and there is nothing to nudge
        // (note that vertices were themselves nudged so the facet is nudged).
        f->init_plane();
        Transformation::normalize_plane_coefficients(f);

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_3.OFF", f);
#endif

        // Here we do not need to add the fixed facet to incident determined vertices
        // because all vertices are already fully determined
        return true;
      }

      std::vector<const Point_3*> fixed_points;
      for (const VertexSPtr& fixed_v : fixing_vertices[f]) {
        fixed_points.push_back(&(fixed_v->point()));
      }

      perturbPlaneCoefficientsFixedPoints(f, range, fixed_points);

      CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " is now fixed at " << f->get_plane() << " [measure=" << Size_shenanigans::length(f->get_plane()) << "]");

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_" + std::to_string(fixed_points.size()) + ".OFF", f);
#endif

      CGAL_SS3_TRANSF_TRACE_V(64, "Newly fixed facet F" << f->id() << " determines its high-degree incident vertices...");

      // Need to now tag the vertices of the facet
      for (const VertexSPtr& v : f->vertices()) {
        CGAL_SS3_TRANSF_TRACE_V(64, "incident: " << v->id() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)");
        if (!is_vertex_determined(v)) {
          determining_facets[v].insert(f);
          CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is determined by F" << f->id() << " (b)");
        }
      }

      return true;
    };

    determine_vertex = [&](const VertexSPtr& v)
    {
      CGAL_precondition(is_high_degree(v));
      CGAL_precondition(is_vertex_determined(v));

      CGAL_SS3_TRANSF_TRACE_CODE(auto it = determining_facets[v].begin();)
      CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is now fully determined by"
                                  << " F" << (*it)->id() << " [measure=" << Size_shenanigans::length((*it)->get_plane())
                                  << "] F" << (*std::next(it))->id() << " [measure=" << Size_shenanigans::length((*std::next(it))->get_plane())
                                  << "] F" << (*std::next(it, 2))->id() << " [measure=" << Size_shenanigans::length((*std::next(it, 2))->get_plane()) << "]");

      // set the nudged position for the vertex: a nudge constrained by already fixed incident facets
      nudge_constrained_vertex(v);

      CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is now determined at " << v->point() << " [measure=" << Size_shenanigans::length(v->point()) << "]");

      // compute the plane coefficients of any incident facet that becomes fixed
      // by this vertex becoming determined
      for (FacetWPtr wf : v->facets()) {
        if (FacetSPtr f = wf.lock()) {
          add_fixing_vertex(f, v);
        }
      }
    };

    CGAL_SS3_TRANSF_TRACE_V(8, "== Main facet flood... ==");

    while (!facets_to_process.empty()) {
      facets_to_process.sort(facet_sorter); // @todo priority queue...
      FacetSPtr facet = facets_to_process.front();
      facets_to_process.pop_front();

      CGAL_SS3_TRANSF_TRACE_V(16, "Pop F" << facet->id());

      CGAL_assertion(!facet->is_triangle());
      CGAL_assertion(fixing_vertices[facet].size() <= 2);

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "  Fixing vertices:";)
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : fixing_vertices[facet]))
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->id();)
      CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/visited_face_" + std::to_string(visited_face_id++) + ".OFF", facet);
#endif

      CGAL_assertion(facet->vertices().size() >= 3);

      if (is_facet_overconstrained(facet) || should_triangulate_facet(facet)) {
        triangulate_facet(facet);
        continue;
      }

      // Now, adding the facet to the high-degree vertices will not over constrain the facet, so do it:
      for (const VertexSPtr& v : facet->vertices()) {
        if (!is_vertex_determined(v)) {
          determining_facets[v].insert(facet);
          CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is determined by F" << facet->id() << " (d)");
          if (is_high_degree(v) && is_vertex_determined(v)) {
            // When the vertex becomes fixed (its 3 determining facets become known), we need:
            // - to perturb the position of the vertex
            // - to update all incident facets to check if they are now fixed and in that case,
            //   compute their plane coefficients
            determine_vertex(v);
          }
        }
      }
    }

    // Some facets might have high-degree vertices, but still some freedom of movement
    // after the flooding, fix them
    //
    // @todo could we not simply nudge high-degree vertices and fix everything left (triangle or not)?
    CGAL_SS3_TRANSF_TRACE_V(16, "== Deal with remaining facets with high degree vertices... ==");

    for (const FacetSPtr& f : polyhedron->facets()) {
      if (f->is_triangle() || is_facet_fixed(f)) {
        continue;
      }

      CGAL_SS3_TRANSF_TRACE_V(32, "Nudge and fix F" << f->id() << " [remaining]");

      std::vector<const Point_3*> fixed_points;
      for (const VertexSPtr& v : fixing_vertices[f]) {
        CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is a fixing vertex");
        fixed_points.push_back(&(v->point()));
      }

      perturbPlaneCoefficientsFixedPoints(f, range, fixed_points);

      // fixing the facet cannot determine a vertex because the facet has already been visited

      // add random vertices to mark the facet as fixed
      static int dummy_id = -1;
      while (!is_facet_fixed(f)) {
        VertexSPtr dummy_v = Vertex::create(CGAL::ORIGIN);
        dummy_v->set_id(dummy_id--);
        fixing_vertices[f].insert(dummy_v);
      }

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_remaining.OFF", f);
#endif

      CGAL_postcondition(is_facet_fixed(f));
    }

    // At this point, everything that is not a high-degree triangular facet should be fixed
    for (const FacetSPtr& f : polyhedron->facets()) {
      if (f->is_triangle() || !has_high_degree_vertices(f)) {
        continue;
      }
      CGAL_assertion(is_facet_fixed(f));
    }

    // Nudge vertices that can still be nudged, for randomness
    CGAL_SS3_TRANSF_TRACE_V(16, "== Nudge undetermined high-degree vertices... ==");

    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (is_high_degree(v) && !is_vertex_determined(v)) {
        CGAL_SS3_TRANSF_TRACE_V(32, "  V" << v->id() << " is high degree and not fully determined, nudge it");
        nudge_constrained_vertex(v);

        // determine the vertex
        // since we know only triangle facets are left, we don't need to cascade and check
        // if incident facets become fixed
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (!is_facet_fixed(f)) {
              // the facet cannot be without high-degree vertices since v is high degree
              CGAL_assertion(f->is_triangle());
              fixing_vertices[f].insert(v);
            }
          }
        }

        // add dummy facets to mark the vertex as determined
        static int dummy_id = -1;
        while (!is_vertex_determined(v)) {
          FacetSPtr dummy_f = Facet::create();
          dummy_f->set_id(dummy_id--);
          determining_facets[v].insert(dummy_f);
        }

        CGAL_postcondition(is_vertex_determined(v));
      }
    }

    // Now handle triangle faces with high degrees
    CGAL_SS3_TRANSF_TRACE_V(16, "== Deal with remaining triangles... ==");

    for (const FacetSPtr& f : polyhedron->facets()) {
      if (!f->is_triangle() || !has_high_degree_vertices(f)) {
        continue;
      }

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "Fix triangle F" << f->id() << " [");
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : f->vertices()) {)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << v->id());
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " (" << v->degree() << ")");
      CGAL_SS3_TRANSF_TRACE_CODE(if (is_vertex_determined(v)) { ss << "*"; })
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " "; } ss << "]";)
      CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

      CGAL_SS3_TRANSF_TRACE_V(32, "Nudge and fix F" << f->id() << " [triangle]");

      std::vector<const Point_3*> fixed_points;
      for (const VertexSPtr& v : fixing_vertices[f]) {
        CGAL_SS3_TRANSF_TRACE_V(32, "  V" << v->id() << " is a fixing vertex");
        fixed_points.push_back(&(v->point()));
      }

      if (fixed_points.size() == 3) {
        f->init_plane();
        Transformation::normalize_plane_coefficients(f);
      } else {
        perturbPlaneCoefficientsFixedPoints(f, range, fixed_points);
      }

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_triangle.OFF", f);
#endif

      // We still need to update the determining facets because some neighboring
      // facets could be an unfixed high-degree triangle
      for (const VertexSPtr& v : f->vertices()) {
        if (!is_vertex_determined(v)) {
          determining_facets[v].insert(f);
          CGAL_SS3_TRANSF_TRACE_V(32, "  V" << v->id() << " is determined by F" << f->id() << " (f)");
          // no need to cascade here, we know only triangles are left
        }
      }

      for (const VertexSPtr& v : f->vertices()) {
        fixing_vertices[f].insert(v);
      }

      CGAL_postcondition(is_facet_fixed(f));
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/tilt_v3-pre_reset.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

    CGAL_SS3_TRANSF_TRACE_V(16, "Reset the position of not-fully-constrained vertices...");

    // Recompute all points which were not fixed (degree 3 vertices)
    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (!is_high_degree(v)) {
        Transformation::reset_point(v);

        // determine (without cascading)
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            determining_facets[v].insert(f);
          }
        }
      } else {
        CGAL_assertion(is_vertex_determined(v)); // high degree vertices have already been determined
      }
    }

    CGAL_SS3_TRANSF_TRACE_V(8, "All facets processed");

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/tilt_v3.obj", polyhedron, parameters::do_not_triangulate_faces(true));
    IO::write_OBJ("results/tilt_v3-triangulated.obj", polyhedron, parameters::do_not_triangulate_faces(false));
#endif

    CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
    CGAL_assertion(is_vertex_determined(v));
    CGAL_assertion_code(})

    CGAL_assertion_code(for (const FacetSPtr& f : polyhedron->facets()) {)
    CGAL_assertion(fixing_vertices[f].size() <= 3);
    CGAL_assertion_code(})

    CGAL_assertion_code(for (const FacetSPtr& facet : polyhedron->facets()) {)
    CGAL_assertion_code(for (const VertexSPtr& v : facet->vertices()) {)
    CGAL_assertion(facet->get_plane().has_on(v->point()));
    CGAL_assertion_code(})
    CGAL_assertion_code(})

    CGAL_SS3_TRANSF_TRACE_V(8, "Had to triangulate " << had_to_triangulate_n << " facets");

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " has depth " << CGAL::depth(v->point()));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " has depth " << CGAL::depth(f->get_plane()));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " has length " << Size_shenanigans::length(v->point()));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " has length " << Size_shenanigans::length(f->get_plane()));
  }

#if 1
  // duplicated because all faces have their source here
  static void get_clipped_plane_faces(const VertexSPtr vertex,
                                      const Iso_cuboid_3& bbox,
                                      std::vector<Point_3>& points,
                                      std::vector<std::vector<std::size_t> >& triangles,
                                      std::vector<FacetSPtr>& triangle_to_facet)
  {
    using Vector_3 = typename GeomTraits::Vector_3;

    for(const auto& facet_wptr : vertex->facets())
    {
      if(FacetSPtr facet = facet_wptr.lock())
      {
        const Plane_3& plane = facet->get_plane();

        std::vector<Point_3> local_range;
        auto res = CGAL::intersection(bbox, plane);
        if (!res) {
          // Should not happen
          CGAL_SS3_TRANSF_TRACE_V(1, "no intersection between plane and bbox?!");
          CGAL_assertion(false);
          std::exit(1);
        } else if (const Triangle_3* itr = std::get_if<Triangle_3>(&*res)) {
          for (int i=0; i<3; ++i) {
            local_range.push_back((*itr)[i]);
          }
        } else if (const std::vector<Point_3>* ir = std::get_if<std::vector<Point_3> >(&*res)) {
          for (const Point_3& p : *ir) {
            local_range.push_back(p);
          }
        } else {
          CGAL_SS3_TRANSF_TRACE_V(1, "plane/bbox intersection is not a polygon");
          CGAL_assertion(false);
          std::exit(1);
        }

        // Ensure orientation: normal of local_range must match plane's orientation
        if(local_range.size() >= 3) {
          Vector_3 plane_normal = plane.orthogonal_vector();
          Vector_3 tri_normal = CGAL::cross_product(local_range[1] - local_range[0], local_range[2] - local_range[1]);
          if(tri_normal * plane_normal < 0) {
            std::reverse(local_range.begin(), local_range.end());
          }
        }

        // Triangulate by fanning from the first point
        std::size_t base_idx = points.size();
        for(const Point_3& p : local_range) {
          points.push_back(p);
        }

        if(local_range.size() >= 3) {
          // Build a single polygon for triangulation
          std::vector<std::vector<std::size_t> > polygons(1);
          for(std::size_t i = 0; i < local_range.size(); ++i)
            polygons.back().push_back(base_idx + i);

          CGAL::Polygon_mesh_processing::triangulate_polygons(points, polygons);

          for(const auto& tri : polygons) {
            CGAL_assertion(tri.size() == 3);
            triangles.push_back(tri);
            triangle_to_facet.push_back(facet);
          }
        }
      }
    }
  }
#else
  static void get_clipped_plane_faces(const VertexSPtr vertex,
                                      const Iso_cuboid_3& bbox,
                                      std::vector<Point_3>& points,
                                      std::vector<std::vector<std::size_t> >& triangles,
                                      std::vector<FacetSPtr>& triangle_2_sptr)
  {
    using Vector_3 = typename GeomTraits::Vector_3;
    const Point_3& center = vertex->point();

    for(const auto& facet_wptr : vertex->facets())
    {
      if(FacetSPtr facet = facet_wptr.lock())
      {
        const Plane_3& plane = facet->get_plane();

        std::vector<Point_3> local_range;
        auto res = CGAL::intersection(bbox, plane);
        if (!res) {
          // Should not happen, as bbox is constructed to contain all intersections
          CGAL_SS3_TRANSF_TRACE_V(1, "no intersection between plane and bbox");
          CGAL_assertion(false);
          std::exit(1);
        } else if (const Triangle_3* itr = std::get_if<Triangle_3>(&*res)) {
          for (int i=0; i<3; ++i) {
            local_range.push_back((*itr)[i]);
          }
        } else if (const std::vector<Point_3>* ir = std::get_if<std::vector<Point_3> >(&*res)) {
          for (const Point_3& p : *ir) {
            local_range.push_back(p);
          }
        } else {
          CGAL_SS3_TRANSF_TRACE_V(1, "plane/bbox intersection is not a polygon");
          CGAL_assertion(false);
          std::exit(1);
        }

        // Ensure orientation: normal of local_range must match plane's orientation
        if(local_range.size() >= 3) {
          Vector_3 plane_normal = plane.orthogonal_vector();
          Vector_3 tri_normal = CGAL::cross_product(local_range[1] - local_range[0], local_range[2] - local_range[1]);
          if(tri_normal * plane_normal < 0) {
            std::reverse(local_range.begin(), local_range.end());
          }
        }

        // Insert extremities (projections) into local_range at the correct place
        const Point_3& prev_pt = vertex->prev(facet)->point();
        const Point_3& v_pt = vertex->point();
        const Point_3& next_pt = vertex->next(facet)->point();

        Vector_3 prev_dir = prev_pt - v_pt;
        Vector_3 next_dir = next_pt - v_pt;
        auto res_prev = CGAL::intersection(Ray_3{v_pt, prev_dir}, bbox);
        auto res_next = CGAL::intersection(Ray_3{v_pt, next_dir}, bbox);
        CGAL_assertion(res_prev && res_next);

        auto get_ray_bbox_extremity = [&](const auto& res, const Point_3& src) -> std::optional<Point_3> {
          if(const Segment_3* seg = std::get_if<Segment_3>(&*res)) {
            if(seg->source() == src)
              return seg->target();
            else if(seg->target() == src)
              return seg->source();
            else
              return std::nullopt;
          }
          return std::nullopt;
        };

        std::optional<Point_3> opt_prev = get_ray_bbox_extremity(res_prev, v_pt);
        std::optional<Point_3> opt_next = get_ray_bbox_extremity(res_next, v_pt);
        CGAL_assertion(opt_prev && opt_next);

        // check linearly to find in which segment the point belongs
        // (it belongs by construction)
        auto insert_in_order = [&local_range](const Point_3& new_p)
        {
          for(auto it = local_range.begin(); it != local_range.end(); ++it) {
            auto it_next = std::next(it);
            if(it_next == local_range.end()) {
              it_next = local_range.begin();
            }
            if(new_p != *it && new_p != *it_next &&
               CGAL::collinear(*it, new_p, *it_next) &&
               CGAL::collinear_are_strictly_ordered_along_line(*it, new_p, *it_next))
            {
              // std::cout << "insert " << new_p << " between " << *it << " and " << *it_next << std::endl;
              local_range.insert(it_next, new_p);
              break;
            }
          }
        };

        insert_in_order(*opt_prev);
        insert_in_order(*opt_next);

        // Triangulate by fanning from the center vertex (vertex->point())
        std::size_t center_idx = points.size();
        points.push_back(center);

        std::size_t base_idx = points.size();
        for(const Point_3& p : local_range) {
          points.push_back(p);
        }

        // --- Sector logic ---
        // Get prev/next points for this facet at the center vertex
        Vector_3 n = plane.orthogonal_vector();

        // Orientation planes: through center, normal is n x (prev-center) and n x (next-center)
        Vector_3 v_prev = prev_pt - center;
        Vector_3 v_next = next_pt - center;
        Vector_3 n_prev = CGAL::cross_product(n, v_prev);
        Vector_3 n_next = CGAL::cross_product(n, v_next);
        Plane_3 plane_prev(center, n_prev);
        Plane_3 plane_next(center, n_next);

        // Determine if angle at center is > 180°
        // If next is to the left of prev (in the facet's orientation), angle < 180°
        // If next is to the right of prev, angle > 180°
        // Use orientation of (center, prev, next) with normal n
        bool angle_gt_180 = (CGAL::orientation(center, next_pt, prev_pt, center + n) == CGAL::NEGATIVE);
        // std::cout << "center = " << center << std::endl;
        // std::cout << "prev_pt = " << prev_pt << std::endl;
        // std::cout << "next_pt = " << next_pt << std::endl;
        // std::cout << "center + n = " << center + n << std::endl;
        // std::cout << "angle_gt_180 = " << angle_gt_180 << std::endl;

        for(std::size_t i=0; i<local_range.size(); ++i) {
          std::size_t i1 = base_idx + i;
          std::size_t i2 = base_idx + ((i+1)%local_range.size());

          // Compute midpoint of the two extremities
          const Point_3& p1 = points[i1];
          const Point_3& p2 = points[i2];
          Point_3 mid = CGAL::midpoint(p1, p2);
          // std::cout << "test: " << mid << std::endl;

          // Test if mid is between the two orientation planes
          bool on_pos_side_prev = (plane_prev.oriented_side(mid) == CGAL::ON_NEGATIVE_SIDE);
          bool on_pos_side_next = (plane_next.oriented_side(mid) == CGAL::ON_POSITIVE_SIDE);
          // std::cout << "on_pos_side_prev = " << on_pos_side_prev << std::endl;
          // std::cout << "on_pos_side_next = " << on_pos_side_next << std::endl;

          bool in_sector = (on_pos_side_prev && on_pos_side_next);
          // If angle > 180°, invert logic: triangles inside are NOT part of the facet
          bool is_facet_triangle = angle_gt_180 ? !in_sector : in_sector;

          triangles.push_back({center_idx, i1, i2});
          triangle_2_sptr.push_back(is_facet_triangle ? facet : nullptr);
        }
      }
    }
  }
#endif

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
  static void apply_rand_plane_tilts_V4(const PolyhedronSPtr& polyhedron)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    CGAL_SS3_TRANSF_TRACE_V(4, "Random Plane Tilt (v4)");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    // @todo is that even needed?
    Transformation::normalize_facet_planes(polyhedron);

    // The approach here is to perform plane perturbations, followed by extracting of a solid,
    // yielding a general position polyhedral surface.
    //
    // As long as vertices are stable under the perturbation, surface extraction is reduced
    // to the selection of a local set of facets within a local arrangement around vertices
    // that have a degree > 3.
    //
    // At each input high-degree vertex, we must recover a valid split from the arrangement
    // of the perturbed planes.

    ConfigurationSPtr config = Configuration::get_instance();
    double nudge_range = config->get_double("Preprocessing", "perturbation_epsilon");

    // -- PART 1 --
    // Compute *perturbed* planes such that at a given vertex, any 3-intersection of incident
    // planes intersect closely to the input vertex.

    CGAL_SS3_TRANSF_TRACE_V(8, "Part 1: preprocess planes");

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/V4_input.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

    // Here are the rough steps of the algorithm:
    // - For all vertices, compute the max displacement (as in perturbPlanesCoefficientsControlled())
    // - For all faces, compute a naive normalized perturbed plane (perturbPlaneCoefficientsNudge())
    // - For all vertices, check if the vertex is stable, i.e. any 3-intersection
    //   is close to the vertex base position (is_stable_vertex()).
    // - If a vertex is not stable, mark it as a new "anchor" for its incident faces.
    // - If some vertices are new anchors, fix their position with a random nudge (with rand_vec)
    //   from the base position. Then, recompute incident planes as perturbation of the input planes,
    //   but with the constraint of passing through this anchor (with perturbPlaneCoefficientsFixedPoints()).
    // - If a facet has strictly more than 3 anchors, triangulate it.
    // - Loop the previous steps until there are no more new anchors. This must finish
    //   because at the worst we have triangulated all faces.

    {
      std::unordered_map<FacetSPtr, Plane_3> original_planes;
      std::unordered_map<VertexSPtr, FT> sq_max_displacements;
      std::unordered_set<VertexSPtr> anchors;

      for (const FacetSPtr& facet : polyhedron->facets()) {
        original_planes[facet] = facet->get_plane();
      }

      FT global_sq_max = CGAL::square(1e-5); // @tmp hardcoded
      for (const VertexSPtr& v : polyhedron->vertices()) {
        FT local_sq_max = global_sq_max;
        // @fixme this local bound ought to be LFS-based rather than incident edge length-based
        for (const EdgeWPtr& we : v->edges()) {
          if (const EdgeSPtr e = we.lock()) {
            VertexSPtr ov = e->other(v);
            const FT sq_dist = CGAL::squared_distance(v->point(), ov->point());
            const FT limit_sq = FT(0.01) * sq_dist;
            if (limit_sq < local_sq_max) {
              local_sq_max = limit_sq;
            }
          }
        }
        sq_max_displacements[v] = local_sq_max;
      }

      std::unordered_set<FacetSPtr> facets_to_recompute;
      for (const FacetSPtr& f : polyhedron->facets()) {
        facets_to_recompute.insert(f);
      }

      // Apply smaller and smaller perturbations, until all anchors are stable.
      // Since 0-perturbation has planes passing through their anchors, this process must converge.
      auto recompute_facet_planes = [&](const auto& facets) {
        // Start by computing the base plane, which passes through the anchors.
        // This matters because otherwise a perturbation going to 0 might not mean
        // that the plane passes through the anchors.

        std::unordered_map<FacetSPtr, Plane_3> pre_nudge_planes;
        std::unordered_map<FacetSPtr, FT> nudge_ranges;

        for (const FacetSPtr& f : facets) {
          CGAL_precondition(Transformation::has_normalized_plane(f));
          nudge_ranges[f] = nudge_range;

          // @fixme? do we need some kind of "link stability" to ensure the link does not tangle?...
          // --> maybe we don't because we have stability of all points, hence the edges are stable too
          CGAL_SS3_TRANSF_TRACE_V(16, "Recompute F" << f->id());
          CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                              << f->get_plane().c() << " " << f->get_plane().d() << "]");

          boost::container::small_vector<const Point_3*, 3> anchor_points;
          for (const VertexSPtr& v : f->vertices()) {
            if (anchors.count(v)) {
              anchor_points.push_back(&(v->point()));
              CGAL_SS3_TRANSF_TRACE_V(16, "    V" << v->id() << " is an anchor point at " << v->point());
            }
          }

          CGAL_assertion(anchor_points.size() <= 3);

          if (anchor_points.size() == 0) {
            perturbPlaneCoefficientsNudge(f, nudge_range);
          } else if (anchor_points.size() == 1) {
            const Point_3& p0 = *anchor_points[0];
            const Vector_3& n = f->get_plane().orthogonal_vector();
            const FT a = p0.x() * n.x();
            const FT b = p0.y() * n.y();
            const FT c = p0.z() * n.z();
            const FT d = -(a + b + c);
            f->set_plane(Plane_3(n.x(), n.y(), n.z(), d));
            CGAL_postcondition(f->get_plane().has_on(p0));
          } else if (anchor_points.size() == 2) {
            const Point_3& p0 = *anchor_points[0];
            const Point_3& p1 = *anchor_points[1];
            const FT ux = p1.x() - p0.x();
            const FT uy = p1.y() - p0.y();
            const FT uz = p1.z() - p0.z();
            const FT uu = ux * ux + uy * uy + uz * uz;

            const Vector_3& n = f->get_plane().orthogonal_vector();
            const FT dot = n.x() * ux + n.y() * uy + n.z() * uz;
            const FT ab = n.x() - (dot * ux / uu);
            const FT bb = n.y() - (dot * uy / uu);
            const FT cb = n.z() - (dot * uz / uu);

            const FT vx = uy * cb - uz * bb;
            const FT vy = uz * ab - ux * cb;
            const FT vz = ux * bb - uy * ab;

            static std::random_device rd;
            unsigned int s = 0; // rd()
            // CGAL_SS3_TRANSF_TRACE("seed = " << s);
            static std::mt19937 gen(s);
            static std::uniform_real_distribution<> rdist(-nudge_range, nudge_range);

            double eps = rdist(gen);

            FT a1 = ab + eps * vx;
            FT b1 = bb + eps * vy;
            FT c1 = cb + eps * vz;
            FT d1 = -(a1 * p0.x() + b1 * p0.y() + c1 * p0.z());
            f->set_plane(Plane_3(a1, b1, c1, d1));
            CGAL_postcondition(f->get_plane().has_on(p0));
            CGAL_postcondition(f->get_plane().has_on(p1));
          } else if (anchor_points.size() == 3) {
            const Point_3& p0 = *(anchor_points[0]);
            const Point_3& p1 = *(anchor_points[1]);
            const Point_3& p2 = *(anchor_points[2]);

            // @fixme should a large deviation of normal be understood as an unstable facet
            // that ought to be triangulated? If selected anchors were aligned, a small
            // nudge is a large normal change.
            // Maybe it's not needed because if the normal changed a lot, then it's likely
            // other vertices won't be stable and the facet will be triangulated.
            Plane_3 new_pl(p0, p1, p2);

            // We don't know about the order of the fixed points, so trust the input normal for orientation
            //
            // @fixme robustness when there are very thin triangles?... the nudge is small so it should be fine...
            if (new_pl.orthogonal_vector() * original_planes.at(f).orthogonal_vector() < 0) {
              new_pl = new_pl.opposite();
            }

            f->set_plane(new_pl);
            Transformation::normalize_plane_coefficients(f); // doesn't pass through the points anymore...

            // this should fail due to inexact normalization?...
            CGAL_postcondition(f->get_plane().has_on(p0));
            CGAL_postcondition(f->get_plane().has_on(p1));
            CGAL_postcondition(f->get_plane().has_on(p2));
          }

          CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");

          pre_nudge_planes[f] = f->get_plane();
        }

        for (;;) {
          // @todo don't recompute if the nudge range hasn't changed
          for (const FacetSPtr& f : facets) {
            // Now, we compute a nudge, using the fixed nudge direction
            const Plane_3& base_plane = pre_nudge_planes[f];
            const FT& nudge_coeff = nudge_ranges[f];
            Plane_3 nudged_plane { base_plane.a(), base_plane.b(), base_plane.c(), base_plane.d() + nudge_coeff };
            f->set_plane(nudged_plane);
            CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                              << f->get_plane().c() << " " << f->get_plane().d() << "]");

            // std::cout << "distance to anchors:" << std::endl;
            // for (const VertexSPtr& v : f->vertices()) {
            //   if (anchors.count(v)) {
            //     FT dist = CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), f->get_plane()));
            //     std::cout << "  V" << v->id() << " dist: " << dist << std::endl;
            //   }
            // }
          }

          // Now, all facets have updated (nudged) planes, so check the stability of the anchors
          bool stable_anchors = true;
          std::unordered_map<FacetSPtr, FT> new_nudge_ranges = nudge_ranges;
          for (const VertexSPtr& anchor : polyhedron->vertices()) { // @local
            if (!anchors.count(anchor)) {
              continue;
            }
            if (!is_stable_vertex(anchor, polyhedron, sq_max_displacements[anchor])) {
              CGAL_SS3_TRANSF_TRACE_V(16, "Anchor V" << anchor->id() << " is not stable, reducing incident plane nudges");
              stable_anchors = false;
              for (const FacetWPtr& wf : anchor->facets()) {
                if (FacetSPtr f = wf.lock()) {
                  new_nudge_ranges[f] = FT(0.01) * nudge_ranges.at(f);
                  CGAL_SS3_TRANSF_TRACE_V(16, "  Nudge range of F" << f->id() << " to " << new_nudge_ranges[f]);
                }
              }
            }
          }

          if (stable_anchors) {
            for (const FacetSPtr& f : facets) {
              Transformation::normalize_plane_coefficients(f);
            }
            break;
          }
          nudge_ranges.swap(new_nudge_ranges);
        }
      };

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream anchor_out("results/anchors.xyz");
      anchor_out.precision(17);
#endif

      auto make_random_anchor = [&](const VertexSPtr& v) {
        static std::random_device rd;
        unsigned int s = 0; // rd()
        // CGAL_SS3_TRANSF_TRACE("seed = " << s);
        static std::mt19937 gen(s);
        static std::uniform_real_distribution<> rdist(-nudge_range, nudge_range);

        double eps_x = rdist(gen);
        double eps_y = rdist(gen);
        double eps_z = rdist(gen);

        Point_3 new_p(v->point().x() + eps_x, v->point().y() + eps_y, v->point().z() + eps_z);
        CGAL_SS3_TRANSF_TRACE_V(8, "Make V" << v->id() << " an anchor at " << new_p);
        CGAL_SS3_TRANSF_TRACE_V(8, "  Distance from original: " << CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), new_p)));
        v->set_point(new_p);
      };

      auto has_dirty_incident_face = [&](const VertexSPtr& v) -> bool {
        for (const FacetWPtr& wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (facets_to_recompute.count(f)) {
              return true;
            }
          }
        }
        return false;
      };

      for (;;) {
        // Apply perturbation with given anchors, and note faces which have unstable vertices

#if 1
        // The smart algorithm does not yet work, but this complete range is a massive amount of useless work, recomputing facet planes over and over.
        recompute_facet_planes(polyhedron->facets());
        facets_to_recompute.clear();
#else
        recompute_facet_planes(facets_to_recompute);
        facets_to_recompute.clear();
#endif

        // facet perturbations have been recomputed with their temporary anchors
        // check stability of other vertices
        // @todo it's pointless to check vertices for which the incident facets have no changed.
        std::vector<VertexSPtr> new_anchors;
        for (const VertexSPtr& v : polyhedron->vertices()) {
          // if the facet is already marked for recomputation due to a new anchor, it's pointless
          // to add more anchors already. Instead, give it a chance to be more stable with the new anchor.
          if (has_dirty_incident_face(v)) {
            continue;
          }

          if (!is_stable_vertex(v, polyhedron, sq_max_displacements[v]) && !anchors.count(v)) {
            make_random_anchor(v);
            anchors.insert(v);
            new_anchors.push_back(v);
#ifdef CGAL_SS3_DUMP_FILES
            anchor_out << v->point() << "\n";
#endif
            for (FacetWPtr wf : v->facets()) {
              if (FacetSPtr f = wf.lock()) {
                facets_to_recompute.insert(f);
              }
            }
          }
        }

        if (new_anchors.empty()) {
          break;
        }

        // Some facets may now be overconstrained and must be triangulated
        std::unordered_set<FacetSPtr> facets_to_consider;

        for (const VertexSPtr& v : new_anchors) {
          for (const FacetWPtr& wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
              facets_to_consider.insert(f);
            }
          }
        }

        while (!facets_to_consider.empty()) {
          FacetSPtr f = *facets_to_consider.begin();
          facets_to_consider.erase(facets_to_consider.begin());

          std::size_t anchor_count = 0;
          for (const VertexSPtr& fv : f->vertices()) {
            if (anchors.count(fv)) {
              ++anchor_count;
            }
          }

          if (anchor_count > 3) {
            CGAL_SS3_TRANSF_TRACE_V(8, "  Must triangulate F" << f->id());
            Plane_3 original_plane = original_planes.at(f);
            auto [local_vertices, new_facets] = Transformation::triangulate_facet(f, polyhedron);
            CGAL_postcondition(polyhedron->is_consistent());

            original_planes.erase(f);
            facets_to_recompute.erase(f);
            for (const FacetSPtr& nf : new_facets) {
              original_planes[nf] = original_plane;
              facets_to_consider.insert(nf);
              facets_to_recompute.insert(nf);
            }
          }
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      anchor_out.close();
#endif

#ifdef CGAL_SS3_DUMP_FILES
      IO::write_OBJ("results/V4_preprocessed.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

      CGAL_postcondition_code(for (const VertexSPtr v : polyhedron->vertices()))
      CGAL_postcondition(is_stable_vertex(v, polyhedron, sq_max_displacements[v]));
    }

    CGAL_postcondition(do_all_plane_pairs_intersect(polyhedron));
    CGAL_postcondition(do_all_plane_triplets_intersect(polyhedron));

    // -- PART 2 --
    // The perturbation planes are set up, recompute all vertex positions, and split high-degree
    // vertices when needed.

    CGAL_SS3_TRANSF_TRACE_V(8, "Part 2: split high-degree vertices");

    std::list<VertexSPtr> vertices_tosplit;
    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (v->degree() > 3) {
        vertices_tosplit.push_back(v);
      }
    }

    CGAL_SS3_CORE_TRACE_V(4, vertices_tosplit.size() << " vertices to split");

    for (const VertexSPtr& vertex : vertices_tosplit) {
      CGAL_SS3_CORE_TRACE_V(8, "Splitting " << vertex->to_string());

      vertex->sort();

      // soup to be used in the arrangement
      std::vector<Point_3> points;
      std::vector<std::vector<std::size_t> > triangles;
      std::vector<FacetSPtr> triangle_to_facet;

      // Create a sufficiently large bounding box containing all plane intersections
      // @todo obviously !slightly! suboptimal
      CGAL::Bbox_3 bb = vertex->point().bbox();

      FT max_sq_displacement = 0;

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream inter_out("results/3-inter.xyz");
      inter_out.precision(17);
#endif

      // @todo no need to recompute this since we have stable vertices
      for (auto it_wf1 = vertex->facets().begin() ; it_wf1 != vertex->facets().end(); ++it_wf1) {
        if (FacetSPtr f1 = it_wf1->lock()) {
          for (auto it_wf2 = std::next(it_wf1) ; it_wf2 != vertex->facets().end(); ++it_wf2) {
            if (FacetSPtr f2 = it_wf2->lock()) {
              for (auto it_wf3 = std::next(it_wf2) ; it_wf3 != vertex->facets().end(); ++it_wf3) {
                if (FacetSPtr f3 = it_wf3->lock()) {
                  std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                  if (p_new.has_value()) {
                    FT sqd = CGAL::squared_distance(vertex->point(), *p_new);
                    max_sq_displacement = std::max(max_sq_displacement, sqd);

                    bb += p_new->bbox();
#ifdef CGAL_SS3_DUMP_FILES
                    inter_out << *(p_new) << "\n";
#endif
                  }
                }
              }
            }
          }
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(64, "max_sq_displacement = " << max_sq_displacement);

      bb.scale(1.5);

      Iso_cuboid_3 bbox { bb };

      CGAL_SS3_TRANSF_TRACE_V(64, "splitting bounding box: " << bbox);
      CGAL_SS3_TRANSF_TRACE_V(64, "x span " << bb.x_span() << ", y span " << bb.y_span() << ", z span " << bb.z_span());

      // Get the triangles from the base planes, with tagged faces
      get_clipped_plane_faces(vertex, bbox, points, triangles, triangle_to_facet);
      CGAL_SS3_TRANSF_TRACE_V(64, points.size() << " points, " << triangles.size() << " triangles [base]");

#ifdef CGAL_SS3_DUMP_FILES
      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            id_to_triangles[fsptr->id()].push_back(triangles[i]);
          }
        }

        for (const auto& kv : id_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_base_" << kv.first << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
        }

        CGAL_assertion(id_to_triangles.size() == vertex->degree());
      }
#endif

      // Add the bounding box faces
      std::size_t base_idx = points.size();
      points.push_back(Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin())); // v0
      points.push_back(Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin())); // v1
      points.push_back(Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin())); // v2
      points.push_back(Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin())); // v3
      points.push_back(Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax())); // v4
      points.push_back(Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax())); // v5
      points.push_back(Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax())); // v6
      points.push_back(Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax())); // v7

      // bottom face
      triangles.push_back({base_idx+0, base_idx+2, base_idx+1});
      triangles.push_back({base_idx+0, base_idx+3, base_idx+2});
      // top face
      triangles.push_back({base_idx+4, base_idx+5, base_idx+6});
      triangles.push_back({base_idx+4, base_idx+6, base_idx+7});
      // front face
      triangles.push_back({base_idx+0, base_idx+1, base_idx+5});
      triangles.push_back({base_idx+0, base_idx+5, base_idx+4});
      // back face
      triangles.push_back({base_idx+2, base_idx+3, base_idx+7});
      triangles.push_back({base_idx+2, base_idx+7, base_idx+6});
      // left face
      triangles.push_back({base_idx+0, base_idx+4, base_idx+7});
      triangles.push_back({base_idx+0, base_idx+7, base_idx+3});
      // right face
      triangles.push_back({base_idx+1, base_idx+2, base_idx+6});
      triangles.push_back({base_idx+1, base_idx+6, base_idx+5});

      triangle_to_facet.resize(triangles.size(), nullptr);

      CGAL_SS3_TRANSF_TRACE_V(64, points.size() << " points, " << triangles.size() << " triangles [w/ bbox]");

      for (std::size_t i=0; i<triangle_to_facet.size(); ++i) {
        CGAL_SS3_TRANSF_TRACE_V(64, "triangle " << i << " ptr: "
                                << (triangle_to_facet[i] ? triangle_to_facet[i]->id() : -1));
      }

      PMP::merge_duplicate_points_in_polygon_soup(points, triangles);

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_OFF("results/arr_soup.off", points, triangles, CGAL::parameters::stream_precision(17));
#endif

      CGAL_assertion(triangles.size() == triangle_to_facet.size());

      CGAL_SS3_TRANSF_TRACE_V(1, "autorefining...");
      std::vector<FacetSPtr> updated_triangle_to_facet;
      Range_updating_autoref_visitor<FacetSPtr> autoref_visitor(triangle_to_facet, updated_triangle_to_facet);

      PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::visitor(autoref_visitor));
      triangle_to_facet = std::move(updated_triangle_to_facet);
      CGAL_assertion(triangles.size() == triangle_to_facet.size());

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_OFF("results/arr_autorefined.off", points, triangles, CGAL::parameters::stream_precision(17));
#endif

      // dump the arrangement, with normalized points fitting in a unit cube, for easier visualization
      std::vector<Point_3> normalized_points;
      for (const Point_3& p : points) {
        normalized_points.push_back(Point_3((p.x() - bbox.xmin()) / bb.x_span(),
                                            (p.y() - bbox.ymin()) / bb.y_span(),
                                            (p.z() - bbox.zmin()) / bb.z_span()));
      }

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_OFF("results/arr_autorefined_normalized.off", normalized_points, triangles, CGAL::parameters::stream_precision(17));
#endif

      for (std::size_t i=0; i<triangles.size(); ++i) {
        CGAL_SS3_TRANSF_TRACE_V(64, "[4] triangle " << i << " ptr: " << (triangle_to_facet[i] ? triangle_to_facet[i]->id() : -1));
      }

#ifdef CGAL_SS3_DUMP_FILES
      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            id_to_triangles[fsptr->id()].push_back(triangles[i]);
          }
        }

        for (const auto& kv : id_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_autorefined_" << kv.first << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
          std::ostringstream oss_n;
          oss_n << "results/arr_autorefined_" << kv.first << "_normalized.off";
          CGAL::IO::write_OFF(oss_n.str(), normalized_points, kv.second, CGAL::parameters::stream_precision(17));
        }

        CGAL_assertion(id_to_triangles.size() == vertex->degree());
      }
#endif

      CGAL_SS3_TRANSF_TRACE_V(64, "repairing...");
      Range_updating_repair_PS_visitor<FacetSPtr> repair_ps_visitor(triangle_to_facet);
      PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                    CGAL::parameters::visitor(repair_ps_visitor)
                                                                    .erase_all_duplicates(false) /*keep one*/
                                                                    .require_same_orientation(false)
                                                                    .verbose(true));

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_OFF("results/arr_repaired.off", points, triangles, CGAL::parameters::stream_precision(17));
#endif

#ifdef CGAL_SS3_DUMP_FILES
      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            id_to_triangles[fsptr->id()].push_back(triangles[i]);
          }
        }

        for (const auto& kv : id_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_repaired_" << kv.first << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
        }
      }
#endif

      CGAL_assertion(triangles.size() == triangle_to_facet.size());

#ifdef CGAL_SS3_DUMP_FILES
      // dump only the triangles with a non nullptr
      {
        std::vector<std::vector<std::size_t> > input_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          if (triangle_to_facet[i]) {
            input_triangles.push_back(triangles[i]);
          }
        }

        CGAL::IO::write_OFF("results/arr_recovered_input.off", points, input_triangles, CGAL::parameters::stream_precision(17));
      }
#endif

      using PID = std::size_t;
      using TID = std::size_t;
      using VID = std::size_t;

      auto fill_edge_map = [](const std::vector<Point_3>& points,
                              const std::vector<std::vector<PID> >& triangles,
                              std::vector<std::unordered_map<PID, std::vector<TID> > >& edge_map)
      {
        CGAL_precondition(edge_map.size() == points.size());

        // collect duplicated edges
        for (TID ti=0; ti<triangles.size(); ++ti) {
          for (std::size_t j=0; j<3; ++j) {
            std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[ti][j],
                                                                triangles[ti][(j+1)%3]);
            edge_map[e_pids.first][e_pids.second].push_back(ti);
          }
        }

        namespace pred = CGAL::Polygon_mesh_processing::Corefinement;

        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
          for (auto& pid1_and_edges : edge_map[pid0]) {
            std::vector<TID>& inc_triangles = pid1_and_edges.second;

            if (inc_triangles.size() == 2) { // edge is only incident to a single SS3 face
              continue;
            }

            const PID pid1 = pid1_and_edges.first;
            CGAL_assertion(pid0 != pid1);

            CGAL_SS3_TRANSF_TRACE_V(64, "processing NM edge: " << points[pid0] << " -- " << points[pid1]);

            auto get_third_point_id = [&triangles, pid0, pid1](TID tid) -> PID
            {
              std::size_t third;

              // need to be careful that the orientation of the edge might not match the orientation of the triangle
              if (triangles[tid][0] == pid0 || triangles[tid][0] == pid1) {
                if (triangles[tid][1] == pid0 || triangles[tid][1] == pid1) {
                  third = triangles[tid][2];
                } else {
                  third = triangles[tid][1];
                }
              } else {
                third = triangles[tid][0];
              }

              CGAL_postcondition(third != pid0 && third != pid1);
              return third;
            };

            const Point_3& ref_pt = points.at(get_third_point_id(inc_triangles[0]));
            auto less = [&ref_pt, &points, pid0, pid1, get_third_point_id](TID tid1, TID tid2)
            {
              return pred::sorted_around_edge<GeomTraits>(points.at(pid0), points.at(pid1),
                                                          ref_pt,
                                                          points.at(get_third_point_id(tid1)),
                                                          points.at(get_third_point_id(tid2)));
            };

            std::sort(inc_triangles.begin()+1, inc_triangles.end(), less);

            // std::cout << "Around edge [" << pid0 << " " << pid1 << "], faces are sorted: ";
            // for(TID tid : inc_triangles)
            //   std::cout << " " << tid;
            // std::cout << std::endl;
          }
        }
      }; // lambda 'fill_edge_map'

      std::vector<std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
      fill_edge_map(points, triangles, edge_map);

      auto build_volume_CC = [](const TID seed_tid,
                                const VID CC_ID,
                                const bool start_from_inverted_face,
                                const std::vector<Point_3>& points,
                                const std::vector<std::vector<PID> >& triangles,
                                const auto& edge_map,
                                auto& volume_CCs,
                                auto& face_volume_IDs)
      {
        CGAL_SS3_TRANSF_TRACE_V(64, "Building volume #" << CC_ID << " from seed face " << seed_tid);

        volume_CCs.emplace_back();

        std::stack<std::pair<TID, bool> > to_visit;
        to_visit.emplace(seed_tid, start_from_inverted_face);

        while (!to_visit.empty())
        {
          TID current_tid;
          bool invert_face;
          std::tie(current_tid, invert_face) = to_visit.top();
          to_visit.pop();

          CGAL_SS3_TRANSF_TRACE_V(64,
            "At face " << current_tid << " [" << triangles[current_tid][0]
                                      << ", " << triangles[current_tid][1]
                                      << ", " << triangles[current_tid][2] << "]");
          CGAL_SS3_TRANSF_TRACE_V(64, "  invert: " << invert_face);
          CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS: " << face_volume_IDs[current_tid][0] << " " << face_volume_IDs[current_tid][1]);

          std::size_t pos = invert_face ? 0 : 1;
          if (face_volume_IDs[current_tid][pos] == CC_ID) {
            // already visited this facet during the flooding of this volume's boundary
            continue;
          }

          CGAL_assertion(face_volume_IDs[current_tid][pos] == VID(-1)); // triangle should only be encountered once
          CGAL_warning(face_volume_IDs[current_tid][(pos+1)%2] != CC_ID); // Moebius shenanigans should be an instance of a bug

          volume_CCs.back().push_back(current_tid);

          // mark face as visited
          face_volume_IDs[current_tid][pos] = CC_ID;

          // flood through the edges
          for (int j=0; j<3; ++j) {
            std::pair<PID, PID> e_pids = CGAL::make_sorted_pair(triangles[current_tid][j],
                                                                triangles[current_tid][(j+1)%3]);
            const std::vector<TID>& inc_triangles = edge_map.at(e_pids.first).at(e_pids.second);
            CGAL_assertion(!inc_triangles.empty());

            CGAL_SS3_TRANSF_TRACE_V(64, "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]");
            CGAL_SS3_TRANSF_TRACE_V(64, "    pos: " << points[e_pids.first] << " " << points[e_pids.second]);

            // The faces are ordered CCW while looking from pid0.
            // So the walking while looking from [j] depends on whether [j] is pid0 or not
            int iter_direction = (e_pids.first == triangles[current_tid][j]) ? 1 : -1;
            CGAL_SS3_TRANSF_TRACE_V(64, "    iter_direction = " << iter_direction);

            // and it also depends on whether we are walking above or below the face
            iter_direction *= invert_face ? 1 : -1;
            CGAL_SS3_TRANSF_TRACE_V(64, "    invert_face = " << invert_face);

            TID next_tid = current_tid;
            for (;;) {
              if (inc_triangles.size() == 1) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Warning: dangling triangle...");
                CGAL_SS3_TRANSF_TRACE_V(1, "    over the edge, the triangle is ITSELF " << current_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "]");
                to_visit.emplace(current_tid, !invert_face);
                break;
              } else if (inc_triangles.size() == 2) {
                // we should only be there once, meaning if we do not ignore orientations,
                // then the faces MUST be compatible
                CGAL_assertion(next_tid == current_tid);

                next_tid = (inc_triangles[0] == current_tid) ? inc_triangles[1] : inc_triangles[0];
                CGAL_SS3_TRANSF_TRACE_V(64, "    over the edge, the triangle is TRIVIALLY " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "]");
                CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1]);
                CGAL_assertion(next_tid != current_tid);
              } else {
                // tricky part, now
                auto tid_it = std::find(std::begin(inc_triangles), std::end(inc_triangles), next_tid /*updates on every iteration*/);
                CGAL_assertion(tid_it != inc_triangles.end());

                if (iter_direction == 1) { // CCW
                  CGAL_SS3_TRANSF_TRACE_V(64, "    CCW walk");
                  auto next_it = std::next(tid_it);
                  next_tid = (next_it == inc_triangles.end()) ? inc_triangles[0] : *next_it;
                } else { // CW
                  CGAL_SS3_TRANSF_TRACE_V(64, "    CW walk");
                  next_tid = (tid_it == inc_triangles.begin()) ? inc_triangles.back() : *(std::prev(tid_it));
                }

                CGAL_SS3_TRANSF_TRACE_V(64, "    over the edge, the triangle is " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "]");
                CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1]);
                CGAL_assertion(next_tid != current_tid);
              }

              // If the edge has the same direction in both faces (aka, the orientation changes),
              // then we have to flip the direction of turning around the edge)
              const auto j_it = std::find(std::begin(triangles[next_tid]),
                                          std::end(triangles[next_tid]),
                                          triangles[current_tid][j]);
              CGAL_assertion(j_it != std::end(triangles[next_tid]));
              const std::size_t pos = std::distance(std::begin(triangles[next_tid]), j_it);
              CGAL_assertion(triangles[next_tid][pos] == triangles[current_tid][j]);

              const bool flip_side = (triangles[next_tid][(pos+1)%3] == triangles[current_tid][(j+1)%3]);
              CGAL_SS3_TRANSF_TRACE_V(64, "    flipping? " << flip_side << " (N: " << triangles[next_tid][(pos+1)%3] << " C: " << triangles[current_tid][(j+1)%3] << ")");

              CGAL_SS3_TRANSF_TRACE_V(64, "Final TID = " << next_tid);
              if (flip_side) {
                to_visit.emplace(next_tid, !invert_face);
              } else {
                to_visit.emplace(next_tid, invert_face);
              }

              break;
            }
          }
        }
      }; // lambda 'build_volume_CC'

      CGAL_SS3_TRANSF_TRACE_V(1, "building volumes...");

      // identify volumes in the arrangement, and tag faces of the volumes
      // that are incident to the base face(s)
      std::vector<std::vector<TID> > volume_CCs; // range of ranges (volumes) of triangle IDs
      std::vector<std::array<VID, 2> > face_volume_IDs(triangles.size(),
                                                       // [0] is down, [1] is up
                                                       std::array<VID, 2>{VID(-1), VID(-1)});

      VID vid = 0;
      for(std::size_t i=0; i<triangles.size(); ++i) {
        // do not start from bbox faces: we will visit them anyway
        if (!triangle_to_facet[i])
          continue;
        if (face_volume_IDs[i][0] == VID(-1))
          build_volume_CC(i, vid++, true, points, triangles, edge_map, volume_CCs, face_volume_IDs);
        if (face_volume_IDs[i][1] == VID(-1))
          build_volume_CC(i, vid++, false, points, triangles, edge_map, volume_CCs, face_volume_IDs);
      }

      CGAL_SS3_TRANSF_TRACE_V(1, volume_CCs.size() << " volume CCs");

      for (std::size_t i=0; i<volume_CCs.size(); ++i) {
        // build a mesh from the soup
        std::vector<Point_3> cc_points = points;
        std::vector<std::vector<PID> > cc_triangles;
        for (TID tid : volume_CCs[i]) {
          cc_triangles.push_back(triangles[tid]);
        }

#ifdef CGAL_SS3_DUMP_FILES
        std::ostringstream oss;
        oss << "results/volume_cc_" << i << ".off";
        CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
#endif
      }

      enum class CC_in_out_flag
      {
        UNINITIALIZED = 0,
        TBD,
        INSIDE,
        OUTSIDE
      };

      std::vector<CC_in_out_flag> in_out_flags(volume_CCs.size(), CC_in_out_flag::UNINITIALIZED);
      std::size_t classified_CCs = 0;

      // Classify some trivial CCs:
      // - if every non-bbox face points outwards, the CC is necessarily in
      // - if every non-bbox face points inwards, the CC is necessarily out
      for(std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        // bottom
        VID bot_vid = face_volume_IDs[i][0];
        CGAL_assertion(bot_vid != VID(-1));
        if (in_out_flags[bot_vid] == CC_in_out_flag::UNINITIALIZED)
          in_out_flags[bot_vid] = CC_in_out_flag::INSIDE; // [0], bottom, face points outwards
        else if (in_out_flags[bot_vid] == CC_in_out_flag::OUTSIDE)
          in_out_flags[bot_vid] = CC_in_out_flag::TBD;

        // top
        VID top_vid = face_volume_IDs[i][1];
        CGAL_assertion(top_vid != VID(-1));
        if (in_out_flags[top_vid] == CC_in_out_flag::UNINITIALIZED)
          in_out_flags[top_vid] = CC_in_out_flag::OUTSIDE; // [1], top, face points inwards
        else if (in_out_flags[top_vid] == CC_in_out_flag::INSIDE)
          in_out_flags[top_vid] = CC_in_out_flag::TBD;
      }

#ifdef CGAL_SS3_DUMP_FILES
      // base info dump
      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(in_out_flags[i] != CC_in_out_flag::UNINITIALIZED);
          if (in_out_flags[i] == CC_in_out_flag::INSIDE) {
            CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known inside (base)");
          } else if (in_out_flags[i] == CC_in_out_flag::OUTSIDE) {
            CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known outside (base)");
          } else {
            ++undetermined_n;
          }
        }
        CGAL_SS3_TRANSF_TRACE_V(64, undetermined_n << " undetermined cells");

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps) {
          std::vector<Point_3> cc_points = points;
          std::vector<std::vector<PID> > cc_triangles;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            CGAL_assertion(e.first != CC_in_out_flag::TBD);
            CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

            if (in_out_flags[i] != e.first)
              continue;

            for (TID tid : volume_CCs[i])
              cc_triangles.push_back(triangles[tid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_base.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_base_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_triangles, CGAL::parameters::stream_precision(17));
          }
        }
      }
#endif

      // Set up the non obvious known volumes

      boost::dynamic_bitset<> is_boundary_point(points.size(), 0);
      for (std::size_t pid=0; pid<points.size(); ++pid) {
        const Point_3& p = points[pid];
        if (p.x() == bbox.xmin() || p.x() == bbox.xmax() ||
            p.y() == bbox.ymin() || p.y() == bbox.ymax() ||
            p.z() == bbox.zmin() || p.z() == bbox.zmax()) {
          is_boundary_point.set(pid);
        }
      }

      // Generalization of above: for the complete trace of the star on the bounding box,
      // we know the cell above is OUTSIDE and the cell beneath is INSIDE

      // We need to have the border incident edges for any type of boundary point
      std::unordered_map<PID, std::unordered_map<PID, std::vector<TID> > > symmetrical_border_edge_map;
      for (TID ti=0; ti<triangles.size(); ++ti) {
        for (std::size_t j=0; j<3; ++j) {
          const PID pid0 = triangles[ti][j];
          const PID pid1 = triangles[ti][(j+1)%3];
          if (!is_boundary_point[pid0] || !is_boundary_point[pid1]) {
            continue;
          }
          // avoid edges across the box
          const Point_3 m = CGAL::midpoint(points[pid0], points[pid1]);
          if (m.x() == bbox.xmin() || m.x() == bbox.xmax() ||
              m.y() == bbox.ymin() || m.y() == bbox.ymax() ||
              m.z() == bbox.zmin() || m.z() == bbox.zmax()) {
            symmetrical_border_edge_map[pid0][pid1].push_back(ti);
            symmetrical_border_edge_map[pid1][pid0].push_back(ti);
          }
        }
      }

      // Helper to mark cell directly above and below a facet incident to an edge
      auto mark_cells = [&](const PID ipid0, const PID ipid1, const FacetSPtr& current_facet)
      {
        const auto [pid0, pid1] = CGAL::make_sorted_pair(ipid0, ipid1);
        const auto edge_it = symmetrical_border_edge_map.find(pid0);
        CGAL_assertion(edge_it != symmetrical_border_edge_map.end());
        const auto pid1_it = edge_it->second.find(pid1);
        CGAL_assertion(pid1_it != edge_it->second.end());

        const auto& tids = pid1_it->second;
        for (TID tid : tids) {
          if (triangle_to_facet[tid] != current_facet) {
            continue;
          }

          in_out_flags[face_volume_IDs[tid][1]] = CC_in_out_flag::OUTSIDE;
          in_out_flags[face_volume_IDs[tid][0]] = CC_in_out_flag::INSIDE;

          return;
        }
      };

      // Helper to find next boundary point on current_facet in the direction of next_facet
      auto find_next_boundary_point = [&](const PID prev_pid,
                                          const PID current_pid,
                                          const FacetSPtr& prev_facet,
                                          const FacetSPtr& current_facet,
                                          const FacetSPtr& next_facet) -> std::pair<PID, bool>
      {
        std::set<PID> incident_boundary_points;

        auto edge_it = symmetrical_border_edge_map.find(current_pid);
        CGAL_assertion(edge_it != symmetrical_border_edge_map.end());

        for (const auto& [other_pid, tids] : edge_it->second) {
          CGAL_assertion(is_boundary_point[current_pid] && is_boundary_point[other_pid]);

          bool on_current_facet = false;
          for (TID tid : tids) {
            if (triangle_to_facet[tid] == current_facet) {
              on_current_facet = true;
              break;
            }
          }

          if (on_current_facet) {
            incident_boundary_points.insert(other_pid);
          }
        }

        CGAL_SS3_TRANSF_TRACE_V(64, incident_boundary_points.size() << " incident boundary points");
        CGAL_assertion(incident_boundary_points.size() == 2);

        PID first_pid = *(incident_boundary_points.begin());
        PID second_pid = *(std::next(incident_boundary_points.begin()));

        if (prev_pid != PID(-1)) {
          return {(prev_pid == first_pid ? second_pid : first_pid), true};
        }

        // If we are here, candidate_pid is the point that is at the intersection
        // of prev_facet and facet and we need to know in which direction to walk
        // along facet.

        const Plane_3& pl = current_facet->get_plane();

        // determine convexity of prev/current
        std::list<EdgeSPtr> common_edges = prev_facet->find_edges(current_facet);
        CGAL_assertion(!common_edges.empty());

        // can't trust the edge is_reflex() function because at that point because the facets
        // have been perturbed and the convexity might have changed...
        auto is_reflex = [&](const EdgeSPtr& e) -> bool {
          bool result = false;
          const FacetSPtr facet_l = e->get_facet_L();
          const FacetSPtr facet_r = e->get_facet_R();
          CGAL_SS3_DEBUG_SPTR(facet_l);
          CGAL_SS3_DEBUG_SPTR(facet_r);
          const Plane_3& plane_l = facet_l->get_plane();
          const Plane_3& plane_r = facet_r->get_plane();
          std::optional<Line_3> oline = Kernel_wrapper::intersection(plane_l, plane_r);
          CGAL_assertion(bool(oline));
          Vector_3 dir = oline->to_vector();
          CGAL_assertion(dir != CGAL::NULL_VECTOR);
          // possibly reorient 'dir' to align with the direction of the edge
          // note that the edge is stable since its vertices are stable
          if (dir * Vector_3(e->source()->point(), e->target()->point()) < 0) {
            dir = -dir;
          }
          const Point_3 p_src = oline->point();
          const Vector_3 normal_l = plane_l.orthogonal_vector();
          CGAL_assertion(normal_l != CGAL::NULL_VECTOR);
          Point_3 p = p_src + CGAL::cross_product(normal_l, dir);
          if (plane_r.oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
            result = true;
          }
          return result;
        };


        bool is_convex = !(is_reflex(common_edges.front()));
        CGAL_SS3_TRANSF_TRACE_V(64, "is_convex = " << is_convex);

        // now, if the edge is convex in the input polyhedron, we must turn "right", meaning,
        // the next point is on the negative side of the plane
        if (prev_facet->get_plane().oriented_side(points[first_pid]) == CGAL::NEGATIVE)
          return {(is_convex ? first_pid : second_pid), true};
        else
          return {(is_convex ? second_pid : first_pid), true};

        return {PID(-1), false};
      };

      // Walk the trace of the star on the bbox ("box link")
      EdgeSPtr start_edge = vertex->first_edge();
      EdgeSPtr edge = start_edge;

      std::map<EdgeSPtr, PID> edge_start_pid;
      auto compute_start_pid = [&](const EdgeSPtr& edge) -> PID {
        FacetSPtr prev_facet = edge->right(vertex);
        FacetSPtr facet = edge->left(vertex);

        VertexSPtr other_v = edge->other(vertex);
        Vector_3 orig_v { vertex->point(), other_v->point() };
        CGAL_assertion(orig_v != CGAL::NULL_VECTOR);

        CGAL_SS3_TRANSF_TRACE_V(64, "compute start of " << prev_facet->id() << " " << facet->id());

        for (PID pid0 = 0; pid0 < points.size(); ++pid0) {
          for (const auto& [pid1, tids] : edge_map[pid0]) {
            // The edge should be incident to only these two planes.
            if (tids.size() != 4)
              continue;

            // The edge should have one point on the boundary, one point inside.
            // It cannot cross the inner box because another plane will intersect the line.
            if (is_boundary_point[pid0] == is_boundary_point[pid1]) {
              continue;
            }

            std::set<FacetSPtr> incident_faces;
            for (TID tid : tids)
              incident_faces.insert(triangle_to_facet[tid]);
            CGAL_assertion(incident_faces.size() == 2);

            if (!incident_faces.count(prev_facet) || !incident_faces.count(facet)) {
              continue;
            }

            bool flip = is_boundary_point[pid0];
            PID test_spid0 = flip ? pid1 : pid0;
            PID test_spid1 = flip ? pid0 : pid1;
            CGAL_assertion(is_boundary_point[test_spid1]);

            CGAL_SS3_TRANSF_TRACE_V(64, "flip? " << flip);
            CGAL_SS3_TRANSF_TRACE_V(64, "candidate 0 " << normalized_points[test_spid0]);
            CGAL_SS3_TRANSF_TRACE_V(64, "candidate 1 " << normalized_points[test_spid1]);

            Vector_3 new_v { points[test_spid0], points[test_spid1] };
            CGAL_SS3_TRANSF_TRACE_V(64, "sp = " << CGAL::scalar_product(orig_v, new_v));

            if (CGAL::scalar_product(orig_v, new_v) >= 0) {
              CGAL_SS3_TRANSF_TRACE_V(64, "start is " << normalized_points[test_spid1]);
              return test_spid1;
            }
          }
        }

        CGAL_assertion(false);
        return PID(-1);
      };

      do {
        edge_start_pid[edge] = compute_start_pid(edge);
        edge = edge->next(vertex);
      } while(edge != start_edge);

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream link_out("results/link.txt");
      link_out.precision(17);
#endif

      edge = start_edge;
      do {
        FacetSPtr prev_facet = edge->right(vertex);
        FacetSPtr facet = edge->left(vertex);
        EdgeSPtr next_edge = edge->next(vertex);
        FacetSPtr next_facet = next_edge->left(vertex);

        CGAL_SS3_TRANSF_TRACE_V(64, "walking on " << facet->id() << " (prev: " << prev_facet->id() << ")");

        PID current_pid = edge_start_pid[edge];
        CGAL_assertion(current_pid != PID(-1));
        CGAL_SS3_TRANSF_TRACE_V(64, "Start at " << current_pid << " pos: " << normalized_points[current_pid]);

        // Now walk the star link on arrangement edges

        PID prev_pid (-1);
        for(;;) {
          CGAL_assertion(facet->get_plane().has_on(points[current_pid]));

          // find the next point on 'facet' in the direction of 'next_facet'
          auto [next_pid, valid] = find_next_boundary_point(prev_pid, current_pid, prev_facet, facet, next_facet);
          CGAL_assertion(valid);
          CGAL_assertion(is_boundary_point[current_pid] && is_boundary_point[next_pid]);

#ifdef CGAL_SS3_DUMP_FILES
          link_out << "2 " << normalized_points[current_pid] << " " << normalized_points[next_pid] << "\n";
#endif
          CGAL_SS3_TRANSF_TRACE_V(64, "now at " << next_pid << " pos: " << normalized_points[next_pid]);

          mark_cells(current_pid, next_pid, facet);

          prev_pid = current_pid;
          current_pid = next_pid;

          if (current_pid == edge_start_pid[next_edge]) {
            break;
          }
        }

        edge = next_edge;
      }
      while (edge != start_edge);

#ifdef CGAL_SS3_DUMP_FILES
      link_out.close();
#endif

      CGAL_SS3_TRANSF_TRACE_V(64, "Known cells marked from star link walk");

#ifdef CGAL_SS3_DUMP_FILES
      // intermediate info dump
      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(in_out_flags[i] != CC_in_out_flag::UNINITIALIZED);
          if (in_out_flags[i] == CC_in_out_flag::INSIDE) {
             CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known inside (boundary edges)");
          } else if (in_out_flags[i] == CC_in_out_flag::OUTSIDE) {
            CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known outside (boundary edges)");
          } else {
            ++undetermined_n;
          }
        }
        CGAL_SS3_TRANSF_TRACE_V(64, undetermined_n << " undetermined cells");

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps) {
          std::vector<Point_3> cc_points = points;
          std::vector<std::vector<PID> > cc_triangles;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            if (in_out_flags[i] != e.first)
              continue;

            for (TID tid : volume_CCs[i])
              cc_triangles.push_back(triangles[tid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_intermediate.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));

            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_intermediate_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_triangles, CGAL::parameters::stream_precision(17));
          }
        }

        // dump boundary facets
        {
          CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_triangles;

          for (std::size_t i=0; i<triangles.size(); ++i) {
            FacetSPtr fsptr = triangle_to_facet[i];
            if (!fsptr)
              continue;

            VID bot_vid = face_volume_IDs[i][0];
            VID top_vid = face_volume_IDs[i][1];
            CGAL_assertion(bot_vid != VID(-1) && top_vid != VID(-1));

            if (in_out_flags[top_vid] == CC_in_out_flag::OUTSIDE && in_out_flags[bot_vid] == CC_in_out_flag::INSIDE) {
              CGAL_SS3_TRANSF_TRACE_V(64, "Boundary facet " << i << " with top/bottom volumes " << top_vid << "/" << bot_vid);
              id_to_triangles[fsptr->id()].push_back(triangles[i]);
            }
          }

          for (const auto& kv : id_to_triangles) {
            std::ostringstream oss;
            oss << "results/intermediate_boundary_facet_" << kv.first << ".off";
            CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/intermediate_boundary_facet_" << kv.first << "_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, kv.second, CGAL::parameters::stream_precision(17));
          }

          CGAL_assertion(id_to_triangles.size() == vertex->degree());
        }
      }
#endif

      // Sanity checks:
      // - a boundary facet cannot have an OUTSIDE volume on its bottom and an INSIDE volume on its top
      for (std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        VID bot_vid = face_volume_IDs[i][0];
        VID top_vid = face_volume_IDs[i][1];
        CGAL_assertion(!(in_out_flags[top_vid] == CC_in_out_flag::INSIDE && in_out_flags[bot_vid] == CC_in_out_flag::OUTSIDE));
      }

      // Now, we want to deal with the tentative cells

      using operations_research::Domain;
      using operations_research::sat::BoolVar;
      using operations_research::sat::IntVar;
      using operations_research::sat::CpModelBuilder;
      using operations_research::sat::CpSolverResponse;
      using operations_research::sat::CpSolverStatus;
      using operations_research::sat::LinearExpr;
      using operations_research::sat::SatParameters;
      using operations_research::sat::SolutionBooleanValue;
      using operations_research::sat::SolveWithParameters;

      CpModelBuilder model;

      const std::size_t C = volume_CCs.size();

      std::vector<BoolVar> x;
      x.reserve(C);
      for (int c = 0; c < C; ++c)
        x.push_back(model.NewBoolVar());

      // =====================================================================
      //  CONSTRAINT — Fixed cells
      // =====================================================================

      // cells that are known, are known
      unsigned int x_unknowns = 0;
      for (int c = 0; c < C; ++c) {
        if (in_out_flags[c] == CC_in_out_flag::OUTSIDE) {
          model.FixVariable(x[c], false);
        } else if (in_out_flags[c] == CC_in_out_flag::INSIDE) {
          model.FixVariable(x[c], true);
        } else {
          ++x_unknowns;
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(8, "Cell unknowns: " << x_unknowns << " / " << C);

      // =====================================================================
      //  CONSTRAINT — Boundary facets point outside
      // =====================================================================

      // The boundary must be well oriented, hence a boundary facet must have the 'INSIDE' marker
      // on its negative side, and the 'OUTSIDE' marker on its positive side.
      // In model terms, this is:
      //   x[pos] => x[neg]. (if above is true, below is true)
      for(std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        model.AddImplication(x[face_volume_IDs[i][1]], x[face_volume_IDs[i][0]]);
      }

      // @speed could use hints using intersection of volumes?

      // Variables for the facets
      std::vector<BoolVar> b(triangles.size());
      unsigned int b_unknowns = 0;
      for (std::size_t tid=0; tid<triangles.size(); ++tid) {
        b[tid] = model.NewBoolVar();
        const VID bot_vid = face_volume_IDs[tid][0];
        const VID top_vid = face_volume_IDs[tid][1];

        if (bot_vid == VID(-1) || top_vid == VID(-1)) {
          model.FixVariable(b[tid], false);
        } else if (in_out_flags[bot_vid] == CC_in_out_flag::INSIDE && in_out_flags[top_vid] == CC_in_out_flag::OUTSIDE) {
          model.FixVariable(b[tid], true);
        } else { // @todo we can fix the facet that have incident fixed cells
          // A facet is on the boundary if and only if its two incident volumes are different
          model.AddNotEqual(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[tid]);
          model.AddEquality(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[tid].Not());
          ++b_unknowns;
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(8, "Face unknowns: " << b_unknowns << " / " << triangles.size());

      // =====================================================================
      //  CONSTRAINT — Global Euler characteristic (V_b − E_b + F_b = 1)
      // =====================================================================

      std::unordered_map<PID, std::unordered_set<TID> > vertex_incident_facets;
      for (TID ti=0; ti<triangles.size(); ++ti) {
        for (int j=0; j<3; ++j)
          vertex_incident_facets[triangles[ti][j]].insert(ti);
      }

      // F_b
      LinearExpr F_b = LinearExpr::Sum(b);

      // E_b : edge on surface iff >= 1 incident boundary facet
      std::vector<BoolVar> e_b;
      e_b.reserve(edge_map.size());

      for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
        for (auto& pid1_and_edges : edge_map[pid0]) {
          e_b.push_back(model.NewBoolVar());
          std::vector<TID>& inc_triangles = pid1_and_edges.second;
          CGAL_assertion(!inc_triangles.empty());
          std::vector<IntVar> inc;
          inc.reserve(inc_triangles.size());
          for (TID tid : inc_triangles)
            inc.emplace_back(b[tid]);
          model.AddMaxEquality(e_b.back(), inc);
        }
      }
      LinearExpr E_b = LinearExpr::Sum(e_b);

      // V_b : vertex on surface iff >= 1 incident boundary facet
      std::vector<BoolVar> v_b(points.size());
      for (int v=0; v<points.size(); ++v) {
        v_b[v] = model.NewBoolVar();
        std::unordered_set<TID>& inc_triangles = vertex_incident_facets[v];
        CGAL_assertion(!inc_triangles.empty());
        std::vector<IntVar> inc;
        inc.reserve(inc_triangles.size());
        for (TID tid : vertex_incident_facets[v])
          inc.emplace_back(b[tid]);
        model.AddMaxEquality(v_b[v], inc);
      }
      LinearExpr V_b = LinearExpr::Sum(v_b);

      model.AddEquality(V_b - E_b + F_b, 1);

      // =====================================================================
      //  CONSTRAINT — Per-facet Euler characteristic (V_c − E_c + F_c = 1)
      // =====================================================================

      for (FacetWPtr wf : vertex->facets()) {
        if (FacetSPtr f = wf.lock()) {
          std::vector<BoolVar> b_c;
          for (TID tid=0; tid<triangles.size(); ++tid) {
            if (triangle_to_facet[tid] == f)
              b_c.push_back(b[tid]);
          }
          CGAL_assertion(!b_c.empty());

          LinearExpr F_c = LinearExpr::Sum(b_c);

          // E_c
          std::vector<BoolVar> e_c;
          e_c.reserve(std::pow(vertex->degree(), 2));

          int e = 0;
          for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
              std::vector<TID>& inc_triangles = pid1_and_edges.second;
              std::vector<IntVar> inc;
              for (TID tid : inc_triangles) {
                if (triangle_to_facet[tid] == f)
                  inc.emplace_back(b[tid]);
              }
              if (inc.empty())
                continue;
              e_c.push_back(model.NewBoolVar());
              model.AddMaxEquality(e_c.back(), inc);
            }
          }

          CGAL_assertion(!e_c.empty());
          LinearExpr E_c = LinearExpr::Sum(e_c);

          // V_c
          std::vector<BoolVar> v_c;
          v_c.reserve(std::pow(vertex->degree(), 2));

          for (int v=0; v<points.size(); ++v) {
            std::vector<IntVar> inc;
            std::unordered_set<TID>& inc_triangles = vertex_incident_facets[v];
            for (TID tid : inc_triangles) {
              if (triangle_to_facet[tid] == f)
                inc.emplace_back(b[tid]);
            }
            if (inc.empty())
              continue;
            v_c.push_back(model.NewBoolVar());
            model.AddMaxEquality(v_c.back(), inc);
          }
          CGAL_assertion(!v_c.empty());
          LinearExpr V_c = LinearExpr::Sum(v_c);

          // Only enforce if color is active on the boundary
          model.AddEquality(V_c - E_c + F_c, 1);
        }
      }

      // =====================================================================
      //  CONSTRAINT — Global edge manifoldness
      // =====================================================================

      // @todo we can ignore entirely-boundary edges (pid0/pid1 on boundary **AND** midpoint not on bbox)

      for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
        for (auto& pid1_and_edges : edge_map[pid0]) {
          LinearExpr inc_n = 0;
          std::vector<TID>& inc_triangles = pid1_and_edges.second;
          for (TID tid : inc_triangles)
            inc_n += b[tid];
          model.AddLessOrEqual(inc_n, 2);
        }
      }

#if 0 // this cannot happen because all the facets of a same color live in the same plane
      // =====================================================================
      //  CONSTRAINT — Local edge manifoldness
      // =====================================================================

      for (FacetWPtr wf : vertex->facets()) {
        if (FacetSPtr f = wf.lock()) {
          for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
              LinearExpr inc_n = 0;
              std::vector<TID>& inc_triangles = pid1_and_edges.second;
              for (TID tid : inc_triangles) {
                if (triangle_to_facet[tid] == f)
                  inc_n += b[tid];
              }
              model.AddLessOrEqual(inc_n, 2);
            }
          }
        }
      }
#endif

      // =====================================================================
      //  CONSTRAINT — Global vertex manifoldness
      // =====================================================================

      // @todo

      // =====================================================================
      //  CONSTRAINT — Local vertex manifoldness (per facet color)
      // =====================================================================

      // @fixme seems wrong to me: we can still pinch within the same CC
      // would have to count the number of borders around vertices instead...

#define CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
# if 0 // Flow-based
      for (FacetWPtr wf : vertex->facets()) {
        if (FacetSPtr f = wf.lock()) {
          std::vector<TID> tids_of_f;
          for (TID tid = 0; tid < triangles.size(); ++tid) {
            if (triangle_to_facet[tid] == f)
              tids_of_f.push_back(tid);
          }
          CGAL_assertion(!tids_of_f.empty());

          const int max_flow = static_cast<int>(tids_of_f.size());

          // set the fixed root
          TID root_tid = TID(-1);
          for (TID tid : tids_of_f) {
            const VID bot = face_volume_IDs[tid][0];
            const VID top = face_volume_IDs[tid][1];
            if (in_out_flags[bot] == CC_in_out_flag::INSIDE &&
                in_out_flags[top] == CC_in_out_flag::OUTSIDE) {
              root_tid = tid;
              break;
            }
          }
          CGAL_SS3_TRANSF_TRACE_V(8, "Root of F" << f->id() << " is " << root_tid);
          CGAL_assertion(root_tid != TID(-1));

          std::unordered_map<TID, std::vector<TID>> adj;
          for (std::size_t pid0 = 0; pid0 < points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
              std::vector<TID>& inc = pid1_and_edges.second;
              for (std::size_t i = 0; i < inc.size(); ++i) {
                for (std::size_t j = i + 1; j < inc.size(); ++j) {
                  TID t1 = inc[i], t2 = inc[j];
                  if (triangle_to_facet[t1] == f && triangle_to_facet[t2] == f) {
                    adj[t1].push_back(t2);
                    adj[t2].push_back(t1);
                  }
                }
              }
            }
          }

          // directed flow variables for each undirected edge
          std::unordered_map<TID, std::unordered_map<TID, IntVar>> flow;
          for (TID t1 : tids_of_f) {
            for (TID t2 : adj[t1]) {
              if (t1 < t2) { // create once per edge
                IntVar fwd = model.NewIntVar(Domain(0, max_flow));
                IntVar bwd = model.NewIntVar(Domain(0, max_flow));
                flow[t1][t2] = fwd;
                flow[t2][t1] = bwd;
                // Capacity: no flow if either facet is inactive
                model.AddEquality(fwd, 0).OnlyEnforceIf(b[t1].Not());
                model.AddEquality(fwd, 0).OnlyEnforceIf(b[t2].Not());
                model.AddEquality(bwd, 0).OnlyEnforceIf(b[t1].Not());
                model.AddEquality(bwd, 0).OnlyEnforceIf(b[t2].Not());
              }
            }
          }

          // flow conservation
          LinearExpr b_f; // sum of active facets of this color
          for (TID t : tids_of_f)
            b_f += b[t];

          for (TID t : tids_of_f) {
            LinearExpr inflow, outflow;
            for (TID t2 : adj[t]) {
              inflow  += flow[t2][t]; // flow coming from t2 into t
              outflow += flow[t][t2]; // flow going from t to t2
            }
            if (t == root_tid) {
              // Root supplies (total_active - 1) units
              model.AddEquality(outflow - inflow + 1, b_f);
            } else {
              // Every other node consumes exactly b[t] (1 if active, 0 if not)
              model.AddEquality(inflow - outflow, b[t]);
            }
          }
        }
      }
# else // spanning tree-based
      for (FacetWPtr wf : vertex->facets()) {
        if (FacetSPtr f = wf.lock()) {
          std::vector<TID> tids_of_f;
          for (TID tid = 0; tid < triangles.size(); ++tid) {
            if (triangle_to_facet[tid] == f)
              tids_of_f.push_back(tid);
          }
          CGAL_assertion(!tids_of_f.empty());
          const int N = tids_of_f.size();

          // Pick the fixed root
          TID root_tid = TID(-1);
          for (TID tid : tids_of_f) {
            const VID bot = face_volume_IDs[tid][0];
            const VID top = face_volume_IDs[tid][1];
            // If both cells are fixed and different => b[tid] is fixed true
            const bool bot_fixed = (in_out_flags[bot] == CC_in_out_flag::INSIDE);
            const bool top_fixed = (in_out_flags[top] == CC_in_out_flag::OUTSIDE);
            if (bot_fixed && top_fixed) {
              root_tid = tid;
              break;
            }
          }
          std::cout << "root of F" << f->id() << " is " << root_tid << std::endl;
          CGAL_assertion(root_tid != TID(-1));

          // Distance variables: d[tid] is the distance from the root.
          std::unordered_map<TID, IntVar> d;
          for (TID tid : tids_of_f) {
            d[tid] = model.NewIntVar({0, N}); // Max distance is the number of triangles
          }

          // The root is strictly at distance 0
          model.AddEquality(d[root_tid], 0);

          // Build local adjacency map for this specific color
          std::unordered_map<TID, std::vector<TID>> adj;
          for (std::size_t pid0 = 0; pid0 < points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
              std::vector<TID>& inc = pid1_and_edges.second;
              for (std::size_t i = 0; i < inc.size(); ++i) {
                for (std::size_t j = i + 1; j < inc.size(); ++j) {
                  TID t1 = inc[i], t2 = inc[j];
                  if (triangle_to_facet[t1] == f && triangle_to_facet[t2] == f) {
                    adj[t1].push_back(t2);
                    adj[t2].push_back(t1);
                  }
                }
              }
            }
          }

          for (TID tid : tids_of_f) {
            if (tid == root_tid)
              continue;

            std::vector<BoolVar> parent_vars;

            for (TID nbr : adj[tid]) {
              BoolVar is_parent = model.NewBoolVar();
              parent_vars.push_back(is_parent);

              // If 'nbr' is the parent, 'nbr' MUST also be active on the boundary
              model.AddImplication(is_parent, b[nbr]);
              // If 'nbr' is the parent, distance strictly increases by 1
              model.AddEquality(d[tid], d[nbr] + 1).OnlyEnforceIf(is_parent);
            }

            // Connectivity Constraint:
            // If this facet is on the boundary (b[tid] == 1), it must have EXACTLY 1 parent.
            // If it is not on the boundary (b[tid] == 0), it has 0 parents.
            model.AddEquality(LinearExpr::Sum(parent_vars), b[tid]);
          }
        }
      }
# endif // flow or tree
#endif

      // =====================================================================
      //  SOLVE
      // =====================================================================

      std::vector<CC_in_out_flag> solution(volume_CCs.size(), CC_in_out_flag::UNINITIALIZED);

#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
      const auto& proto = model.Build();
      CGAL_SS3_TRANSF_TRACE_V(8, "Variables: " << proto.variables_size());
      CGAL_SS3_TRANSF_TRACE_V(8, "Constraints: " << proto.constraints_size());

      SatParameters params;
      params.set_stop_after_first_solution(true);
      params.set_enumerate_all_solutions(false);

      std::vector<bool> current_b_values(triangles.size());

      bool ok = false;
      const CpSolverResponse response = SolveWithParameters(proto, params);
      if (response.status() == CpSolverStatus::OPTIMAL ||
          response.status() == CpSolverStatus::FEASIBLE) {
        ok = true;

        for (int c = 0; c < C; ++c) {
          solution[c] = SolutionBooleanValue(response, x[c]) ? CC_in_out_flag::INSIDE
                                                             : CC_in_out_flag::OUTSIDE;
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(32, "OK is " << ok);
      if (!ok) {
        CGAL_SS3_TRANSF_TRACE_V(1, "ERROR: failed to find a solution [" << response.status() << "]");
        std::abort();
      }
#else
      // Since vertex manifoldness is difficult to express with constraints, try and try
      // till we succeed
      bool valid_solution_found = false;
      int iteration = -1;
      const int max_iterations = 1000;

      while (iteration < max_iterations) {
        ++iteration;
        std::cout << "--- Starting solver iteration " << iteration << " ---";

        const auto& proto = model.Build();
        CGAL_SS3_TRANSF_TRACE_V(8, "Variables: " << proto.variables_size());
        CGAL_SS3_TRANSF_TRACE_V(8, "Constraints: " << proto.constraints_size());

        SatParameters params;
        params.set_stop_after_first_solution(true);
        params.set_enumerate_all_solutions(false);

        std::vector<bool> current_b_values(triangles.size());

        bool ok = false;
        const CpSolverResponse response = SolveWithParameters(proto, params);
        if (response.status() == CpSolverStatus::OPTIMAL ||
            response.status() == CpSolverStatus::FEASIBLE) {
          ok = true;

          for (int c = 0; c < C; ++c) {
            solution[c] = SolutionBooleanValue(response, x[c]) ? CC_in_out_flag::INSIDE
                                                               : CC_in_out_flag::OUTSIDE;
          }
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "OK is " << ok);
        if (!ok) {
          CGAL_SS3_TRANSF_TRACE_V(1, "ERROR: failed to find a solution [" << response.status() << "]");
          std::abort();
        }

#ifdef CGAL_SS3_DUMP_FILES
        // Tentative dump
        {
          unsigned int undetermined_n = 0;
          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            CGAL_assertion(solution[i] != CC_in_out_flag::UNINITIALIZED);
            if (solution[i] == CC_in_out_flag::INSIDE) {
              CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known inside (tentative)");
            } else if (solution[i] == CC_in_out_flag::OUTSIDE) {
              CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known outside (tentative)");
            } else {
              ++undetermined_n;
            }
          }
          CGAL_SS3_TRANSF_TRACE_V(64, undetermined_n << " undetermined cells");

          std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
            {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
          for (auto e : dumps) {
            std::vector<Point_3> cc_points = points;
            std::vector<std::vector<PID> > cc_triangles;

            for (std::size_t i=0; i<volume_CCs.size(); ++i) {
              CGAL_assertion(e.first != CC_in_out_flag::TBD);
              CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

              if (solution[i] != e.first)
                continue;

              for (TID tid : volume_CCs[i])
                cc_triangles.push_back(triangles[tid]);

              std::ostringstream oss;
              oss << "results/volumes_" << e.second << "_tentative.off";
              CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
              std::ostringstream oss_n;
              oss_n << "results/volumes_" << e.second << "_tentative_normalized.off";
              CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_triangles, CGAL::parameters::stream_precision(17));
            }
          }
        }
#endif

        // Now, check for vertex manifoldness both in the global and local scale
        PID nm_vertex_id = -1;
        std::optional<FacetSPtr> failed_f = std::nullopt; // track which facet failed

        auto check_CC = [&](std::optional<FacetSPtr> of = std::nullopt) -> bool
        {
          std::vector<std::vector<PID> > local_triangles;
          for (std::size_t i=0; i<triangles.size(); ++i) {
            if (!triangle_to_facet[i])
              continue;

            if (solution[face_volume_IDs[i][0]] == solution[face_volume_IDs[i][1]])
              continue;

            if (of.has_value() && triangle_to_facet[i] != of.value())
              continue;

            local_triangles.push_back(triangles[i]);
          }

          CGAL_assertion(!local_triangles.empty());

          typedef std::vector<PID>                                            PointRange;
          typedef std::vector<std::vector<PID> >                              PolygonRange;
          typedef Polygon_mesh_processing::internal::Polygon_soup_orienter<PointRange, PolygonRange>   Orienter;

          typename Orienter::Edge_map edges(points.size());
          typename Orienter::Marked_edges marked_edges;
          Orienter::fill_edge_map(edges, marked_edges, local_triangles);
          CGAL_assertion(marked_edges.empty()); // no NM edges is part of the constraints

          Orienter::has_singular_vertices(points.size(), local_triangles, edges, marked_edges, nm_vertex_id);

          return (nm_vertex_id == -1);
        };

        if (!check_CC()) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Warning: issue with global boundary");
        } else {
          for (FacetWPtr wf : vertex->facets()) {
            if (FacetSPtr f = wf.lock()) {
              if (!check_CC(f)) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Warning: issue with CC of F" << f->id());
                failed_f = f;
                break; // nm_vertex_id is populated, and we know exactly which color caused it
              }
            }
          }
        }

        if (nm_vertex_id != PID(-1)) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Non-manifold vertex " << nm_vertex_id << " found");
          std::cout << "at position " << normalized_points[nm_vertex_id] << std::endl;

          for (TID tid : vertex_incident_facets[nm_vertex_id]) {
            if (!triangle_to_facet[tid])
              continue;
            VID bot = face_volume_IDs[tid][0], top = face_volume_IDs[tid][1];
            std::cout << "tid " << tid << " bot_flag=" << (int)in_out_flags[bot]
                      << " top_flag=" << (int)in_out_flags[top]
                      << " b=" << (solution[bot] != solution[top]) << "\n";
          }

          std::vector<BoolVar> nogood_terms;
          bool some_not_fixed = false;

          for (TID tid : vertex_incident_facets[nm_vertex_id]) {
            if (!triangle_to_facet[tid])
              continue;

            if (failed_f.has_value() && triangle_to_facet[tid] != failed_f.value()) {
              continue;
            }

            VID bot = face_volume_IDs[tid][0];
            VID top = face_volume_IDs[tid][1];
            bool bot_fixed = (in_out_flags[bot] == CC_in_out_flag::INSIDE);
            bool top_fixed = (in_out_flags[top] == CC_in_out_flag::OUTSIDE);
            if (!bot_fixed || !top_fixed) {
              some_not_fixed = true; // At least one cell can be flipped!
            }

            // Reconstruct whether this facet was on the boundary in the rejected solution
            bool is_boundary = (solution[bot] != solution[top]);
            nogood_terms.push_back(is_boundary ? b[tid].Not() : b[tid]);
          }

          // Safety abort if the input geometry strictly forces a non-manifold pinch
          if (!some_not_fixed) {
             CGAL_SS3_TRANSF_TRACE_V(1, "Error: Non-manifold vertex " << nm_vertex_id
                << " is completely surrounded by fixed input cells.");
             std::abort();
          }

          if (!nogood_terms.empty()) {
            model.AddBoolOr(nogood_terms);
          }
        } else {
          CGAL_SS3_TRANSF_TRACE_V(1, "Valid solution found");
          valid_solution_found = true;
          break;
        }
      }

      CGAL_assertion(valid_solution_found);
#endif // CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS


#ifdef CGAL_SS3_DUMP_FILES
      // Final dump
      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(solution[i] != CC_in_out_flag::UNINITIALIZED);
          if (solution[i] == CC_in_out_flag::INSIDE) {
            CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known inside (final)");
          } else if (solution[i] == CC_in_out_flag::OUTSIDE) {
            CGAL_SS3_TRANSF_TRACE_V(64, "volume " << i << " is known outside (final)");
          } else {
            ++undetermined_n;
          }
        }
        CGAL_SS3_TRANSF_TRACE_V(64, undetermined_n << " undetermined cells");

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps) {
          std::vector<Point_3> cc_points = points;
          std::vector<std::vector<PID> > cc_triangles;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            CGAL_assertion(e.first != CC_in_out_flag::TBD);
            CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

            if (solution[i] != e.first)
              continue;

            for (TID tid : volume_CCs[i])
              cc_triangles.push_back(triangles[tid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_final.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_final_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_triangles, CGAL::parameters::stream_precision(17));
          }
        }
      }
#else // CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

#endif // CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

      // a valid partition has:
      // - C1) all cells are either inside or outside
      // - C2) boundary is manifold
      // - C3) only 1 "inside" face-connected component & 1 "outside" face-CC
      // - C4) only 1 simply connected component per input face
      auto is_valid_partition = [&](const std::vector<std::vector<TID> >& volume_CCs,
                                    const std::vector<std::array<VID, 2>>& face_volume_IDs,
                                    const std::vector<CC_in_out_flag>& in_out_flags,
                                    const std::vector<std::vector<PID> >& triangles,
                                    const std::vector<Point_3>& points,
                                    const std::vector<std::unordered_map<PID, std::vector<TID> > >& edge_map) -> bool
      {
        namespace PMP = CGAL::Polygon_mesh_processing;

        // C1
        for(const auto& flag : in_out_flags) {
          if(flag == CC_in_out_flag::UNINITIALIZED || flag == CC_in_out_flag::TBD) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: failed C1");
            return false;
          }
        }

        // C2+C3+C4
        auto check_CC = [&](std::optional<FacetSPtr> of = std::nullopt) -> bool {
          std::vector<std::vector<PID> > local_triangles;
          for (std::size_t i=0; i<triangles.size(); ++i) {
            if (!triangle_to_facet[i])
              continue;

            if (in_out_flags[face_volume_IDs[i][0]] == in_out_flags[face_volume_IDs[i][1]])
              continue;

            if (of.has_value() && triangle_to_facet[i] != of.value())
              continue;

            local_triangles.push_back(triangles[i]);
          }

          CGAL_assertion(!local_triangles.empty());

          if (!PMP::is_polygon_soup_a_polygon_mesh(local_triangles)) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: failed C34-PS_PM");
#ifdef CGAL_SS3_DUMP_FILES
            CGAL::IO::write_polygon_soup("results/bad_CC_normalized.off", normalized_points, local_triangles, CGAL::parameters::stream_precision(17));
#endif
            return false;
          }

          // @todo avoid using a Surface_mesh
          CGAL::Surface_mesh<Point_3> sm;
          PMP::polygon_soup_to_polygon_mesh(points, local_triangles, sm);
          CGAL_assertion(faces(sm).size() == local_triangles.size());

          const unsigned int nb = number_of_borders(sm);
          if (nb != 1) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: failed C34-borders (" << nb << ")");
            return false;
          }

          if (PMP::does_self_intersect(sm)) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: CC self intersects");
            return false;
          }

          return true;
        };

        if (!check_CC()) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Error: issue with global boundary");
          return false;
        }

        for (FacetWPtr wf : vertex->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (!check_CC(f)) {
              CGAL_SS3_TRANSF_TRACE_V(1, "Error: issue with CC of F" << f->id());
              return false;
            }
          }
        }

        return true;
      };

      bool is_valid = is_valid_partition(volume_CCs, face_volume_IDs, solution, triangles, points, edge_map);
      CGAL_SS3_TRANSF_TRACE_V(32, "valid partition? " << is_valid);
      if (!is_valid)
        std::abort();

      // Now, extract the connectivity from the arrangement and plug it back into the polyhedron itself
      std::vector<std::vector<PID> > boundary_triangles;
      std::vector<FacetSPtr> boundary_triangle_to_facet;
      for (std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        if (solution[face_volume_IDs[i][0]] == solution[face_volume_IDs[i][1]])
          continue;

        boundary_triangles.push_back(triangles[i]);
        boundary_triangle_to_facet.push_back(triangle_to_facet[i]);
      }

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_polygon_soup("results/boundary.off", points, boundary_triangles);
#endif

      CGAL_assertion(!boundary_triangles.empty());
      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(boundary_triangles));

      // @todo avoid using a CGAL::Surface_mesh
      using Mesh = CGAL::Surface_mesh<Point_3>;
      using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

      Mesh bsm;
      auto ifpm = get(CGAL::dynamic_face_property_t<FacetSPtr>{}, bsm);

      auto out = boost::make_function_output_iterator(
        [ifpm, &boundary_triangle_to_facet](const std::pair<std::size_t, face_descriptor>& pid_f) {
          // we shouldn't have triangles that are not from a face here
          CGAL_assertion(boundary_triangle_to_facet[pid_f.first] != FacetSPtr());
          put(ifpm, pid_f.second, boundary_triangle_to_facet[pid_f.first]);
        });

      PMP::polygon_soup_to_polygon_mesh(points, boundary_triangles, bsm,
                                        CGAL::parameters::polygon_to_face_output_iterator(out));
      CGAL_assertion(faces(bsm).size() == boundary_triangles.size());

      // ordered facets in the polyhedron
      CGAL::unordered_flat_map<FacetSPtr, int> local_facet_indices;
      edge = start_edge;
      FacetSPtr current_facet = edge->left(vertex);

      int curr_index = 0;
      do {
        local_facet_indices[current_facet] = curr_index++;
        edge = edge->next(vertex);
        current_facet = edge->other(current_facet);
      } while (edge != start_edge);

      using vec2i = boost::shared_array<int>;
      using combi = std::vector<vec2i>;

      auto reconstruct_combi = [&]() -> combi {
        auto create_split = [](int begin, int end) -> vec2i {
          vec2i result(new int[2]);
          result[0] = begin;
          result[1] = end;
          return result;
        };

        std::set<std::pair<int, int> > unique_splits;
        for (edge_descriptor e : edges(bsm))
        {
          if (is_border(e, bsm))
            continue;

          halfedge_descriptor h = halfedge(e, bsm);
          face_descriptor f1 = face(h, bsm);
          face_descriptor f2 = face(opposite(h, bsm), bsm);

          // map mesh faces -> original input facets
          // (use the type produced by ifpm; here assumed comparable + indexable)
          FacetSPtr input_facet_1 = get(ifpm, f1);
          FacetSPtr input_facet_2 = get(ifpm, f2);
          CGAL_assertion(input_facet_1 != FacetSPtr() && input_facet_2 != FacetSPtr());

          if (input_facet_1 == input_facet_2)
            continue;

          // if the facets are neighbors, it's not a split edge
          EdgeSPtr common_e = input_facet_1->find_edge(input_facet_2);
          if (common_e != EdgeSPtr())
            continue;

          // map original facets -> fan positions
          int a = local_facet_indices.at(input_facet_1);
          int b = local_facet_indices.at(input_facet_2);

          // canonicalize so begin < end (matches create_single_split_combinations)
          int begin = (std::min)(a, b);
          int end = (std::max)(a, b);

          CGAL_SS3_TRANSF_TRACE_V(64, "create split between F" << input_facet_1->id() << " and F" << input_facet_2->id());
          unique_splits.emplace(begin, end);
        }

        combi result;
        result.reserve(unique_splits.size());
        for (const std::pair<int, int>& s : unique_splits)
          result.push_back(create_split(s.first, s.second));

        // sort into the canonical order used by the generator
        // compare_splits returns +1 when split1 < split2 lexicographically.
        auto sorter = [](const vec2i& s1, const vec2i& s2)
        {
          auto compare_splits = [](const vec2i& split1, const vec2i& split2) -> int {
            int result = 0;
            if (split1[0] < split2[0] || (split1[0] == split2[0] && split1[1] < split2[1]))
            result = 1;
            else if (split1[0] > split2[0] || (split1[0] == split2[0] && split1[1] > split2[1]))
            result = -1;
            return result;
          };

          return (compare_splits(s1, s2) > 0);
        };

        std::sort(result.begin(), result.end(), sorter);

        // a fully-split deg N vertex requires N-3 splits
        CGAL_SS3_TRANSF_TRACE_V(64, "result size: " << result.size());
        CGAL_SS3_TRANSF_TRACE_V(64, "degree: " << vertex->degree());
        CGAL_postcondition(result.size() == vertex->degree() - 3);

        return result;
      };

      // @todo just split_vertex(vertex, split_combi)...?
      using Combi_vertex_splitter = algorithm::Combi_vertex_splitter<GeomTraits>;

      combi split_combi = reconstruct_combi();
      PolyhedronSPtr poly_c = Combi_vertex_splitter::copy_vertex(vertex);
      VertexSPtr vertex_c = poly_c->vertices().front();
      Combi_vertex_splitter::split_vertex(vertex_c, split_combi);
      Combi_vertex_splitter::apply(poly_c, vertex);
      CGAL_postcondition(polyhedron->is_consistent());
    }

    CGAL_postcondition_code(bool success =)
      Transformation::reset_points(polyhedron);

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/perturbed_V4.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

    CGAL_postcondition(success);
    CGAL_postcondition(polyhedron && polyhedron->is_consistent());
    CGAL_postcondition(!Self_intersection::has_self_intersecting_surface(polyhedron));
  }
#endif // CGAL_SPS3_USE_V4_PERTURBATION

  // Perturbation to ensure generic configuration.
  // We always need to ensure that points are exactly on the planes of their incident facets.
  static void apply_rand_perturbation(PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Applying random perturbation to the polyhedron...");

    const FT sq_tolerance = 0;
    for (const VertexSPtr vertex : polyhedron->vertices()) {
      is_stable_vertex(vertex, polyhedron, sq_tolerance);
    }

    // if the input is all triangles, simply perturb points directly
#ifndef CGAL_SPS3_USE_V4_PERTURBATION
    // don't do this because we don't split vertices with V4
    if (is_triangle_polyhedron(polyhedron))
      return rand_move_points(polyhedron);
#endif

    // Generic approach
    Transformation::normalize_facet_planes(polyhedron);

    ConfigurationSPtr config = Configuration::get_instance();
    bool safe_mode = config->get_Boolean("Preprocessing", "check_degenerate_configuration");

    PolyhedronSPtr p_mem = polyhedron->clone();

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
    apply_rand_plane_tilts_V4(polyhedron);
#else
    apply_rand_plane_tilts_V3(polyhedron);
#endif

    if (safe_mode) {
      CGAL_SS3_TRANSF_TRACE_V(16, "Safe mode is enabled, checking validity of the perturbation...");

      for (;;) {
        // IO::write_OBJ("results/last_perturbation.obj", polyhedron, parameters::stream_precision(17).do_not_triangulate_faces(true));

        if (check_perturbed_positions_proximity(polyhedron, p_mem) &&
            do_all_plane_pairs_intersect(polyhedron) &&
            !Self_intersection::has_self_intersecting_surface(polyhedron) &&
            do_all_plane_triplets_intersect(polyhedron)) {
          CGAL_SS3_TRANSF_TRACE_V(16, "Safe mode is enabled, checking validity of the perturbation...");
          break;
        }

        CGAL_SS3_TRANSF_TRACE_V(4, "Perturbation failed, retrying...");

        polyhedron = p_mem->clone();
        if (!is_triangle_polyhedron(polyhedron)) {
          Transformation::triangulate_facets(polyhedron);
          p_mem = polyhedron->clone();
        }

        rand_move_points(polyhedron);
      }
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/perturbed.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H */
