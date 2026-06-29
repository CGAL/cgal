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
    unsigned int s = rd();
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(min, max);

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
          CGAL_SS3_TRANSF_TRACE_V(32, "Degenerate facet pair:");
          CGAL_SS3_TRANSF_TRACE_V(32, "  " << facet1->to_string());
          CGAL_SS3_TRANSF_TRACE_V(32, "  " << facet2->to_string());
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
            CGAL_SS3_TRANSF_TRACE("Degenerate facet triplet:");
            CGAL_SS3_TRANSF_TRACE("  " << facet1->to_string());
            CGAL_SS3_TRANSF_TRACE("  " << facet2->to_string());
            CGAL_SS3_TRANSF_TRACE("  " << facet3->to_string());
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
      double rx = CGAL::to_double(p.x()) + v_r[0];
      double ry = CGAL::to_double(p.y()) + v_r[1];
      double rz = CGAL::to_double(p.z()) + v_r[2];
      vertex->set_point(Point_3{rx, ry, rz});
    }

    // recompute normalized planes to ensure points are on the supporting planes
    polyhedron->init_planes();
    Transformation::normalize_facet_planes(polyhedron);
    CGAL_assertion_code(bool success =)
      Transformation::reset_points(polyhedron);
    CGAL_assertion(success);
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
          if (limit_sq < local_sq_max)
            local_sq_max = limit_sq;
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
        std::cout << "perturb F" << f->id() << " with scale " << scale << std::endl;

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

      if (all_ok)
      {
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
          std::cerr << "no intersection between plane and bbox?!" << std::endl;
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
          std::cerr << "plane/bbox intersection is not a polygon" << std::endl;
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
          std::cerr << "no intersection between plane and bbox" << std::endl;
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
          std::cerr << "plane/bbox intersection is not a polygon" << std::endl;
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
        std::cout << "center = " << center << std::endl;
        std::cout << "prev_pt = " << prev_pt << std::endl;
        std::cout << "next_pt = " << next_pt << std::endl;
        std::cout << "center + n = " << center + n << std::endl;
        std::cout << "angle_gt_180 = " << angle_gt_180 << std::endl;

        for(std::size_t i=0; i<local_range.size(); ++i) {
          std::size_t i1 = base_idx + i;
          std::size_t i2 = base_idx + ((i+1)%local_range.size());

          // Compute midpoint of the two extremities
          const Point_3& p1 = points[i1];
          const Point_3& p2 = points[i2];
          Point_3 mid = CGAL::midpoint(p1, p2);
          std::cout << "test: " << mid << std::endl;

          // Test if mid is between the two orientation planes
          bool on_pos_side_prev = (plane_prev.oriented_side(mid) == CGAL::ON_NEGATIVE_SIDE);
          bool on_pos_side_next = (plane_next.oriented_side(mid) == CGAL::ON_POSITIVE_SIDE);
          std::cout << "on_pos_side_prev = " << on_pos_side_prev << std::endl;
          std::cout << "on_pos_side_next = " << on_pos_side_next << std::endl;

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

    // The approach here is to do a subtle plan perturbation, and then to reconstruct a solid
    // As long as the perturbation is small enough and the input is manifold,
    // the only topological change that happens is a split of vertices, resulting in
    // a degree-3-vertices-only, perturbed polyhedron.
    //
    // At each input high-degree vertex, we must recover a valid split from the arrangement
    // of the perturbed planes.

    // Generic approach
    Transformation::normalize_facet_planes(polyhedron);

    ConfigurationSPtr config = Configuration::get_instance();
    double range = config->get_double("Preprocessing", "perturbation_epsilon");

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
    perturbPlanesCoefficientsControlled(polyhedron, range, 1. /*max displacement*/);
#else
    // @todo apply a garanteed small-enough nudge
    for (const FacetSPtr& facet : polyhedron->facets()) {
      CGAL_SS3_TRANSF_TRACE_V(32, "Nudge F" << facet->id());
      perturbPlaneCoefficientsNudge(facet, range);
    }
#endif

    std::list<VertexSPtr> vertices_tosplit;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex->degree() > 3) {
        vertices_tosplit.push_back(vertex);
      }
    }

    CGAL_SS3_CORE_TRACE_V(4, vertices_tosplit.size() << " vertices to split");

    using Arr_vertex_splitter = algorithm::Arr_vertex_splitter<GeomTraits>;
    using Arr_vertex_splitter_sptr = std::shared_ptr<Arr_vertex_splitter>;
    Arr_vertex_splitter_sptr vertex_splitter = Arr_vertex_splitter::create();

    for (const VertexSPtr& vertex : vertices_tosplit) {
      CGAL_SS3_CORE_TRACE_V(8, "Splitting vertex at " << vertex->point());

      vertex->sort();

      // soup to be used in the arrangement
      std::vector<Point_3> points;
      std::vector<std::vector<std::size_t> > triangles;
      std::vector<FacetSPtr> triangle_to_facet;

      // Compute bounding box of intersections
      Iso_cuboid_3 bbox = Arr_vertex_splitter::compute_intersection_bbox(vertex);

      // Get the triangles from the base planes, with tagged faces
      get_clipped_plane_faces(vertex, bbox, points, triangles, triangle_to_facet);
      std::cout << points.size() << " points, " << triangles.size() << " triangles [base]" << std::endl;

      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > facet_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            facet_to_triangles[fsptr].push_back(triangles[i]);
          }
        }
        std::size_t idx = 0;
        for (const auto& kv : facet_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_base_" << idx << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
          ++idx;
        }
      }

      // Add the bounding box faces
#if 1
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

      std::cout << points.size() << " points, " << triangles.size() << " triangles [w/ bbox]" << std::endl;
#endif

      for (std::size_t i=0; i<triangle_to_facet.size(); ++i) {
        std::cout << "triangle " << i << " ==> ptr: " << triangle_to_facet[i] << std::endl;
      }

      PMP::merge_duplicate_points_in_polygon_soup(points, triangles);
      CGAL::IO::write_OFF("results/arr_soup.off", points, triangles, CGAL::parameters::stream_precision(17));

      CGAL_assertion(triangles.size() == triangle_to_facet.size());

      std::cout << "autorefining..." << std::endl;
      std::vector<FacetSPtr> updated_triangle_to_facet;
      Range_updating_autoref_visitor<FacetSPtr> autoref_visitor(triangle_to_facet, updated_triangle_to_facet);

      PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::visitor(autoref_visitor));
      triangle_to_facet = std::move(updated_triangle_to_facet);
      CGAL_assertion(triangles.size() == triangle_to_facet.size());

      CGAL::IO::write_OFF("results/arr_autorefined.off", points, triangles, CGAL::parameters::stream_precision(17));

      for (std::size_t i=0; i<triangles.size(); ++i) {
        std::cout << "[4] triangle " << i << " ptr: " << triangle_to_facet[i] << std::endl;
      }

      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > facet_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            facet_to_triangles[fsptr].push_back(triangles[i]);
          }
        }
        std::size_t idx = 0;
        for (const auto& kv : facet_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_autorefined_" << idx << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
          ++idx;
        }
      }

      std::cout << "repairing..." << std::endl;
      Range_updating_repair_PS_visitor<FacetSPtr> repair_ps_visitor(triangle_to_facet);
      PMP::merge_duplicate_polygons_in_polygon_soup(points, triangles,
                                                    CGAL::parameters::visitor(repair_ps_visitor)
                                                                    .erase_all_duplicates(false) /*keep one*/
                                                                    .require_same_orientation(false)
                                                                    .verbose(true));
      CGAL::IO::write_OFF("results/arr_repaired.off", points, triangles, CGAL::parameters::stream_precision(17));

      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        std::map<FacetSPtr, std::vector<std::vector<std::size_t> > > facet_to_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          FacetSPtr fsptr = triangle_to_facet[i];
          if (fsptr) {
            facet_to_triangles[fsptr].push_back(triangles[i]);
          }
        }
        std::size_t idx = 0;
        for (const auto& kv : facet_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_repaired_" << idx << ".off";
          CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
          ++idx;
        }
      }

      CGAL_assertion(triangles.size() == triangle_to_facet.size());

      // dump only the triangles with a non nullptr
      {
        std::vector<std::vector<std::size_t> > input_triangles;
        for (std::size_t i=0; i<triangles.size(); ++i) {
          if (triangle_to_facet[i]) {
            input_triangles.push_back(triangles[i]);
          }
        }

        std::cout << input_triangles.size() << " input triangles" << std::endl;
        CGAL::IO::write_OFF("results/arr_recovered_input.off", points, input_triangles, CGAL::parameters::stream_precision(17));
      }

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

            std::cout << "processing a non-manifold edge: " << points[pid0] << " -- " << points[pid1] << std::endl;

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
        std::cout << "Building volume #" << CC_ID << " from seed face " << seed_tid << std::endl;

        volume_CCs.emplace_back();

        std::stack<std::pair<TID, bool> > to_visit;
        to_visit.emplace(seed_tid, start_from_inverted_face);

        while (!to_visit.empty())
        {
          TID current_tid;
          bool invert_face;
          std::tie(current_tid, invert_face) = to_visit.top();
          to_visit.pop();

          std::cout << "At face " << current_tid << " [" << triangles[current_tid][0]
                                                 << ", " << triangles[current_tid][1]
                                                 << ", " << triangles[current_tid][2] << "], ";
          std::cout << "invert: " << invert_face << ", ";
          std::cout << "VIDS: " << face_volume_IDs[current_tid][0] << " " << face_volume_IDs[current_tid][1] << std::endl;

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

            std::cout << "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]" << std::endl;
            std::cout << "    pos: " << points[e_pids.first] << " " << points[e_pids.second] << std::endl;

            // The faces are ordered CCW while looking from pid0.
            // So the walking while looking from [j] depends on whether [j] is pid0 or not
            int iter_direction = (e_pids.first == triangles[current_tid][j]) ? 1 : -1;
            std::cout << "    iter_direction = " << iter_direction << std::endl;

            // and it also depends on whether we are walking above or below the face
            iter_direction *= invert_face ? 1 : -1;
            std::cout << "    invert_face = " << invert_face << std::endl;

            TID next_tid = current_tid;
            for (;;) {
              if (inc_triangles.size() == 1) {
                std::cerr << "Warning: dangling triangle..." << std::endl;
                std::cout << "    over the edge, the triangle is ITSELF " << current_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                to_visit.emplace(current_tid, !invert_face);
                break;
              } else if (inc_triangles.size() == 2) {
                // we should only be there once, meaning if we do not ignore orientations,
                // then the faces MUST be compatible
                CGAL_assertion(next_tid == current_tid);

                next_tid = (inc_triangles[0] == current_tid) ? inc_triangles[1] : inc_triangles[0];
                std::cout << "    over the edge, the triangle is TRIVIALLY " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
                CGAL_assertion(next_tid != current_tid);
              } else {
                // tricky part, now
                auto tid_it = std::find(std::begin(inc_triangles), std::end(inc_triangles), next_tid /*updates on every iteration*/);
                CGAL_assertion(tid_it != inc_triangles.end());

                if (iter_direction == 1) { // CCW
                  std::cout << "    CCW walk" << std::endl;
                  auto next_it = std::next(tid_it);
                  next_tid = (next_it == inc_triangles.end()) ? inc_triangles[0] : *next_it;
                } else { // CW
                  std::cout << "    CW walk" << std::endl;
                  next_tid = (tid_it == inc_triangles.begin()) ? inc_triangles.back() : *(std::prev(tid_it));
                }

                std::cout << "    over the edge, the triangle is " << next_tid << " [" << triangles[next_tid][0] << ", " << triangles[next_tid][1] << ", " << triangles[next_tid][2] << "], ";
                std::cout << "VIDS " << face_volume_IDs[next_tid][0] << " " << face_volume_IDs[next_tid][1] << std::endl;
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
              std::cout << "    flipping? " << flip_side
                        << " (N: " << triangles[next_tid][(pos+1)%3] << " C: " << triangles[current_tid][(j+1)%3] << ")" << std::endl;

              std::cout << "Final TID = " << next_tid << std::endl;
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

      std::cout << volume_CCs.size() << " volume CCs" << std::endl;

      for (std::size_t i=0; i<volume_CCs.size(); ++i) {
        // build a mesh from the soup
        std::vector<Point_3> cc_points = points;
        std::vector<std::vector<PID> > cc_triangles;
        for (TID tid : volume_CCs[i]) {
          cc_triangles.push_back(triangles[tid]);
        }

        std::ostringstream oss;
        oss << "results/volume_cc_" << i << ".off";
        CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
        std::cout << "Wrote volume CC " << i << " with " << cc_triangles.size() << " triangles." << std::endl;
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

      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(in_out_flags[i] != CC_in_out_flag::UNINITIALIZED);
          if (in_out_flags[i] == CC_in_out_flag::INSIDE)
            std::cout << "volume " << i << " is known inside (base)" << std::endl;
          else if (in_out_flags[i] == CC_in_out_flag::OUTSIDE)
            std::cout << "volume " << i << " is known outside (base)" << std::endl;
          else
            ++undetermined_n;
        }
        std::cout << undetermined_n << " undetermined cells" << std::endl;

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps)
        {
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
          }
        }
      }

      // Set up the non obvious known volumes

      boost::dynamic_bitset<> is_boundary_point(points.size(), 0);
      for (std::size_t pid=0; pid<points.size(); ++pid)
      {
        const Point_3& p = points[pid];
        if (p.x() == bbox.xmin() || p.x() == bbox.xmax() ||
            p.y() == bbox.ymin() || p.y() == bbox.ymax() ||
            p.z() == bbox.zmin() || p.z() == bbox.zmax())
          is_boundary_point.set(pid);
      }

#if 0
      // For each edge of the mesh, we have 2 pairs of 4 volume cells incident to the bounding box
      // We keep the four volume cells that are on the same side as the edge's endpoint opposite
      // of the vertex.
      // In these four volume cells, assign inside/outside: if the edge is convex, then it's 3
      // volumes outside and 1 volume inside. If the edge is concave, then it's 1 volume inside
      // outside and 3 volumes inside.

      for(PID pid0 = 0; pid0 < points.size(); ++pid0) {
        for(const auto& [pid1, tids] : edge_map[pid0]) {
          CGAL_assertion(tids.size() <= 4);
          if (tids.size() != 4)
            continue;
          if (is_boundary_point[pid0] == is_boundary_point[pid1])  // internal or on border edge
            continue;

          // now, all incident volumes will be on the boundary since the edge is on the boundary
          std::set<FacetSPtr> orig_facets;
          for(TID tid : tids) {
            CGAL_assertion(face_volume_IDs[tid][0] != VID(-1) && face_volume_IDs[tid][1] != VID(-1));
            orig_facets.insert(triangle_to_facet[tid]);
          }
          CGAL_assertion(orig_facets.size() == 2);

          FacetSPtr f0 = *(orig_facets.begin());
          FacetSPtr f1 = *(std::next(orig_facets.begin(), 1));
          EdgeSPtr common_e = f0->find_edge(f1);
          if (common_e == EdgeSPtr()) // these faces are not adjacent in the input
            continue;

          // get the point that is *not* the center of the star
          VertexSPtr other_v = common_e->other(vertex);

          Vector_3 orig_v { vertex->point(), other_v->point() };

          // sort by distance to the center vertex to get the oriented edge
          // @fixme sounds like an inconsistency because we are going to check scalar products
          // between vectors that do NOT correspond to post-perturbation state
          bool flip = (CGAL::compare_distance(vertex->point(), points[pid1], points[pid0]) == CGAL::SMALLER);
          PID spid0 = flip ? pid1 : pid0;
          PID spid1 = flip ? pid0 : pid1;

          Vector_3 new_v { points[spid0], points[spid1] };

          std::cout << "--" << std::endl;
          std::cout << points[pid0] << " " << points[pid1] << " is a tentative edge" << std::endl;
          std::cout << points[spid0] << " " << points[spid1] << " (oriented)" << std::endl;
          std::cout << "vertex = " << vertex->point() << std::endl;
          std::cout << "other = " << other_v->point() << std::endl;
          std::cout << CGAL::scalar_product(orig_v, new_v) << std::endl;

          if (CGAL::scalar_product(orig_v, new_v) < 0)
            continue;

          std::cout << points[pid0] << " " << points[pid1] << " is a good edge" << std::endl;

          // Now, we have an edge touching the bbox (but not on the bbox), equivalent to (part of)
          // the edge in the initial star
          // -> mark the volumes as inside and outside depending on convexity

          // @fixme same inconsistency here: planes have been tilted and so we use the new plane
          // but third_v->point() is an OLD position...
          VertexSPtr third_v = *(f1->vertices().begin());
          const Plane_3& pl0 = f0->get_plane();
          for (VertexSPtr v1 : f1->vertices()) {
            if (CGAL::compare_distance(pl0, v1->point(), third_v->point()) == CGAL::LARGER)
              third_v = v1;
          }

          std::map<VID, unsigned int> positive_hits;
          for(TID tid : tids)
            ++(positive_hits[face_volume_IDs[tid][1]]);

          bool is_convex = (f0->get_plane().oriented_side(third_v->point()) == CGAL::NEGATIVE);
          std::cout << "it is a " << (is_convex ? "convex" : "concave") << " edge" << std::endl;
          std::cout << "f0f1: " << f0->id() << " " << f1->id() << std::endl;
          Vector_3 n = f0->get_plane().orthogonal_vector();
          n = n/CGAL::approximate_sqrt(n.squared_length());
          std::cout << "normal " << vertex->point() << " " << vertex->point() + n << std::endl;
          std::cout << "third point was " << third_v->point() << std::endl;

          for (const auto& e : positive_hits) {
            if (e.second == 0)
              in_out_flags[e.first] = CC_in_out_flag::INSIDE;
            else if (e.second == 2)
              in_out_flags[e.first] = CC_in_out_flag::OUTSIDE;
            else
              in_out_flags[e.first] = (is_convex ? CC_in_out_flag::OUTSIDE : CC_in_out_flag::INSIDE);
          }
        }
      }
#else
      // Generalization of above: for the complete trace of the star on the bounding box,
      // we know the cell above is OUTSIDE and the cell beneath is INSIDE

      // We need to have the border incident edges for any type of boundary point
      std::unordered_map<PID, std::unordered_map<PID, std::vector<TID> > > symmetrical_border_edge_map;
      for (TID ti=0; ti<triangles.size(); ++ti) {
        for (std::size_t j=0; j<3; ++j) {
          PID pid0 = triangles[ti][j];
          PID pid1 = triangles[ti][(j+1)%3];
          if (!is_boundary_point[pid0] || !is_boundary_point[pid1])
            continue;
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
      auto mark_cells = [&](PID ipid0, PID ipid1, FacetSPtr current_facet)
      {
        auto [pid0, pid1] = CGAL::make_sorted_pair(ipid0, ipid1);
        auto edge_it = symmetrical_border_edge_map.find(pid0);
        CGAL_assertion(edge_it != symmetrical_border_edge_map.end());
        auto pid1_it = edge_it->second.find(pid1);
        CGAL_assertion(pid1_it != edge_it->second.end());

        const auto& tids = pid1_it->second;
        for (TID tid : tids) {
          if (triangle_to_facet[tid] != current_facet)
            continue;

          in_out_flags[face_volume_IDs[tid][1]] = CC_in_out_flag::OUTSIDE;
          in_out_flags[face_volume_IDs[tid][0]] = CC_in_out_flag::INSIDE;

          return;
        }
      };

      // Helper to find next boundary point on current_facet in the direction of next_facet
      auto find_next_boundary_point = [&](PID prev_pid,
                                          PID current_pid,
                                          FacetSPtr prev_facet,
                                          FacetSPtr current_facet,
                                          FacetSPtr next_facet) -> std::pair<PID, bool>
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

          if (on_current_facet)
            incident_boundary_points.insert(other_pid);
        }

        std::cout << incident_boundary_points.size() << " incident boundary points" << std::endl;
        CGAL_assertion(incident_boundary_points.size() == 2);

        PID first_pid = *(incident_boundary_points.begin());
        PID second_pid = *(std::next(incident_boundary_points.begin()));

        if (prev_pid != PID(-1))
          return {(prev_pid == first_pid ? second_pid : first_pid), true};

        // If we are here, candidate_pid is the point that is at the intersection
        // of prev_facet and facet and we need to know in which direction to walk
        // along facet. This is deduced from the convexity of the edge in the input polyhedron.
        //
        // @fixme this assumes that perturbation does not modify the convexity of edges...

        const Plane_3& pl = current_facet->get_plane();

        // determine convexity of prev/current
        VertexSPtr third_v = *(current_facet->vertices().begin());
        const Plane_3& prev_pl = prev_facet->get_plane();
        for (VertexSPtr v : current_facet->vertices()) {
          if (CGAL::compare_distance(prev_pl, v->point(), third_v->point()) == CGAL::LARGER)
            third_v = v;
        }

        bool is_convex = (prev_facet->get_plane().oriented_side(third_v->point()) == CGAL::NEGATIVE);

        // now, if the edge is convex in the input polyhedron, we must turn "right", meaning,
        // the next point is on the negative side of the plane too
        if (prev_facet->get_plane().oriented_side(points[first_pid]) == CGAL::NEGATIVE)
          return {(is_convex ? first_pid : second_pid), true};
        else
          return {(is_convex ? second_pid : first_pid), true};

        return {PID(-1), false};
      };

      // Walk the trace of the star on the bbox ("box link")
      FacetSPtr start_facet = vertex->first_facet();
      FacetSPtr facet = start_facet;

      std::ofstream link_out("results/link.txt");
      link_out.precision(17);

      do
      {
        FacetSPtr prev_facet = facet->prev(vertex);
        std::cout << "walking on " << facet->id() << " (prev: " << prev_facet->id() << ")" << std::endl;

        PID current_pid = PID(-1);

      // Find an arrangement edge that:
      // 1. Has prev_facet and facet as the two incident input faces
      // 2. Has one endpoint on the boundary, one not
      // 3. Satisfies the orientation check

        bool found_start_point = false;

        for (PID pid0 = 0; pid0 < points.size(); ++pid0) {
          // (non-border) edge map here because we seek a non-border edge
          for (const auto& [pid1, tids] : edge_map[pid0]) {
            if (tids.size() != 4)
              continue;

            if (is_boundary_point[pid0] == is_boundary_point[pid1]) // internal or on border edge
              continue;

            std::set<FacetSPtr> incident_faces;
            for (TID tid : tids)
              incident_faces.insert(triangle_to_facet[tid]);
            CGAL_assertion(incident_faces.size() == 2);

            if (!incident_faces.count(prev_facet) || !incident_faces.count(facet))
              continue;

            std::cout << points[pid0] << " " << points[pid1] << " is a tentative edge" << std::endl;

            // Found a matching edge, now check orientation
            EdgeSPtr common_e = prev_facet->find_edge(facet);
            CGAL_assertion(common_e != EdgeSPtr());

            VertexSPtr other_v = common_e->other(vertex);
            Vector_3 orig_v{vertex->point(), other_v->point()};

            bool flip = (CGAL::compare_distance(vertex->point(), points[pid1], points[pid0]) == CGAL::SMALLER);
            PID test_spid0 = flip ? pid1 : pid0;
            PID test_spid1 = flip ? pid0 : pid1;

            Vector_3 new_v { points[test_spid0], points[test_spid1] };

            if (CGAL::scalar_product(orig_v, new_v) >= 0) {
              found_start_point = true;
              current_pid = test_spid1;
              break;
            }
          }

          if (found_start_point)
            break;
        }

        CGAL_assertion(found_start_point);
        std::cout << "Start at " << current_pid << " pos: " << points[current_pid] << std::endl;

        // Now walk the star link on arrangement edges

        FacetSPtr next_facet = facet->next(vertex);
        const Plane_3& next_pl = next_facet->get_plane();

        PID prev_pid (-1);
        for(;;) {
          CGAL_assertion(!next_pl.has_on(points[current_pid]));

          // find the next point on 'facet' in the direction of 'next_facet'
          auto [next_pid, valid] = find_next_boundary_point(prev_pid, current_pid, prev_facet, facet, next_facet);
          CGAL_assertion(valid);
          CGAL_assertion(is_boundary_point[current_pid] && is_boundary_point[next_pid]);

          link_out << "2 " << points[current_pid] << " " << points[next_pid] << std::endl;

          mark_cells(current_pid, next_pid, facet);

          prev_pid = current_pid;
          current_pid = next_pid;

          // Check if we've reached the boundary between current_facet and next_facet
          if (next_pl.has_on(points[current_pid]))
            break;
        }

        facet = next_facet;
      }
      while (facet != start_facet);

      std::cout << "Known cells marked from star link walk" << std::endl;
#endif

      // Intermediate dump
      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(in_out_flags[i] != CC_in_out_flag::UNINITIALIZED);
          if (in_out_flags[i] == CC_in_out_flag::INSIDE)
            std::cout << "volume " << i << " is known inside (boundary edges)" << std::endl;
          else if (in_out_flags[i] == CC_in_out_flag::OUTSIDE)
            std::cout << "volume " << i << " is known outside (boundary edges)" << std::endl;
          else
            ++undetermined_n;
        }
        std::cout << undetermined_n << " undetermined cells" << std::endl;

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps)
        {
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
          }
        }
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

      // cells that are known, are known
      for (int c = 0; c < C; ++c) {
        if (in_out_flags[c] == CC_in_out_flag::OUTSIDE)
          model.FixVariable(x[c], false);
        else if (in_out_flags[c] == CC_in_out_flag::INSIDE)
          model.FixVariable(x[c], true);
      }

      // The boundary must be well oriented, so the boundary can only have an INSIDE volume
      // on its negative side, and an OUTSIDE volume on its positive side.
      // In model terms, this is:
      //   x[pos] => x[neg]. (if above is true, below is true)
      for(std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        model.AddImplication(x[face_volume_IDs[i][1]], x[face_volume_IDs[i][0]]);
      }

      // missing constraints:
      // @todo manifoldness

      // @speed could use hints using intersection of volumes?

      // Variables for the facets
      std::vector<BoolVar> b(triangles.size());
      for (std::size_t tid=0; tid<triangles.size(); ++tid) {
        b[tid] = model.NewBoolVar();
        const VID bot_vid = face_volume_IDs[tid][0];
        const VID top_vid = face_volume_IDs[tid][1];

        if (bot_vid == VID(-1) || top_vid == VID(-1)) {
          model.FixVariable(b[tid], false);
        } else {
          model.AddNotEqual(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[tid]);
          model.AddEquality(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[tid].Not());
        }
      }

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
      e_b.reserve(std::pow(vertex->degree(), 3));

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

      // =====================================================================
      //  CONSTRAINT — Global vertex manifoldness
      // =====================================================================

      // check manifoldness on the border: only 2 boundary edges

      // @todo

      // =====================================================================
      //  CONSTRAINT — Local vertex manifoldness
      // =====================================================================


      // @todo




      // Solve
      const auto& proto = model.Build();
      LOG(INFO) << "Variables: " << proto.variables_size();
      LOG(INFO) << "Constraints: " << proto.constraints_size();

      SatParameters params;
      params.set_stop_after_first_solution(true);
      params.set_enumerate_all_solutions(false);

      bool ok = false;
      const CpSolverResponse response = SolveWithParameters(proto, params);
      if (response.status() == CpSolverStatus::OPTIMAL ||
          response.status() == CpSolverStatus::FEASIBLE) {
        ok = true;

        for (int c = 0; c < C; ++c)
          in_out_flags[c] = SolutionBooleanValue(response, x[c]) ? CC_in_out_flag::INSIDE
                                                                 : CC_in_out_flag::OUTSIDE;
      }

      std::cout << "OK is " << ok << std::endl;

      // Final dump
      {
        unsigned int undetermined_n = 0;
        for (std::size_t i=0; i<volume_CCs.size(); ++i) {
          CGAL_assertion(in_out_flags[i] != CC_in_out_flag::UNINITIALIZED);
          if (in_out_flags[i] == CC_in_out_flag::INSIDE)
            std::cout << "volume " << i << " is known inside (final)" << std::endl;
          else if (in_out_flags[i] == CC_in_out_flag::OUTSIDE)
            std::cout << "volume " << i << " is known outside (final)" << std::endl;
          else
            ++undetermined_n;
        }
        std::cout << undetermined_n << " undetermined cells" << std::endl;

        std::vector<std::pair<CC_in_out_flag, std::string> > dumps =
          {{CC_in_out_flag::INSIDE, "INSIDE"}, {CC_in_out_flag::OUTSIDE, "OUTSIDE"}};
        for (auto e : dumps)
        {
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
            oss << "results/volumes_" << e.second << "_final.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_triangles, CGAL::parameters::stream_precision(17));
          }
        }
      }

      // a valid partition has:
      // - C1) all cells are either inside or outside
      // - C2) boundary is manifold
      // - C3) only 1 "inside" face-connected component & 1 "outside" face-CC
      // - C4) only 1 simply connected component per input face
      auto is_valid_partition = [&](const std::vector<std::vector<TID>>& volume_CCs,
                                   const std::vector<std::array<VID, 2>>& face_volume_IDs,
                                   const std::vector<CC_in_out_flag>& in_out_flags,
                                   const std::vector<std::vector<PID>>& triangles,
                                   const std::vector<Point_3>& points,
                                   const std::vector<std::unordered_map<PID, std::vector<TID>>>& edge_map) -> bool
      {
        namespace PMP = CGAL::Polygon_mesh_processing;

        // C1
        for(const auto& flag : in_out_flags) {
          if(flag == CC_in_out_flag::UNINITIALIZED || flag == CC_in_out_flag::TBD) {
            std::cerr << "Error: failed C1" << std::endl;
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
            std::cerr << "Error: failed C34-PS_PM" << std::endl;
            CGAL::IO::write_polygon_soup("results/bad_CC.off", points, local_triangles, CGAL::parameters::stream_precision(17));
            return false;
          }

          // @todo avoid using a Surface_mesh
          CGAL::Surface_mesh<Point_3> sm;
          PMP::polygon_soup_to_polygon_mesh(points, local_triangles, sm);
          CGAL_assertion(faces(sm).size() == local_triangles.size());

          const unsigned int nb = number_of_borders(sm);
          if (nb != 1) {
            std::cerr << "Error: failed C34-borders (" << nb << ")" << std::endl;
            return false;
          }

          if (PMP::does_self_intersect(sm)) {
            std::cerr << "Error: CC self intersects" << std::endl;
            return false;
          }

          return true;
        };

        if (!check_CC()) {
          std::cerr << "Error: issue with global boundary" << std::endl;
          return false;
        }

        for (FacetWPtr wf : vertex->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (!check_CC(f)) {
              std::cerr << "Error: issue with CC of F" << f->id() << std::endl;
              return false;
            }
          }
        }

        return true;
      };

      bool is_valid = is_valid_partition(volume_CCs, face_volume_IDs, in_out_flags, triangles, points, edge_map);
      std::cout << "valid partition? " << is_valid << std::endl;
      if (!is_valid)
        std::exit(1);

      // Now, extract the connectivity from the arrangement and plug it back into the polyhedron itself
      std::vector<std::vector<PID> > boundary_triangles;
      std::vector<FacetSPtr> boundary_triangle_to_facet;
      for (std::size_t i=0; i<triangles.size(); ++i) {
        if (!triangle_to_facet[i])
          continue;

        if (in_out_flags[face_volume_IDs[i][0]] == in_out_flags[face_volume_IDs[i][1]])
          continue;

        boundary_triangles.push_back(triangles[i]);
        boundary_triangle_to_facet.push_back(triangle_to_facet[i]);
      }

      CGAL::IO::write_polygon_soup("results/boundary.off", points, boundary_triangles);

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
      EdgeSPtr edge = vertex->find_edge(start_facet);
      EdgeSPtr edge_first = edge;
      FacetSPtr current_facet = start_facet;

      int curr_index = 0;
      do {
        local_facet_indices[current_facet] = curr_index++;
        edge = edge->next(vertex);
        current_facet = edge->other(current_facet);
      } while (edge != edge_first);

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

          std::cout << "create split between F" << input_facet_1->id() << " and F" << input_facet_2->id() << std::endl;
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
        std::cout << "result size: " << result.size() << std::endl;
        std::cout << "degree: " << vertex->degree() << std::endl;
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

    CGAL_assertion_code(bool success =)
      Transformation::reset_points(polyhedron);
    CGAL_assertion(success);
    CGAL_postcondition(polyhedron && polyhedron->is_consistent());

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/perturbed_V4.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif
  }
#endif // CGAL_SPS3_USE_V4_PERTURBATION

  // Perturbation to ensure generic configuration.
  // We always need to ensure that points are exactly on the planes of their incident facets.
  static void apply_rand_perturbation(PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Applying random perturbation to the polyhedron...");

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
