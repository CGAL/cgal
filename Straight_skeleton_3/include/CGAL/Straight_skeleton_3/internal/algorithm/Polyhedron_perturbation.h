// Copyright (c) 2024-2026 GeometryFactory (France)
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
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/enum.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
# include "ortools/sat/cp_model.h"
# include "ortools/sat/sat_parameters.pb.h"
#endif

#include <algorithm>
#ifdef CGAL_SS3_DUMP_FILES
# include <fstream>
#endif
#include <iterator>
#include <limits>
#include <list>
#include <random>

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

  struct Stability_failure
  {
    VertexSPtr v;
    FacetSPtr f0;
    FacetSPtr f1;
    FacetSPtr f2;
    FT sq_distance;
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
    CGAL_SS3_TRANSF_TRACE_V(4, "Check all 2-, 3-combinations of planes");

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
    * checks for degeneracies:
    *  - all pairs of planes should intersect in a non-degenerate line.
    *  - all triplets of planes should intersect in a point.
    * avoiding the O(n^3) naive complexity
    */
  static bool are_planes_in_general_position(const PolyhedronSPtr& polyhedron,
                                             std::vector<Stability_failure>* failures = nullptr)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Check if planes are in general position...");
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    CGAL_precondition(polyhedron->facets().size() >= 3);

    // Three planes meet in a single point iff their normals are linearly independent.
    // With c_ij = n_i x n_j, we have det(n_i, n_j, n_k) = 0 iff c_ij and c_ik are parallel,
    // and thus c_ij = 0 detects parallel normals directly.

    // Two canonicalized projective normals are the same direction iff proportional with a positive factor
    struct Projective_normal { FT x, y, z; };
    using Direction = std::pair<Projective_normal, std::size_t>;

    // canonical := first non-zero coordinate is positive
    auto canonicalize = [](const Vector_3& c) -> Projective_normal
    {
      const FT& cx = c.x();
      const FT& cy = c.y();
      const FT& cz = c.z();

      if (cx < 0 || (cx == 0 && cy < 0) || (cx == 0 && cy == 0 && cz < 0)) {
        return { -cx, -cy, -cz };
      } else {
        return { cx, cy, cz };
      }
    };

    // Strict weak ordering on sign-canonicalized non-zero triples.
    //
    // Precondition (what makes the cross-multiplications below valid): the
    // first non-zero coordinate of each input is > 0, so every denominator the
    // ratios are taken over is strictly positive and multiplying by it
    // preserves order.
    //
    // Category = position of the first non-zero coordinate. Categories are
    // disjoint and projectively separated, so ordering by category first is
    // consistent; within a category the order is lexicographic on the
    // coordinate ratios, compared without division.
    auto projective_compare = [](const Projective_normal& a,
                                 const Projective_normal& b) -> bool
    {
      auto get_category = [](const Projective_normal& v) -> int {
        if (v.x > 0) return 0;
        if (v.y > 0) return 1;
        return 2;
      };

      const int cat_a = get_category(a);
      const int cat_b = get_category(b);

      if (cat_a != cat_b) {
        return cat_a < cat_b;
      }

      if (cat_a == 0) {
        // Both x > 0, so compare y_a/x_a with y_b/x_b, then z_a/x_a with z_b/x_b.
        const FT ya_xb = a.y * b.x;
        const FT yb_xa = b.y * a.x;
        if (ya_xb != yb_xa) {
          return ya_xb < yb_xa;
        }
        return a.z * b.x < b.z * a.x;
      }

      if (cat_a == 1) {
        // Both x == 0, y > 0, so compare z_a/y_a with z_b/y_b.
        return a.z * b.y < b.z * a.y;
      }

      // both are (0, 0, z) with z > 0, hence identical directions.
      return false;
    };

    // the comparator refuses both orders <=> same projective direction
    auto projective_equal = [&projective_compare](const Projective_normal& a,
                                                  const Projective_normal& b) -> bool
    {
      return (!projective_compare(a, b) && !projective_compare(b, a));
    };

    std::vector<std::pair<FacetSPtr, Vector_3> > normals;
    const std::size_t nf = polyhedron->facets().size();
    normals.reserve(nf);

    for (const FacetSPtr& facet : polyhedron->facets()) {
      const Vector_3 n = facet->get_plane().orthogonal_vector();
      normals.emplace_back(facet, n);
    }

    std::vector<Direction> directions;
    directions.reserve(nf - 1);

    for (std::size_t i=0; i<nf; ++i) {
      directions.clear();
      const Vector_3& n1 = normals[i].second;
      for (std::size_t j=i+1; j<nf; ++j) {
        const Vector_3& n2 = normals[j].second;
        const Vector_3 c_ij = CGAL::cross_product(n1, n2);
        if (c_ij == CGAL::NULL_VECTOR) { // Parallel normals: *every* triplet containing both facets is degenerate
          std::size_t k = 0;
          while (k == i || k == j) {
            ++k;
          }
          CGAL_SS3_TRANSF_TRACE_V(1, "Degenerate facet triplet: " << normals[i].first->id() << " " << normals[j].first->id() << " -1");
          if (failures) {
            failures->push_back({VertexSPtr(), normals[i].first, normals[j].first, normals[k].first, FT(-1)});
            continue;
          } else {
            return false;
          }
        }

        directions.push_back(std::make_pair(canonicalize(c_ij), j));
      }

      // Duplicate direction <=> c_ij || c_ik <=> det(n_i, n_j, n_k) == 0.
      std::sort(directions.begin(), directions.end(),
                [&projective_compare](const Direction& d1, const Direction& d2) -> bool {
                  return projective_compare(d1.first, d2.first);
                });

      const auto dup = std::adjacent_find(
          directions.begin(), directions.end(),
          [&projective_equal](const Direction& d1, const Direction& d2) -> bool {
            return projective_equal(d1.first, d2.first);
          });

      if (dup != directions.end()) {
        const std::size_t jidx = dup->second;
        const std::size_t kidx = (std::next(dup))->second;
        CGAL_SS3_TRANSF_TRACE_V(1, "Degenerate facet triplet: " << normals[i].first->id() << " " << normals[jidx].first->id() << " " << normals[kidx].first->id());
        if (failures) {
          failures->push_back({VertexSPtr(), normals[i].first, normals[jidx].first, normals[kidx].first, FT(-2)});
          continue;
        } else {
          return false;
        }
      }
    }

    if (failures) {
      return failures->empty();
    } else {
      return true;
    }
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
      FT rx = p.x() + FT(v_r[0]);
      FT ry = p.y() + FT(v_r[1]);
      FT rz = p.z() + FT(v_r[2]);
      Point_3 new_pos{rx, ry, rz};
      CGAL_SS3_TRANSF_TRACE_V(32, "Point from " << vertex->point() << " to " << new_pos);
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

    CGAL_SS3_TRANSF_TRACE_V(32, "Perturb (Nudge) Facet " << facet->id());
    CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                        << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    auto nudge = [&](const FT& v)
    {
      static std::random_device rd;
      unsigned int s = 0; // rd();
      // CGAL_SS3_TRANSF_TRACE("seed = " << s);
      static std::mt19937 gen(s);
      static std::uniform_real_distribution<> rdist(-range, range);

      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // The point is also that once 'nv' is a double, its interval will be a singleton,
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

    CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                      << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    CGAL_postcondition(Transformation::has_normalized_plane(facet));
  }

  /**
  * nudges the plane coefficients but ensure that the perturbed plane goes through 0, 1, or 2 fixed points.
  * If 0 points: nudge all coefficients independently.
  * If 1 point: nudge (a, b, c), recompute d so the plane passes through the point.
  * If 2 points: nudge (a, b, c) with the constraint that the new plane passes through both points.
  */
  template <typename VertexRange>
  static void perturbPlaneCoefficientsFixedPoints(const FacetSPtr& facet,
                                                  const double range,
                                                  const VertexRange& fixed_vertices)
  {
    CGAL_SS3_TRANSF_TRACE_V(32, "Perturb (Fixed) F" << facet->id() << " [" << facet->vertices().size() << " vs]");
    CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                        << facet->get_plane().c() << " " << facet->get_plane().d() << "]");
    CGAL_SS3_TRANSF_TRACE_V(32, "  with " << fixed_vertices.size() << " fixed vertices");
    for (const VertexSPtr& fixed_vertex : fixed_vertices) {
      CGAL_SS3_TRANSF_TRACE_V(32, "    fixed V" << fixed_vertex->id() << " at " << fixed_vertex->point());
    }

    CGAL_precondition(Transformation::has_normalized_plane(facet));

    boost::container::small_vector<const Point_3*, 3> fixed_points;
    for (const VertexSPtr& fixed_vertex : fixed_vertices) {
      fixed_points.push_back(&(fixed_vertex->point()));
    }

    const std::size_t fixed_points_count = fixed_points.size();
    CGAL_assertion(fixed_points_count < 3);

    static std::random_device rd;
    unsigned int s = 0; // rd()
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(-range, range);

    auto nudge = [&](const FT& v)
    {
      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // The point is also that once 'nv' is a double, its interval will be a singleton,
      // and we will have access to static filters
      double step = rdist(gen);
      double nv = CGAL::to_double(v) + step;
      return nv;
    };

    // 0 fixed points: nudge (a, b, c), recompute d so the plane passes through the facet centroid
    if (fixed_points_count == 0) {
      double cx = 0.;
      double cy = 0.;
      double cz = 0.;
      unsigned int vertex_count = 0;
      for (const VertexSPtr& vertex : facet->vertices()) {
        const Point_3& p = vertex->point();
        cx += CGAL::to_double(p.x());
        cy += CGAL::to_double(p.y());
        cz += CGAL::to_double(p.z());
        ++vertex_count;
      }
      CGAL_assertion(vertex_count > 0);
      const double x0 = cx / double(vertex_count);
      const double y0 = cy / double(vertex_count);
      const double z0 = cz / double(vertex_count);

      const double na = nudge(facet->get_plane().a());
      const double nb = nudge(facet->get_plane().b());
      const double nc = nudge(facet->get_plane().c());
      const double nd = -(na * x0 + nb * y0 + nc * z0);

      facet->set_plane(Plane_3{na, nb, nc, nd});
    } else if (fixed_points_count == 1) {
      // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
      const Point_3& p0 = *(fixed_points[0]);

      const double na = nudge(facet->get_plane().a());
      const double nb = nudge(facet->get_plane().b());
      const double nc = nudge(facet->get_plane().c());

      const FT& x0 = p0.x();
      const FT& y0 = p0.y();
      const FT& z0 = p0.z();
      FT d = - (na * x0 + nb * y0 + nc * z0);
      facet->set_plane(Plane_3{na, nb, nc, d});
      CGAL_postcondition(facet->get_plane().has_on(p0));
    } else if (fixed_points_count == 2) {
#if 1
      // The planes containing p0 and p1 form a 2-dimensional linear subspace (a pencil)
      // of the 4-dimensional coefficient space (a,b,c,d).
      // We build an exact (division-free) basis {(N1,D1), (N2,D2)} of that pencil, and set
      //  (a1, b1, c1, d1) = alpha * (N1,D1) + beta * (N2,D2)
      // giving a constrained, nudged plane close to the current one.

      const Point_3& p0 = *(fixed_points[0]);
      const Point_3& p1 = *(fixed_points[1]);
      CGAL_assertion(p0 != p1);

      auto dot3 = [](const double a[3], const double b[3]) -> double {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      };
      auto norm3 = [&dot3](const double a[3]) -> double {
        return std::sqrt(dot3(a, a));
      };

      const FT ux = p1.x() - p0.x();
      const FT uy = p1.y() - p0.y();
      const FT uz = p1.z() - p0.z();

      const double u_d[3] = { CGAL::to_double(ux),
                              CGAL::to_double(uy),
                              CGAL::to_double(uz) };
      const double au[3] = { std::abs(u_d[0]), std::abs(u_d[1]), std::abs(u_d[2]) };
      const double ul = norm3(u_d);
      CGAL_assertion(ul > 0.); // means we don't have p0 and p1 absurdly close

      // Pivot on the largest component of u.
      // Then |N1 x N2| = |u_pivot| * |u|, and the angle between N1 and N2 is >= 60 degrees
      // Indeed, if a = |u_pivot| is the largest and b, c are the others,
      // |cos| = bc / (sqrt(b^2+a^2) * sqrt(c^2+a^2)) <= bc / (b*sqrt(2) * c*sqrt(2)) = 1/2).
      // The 2x2 solve below is therefore well-conditioned.
      int pivot = 0;
      if (au[1] >= au[0] && au[1] >= au[2]) {
        pivot = 1;
      } else if (au[2] >= au[0] && au[2] >= au[1]) {
        pivot = 2;
      }

      // Compute the basis of the pencil
      // D_i = -(N_i . p0) is written expanded such that every product has magnitude ~|p0|*|u|
      // (no |p0|^2 terms): this avoids the cancellation of the Pluecker form p0 x p1
      // and preserves interval-filter margin.
      FT N1x, N1y, N1z, D1;
      FT N2x, N2y, N2z, D2;

      switch (pivot) {
        case 0: // |ux| largest
          N1x = uy;
          N1y = - ux;
          N1z = FT(0);
          N2x = uz;
          N2y = FT(0);
          N2z = - ux;
          D1 = ux * p0.y() - uy * p0.x();
          D2 = ux * p0.z() - uz * p0.x();
          break;
        case 1: // |uy| largest
          N1x = - uy;
          N1y =  ux;
          N1z = FT(0);
          N2x = FT(0);
          N2y =  uz;
          N2z = - uy;
          D1 = uy * p0.x() - ux * p0.y();
          D2 = uy * p0.z() - uz * p0.y();
          break;
        default: // |uz| largest
          N1x = - uz;
          N1y = FT(0);
          N1z =  ux;
          N2x = FT(0);
          N2y = - uz;
          N2z =  uy;
          D1 = uz * p0.x() - ux * p0.z();
          D2 = uz * p0.y() - uy * p0.z();
          break;
      }

      // select a plane in the pencil

      const double n1[3] = { CGAL::to_double(N1x), CGAL::to_double(N1y), CGAL::to_double(N1z) };
      const double n2[3] = { CGAL::to_double(N2x), CGAL::to_double(N2y), CGAL::to_double(N2z) };
      const double l1 = norm3(n1);
      const double l2 = norm3(n2);
      CGAL_assertion(l1 > 0. && l2 > 0.);
      const double e1[3] = { n1[0]/l1, n1[1]/l1, n1[2]/l1 };
      const double e2[3] = { n2[0]/l2, n2[1]/l2, n2[2]/l2 };

      // Target direction: cast to double to drop the DAG of the previous plane
      double t[3] = { CGAL::to_double(facet->get_plane().a()),
                      CGAL::to_double(facet->get_plane().b()),
                      CGAL::to_double(facet->get_plane().c()) };
      const double tn = norm3(t);
      CGAL_assertion(tn > 0.);
      t[0] /= tn;
      t[1] /= tn;
      t[2] /= tn;

      // Only directions orthogonal to u are reachable.
      // |t.uhat| = sin(phi), phi = the angle by which the two anchors force the normal to swing.
      const double uh[3] = { u_d[0]/ul, u_d[1]/ul, u_d[2]/ul };
      const double tu = dot3(t, uh);
      double tp[3] = { t[0] - tu*uh[0], t[1] - tu*uh[1], t[2] - tu*uh[2] };
      const double tpn = norm3(tp); // == cos(phi)
      double phi = std::atan2(std::abs(tu), tpn);

      // if the new facet's normal is (nearly) parallel to the anchor edge direction,
      // then the resulting plane through p0,p1 is not close to the current one.
      const double threshold = 1e-8;
      const bool degenerate = (tpn <= threshold);
      if (degenerate) {
        // Return a valid plane, and let other mechanisms downstream deal with this bad perturbation...
        CGAL_SS3_TRANSF_TRACE_V(32, "  Warning: F" << facet->id() << ": normal nearly parallel to the anchor edge (|cos| = " << tpn << ")");
        tp[0] = e1[0];
        tp[1] = e1[1];
        tp[2] = e1[2];
        phi = 0.5 * CGAL_PI;
      } else {
        tp[0] /= tpn;
        tp[1] /= tpn;
        tp[2] /= tpn;
      }

      // Rotate around 'u' by a random angle
      using rdist_param = std::uniform_real_distribution<>::param_type;
      const double theta = rdist(gen, rdist_param(-range, range));
      const double ct = std::cos(theta), st = std::sin(theta);

      const double w[3] = { uh[1]*tp[2] - uh[2]*tp[1],
                            uh[2]*tp[0] - uh[0]*tp[2],
                            uh[0]*tp[1] - uh[1]*tp[0] };

      const double tj[3] = { ct*tp[0] + st*w[0],
                             ct*tp[1] + st*w[1],
                             ct*tp[2] + st*w[2] };

      // Express 'tj' in the (unit) basis (e1,e2) of u^perp
      const double g = dot3(e1, e2);
      const double r1 = dot3(e1, tj);
      const double r2 = dot3(e2, tj);
      const double det = 1.0 - g*g;
      CGAL_assertion(det > 0.5);
      const double alpha = (r1 - g*r2) / det / l1; // coefficient of N1
      const double beta  = (r2 - g*r1) / det / l2; // coefficient of N2

      // Round alpha,beta onto a common dyadic exponent: one shared power-of-two
      // denominator instead of two unrelated ones, both numerators <= 2^52.
      const double amax = std::max(std::abs(alpha), std::abs(beta));
      CGAL_assertion(amax > 0.);
      int aexp = 0;
      std::frexp(amax, &aexp); // amax in [2^(aexp-1), 2^aexp)
      const int shift = 52 - aexp;
      const double scl     = std::ldexp(1.0, shift);
      const double inv_scl = std::ldexp(1.0, -shift);
      CGAL_assertion(scl > 0. && std::isfinite(scl) && inv_scl > 0.);

      const double lam_d = std::round(alpha * scl); // integral, |.| <= 2^52
      const double mu_d  = std::round(beta  * scl);
      CGAL_assertion(lam_d != 0. || mu_d != 0.);

      const FT inv_scl_ft(inv_scl); // exact power of two
      const FT alpha_ft = FT(lam_d) * inv_scl_ft;
      const FT beta_ft  = FT(mu_d)  * inv_scl_ft;

      // Back to the exact world
      const FT a1 = alpha_ft * N1x + beta_ft * N2x;
      const FT b1 = alpha_ft * N1y + beta_ft * N2y;
      const FT c1 = alpha_ft * N1z + beta_ft * N2z;
      const FT d1 = alpha_ft * D1  + beta_ft * D2;

      facet->set_plane(Plane_3{a1, b1, c1, d1});

      // By construction, for any 'alpha' and 'beta':
      //   a1*p0 + d1 = alpha*(N1.p0 + D1) + beta*(N2.p0 + D2) = 0
      //   a1*p1 + d1 = alpha*(N1.u)       + beta*(N2.u)       = 0
      CGAL_postcondition(facet->get_plane().has_on(p0));
      CGAL_postcondition(facet->get_plane().has_on(p1));

      // (a1,b1,c1) reproduces the UNIT vector tj, so the plane is already
      // normalized to ~1 ulp: no rescaling (and no extra 53-bit factor) needed.
      // The new normal is within theta + phi of the old one.
      const double nd[3] = { CGAL::to_double(a1),
                             CGAL::to_double(b1),
                             CGAL::to_double(c1) };
      const double ndn = norm3(nd);
      const double nu[3] = { nd[0]/ndn, nd[1]/ndn, nd[2]/ndn };
      const double cosang = dot3(nu, t);
      CGAL_assertion(ndn > 0.);
      CGAL_assertion(std::abs(ndn - 1.0) < 1e-9);
      CGAL_assertion(degenerate || cosang >= std::cos(std::abs(theta) + phi) - 1e-9);
#else
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
      const FT ux = p1x - p0x;
      const FT uy = p1y - p0y;
      const FT uz = p1z - p0z;
      const FT uu = ux*ux + uy*uy + uz*uz;

      // Original normal
      const FT& a0 = facet->get_plane().a();
      const FT& b0 = facet->get_plane().b();
      const FT& c0 = facet->get_plane().c();

      // Project original normal onto plane orthogonal to u
      const FT dot = a0*ux + b0*uy + c0*uz;
      const FT scale = dot / uu;
      const FT ab = a0 - scale * ux;
      const FT bb = b0 - scale * uy;
      const FT cb = c0 - scale * uz;

      // Find a direction to nudge (cross product)
      // v = u x B = u x (N - (dot / uu) u) = u x N - (dot / uu) u x u = u x N
      // const FT vx = uy * cb - uz * bb;
      // const FT vy = uz * ab - ux * cb;
      // const FT vz = ux * bb - uy * ab;
      const FT vx = uy * c0 - uz * b0;
      const FT vy = uz * a0 - ux * c0;
      const FT vz = ux * b0 - uy * a0;

      // Nudge the normal
      const double epsilon = rdist(gen);
      const FT a1 = ab + epsilon * vx;
      const FT b1 = bb + epsilon * vy;
      const FT c1 = cb + epsilon * vz;
      const FT d1 = - (a1 * p0x + b1 * p0y + c1 * p0z);
      facet->set_plane(Plane_3{a1, b1, c1, d1});
      CGAL_postcondition(facet->get_plane().has_on(p0));
      CGAL_postcondition(facet->get_plane().has_on(p1));
#endif
    }

    CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                      << facet->get_plane().c() << " " << facet->get_plane().d() << "]");

    CGAL_postcondition(Transformation::has_normalized_plane(facet));
  }

  static void perturbPlaneCoefficientsHighDegrees(const FacetSPtr& facet,
                                                  const double range)
  {
    std::set<VertexSPtr> high_degree_vertices;
    for (const VertexSPtr& v : facet->vertices()) {
      if (v->degree() > 3) {
        high_degree_vertices.insert(v);
      }
    }

    return perturbPlaneCoefficientsFixedPoints(facet, range, high_degree_vertices);
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

    for (const FacetSPtr& facet : polyhedron->facets()) {
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

      CGAL::unordered_flat_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<CGAL::unordered_flat_map<PCDT_FH, bool>> in_domain(in_domain_map);
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
        // Give priority to the one with the least amount of determined vertices
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

    // sometimes we could fix as a polygon, but we need triangulate for other reasons
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

      // Here, the facet has now received enough determined vertices to become fixed.
      // So, fix it: compute its random perturbation, and add the facet ID to its vertices.

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "F" << f->id() << " is now fixed by");
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : fixing_vertices[f]))
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->id() << " [measure=" << Size_shenanigans::length(fv->point()) << "]");
      CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

      if (f->is_triangle()) {
        CGAL_assertion(fixing_vertices[f].size() == 3); // just to be clear

        // for triangles, all vertices are determined, and there is nothing to nudge
        // (note that vertices were themselves nudged, thus the facet is nudged).
        f->init_plane();
        Transformation::normalize_facet_plane(f);

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_3.OFF", f);
#endif

        // Here we do not need to add the fixed facet to incident determined vertices
        // because all vertices are already fully determined
        return true;
      }

      perturbPlaneCoefficientsFixedPoints(f, range, fixing_vertices[f]);

      CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " is now fixed at " << f->get_plane() << " [measure=" << Size_shenanigans::length(f->get_plane()) << "]");

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_" + std::to_string(fixing_vertices[f].size()) + ".OFF", f);
#endif

      CGAL_SS3_TRANSF_TRACE_V(64, "Newly fixed facet F" << f->id() << " determines its high-degree incident vertices...");

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

      perturbPlaneCoefficientsFixedPoints(f, range, fixing_vertices[f]);

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

      if (fixing_vertices[f].size() == 3) {
        f->init_plane();
        Transformation::normalize_facet_plane(f);
      } else {
        perturbPlaneCoefficientsFixedPoints(f, range, fixing_vertices[f]);
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
                                      std::vector<FacetSPtr>& polygon_to_facet)
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
          std::abort();
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
          std::abort();
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
            polygon_to_facet.push_back(facet);
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
          std::abort();
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
          std::abort();
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
    namespace pred = PMP::Corefinement;

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

    CGAL::unordered_flat_map<VertexSPtr, Point_3> original_points;
    for (const VertexSPtr& v : polyhedron->vertices()) {
      original_points[v] = v->point();
    }

    CGAL::unordered_flat_map<FacetSPtr, Plane_3> original_planes;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      original_planes[facet] = facet->get_plane();
    }

    FT global_sq_max = CGAL::square(1e-7); // @tmp hardcoded
    CGAL::unordered_flat_map<VertexSPtr, FT> sq_max_displacements;

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

    auto is_stable = [&](const VertexSPtr& v,
                         const FT& max_sq_displacement) -> bool
    {
      for (auto it_wf1 = v->facets().begin(); it_wf1 != v->facets().end(); ++it_wf1) {
        if (FacetSPtr f1 = it_wf1->lock()) {
          for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
            if (FacetSPtr f2 = it_wf2->lock()) {
              for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                if (FacetSPtr f3 = it_wf3->lock()) {
                  std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                  if (!p_new.has_value()) {
                    CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
                    CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                    continue;
                  }

                  CGAL_assertion(f1->get_plane().has_on(*p_new));
                  CGAL_assertion(f2->get_plane().has_on(*p_new));
                  CGAL_assertion(f3->get_plane().has_on(*p_new));

                  const FT sqd = CGAL::squared_distance(original_points[v], p_new.value());

                  if (sqd > max_sq_displacement) {
                    CGAL_SS3_TRANSF_TRACE_V(32, "  TOO FAR: " << v->to_string());
                    CGAL_SS3_TRANSF_TRACE_V(32, "  from " << original_points[v] << " to " << p_new.value());
                    CGAL_SS3_TRANSF_TRACE_V(32, "  sq dist " << sqd << " VS " << max_sq_displacement);
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f1->id() << " [" << f1->get_plane().a() << " " << f1->get_plane().b() << " "
                                                                          << f1->get_plane().c() << " " << f1->get_plane().d() << "] "
                                                                          << CGAL::squared_distance(f1->get_plane(), original_points[v]));
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f2->id() << " [" << f2->get_plane().a() << " " << f2->get_plane().b() << " "
                                                                          << f2->get_plane().c() << " " << f2->get_plane().d() << "] "
                                                                          << CGAL::squared_distance(f2->get_plane(), original_points[v]));
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f3->id() << " [" << f3->get_plane().a() << " " << f3->get_plane().b() << " "
                                                                          << f3->get_plane().c() << " " << f3->get_plane().d() << "] "
                                                                          << CGAL::squared_distance(f3->get_plane(), original_points[v]));
                    return false;
                  }
                }
              }
            }
          }
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is stable");
      return true;
    };

    CGAL_SS3_TRANSF_TRACE_V(8, "BEFORE");

    // @todo if we have general position + stability + no high degree, we could skip perturbation
    CGAL_SS3_TRANSF_TRACE_CODE(std::vector<Stability_failure> failures;)
    CGAL_SS3_TRANSF_TRACE_V(32, "Unperturbed is in general position?\n" << are_planes_in_general_position(polyhedron, &failures));

#ifdef CGAL_SS3_DUMP_FILES
    std::ofstream out_unstable_before("results/unstable_vertices_before.xyz");
    out_unstable_before.precision(17);
#endif

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE_CODE(if(!is_stable(v, sq_max_displacements[v])) {);
#ifdef CGAL_SS3_DUMP_FILES
    CGAL_SS3_TRANSF_TRACE_CODE(out_unstable_before << v->point() << "\n");
#endif
    CGAL_SS3_TRANSF_TRACE_CODE(})

#ifdef CGAL_SS3_DUMP_FILES
    out_unstable_before.close();
#endif

    CGAL_SS3_TRANSF_TRACE_V(8, "START");

    {
      // V1: use unstable vertices as anchors
      // V2: fix vertices and facets (same as perturb_v3 but unstable vertices instead of high-degree vertices)
#define CGAL_SS3_PERTURB_V4_PART1_V1_QUATER

#ifdef CGAL_SS3_PERTURB_V4_PART1_V1
      // Part 1 V2
      // =========

      // - (1) Perturbations on all faces
      // - (2) Priority queue of unstable vertices by displacement (store vertex+facet+displacement)
      // - (3) Get most unstable anchor in incident faces of highest displacement (or anchor in face
      //       that appears in the largest number of high displacement vertices)
      // - (4) recompute the facet plane, taking into account the anchors
      // - (5) update the stability cost of the vertices of the just-anchored face

      CGAL::unordered_flat_map<FacetSPtr, boost::container::small_vector<VertexSPtr, 3> > anchors;

      struct Vertex_stability_info
      {
        bool is_stable = true;
        FT max_sq_displacement = FT(0);
        FacetSPtr worst_facet;
      };

      struct Vertex_stability_record
      {
        VertexSPtr vertex;
        FacetSPtr facet;
        FT max_sq_displacement;
        bool is_active = true;
        std::size_t vertex_id = -1;
      };

      struct Vertex_stability_record_less
      {
        explicit Vertex_stability_record_less(const std::vector<Vertex_stability_record>* records = nullptr)
          : records(records)
        {}

        bool operator()(std::size_t lhs, std::size_t rhs) const
        {
          CGAL_precondition(records != nullptr);
          return (*records)[lhs].max_sq_displacement > (*records)[rhs].max_sq_displacement;
        }

        const std::vector<Vertex_stability_record>* records;
      };

      auto evaluate_vertex_stability = [&](const VertexSPtr& v,
                                           const FT& sq_displacement_bound) -> Vertex_stability_info
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Check stability of V" << v->id() << " (deg: " << v->degree() << ")");

        Vertex_stability_info info;
        info.max_sq_displacement = sq_displacement_bound;

        for (auto it_wf1 = v->facets().begin(); it_wf1 != v->facets().end(); ++it_wf1) {
          if (FacetSPtr f1 = it_wf1->lock()) {
            for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
              if (FacetSPtr f2 = it_wf2->lock()) {
                for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                  if (FacetSPtr f3 = it_wf3->lock()) {
                    std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                    if (!p_new.has_value()) {
                      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
                      CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                      continue;
                    }

                    const FT sqd = CGAL::squared_distance(v->point(), p_new.value());
                    if (sqd > info.max_sq_displacement) {
                      auto dump_facet = [&](const FacetSPtr& f)
                      {
                        std::stringstream oss;
                        oss << "F" << f->id() << " " << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                     << f->get_plane().c() << " " << f->get_plane().d() << " ";
                        oss << "[dist = " << CGAL::approximate_sqrt(CGAL::squared_distance(f->get_plane(), v->point())) << "] ";
                        oss << "[#v = " << f->vertices().size() << "] ";
                        oss << "[anchors {";
                        const auto it_a = anchors.find(f);
                        if (it_a != anchors.end()) {
                          for (const VertexSPtr& av : it_a->second)
                            oss << " V" << av->id();
                        }
                        oss << " }]";
                        oss << " }]";
                        return oss.str();
                      };

                      CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is too far");
                      CGAL_SS3_TRANSF_TRACE_V(64, "  from " << v->point() << " to " << p_new.value());
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f1));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f2));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f3));

                      info.max_sq_displacement = sqd;

                      // pick the facet that has the least amount of anchors, hoping that anchoring that facet
                      // will solve the stability of possible other unstable vertices
                      info.worst_facet = (std::min)({f1, f2, f3}, [&](const FacetSPtr& lf, const FacetSPtr& rf) -> bool {
                          auto lit = anchors.find(lf);
                          auto rit = anchors.find(rf);

                          bool l_has_v = (lit != anchors.end()) &&
                                        (std::find(lit->second.begin(), lit->second.end(), v) != lit->second.end());
                          bool r_has_v = (rit != anchors.end()) &&
                                        (std::find(rit->second.begin(), rit->second.end(), v) != rit->second.end());

                          // prefer the facet that does not already have 'v' as anchor
                          if (l_has_v != r_has_v)
                              return !l_has_v;

#if 0
                          // otherwise, pick the one with the least amount of anchors
                          const std::size_t l_count = (lit != anchors.end()) ? lit->second.size() : 0;
                          const std::size_t r_count = (rit != anchors.end()) ? rit->second.size() : 0;
                          return l_count < r_count;
#else
                          // prioritize small faces
                          const std::size_t ln = lf->vertices().size();
                          const std::size_t rn = rf->vertices().size();
                          return ln < rn;
#endif
                      });

                      CGAL_SS3_TRANSF_TRACE_V(64, "  worst is F" << info.worst_facet->id());
                    }
                  }
                }
              }
            }
          }
        }

        info.is_stable = (info.max_sq_displacement <= sq_displacement_bound);

        if (info.is_stable) {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is stable (max displacement = " << CGAL::approximate_sqrt(info.max_sq_displacement) << ")");
        } else {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is unstable and would move at max "
              << CGAL::approximate_sqrt(info.max_sq_displacement) << " (tolerance "
              << CGAL::approximate_sqrt(sq_displacement_bound) << ")");
        }

        return info;
      };

      auto recompute_facet_plane = [&](const FacetSPtr& f)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Recompute F" << f->id());
        CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");

        static const boost::container::small_vector<VertexSPtr, 3> s_no_anchors;
        const auto it_a = anchors.find(f);
        const auto& anchor_vertices = (it_a != anchors.end()) ? it_a->second : s_no_anchors;

        for (const VertexSPtr& v : anchor_vertices) {
          CGAL_SS3_TRANSF_TRACE_V(32, "    V" << v->id() << " is an anchor of F" << f->id() << " at " << v->point());
        }

        if (anchor_vertices.size() < 3) {
          perturbPlaneCoefficientsFixedPoints(f, nudge_range, anchor_vertices);
        } else { // anchor_vertices.size() >= 3
          const Point_3& p0 = anchor_vertices[0]->point();
          const Point_3& p1 = anchor_vertices[1]->point();
          const Point_3& p2 = anchor_vertices[2]->point();

          // @fixme should a large deviation of normal be understood as an unstable facet
          // that ought to be triangulated? If selected anchors were aligned, a small
          // nudge is a large normal change.
          // Maybe it's not needed because if the normal changed a lot, then it's likely
          // other vertices won't be stable and the facet will be triangulated.
          Plane_3 new_pl(p0, p1, p2);

          // We don't know about the order of the fixed points, so trust the input normal for orientation
          // @fixme robustness with thin triangles... the nudge is small so it should be fine...
          if (new_pl.orthogonal_vector() * original_planes.at(f).orthogonal_vector() < 0) {
            new_pl = new_pl.opposite();
          }

          f->set_plane(new_pl);
          Transformation::normalize_facet_plane(f); // doesn't pass through the points anymore...

          // this should fail due to inexact normalization?...
          CGAL_postcondition(f->get_plane().has_on(p0));
          CGAL_postcondition(f->get_plane().has_on(p1));
          CGAL_postcondition(f->get_plane().has_on(p2));

          CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "] [3 anchors]");
        }
      };

      // (1)
      CGAL_SS3_TRANSF_TRACE_V(16, "Initial plane perturbations");
      for (const FacetSPtr& f : polyhedron->facets()) {
        recompute_facet_plane(f);
      }

      // (2)
      std::vector<Vertex_stability_record> unstable_vertices;
      CGAL::unordered_flat_map<VertexSPtr, std::size_t> unstable_record_indices;

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream out_unstable_base("results/base_unstable_vertices.xyz");
      out_unstable_base.precision(17);
#endif

      // initial fill
      CGAL_SS3_TRANSF_TRACE_V(16, "Initial stable/unstable classification");
      for (const VertexSPtr& v : polyhedron->vertices()) {
        const Vertex_stability_info info = evaluate_vertex_stability(v, sq_max_displacements[v]);
        if (!info.is_stable) {
#ifdef CGAL_SS3_DUMP_FILES
          out_unstable_base << v->point() << "\n";
#endif
          unstable_record_indices[v] = unstable_vertices.size();
          unstable_vertices.push_back({v, info.worst_facet, info.max_sq_displacement, true, 0});
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      out_unstable_base.close();
#endif

      using Stability_queue = Modifiable_priority_queue<std::size_t,
                                                        Vertex_stability_record_less,
                                                        boost::identity_property_map,
                                                        CGAL_BOOST_PAIRING_HEAP>;
      Stability_queue pq(polyhedron->vertices().size(), Vertex_stability_record_less(&unstable_vertices));
      for (std::size_t i = 0; i < unstable_vertices.size(); ++i) {
        pq.push(i);
      }

      auto update_vertex_record = [&](const VertexSPtr& v)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "update record of V" << v->id());
        const Vertex_stability_info info = evaluate_vertex_stability(v, sq_max_displacements[v]);

        const auto it = unstable_record_indices.find(v);
        if (it != unstable_record_indices.end()) { // already exists
          CGAL_SS3_TRANSF_TRACE_V(64, "record exists");
          Vertex_stability_record& record = unstable_vertices[it->second];
          record.vertex = v;
          record.vertex_id = v->id();
          record.facet = info.worst_facet;
          record.max_sq_displacement = info.max_sq_displacement;
          record.is_active = !info.is_stable;
          CGAL_SS3_TRANSF_TRACE_V(64, "activity: " << record.is_active);
          if (!record.is_active) {
            CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is now stable");
            record.facet = nullptr;
            record.max_sq_displacement = FT(0);
            // it could be purged from the PQ here, but it is simpler to perform it lazily, at pop time
          } else {
            if (pq.contains(it->second)) {
              CGAL_SS3_TRANSF_TRACE_V(64, "queue contains");
              pq.update(it->second);
            } else {
              pq.push(it->second);
            }
          }
        } else {
          if (!info.is_stable) {
            const std::size_t index = unstable_vertices.size();
            unstable_vertices.push_back({v, info.worst_facet, info.max_sq_displacement, true, 0});
            unstable_record_indices[v] = index;
            CGAL_SS3_TRANSF_TRACE_V(64, "new unstable vertex");
            pq.push(index);
          }
        }
      };

      std::unordered_set<VertexSPtr> all_anchors;

      // main loop
      while (!pq.empty()) {
        // (3) get the most unstable vertex and anchor it on the facet that caused the worst displacement
        CGAL_SS3_TRANSF_TRACE_V(16, "stable/unstable: main loop (" << pq.size() << ")");

        const std::size_t unstable_index = pq.top_and_pop();
        Vertex_stability_record record = unstable_vertices[unstable_index];
        if (!record.is_active) {
          continue;
        }

        CGAL_SS3_TRANSF_TRACE_V(16, "pop V" << record.vertex->id());
        CGAL_SS3_TRANSF_TRACE_V(16, "pop F" << record.facet->id());

        const FacetSPtr& f = record.facet;
        auto& facet_anchors = anchors[f];
        CGAL_assertion(std::find(facet_anchors.begin(), facet_anchors.end(), record.vertex) == facet_anchors.end());
        facet_anchors.push_back(record.vertex);

        // if the vertex is not already an anchor in any incident face, we fix its position
        auto make_random_anchor = [&](const VertexSPtr& v) {
          static std::random_device rd;
          unsigned int s = 0; // rd()
          // CGAL_SS3_TRANSF_TRACE("seed = " << s);
          static std::mt19937 gen(s);
          static std::uniform_real_distribution<> rdist(-nudge_range, nudge_range);

          double eps_x = rdist(gen);
          double eps_y = rdist(gen);
          double eps_z = rdist(gen);

#if 1
          const Point_3 new_p(CGAL::to_double(v->point().x() + eps_x),
                              CGAL::to_double(v->point().y() + eps_y),
                              CGAL::to_double(v->point().z() + eps_z));
#else
          Point_3 new_p(v->point().x() + eps_x,
                        v->point().y() + eps_y,
                        v->point().z() + eps_z);
#endif

          CGAL_SS3_TRANSF_TRACE_V(8, "New anchor V" << v->id() << " at " << new_p);
          CGAL_SS3_TRANSF_TRACE_V(8, "  Distance from original: " << CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), new_p)));

          v->set_point(new_p);
        };

        const VertexSPtr& v = record.vertex;
        auto ires = all_anchors.insert(v);
        if (ires.second) {
          make_random_anchor(v);
        }

        std::vector<FacetSPtr> facets_to_recompute;

        if (facet_anchors.size() > 3) {
          CGAL_SS3_TRANSF_TRACE_V(8, "  Must triangulate F" << f->id());
          auto [_, new_facets] = Transformation::triangulate_facet(f, polyhedron);
          const Plane_3 original_plane = original_planes.at(f);
          for (const FacetSPtr& nf : new_facets) {
            original_planes[nf] = original_plane;
            facets_to_recompute.push_back(nf);
          }
          original_planes.erase(f);

          // @todo we could backport the anchors to the new faces
          anchors.erase(f);
        } else {
          facets_to_recompute.push_back(f);
        }

        for (const FacetSPtr& nf : facets_to_recompute) {
          // (4) recompute the facet equation
          recompute_facet_plane(nf);

          // (5) update the facet's vertices stability in the queue
          for (const VertexSPtr& v : nf->vertices()) {
            update_vertex_record(v);
          }
        }
      }
#endif // CGAL_SS3_PERTURB_V4_PART1_V1

#ifdef CGAL_SS3_PERTURB_V4_PART1_V1_BIS
      // Some preprocessing: identify edges of the polyhedron whose vertices have more than
      //  three (or more) common incident facets, i.e. the set intersection of the set
      // of facets incident to each edge extremity has size > 3.
      // In such a configuration, we want to identify the facet(s) which are not incident
      // to the edge, i.e. the facet touches the source 'sv' and target 'tv' vertices, but there is
      // no edge 'sv-tv'. We split that facet such that 'sv' and 'tv' are on different facets
      // after the split.
      // This is to be performed otherwise we cannot guarantee general position if these vertices
      // were to become anchors.
      std::list<FacetSPtr> facets_to_preprocess;
      for (const FacetSPtr& f : polyhedron->facets()) {
        facets_to_preprocess.push_back(f);
      }

      for (;;) {
        bool did_something = false;
        std::list<FacetSPtr> facets_to_exclude;

        for (const FacetSPtr& f : facets_to_preprocess) {
          for (const EdgeSPtr& e : f->edges()) {
            const VertexSPtr& sv = e->src(f);
            const VertexSPtr& tv = e->tgt(f);

            // If we have not simplified exact coplanarity, we could create such a configuration
            // with degree 3 vertices...
            if (sv->degree() <= 3 || tv->degree() <= 3) {
              continue;
            }

            std::set<FacetSPtr> sv_facets;
            std::set<FacetSPtr> tv_facets;
            std::set<FacetSPtr> common_facets;

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

            std::set_intersection(sv_facets.begin(), sv_facets.end(),
                                  tv_facets.begin(), tv_facets.end(),
                                  std::inserter(common_facets, common_facets.begin()));

            for (const FacetSPtr& fprime : common_facets) {
              if (fprime == f) {
                continue;
              }

              if (tv->next(fprime) == sv || sv->next(fprime) == tv) {
                continue;
              }

              facets_to_exclude.push_back(fprime);
              CGAL_SS3_TRANSF_TRACE_V(32, "Facet F" << fprime->id() << " needs triangulating due to missing edge between V" << sv->id() << " and V" << tv->id());

              // @todo Transformation::split() with minimal split to create only two facets?
              const Plane_3 original_plane = original_planes.at(f);
              auto [_, new_facets] = Transformation::triangulate_facet(fprime, polyhedron);
              for (const FacetSPtr& nf : new_facets) {
                original_planes[nf] = original_plane;
              }
              original_planes.erase(f);

              did_something = true;
              break;
            }

            if (did_something) {
              break;
            }
          }

          if (did_something) {
            break;
          }
        }

        for (const FacetSPtr& f : facets_to_exclude) {
          facets_to_preprocess.remove(f);
        }

        if (!did_something) {
          break;
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      IO::write_OBJ("results/v4_preprocessed.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

      // Nudge all facets randomly
      for (const FacetSPtr& f : polyhedron->facets()) {
        perturbPlaneCoefficientsFixedPoints(f, nudge_range, std::vector<VertexSPtr>());
      }

      // Compute the initial set of unstable vertices
      CGAL::unordered_flat_set<VertexSPtr> unstable_vertices;
      for (const VertexSPtr& v : polyhedron->vertices()) {
        if (!is_stable(v, sq_max_displacements[v])) {
          unstable_vertices.insert(v);
        }
      }
      CGAL_SS3_TRANSF_TRACE_V(8, "initial unstable vertices = " << unstable_vertices.size());

      auto unstable_vertex_count = [&](const FacetSPtr& f) -> std::size_t
      {
        std::size_t count = 0;
        for (const VertexSPtr& v : f->vertices()) {
          if (unstable_vertices.count(v) != 0) {
            ++count;
          }
        }
        return count;
      };

      auto has_unstable_vertices = [&](const FacetSPtr& f) -> bool
      {
        return unstable_vertex_count(f) != 0;
      };

      CGAL::unordered_flat_set<VertexSPtr> fixed_vertices;
      CGAL::unordered_flat_set<FacetSPtr> fixed_facets;

      struct Facet_stability_record
      {
        FacetSPtr facet;
        std::size_t unstable_count = 0;
        bool is_active = true;
      };

      struct Facet_stability_record_less
      {
        explicit Facet_stability_record_less(const std::vector<Facet_stability_record>* records = nullptr)
          : records(records)
        {}

        bool operator()(std::size_t lhs, std::size_t rhs) const
        {
          CGAL_precondition(records != nullptr);
          if ((*records)[lhs].unstable_count != (*records)[rhs].unstable_count) {
            return (*records)[lhs].unstable_count > (*records)[rhs].unstable_count;
          }
          return (*records)[lhs].facet->vertices().size() > (*records)[rhs].facet->vertices().size();
        }

        const std::vector<Facet_stability_record>* records;
      };

      std::vector<Facet_stability_record> facet_records;
      CGAL::unordered_flat_map<FacetSPtr, std::size_t> facet_record_indices;
      using Facet_queue = Modifiable_priority_queue<std::size_t,
                                                    Facet_stability_record_less,
                                                    boost::identity_property_map,
                                                    CGAL_BOOST_PAIRING_HEAP>;
      Facet_queue active_queue(polyhedron->facets().size(), Facet_stability_record_less(&facet_records));

      // Records & queue for facets deferred for phase 2
      std::vector<Facet_stability_record> deferred_records;
      CGAL::unordered_flat_map<FacetSPtr, std::size_t> deferred_record_indices;
      // @fixme this is terrible and does not guarantee anything
      const std::size_t facet_queue_capacity = std::max<std::size_t>(polyhedron->facets().size() * 4U, std::size_t(16));
      Facet_queue deferred_queue(facet_queue_capacity, Facet_stability_record_less(&deferred_records));

      auto is_facet_fixed = [&](const FacetSPtr& f) -> bool
      {
        return fixed_facets.find(f) != fixed_facets.end();
      };

      auto count_fixed_incident_facets = [&](const VertexSPtr& v) -> std::size_t
      {
        std::size_t count = 0;
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (is_facet_fixed(f)) {
              ++count;
            }
          }
        }
        return count;
      };

      auto update_facet_record = [&](const FacetSPtr& f,
                                     std::vector<Facet_stability_record>& records,
                                     CGAL::unordered_flat_map<FacetSPtr, std::size_t>& record_indices,
                                     Modifiable_priority_queue<std::size_t,
                                                               Facet_stability_record_less,
                                                               boost::identity_property_map,
                                                               CGAL_BOOST_PAIRING_HEAP>& queue)
      {
        const std::size_t unstable_count = unstable_vertex_count(f);
        const auto it = record_indices.find(f);
        if (it != record_indices.end()) {
          Facet_stability_record& rec = records[it->second];
          rec.unstable_count = unstable_count;
          rec.is_active = (unstable_count != 0);
          if (rec.is_active) {
            if (queue.contains(it->second)) {
              queue.update(it->second);
            } else {
              queue.push(it->second);
            }
          }
        } else if (unstable_count != 0) {
          const std::size_t index = records.size();
          records.push_back({f, unstable_count, true});
          record_indices[f] = index;
          queue.push(index);
        }
      };

      auto update_vertex_stability = [&](const VertexSPtr& v,
                                        const bool defer_to_deferred_queue)
      {
        const bool was_unstable = (unstable_vertices.count(v) != 0);
        const bool is_unstable = !is_stable(v, sq_max_displacements[v]);
        if (is_unstable) {
          unstable_vertices.insert(v);
        } else {
          unstable_vertices.erase(v);
        }

        if (was_unstable != is_unstable) {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " stability changed -> " << (is_unstable ? "unstable" : "stable"));
        }

        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (is_facet_fixed(f)) {
              continue;
            }

            if (defer_to_deferred_queue) {
              update_facet_record(f, deferred_records, deferred_record_indices, deferred_queue);
            } else {
              update_facet_record(f, facet_records, facet_record_indices, active_queue);
            }
          }
        }
      };

      auto fix_vertex = [&](const VertexSPtr& v)
      {
        std::vector<const Plane_3*> constraining_planes;
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (is_facet_fixed(f)) {
              constraining_planes.push_back(&(f->get_plane()));
            }
          }
        }

        CGAL_assertion(constraining_planes.size() == 3);

        std::optional<Point_3> point = Kernel_wrapper::intersection(*(constraining_planes[0]),
                                                                    *(constraining_planes[1]),
                                                                    *(constraining_planes[2]));
        if (!point) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Error: triplet of planes does not define a point!");
          std::abort();
        }

        v->set_point(*point);
        fixed_vertices.insert(v);
      };

      auto is_vertex_partially_stable = [&](const VertexSPtr& v,
                                            const FacetSPtr& candidate_facet) -> bool
      {
        std::vector<FacetSPtr> previous_fixed_incident;
        previous_fixed_incident.reserve(3);
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (f != candidate_facet && is_facet_fixed(f)) {
              previous_fixed_incident.push_back(f);
            }
          }
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " partial stability with " << previous_fixed_incident.size() << " previous fixed facet(s) (" << v->degree() << " facets)");
        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream fixed_facets_ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : previous_fixed_incident) { fixed_facets_ss << " F" << f->id(); })
        CGAL_SS3_TRANSF_TRACE_V(32, "  fixed facets:" << fixed_facets_ss.str());

        if (previous_fixed_incident.empty()) {
          const FT sqd = CGAL::squared_distance(original_points[v], candidate_facet->get_plane());
          const bool stable = (sqd <= sq_max_displacements[v]);
          CGAL_SS3_TRANSF_TRACE_V(32, "Stability [0]: " << stable << " (sq distance = " << sqd << ", bound = " << sq_max_displacements[v] << ")");
          return stable;
        }

        if (previous_fixed_incident.size() == 1) {
          const FacetSPtr& g = previous_fixed_incident.front();
          std::optional<Line_3> line = Kernel_wrapper::intersection(candidate_facet->get_plane(), g->get_plane());
          if (!line) {
            CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " not partially stable [1]: candidate facet and F" << g->id() << " are parallel");
            return false;
          }
          const FT sqd = CGAL::squared_distance(original_points[v], line->projection(original_points[v]));
          const bool stable = (sqd <= sq_max_displacements[v]);
          CGAL_SS3_TRANSF_TRACE_V(32, "Stability [1]: " << stable << " (sq distance = " << sqd << ", bound = " << sq_max_displacements[v] << ")");
          return stable;
        }

        if (previous_fixed_incident.size() == 2) {
          const FacetSPtr& g0 = previous_fixed_incident[0];
          const FacetSPtr& g1 = previous_fixed_incident[1];
          std::optional<Point_3> p = Kernel_wrapper::intersection(candidate_facet->get_plane(), g0->get_plane(), g1->get_plane());
          if (!p) {
            CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " not partially stable [2]: candidate + F" << g0->id() << " + F" << g1->id() << " do not intersect");
            return false;
          }
          const FT sqd = CGAL::squared_distance(original_points[v], p.value());
          const bool stable = (sqd <= sq_max_displacements[v]);
          CGAL_SS3_TRANSF_TRACE_V(32, "Stability [2]: " << stable << " w/ F" << candidate_facet->id() << "(candidate) + F" << g0->id() << " + F" << g1->id() << " (sq distance = " << sqd << ", bound = " << sq_max_displacements[v] << ")");
          return stable;
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " has more than 2 fixed neighbors; only candidate-pair intersections are asserted");

        CGAL_assertion(previous_fixed_incident.size() > 2);
        for (std::size_t i = 0; i < previous_fixed_incident.size(); ++i) {
          for (std::size_t j = i + 1; j < previous_fixed_incident.size(); ++j) {
            std::optional<Point_3> p = Kernel_wrapper::intersection(candidate_facet->get_plane(),
                                                                    previous_fixed_incident[i]->get_plane(),
                                                                    previous_fixed_incident[j]->get_plane());
            if (!p) { //@fixme should be an assert? (others as well)
              CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " not partially stable·[3]: candidate + F" << previous_fixed_incident[i]->id() << " + F" << previous_fixed_incident[j]->id() << " do not intersect");
              return false;
            }

            const FT sqd = CGAL::squared_distance(original_points[v], p.value());
            const bool stable = (sqd <= sq_max_displacements[v]);
            CGAL_SS3_TRANSF_TRACE_V(32, "Stability [3]: " << stable << " w/ F" << candidate_facet->id() << "(candidate) + F" << previous_fixed_incident[i]->id() << " + F" << previous_fixed_incident[j]->id() << " (sq distance = " << sqd << ", bound = " << sq_max_displacements[v] << ")");
            if (!stable) {
              return stable;
            }
          }
        }

        return true;
      };

      // check partial stability for all triangle vertices
      auto all_vertices_partially_stable = [&](const FacetSPtr& f) -> bool {
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_vertex_partially_stable(v, f)) {
            return false;
          }
        }
        return true;
      };

      auto recompute_facet_plane = [&](const FacetSPtr& f)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Recompute F" << f->id());
        CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");

        std::vector<VertexSPtr> anchor_vertices;
        for (const VertexSPtr& v : f->vertices()) {
          if (fixed_vertices.count(v)) {
            anchor_vertices.push_back(v);
          }
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "try to fix F" << f->id() << " with " << anchor_vertices.size() << " constraints...");
        for (const VertexSPtr& v : anchor_vertices) {
          CGAL_SS3_TRANSF_TRACE_V(32, "    V" << v->id() << " is an anchor of F" << f->id() << " at " << v->point());
        }

        // If there is no anchor, there is no point adjusting the perturbation yet,
        // since an initial one has already been performed before the first round
        // of stability evaluations.
        if (anchor_vertices.size() > 0 && anchor_vertices.size() < 3) {
          perturbPlaneCoefficientsFixedPoints(f, nudge_range, anchor_vertices);
        } else if (anchor_vertices.size() >= 3) {
          const Point_3& p0 = anchor_vertices[0]->point();
          const Point_3& p1 = anchor_vertices[1]->point();
          const Point_3& p2 = anchor_vertices[2]->point();

          // @fixme should a large deviation of normal be understood as an unstable facet
          // that ought to be triangulated? If selected anchors were aligned, a small
          // nudge is a large normal change.
          // Maybe it's not needed because if the normal changed a lot, then it's likely
          // other vertices won't be stable and the facet will be triangulated.
          Plane_3 new_pl(p0, p1, p2);

          // We don't know about the order of the fixed points, so trust the input normal for orientation
          // @fixme robustness with thin triangles... the nudge is small so it should be fine...
          if (new_pl.orthogonal_vector() * original_planes.at(f).orthogonal_vector() < 0) {
            new_pl = new_pl.opposite();
          }

          f->set_plane(new_pl);
          Transformation::normalize_facet_plane(f); // doesn't pass through the points anymore...

          // this should fail due to inexact normalization?...
          CGAL_postcondition(f->get_plane().has_on(p0));
          CGAL_postcondition(f->get_plane().has_on(p1));
          CGAL_postcondition(f->get_plane().has_on(p2));

          CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "] [3 anchors]");
        }
      };

      auto try_to_fix_facet = [&](const FacetSPtr& f) -> bool
      {
        CGAL_precondition(!is_facet_fixed(f));

        // This isn't good because we constrain a lot here for stability and triangulating
        // creates more unstable points and cascading. Thus, we get more cascading, not less
        // when we triangulate (contrary to the V3 perturbation).
        //
        // const std::size_t max_length = 100;
        // std::size_t l = Size_shenanigans::length(f->get_plane());
        // if (l > max_length) {
        //   CGAL_SS3_TRANSF_TRACE_V(32, "Facet F" << f->id() << " would be too long");
        //   CGAL_SS3_TRANSF_TRACE_V(32, CGAL::exact(f->get_plane()) << " (l=" << l << ")");
        //   CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " deferred: too large representation");
        //   return false;
        // }

        const Plane_3 old_plane = f->get_plane(); // intentional copy

        recompute_facet_plane(f);

        bool good_nudge = true;
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_vertex_partially_stable(v, f)) {
            good_nudge = false;
            f->set_plane(old_plane);
            break;
          }
        }

        if (!good_nudge) {
          CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " deferred: not partially stable");
          return false;
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "fix F" << f->id());

        // Found a good nudge with all vertices partially stable, fix the facet
        fixed_facets.insert(f);
        return true;
      };

      // we only put facets with unstable vertices in queue because we start by
      // nudging all facets, so general position is already acquired and we don't need
      // to come back to stable facets.
      for (const FacetSPtr& f : polyhedron->facets()) {
        if (!has_unstable_vertices(f)) {
          continue;
        }
        if (f->is_triangle()) {
          update_facet_record(f, deferred_records, deferred_record_indices, deferred_queue);
        } else {
          update_facet_record(f, facet_records, facet_record_indices, active_queue);
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(32, "Initial active queue: " << active_queue.size());
      CGAL_SS3_TRANSF_TRACE_V(32, "Initial deferred queue: " << deferred_queue.size());

      while (!active_queue.empty()) {
        const std::size_t unstable_index = active_queue.top_and_pop();
        Facet_stability_record record = facet_records[unstable_index];
        CGAL_SS3_TRANSF_TRACE_V(16, "pop active " << record.facet->to_string() << " (unstable_count = " << record.unstable_count << ")");
        if (!record.is_active) {
          continue;
        }

        CGAL_assertion(record.facet != nullptr);
        CGAL_assertion(!record.facet->is_triangle());

        if (!try_to_fix_facet(record.facet)) {
          CGAL_SS3_TRANSF_TRACE_V(32, "defer F" << record.facet->id());
          update_facet_record(record.facet, deferred_records, deferred_record_indices, deferred_queue);
          continue;
        }

        for (const VertexSPtr& v : record.facet->vertices()) {
          std::size_t fixed_inc = count_fixed_incident_facets(v);
          if (fixed_inc == 3) {
            fix_vertex(v);
          }
          update_vertex_stability(v, true);
        }
      }

      // Phase 2: process deferred facets. For polygonal facets we first try again to fix them;
      // only if they remain unfixable do we then triangulate and treat the resulting triangles.
      while (!deferred_queue.empty()) {
        const std::size_t def_index = deferred_queue.top_and_pop();
        Facet_stability_record def_record = deferred_records[def_index];
        if (!def_record.is_active || def_record.facet == nullptr || is_facet_fixed(def_record.facet)) {
          continue;
        }

        FacetSPtr f = def_record.facet;
        CGAL_SS3_TRANSF_TRACE_V(32, "pop deferred F" << f->id());

        if (!f->is_triangle()) {
          if (try_to_fix_facet(f)) {
            for (const VertexSPtr& v : f->vertices()) {
              std::size_t fixed_inc = count_fixed_incident_facets(v);
              if (fixed_inc == 3) {
                fix_vertex(v);
              }
              update_vertex_stability(v, true);
            }
            continue;
          }

          CGAL_SS3_TRANSF_TRACE_V(32, "Deferred F" << f->id() << " still unfixable, triangulate it");
          const Plane_3 original_plane = original_planes.at(f);

          auto [local_vertices, new_facets] = Transformation::triangulate_facet(f, polyhedron);
          for (const FacetSPtr& nf : new_facets) {
            CGAL_SS3_TRANSF_TRACE_V(32, "new facet: F" << nf->id());
            original_planes[nf] = original_plane;
            recompute_facet_plane(nf);
            update_facet_record(nf, deferred_records, deferred_record_indices, deferred_queue);
          }
          original_planes.erase(f);

          for (const VertexSPtr& v : local_vertices) {
            update_vertex_stability(v, true);
          }

          continue;
        }

        // f is a triangle: proceed with triangle handling
        FacetSPtr tri_f = f;
        CGAL_assertion(tri_f != FacetSPtr() && tri_f->is_triangle());

        CGAL_SS3_TRANSF_TRACE_V(32, "Process triangle F" << tri_f->id());

        // initial constraints: incident vertices that are already fixed
        boost::container::small_vector<VertexSPtr, 3> constraints;
        for (const VertexSPtr& v : tri_f->vertices()) {
          if (fixed_vertices.count(v)) {
            constraints.push_back(v);
          }
        }

        // try initial nudge with existing fixed vertices
        recompute_facet_plane(tri_f);

        // otherwise iteratively create up to 3 new anchors (per-triangle) and retry
        if (all_vertices_partially_stable(tri_f)) {
          fixed_facets.insert(tri_f);
          // notify any queues about this change
          for (const VertexSPtr& v : tri_f->vertices()) {
            update_vertex_stability(v, true);
          }
          continue;
        }

        auto compute_anchor_pos = [&](const VertexSPtr& v) -> Point_3 {
          // build list of constraining planes from already-fixed incident facets
          std::vector<const Plane_3*> constraining_planes;
          for (FacetWPtr wf : v->facets()) {
            if (FacetSPtr ff = wf.lock()) {
              if (fixed_facets.count(ff)) {
                constraining_planes.push_back(&ff->get_plane());
              }
            }
          }

          CGAL_SS3_TRANSF_TRACE_V(32, "Anchor V" << v->id() << " with " << constraining_planes.size() << " incident fixed facets");

          const Point_3& p = original_points[v];
          std::array<double, 3> r = rand_vec(-nudge_range/2.0, nudge_range/2.0);
          const double px_t = CGAL::to_double(p.x()) + r[0];
          const double py_t = CGAL::to_double(p.y()) + r[1];
          const double pz_t = CGAL::to_double(p.z()) + r[2];

          if (constraining_planes.empty()) {
            return Point_3{ FT(px_t), FT(py_t), FT(pz_t) };
          } else if (constraining_planes.size() == 1) {
#if 0
            const Plane_3& plane = *(constraining_planes[0]);
            Point_3 pp = plane.projection(Point_3{ FT(px_t), FT(py_t), FT(pz_t) });
            CGAL_assertion(constraining_planes[0]->has_on(pp));
            return pp;
#else
            const Plane_3& plane = *(constraining_planes[0]);
            const double ad = CGAL::to_double(plane.a());
            const double bd = CGAL::to_double(plane.b());
            const double cd = CGAL::to_double(plane.c());
            const double dd = CGAL::to_double(plane.d());

            const double n2 = ad*ad + bd*bd + cd*cd;
            CGAL_assertion(n2 > 0.0);

            // project nudged target point onto plane in double precision
            const double dist = (ad * px_t + bd * py_t + cd * pz_t + dd) / n2;
            const double proj_x = px_t - dist * ad;
            const double proj_y = py_t - dist * bd;
            const double proj_z = pz_t - dist * cd;

            // pivot on largest normal component to guarantee |slope| <= 1.0
            const double abs_a = std::abs(ad);
            const double abs_b = std::abs(bd);
            const double abs_c = std::abs(cd);

            const FT& a = plane.a();
            const FT& b = plane.b();
            const FT& c = plane.c();
            const FT& d = plane.d();

            Point_3 pp;
            if (abs_c >= abs_a && abs_c >= abs_b) {
              // Free variables: x and y
              FT x(proj_x);
              FT y(proj_y);
              FT z = -(a * x + b * y + d) / c;
              pp = Point_3(x, y, z);
            } else if (abs_b >= abs_a && abs_b >= abs_c) {
              // Free variables: x and z
              FT x(proj_x);
              FT z(proj_z);
              FT y = -(a * x + c * z + d) / b;
              pp = Point_3(x, y, z);
            } else {
              // Free variables: y and z
              FT y(proj_y);
              FT z(proj_z);
              FT x = -(b * y + c * z + d) / a;
              pp = Point_3(x, y, z);
            }

            CGAL_assertion(constraining_planes[0]->has_on(pp));
            return pp;
#endif
          } else if (constraining_planes.size() == 2) {
#if 0
            const Plane_3& plane1 = *(constraining_planes[0]);
            const Plane_3& plane2 = *(constraining_planes[1]);
            std::optional<Line_3> line = Kernel_wrapper::intersection(plane1, plane2);
            // std::cout << "  Constraint Line" << std::endl;
            // std::cout << "    " << line->point(0) << std::endl;
            // std::cout << "    " << line->point(1) << std::endl;
            Point_3 pp = line->projection(Point_3{ FT(px_t), FT(py_t), FT(pz_t) });
            CGAL_assertion(constraining_planes[0]->has_on(pp));
            CGAL_assertion(constraining_planes[1]->has_on(pp));
            return pp;
#else
            const Plane_3& plane1 = *(constraining_planes[0]);
            const Plane_3& plane2 = *(constraining_planes[1]);

            const double a1d = CGAL::to_double(plane1.a());
            const double b1d = CGAL::to_double(plane1.b());
            const double c1d = CGAL::to_double(plane1.c());
            const double d1d = CGAL::to_double(plane1.d());

            const double a2d = CGAL::to_double(plane2.a());
            const double b2d = CGAL::to_double(plane2.b());
            const double c2d = CGAL::to_double(plane2.c());
            const double d2d = CGAL::to_double(plane2.d());

            // Line direction in double: u = n1 x n2
            const double uxd = b1d * c2d - c1d * b2d;
            const double uyd = c1d * a2d - a1d * c2d;
            const double uzd = a1d * b2d - b1d * a2d;
            const double detG = uxd * uxd + uyd * uyd + uzd * uzd;
            CGAL_assertion(detG > 0.0); // Planes are not parallel

            // Clean 2-plane orthogonal projection in double precision
            const double dot11 = a1d*a1d + b1d*b1d + c1d*c1d;
            const double dot22 = a2d*a2d + b2d*b2d + c2d*c2d;
            const double dot12 = a1d*a2d + b1d*b2d + c1d*c2d;

            const double res1 = a1d * px_t + b1d * py_t + c1d * pz_t + d1d;
            const double res2 = a2d * px_t + b2d * py_t + c2d * pz_t + d2d;

            const double alpha = (res1 * dot22 - res2 * dot12) / detG;
            const double beta  = (res2 * dot11 - res1 * dot12) / detG;

            const double proj_x = px_t - alpha * a1d - beta * a2d;
            const double proj_y = py_t - alpha * b1d - beta * b2d;
            const double proj_z = pz_t - alpha * c1d - beta * c2d;

            // Pivot on largest line direction component
            const double abs_ux = std::abs(uxd);
            const double abs_uy = std::abs(uyd);
            const double abs_uz = std::abs(uzd);

            const FT& a1 = plane1.a(); const FT& b1 = plane1.b(); const FT& c1 = plane1.c(); const FT& d1 = plane1.d();
            const FT& a2 = plane2.a(); const FT& b2 = plane2.b(); const FT& c2 = plane2.c(); const FT& d2 = plane2.d();

            const FT ux = b1 * c2 - c1 * b2;
            const FT uy = c1 * a2 - a1 * c2;
            const FT uz = a1 * b2 - b1 * a2;

            Point_3 pp;
            if (abs_uz >= abs_ux && abs_uz >= abs_uy) {
              // Free variable: z (exact 53-bit dyadic, 1 limb)
              FT z(proj_z);
              FT x = (ux * z + (b1 * d2 - d1 * b2)) / uz;
              FT y = (uy * z + (d1 * a2 - a1 * d2)) / uz;
              pp = Point_3(x, y, z);
            } else if (abs_uy >= abs_ux && abs_uy >= abs_uz) {
              // Free variable: y
              FT y(proj_y);
              FT x = (ux * y + (d1 * c2 - c1 * d2)) / uy;
              FT z = (uz * y + (a1 * d2 - d1 * a2)) / uy;
              pp = Point_3(x, y, z);
            } else {
              // Free variable: x
              FT x(proj_x);
              FT y = (uy * x + (c1 * d2 - d1 * c2)) / ux;
              FT z = (uz * x + (d1 * b2 - b1 * d2)) / ux;
              pp = Point_3(x, y, z);
            }

            CGAL_assertion(constraining_planes[0]->has_on(pp));
            CGAL_assertion(constraining_planes[1]->has_on(pp));
            return pp;
#endif
          } else {
            std::optional<Point_3> point = Kernel_wrapper::intersection(*constraining_planes[0], *constraining_planes[1], *constraining_planes[2]);
            CGAL_assertion(point.has_value());
            return point.value();
          }
        };

        // add anchors until the triangle becomes partially stable or until we have exhausted anchors
        // @todo avoid recomputing partial stability over and over
        while (constraints.size() < 3 && !all_vertices_partially_stable(tri_f)) {
          CGAL_SS3_TRANSF_TRACE_V(32, "Must stabilize F" << tri_f->id() << " [" << constraints.size() << " anchors]...");

          // find a vertex that is not partially stable and not yet used as a constraint
          VertexSPtr target_v;
          for (const VertexSPtr& v : tri_f->vertices()) {
            // ignore if already fixed
            if (std::find(constraints.begin(), constraints.end(), v) != constraints.end()) {
              continue;
            }
            if (!is_vertex_partially_stable(v, tri_f)) {
              target_v = v;
              break;
            }
          }

          CGAL_assertion (target_v != VertexSPtr());

          const Point_3 anchor = compute_anchor_pos(target_v);
          CGAL_SS3_TRANSF_TRACE_V(32, "Anchor V" << target_v->id() << " at " << anchor);

          // @fixme
          // 1. is that guaranteed?
          // 2. is stability as a whole guaranteed? if not, we could have to do something really nasty
          // like backtracking (or restarting) and fixing that vertex in all facets from the beginning...
          CGAL_assertion_code(const FT sqd = CGAL::squared_distance(anchor, original_points[target_v]));
          CGAL_assertion(sqd < sq_max_displacements[target_v]);

          // set the original vertex to the anchor position and mark it fixed so subsequent triangles use it
          target_v->set_point(anchor);
          fixed_vertices.insert(target_v);

          // recompute plane constrained to current anchors (temporary anchors are only for this facet)
          constraints.push_back(target_v);
          recompute_facet_plane(tri_f);

          // notify queues that this vertex changed stability
          update_vertex_stability(target_v, true);
        }

        if (all_vertices_partially_stable(tri_f)) {
          fixed_facets.insert(tri_f);
          for (const VertexSPtr& v : tri_f->vertices()) {
            update_vertex_stability(v, true);
          }
        } else {
          CGAL_SS3_TRANSF_TRACE_V(1, "Warning: Triangle F" << tri_f->id() << " remains partially unstable after anchoring");
          std::abort();
        }
      }
#endif // CGAL_SS3_PERTURB_V4_PART1_V1_BIS

#ifdef CGAL_SS3_PERTURB_V4_PART1_V1_QUATER
      // Some preprocessing: identify edges of the polyhedron whose vertices have more than
      //  three (or more) common incident facets, i.e. the set intersection of the set
      // of facets incident to each edge extremity has size > 3.
      // In such a configuration, we want to identify the facet(s) which are not incident
      // to the edge, i.e. the facet touches the source 'sv' and target 'tv' vertices, but there is
      // no edge 'sv-tv'. We split that facet such that 'sv' and 'tv' are on different facets
      // after the split.
      // This is to be performed otherwise we cannot guarantee general position if these vertices
      // were to become anchors.
      std::list<FacetSPtr> facets_to_preprocess;
      for (const FacetSPtr& f : polyhedron->facets()) {
        facets_to_preprocess.push_back(f);
      }

      for (;;) {
        bool did_something = false;
        std::list<FacetSPtr> facets_to_exclude;

        for (const FacetSPtr& f : facets_to_preprocess) {
          for (const EdgeSPtr& e : f->edges()) {
            const VertexSPtr& sv = e->src(f);
            const VertexSPtr& tv = e->tgt(f);

            // If we have not simplified exact coplanarity, we could create such a configuration
            // with degree 3 vertices...
            if (sv->degree() <= 3 || tv->degree() <= 3) {
              continue;
            }

            std::set<FacetSPtr> sv_facets;
            std::set<FacetSPtr> tv_facets;
            std::set<FacetSPtr> common_facets;

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

            std::set_intersection(sv_facets.begin(), sv_facets.end(),
                                  tv_facets.begin(), tv_facets.end(),
                                  std::inserter(common_facets, common_facets.begin()));

            for (const FacetSPtr& fprime : common_facets) {
              if (fprime == f) {
                continue;
              }

              if (tv->next(fprime) == sv || sv->next(fprime) == tv) {
                continue;
              }

              facets_to_exclude.push_back(fprime);
              CGAL_SS3_TRANSF_TRACE_V(32, "Facet F" << fprime->id() << " needs triangulating due to missing edge between V" << sv->id() << " and V" << tv->id());

              // @todo Transformation::split() with minimal split to create only two facets?
              const Plane_3 original_plane = original_planes.at(f);
              auto [_, new_facets] = Transformation::triangulate_facet(fprime, polyhedron);
              for (const FacetSPtr& nf : new_facets) {
                original_planes[nf] = original_plane;
              }
              original_planes.erase(f);

              did_something = true;
              break;
            }

            if (did_something) {
              break;
            }
          }

          if (did_something) {
            break;
          }
        }

        for (const FacetSPtr& f : facets_to_exclude) {
          facets_to_preprocess.remove(f);
        }

        if (!did_something) {
          break;
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      IO::write_OBJ("results/v4_preprocessed.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

      // - (0) Perturbations on all faces
      // - (A) Freeze large polygons at their perturbed plane, as long as every vertex keeps
      //       <= 3 frozen incident facets whose partial intersection is within kappa*delta
      // - (B) Priority queue of unstable vertices by displacement (vertex+facet+displacement)
      //   (3) pop the most unstable vertex; anchor it in the offending facet that is neither
      //       frozen nor already anchored at that vertex (triangle > polygon w/ <3 anchors > polygon w/ 3)
      //   (4) recompute the facet plane, taking into account the anchors
      //   (5) update the stability cost of the vertices of the just-anchored face
      //
      // Anchors are placed ON the frozen structure of the vertex (point / line / plane), so
      // frozen facets are always concurrent with the anchor and never block progress.

      const double kappa = 0.5; // fraction of delta that the frozen structure may consume

      CGAL::unordered_flat_map<FacetSPtr, boost::container::small_vector<VertexSPtr, 3> > anchors;
      CGAL::unordered_flat_set<FacetSPtr> frozen;
      CGAL::unordered_flat_map<VertexSPtr, boost::container::small_vector<FacetSPtr, 3> > frozen_at;
      std::unordered_set<VertexSPtr> all_anchors;

      // remaining tolerance around v->point() (shrinks once the vertex is anchored)
      CGAL::unordered_flat_map<VertexSPtr, FT> sq_budget;
      for (const VertexSPtr& v : polyhedron->vertices()) {
        sq_budget[v] = sq_max_displacements[v];
      }

      std::mt19937 gen(0);

      auto is_anchor_of = [&](const VertexSPtr& v, const FacetSPtr& f) -> bool
      {
        const auto it = anchors.find(f);
        return it != anchors.end() && std::find(it->second.begin(), it->second.end(), v) != it->second.end();
      };

      auto anchor_count = [&](const FacetSPtr& f) -> std::size_t
      {
        const auto it = anchors.find(f);
        return (it != anchors.end()) ? it->second.size() : 0;
      };

      struct Vertex_stability_info
      {
        bool is_stable = true;
        FT max_sq_displacement = FT(0);
        FacetSPtr worst_facet;
      };

      struct Vertex_stability_record
      {
        VertexSPtr vertex;
        FacetSPtr facet;
        FT max_sq_displacement;
        bool is_active = true;
        std::size_t vertex_id = -1;
      };

      struct Vertex_stability_record_less
      {
        explicit Vertex_stability_record_less(const std::vector<Vertex_stability_record>* records = nullptr)
          : records(records)
        {}

        bool operator()(std::size_t lhs, std::size_t rhs) const
        {
          CGAL_precondition(records != nullptr);
          const Vertex_stability_record& l = (*records)[lhs];
          const Vertex_stability_record& r = (*records)[rhs];
          if (l.max_sq_displacement != r.max_sq_displacement)
            return l.max_sq_displacement > r.max_sq_displacement;
          return l.vertex_id < r.vertex_id;
        }

        const std::vector<Vertex_stability_record>* records;
      };

      auto evaluate_vertex_stability = [&](const VertexSPtr& v,
                                           const FT& sq_displacement_bound) -> Vertex_stability_info
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Check stability of V" << v->id() << " (deg: " << v->degree() << ")");

        Vertex_stability_info info;
        info.max_sq_displacement = sq_displacement_bound;

        for (auto it_wf1 = v->facets().begin(); it_wf1 != v->facets().end(); ++it_wf1) {
          if (FacetSPtr f1 = it_wf1->lock()) {
            for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
              if (FacetSPtr f2 = it_wf2->lock()) {
                for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                  if (FacetSPtr f3 = it_wf3->lock()) {
                    std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                    if (!p_new.has_value()) {
                      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
                      CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                      continue;
                    }

                    const FT sqd = CGAL::squared_distance(v->point(), p_new.value());
                    if (sqd > info.max_sq_displacement) {
                      auto dump_facet = [&](const FacetSPtr& f)
                      {
                        std::stringstream oss;
                        oss << "F" << f->id() << " " << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                     << f->get_plane().c() << " " << f->get_plane().d() << " ";
                        oss << "[dist = " << CGAL::approximate_sqrt(CGAL::squared_distance(f->get_plane(), v->point())) << "] ";
                        oss << "[#v = " << f->vertices().size() << "] ";
                        oss << (frozen.count(f) ? "[frozen] " : "");
                        oss << "[anchors {";
                        for (const VertexSPtr& av : anchors[f])
                          oss << " V" << av->id();
                        oss << " }]";
                        return oss.str();
                      };

                      CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is too far");
                      CGAL_SS3_TRANSF_TRACE_V(64, "  from " << v->point() << " to " << p_new.value());
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f1));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f2));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f3));

                      info.max_sq_displacement = sqd;

                      // Candidates: facets of the triplet that are neither frozen nor already anchored at v.
                      // Progress argument: frozen facets and facets anchored at v all pass through the anchor,
                      // so a displaced triplet always contains at least one candidate.
                      boost::container::small_vector<FacetSPtr, 3> cands;
                      for (const FacetSPtr& g : {f1, f2, f3}) {
                        if (frozen.count(g) != 0) continue;
                        if (is_anchor_of(v, g)) continue;
                        cands.push_back(g);
                      }

                      if (cands.empty()) {
                        // Only possible before v is anchored (v->point() == original point) and with a
                        // fully frozen triplet, which Phase A guarantees to be within kappa*delta.
                        CGAL_SS3_TRANSF_TRACE_V(1, "Error: unstable triplet with no absorber at V" << v->id());
                        std::abort();
                      }

                      info.worst_facet = *std::min_element(cands.begin(), cands.end(),
                        [&](const FacetSPtr& lf, const FacetSPtr& rf) -> bool {
                          // free absorber first
                          if (lf->is_triangle() != rf->is_triangle()) return lf->is_triangle();
                          // avoid facets that would need a triangulation
                          const std::size_t la = anchor_count(lf);
                          const std::size_t ra = anchor_count(rf);
                          if ((la == 3) != (ra == 3)) return la < ra;
                          // prioritize small faces
                          const std::size_t sl = lf->vertices().size();
                          const std::size_t sr = rf->vertices().size();
                          if (sl != sr) return sl < sr;
                          // total order tie-breaker
                          return lf->id() < rf->id();
                        });

                      CGAL_SS3_TRANSF_TRACE_V(64, "  worst is F" << info.worst_facet->id());
                    }
                  }
                }
              }
            }
          }
        }

        info.is_stable = (info.max_sq_displacement <= sq_displacement_bound);

        if (info.is_stable) {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is stable (max displacement = " << CGAL::approximate_sqrt(info.max_sq_displacement) << ")");
        } else {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is unstable and would move at max "
              << CGAL::approximate_sqrt(info.max_sq_displacement) << " (tolerance "
              << CGAL::approximate_sqrt(sq_displacement_bound) << ")");
        }

        return info;
      };

      auto recompute_facet_plane = [&](const FacetSPtr& f)
      {
        CGAL_precondition(frozen.count(f) == 0);

        CGAL_SS3_TRANSF_TRACE_V(32, "Recompute F" << f->id());
        CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");

        const auto& anchor_vertices = anchors[f];
        for (const VertexSPtr& v : anchor_vertices) {
          CGAL_SS3_TRANSF_TRACE_V(32, "    V" << v->id() << " is an anchor of F" << f->id() << " at " << v->point());
        }

        if (anchor_vertices.size() < 3) {
          perturbPlaneCoefficientsFixedPoints(f, nudge_range, anchor_vertices);
        } else { // anchor_vertices.size() == 3
          const Point_3& p0 = anchor_vertices[0]->point();
          const Point_3& p1 = anchor_vertices[1]->point();
          const Point_3& p2 = anchor_vertices[2]->point();

          Plane_3 new_pl(p0, p1, p2);

          // We don't know about the order of the fixed points, so trust the input normal for orientation
          if (new_pl.orthogonal_vector() * original_planes.at(f).orthogonal_vector() < 0) {
            new_pl = new_pl.opposite();
          }

          f->set_plane(new_pl);
          Transformation::normalize_facet_plane(f);

          CGAL_postcondition(f->get_plane().has_on(p0));
          CGAL_postcondition(f->get_plane().has_on(p1));
          CGAL_postcondition(f->get_plane().has_on(p2));

          CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "] [3 anchors]");
        }
      };

      // (0) initial perturbations
      CGAL_SS3_TRANSF_TRACE_V(16, "Initial plane perturbations");
      for (const FacetSPtr& f : polyhedron->facets()) {
        recompute_facet_plane(f);
      }

      // -----------------------------------------------------------------------
      // (A) Freeze large polygons
      // -----------------------------------------------------------------------

      // Is the frozen structure at v (+ candidate) still within kappa*delta of the original point?
      auto frozen_structure_ok = [&](const VertexSPtr& v, const FacetSPtr& cand) -> bool
      {
        boost::container::small_vector<const Plane_3*, 4> pl;
        for (const FacetSPtr& g : frozen_at[v]) {
          pl.push_back(&(g->get_plane()));
        }
        if (cand) {
          pl.push_back(&(cand->get_plane()));
        }
        if (pl.size() > 3) {
          return false;
        }

        const Point_3& p = original_points[v];
        const FT bound = FT(kappa * kappa) * sq_max_displacements[v];

        if (pl.size() == 1) {
          return CGAL::squared_distance(p, *pl[0]) <= bound;
        }
        if (pl.size() == 2) {
          std::optional<Line_3> L = Kernel_wrapper::intersection(*pl[0], *pl[1]);
          return L && CGAL::squared_distance(p, L->projection(p)) <= bound;
        }
        if (pl.size() == 3) {
          std::optional<Point_3> q = Kernel_wrapper::intersection(*pl[0], *pl[1], *pl[2]);
          return q && CGAL::squared_distance(p, *q) <= bound;
        }
        return true;
      };

      {
        std::vector<FacetSPtr> order(polyhedron->facets().begin(), polyhedron->facets().end());
        std::sort(order.begin(), order.end(), [](const FacetSPtr& a, const FacetSPtr& b) {
          if (a->vertices().size() != b->vertices().size())
            return a->vertices().size() > b->vertices().size();
          return a->id() < b->id();
        });

        for (const FacetSPtr& f : order) {
          if (f->is_triangle()) {
            // triangles are the absorbers, never freeze them
            continue;
          }

          bool ok = true;
          for (const VertexSPtr& v : f->vertices()) {
            if (!frozen_structure_ok(v, f)) {
              ok = false;
              break;
            }
          }
          if (!ok) {
            CGAL_SS3_TRANSF_TRACE_V(32, "cannot freeze F" << f->id());
            continue;
          }

          frozen.insert(f);
          for (const VertexSPtr& v : f->vertices()) frozen_at[v].push_back(f);
          CGAL_SS3_TRANSF_TRACE_V(16, "freeze F" << f->id() << " (" << f->vertices().size() << " vertices)");
        }

        CGAL_SS3_TRANSF_TRACE_V(8, "frozen facets: " << frozen.size() << " / " << polyhedron->facets().size());
      }

      // -----------------------------------------------------------------------
      // (B) Vertex-driven anchoring on the frozen structure
      // -----------------------------------------------------------------------

      // Anchor placement: on the intersection of the frozen facets at v (point/line/plane),
      // random otherwise. Remaining budget is measured from the anchor.
#if 1
      // Anchor placement on the frozen structure at v:
      //   0 frozen: pure double point (all coordinates dyadic)
      //   1 frozen: 2 free double coordinates, 1 exact dependent (pivot = largest |normal| comp.)
      //   2 frozen: 1 free double coordinate, 2 exact dependent (pivot = largest |direction| comp.)
      //   3 frozen: exact intersection (unavoidable)
      // Only the dependent coordinates carry the frozen planes' coefficients; growth is minimal.
      auto make_anchor = [&](const VertexSPtr& v, const FacetSPtr& f_hint)
      {
        const Point_3& p = original_points[v];

        boost::container::small_vector<FacetSPtr, 3> S(frozen_at[v].begin(), frozen_at[v].end());
        std::sort(S.begin(), S.end(), [](const FacetSPtr& a, const FacetSPtr& b) { return a->id() < b->id(); });

        const double delta = std::sqrt(CGAL::to_double(sq_max_displacements[v]));
        const double slide = std::min(nudge_range, 0.25 * (1.0 - kappa) * delta);
        std::uniform_real_distribution<> sdist(-slide, slide);

        // double target: original point + random slide
        double tx = CGAL::to_double(p.x()) + sdist(gen);
        double ty = CGAL::to_double(p.y()) + sdist(gen);
        double tz = CGAL::to_double(p.z()) + sdist(gen);

        auto coeffs_d = [](const Plane_3& pl, double& a, double& b, double& c, double& d) {
          a = CGAL::to_double(pl.a()); b = CGAL::to_double(pl.b());
          c = CGAL::to_double(pl.c()); d = CGAL::to_double(pl.d());
        };

        Point_3 q;

        if (S.size() == 0) {
          // stay close to the requesting facet's current plane (limits cascades), in double only
          double a, b, c, d; coeffs_d(f_hint->get_plane(), a, b, c, d);
          const double n2 = a*a + b*b + c*c;
          const double px = CGAL::to_double(p.x()), py = CGAL::to_double(p.y()), pz = CGAL::to_double(p.z());
          const double dist = (a*px + b*py + c*pz + d) / n2;
          if (dist * dist * n2 <= kappa * kappa * delta * delta) {
            tx -= dist * a; ty -= dist * b; tz -= dist * c;
          }
          q = Point_3(FT(tx), FT(ty), FT(tz));
        }
        else if (S.size() == 1) {
          const Plane_3& pl = S[0]->get_plane();
          double ad, bd, cd, dd; coeffs_d(pl, ad, bd, cd, dd);
          const double n2 = ad*ad + bd*bd + cd*cd;
          CGAL_assertion(n2 > 0.0);

          // project the double target onto the plane in double
          const double dist = (ad*tx + bd*ty + cd*tz + dd) / n2;
          const double px = tx - dist*ad, py = ty - dist*bd, pz = tz - dist*cd;

          const FT& a = pl.a(); const FT& b = pl.b(); const FT& c = pl.c(); const FT& d = pl.d();
          const double aa = std::abs(ad), ab = std::abs(bd), ac = std::abs(cd);

          if (ac >= aa && ac >= ab) {
            FT x(px), y(py);
            q = Point_3(x, y, -(a*x + b*y + d) / c);
          } else if (ab >= aa && ab >= ac) {
            FT x(px), z(pz);
            q = Point_3(x, -(a*x + c*z + d) / b, z);
          } else {
            FT y(py), z(pz);
            q = Point_3(-(b*y + c*z + d) / a, y, z);
          }
          CGAL_assertion(pl.has_on(q));
        }
        else if (S.size() == 2) {
          const Plane_3& p1 = S[0]->get_plane();
          const Plane_3& p2 = S[1]->get_plane();
          double a1d, b1d, c1d, d1d; coeffs_d(p1, a1d, b1d, c1d, d1d);
          double a2d, b2d, c2d, d2d; coeffs_d(p2, a2d, b2d, c2d, d2d);

          // line direction u = n1 x n2 (double)
          const double uxd = b1d*c2d - c1d*b2d;
          const double uyd = c1d*a2d - a1d*c2d;
          const double uzd = a1d*b2d - b1d*a2d;
          const double detG = uxd*uxd + uyd*uyd + uzd*uzd;
          CGAL_assertion(detG > 0.0);

          // orthogonal projection of the double target onto the line, in double
          const double g11 = a1d*a1d + b1d*b1d + c1d*c1d;
          const double g22 = a2d*a2d + b2d*b2d + c2d*c2d;
          const double g12 = a1d*a2d + b1d*b2d + c1d*c2d;
          const double r1 = a1d*tx + b1d*ty + c1d*tz + d1d;
          const double r2 = a2d*tx + b2d*ty + c2d*tz + d2d;
          const double alpha = (r1*g22 - r2*g12) / detG;
          const double beta  = (r2*g11 - r1*g12) / detG;
          const double px = tx - alpha*a1d - beta*a2d;
          const double py = ty - alpha*b1d - beta*b2d;
          const double pz = tz - alpha*c1d - beta*c2d;

          const FT& a1 = p1.a(); const FT& b1 = p1.b(); const FT& c1 = p1.c(); const FT& d1 = p1.d();
          const FT& a2 = p2.a(); const FT& b2 = p2.b(); const FT& c2 = p2.c(); const FT& d2 = p2.d();
          const FT ux = b1*c2 - c1*b2;
          const FT uy = c1*a2 - a1*c2;
          const FT uz = a1*b2 - b1*a2;

          const double aux = std::abs(uxd), auy = std::abs(uyd), auz = std::abs(uzd);

          if (auz >= aux && auz >= auy) {
            FT z(pz);
            q = Point_3((ux*z + (b1*d2 - d1*b2)) / uz,
                        (uy*z + (d1*a2 - a1*d2)) / uz,
                        z);
          } else if (auy >= aux && auy >= auz) {
            FT y(py);
            q = Point_3((ux*y + (d1*c2 - c1*d2)) / uy,
                        y,
                        (uz*y + (a1*d2 - d1*a2)) / uy);
          } else {
            FT x(px);
            q = Point_3(x,
                        (uy*x + (c1*d2 - d1*c2)) / ux,
                        (uz*x + (d1*b2 - b1*d2)) / ux);
          }
          CGAL_assertion(p1.has_on(q));
          CGAL_assertion(p2.has_on(q));
        }
        else { // S.size() == 3
          std::optional<Point_3> x = Kernel_wrapper::intersection(S[0]->get_plane(),
                                                                  S[1]->get_plane(),
                                                                  S[2]->get_plane());
          CGAL_assertion(x.has_value());
          q = *x;
        }

        const FT sqd = CGAL::squared_distance(q, p);
        CGAL_assertion(sqd < sq_max_displacements[v]);

        const double used = std::sqrt(CGAL::to_double(sqd));
        const double rem = 0.999 * (delta - used);
        CGAL_assertion(rem > 0.0);
        sq_budget[v] = FT(rem * rem);

        CGAL_SS3_TRANSF_TRACE_V(8, "New anchor V" << v->id() << " at " << q << " [" << S.size() << " frozen]");
        CGAL_SS3_TRANSF_TRACE_V(8, "  Distance from original: " << used << ", remaining budget: " << rem);

        v->set_point(q);
      };
#else
      auto make_anchor = [&](const VertexSPtr& v, const FacetSPtr& f_hint)
      {
        const Point_3& p = original_points[v];

        boost::container::small_vector<FacetSPtr, 3> S(frozen_at[v].begin(), frozen_at[v].end());
        std::sort(S.begin(), S.end(), [](const FacetSPtr& a, const FacetSPtr& b) { return a->id() < b->id(); });

        const double delta = std::sqrt(CGAL::to_double(sq_max_displacements[v]));
        // slack left after the frozen structure (kappa*delta); keep the random part well inside
        const double slide = std::min(nudge_range, 0.25 * (1.0 - kappa) * delta);
        std::uniform_real_distribution<> sdist(-slide, slide);

        Point_3 q;
        if (S.size() == 3) {
          q = *Kernel_wrapper::intersection(S[0]->get_plane(), S[1]->get_plane(), S[2]->get_plane());
        } else if (S.size() == 2) {
          Line_3 L = *Kernel_wrapper::intersection(S[0]->get_plane(), S[1]->get_plane());
          const Vector_3 dir = L.to_vector();
          const double inv_len = 1.0 / std::sqrt(CGAL::to_double(dir.squared_length()));
          q = L.projection(p) + FT(sdist(gen) * inv_len) * dir; // exactly on the line
        } else if (S.size() == 1) {
          const Plane_3& P = S[0]->get_plane();
          q = P.projection(p + Vector_3(FT(sdist(gen)), FT(sdist(gen)), FT(sdist(gen)))); // exactly on the plane
        } else {
          // No frozen structure: stay close to the current plane of the requesting facet
          // (limits the cascade), fall back to a random point near p.
          const Plane_3& P = f_hint->get_plane();
          const Point_3 proj = P.projection(p);
          if (CGAL::squared_distance(proj, p) <= FT(kappa * kappa) * sq_max_displacements[v]) {
            q = proj + Vector_3(FT(sdist(gen)), FT(sdist(gen)), FT(sdist(gen)));
          } else {
            q = p + Vector_3(FT(sdist(gen)), FT(sdist(gen)), FT(sdist(gen)));
          }
        }

        const FT sqd = CGAL::squared_distance(q, p);
        CGAL_assertion(sqd < sq_max_displacements[v]);

        // remaining budget from the anchor (conservative, double)
        const double used = std::sqrt(CGAL::to_double(sqd));
        const double rem = 0.999 * (delta - used);
        CGAL_assertion(rem > 0.0);
        sq_budget[v] = FT(rem * rem);

        CGAL_SS3_TRANSF_TRACE_V(8, "New anchor V" << v->id() << " at " << q << " [" << S.size() << " frozen]");
        CGAL_SS3_TRANSF_TRACE_V(8, "  Distance from original: " << used << ", remaining budget: " << rem);

        v->set_point(q);
      };
#endif

      std::vector<Vertex_stability_record> unstable_vertices;
      CGAL::unordered_flat_map<VertexSPtr, std::size_t> unstable_record_indices;

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream out_unstable_base("results/base_unstable_vertices.xyz");
      out_unstable_base.precision(17);
#endif

      CGAL_SS3_TRANSF_TRACE_V(16, "Initial stable/unstable classification");
      for (const VertexSPtr& v : polyhedron->vertices()) {
        const Vertex_stability_info info = evaluate_vertex_stability(v, sq_budget[v]);
        if (!info.is_stable) {
#ifdef CGAL_SS3_DUMP_FILES
          out_unstable_base << v->point() << "\n";
#endif
          unstable_record_indices[v] = unstable_vertices.size();
          unstable_vertices.push_back({v, info.worst_facet, info.max_sq_displacement, true, static_cast<std::size_t>(v->id())});
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      out_unstable_base.close();
#endif

      using Stability_queue = Modifiable_priority_queue<std::size_t,
                                                        Vertex_stability_record_less,
                                                        boost::identity_property_map,
                                                        CGAL_BOOST_PAIRING_HEAP>;
      Stability_queue pq(polyhedron->vertices().size(), Vertex_stability_record_less(&unstable_vertices));
      for (std::size_t i = 0; i < unstable_vertices.size(); ++i) {
        pq.push(i);
      }

      auto update_vertex_record = [&](const VertexSPtr& v)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "update record of V" << v->id());
        const Vertex_stability_info info = evaluate_vertex_stability(v, sq_budget[v]);

        const auto it = unstable_record_indices.find(v);
        if (it != unstable_record_indices.end()) {
          Vertex_stability_record& record = unstable_vertices[it->second];
          record.vertex = v;
          record.vertex_id = v->id();
          record.facet = info.worst_facet;
          record.max_sq_displacement = info.max_sq_displacement;
          record.is_active = !info.is_stable;
          if (info.is_stable) {
            CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is now stable");
            if (pq.contains(it->second)) {
              pq.erase(it->second);
            }
            record.is_active = false;
            record.facet = nullptr;
            record.max_sq_displacement = FT(0);
          } else {
            record.is_active = true;
            record.facet = info.worst_facet;
            record.max_sq_displacement = info.max_sq_displacement;
            if (pq.contains(it->second)) {
              pq.update(it->second);
            } else {
              pq.push(it->second);
            }
          }
        } else if (!info.is_stable) {
          const std::size_t index = unstable_vertices.size();
          CGAL_assertion_code(for (const auto& r : unstable_vertices))
          CGAL_assertion(r.vertex == v || r.vertex_id != static_cast<std::size_t>(v->id()));
          unstable_vertices.push_back({v, info.worst_facet, info.max_sq_displacement, true, static_cast<std::size_t>(v->id())});
          unstable_record_indices[v] = index;
          CGAL_SS3_TRANSF_TRACE_V(64, "new unstable vertex");
          pq.push(index);
        }
      };

      // main loop
      while (!pq.empty()) {
        CGAL_SS3_TRANSF_TRACE_V(16, "stable/unstable: main loop (" << pq.size() << ")");

        const std::size_t unstable_index = pq.top_and_pop();
        Vertex_stability_record record = unstable_vertices[unstable_index];
        if (!record.is_active) {
          continue;
        }

        const VertexSPtr& v = record.vertex;
        const FacetSPtr& f = record.facet;
        CGAL_SS3_TRANSF_TRACE_V(16, "pop V" << v->id() << " / F" << f->id());

        CGAL_assertion(frozen.count(f) == 0);
        CGAL_assertion(!is_anchor_of(v, f));

        // (3) place the anchor (once per vertex), register it in f
        if (all_anchors.insert(v).second) {
          make_anchor(v, f);
        }
        auto& facet_anchors = anchors[f];
        facet_anchors.push_back(v);

        std::vector<FacetSPtr> facets_to_recompute;

        if (facet_anchors.size() > 3) {
          CGAL_assertion(!f->is_triangle());
          CGAL_SS3_TRANSF_TRACE_V(8, "  Must triangulate F" << f->id());

          const boost::container::small_vector<VertexSPtr, 3> old_anchors(facet_anchors.begin(), facet_anchors.end());
          const Plane_3 original_plane = original_planes.at(f);

          auto [_, new_facets] = Transformation::triangulate_facet(f, polyhedron);
          original_planes.erase(f);
          anchors.erase(f);

          for (const FacetSPtr& nf : new_facets) {
            original_planes[nf] = original_plane;
            // backport anchors: a triangle inherits the anchors among its vertices
            for (const VertexSPtr& w : nf->vertices()) {
              if (std::find(old_anchors.begin(), old_anchors.end(), w) != old_anchors.end()) {
                anchors[nf].push_back(w);
              }
            }
            facets_to_recompute.push_back(nf);
          }
        } else {
          facets_to_recompute.push_back(f);
        }

        for (const FacetSPtr& nf : facets_to_recompute) {
          // (4) recompute the facet equation
          recompute_facet_plane(nf);

          // (5) update the facet's vertices stability in the queue
          for (const VertexSPtr& w : nf->vertices()) {
            update_vertex_record(w);
          }
        }
      }

      // -----------------------------------------------------------------------
      // Exact certification: stability w.r.t. the ORIGINAL points + general position
      // -----------------------------------------------------------------------
      for (const VertexSPtr& v : polyhedron->vertices()) {
        std::vector<FacetSPtr> ifs;
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            ifs.push_back(f);
          }
        }
        for (std::size_t i=0; i<ifs.size(); ++i) {
          for (std::size_t j=i+1; j<ifs.size(); ++j) {
            if (CGAL::cross_product(ifs[i]->get_plane().orthogonal_vector(),
                                    ifs[j]->get_plane().orthogonal_vector()) == CGAL::NULL_VECTOR) {
              CGAL_SS3_TRANSF_TRACE_V(1, "Error: parallel planes at V" << v->id());
              CGAL_SS3_TRANSF_TRACE_V(1, "Facets: F" << ifs[i]->id() << " and F" << ifs[j]->id());
              std::abort();
            }
            for (std::size_t k=j+1; k<ifs.size(); ++k) {
              std::optional<Point_3> x = Kernel_wrapper::intersection(ifs[i]->get_plane(), ifs[j]->get_plane(), ifs[k]->get_plane());
              if (!x) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Error: singular triplet at V" << v->id());
                CGAL_SS3_TRANSF_TRACE_V(1, "Facets: F" << ifs[i]->id() << ", F" << ifs[j]->id() << ", and F" << ifs[k]->id());
                std::abort();
              }
              if (CGAL::squared_distance(*x, original_points[v]) > sq_max_displacements[v]) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Error: stability violated at V" << v->id());
                std::abort();
              }
            }
          }
        }
      }
#endif

#ifdef CGAL_SS3_PERTURB_V4_PART1_V2
      CGAL_SS3_TRANSF_TRACE_V(16, "Initial plane perturbations");
      for (const FacetSPtr& f : polyhedron->facets()) {
        perturbPlaneCoefficientsFixedPoints(f, nudge_range, std::vector<VertexSPtr>());
      }

      // unstable vertex --> first 3 incident facets determining the position of the vertex
      CGAL::unordered_flat_map<VertexSPtr, CGAL::unordered_flat_set<FacetSPtr> > determining_facets;

      // fixed vertices/facets
      std::unordered_set<VertexSPtr> fixed_vertices;
      std::unordered_set<FacetSPtr> fixed_facets;

      // facet --> (determined) unstable vertices
      //
      // The facet becomes fixed at 2 vertices and not 3 vertices despite the vertices being perturbed
      // because if we do 3 random perturbations of vertices, the normal can vary wildly.
      //
      // @todo Ideally, it could be fixed with 3 unstable vertices and a smarter perturbation
      // (which takes into account all incident facets of these 3 fixing vertices...)
      CGAL::unordered_flat_map<FacetSPtr, CGAL::unordered_flat_set<VertexSPtr> > determining_vertices;

      // -------------------------------------------------------------------------------------------

      auto has_unstable_vertices = [&](const FacetSPtr& f) -> bool
      {
        CGAL_SS3_TRANSF_TRACE_V(64, "Checking if F" << f->id() << " (" << f->vertices().size() << " nv) has unstable vertices");
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_stable(v, sq_max_displacements[v])) {
            return true;
          }
        }
        return false;
      };

      // -------------------------------------------------------------------------------------------

      auto is_vertex_fully_determined = [&](const VertexSPtr& v) -> bool
      {
        CGAL_SS3_TRANSF_TRACE_V(64, "Checking if V" << v->id() << " (deg: " << v->degree() << ") is fully determined [" << determining_facets[v].size() << "]");
        auto it = determining_facets.find(v);
        return (it != determining_facets.end() && it->second.size() == 3);
      };

      // -------------------------------------------------------------------------------------------

      auto is_facet_fully_determined = [&](const FacetSPtr& f) -> bool
      {
        CGAL_SS3_TRANSF_TRACE_V(64, "Checking if F" << f->id() << " (" << f->vertices().size() << " nv) is fully determined");
        CGAL_SS3_TRANSF_TRACE_V(64, "  number of determining vertices: " << determining_vertices[f].size());
        CGAL_assertion(determining_vertices[f].size() <= 3);
        return (f->is_triangle() && determining_vertices[f].size() == 3) ||
                (!f->is_triangle() && determining_vertices[f].size() == 2);
      };

      // -------------------------------------------------------------------------------------------

      auto add_determining_facet = [&](const VertexSPtr& v, const FacetSPtr& f)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "  Determine V" << v->id() << " with F" << f->id());
        determining_facets[v].insert(f);
        CGAL_precondition(determining_facets[v].size() <= 3);
      };

      // -------------------------------------------------------------------------------------------

      auto is_facet_fixed = [&](const FacetSPtr& f) -> bool
      {
        return fixed_facets.find(f) != fixed_facets.end();
      };

#if 0
      // -------------------------------------------------------------------------------------------
      // replaced by 'fix_facet_and_check' because we need to check for newly unstable vertices after perturbation

      auto fix_facet = [&](const FacetSPtr& f)
      {
        CGAL_precondition(is_facet_fully_determined(f));
        CGAL_precondition(!is_facet_fixed(f));

        if (f->is_triangle()) {
          CGAL_assertion(determining_vertices[f].size() == 3); // just to be clear

          fixed_facets.insert(f);

          // for triangles, all vertices are determined, and there is nothing to nudge
          // (note that vertices were themselves nudged, thus the facet is nudged).
          f->init_plane();
          Transformation::normalize_facet_plane(f);
          // Here we do not need to add the fixed facet to incident determined vertices
          // because all vertices are already fully determined
          return;
        }

        perturbPlaneCoefficientsFixedPoints(f, nudge_range, determining_vertices[f]);

        CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " now has plane " << f->get_plane() << " [measure=" << Size_shenanigans::length(f->get_plane()) << "]");

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "F" << f->id() << " is now a fixed face, by");
        CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : determining_vertices[f]))
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->id() << " [measure=" << Size_shenanigans::length(fv->point()) << "]");
        CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

        CGAL_SS3_TRANSF_TRACE_V(32, "Newly fixed facet F" << f->id() << " determines its unstable incident vertices...");

        for (const VertexSPtr& v : f->vertices()) {
          CGAL_SS3_TRANSF_TRACE_V(32, "incident V" << v->id() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facet(s))");
          if (!is_vertex_fully_determined(v)) {
            add_determining_facet(v, f);
          }
        }

        fixed_facets.insert(f);
      };
#endif

      // -------------------------------------------------------------------------------------------

      // This is the main list of facets that we will process
      std::list<FacetSPtr> facets_to_process;

      auto fix_facet_and_check = [&](const FacetSPtr& f)
      {
        CGAL::unordered_flat_map<VertexSPtr, bool> previously_stable;
        CGAL::unordered_flat_set<VertexSPtr> f_vertices;
        for (const VertexSPtr& v : f->vertices()) {
          previously_stable[v] = is_stable(v, sq_max_displacements[v]);
          f_vertices.insert(v);
        }

        perturbPlaneCoefficientsFixedPoints(f, nudge_range, determining_vertices[f]);

        bool should_triangulate = false;
        for (const VertexSPtr& v : f->vertices()) {
          bool was_stable_before = previously_stable[v];
          bool is_stable_now = is_stable(v, sq_max_displacements[v]);

          if (was_stable_before && !is_stable_now) {
            CGAL_SS3_TRANSF_TRACE_V(32, "  Vertex V" << v->id() << " became unstable due to perturbation of F" << f->id());
            should_triangulate = true;

            // Check all incident facets of this newly unstable vertex
            for (FacetWPtr wf : v->facets()) {
              if (FacetSPtr f_inc = wf.lock()) {
                if (f_inc == f) {
                  continue; // Skip the current facet being perturbed
                }

                if (is_facet_fixed(f_inc)) {
                  // nothing else to do, and the vertex should already be determined by the facet
                  CGAL_assertion(determining_facets[v].count(f_inc));
                } else {
                  // If an incident facet did not have any unstable vertices, then it must be added
                  // to facets being processed
                  bool had_no_unstable_vertices = true;
                  for (const VertexSPtr& u : f_inc->vertices()) {
                    bool u_was_stable = false;
                    if (f_vertices.count(u)) {
                      u_was_stable = previously_stable[u];
                    } else {
                      u_was_stable = is_stable(u, sq_max_displacements[u]);
                    }

                    if (!u_was_stable) {
                      had_no_unstable_vertices = false;
                      break;
                    }
                  }

                  // If that incident facet had unstable vertices and is not fixed, then it should
                  // already be in the queue
                  if (had_no_unstable_vertices) {
                    // The main loop queue only processes non-triangles (it asserts !is_triangle())
                    if (!f_inc->is_triangle()) {
                      CGAL_assertion(std::find(facets_to_process.begin(), facets_to_process.end(), f_inc) == facets_to_process.end());

                      CGAL_SS3_TRANSF_TRACE_V(32, "  Adding F" << f_inc->id() << " to queue because of newly unstable V" << v->id());
                      facets_to_process.push_back(f_inc);
                    }
                  }
                }
              }
            }
          }
        }

        if (should_triangulate) {
          CGAL_SS3_TRANSF_TRACE_V(32, "  Triangulating F" << f->id() << " because a vertex became unstable");
          triangulate_facet(f, polyhedron);
        }
      };

      // -------------------------------------------------------------------------------------------

      auto add_determining_vertex = [&](const FacetSPtr& f, const VertexSPtr& v)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "  Determine F" << f->id() << " with V" << v->id());

        if (is_facet_fully_determined(f)) {
          CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f->id() << " is already fixed");
          return;
        }

        determining_vertices[f].insert(v);

        if (is_facet_fully_determined(f)) {
          fix_facet(f);
        }
      };

      // -------------------------------------------------------------------------------------------

      // @fixme can this move a vertex 'v' too far away from its original position when
      // it is constrained by 1 or 2 fixed facets?
      auto nudge_constrained_vertex = [&](const VertexSPtr& v)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "  Nudging V" << v->id() << " from " << v->point());

        std::vector<const Plane_3*> constraining_planes;
        for (const FacetSPtr& df : determining_facets[v]) {
          if (is_facet_fully_determined(df)) {
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

        std::array<double, 3> v_r = rand_vec(-nudge_range/2.0, nudge_range/2.0);
        double x = CGAL::to_double(p.x()) + v_r[0];
        double y = CGAL::to_double(p.y()) + v_r[1];
        double z = CGAL::to_double(p.z()) + v_r[2];

        Point_3 p_nudged { x, y, z };
        CGAL_SS3_TRANSF_TRACE_V(32, "base nudge: " << x << " " << y << " " << z);

        Point_3 p_new;
        if (n_fixed == 0) {
          p_new = p_nudged;
        } else if (n_fixed == 1) {
          const Plane_3& plane = *(constraining_planes[0]);
          p_new = plane.projection(p_nudged);
        } else if (n_fixed == 2) {
          const Plane_3& plane1 = *(constraining_planes[0]);
          const Plane_3& plane2 = *(constraining_planes[1]);
          std::optional<Line_3> line = Kernel_wrapper::intersection(plane1, plane2);
          p_new = line->projection(p_nudged);
        }

        CGAL_SS3_TRANSF_TRACE_V(32, "  Nudged V" << v->id() << " to " << p_new);

        v->set_point(p_new);
      };

      // -------------------------------------------------------------------------------------------

      auto is_vertex_fixed = [&](const VertexSPtr& v) -> bool
      {
        return fixed_vertices.find(v) != fixed_vertices.end();
      };

      // -------------------------------------------------------------------------------------------

      auto fix_vertex = [&](const VertexSPtr& v)
      {
        CGAL_precondition(!is_stable(v, sq_max_displacements[v]));
        CGAL_precondition(is_vertex_fully_determined(v));
        CGAL_precondition(!is_vertex_fixed(v));

        CGAL_SS3_TRANSF_TRACE_CODE(auto it = determining_facets[v].begin();)
        CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is now fully determined by"
                                    << " F" << (*it)->id() << " [measure=" << Size_shenanigans::length((*it)->get_plane())
                                    << "] F" << (*std::next(it))->id() << " [measure=" << Size_shenanigans::length((*std::next(it))->get_plane())
                                    << "] F" << (*std::next(it, 2))->id() << " [measure=" << Size_shenanigans::length((*std::next(it, 2))->get_plane()) << "]");

        // set the nudged position for the vertex: a nudge constrained by the fixed incident facets
        nudge_constrained_vertex(v);

        CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is now a fixed vertex at " << v->point() << " [measure=" << Size_shenanigans::length(v->point()) << "]");

        fixed_vertices.insert(v);

        // a fixed vertex determines its incident facets
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            add_determining_vertex(f, v);
          }
        }
      };

      // -------------------------------------------------------------------------------------------

      auto is_facet_overconstrained = [&](const FacetSPtr& f) -> bool
      {
        if (f->is_triangle() || is_facet_fully_determined(f)) {
          return false;
        }

        // @todo
        // the facet is overconstrained if it has more than 2 unstable vertices that are either:
        // - already fixed, or
        // - would be fixed if we were to fix this facet (i.e. they have 2 determining facets and are incident to this facet)

        // @fixme this isn't not a good test here, because the facet's plane changes when it is determined.
        // Secondly, can we do subset-of-facets-stability?
        // Anyhow, maybe the best way to proceed is to not have this function in this pipeline.
        // Instead, simply do fix_facet(), and then check afterwards: have we over-constrained
        // the adjacent facets (and the facet itself)? If so, backtrack and triangulate the facet
        // instead of fixing it.
        //
        // we cannot fix that facet if determining unstable vertices of this facet
        // would create too many fixed vertices in *any* unfixed facet incident to
        // the determined (unstable) vertices of this facet
        CGAL::unordered_flat_map<FacetSPtr, unsigned int> facets_to_test; // facets --> number of appearances
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_stable(v, sq_max_displacements[v])) {
            for (FacetWPtr inc_f : v->facets()) {
              if (FacetSPtr f = inc_f.lock()) {
                if (!is_facet_fully_determined(f)) {
                  ++facets_to_test[f];
                }
              }
            }
          }
        }

        for (const auto& [ft, count] : facets_to_test) {
          // Count the number of unstable vertices with either:
          // - 3 determining facets
          // - 2 determining facets and incident to 'facet'
          // These are vertices that are fixed, or would be fixed once we 'add' the facet
          // to its unstable vertices.
          unsigned int constrain_n = 0;
          for (const VertexSPtr& v : f->vertices()) {
            if (!is_stable(v, sq_max_displacements[v])) {
              if (is_vertex_fixed(v)) {
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

      // -------------------------------------------------------------------------------------------

      // sometimes we could fix as a polygon, but we need triangulate for other reasons
      auto should_triangulate_facet = [&](const FacetSPtr& f) -> bool
      {
        if (f->is_triangle() || is_facet_fully_determined(f)) {
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

      // -------------------------------------------------------------------------------------------

      auto triangulate_facet = [&](const FacetSPtr& facet_tt)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Triangulate F" << facet_tt->id());

        CGAL_assertion(!is_facet_fully_determined(facet_tt));

        // the facet is not yet fixed, so no vertex can have it as determining facet
        CGAL_assertion_code(for (const VertexSPtr& v : facet_tt->vertices()) {)
        CGAL_assertion(determining_facets[v].size() <= 3);
        CGAL_assertion(determining_facets[v].count(facet_tt) == 0);
        CGAL_assertion_code(})

        auto [local_vertices, new_facets] = Transformation::triangulate_facet(facet_tt, polyhedron);

        for (const VertexSPtr& v : local_vertices) {
          CGAL_SS3_TRANSF_TRACE_V(64, "local vertex " << v->id() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)");

          if (is_vertex_fixed(v)) {
            CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is already fixed, skipping");
            continue;
          }

          for (FacetWPtr wf : v->facets()) {
            if (FacetSPtr fptr = wf.lock()) {
              if (is_facet_fully_determined(fptr)) {
                CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is determined by F" << fptr->id() << " (c)");
                add_determining_facet(v, fptr);

                if (is_vertex_fully_determined(v)) {
                  fix_vertex(v);
                  break;
                }
              }
            }
          }
        }

        // already-fixed vertices are fixed points for the new facets
        for (const FacetSPtr& nf : new_facets) {
          CGAL_SS3_TRANSF_TRACE_V(32, "spawned F" << nf->id());

          for (const VertexSPtr& iv : nf->vertices()) {
            if (is_vertex_fixed(iv)) {
              CGAL_SS3_TRANSF_TRACE_V(64, "newborn F" << nf->id() << " is constrained by V" << iv->id());
              add_determining_vertex(nf, iv);
            }
          }
        }
      };

      // -------------------------------------------------------------------------------------------

      auto facet_sorter = [&](const FacetSPtr& a, const FacetSPtr& b)
      {
        auto uv_count = [&](const FacetSPtr& f) -> unsigned int {
          unsigned int uv_n = 0;
          for (const VertexSPtr& v : f->vertices()) {
            if (!is_stable(v, sq_max_displacements[v])) {
              ++uv_n;
            }
          }
          return uv_n;
        };

        // Give priority to facets with no fixed vertices as to avoid having to triangulate.
        // If both or neither have constrained vertices, give priority to the largest unstable count.
        //
        // The point is to avoid cascading exact number types, even if we have to triangulate a little more
        auto get_fixed_count = [&](const FacetSPtr& f) -> unsigned int
        {
          unsigned int res = 0;
          for (const VertexSPtr& v : f->vertices()) {
            if (is_vertex_fixed(v)) {
              ++res;
            }
          }
          return res;
        };

        unsigned int adn = get_fixed_count(a);
        unsigned int bdn = get_fixed_count(b);

        CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->id() << " has " << adn << " fully determined vertices");
        CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->id() << " has " << bdn << " fully determined vertices");

        if (adn != bdn) {
          // Give priority to the one with the least amount of fully determined vertices
          return adn < bdn;
        }

        // same number of fully determined vertices, give priority to the facet with the most unstable vertices
        unsigned int a_hdv_n = uv_count(a);
        unsigned int b_hdv_n = uv_count(b);

        CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->id() << " has " << a_hdv_n << " unstable vertice(s)");
        CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->id() << " has " << b_hdv_n << " unstable vertice(s)");

        if (a_hdv_n != b_hdv_n) {
          // Give priority to the one with the most unstable vertices
          return a_hdv_n > b_hdv_n;
        }

        // same number of fully determined vertices and unstable vertices, give priority to the largest facet
        return a->vertices().size() > b->vertices().size();
      };

      // -------------------------------------------------------------------------------------------

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream unst_out("results/unstable_vertices.xyz");
      unst_out.precision(17);
      for (const VertexSPtr& v : polyhedron->vertices()) {
        if (!is_stable(v, sq_max_displacements[v])) {
          unst_out << v->point() << "\n";
        }
      }
      unst_out.close();
#endif

      for (const FacetSPtr& f : polyhedron->facets()) {
        if (f->is_triangle() || !has_unstable_vertices(f)) {
          continue;
        }
        facets_to_process.push_back(f);
      }

      CGAL_SS3_TRANSF_TRACE_V(16, "== Main queue... ==");
      CGAL_SS3_TRANSF_TRACE_V(16, "  " << facets_to_process.size() << " initial facets to process");

      while (!facets_to_process.empty()) {
        CGAL_SS3_TRANSF_TRACE_V(16, "Sort again...");
        facets_to_process.sort(facet_sorter); // @todo use a priority queue

        FacetSPtr f = facets_to_process.front();
        facets_to_process.pop_front();
        CGAL_SS3_TRANSF_TRACE_V(16, "Pop F" << f->id());

        CGAL_assertion(!f->is_triangle());
        CGAL_assertion(!is_facet_fixed(f));

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "  " << determining_vertices[f].size() << " determining vertices:";)
        CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : determining_vertices[f]))
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->id();)
        CGAL_SS3_TRANSF_TRACE_V(32, ss.str());
        CGAL_assertion(determining_vertices[f].size() <= 2);

        if (is_facet_overconstrained(f) || should_triangulate_facet(f)) {
          CGAL_SS3_TRANSF_TRACE_V(32, "  Must triangulate F" << f->id());
          triangulate_facet(f);
          continue;
        }

        // Now, adding the facet to the unstable vertices will not over-constrain the facet, so do it:
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_vertex_fixed(v)) {
            add_determining_facet(v, f);
            CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " is determined by F" << f->id() << " (d)");
            if (is_vertex_fully_determined(v) && !is_stable(v, sq_max_displacements[v])) {
              // When the vertex becomes fixed (its 3 determining facets are known), we need to:
              // - perturb the position of the vertex
              // - update all incident facets to check if they are now fully determined and if it is the case
              //   recompute their plane coefficients
              fix_vertex(v);
            }
          }
        }
      }

      // Some facets might be unfixed and with unstable vertices, fix them now
      CGAL_SS3_TRANSF_TRACE_V(16, "== Deal with remaining facets with unstable vertices... ==");

      for (const FacetSPtr& f : polyhedron->facets()) {
        if (f->is_triangle() || is_facet_fixed(f) || !has_unstable_vertices(f)) {
          continue;
        }

        CGAL_assertion(!is_facet_fully_determined(f));

        CGAL_SS3_TRANSF_TRACE_V(32, "Nudge and fix F" << f->id() << " [remaining]");

        // @fixme this can create facets with unstable vertices, which will need to be treated
        // in the main loop or in this loop... Some kind of for(;;) is necessary...
        fix_facet_and_check(f);

        // fixing the facet cannot determine a vertex because the facet has already been visited in the main loop
        fixed_facets.insert(f);
      }

      // At this point, only unstable triangles remain
      for (const FacetSPtr& f : polyhedron->facets()) {
        if (f->is_triangle() || !has_unstable_vertices(f)) {
          continue;
        }
        CGAL_assertion(is_facet_fixed(f));
      }

      // Nudge vertices that can still be nudged, for randomness
      CGAL_SS3_TRANSF_TRACE_V(16, "== Nudge non-fixed, unstable vertices... ==");

      for (const VertexSPtr& v : polyhedron->vertices()) {
        if (!is_vertex_fixed(v) && !is_stable(v, sq_max_displacements[v])) {
          CGAL_SS3_TRANSF_TRACE_V(32, "  V" << v->id() << " is unstable and not fixed, nudge it");
          nudge_constrained_vertex(v);

          fixed_vertices.insert(v);

          // since we know only triangle facets are not yet fixed, we don't need to cascade and check
          // if incident facets become fixed
          for (FacetWPtr wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
              if (!is_facet_fully_determined(f)) {
                // the facet cannot be without unstable vertices since v is unstable
                CGAL_assertion(f->is_triangle());
                determining_vertices[f].insert(v);
              }
            }
          }
        }
      }

      // Now handle triangle faces with unstable vertices
      CGAL_SS3_TRANSF_TRACE_V(16, "== Deal with remaining triangles... ==");

      for (const FacetSPtr& f : polyhedron->facets()) {
        if (!f->is_triangle() || !has_unstable_vertices(f)) {
          continue;
        }

        CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "Fix F" << f->id() << " [");
        CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : f->vertices()) {)
        CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << v->id());
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " (" << v->degree() << ")");
        CGAL_SS3_TRANSF_TRACE_CODE(if (is_vertex_fully_determined(v)) { ss << "*"; })
        CGAL_SS3_TRANSF_TRACE_CODE(ss << " "; } ss << "]";)
        CGAL_SS3_TRANSF_TRACE_V(32, ss.str());

        CGAL_SS3_TRANSF_TRACE_V(32, "Nudge and fix F" << f->id() << " [last]");

        if (determining_vertices[f].size() == 3) {
          f->init_plane();
          Transformation::normalize_facet_plane(f);
        } else {
          fix_facet_and_check(f);
        }

        fixed_facets.insert(f);

        // We still need to update the determining facets because some neighboring
        // facet could be an unfixed unstable face
        for (const VertexSPtr& v : f->vertices()) {
          if (!is_vertex_fully_determined(v)) {
            add_determining_facet(v, f);
            CGAL_SS3_TRANSF_TRACE_V(32, "  V" << v->id() << " is determined by F" << f->id() << " (f)");
            // no need to cascade here, because we know only stable vertices are left
          }
        }

        for (const VertexSPtr& v : f->vertices()) {
          if (!is_facet_fully_determined(f)) // should be just for debugging
            determining_vertices[f].insert(v);
        }

        CGAL_postcondition(is_facet_fully_determined(f));
      }

      CGAL_SS3_TRANSF_TRACE_V(16, "Reset the position of degree 3 vertices...");

      for (const VertexSPtr& v : polyhedron->vertices()) {
        // At this point, high-degree vertices do not live in a single position, but will be split
        if (v->degree() == 3 && !is_vertex_fixed(v)) {
          Transformation::reset_point(v);
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(8, "All facets perturbed");

# if 0
      CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
      CGAL_assertion(is_vertex_fixed(v));
      CGAL_assertion_code(})

      CGAL_assertion_code(for (const FacetSPtr& f : polyhedron->facets()) {)
      CGAL_assertion(is_facet_fixed(f));
      CGAL_assertion_code(})

      CGAL_assertion_code(for (const FacetSPtr& facet : polyhedron->facets()) {)
      CGAL_assertion_code(for (const VertexSPtr& v : facet->vertices()) {)
      CGAL_assertion(facet->get_plane().has_on(v->point()));
      CGAL_assertion_code(})
      CGAL_assertion_code(})
# endif

#endif // CGAL_SS3_PERTURB_V4_PART1_V2

  #ifdef CGAL_SS3_DUMP_FILES
      IO::write_OBJ("results/V4_general_position.obj", polyhedron, parameters::do_not_triangulate_faces(true));
      IO::write_OBJ("results/V4_general_position-triangulated.obj", polyhedron, parameters::do_not_triangulate_faces(false));
  #endif

      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
      CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " has length " << Size_shenanigans::length(v->point()));

      CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
      CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " has length " << Size_shenanigans::length(f->get_plane()));

      CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
      CGAL_assertion(is_stable(v, sq_max_displacements[v]));
      CGAL_assertion_code(})

      CGAL_SS3_TRANSF_TRACE_V(8, "Shift facets incident to high-degree vertices");

      // Now, we have stability and enforced planes go through anchors where needed.
      // Since we have stability and general position, we can now split high-degree vertices
      // by slighlty shifting planes, thus creating an arrangement of planes within which
      // we select the new surface.

      // Only facets incident to a high-degree vertex need this random translation.
      CGAL_SS3_TRANSF_TRACE_CODE(std::size_t hdv = 0;)

      CGAL::unordered_flat_set<FacetSPtr> facets_to_translate;
      for (const VertexSPtr& v : polyhedron->vertices()) {
        if (v->degree() > 3) {
          CGAL_SS3_TRANSF_TRACE_CODE(++hdv;)
          for (const FacetWPtr& wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
              facets_to_translate.insert(f);
            }
          }
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(8, "Shift " << facets_to_translate.size() << " facets incident to " << hdv << " high-degree vertices");

      CGAL::unordered_flat_map<FacetSPtr, Plane_3> pre_nudge_planes;
      CGAL::unordered_flat_map<FacetSPtr, FT> nudge_ranges;

      for (const FacetSPtr& f : facets_to_translate) {
        nudge_ranges[f] = nudge_range;
        pre_nudge_planes[f] = f->get_plane();
      }

      CGAL::unordered_flat_set<FacetSPtr> active_facets = facets_to_translate;
      for (;;) {
        if (active_facets.empty()) {
          CGAL_SS3_TRANSF_TRACE_V(16, "  No facet needs translation reduction; stopping");
          break;
        }

        // compute local set of vertices affected by the active facets this iteration
        CGAL::unordered_flat_set<VertexSPtr> vertices_to_check;
        for (const FacetSPtr& f : active_facets) {
          for (const VertexSPtr& fv : f->vertices()) {
            vertices_to_check.insert(fv);
          }
        }

        CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
        CGAL_assertion_code(if (vertices_to_check.count(v)) continue;)
        CGAL_assertion(is_stable(v, sq_max_displacements[v]));
        CGAL_assertion_code(})

        for (const FacetSPtr& f : active_facets) {
          const Plane_3& base_plane = pre_nudge_planes[f];
          const FT& nudge_coeff = nudge_ranges[f];

          CGAL_SS3_TRANSF_TRACE_V(32, "Perturb (Translate) Facet " << f->id());
          CGAL_SS3_TRANSF_TRACE_V(32, "  Nudge coefficient: " << nudge_coeff);
          CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                              << f->get_plane().c() << " " << f->get_plane().d() << "]");

          Plane_3 nudged_plane { base_plane.a(), base_plane.b(), base_plane.c(), base_plane.d() + nudge_coeff };
          f->set_plane(nudged_plane);

          CGAL_SS3_TRANSF_TRACE_V(32, "  To coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");
        }

        // Now, facets have updated with nudged planes -> check stability of the
        // vertices affected by those facets to ensure the nudge was small enough.
        // Note that once a vertex has become stable again after a reduced nudge,
        // it cannot become unstable again.
        bool all_stable = true;
        CGAL::unordered_flat_set<FacetSPtr> next_active_facets;
        for (const VertexSPtr& v : vertices_to_check) {
          if (is_stable(v, sq_max_displacements[v])) {
            CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is stable => no longer checked");
            continue;
          }

          // @todo we could reduce only the nudge of one facet at a time (one of the facets responsible
          // for stability failure).
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is unstable => reduce nudging of incident facets");
          all_stable = false;

          for (const FacetWPtr& wf : v->facets()) {
            if (FacetSPtr f = wf.lock()) {
              if (!nudge_ranges.count(f)) {
                continue;
              }

              CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f->id() << " selected for nudge reduction");
              next_active_facets.insert(f);
            }
          }
        }

        if (all_stable) {
          CGAL_SS3_TRANSF_TRACE_V(8, "  All vertices stable");
          for (const FacetSPtr& f : facets_to_translate) {
            Transformation::normalize_facet_plane(f); // @fixme is this really needed?...
          }
          break;
        }

        // apply a single reduction per facet now that next_active_facets is complete
        for (const FacetSPtr& f : next_active_facets) {
          auto it = nudge_ranges.find(f);
          CGAL_assertion(it != nudge_ranges.end());
          it->second *= FT(0.01);
          CGAL_SS3_TRANSF_TRACE_V(32, "  Nudge range of F" << f->id() << " set to " << it->second);
        }

        // next iteration only needs to consider facets that were reduced/are still active
        active_facets = next_active_facets;
      }

#ifdef CGAL_SS3_DUMP_FILES
        IO::write_OBJ("results/V4_anchored.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

      for (const FacetSPtr& f : polyhedron->facets()) {
        CGAL_SS3_TRANSF_TRACE_V(32, "  Facet #" << f->id() << " ["
                              << f->get_plane().a() << " " << f->get_plane().b() << " "
                              << f->get_plane().c() << " " << f->get_plane().d() << "]");
      }

      CGAL_postcondition(are_planes_in_general_position(polyhedron));

      CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()))
      CGAL_assertion(is_stable(v, sq_max_displacements[v]));
    }

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
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<FacetSPtr> polygon_to_facet;

      // Create a sufficiently large bounding box containing all plane intersections
      CGAL_assertion(is_stable(vertex, sq_max_displacements[vertex]));

      CGAL::Bbox_3 vbb = vertex->point().bbox();
      const double md = approx(approximate_sqrt(sq_max_displacements[vertex])).sup();
      CGAL::Bbox_3 bb = { vbb.xmin() - md, vbb.ymin() - md, vbb.zmin() - md,
                          vbb.xmax() + md, vbb.ymax() + md, vbb.zmax() + md };

      bb.scale(1.5);
      Iso_cuboid_3 bbox { bb };

      CGAL_SS3_TRANSF_TRACE_V(64, "splitting bounding box: " << bbox);
      CGAL_SS3_TRANSF_TRACE_V(64, "x span " << bb.x_span() << ", y span " << bb.y_span() << ", z span " << bb.z_span());

      CGAL_assertion_code(for (auto it_wf1 = vertex->facets().begin() ; it_wf1 != vertex->facets().end(); ++it_wf1) {)
      CGAL_assertion_code(if (FacetSPtr f1 = it_wf1->lock()) {)
      CGAL_assertion_code(for (auto it_wf2 = std::next(it_wf1) ; it_wf2 != vertex->facets().end(); ++it_wf2) {)
      CGAL_assertion_code(if (FacetSPtr f2 = it_wf2->lock()) {)
      CGAL_assertion_code(for (auto it_wf3 = std::next(it_wf2) ; it_wf3 != vertex->facets().end(); ++it_wf3) {)
      CGAL_assertion_code(if (FacetSPtr f3 = it_wf3->lock()) {)
      CGAL_assertion_code(std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());)
      CGAL_assertion_code(if (p_new.has_value()) {)
      CGAL_assertion(bb.xmin() < p_new->x() && p_new->x() < bb.xmax());
      CGAL_assertion(bb.ymin() < p_new->y() && p_new->y() < bb.ymax());
      CGAL_assertion(bb.zmin() < p_new->z() && p_new->z() < bb.zmax());
      CGAL_assertion_code(}}}}}}})

      // Compute the arrangement using PMP::split() and the supporting planes
      // of the facets incident to the vertex.
      //
      // @todo most of what is below could be removed if we had a function that performed
      // split --> LCC. But in addition to this, we need to be able to walk a path of edges
      // during link recovery in the constrain problem setup.
      using Mesh = CGAL::Surface_mesh<Point_3>;
      using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

      Mesh mesh;
      CGAL::make_hexahedron(bb, mesh);
      auto fpm = mesh.template add_property_map<face_descriptor, FacetSPtr>("f:fsptr").first;

      // Start from the splitting box itself; after all the splits, we convert the mesh back to a soup
      // while keeping, on each face, the source facet that generated it.

      for (face_descriptor fd : faces(mesh)) {
        put(fpm, fd, FacetSPtr());
      }

      Face_property_map_updating_coref_visitor<Mesh, decltype(fpm)> visitor(fpm);

      for (const auto& facet_wptr : vertex->facets()) {
        if (FacetSPtr facet = facet_wptr.lock()) {
          const Plane_3& plane = facet->get_plane();
          CGAL::Polygon_mesh_processing::split(mesh, plane, CGAL::parameters::do_not_triangulate_faces(true).visitor(visitor));
          std::vector<halfedge_descriptor> hedges;
          CGAL::extract_boundary_cycles(mesh, std::back_inserter(hedges));
          for (halfedge_descriptor h : hedges) {
            CGAL::Euler::fill_hole(h, mesh);
            const face_descriptor fd = face(h, mesh);
            put(fpm, fd, facet);
          }
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      auto nvpm = get(CGAL::dynamic_vertex_property_t<Point_3>{}, mesh);
      for (vertex_descriptor vd : vertices(mesh)) {
        put(nvpm, vd, Point_3((mesh.point(vd).x() - bbox.xmin()) / bb.x_span(),
                              (mesh.point(vd).y() - bbox.ymin()) / bb.y_span(),
                              (mesh.point(vd).z() - bbox.zmin()) / bb.z_span()));
      }

      CGAL::IO::write_polygon_mesh("results/arr_split_normalized.off", mesh,
                                   CGAL::parameters::stream_precision(17).vertex_point_map(nvpm));
#endif

      CGAL_SS3_TRANSF_TRACE_V(16, "split mesh -> " << vertices(mesh).size() << " points, " << faces(mesh).size() << " faces");

      CGAL_postcondition(is_closed(mesh));
      CGAL_postcondition(is_valid_polygon_mesh(mesh));

      // Convert the split mesh back to a polygon soup while preserving the facet association
      points.clear();
      polygons.clear();
      polygons.reserve(num_faces(mesh));
      PMP::polygon_mesh_to_polygon_soup(mesh, points, polygons);

#ifdef CGAL_SS3_DUMP_FILES
      std::vector<Point_3> normalized_points;
      for (const Point_3& p : points) {
        normalized_points.push_back(Point_3((p.x() - bbox.xmin()) / bb.x_span(),
                                            (p.y() - bbox.ymin()) / bb.y_span(),
                                            (p.z() - bbox.zmin()) / bb.z_span()));
      }
#endif

      polygon_to_facet.resize(polygons.size(), FacetSPtr());

      std::size_t face_index = 0;
      for (face_descriptor fd : faces(mesh)) {
        if (face_index >= polygon_to_facet.size()) {
          break;
        }
        polygon_to_facet[face_index++] = get(fpm, fd);
      }

      CGAL_assertion(polygons.size() == polygon_to_facet.size());

      CGAL_SS3_TRANSF_TRACE_V(16, "split soup -> " << points.size() << " points, " << polygons.size() << " polygons");
      for (std::size_t i=0; i<polygon_to_facet.size(); ++i) {
        CGAL_SS3_TRANSF_TRACE_V(64, "polygon " << i << " ptr: " << (polygon_to_facet[i] ? polygon_to_facet[i]->id() : -1));
      }

      PMP::merge_duplicate_points_in_polygon_soup(points, polygons);

#ifdef CGAL_SS3_DUMP_FILES
      normalized_points.clear();
      for (const Point_3& p : points) {
        normalized_points.push_back(Point_3((p.x() - bbox.xmin()) / bb.x_span(),
                                            (p.y() - bbox.ymin()) / bb.y_span(),
                                            (p.z() - bbox.zmin()) / bb.z_span()));
      }
#endif

      Range_updating_repair_PS_visitor<FacetSPtr> repair_ps_visitor(polygon_to_facet);
      PMP::merge_duplicate_polygons_in_polygon_soup(points, polygons,
                                                    CGAL::parameters::visitor(repair_ps_visitor)
                                                                     .erase_all_duplicates(false)
                                                                     .require_same_orientation(false)
                                                                     .verbose(true));

      CGAL_SS3_TRANSF_TRACE_V(16, "repaired soup -> " << points.size() << " points, " << polygons.size() << " polygons");

      CGAL_assertion(polygons.size() == polygon_to_facet.size());

      // Check and possibly flip the orientation such that the polygon soup is consistently oriented
      // with respect to the original facets
      for (std::size_t i=0; i<polygons.size(); ++i) {
        const FacetSPtr& fsptr = polygon_to_facet[i];
        if (fsptr) {
          const Plane_3& plane = fsptr->get_plane();
          const Point_3& p0 = points[polygons[i][0]];
          const Point_3& p1 = points[polygons[i][1]];
          const Point_3& p2 = points[polygons[i][2]];
          Vector_3 v01 = p1 - p0;
          Vector_3 v02 = p2 - p0;
          Vector_3 normal = CGAL::cross_product(v01, v02);
          if (plane.a() * normal.x() + plane.b() * normal.y() + plane.c() * normal.z() < 0) {
            // The orientation is opposite to the original facet's plane, flip the polygon
            std::reverse(polygons[i].begin(), polygons[i].end());
          }
        }
      }

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_polygon_soup("results/arr_split_normalized.off", normalized_points, polygons, CGAL::parameters::stream_precision(17));
#endif

#ifdef CGAL_SS3_DUMP_FILES
      // Dump a polygon soup for each set of triangle faces associated to a specific fsptr value
      {
        CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_triangles;
        for (std::size_t i=0; i<polygons.size(); ++i) {
          FacetSPtr fsptr = polygon_to_facet[i];
          if (fsptr) {
            id_to_triangles[fsptr->id()].push_back(polygons[i]);
          }
        }

        for (const auto& kv : id_to_triangles) {
          std::ostringstream oss;
          oss << "results/arr_final_" << kv.first << "_normalized.off";
          CGAL::IO::write_OFF(oss.str(), normalized_points, kv.second, CGAL::parameters::stream_precision(17));
        }
      }
#endif

      CGAL_assertion(polygons.size() == polygon_to_facet.size());

#ifdef CGAL_SS3_DUMP_FILES
      // dump only the triangles with a non nullptr
      {
        std::vector<std::vector<std::size_t> > input_triangles;
        for (std::size_t i=0; i<polygons.size(); ++i) {
          if (polygon_to_facet[i]) {
            input_triangles.push_back(polygons[i]);
          }
        }

        CGAL::IO::write_OFF("results/arr_recovered_input.off", normalized_points, input_triangles, CGAL::parameters::stream_precision(17));
      }
#endif

      using PID = std::size_t;
      using FID = std::size_t;
      using VID = std::size_t;

      CGAL_SS3_TRANSF_TRACE_V(16, "fill edge map...");

      std::vector<CGAL::unordered_flat_map<PID, boost::container::small_vector<FID, 4> > > edge_map(points.size());

#define CGAL_SS3_FILL_EDGE_MAP_V2
#ifdef CGAL_SS3_FILL_EDGE_MAP_V2
      // Build edge incidence map directly from the polygon soup (polygons vector).
      // This collects for each canonical edge (min,max) the list of incident polygon ids.
      // Afterwards, handle the common 4-face case by interleaving polygons belonging
      // to two distinct input facets into an a-b-a-b order using `polygon_to_facet`.
      {
        // Collect edges
        for (FID fid = 0; fid < static_cast<FID>(polygons.size()); ++fid) {
          const auto& poly = polygons[fid];
          const std::size_t s = poly.size();
          for (std::size_t j = 0; j < s; ++j) {
            PID a = poly[j];
            PID b = poly[(j+1)%s];
            if (b < a) std::swap(a, b);
            auto& m = edge_map[a];
            auto [it, inserted] = m.try_emplace(b);
            it->second.push_back(fid);
            CGAL_assertion(it->second.size() <= 4);
          }
        }

        // Sort incident faces in CCW order around the edge seen from pid0
        for (PID pid0=0; pid0<static_cast<PID>(edge_map.size()); ++pid0) {
          for (auto& pid1_and_edges : edge_map[pid0]) {
            const PID pid1 = pid1_and_edges.first;
            auto& inc_polygons = pid1_and_edges.second;
            const std::size_t ninc = inc_polygons.size();
            CGAL_assertion(ninc <= 4);
            if (ninc <= 2) {
              continue;
            }

            auto get_third_point_id = [&polygons, pid0, pid1](const FID fid) -> PID
            {
              PID third;
              const auto& poly = polygons[fid];
              if (poly[0] == pid0 || poly[0] == pid1) {
                if (poly[1] == pid0 || poly[1] == pid1) {
                  third = poly[2];
                } else {
                  third = poly[1];
                }
              } else {
                third = poly[0];
              }

              CGAL_postcondition(third != pid0 && third != pid1);
              return third;
            };

            if (ninc == 3) {
#if 0
              // Two of these polygons must be from the bounding box (no associated facet ptr).
              // The third is a real facet. We need to order them around the edge
              // such that the sequence is H1-A-H2 or H2-A-H1 (CCW around pid0).
              FID real_fid = -1;
              boost::container::small_vector<FID, 2> hfids;
              for (int k=0; k<3; ++k) {
                if (polygon_to_facet[inc_polygons[k]]) {
                  real_fid = inc_polygons[k];
                } else {
                  hfids.push_back(inc_polygons[k]);
                }
              }

              CGAL_postcondition(real_fid != -1);
              CGAL_postcondition(hfids.size() == 2);
              CGAL_postcondition(real_fid != hfids[0] && real_fid != hfids[1]);

              const Point_3& ref_pt = points[get_third_point_id(real_fid)];

              bool h0_before_h1 = pred::sorted_around_edge<GeomTraits>(
                                    points[pid0], points[pid1],
                                    ref_pt,
                                    points[get_third_point_id(hfids[0])],
                                    points[get_third_point_id(hfids[1])]);

              inc_polygons[0] = (h0_before_h1 ? hfids[1] : hfids[0]);
              inc_polygons[1] = real_fid;
              inc_polygons[2] = (h0_before_h1 ? hfids[0] : hfids[1]);
#else
              const Point_3& ref_pt = points[get_third_point_id(inc_polygons[0])];
              bool h1_before_h2 = pred::sorted_around_edge<GeomTraits>(
                                    points[pid0], points[pid1],
                                    ref_pt,
                                    points[get_third_point_id(inc_polygons[1])],
                                    points[get_third_point_id(inc_polygons[2])]);

              if (!h1_before_h2) {
                std::swap(inc_polygons[1], inc_polygons[2]);
              }
#endif

            } else {
              // in the case '4', we can apply an A-B-A-B pattern since we know the edge
              // is incident to exactly 2x2 faces issued from 2 planes.
              std::array<FacetSPtr,4> fptrs;
              for (std::size_t k=0; k<4; ++k) {
                fptrs[k] = polygon_to_facet[inc_polygons[k]];
                CGAL_assertion(fptrs[k] != FacetSPtr());
              }

              // Identify the two distinct facet pointers and group polygon ids by facet
              FacetSPtr f0 = fptrs[0];
              FacetSPtr f1 = nullptr;
              boost::container::small_vector<FID, 2> group0, group1;
              for (std::size_t k=0; k<4; ++k) {
                if (fptrs[k] == f0) {
                  group0.push_back(inc_polygons[k]);
                } else {
                  if (!f1) {
                    f1 = fptrs[k];
                  }
                  CGAL_assertion(fptrs[k] == f1);
                  group1.push_back(inc_polygons[k]);
                }
              }

              CGAL_assertion(group0.size() == 2 && group1.size() == 2);

              // Determine ordering around the edge so faces are CCW when looking from pid0.
              // We pick an arbitrary face from group0 (A1) as reference and decide which
              // of the two group1 polygons (B candidates) comes next in CCW order.
              const Point_3& ref_pt = points[get_third_point_id(group0[0])];

              bool b0_before_b1 = pred::sorted_around_edge<GeomTraits>(
                                    points[pid0], points[pid1],
                                    ref_pt,
                                    points[get_third_point_id(group1[0])],
                                    points[get_third_point_id(group1[1])]);

              inc_polygons[0] = group0[0];
              inc_polygons[1] = (b0_before_b1 ? group1[0] : group1[1]);
              inc_polygons[2] = group0[1];
              inc_polygons[3] = (b0_before_b1 ? group1[1] : group1[0]);
            }
          }
        }
      }
#else // CGAL_SS3_FILL_EDGE_MAP_V2
      auto fill_edge_map = [](const std::vector<Point_3>& points,
                              const std::vector<std::vector<PID> >& polygons,
                              std::vector<CGAL::unordered_flat_map<PID, boost::container::small_vector<FID, 4> > >& edge_map)
      {
        CGAL_precondition(edge_map.size() == points.size());

        // collect duplicated edges in polygons
        auto add_edge = [&](PID a, PID b, FID fid)
        {
          if (b < a)
            std::swap(a, b);

          CGAL_SS3_TRANSF_TRACE_V(64, "collect edge: " << points[a] << " -- " << points[b] << " from polygon " << fid);

          auto& m = edge_map[a];
          auto [it, _] = m.try_emplace(b);
          it->second.push_back(fid);
        };

        for (FID fid=0; fid<polygons.size(); ++fid) {
          const auto& poly = polygons[fid];
          const std::size_t s = poly.size();
          for (std::size_t j=0; j+1<s; ++j) {
            add_edge(poly[j], poly[j + 1], fid);
          }
          add_edge(poly[s-1], poly[0], fid);
        }

        for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
          for (auto& pid1_and_edges : edge_map[pid0]) {
            const PID pid1 = pid1_and_edges.first;
            CGAL_SS3_TRANSF_TRACE_V(64, "processing NM edge: " << points[pid0] << " -- " << points[pid1]);
            CGAL_assertion(pid0 != pid1);

            auto& inc_polygons = pid1_and_edges.second;
            CGAL_assertion(inc_polygons.size() > 1);

            if (inc_polygons.size() <= 2) { // edge is only incident to a single SS3 face
              continue;
            }

            // the orientation of the edge might not match the orientation of the face
            // note that we know that by construction, we have strictly convex faces
            auto get_third_point_id = [&polygons, pid0, pid1](const FID fid) -> PID
            {
              PID third;
              const auto& poly = polygons[fid];
              if (poly[0] == pid0 || poly[0] == pid1) {
                if (poly[1] == pid0 || poly[1] == pid1) {
                  third = poly[2];
                } else {
                  third = poly[1];
                }
              } else {
                third = poly[0];
              }

              CGAL_postcondition(third != pid0 && third != pid1);
              return third;
            };

            const Point_3& ref_pt = points[get_third_point_id(inc_polygons[0])];

            auto less = [&ref_pt, &points, pid0, pid1, get_third_point_id](const FID fid1, const FID fid2)
            {
              return pred::sorted_around_edge<GeomTraits>(points[pid0], points[pid1],
                                                          ref_pt,
                                                          points[get_third_point_id(fid1)],
                                                          points[get_third_point_id(fid2)]);
            };

            std::sort(inc_polygons.begin()+1, inc_polygons.end(), less);

            // std::cout << "Around edge [" << pid0 << " " << pid1 << "], faces are sorted: ";
            // for(FID fid : inc_polygons) {
            //   std::cout << " " << fid;
            // }
            // std::cout << std::endl;
          }
        }
      }; // lambda 'fill_edge_map'

      fill_edge_map(points, polygons, edge_map);
#endif // CGAL_SS3_FILL_EDGE_MAP_V2

      CGAL_SS3_TRANSF_TRACE_V(16, "edge map size: " << edge_map.size());

      // for (PID pid0=0; pid0<static_cast<PID>(edge_map.size()); ++pid0) {
      //   for (const auto& pid1_and_edges : edge_map[pid0]) {
      //     const auto& inc_polygons = pid1_and_edges.second;
      //     std::cout << pid0 << " " << pid1_and_edges.first;
      //     for (const auto& fid : inc_polygons)
      //       std::cout << " " << fid;
      //     std::cout << std::endl;
      //   }
      // }

      auto build_volume_CC = [](const FID seed_fid,
                                const VID CC_ID,
                                const bool start_from_inverted_face,
                                const std::vector<Point_3>& points,
                                const std::vector<std::vector<PID> >& polygons,
                                const auto& edge_map,
                                auto& volume_CCs,
                                auto& face_volume_IDs)
      {
        CGAL_SS3_TRANSF_TRACE_V(64, "Building volume #" << CC_ID << " from seed face " << seed_fid);

        volume_CCs.emplace_back();

        std::stack<std::pair<FID, bool> > to_visit;
        to_visit.emplace(seed_fid, start_from_inverted_face);

        while (!to_visit.empty())
        {
          FID current_fid;
          bool invert_face;
          std::tie(current_fid, invert_face) = to_visit.top();
          to_visit.pop();

          CGAL_SS3_TRANSF_TRACE_V(64, "At face " << current_fid);
          CGAL_SS3_TRANSF_TRACE_V(64, "  invert: " << invert_face);
          CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS: " << face_volume_IDs[current_fid][0] << " " << face_volume_IDs[current_fid][1]);

          std::size_t pos = invert_face ? 0 : 1;
          if (face_volume_IDs[current_fid][pos] == CC_ID) {
            // already visited this facet during the flooding of this volume's boundary
            continue;
          }

          CGAL_assertion(face_volume_IDs[current_fid][pos] == VID(-1)); // polygon should only be encountered once
          CGAL_warning(face_volume_IDs[current_fid][(pos+1)%2] != CC_ID); // Moebius shenanigans should be an instance of a bug

          volume_CCs.back().push_back(current_fid);

          // mark face as visited
          face_volume_IDs[current_fid][pos] = CC_ID;

          // flood through the edges
          const std::size_t s = polygons[current_fid].size();
          for (std::size_t j=0; j<s; ++j) {
            PID e_pid0 = polygons[current_fid][j];
            PID e_pid1 = polygons[current_fid][(j+1)%s];
            if (e_pid1 < e_pid0)
              std::swap(e_pid0, e_pid1);

            std::pair<PID, PID> e_pids(e_pid0, e_pid1);
            const auto& inc_polygons = edge_map.at(e_pids.first).at(e_pids.second);
            CGAL_assertion(!inc_polygons.empty());

            CGAL_SS3_TRANSF_TRACE_V(64, "  ~~ Crossing edge [" << e_pids.first << ", " << e_pids.second << "]");
            CGAL_SS3_TRANSF_TRACE_V(64, "    pos: " << points[e_pids.first] << " " << points[e_pids.second]);

            // The faces are ordered CCW while looking from pid0.
            // So the walking while looking from [j] depends on whether [j] is pid0 or not
            int iter_direction = (e_pids.first == polygons[current_fid][j]) ? 1 : -1;
            CGAL_SS3_TRANSF_TRACE_V(64, "    iter_direction = " << iter_direction);

            // and it also depends on whether we are walking above or below the face
            iter_direction *= invert_face ? 1 : -1;
            CGAL_SS3_TRANSF_TRACE_V(64, "    invert_face = " << invert_face);

            FID next_fid = current_fid;
            for (;;) {
              if (inc_polygons.size() == 1) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Warning: dangling polygon...");
                CGAL_SS3_TRANSF_TRACE_V(1, "    over the edge, the polygon is ITSELF " << current_fid);
                to_visit.emplace(current_fid, !invert_face);
                break;
              } else if (inc_polygons.size() == 2) {
                // we should only be there once, meaning if we do not ignore orientations,
                // then the faces MUST be compatible
                CGAL_assertion(next_fid == current_fid);

                next_fid = (inc_polygons[0] == current_fid) ? inc_polygons[1] : inc_polygons[0];
                CGAL_SS3_TRANSF_TRACE_V(64, "    over the edge, the polygon is TRIVIALLY " << next_fid);
                CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS " << face_volume_IDs[next_fid][0] << " " << face_volume_IDs[next_fid][1]);
                CGAL_assertion(next_fid != current_fid);
              } else {
                // tricky part, now
                auto fid_it = std::find(std::begin(inc_polygons), std::end(inc_polygons), next_fid /*updates on every iteration*/);
                CGAL_assertion(fid_it != inc_polygons.end());

                if (iter_direction == 1) { // CCW
                  CGAL_SS3_TRANSF_TRACE_V(64, "    CCW walk");
                  auto next_it = std::next(fid_it);
                  next_fid = (next_it == inc_polygons.end()) ? inc_polygons[0] : *next_it;
                } else { // CW
                  CGAL_SS3_TRANSF_TRACE_V(64, "    CW walk");
                  next_fid = (fid_it == inc_polygons.begin()) ? inc_polygons.back() : *(std::prev(fid_it));
                }

                CGAL_SS3_TRANSF_TRACE_V(64, "    over the edge, the polygon is " << next_fid);
                CGAL_SS3_TRANSF_TRACE_V(64, "  VIDS " << face_volume_IDs[next_fid][0] << " " << face_volume_IDs[next_fid][1]);
                CGAL_assertion(next_fid != current_fid);
              }

              // If the edge has the same direction in both faces (aka, the orientation changes),
              // then we have to flip the direction of turning around the edge)
              const auto j_it = std::find(std::begin(polygons[next_fid]),
                                          std::end(polygons[next_fid]),
                                          polygons[current_fid][j]);
              CGAL_assertion(j_it != std::end(polygons[next_fid]));
              const std::size_t pos = std::distance(std::begin(polygons[next_fid]), j_it);
              CGAL_assertion(polygons[next_fid][pos] == polygons[current_fid][j]);

              const std::size_t next_s = polygons[next_fid].size();
              const bool flip_side = (polygons[next_fid][(pos+1)%next_s] == polygons[current_fid][(j+1)%s]);
              CGAL_SS3_TRANSF_TRACE_V(64, "    flipping? " << flip_side << " (N: " << polygons[next_fid][(pos+1)%next_s] << " C: " << polygons[current_fid][(j+1)%s] << ")");

              CGAL_SS3_TRANSF_TRACE_V(64, "Final FID = " << next_fid);
              if (flip_side) {
                to_visit.emplace(next_fid, !invert_face);
              } else {
                to_visit.emplace(next_fid, invert_face);
              }

              break;
            }
          }
        }
      }; // lambda 'build_volume_CC'

      CGAL_SS3_TRANSF_TRACE_V(16, "building volumes...");

      // identify volumes in the arrangement, and tag faces of the volumes
      // that are incident to the base face(s)
      std::vector<std::vector<FID> > volume_CCs; // range of ranges (volumes) of polygon IDs
      std::vector<std::array<VID, 2> > face_volume_IDs(polygons.size(),
                                                       // [0] is down, [1] is up
                                                       std::array<VID, 2>{VID(-1), VID(-1)});

      VID vid = 0;
      for(std::size_t i=0; i<polygons.size(); ++i) {
        // do not start from bbox faces: we will visit them anyway
        if (!polygon_to_facet[i])
          continue;
        if (face_volume_IDs[i][0] == VID(-1))
          build_volume_CC(i, vid++, true /*up*/, points, polygons, edge_map, volume_CCs, face_volume_IDs);
        if (face_volume_IDs[i][1] == VID(-1))
          build_volume_CC(i, vid++, false /*down*/, points, polygons, edge_map, volume_CCs, face_volume_IDs);
      }

      CGAL_SS3_TRANSF_TRACE_V(16, volume_CCs.size() << " volume CCs");

      for (std::size_t i=0; i<volume_CCs.size(); ++i) {
        // build a mesh from the soup
        std::vector<Point_3> cc_points = points;
        std::vector<std::vector<PID> > cc_polygons;
        for (FID fid : volume_CCs[i]) {
          cc_polygons.push_back(polygons[fid]);
        }

#ifdef CGAL_SS3_DUMP_FILES
        std::ostringstream oss;
        oss << "results/volume_cc_" << i << ".off";
        CGAL::IO::write_OFF(oss.str(), cc_points, cc_polygons, CGAL::parameters::stream_precision(17));
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

      // Classify some trivial CCs:
      // - if every non-bbox face points outwards, the CC is necessarily in
      // - if every non-bbox face points inwards, the CC is necessarily out
      for(std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
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
          std::vector<std::vector<PID> > cc_polygons;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            CGAL_assertion(e.first != CC_in_out_flag::TBD);
            CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

            if (in_out_flags[i] != e.first)
              continue;

            for (FID fid : volume_CCs[i])
              cc_polygons.push_back(polygons[fid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_base.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_polygons, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_base_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_polygons, CGAL::parameters::stream_precision(17));
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
      CGAL::unordered_flat_map<PID, CGAL::unordered_flat_map<PID, std::vector<FID> > > symmetrical_border_edge_map;
      for (FID fid=0; fid<polygons.size(); ++fid) {
        const std::size_t s = polygons[fid].size();
        for (std::size_t j=0; j<s; ++j) {
          const PID pid0 = polygons[fid][j];
          const PID pid1 = polygons[fid][(j+1)%s];
          if (!is_boundary_point[pid0] || !is_boundary_point[pid1]) {
            continue;
          }
          // avoid edges across the box
          const Point_3 m = CGAL::midpoint(points[pid0], points[pid1]);
          if (m.x() == bbox.xmin() || m.x() == bbox.xmax() ||
              m.y() == bbox.ymin() || m.y() == bbox.ymax() ||
              m.z() == bbox.zmin() || m.z() == bbox.zmax()) {
            symmetrical_border_edge_map[pid0][pid1].push_back(fid);
            symmetrical_border_edge_map[pid1][pid0].push_back(fid);
          }
        }
      }

      // Helper to mark cell directly above and below a facet incident to an edge
      auto mark_cells = [&](const PID ipid0, const PID ipid1, const FacetSPtr& current_facet)
      {
        PID pid0 = ipid0;
        PID pid1 = ipid1;
        if (pid1 < pid0)
          std::swap(pid0, pid1);

        const auto edge_it = symmetrical_border_edge_map.find(pid0);
        CGAL_assertion(edge_it != symmetrical_border_edge_map.end());
        const auto pid1_it = edge_it->second.find(pid1);
        CGAL_assertion(pid1_it != edge_it->second.end());

        const auto& fids = pid1_it->second;
        for (FID fid : fids) {
          if (polygon_to_facet[fid] != current_facet) {
            continue;
          }

          in_out_flags[face_volume_IDs[fid][1]] = CC_in_out_flag::OUTSIDE;
          in_out_flags[face_volume_IDs[fid][0]] = CC_in_out_flag::INSIDE;

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

        for (const auto& [other_pid, fids] : edge_it->second) {
          CGAL_assertion(is_boundary_point[current_pid] && is_boundary_point[other_pid]);

          bool on_current_facet = false;
          for (FID fid : fids) {
            if (polygon_to_facet[fid] == current_facet) {
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

        // determine convexity of prev/current
        std::list<EdgeSPtr> common_edges = prev_facet->find_edges(current_facet);
        CGAL_assertion(!common_edges.empty());

        // can't trust the edge is_reflex() function at that point, because the facets
        // have been perturbed and the convexity might have changed
        auto is_reflex = [&](const EdgeSPtr& e) -> bool {
          CGAL_SS3_TRANSF_TRACE_V(64, "is_reflex = " << e->source()->point() << " " << e->target()->point());
          CGAL_assertion(e->source()->point() != e->target()->point());

          bool result = false;
          const FacetSPtr facet_l = e->get_facet_L();
          const FacetSPtr facet_r = e->get_facet_R();
          CGAL_SS3_DEBUG_SPTR(facet_l);
          CGAL_SS3_DEBUG_SPTR(facet_r);
          const Plane_3& plane_l = facet_l->get_plane();
          const Plane_3& plane_r = facet_r->get_plane();
          CGAL_SS3_TRANSF_TRACE_V(64, "planes = " << plane_l << " " << plane_r);
          std::optional<Line_3> oline = Kernel_wrapper::intersection(plane_l, plane_r);
          CGAL_assertion(bool(oline));
          Vector_3 dir = oline->to_vector();
          CGAL_SS3_TRANSF_TRACE_V(64, "dir = " << dir);
          CGAL_assertion(dir != CGAL::NULL_VECTOR);
          // possibly reorient 'dir' to align with the direction of the edge
          // note that the edge is stable since its vertices are stable
          if (dir * Vector_3(e->source()->point(), e->target()->point()) < 0) {
            dir = -dir;
          }
          CGAL_SS3_TRANSF_TRACE_V(64, "canonical dir = " << dir);
          const Point_3 p_src = oline->point();
          const Vector_3 normal_l = plane_l.orthogonal_vector();
          CGAL_assertion(normal_l != CGAL::NULL_VECTOR);
          Point_3 p = p_src + CGAL::cross_product(normal_l, dir);
          if (plane_r.oriented_side(p) == CGAL::ON_POSITIVE_SIDE) {
            result = true;
          }
          CGAL_SS3_TRANSF_TRACE_V(64, "result = " << result);
          return result;
        };

        // find the correct common edge
        EdgeSPtr common_edge;
        for (const EdgeSPtr& e : common_edges) {
          if (e->source() == vertex || e->target() == vertex) {
            common_edge = e;
            break;
          }
        }
        CGAL_SS3_DEBUG_SPTR(common_edge);

        bool is_convex = !(is_reflex(common_edge));
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
          for (const auto& [pid1, fids] : edge_map[pid0]) {
            // The edge should be incident to only these two planes.
            if (fids.size() != 4)
              continue;

            // The edge should have one point on the boundary, one point inside.
            // It cannot cross the inner box because another plane will intersect the line.
            if (is_boundary_point[pid0] == is_boundary_point[pid1]) {
              continue;
            }

            std::set<FacetSPtr> incident_faces;
            for (FID fid : fids)
              incident_faces.insert(polygon_to_facet[fid]);
            CGAL_assertion(incident_faces.size() == 2);

            if (!incident_faces.count(prev_facet) || !incident_faces.count(facet)) {
              continue;
            }

            bool flip = is_boundary_point[pid0];
            PID test_spid0 = flip ? pid1 : pid0;
            PID test_spid1 = flip ? pid0 : pid1;
            CGAL_assertion(is_boundary_point[test_spid1]);

            Vector_3 new_v { points[test_spid0], points[test_spid1] };

            if (CGAL::scalar_product(orig_v, new_v) >= 0) {
              CGAL_SS3_TRANSF_TRACE_V(64, "start is " << points[test_spid1]);
              CGAL_SS3_TRANSF_TRACE_V(64, "start is " << normalized_points[test_spid1] << " (normalized)");
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

        // Now walk the star link on arrangement edges

        PID prev_pid (-1);
        for(;;) {
          CGAL_assertion(facet->get_plane().has_on(points[current_pid]));

          CGAL_SS3_TRANSF_TRACE_V(64, "at " << points[current_pid]);
          CGAL_SS3_TRANSF_TRACE_V(64, "at " << normalized_points[current_pid] << " (normalized)");

          // find the next point on 'facet' in the direction of 'next_facet'
          auto [next_pid, valid] = find_next_boundary_point(prev_pid, current_pid, prev_facet, facet, next_facet);
          CGAL_assertion(valid);
          CGAL_assertion(is_boundary_point[current_pid] && is_boundary_point[next_pid]);

#ifdef CGAL_SS3_DUMP_FILES
          link_out << "2 " << normalized_points[current_pid] << " " << normalized_points[next_pid] << "\n";
#endif
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
          std::vector<std::vector<PID> > cc_polygons;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            if (in_out_flags[i] != e.first)
              continue;

            for (FID fid : volume_CCs[i])
              cc_polygons.push_back(polygons[fid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_intermediate.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_polygons, CGAL::parameters::stream_precision(17));

            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_intermediate_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_polygons, CGAL::parameters::stream_precision(17));
          }
        }

        // dump boundary facets
        {
          CGAL::unordered_flat_map<std::size_t, std::vector<std::vector<std::size_t> > > id_to_polygons;

          for (std::size_t i=0; i<polygons.size(); ++i) {
            FacetSPtr fsptr = polygon_to_facet[i];
            if (!fsptr)
              continue;

            VID bot_vid = face_volume_IDs[i][0];
            VID top_vid = face_volume_IDs[i][1];
            CGAL_assertion(bot_vid != VID(-1) && top_vid != VID(-1));

            if (in_out_flags[top_vid] == CC_in_out_flag::OUTSIDE && in_out_flags[bot_vid] == CC_in_out_flag::INSIDE) {
              CGAL_SS3_TRANSF_TRACE_V(64, "Boundary facet " << i << " with top/bottom volumes " << top_vid << "/" << bot_vid);
              id_to_polygons[fsptr->id()].push_back(polygons[i]);
            }
          }

          for (const auto& kv : id_to_polygons) {
            std::ostringstream oss;
            oss << "results/intermediate_boundary_facet_" << kv.first << ".off";
            CGAL::IO::write_OFF(oss.str(), points, kv.second, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/intermediate_boundary_facet_" << kv.first << "_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, kv.second, CGAL::parameters::stream_precision(17));
          }

          CGAL_assertion(id_to_polygons.size() == vertex->degree());
        }
      }
#endif

      // Sanity checks:
      // - a boundary facet cannot have an OUTSIDE volume on its bottom and an INSIDE volume on its top
      for (std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
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

      CGAL_SS3_TRANSF_TRACE_V(16, "Cell unknowns: " << x_unknowns << " / " << C);

      // =====================================================================
      //  CONSTRAINT — Boundary facets point outside
      // =====================================================================

      // The boundary must be well oriented, hence a boundary facet must have the 'INSIDE' marker
      // on its negative side, and the 'OUTSIDE' marker on its positive side.
      // In model terms, this is:
      //   x[pos] => x[neg]. (if above is true, below is true)
      for(std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        model.AddImplication(x[face_volume_IDs[i][1]], x[face_volume_IDs[i][0]]);
      }

      // @speed could use hints using intersection of volumes?

      // Variables for the facets
      std::vector<BoolVar> b(polygons.size());
      unsigned int b_unknowns = 0;
      for (std::size_t fid=0; fid<polygons.size(); ++fid) {
        b[fid] = model.NewBoolVar();
        const VID bot_vid = face_volume_IDs[fid][0];
        const VID top_vid = face_volume_IDs[fid][1];

        if (bot_vid == VID(-1) || top_vid == VID(-1)) {
          model.FixVariable(b[fid], false);
        } else if (in_out_flags[bot_vid] == CC_in_out_flag::INSIDE && in_out_flags[top_vid] == CC_in_out_flag::OUTSIDE) {
          model.FixVariable(b[fid], true);
        } else { // @todo we can fix the facet that have incident fixed cells
          // A facet is on the boundary if and only if its two incident volumes are different
          model.AddNotEqual(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[fid]);
          model.AddEquality(x[bot_vid], x[top_vid]).OnlyEnforceIf(b[fid].Not());
          ++b_unknowns;
        }
      }

      CGAL_SS3_TRANSF_TRACE_V(16, "Face unknowns: " << b_unknowns << " / " << polygons.size());

      // =====================================================================
      //  CONSTRAINT — Global Euler characteristic (V_b − E_b + F_b = 1)
      // =====================================================================

      CGAL::unordered_flat_map<PID, std::unordered_set<FID> > vertex_incident_facets;
      for (FID fid=0; fid<polygons.size(); ++fid) {
        const std::size_t s = polygons[fid].size();
        for (int j=0; j<s; ++j) {
          vertex_incident_facets[polygons[fid][j]].insert(fid);
        }
      }

      // F_b
      LinearExpr F_b = LinearExpr::Sum(b);

      // E_b : edge on surface iff >= 1 incident boundary facet
      std::vector<BoolVar> e_b;
      e_b.reserve(edge_map.size());

      for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
        for (auto& pid1_and_edges : edge_map[pid0]) {
          e_b.push_back(model.NewBoolVar());
          auto& inc_polygons = pid1_and_edges.second;
          CGAL_assertion(!inc_polygons.empty());
          std::vector<IntVar> inc;
          inc.reserve(inc_polygons.size());
          for (FID fid : inc_polygons)
            inc.emplace_back(b[fid]);
          model.AddMaxEquality(e_b.back(), inc);
        }
      }
      LinearExpr E_b = LinearExpr::Sum(e_b);

      // V_b : vertex on surface iff >= 1 incident boundary facet
      std::vector<BoolVar> v_b(points.size());
      for (int v=0; v<points.size(); ++v) {
        v_b[v] = model.NewBoolVar();
        std::unordered_set<FID>& inc_polygons = vertex_incident_facets[v];
        CGAL_assertion(!inc_polygons.empty());
        std::vector<IntVar> inc;
        inc.reserve(inc_polygons.size());
        for (FID fid : vertex_incident_facets[v])
          inc.emplace_back(b[fid]);
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
          for (FID fid=0; fid<polygons.size(); ++fid) {
            if (polygon_to_facet[fid] == f)
              b_c.push_back(b[fid]);
          }
          CGAL_assertion(!b_c.empty());

          LinearExpr F_c = LinearExpr::Sum(b_c);

          // E_c
          std::vector<BoolVar> e_c;
          e_c.reserve(std::pow(vertex->degree(), 2));

          for (std::size_t pid0=0; pid0<points.size(); ++pid0) {
            for (auto& pid1_and_edges : edge_map[pid0]) {
              auto& inc_polygons = pid1_and_edges.second;
              std::vector<IntVar> inc;
              for (FID fid : inc_polygons) {
                if (polygon_to_facet[fid] == f)
                  inc.emplace_back(b[fid]);
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
            std::unordered_set<FID>& inc_polygons = vertex_incident_facets[v];
            for (FID fid : inc_polygons) {
              if (polygon_to_facet[fid] == f)
                inc.emplace_back(b[fid]);
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
        for (const auto& pid1_and_edges : edge_map[pid0]) {
          LinearExpr inc_n = 0;
          const auto& inc_polygons = pid1_and_edges.second;
          for (FID fid : inc_polygons)
            inc_n += b[fid];
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
            for (const auto& pid1_and_edges : edge_map[pid0]) {
              LinearExpr inc_n = 0;
              const auto& inc_polygons = pid1_and_edges.second;
              for (FID fid : inc_polygons) {
                if (polygon_to_facet[fid] == f)
                  inc_n += b[fid];
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

      // @fixme the non lazy constraints are more expensive and are probably wrong:
      // we could still pinch within the same CC even with a single root
#define CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
# ifdef CGAL_SLS3_USE_FLOW_BASED_NON_MANIFOLD_VERTEX_CONSTRAINTS
      for (FacetWPtr wf : vertex->facets()) {
        if (FacetSPtr f = wf.lock()) {
          std::vector<FID> fids_of_f;
          for (FID fid = 0; fid < polygons.size(); ++fid) {
            if (polygon_to_facet[fid] == f)
              fids_of_f.push_back(fid);
          }
          CGAL_assertion(!fids_of_f.empty());

          const int max_flow = static_cast<int>(fids_of_f.size());

          // set the fixed root
          FID root_fid = FID(-1);
          for (FID fid : fids_of_f) {
            const VID bot = face_volume_IDs[fid][0];
            const VID top = face_volume_IDs[fid][1];
            if (in_out_flags[bot] == CC_in_out_flag::INSIDE &&
                in_out_flags[top] == CC_in_out_flag::OUTSIDE) {
              root_fid = fid;
              break;
            }
          }
          CGAL_SS3_TRANSF_TRACE_V(32, "Root of F" << f->id() << " is " << root_fid);
          CGAL_assertion(root_fid != FID(-1));

          CGAL::unordered_flat_map<FID, std::vector<FID>> adj;
          for (std::size_t pid0 = 0; pid0 < points.size(); ++pid0) {
            for (const auto& pid1_and_edges : edge_map[pid0]) {
              const auto& inc = pid1_and_edges.second;
              for (std::size_t i = 0; i < inc.size(); ++i) {
                for (std::size_t j = i + 1; j < inc.size(); ++j) {
                  const FID t1 = inc[i], t2 = inc[j];
                  if (polygon_to_facet[t1] == f && polygon_to_facet[t2] == f) {
                    adj[t1].push_back(t2);
                    adj[t2].push_back(t1);
                  }
                }
              }
            }
          }

          // directed flow variables for each undirected edge
          CGAL::unordered_flat_map<FID, CGAL::unordered_flat_map<FID, IntVar>> flow;
          for (FID t1 : fids_of_f) {
            for (FID t2 : adj[t1]) {
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
          for (FID t : fids_of_f)
            b_f += b[t];

          for (FID t : fids_of_f) {
            LinearExpr inflow, outflow;
            for (FID t2 : adj[t]) {
              inflow  += flow[t2][t]; // flow coming from t2 into t
              outflow += flow[t][t2]; // flow going from t to t2
            }
            if (t == root_fid) {
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
          std::vector<FID> fids_of_f;
          for (FID fid = 0; fid < polygons.size(); ++fid) {
            if (polygon_to_facet[fid] == f)
              fids_of_f.push_back(fid);
          }
          CGAL_assertion(!fids_of_f.empty());
          const int N = fids_of_f.size();

          // Pick the fixed root
          FID root_fid = FID(-1);
          for (FID fid : fids_of_f) {
            const VID bot = face_volume_IDs[fid][0];
            const VID top = face_volume_IDs[fid][1];
            // If both cells are fixed and different => b[fid] is fixed true
            const bool bot_fixed = (in_out_flags[bot] == CC_in_out_flag::INSIDE);
            const bool top_fixed = (in_out_flags[top] == CC_in_out_flag::OUTSIDE);
            if (bot_fixed && top_fixed) {
              root_fid = fid;
              break;
            }
          }

          CGAL_assertion(root_fid != FID(-1));

          // Distance variables: d[fid] is the distance from the root.
          CGAL::unordered_flat_map<FID, IntVar> d;
          for (FID fid : fids_of_f) {
            d[fid] = model.NewIntVar({0, N}); // Max distance is the number of polygons
          }

          // The root is strictly at distance 0
          model.AddEquality(d[root_fid], 0);

          // Build local adjacency map for this specific color
          CGAL::unordered_flat_map<FID, std::vector<FID>> adj;
          for (std::size_t pid0 = 0; pid0 < points.size(); ++pid0) {
            for (const auto& pid1_and_edges : edge_map[pid0]) {
              const auto& inc = pid1_and_edges.second;
              for (std::size_t i = 0; i < inc.size(); ++i) {
                for (std::size_t j = i + 1; j < inc.size(); ++j) {
                  FID t1 = inc[i], t2 = inc[j];
                  if (polygon_to_facet[t1] == f && polygon_to_facet[t2] == f) {
                    adj[t1].push_back(t2);
                    adj[t2].push_back(t1);
                  }
                }
              }
            }
          }

          for (FID fid : fids_of_f) {
            if (fid == root_fid)
              continue;

            std::vector<BoolVar> parent_vars;

            for (FID nbr : adj[fid]) {
              BoolVar is_parent = model.NewBoolVar();
              parent_vars.push_back(is_parent);

              // If 'nbr' is the parent, 'nbr' MUST also be active on the boundary
              model.AddImplication(is_parent, b[nbr]);
              // If 'nbr' is the parent, distance strictly increases by 1
              model.AddEquality(d[fid], d[nbr] + 1).OnlyEnforceIf(is_parent);
            }

            // Connectivity Constraint:
            // If this facet is on the boundary (b[fid] == 1), it must have EXACTLY 1 parent.
            // If it is not on the boundary (b[fid] == 0), it has 0 parents.
            model.AddEquality(LinearExpr::Sum(parent_vars), b[fid]);
          }
        }
      }
# endif // CGAL_SLS3_USE_FLOW_BASED_NON_MANIFOLD_VERTEX_CONSTRAINTS
#endif

      // =====================================================================
      //  SOLVE
      // =====================================================================

      std::vector<CC_in_out_flag> solution(volume_CCs.size(), CC_in_out_flag::UNINITIALIZED);

#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
      const auto& proto = model.Build();
      CGAL_SS3_TRANSF_TRACE_V(16, "Variables: " << proto.variables_size());
      CGAL_SS3_TRANSF_TRACE_V(16, "Constraints: " << proto.constraints_size());

      SatParameters params;
      params.set_stop_after_first_solution(true);
      params.set_enumerate_all_solutions(false);

      std::vector<bool> current_b_values(polygons.size());

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
      // Since vertex manifoldness is difficult to express with constraints, and we know a solution exists
      // try and try till we succeed
      int iteration = -1;
      const int max_iterations = 1000; // @fixme could probably calculate that with Catalan numbers

      while (iteration < max_iterations) {
        ++iteration;
        CGAL_SS3_TRANSF_TRACE_V(32, "--- Starting solver iteration " << iteration << " ---");

        const auto& proto = model.Build();
        CGAL_SS3_TRANSF_TRACE_V(16, "Variables: " << proto.variables_size());
        CGAL_SS3_TRANSF_TRACE_V(16, "Constraints: " << proto.constraints_size());

        SatParameters params;
        params.set_stop_after_first_solution(true);
        params.set_enumerate_all_solutions(false);

        std::vector<bool> current_b_values(polygons.size());

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
            std::vector<std::vector<PID> > cc_polygons;

            for (std::size_t i=0; i<volume_CCs.size(); ++i) {
              CGAL_assertion(e.first != CC_in_out_flag::TBD);
              CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

              if (solution[i] != e.first)
                continue;

              for (FID fid : volume_CCs[i])
                cc_polygons.push_back(polygons[fid]);

              std::ostringstream oss;
              oss << "results/volumes_" << e.second << "_tentative.off";
              CGAL::IO::write_OFF(oss.str(), cc_points, cc_polygons, CGAL::parameters::stream_precision(17));
              std::ostringstream oss_n;
              oss_n << "results/volumes_" << e.second << "_tentative_normalized.off";
              CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_polygons, CGAL::parameters::stream_precision(17));
            }
          }
        }
#endif

        // Now, check for vertex manifoldness both in the global and local scale
        PID nm_vertex_id = -1;
        std::optional<FacetSPtr> failed_f = std::nullopt; // track which facet failed

        auto check_CC = [&](std::optional<FacetSPtr> of = std::nullopt) -> bool
        {
          std::vector<std::vector<PID> > local_polygons;
          for (std::size_t i=0; i<polygons.size(); ++i) {
            if (!polygon_to_facet[i])
              continue;

            if (solution[face_volume_IDs[i][0]] == solution[face_volume_IDs[i][1]])
              continue;

            if (of.has_value() && polygon_to_facet[i] != of.value())
              continue;

            local_polygons.push_back(polygons[i]);
          }

          CGAL_assertion(!local_polygons.empty());

          typedef std::vector<PID>                                            PointRange;
          typedef std::vector<std::vector<PID> >                              PolygonRange;
          typedef Polygon_mesh_processing::internal::Polygon_soup_orienter<PointRange, PolygonRange>   Orienter;

          typename Orienter::Edge_map edges(points.size());
          typename Orienter::Marked_edges marked_edges;
          Orienter::fill_edge_map(edges, marked_edges, local_polygons);
          CGAL_assertion(marked_edges.empty()); // no NM edges is part of the constraints

          Orienter::has_singular_vertices(points.size(), local_polygons, edges, marked_edges, nm_vertex_id);

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
          CGAL_SS3_TRANSF_TRACE_V(32, "Non-manifold vertex " << nm_vertex_id << " found");
          CGAL_SS3_TRANSF_TRACE_V(32, "  at position " << points[nm_vertex_id]);

          for (FID fid : vertex_incident_facets[nm_vertex_id]) {
            if (!polygon_to_facet[fid])
              continue;
            VID bot = face_volume_IDs[fid][0], top = face_volume_IDs[fid][1];
            CGAL_SS3_TRANSF_TRACE_V(64, "fid " << fid << " bot_flag=" << (int)in_out_flags[bot]
                                          << " top_flag=" << (int)in_out_flags[top]
                                          << " b=" << (solution[bot] != solution[top]));
          }

          std::vector<BoolVar> nogood_terms;
          bool some_not_fixed = false;

          for (FID fid : vertex_incident_facets[nm_vertex_id]) {
            if (!polygon_to_facet[fid])
              continue;

            if (failed_f.has_value() && polygon_to_facet[fid] != failed_f.value()) {
              continue;
            }

            VID bot = face_volume_IDs[fid][0];
            VID top = face_volume_IDs[fid][1];
            bool bot_fixed = (in_out_flags[bot] == CC_in_out_flag::INSIDE);
            bool top_fixed = (in_out_flags[top] == CC_in_out_flag::OUTSIDE);
            if (!bot_fixed || !top_fixed) {
              some_not_fixed = true; // At least one cell can be flipped!
            }

            // Reconstruct whether this facet was on the boundary in the rejected solution
            bool is_boundary = (solution[bot] != solution[top]);
            nogood_terms.push_back(is_boundary ? b[fid].Not() : b[fid]);
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
          break;
        }
      }

      CGAL_assertion(iteration < max_iterations);
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
          std::vector<std::vector<PID> > cc_polygons;

          for (std::size_t i=0; i<volume_CCs.size(); ++i) {
            CGAL_assertion(e.first != CC_in_out_flag::TBD);
            CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

            if (solution[i] != e.first)
              continue;

            for (FID fid : volume_CCs[i])
              cc_polygons.push_back(polygons[fid]);

            std::ostringstream oss;
            oss << "results/volumes_" << e.second << "_final.off";
            CGAL::IO::write_OFF(oss.str(), cc_points, cc_polygons, CGAL::parameters::stream_precision(17));
            std::ostringstream oss_n;
            oss_n << "results/volumes_" << e.second << "_final_normalized.off";
            CGAL::IO::write_OFF(oss_n.str(), normalized_points, cc_polygons, CGAL::parameters::stream_precision(17));
          }
        }
      }
#endif // CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

      // a valid partition has:
      // - C1) all cells are either inside or outside
      // - C2) boundary is manifold
      // - C3) only 1 "inside" face-connected component & 1 "outside" face-CC
      // - C4) only 1 simply connected component per input face
      auto is_valid_partition = [&](const std::vector<std::vector<FID> >& volume_CCs,
                                    const std::vector<std::array<VID, 2>>& face_volume_IDs,
                                    const std::vector<CC_in_out_flag>& in_out_flags,
                                    const std::vector<std::vector<PID> >& polygons,
                                    const std::vector<Point_3>& points) -> bool
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
          std::vector<std::vector<PID> > local_polygons;
          for (std::size_t i=0; i<polygons.size(); ++i) {
            if (!polygon_to_facet[i])
              continue;

            if (in_out_flags[face_volume_IDs[i][0]] == in_out_flags[face_volume_IDs[i][1]])
              continue;

            if (of.has_value() && polygon_to_facet[i] != of.value())
              continue;

            local_polygons.push_back(polygons[i]);
          }

          CGAL_assertion(!local_polygons.empty());

          if (!PMP::is_polygon_soup_a_polygon_mesh(local_polygons)) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: failed C34-PS_PM");
#ifdef CGAL_SS3_DUMP_FILES
            CGAL::IO::write_polygon_soup("results/bad_CC_normalized.off", normalized_points, local_polygons, CGAL::parameters::stream_precision(17));
#endif
            return false;
          }

          // @todo avoid using a Surface_mesh
          CGAL::Surface_mesh<Point_3> sm;
          PMP::polygon_soup_to_polygon_mesh(points, local_polygons, sm);
          CGAL_assertion(faces(sm).size() == local_polygons.size());

          const unsigned int nb = number_of_borders(sm);
          if (nb != 1) {
            CGAL_SS3_TRANSF_TRACE_V(1, "Error: failed C34-borders (" << nb << ")");
            return false;
          }

          // @tmp not triangulated anymore...
          // if (PMP::does_self_intersect(sm)) {
          //   CGAL_SS3_TRANSF_TRACE_V(1, "Error: CC self intersects");
          //   return false;
          // }

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

      CGAL_postcondition(is_valid_partition(volume_CCs, face_volume_IDs, solution, polygons, points));

      // Now, extract the connectivity from the arrangement and plug it back into the polyhedron itself
      std::vector<std::vector<PID> > boundary_polygons;
      std::vector<FacetSPtr> boundary_polygon_to_facet;
      for (std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
          continue;

        if (solution[face_volume_IDs[i][0]] == solution[face_volume_IDs[i][1]])
          continue;

        boundary_polygons.push_back(polygons[i]);
        boundary_polygon_to_facet.push_back(polygon_to_facet[i]);
      }

#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_polygon_soup("results/boundary.off", points, boundary_polygons);
#endif

      CGAL_assertion(!boundary_polygons.empty());
      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(boundary_polygons));

      // @todo avoid using a CGAL::Surface_mesh
      using Mesh = CGAL::Surface_mesh<Point_3>;
      using edge_descriptor = typename boost::graph_traits<Mesh>::edge_descriptor;
      using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
      using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

      Mesh bsm;
      auto ifpm = get(CGAL::dynamic_face_property_t<FacetSPtr>{}, bsm);

      auto out = boost::make_function_output_iterator(
        [ifpm, &boundary_polygon_to_facet](const std::pair<std::size_t, face_descriptor>& pid_f) {
          // we shouldn't have polygons that are not from a face here
          CGAL_assertion(boundary_polygon_to_facet[pid_f.first] != FacetSPtr());
          put(ifpm, pid_f.second, boundary_polygon_to_facet[pid_f.first]);
        });

      PMP::polygon_soup_to_polygon_mesh(points, boundary_polygons, bsm,
                                        CGAL::parameters::polygon_to_face_output_iterator(out));
      CGAL_assertion(faces(bsm).size() == boundary_polygons.size());

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
    IO::write_OBJ("results/V4_final.obj", polyhedron, parameters::do_not_triangulate_faces(true));
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

    ConfigurationSPtr config = Configuration::get_instance();
    const bool safe_mode = config->get_Boolean("Preprocessing", "check_degenerate_configuration");

    // if the input is all triangles, simply perturb points directly
#ifndef CGAL_SPS3_USE_V4_PERTURBATION
    // don't do this with V4 because we perturb and split at once
    if (is_triangle_polyhedron(polyhedron))
      return rand_move_points(polyhedron);
#endif

    // Generic approach
    Transformation::normalize_facet_planes(polyhedron); // @todo hasn't this already been done before?

    PolyhedronSPtr p_mem;
    if (safe_mode) {
      p_mem = polyhedron->clone();
    }

    CGAL::Real_timer timer;
    timer.start();

#ifdef CGAL_SPS3_USE_V4_PERTURBATION
    apply_rand_plane_tilts_V4(polyhedron);
#else
    apply_rand_plane_tilts_V3(polyhedron);
#endif

    timer.stop();
    std::cout << "perturbation time: " << timer.time() << std::endl;

    if (safe_mode) {
      CGAL_SS3_TRANSF_TRACE_V(8, "Safe mode is enabled, checking validity of the perturbation...");

      for (;;) {
#ifdef CGAL_SS3_DUMP_FILES
        IO::write_OBJ("results/last_perturbation.obj", polyhedron,
                      parameters::stream_precision(17).do_not_triangulate_faces(true));
#endif

        if (are_planes_in_general_position(polyhedron) &&
            !Self_intersection::has_self_intersecting_surface(polyhedron)) {
          CGAL_SS3_TRANSF_TRACE_V(8, "Safe mode is enabled, checking validity of the perturbation...");
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
