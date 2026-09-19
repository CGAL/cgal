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

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_checker.h"
#include "ortools/sat/sat_parameters.pb.h"

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
  static std::mt19937_64& gen()
  {
    // static std::random_device rd;
    constexpr auto s = 0x9E3779B97F4A7C15ull; // rd()
    static thread_local std::mt19937_64 gen(s);
    return gen;
  }

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
    std::uniform_real_distribution<> rdist(min, max);
    return { rdist(gen()), rdist(gen()), rdist(gen()) };
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
    CGAL_SS3_TRANSF_TRACE_V(16, "Check if planes are in general position...");
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
          CGAL_SS3_TRANSF_TRACE_V(1, "Degenerate facet pair: " << normals[i].first->id() << " " << normals[j].first->id() << " -1");
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

  static bool is_triangle_polyhedron(const PolyhedronSPtr& polyhedron)
  {
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->vertices().size() != 3) {
        return false;
      }
    }
    return true;
  }

  // @fixme this approach does not even yield a valid input as the planes
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
    CGAL_SS3_TRANSF_TRACE_V(32, "Perturb (Fixed) F" << facet->id() << " [size: " << facet->vertices().size() << "]");
    CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << facet->get_plane().a() << " " << facet->get_plane().b() << " "
                                                        << facet->get_plane().c() << " " << facet->get_plane().d() << "]");
    CGAL_SS3_TRANSF_TRACE_V(32, "  with " << fixed_vertices.size() << " fixed vertices");
    for (const VertexSPtr& fixed_vertex : fixed_vertices) {
      CGAL_SS3_TRANSF_TRACE_V(32, "    fixed V" << fixed_vertex->id() << " at " << fixed_vertex->point());
    }
    CGAL_SS3_TRANSF_TRACE_V(32, "  Nudge range = " << range);

    CGAL_precondition(Transformation::has_normalized_plane(facet));
    CGAL_precondition(range > 0.);

    boost::container::small_vector<const Point_3*, 3> fixed_points;
    for (const VertexSPtr& fixed_vertex : fixed_vertices) {
      fixed_points.push_back(&(fixed_vertex->point()));
    }

    const std::size_t fixed_points_count = fixed_points.size();
    CGAL_assertion(fixed_points_count < 3);

    // In [r/2, r] to avoid near-zero nudges, and cancellations.
    // Note that planes have normalized coefficients.
    auto signed_magnitude = [&](const double r) -> double
    {
      CGAL_precondition(r > 0.);
      std::uniform_real_distribution<> mdist(0.5 * r, r);
      const double m = mdist(gen());
      return std::bernoulli_distribution(0.5)(gen()) ? m : -m;
    };

    auto nudge = [&](const FT& v)
    {
      // Since we are perturbing, we might as well collapse the DAG of 'v'.
      // The point is also that once 'nv' is a double, its interval will be a singleton,
      // and we will have access to static filters
      const double step = signed_magnitude(range);
      const double nv = CGAL::to_double(v) + step;
      CGAL_assertion(nv != CGAL::to_double(v)); // 'range' too small for this magnitude
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

      // NOTE: 'phi' is typically orders of magnitude larger than 'range' (the anchors force
      // the normal to swing much more than we perturb it). Two facets that share both anchors
      // are therefore projected onto nearly the same element of the pencil, and only 'theta'
      // (plus the jitter below) separates them.

      // Rotate around 'u' by a random angle
      const double theta = signed_magnitude(range);
      const double ct = std::cos(theta), st = std::sin(theta);

      const double w[3] = { uh[1]*tp[2] - uh[2]*tp[1],
                            uh[2]*tp[0] - uh[0]*tp[2],
                            uh[0]*tp[1] - uh[1]*tp[0] };

      const double tj[3] = { ct*tp[0] + st*w[0],
                             ct*tp[1] + st*w[1],
                             ct*tp[2] + st*w[2] };

      // Express 'tj' in the (unit) basis (e1,e2) of u^perp.
      // Caveat: when tj is close to e1 (the common case, since the pencil basis is built
      // from u and the normal is nearly orthogonal to u), 'r2 - g*r1' cancels almost
      // completely and 'beta' retains only a handful of significant bits.
      // That is one of the two reasons for the jitter below.
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

      // Exact sub-ulp jitter.
      //
      // The reachable set of 'beta' across the whole 'theta' range is only a few dozen
      // multiples of the 2^-shift quantum (and 'beta' itself carries only ~8 significant
      // bits because of the cancellation noted above), so two facets sharing the same two
      // anchors can round to the very same (lam_d, mu_d) and end up with planes.
      // Adding a random exact multiple of 2^-(shift + jitter) breaks that tie.
      constexpr int jitter = 24;
      CGAL_assertion(shift + jitter < 1000); // jitter step must stay normal
      const double jitter_step = std::ldexp(1.0, -(shift + jitter));
      std::uniform_int_distribution<int> jdist(-(1 << (jitter - 1)), (1 << (jitter - 1)));
      const int jitter_lam = jdist(gen());
      const int jitter_mu  = jdist(gen());

      const FT alpha_ft = FT(lam_d) * inv_scl_ft + FT(jitter_lam) * FT(jitter_step);
      const FT beta_ft  = FT(mu_d)  * inv_scl_ft + FT(jitter_mu)  * FT(jitter_step);

      CGAL_SS3_TRANSF_TRACE_V(64, "  pencil: phi = " << phi << ", theta = " << theta
                                  << ", lam = " << lam_d << " (+" << jitter_lam << "/2^" << jitter << ")"
                                  << ", mu = "  << mu_d  << " (+" << jitter_mu  << "/2^" << jitter << ")");

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
      const double epsilon = signed_magnitude(range);
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

  static void apply_plane_perturbation_V4(const PolyhedronSPtr& polyhedron)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;
    namespace pred = PMP::Corefinement;

    CGAL_SS3_TRANSF_TRACE_V(4, "Plane Perturbation (v4)");
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
    CGAL_SS3_TRANSF_TRACE_V(8, "Input nudge range: " << nudge_range);

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
    for (const FacetSPtr& f : polyhedron->facets()) {
      CGAL_SS3_TRANSF_TRACE_V(32, "Original plane of F" << f->id() << " is " << f->get_plane());
      original_planes[f] = f->get_plane();
    }

    double characteristic_length = 1.0; // @fixme take bbox's diagonal?
    double relative_tolerance = 1e-7; // @fixme hardcoded
    CGAL::unordered_flat_map<VertexSPtr, std::array<double, 3> > vertex_tolerances;

    // @fixme should also be bounded by the LFS
    for (const VertexSPtr& v : polyhedron->vertices()) {
      auto axis_tol = [&](const FT& c) {
        const double s = std::max(std::abs(CGAL::to_double(c)), characteristic_length);
        return relative_tolerance * s;
      };
      vertex_tolerances[v] = { axis_tol(v->point().x()),
                               axis_tol(v->point().y()),
                               axis_tol(v->point().z()) };
    }

    auto is_stable = [&](const VertexSPtr& v) -> bool
    {
      const Point_3& orig = original_points[v];
      const std::array<double, 3>& tol = vertex_tolerances[v];

      auto within_tolerance_box = [&](const Point_3& p) -> bool {
        return CGAL::abs(p.x() - orig.x()) <= tol[0] &&
               CGAL::abs(p.y() - orig.y()) <= tol[1] &&
               CGAL::abs(p.z() - orig.z()) <= tol[2];
      };

      for (auto it_wf1 = v->facets().begin(); it_wf1 != v->facets().end(); ++it_wf1) {
        if (FacetSPtr f1 = it_wf1->lock()) {
          for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
            if (FacetSPtr f2 = it_wf2->lock()) {
              for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                if (FacetSPtr f3 = it_wf3->lock()) {
                  std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                  if (!p_new.has_value()) {
                    CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point! (a)");
                    CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                    continue;
                  }

                  if (!within_tolerance_box(*p_new)) {
                    CGAL_SS3_TRANSF_TRACE_V(32, "  TOO FAR: " << v->to_string());
                    CGAL_SS3_TRANSF_TRACE_V(32, "  from " << orig << " to " << p_new.value());
                    CGAL_SS3_TRANSF_TRACE_V(32, "  |dx|=" << CGAL::abs(p_new->x() - orig.x()) << " VS tol_x=" << tol[0]);
                    CGAL_SS3_TRANSF_TRACE_V(32, "  |dy|=" << CGAL::abs(p_new->y() - orig.y()) << " VS tol_y=" << tol[1]);
                    CGAL_SS3_TRANSF_TRACE_V(32, "  |dz|=" << CGAL::abs(p_new->z() - orig.z()) << " VS tol_z=" << tol[2]);
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f1->id() << " [" << f1->get_plane().a() << " " << f1->get_plane().b() << " "
                                                                          << f1->get_plane().c() << " " << f1->get_plane().d() << "] || from original: "
                                                                          << CGAL::squared_distance(f1->get_plane(), orig) << " || from new: "
                                                                          << CGAL::squared_distance(f1->get_plane(), p_new.value()));
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f2->id() << " [" << f2->get_plane().a() << " " << f2->get_plane().b() << " "
                                                                          << f2->get_plane().c() << " " << f2->get_plane().d() << "] || from original: "
                                                                          << CGAL::squared_distance(f2->get_plane(), orig) << " || from new: "
                                                                          << CGAL::squared_distance(f2->get_plane(), p_new.value()));
                    CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f3->id() << " [" << f3->get_plane().a() << " " << f3->get_plane().b() << " "
                                                                          << f3->get_plane().c() << " " << f3->get_plane().d() << "] || from original: "
                                                                          << CGAL::squared_distance(f3->get_plane(), orig) << " || from new: "
                                                                          << CGAL::squared_distance(f3->get_plane(), p_new.value()));
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

#if 0
    CGAL_SS3_TRANSF_TRACE_V(8, "BEFORE");

    // @todo if we have general position + stability + no high degree, we could skip perturbation
    CGAL_SS3_TRANSF_TRACE_CODE(std::vector<Stability_failure> failures;)
    CGAL_SS3_TRANSF_TRACE_V(32, "Unperturbed is in general position?\n" << are_planes_in_general_position(polyhedron, &failures));

#ifdef CGAL_SS3_DUMP_FILES
    std::ofstream out_unstable_before("results/unstable_vertices_before.xyz");
    out_unstable_before.precision(17);
#endif

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE_CODE(if (!is_stable(v)) {);
#ifdef CGAL_SS3_DUMP_FILES
    CGAL_SS3_TRANSF_TRACE_CODE(out_unstable_before << v->point() << "\n");
#endif
    CGAL_SS3_TRANSF_TRACE_CODE(})

#ifdef CGAL_SS3_DUMP_FILES
    out_unstable_before.close();
#endif

#endif

    CGAL_SS3_TRANSF_TRACE_V(8, "START");

    {
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
              const Plane_3 original_plane = original_planes.at(fprime); // intentional copy
              original_planes.erase(fprime);

              auto [_, new_facets] = Transformation::triangulate_facet(fprime, polyhedron);
              for (const FacetSPtr& nf : new_facets) {
                original_planes[nf] = original_plane;
              }

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

      auto nudge_point = [&](const VertexSPtr& v,
                             const double lo_frac = 1e-5,
                             const double hi_frac = 1e-4) -> Point_3
      {
        auto nudge_axis = [&](const double c, const double tol) {
          std::uniform_real_distribution<double> mag(lo_frac * tol, hi_frac * tol);
          std::bernoulli_distribution sign(0.5);
          const double eps = mag(gen()) * (sign(gen()) ? 1.0 : -1.0);
          double out = c + eps;
          if (out == c) {
            out = std::nextafter(c, eps >= 0.0 ? +HUGE_VAL : -HUGE_VAL);
        }
          return out;
        };

        const std::array<double, 3>& tolerances = vertex_tolerances.at(v);
        return { nudge_axis(CGAL::to_double(v->point().x()), tolerances[0]),
                 nudge_axis(CGAL::to_double(v->point().y()), tolerances[1]),
                 nudge_axis(CGAL::to_double(v->point().z()), tolerances[2]) };
      };

#define CGAL_SS3_PERTURB_V4_ANCHORING_V2

#ifdef CGAL_SS3_PERTURB_V4_ANCHORING_V1
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
        // Worst, over the incident plane triplets, of max_i |d_i| / tol_i. Dimensionless, hence
        // comparable between vertices with different tolerance boxes; > 1 iff unstable.
        // Stays 0 for a stable vertex.
        double max_violation = 0.;
        FacetSPtr worst_facet;
      };

      struct Vertex_stability_record
      {
        VertexSPtr vertex;
        FacetSPtr facet;
        double max_violation;
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
          return (*records)[lhs].max_violation > (*records)[rhs].max_violation;
        }

        const std::vector<Vertex_stability_record>* records;
      };

      auto evaluate_vertex_stability = [&](const VertexSPtr& v) -> Vertex_stability_info
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Check stability of V" << v->id() << " (deg: " << v->degree() << ")");

        const std::array<double, 3>& tol = vertex_tolerances.at(v);
        CGAL_assertion(tol[0] > 0 && tol[1] > 0 && tol[2] > 0);

        // Reference position: the vertex's current position, i.e. its input position until it
        // gets anchored, and its anchored position afterwards (see make_random_anchor()).
        const Point_3& p = original_points.at(v);

        Vertex_stability_info info;

        for (auto it_wf1 = v->facets().begin(); it_wf1 != v->facets().end(); ++it_wf1) {
          if (FacetSPtr f1 = it_wf1->lock()) {
            for (auto it_wf2 = std::next(it_wf1); it_wf2 != v->facets().end(); ++it_wf2) {
              if (FacetSPtr f2 = it_wf2->lock()) {
                for (auto it_wf3 = std::next(it_wf2); it_wf3 != v->facets().end(); ++it_wf3) {
                  if (FacetSPtr f3 = it_wf3->lock()) {
                    std::optional<Point_3> p_new = Kernel_wrapper::intersection(f1->get_plane(), f2->get_plane(), f3->get_plane());
                    if (!p_new.has_value()) {
                      CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point! (b)");
                      CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                      continue;
                    }

                    const double dx = CGAL::to_double(CGAL::abs(p_new->x() - p.x()));
                    const double dy = CGAL::to_double(CGAL::abs(p_new->y() - p.y()));
                    const double dz = CGAL::to_double(CGAL::abs(p_new->z() - p.z()));
                    const bool violates = (dx > tol[0] || dy > tol[1] || dz > tol[2]);

                    // to order the PQ
                    const double violation = (std::max)({ dx / tol[0], dy / tol[1], dz / tol[2] });

                    if (violates && (!info.worst_facet || violation > info.max_violation)) {
                      auto dump_facet = [&](const FacetSPtr& f)
                      {
                        std::stringstream oss;
                        oss << "F" << f->id() << " " << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                     << f->get_plane().c() << " " << f->get_plane().d() << " ";
                        oss << "[dist = " << CGAL::approximate_sqrt(CGAL::squared_distance(f->get_plane(), p)) << "] ";
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
                      CGAL_SS3_TRANSF_TRACE_V(64, "  from " << p << " to " << p_new.value());
                      CGAL_SS3_TRANSF_TRACE_V(64, "  |dx| " << dx << " VS " << tol[0]);
                      CGAL_SS3_TRANSF_TRACE_V(64, "  |dy| " << dy << " VS " << tol[1]);
                      CGAL_SS3_TRANSF_TRACE_V(64, "  |dz| " << dz << " VS " << tol[2]);
                      CGAL_SS3_TRANSF_TRACE_V(64, "  violation factor " << violation);
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f1));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f2));
                      CGAL_SS3_TRANSF_TRACE_V(64, "  " << dump_facet(f3));

                      info.is_stable = false;
                      info.max_violation = violation;

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

        CGAL_assertion(info.is_stable != bool(info.worst_facet));

        if (info.is_stable) {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is stable");
        } else {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is unstable, exceeding its tolerance box by a factor "
              << info.max_violation << " (tolerances " << tol[0] << " " << tol[1] << " " << tol[2] << ")");
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

          // We don't know about the order of the fixed points, so trust the input normal for orientation
          Plane_3 new_pl(p0, p1, p2);
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

      CGAL_assertion(are_planes_in_general_position(polyhedron));

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
        const Vertex_stability_info info = evaluate_vertex_stability(v);
        if (!info.is_stable) {
#ifdef CGAL_SS3_DUMP_FILES
          out_unstable_base << v->point() << "\n";
#endif
          unstable_record_indices[v] = unstable_vertices.size();
          unstable_vertices.push_back({v, info.worst_facet, info.max_violation, true, 0});
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
        const Vertex_stability_info info = evaluate_vertex_stability(v);

        const auto it = unstable_record_indices.find(v);
        if (it != unstable_record_indices.end()) { // already exists
          CGAL_SS3_TRANSF_TRACE_V(64, "record exists");
          Vertex_stability_record& record = unstable_vertices[it->second];
          record.vertex = v;
          record.vertex_id = v->id();
          record.facet = info.worst_facet;
          record.max_violation = info.max_violation;
          record.is_active = !info.is_stable;
          CGAL_SS3_TRANSF_TRACE_V(64, "activity: " << record.is_active);
          if (!record.is_active) {
            CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is now stable");
            record.facet = nullptr;
            record.max_violation = 0.;
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
            unstable_vertices.push_back({v, info.worst_facet, info.max_violation, true, 0});
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
          const Point_3 new_p = nudge_point(v);

          CGAL_SS3_TRANSF_TRACE_V(8, "New anchor V" << v->id() << " at " << new_p);
          CGAL_SS3_TRANSF_TRACE_V(8, "Original position: " << v->point());
          CGAL_SS3_TRANSF_TRACE_V(8, "  Distance from original: " << CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), new_p)));

          CGAL_assertion(new_p != v->point());

          // The anchor nudge must fit inside the vertex's own tolerance box, otherwise
          // anchoring alone would keep the vertex unstable and the loop would not converge.
          CGAL_assertion(CGAL::abs(new_p.x() - v->point().x()) <= vertex_tolerances.at(v)[0] &&
                         CGAL::abs(new_p.y() - v->point().y()) <= vertex_tolerances.at(v)[1] &&
                         CGAL::abs(new_p.z() - v->point().z()) <= vertex_tolerances.at(v)[2]);

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
#endif // CGAL_SS3_PERTURB_V4_ANCHORING_V1

#ifdef CGAL_SS3_PERTURB_V4_ANCHORING_V2
      // =====================================================================
      //  Perturbation to general position, with bounded vertex displacement
      // =====================================================================
      //
      // Phase 0: perturb every facet plane.
      // Phase 1: greedily freeze the large polygons, then move every vertex to its anchor
      //          position -- a RANDOM point of the intersection of its frozen facets, inside
      //          its tolerance box. A facet is frozen only if that intersection meets the box
      //          for each of its vertices, so afterwards:
      //
      //            (i)  every vertex is inside its tolerance box, and
      //            (ii) every frozen facet's plane passes exactly through all of its vertices.
      //
      // Phase 2: pop the most unstable vertex, anchor it in one of the facets responsible for
      //          the violation, recompute that facet's plane through its anchors, and update
      //          the incident vertices. Repeat until every vertex is stable.
      //
      // A vertex is stable iff every triplet of its incident planes meets inside its tolerance
      // box, measured from original_points[v] -- the single reference used throughout.
      //
      // THE INPUT IS NOT ASSUMED STABLE. If three normals at a vertex are exactly coplanar in
      // the input, perturbing by eps gives a determinant of O(eps) and a triple point displaced
      // by O(eps/eps) = O(1): shrinking the perturbation does not help, so there is no "shrink
      // and retry" anywhere here. Anchoring is the repair, and it is unconditional -- three
      // planes through one common point meet AT that point, exactly, however ill-conditioned
      // they are. Freezing is gated on the same criterion, so a degenerate configuration simply
      // fails to freeze and is left to the anchoring loop.
      //
      // Progress: a frozen facet, and a facet anchored at v, both pass through v->point(). If
      // all three facets of a triplet were such, their common point would be v->point(), which
      // is inside the box by (i), hence not a violation. So a violating triplet always contains
      // a facet that is neither frozen nor anchored at v, i.e. a legal candidate.
      //
      // Termination: each pop consumes one (vertex, facet) incidence for good, and a facet is
      // triangulated when it reaches a 4th anchor -- which a triangle never can -- so each
      // original polygon is triangulated at most once.
      //
      // WHY THE ANCHOR POSITIONS MUST BE RANDOM. General position is a property of ALL facet
      // triplets, including facets that share no vertex and that the stability loop therefore
      // never examines; only the randomization protects those. A facet with 3 anchors has its
      // plane *exactly* determined by them, with no perturbation of its own, so the anchors are
      // its only source of randomness. Placing them deterministically (say, at the closest
      // point of the frozen structure) reproduces the input's own degeneracies: on an
      // axis-aligned model, two opposite faces plus any third face have coplanar normals, and
      // a 3-anchor facet whose anchors sit at their original positions re-creates exactly that.
      //
      // Note: a vertex anchored in more than 3 facets makes those planes concurrent, i.e. it
      // stays a high-degree vertex. That is resolved later by translating planes, which is safe
      // precisely because stability bounds where the resulting degree-3 points can land.

      CGAL::unordered_flat_map<FacetSPtr, boost::container::small_vector<VertexSPtr, 3> > anchors;
      CGAL::unordered_flat_set<FacetSPtr> frozen;
      CGAL::unordered_flat_map<VertexSPtr, boost::container::small_vector<FacetSPtr, 3> > frozen_at;
      std::unordered_set<VertexSPtr> all_anchors;

      auto in_tolerance_box = [&](const VertexSPtr& v, const Point_3& q) -> bool
      {
        const Point_3& p = original_points.at(v);
        const std::array<double, 3>& t = vertex_tolerances.at(v);
        return CGAL::abs(q.x() - p.x()) <= FT(t[0]) &&
               CGAL::abs(q.y() - p.y()) <= FT(t[1]) &&
               CGAL::abs(q.z() - p.z()) <= FT(t[2]);
      };

      // Projection of 't' onto the intersection of the frozen facets at 'v' (plus 'cand', if
      // given). Empty if that intersection is empty, or if there are more than 3 planes --
      // generically they then have no common point at all.
      auto project_on_frozen = [&](const VertexSPtr& v,
                                   const Point_3& t,
                                   const FacetSPtr& cand = FacetSPtr()) -> std::optional<Point_3>
      {
        boost::container::small_vector<const Plane_3*, 4> pl;
        for (const FacetSPtr& g : frozen_at[v]) {
          pl.push_back(&(g->get_plane()));
        }
        if (cand) {
          pl.push_back(&(cand->get_plane()));
        }

        if (pl.size() > 3) {
          return std::nullopt;
        }
        if (pl.empty()) {
          return t;
        }
        if (pl.size() == 1) {
          return pl[0]->projection(t);
        }
        if (pl.size() == 2) {
          std::optional<Line_3> L = Kernel_wrapper::intersection(*pl[0], *pl[1]);
          if (!L.has_value()) {
            return std::nullopt;
          }
          return L->projection(t);
        }
        return Kernel_wrapper::intersection(*pl[0], *pl[1], *pl[2]);
      };

      auto random_point_in_box = [&](const VertexSPtr& v, const double shrink) -> Point_3
      {
        const Point_3& p = original_points.at(v);
        const std::array<double, 3>& tol = vertex_tolerances.at(v);

        // Uniform offset in [-amplitude, amplitude], or exactly 0 when the amplitude is too
        // small for the addition to move the coordinate at all.
        auto random_offset = [&](const double c, const double amplitude) -> double
        {
          constexpr double min_ulps = 8.0; // a few representable values on each side of 0
          if (amplitude < min_ulps * (std::max)(std::abs(c), 1.0) * std::numeric_limits<double>::epsilon()) {
            return 0.;
          }

          std::uniform_real_distribution<> dist(-amplitude, amplitude);
          return dist(gen());
        };

        const double dx = random_offset(CGAL::to_double(p.x()), shrink * tol[0]);
        const double dy = random_offset(CGAL::to_double(p.y()), shrink * tol[1]);
        const double dz = random_offset(CGAL::to_double(p.z()), shrink * tol[2]);

        if (dx == 0. && dy == 0. && dz == 0.) {
          CGAL_SS3_TRANSF_TRACE_V(32, "    V" << v->id() << ": tolerance box too small to randomize the anchor position");
        }

        return { p.x() + FT(dx), p.y() + FT(dy), p.z() + FT(dz) };
      };

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
        double max_violation = 0.; // worst (max_i |d_i| / tol_i) over incident triplets
        FacetSPtr worst_facet;
      };

      struct Vertex_stability_record
      {
        VertexSPtr vertex;
        FacetSPtr facet;
        double max_violation;
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
          if (l.max_violation != r.max_violation)
            return l.max_violation > r.max_violation;
          return l.vertex_id < r.vertex_id;
        }

        const std::vector<Vertex_stability_record>* records;
      };

      auto evaluate_vertex_stability = [&](const VertexSPtr& v) -> Vertex_stability_info
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Check stability of V" << v->id() << " (deg: " << v->degree() << ")");

        const Point_3& p = original_points.at(v);
        const std::array<double, 3>& tol = vertex_tolerances.at(v);
        CGAL_precondition(tol[0] > 0. && tol[1] > 0. && tol[2] > 0.);
        const FT tol_ft[3] = { FT(tol[0]), FT(tol[1]), FT(tol[2]) }; // exact

        Vertex_stability_info info;

        for (auto it1 = v->facets().begin(); it1 != v->facets().end(); ++it1) {
          FacetSPtr f1 = it1->lock();
          if (!f1) continue;
          for (auto it2 = std::next(it1); it2 != v->facets().end(); ++it2) {
            FacetSPtr f2 = it2->lock();
            if (!f2) continue;
            for (auto it3 = std::next(it2); it3 != v->facets().end(); ++it3) {
              FacetSPtr f3 = it3->lock();
              if (!f3) continue;

              std::optional<Point_3> q = Kernel_wrapper::intersection(f1->get_plane(),
                                                                      f2->get_plane(),
                                                                      f3->get_plane());
              if (!q.has_value()) {
                CGAL_SS3_TRANSF_TRACE_V(1, "Warning: triplet of planes does not define a point!");
                CGAL_SS3_TRANSF_TRACE_V(1, "  faces: " << f1->id() << " " << f2->id() << " " << f3->id());
                continue;
              }

              const FT dx = CGAL::abs(q->x() - p.x());
              const FT dy = CGAL::abs(q->y() - p.y());
              const FT dz = CGAL::abs(q->z() - p.z());

              if (dx <= tol_ft[0] && dy <= tol_ft[1] && dz <= tol_ft[2]) {
                continue;
              }

              const double violation = (std::max)({ CGAL::to_double(dx) / tol[0],
                                                    CGAL::to_double(dy) / tol[1],
                                                    CGAL::to_double(dz) / tol[2] });

              if (info.worst_facet && violation <= info.max_violation) {
                continue;
              }

              // candidates: facets of the triplet that are neither frozen nor already anchored at v
              boost::container::small_vector<FacetSPtr, 3> cands;
              for (const FacetSPtr& cf : {f1, f2, f3}) {
                if (frozen.count(cf) != 0 || is_anchor_of(v, cf)) {
                  continue;
                }
                cands.push_back(cf);
              }

              if (cands.empty()) {
                // here means all three planes would pass through v->point(), so their common
                // point would be v->point(), which is inside the box.
                CGAL_SS3_TRANSF_TRACE_V(1, "Error: unstable triplet with no absorber at V" << v->id()
                    << " [F" << f1->id() << " F" << f2->id() << " F" << f3->id() << "]");
                std::abort();
              }

              CGAL_SS3_TRANSF_TRACE_V(64, "  V" << v->id() << " too far by " << violation
                                          << " [F" << f1->id() << " F" << f2->id() << " F" << f3->id() << "]");
              CGAL_SS3_TRANSF_TRACE_V(64, "    from " << p << " to " << *q);
              CGAL_SS3_TRANSF_TRACE_V(64, "    |d| " << dx << " " << dy << " " << dz
                                          << " VS " << tol[0] << " " << tol[1] << " " << tol[2]);

              info.is_stable = false;
              info.max_violation = violation;
              info.worst_facet = *std::min_element(cands.begin(), cands.end(),
                [&](const FacetSPtr& lf, const FacetSPtr& rf) -> bool
                {
                  // a triangle is wiggle room / does not need triangulating
                  if (lf->is_triangle() != rf->is_triangle()) {
                    return lf->is_triangle();
                  }
                  // avoid facets that would give a 4th anchor
                  const std::size_t la = anchor_count(lf), ra = anchor_count(rf);
                  if ((la == 3) != (ra == 3)) {
                    return la < ra;
                  }
                  // prefer small faces
                  const std::size_t sl = lf->vertices().size(), sr = rf->vertices().size();
                  if (sl != sr) {
                    return sl < sr;
                  }
                  return lf->id() < rf->id();
                });

              CGAL_SS3_TRANSF_TRACE_V(64, "    worst is F" << info.worst_facet->id());
            }
          }
        }

        CGAL_assertion(info.is_stable != bool(info.worst_facet));

        if (info.is_stable) {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is stable");
        } else {
          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " is unstable, exceeding its box by a factor "
              << info.max_violation);
        }

        return info;
      };

      auto recompute_facet_plane = [&](const FacetSPtr& f)
      {
        CGAL_SS3_TRANSF_TRACE_V(32, "Recompute F" << f->id());
        CGAL_SS3_TRANSF_TRACE_V(32, "  From coefficients [" << f->get_plane().a() << " " << f->get_plane().b() << " "
                                                            << f->get_plane().c() << " " << f->get_plane().d() << "]");

        CGAL_precondition(frozen.count(f) == 0);
        CGAL_precondition(original_planes.count(f));

        const auto& anchor_vertices = anchors[f];
        for (const VertexSPtr& v : anchor_vertices) {
          CGAL_SS3_TRANSF_TRACE_V(32, "    V" << v->id() << " is an anchor of F" << f->id() << " at " << v->point());
        }

        if (anchor_vertices.size() < 3) {
          perturbPlaneCoefficientsFixedPoints(f, nudge_range, anchor_vertices);
        } else {
          // 3 anchors: the plane is fully determined, with NO perturbation of its own. Its only
          // randomness comes from where the anchors were placed -- see the note at the top.
          const Point_3& p0 = anchor_vertices[0]->point();
          const Point_3& p1 = anchor_vertices[1]->point();
          const Point_3& p2 = anchor_vertices[2]->point();

          // With an unstable input the three anchors can be near-collinear (a sliver facet),
          // in which case the normal below -- and the orientation test against the original --
          // are meaningless.
          CGAL_assertion(!CGAL::collinear(p0, p1, p2));

          Plane_3 new_pl(p0, p1, p2);

          // The order of the anchors is arbitrary, so use the original normal for orientation
          const Plane_3& opl = original_planes.at(f);
          CGAL_assertion(new_pl.orthogonal_vector() * opl.orthogonal_vector() != 0);
          if (new_pl.orthogonal_vector() * opl.orthogonal_vector() < 0) {
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

      // ---------------------------------------------------------------------
      // Phase 0: initial perturbation
      // ---------------------------------------------------------------------

      CGAL_SS3_TRANSF_TRACE_V(16, "Initial plane perturbations");
      for (const FacetSPtr& f : polyhedron->facets()) {
        recompute_facet_plane(f);
      }

      CGAL_assertion(are_planes_in_general_position(polyhedron));

      // ---------------------------------------------------------------------
      // Phase 1a: greedy freezing
      // ---------------------------------------------------------------------
      //
      // Freeze a facet only if every one of its vertices still has an anchor position inside
      // its tolerance box. The test uses the CLOSEST point of the frozen structure, i.e. the
      // most permissive one: if even that is out of the box, no placement can succeed.
      {
        std::vector<FacetSPtr> order(polyhedron->facets().begin(), polyhedron->facets().end());
        std::sort(order.begin(), order.end(), [](const FacetSPtr& a, const FacetSPtr& b) {
          if (a->vertices().size() != b->vertices().size())
            return a->vertices().size() > b->vertices().size();
          return a->id() < b->id();
        });

        for (const FacetSPtr& f : order) {
          if (f->is_triangle()) {
            continue; // triangles are used for wiggle room
          }

          bool ok = true;
          for (const VertexSPtr& v : f->vertices()) {
            const std::optional<Point_3> q = project_on_frozen(v, original_points.at(v), f);
            if (!q.has_value() || !in_tolerance_box(v, *q)) {
              CGAL_SS3_TRANSF_TRACE_V(32, "cannot freeze F" << f->id() << " because of V" << v->id());
              ok = false;
              break;
            }
          }
          if (!ok) {
            continue;
          }

          frozen.insert(f);
          for (const VertexSPtr& v : f->vertices()) {
            frozen_at[v].push_back(f);
          }
          CGAL_SS3_TRANSF_TRACE_V(16, "freeze F" << f->id() << " (" << f->vertices().size() << " vertices)");
        }

        CGAL_SS3_TRANSF_TRACE_V(8, "frozen facets: " << frozen.size() << " / " << polyhedron->facets().size());
      }

      // ---------------------------------------------------------------------
      // Phase 1b: anchor positions
      // ---------------------------------------------------------------------
      //
      // A random point of the tolerance box, projected onto the frozen structure. Projection is
      // 1-Lipschitz but not box-preserving, so the amplitude is halved until the result lands
      // back inside; the closest point, which phase 1a guaranteed to be inside, is the fallback.
      {
        unsigned int n_fallback = 0;

        for (const VertexSPtr& v : polyhedron->vertices()) {
          const Point_3& p = original_points.at(v);
          const std::optional<Point_3> closest = project_on_frozen(v, p);
          CGAL_assertion(closest.has_value());
          CGAL_assertion(in_tolerance_box(v, *closest));

          Point_3 q = *closest;
          bool randomized = false;

          double shrink = 1.0;
          for (int attempt = 0; attempt < 8; ++attempt, shrink *= 0.5) {
            const Point_3 t = random_point_in_box(v, shrink);
            const std::optional<Point_3> r = project_on_frozen(v, t);
            CGAL_assertion(r.has_value());
            if (in_tolerance_box(v, *r)) {
              q = *r;
              randomized = true;
              break;
            }
          }

          if (!randomized) {
            // Only when the frozen structure sits right against the box boundary. The vertex
            // then contributes no randomness, which matters if it ends up as the 3rd anchor
            // of some facet.
            ++n_fallback;
            CGAL_SS3_TRANSF_TRACE_V(8, "V" << v->id() << ": no random anchor position available, "
                                        "falling back on the closest point");
          }

          v->set_point(q);

          CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " anchored at " << q
              << " [" << frozen_at[v].size() << " frozen], distance from original: "
              << CGAL::approximate_sqrt(CGAL::squared_distance(q, p)));
        }

        CGAL_SS3_TRANSF_TRACE_V(8, "vertices with no randomized anchor: " << n_fallback
            << " / " << polyhedron->vertices().size());
      }

      // Invariant (ii): a frozen facet passes exactly through all of its vertices.
      CGAL_assertion_code(
        for (const FacetSPtr& f : frozen)
          for (const VertexSPtr& v : f->vertices())
            CGAL_assertion(f->get_plane().has_on(v->point()));
      )

      // ---------------------------------------------------------------------
      // Phase 2: anchoring
      // ---------------------------------------------------------------------

      std::vector<Vertex_stability_record> unstable_vertices;
      CGAL::unordered_flat_map<VertexSPtr, std::size_t> unstable_record_indices;

#ifdef CGAL_SS3_DUMP_FILES
      std::ofstream out_unstable_base("results/base_unstable_vertices.xyz");
      out_unstable_base.precision(17);
#endif

      CGAL_SS3_TRANSF_TRACE_V(16, "Initial stable/unstable classification");
      for (const VertexSPtr& v : polyhedron->vertices()) {
        const Vertex_stability_info info = evaluate_vertex_stability(v);
        if (!info.is_stable) {
#ifdef CGAL_SS3_DUMP_FILES
          out_unstable_base << v->point() << "\n";
#endif
          unstable_record_indices[v] = unstable_vertices.size();
          unstable_vertices.push_back({v, info.worst_facet, info.max_violation, true,
                                       static_cast<std::size_t>(v->id())});
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
        const Vertex_stability_info info = evaluate_vertex_stability(v);

        const auto it = unstable_record_indices.find(v);
        if (it != unstable_record_indices.end()) {
          Vertex_stability_record& record = unstable_vertices[it->second];
          record.vertex = v;
          record.vertex_id = static_cast<std::size_t>(v->id());
          if (info.is_stable) {
            CGAL_SS3_TRANSF_TRACE_V(64, "V" << v->id() << " is now stable");
            if (pq.contains(it->second)) {
              pq.erase(it->second);
            }
            record.is_active = false;
            record.facet = nullptr;
            record.max_violation = 0.;
          } else {
            record.is_active = true;
            record.facet = info.worst_facet;
            record.max_violation = info.max_violation;
            if (pq.contains(it->second)) {
              pq.update(it->second);
            } else {
              pq.push(it->second);
            }
          }
        } else if (!info.is_stable) {
          const std::size_t index = unstable_vertices.size();
          unstable_vertices.push_back({v, info.worst_facet, info.max_violation, true,
                                       static_cast<std::size_t>(v->id())});
          unstable_record_indices[v] = index;
          CGAL_SS3_TRANSF_TRACE_V(64, "new unstable vertex");
          pq.push(index);
        }
      };

      while (!pq.empty()) {
        CGAL_SS3_TRANSF_TRACE_V(16, "stable/unstable: main loop (" << pq.size() << ")");

        const std::size_t unstable_index = pq.top_and_pop();
        const Vertex_stability_record record = unstable_vertices[unstable_index];
        if (!record.is_active) {
          continue;
        }

        const VertexSPtr v = record.vertex;
        const FacetSPtr f = record.facet;
        CGAL_SS3_TRANSF_TRACE_V(16, "pop V" << v->id() << " / F" << f->id()
                                    << " (violation " << record.max_violation << ")");

        CGAL_assertion(frozen.count(f) == 0);
        CGAL_assertion(!is_anchor_of(v, f));

        // The vertex is already at its anchor position (phase 1b); anchoring is only a matter
        // of pinning 'f' to it.
        all_anchors.insert(v);
        auto& facet_anchors = anchors[f];
        facet_anchors.push_back(v);

        std::vector<FacetSPtr> facets_to_recompute;

        // only a polygon can reach 4 anchors
        if (facet_anchors.size() > 3) {
          CGAL_assertion(!f->is_triangle());
          CGAL_SS3_TRANSF_TRACE_V(8, "  Must triangulate F" << f->id());

          const boost::container::small_vector<VertexSPtr, 3> old_anchors(facet_anchors.begin(),
                                                                          facet_anchors.end());
          const Plane_3 original_plane = original_planes.at(f); // intentional copy

          auto [_, new_facets] = Transformation::triangulate_facet(f, polyhedron);
          original_planes.erase(f);
          anchors.erase(f);

          for (const FacetSPtr& nf : new_facets) {
            original_planes[nf] = original_plane;
            // a triangle inherits the anchors that are among its vertices
            for (const VertexSPtr& w : nf->vertices()) {
              if (std::find(old_anchors.begin(), old_anchors.end(), w) != old_anchors.end()) {
                anchors[nf].push_back(w);
              }
            }
            CGAL_assertion(anchors[nf].size() <= 3);
            facets_to_recompute.push_back(nf);
          }
        } else {
          facets_to_recompute.push_back(f);
        }

        CGAL_SS3_TRANSF_TRACE_V(32, facets_to_recompute.size() << " facet(s) to recompute");
        for (const FacetSPtr& nf : facets_to_recompute) {
          recompute_facet_plane(nf);
          for (const VertexSPtr& iv : nf->vertices()) {
            update_vertex_record(iv);
          }
        }
      }

      for (const FacetSPtr& f : polyhedron->facets()) {
        CGAL_SS3_TRANSF_TRACE_V(32, "  F" << f->id()
                                      << (frozen.count(f) ? " [frozen]" : "")
                                      << " [" << anchor_count(f) << " anchors]"
                                      << " [" << f->get_plane() << "]");
      }
#endif // CGAL_SS3_PERTURB_V4_ANCHORING_V2

  #ifdef CGAL_SS3_DUMP_FILES
      IO::write_OBJ("results/V4_general_position.obj", polyhedron, parameters::do_not_triangulate_faces(true));
      IO::write_OBJ("results/V4_general_position-triangulated.obj", polyhedron, parameters::do_not_triangulate_faces(false));
  #endif

      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
      CGAL_SS3_TRANSF_TRACE_V(32, "V" << v->id() << " has length " << Size_shenanigans::length(v->point()));

      CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
      CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " has length " << Size_shenanigans::length(f->get_plane()));

      CGAL_postcondition_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
      CGAL_postcondition(is_stable(v));
      CGAL_postcondition_code(})

      CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
      CGAL_SS3_TRANSF_TRACE_V(32, "F" << f->id() << " has final plane " << f->get_plane());

      CGAL_postcondition(are_planes_in_general_position(polyhedron));

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
        CGAL_assertion_code(if (vertices_to_check.count(v)) { continue; })
        CGAL_assertion(is_stable(v));
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
          if (is_stable(v)) {
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
          CGAL_SS3_TRANSF_TRACE_V(8, "All vertices stable");
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
      CGAL_assertion(is_stable(v));
    }

    // -- PART 2 --
    // The perturbation planes are now set up, split high-degree vertices

    CGAL_SS3_TRANSF_TRACE_V(8, "Part 2: split high-degree vertices");

    std::list<VertexSPtr> vertices_tosplit;
    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (v->degree() > 3) {
        vertices_tosplit.push_back(v);
      }
    }

    CGAL_SS3_TRANSF_TRACE_V(8, vertices_tosplit.size() << " vertices to split");

    for (const VertexSPtr& vertex : vertices_tosplit) {
      CGAL_SS3_TRANSF_TRACE_V(8, "Splitting " << vertex->to_string());

      vertex->sort();

      // soup to be used in the arrangement
      std::vector<Point_3> points;
      std::vector<std::vector<std::size_t> > polygons;
      std::vector<FacetSPtr> polygon_to_facet;

      // Create a sufficiently large bounding box containing all plane intersections
      CGAL_assertion(is_stable(vertex));

      CGAL::Bbox_3 vbb = original_points[vertex].bbox();
      const std::array<double, 3>& tolerances = vertex_tolerances.at(vertex);
      CGAL::Bbox_3 bb = { vbb.xmin() - tolerances[0], vbb.ymin() - tolerances[1], vbb.zmin() - tolerances[2],
                          vbb.xmax() + tolerances[0], vbb.ymax() + tolerances[1], vbb.zmax() + tolerances[2] };

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
      using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;

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

#if defined(CGAL_SS3_ENABLE_TRACE) || defined(CGAL_SS3_DUMP_FILES)
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

#if defined(CGAL_SS3_ENABLE_TRACE) || defined(CGAL_SS3_DUMP_FILES)
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

      CGAL_SS3_TRANSF_TRACE_V(32, "building volumes...");

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

      CGAL_SS3_TRANSF_TRACE_V(32, volume_CCs.size() << " volume CCs");

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
        UNINITIALIZED,
        TBD,
        INSIDE,
        OUTSIDE
      };

      std::vector<CC_in_out_flag> in_out_flags(volume_CCs.size(), CC_in_out_flag::TBD);

      // Sanity checks: a boundary facet cannot have an OUTSIDE volume on its bottom
      // and an INSIDE volume on its top
      for (std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        VID bot_vid = face_volume_IDs[i][0];
        VID top_vid = face_volume_IDs[i][1];
        CGAL_assertion(!(in_out_flags[top_vid] == CC_in_out_flag::INSIDE && in_out_flags[bot_vid] == CC_in_out_flag::OUTSIDE));
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

      // - for all facets incident to the vertex, there should be at least one polygon
      // that has an INSIDE volume on its bottom and an OUTSIDE volume on its top
      CGAL::unordered_flat_set<FacetSPtr> found_facets;
      for (std::size_t i=0; i<polygons.size(); ++i) {
        if (!polygon_to_facet[i])
          continue;

        CGAL_assertion(face_volume_IDs[i][0] != VID(-1));
        CGAL_assertion(face_volume_IDs[i][1] != VID(-1));
        VID bot_vid = face_volume_IDs[i][0];
        VID top_vid = face_volume_IDs[i][1];
        CGAL_assertion(in_out_flags[bot_vid] != CC_in_out_flag::UNINITIALIZED);
        CGAL_assertion(in_out_flags[top_vid] != CC_in_out_flag::UNINITIALIZED);
        if (in_out_flags[bot_vid] == CC_in_out_flag::INSIDE && in_out_flags[top_vid] == CC_in_out_flag::OUTSIDE) {
          found_facets.insert(polygon_to_facet[i]);
        }
      }
      CGAL_assertion(found_facets.size() == vertex->degree());

      // Now, we want to deal with the tentative cells
      bool tentative_cell_n = 0;

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
      using operations_research::sat::ValidateCpModel;

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

      CGAL_SS3_TRANSF_TRACE_V(32, "Cell unknowns: " << x_unknowns << " / " << C);

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

      CGAL_SS3_TRANSF_TRACE_V(32, "Face unknowns: " << b_unknowns << " / " << polygons.size());

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

#define CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

      // =====================================================================
      //  CONSTRAINT — Global vertex manifoldness
      // =====================================================================

      // @todo if no using non manifold vertex constraints

      // =====================================================================
      //  CONSTRAINT — Local vertex manifoldness (per facet color)
      // =====================================================================

      // @fixme the non lazy constraints are more expensive and are probably wrong:
      // we could still pinch within the same CC even with a single root
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
#endif // CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS

      // =====================================================================
      //  SOLVE
      // =====================================================================

      std::vector<CC_in_out_flag> solution(volume_CCs.size(), CC_in_out_flag::UNINITIALIZED);

#ifndef CGAL_SLS3_USE_LAZY_NON_MANIFOLD_VERTEX_CONSTRAINTS
      const auto& proto = model.Build();
      CGAL_SS3_TRANSF_TRACE_V(32, "Variables: " << proto.variables_size());
      CGAL_SS3_TRANSF_TRACE_V(32, "Constraints: " << proto.constraints_size());

      SatParameters params;
      params.set_stop_after_first_solution(true);
      params.set_enumerate_all_solutions(false);

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
        CGAL_SS3_TRANSF_TRACE_V(8, "--- Starting solver iteration " << iteration << " ---");

        const auto& proto = model.Build();
        CGAL_SS3_TRANSF_TRACE_V(32, "Variables: " << proto.variables_size());
        CGAL_SS3_TRANSF_TRACE_V(32, "Constraints: " << proto.constraints_size());

        CGAL_SS3_TRANSF_TRACE_V(32, "Model:\n" << CpModelStats(proto));

#ifndef NDEBUG
        const std::string validation_error = ValidateCpModel(proto);
        if (!validation_error.empty()) {
          CGAL_SS3_TRANSF_TRACE_V(1, "Error: invalid model (" << validation_error << ")");
          std::abort();
        } else {
          CGAL_SS3_TRANSF_TRACE_V(8, "Model validated successfully");
        }
#endif

        SatParameters params;
        params.set_stop_after_first_solution(true);
        params.set_enumerate_all_solutions(false);

        params.set_random_seed(0);

        // workaround some deadlock issues within absl, and determinism
        params.set_num_workers(1);

        // debug
        // params.set_cp_model_presolve(false);
        // params.set_symmetry_level(0);
        // params.set_max_time_in_seconds(60.0);
        // params.set_log_search_progress(true);
        // params.set_log_to_stdout(true);

        bool ok = false;
        const CpSolverResponse response = SolveWithParameters(proto, params);
        CGAL_SS3_TRANSF_TRACE_V(32, "status = " << CpSolverStatus_Name(response.status()));
        CGAL_SS3_TRANSF_TRACE_V(32, "info   = " << response.solution_info());
        CGAL_SS3_TRANSF_TRACE_V(32, CpSolverResponseStats(response));

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
#ifdef CGAL_SS3_DUMP_FILES
          std::ofstream out("result/failing_model.pb.txt");
          out.precision(17);
          out << proto.DebugString();
#endif
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
          for (const auto& e : dumps) {
            std::vector<Point_3> cc_points = points;
            std::vector<std::vector<PID> > cc_polygons;

            for (std::size_t i=0; i<volume_CCs.size(); ++i) {
              CGAL_assertion(e.first != CC_in_out_flag::TBD);
              CGAL_assertion(e.first != CC_in_out_flag::UNINITIALIZED);

              if (solution[i] != e.first)
                continue;

              for (FID fid : volume_CCs[i]) {
                cc_polygons.push_back(polygons[fid]);
              }

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

          using PointRange = std::vector<PID>;
          using PolygonRange = std::vector<std::vector<PID> >;
          using Orienter = PMP::internal::Polygon_soup_orienter<PointRange, PolygonRange>;

          typename Orienter::Edge_map edges(points.size());
          typename Orienter::Marked_edges marked_edges;
          Orienter::fill_edge_map(edges, marked_edges, local_polygons);
          CGAL_assertion(marked_edges.empty()); // no NM edges is part of the constraints

          Orienter::has_singular_vertices(points.size(), local_polygons, edges, marked_edges, nm_vertex_id);

          return (nm_vertex_id == -1);
        };

        if (!check_CC()) {
          CGAL_SS3_TRANSF_TRACE_V(32, "Issue with global boundary");
        } else {
          for (FacetWPtr wf : vertex->facets()) {
            if (FacetSPtr f = wf.lock()) {
              if (!check_CC(f)) {
                CGAL_SS3_TRANSF_TRACE_V(32, "Issue with CC of F" << f->id());
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
              some_not_fixed = true; // at least one cell can be flipped
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
          CGAL_SS3_TRANSF_TRACE_V(32, "Valid solution found");
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

            if (solution[i] != e.first) {
              continue;
            }

            for (FID fid : volume_CCs[i]) {
              cc_polygons.push_back(polygons[fid]);
            }

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

      auto reconstruct_combi = [&]() -> combi
      {
        auto create_split = [](int begin, int end) -> vec2i
        {
          vec2i result(new int[2]);
          result[0] = begin;
          result[1] = end;
          return result;
        };

        std::set<std::pair<int, int> > unique_splits;
        for (edge_descriptor e : edges(bsm)) {
          if (is_border(e, bsm)) {
            continue;
          }

          halfedge_descriptor h = halfedge(e, bsm);
          face_descriptor f1 = face(h, bsm);
          face_descriptor f2 = face(opposite(h, bsm), bsm);

          // map mesh faces -> original input facets
          // (use the type produced by ifpm; here assumed comparable + indexable)
          FacetSPtr input_facet_1 = get(ifpm, f1);
          FacetSPtr input_facet_2 = get(ifpm, f2);
          CGAL_assertion(input_facet_1 != FacetSPtr() && input_facet_2 != FacetSPtr());

          if (input_facet_1 == input_facet_2) {
            continue;
          }

          // if the facets are neighbors, it's not a split edge
          EdgeSPtr common_e = input_facet_1->find_edge(input_facet_2);
          if (common_e != EdgeSPtr()) {
            continue;
          }

          // map original facets -> fan positions
          int a = local_facet_indices.at(input_facet_1);
          int b = local_facet_indices.at(input_facet_2);

          // canonicalize so begin < end (matches create_single_split_combinations())
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
        // compare_splits returns +1 when split1 < split2, lexicographically
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

      // @todo use split_vertex(vertex, split_combi) directly...?
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
  }

  // Perturbation to ensure generic configuration.
  // We always need to ensure that points are exactly on the planes of their incident facets.
  static void apply_rand_perturbation(PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE_V(4, "Applying random perturbations to the polyhedron...");

    // Generic approach
    Transformation::normalize_facet_planes(polyhedron); // @todo hasn't this already been done before?

    ConfigurationSPtr config = Configuration::get_instance();
    const bool safe_mode = config->get_Boolean("Preprocessing", "check_degenerate_configuration");
    const auto seed = config->get_ull("Preprocessing", "seed");

    CGAL_SS3_TRANSF_TRACE_V(4, "Seed to: " << seed);
    gen().seed(seed);

    PolyhedronSPtr p_mem;
    if (safe_mode) {
      p_mem = polyhedron->clone();
    }

#ifdef CGAL_SS3_RUN_TIMERS
    CGAL::Real_timer timer;
    timer.start();
#endif

    apply_plane_perturbation_V4(polyhedron);

#ifdef CGAL_SS3_RUN_TIMERS
    timer.stop();
    CGAL_SS3_TRANSF_TRACE_V(4, "perturbation time: " << timer.time());
#endif

    if (safe_mode) {
      CGAL_SS3_TRANSF_TRACE_V(8, "Safe mode is enabled, checking validity of the perturbation...");

      for (;;) {
#ifdef CGAL_SS3_DUMP_FILES
        IO::write_OBJ("results/last_perturbation.obj", polyhedron,
                      parameters::stream_precision(17).do_not_triangulate_faces(true));
#endif

        // @fixme
        using G = CGAL::Surface_mesh<Point_3>;

        if (are_planes_in_general_position(polyhedron) &&
            !Self_intersection::template has_self_intersecting_triangulated_surface<G>(polyhedron)) {
          CGAL_SS3_TRANSF_TRACE_V(8, "Successful perturbation");
          break;
        }

        CGAL_SS3_TRANSF_TRACE_V(4, "Warning: perturbation failed, retrying...");

        polyhedron = p_mem->clone();
        if (!is_triangle_polyhedron(polyhedron)) {
          Transformation::triangulate_facets(polyhedron);
          p_mem = polyhedron->clone();
        }

        rand_move_points(polyhedron);
      }
    } else {
      CGAL_assertion(are_planes_in_general_position(polyhedron));
      CGAL_assertion(!Self_intersection::has_self_intersecting_surface(polyhedron));
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
