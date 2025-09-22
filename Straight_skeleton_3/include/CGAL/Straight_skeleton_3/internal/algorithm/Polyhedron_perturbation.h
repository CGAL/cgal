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
 * @file   algo/3d/SelfIntersection.h
 * @author Gernot Walzl
 * @date   2012-07-18
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H

#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_factory.h>
#include <CGAL/Straight_skeleton_3/internal/kernel/Kernel_wrapper.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/enum.h>
#include <CGAL/Random.h>

#include <limits>
#include <list>
#include <random>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class PolyhedronPerturbation
{
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Line_3 = typename Traits::Line_3;
  using Plane_3 = typename Traits::Plane_3;

  using Point3SPtr = std::shared_ptr<Point_3>;
  using Line3SPtr = std::shared_ptr<Line_3>;
  using Plane3SPtr = std::shared_ptr<Plane_3>;

private:
  using Polyhedron = HDS::Polyhedron<Traits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::template Vertex<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::template Facet<Traits>;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using SkelFacetData = typename Polyhedron::SkelFacetData;

private:
  using KernelFactory = kernel::KernelFactory<Traits>;
  using KernelWrapper = kernel::KernelWrapper<Traits>;
  using Transformation = algorithm::PolyhedronTransformation<Traits>;
  using HdsUtils = algorithm::HdsUtils<Traits>;

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
  static std::array<double, 3> randVec(double min, double max)
  {
    static std::random_device rd;
    unsigned int s = 0; // rd()
    // CGAL_SS3_TRANSF_TRACE("seed = " << s);
    static std::mt19937 gen(s);
    static std::uniform_real_distribution<> rdist(min, max);

    return { rdist(gen), rdist(gen), rdist(gen) };
  }

public:
  /**
    * Checking for degeneracies
    */
  static bool doAll2PlanesIntersect(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Check all 2-combinations of planes");
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    typename std::list<FacetSPtr>::iterator it_f1 = polyhedron->facets().begin();
    while (it_f1 != polyhedron->facets().end()) {
      FacetSPtr facet1 = *it_f1++;
      typename std::list<FacetSPtr>::iterator it_f2 = it_f1;
      while (it_f2 != polyhedron->facets().end()) {
        FacetSPtr facet2 = *it_f2++;

        // Do not use CGAL::do_intersect: here we want to check that the result is a point
        if (!KernelWrapper::intersection(facet1->plane(), facet2->plane())) {
          CGAL_SS3_TRANSF_TRACE("Degenerate facet pair:");
          CGAL_SS3_TRANSF_TRACE("  " << facet1->toString());
          CGAL_SS3_TRANSF_TRACE("  " << facet2->toString());
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

  static bool doAll3PlanesIntersect(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Check all 3-combinations of planes");
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
          if (!KernelWrapper::intersection(facet1->plane(), facet2->plane(), facet3->plane())) {
            CGAL_SS3_TRANSF_TRACE("Degenerate facet triplet:");
            CGAL_SS3_TRANSF_TRACE("  " << facet1->toString());
            CGAL_SS3_TRANSF_TRACE("  " << facet2->toString());
            CGAL_SS3_TRANSF_TRACE("  " << facet3->toString());
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

  static bool areAllVerticesDegree3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      if (vertex->degree() != 3) {
        CGAL_SS3_TRANSF_TRACE("High-degree vertex: " << vertex->toString());
        result = false;
        break;
      }
    }
    return result;
  }

  /**
  * Nudge the plane coefficients by a random value in the range [low, high].
  */
  static void perturbPlaneCoefficientsNudge(const FacetSPtr& facet,
                                            const double range)
  {
    CGAL_precondition(Transformation::hasNormalizedPlane(facet));

    CGAL_SS3_TRANSF_TRACE_V(16, "Nudging (Nudge) Facet " << facet->getID());
    CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << facet->plane()->a() << " " << facet->plane()->b() << " "
                                                        << facet->plane()->c() << " " << facet->plane()->d() << "]");

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

    double na = nudge(facet->plane()->a());
    double nb = nudge(facet->plane()->b());
    double nc = nudge(facet->plane()->c());
    double nd = nudge(facet->plane()->d()); // @todo do not nudge 'd'? (mind the 'to_double()')

    double n = CGAL::approximate_sqrt(CGAL::square(na) + CGAL::square(nb) + CGAL::square(nc));
    CGAL_assertion(n != 0); // should not happen since we have normalized and the shift is tiny

    // below doesn't seem to matter? Probably need specific static filters...
#if 0
    facet->setPlane(KernelFactory::createPlane3(na/n, nb/n, nc/n, nd/n));
#else
    // cast to_double() *after* the normalization to have double coordinates in the planes
    // the downside is that we won't have a^2 + b^2 + c^2 == 1,
    // but then again, who does...
    const double a = CGAL::to_double(na/n);
    const double b = CGAL::to_double(nb/n);
    const double c = CGAL::to_double(nc/n);
    const double d = CGAL::to_double(nd/n);
    facet->setPlane(KernelFactory::createPlane3(a, b, c, d));
#endif

    CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << facet->getPlane()->a() << " " << facet->getPlane()->b() << " "
                                                      << facet->getPlane()->c() << " " << facet->getPlane()->d() << "]");

    CGAL_postcondition(Transformation::hasNormalizedPlane(facet));
  }

  /**
  * Nudge the plane coefficients but ensure that the perturbed plane goes through 0, 1, or 2 fixed points.
  * If 0 points: nudge all coefficients independently.
  * If 1 point: nudge (a, b, c), recompute d so the plane passes through the point.
  * If 2 points: nudge (a, b, c) with the constraint that the new plane passes through both points.
  */
  static void perturbPlaneCoefficientsFixedPoints(const FacetSPtr& facet,
                                                  const double range,
                                                  const std::vector<Point3SPtr>& fixed_points)
  {
    CGAL_precondition(Transformation::hasNormalizedPlane(facet));
    CGAL_precondition(fixed_points.size() <= 2);

    CGAL_SS3_TRANSF_TRACE_V(16, "Nudging Facet " << facet->getID());
    CGAL_SS3_TRANSF_TRACE_V(16, "  From coefficients [" << facet->plane()->a() << " " << facet->plane()->b() << " "
                                                        << facet->plane()->c() << " " << facet->plane()->d() << "]");
    CGAL_SS3_TRANSF_TRACE_V(16, "  with " << fixed_points.size() << " fixed points");
    CGAL_SS3_TRANSF_TRACE_CODE(for (Point3SPtr fp : fixed_points))
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
      FT na = nudge_to_simplest_rational_in_interval(facet->plane()->a());
      FT nb = nudge_to_simplest_rational_in_interval(facet->plane()->b());
      FT nc = nudge_to_simplest_rational_in_interval(facet->plane()->c());
      FT nd = nudge_to_simplest_rational_in_interval(facet->plane()->d());
#else
      double na = nudge(facet->plane()->a());
      double nb = nudge(facet->plane()->b());
      double nc = nudge(facet->plane()->c());
      double nd = nudge(facet->plane()->d());
#endif
      facet->setPlane(KernelFactory::createPlane3(na, nb, nc, nd));
    } else if (fixed_points.size() == 1) {
      // 1 fixed point: nudge (a, b, c), recompute d so the plane passes through the point
      Point3SPtr p0 = fixed_points[0];

#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
      FT na = nudge_to_simplest_rational_in_interval(facet->plane()->a());
      FT nb = nudge_to_simplest_rational_in_interval(facet->plane()->b());
      FT nc = nudge_to_simplest_rational_in_interval(facet->plane()->c());
#else
      double na = nudge(facet->plane()->a());
      double nb = nudge(facet->plane()->b());
      double nc = nudge(facet->plane()->c());
#endif
      const FT& x0 = p0->x();
      const FT& y0 = p0->y();
      const FT& z0 = p0->z();
      FT d = - (na * x0 + nb * y0 + nc * z0);
      facet->setPlane(KernelFactory::createPlane3(na, nb, nc, d));
      CGAL_postcondition(facet->getPlane()->has_on(*p0));
    } else if (fixed_points.size() == 2) {
      // 2 fixed points: construct a plane through both points, nudge the normal within the allowed family
      Point3SPtr p0 = fixed_points[0];
      Point3SPtr p1 = fixed_points[1];
      CGAL_assertion(*p0 != *p1);

      const FT& p0x = p0->x();
      const FT& p0y = p0->y();
      const FT& p0z = p0->z();
      const FT& p1x = p1->x();
      const FT& p1y = p1->y();
      const FT& p1z = p1->z();

      // Direction vector between points
      FT ux = p1x - p0x;
      FT uy = p1y - p0y;
      FT uz = p1z - p0z;
      FT uu = ux*ux + uy*uy + uz*uz;

      // Original normal
      const FT& a0 = facet->plane()->a();
      const FT& b0 = facet->plane()->b();
      const FT& c0 = facet->plane()->c();

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
      facet->setPlane(KernelFactory::createPlane3(a1, b1, c1, d1));

      CGAL_postcondition(facet->getPlane()->has_on(*p0));
      CGAL_postcondition(facet->getPlane()->has_on(*p1));
    } else {
      CGAL_SS3_TRANSF_TRACE("Error: called fixed point facet perturbation with > 2 fixed points");
    }

    CGAL_SS3_TRANSF_TRACE_V(16, "  To coefficients [" << facet->getPlane()->a() << " " << facet->getPlane()->b() << " "
                                                      << facet->getPlane()->c() << " " << facet->getPlane()->d() << "]");

    CGAL_postcondition(Transformation::hasNormalizedPlane(facet));
  }

  static void perturbPlaneCoefficientsHighDegrees(const FacetSPtr& facet,
                                                  const double range)
  {
    std::vector<Point3SPtr> high_degree_points;
    for (const VertexSPtr& v : facet->vertices()) {
      if (v->degree() > 3) {
        high_degree_points.push_back(v->getPoint());
      }
    }
    return perturbPlaneCoefficientsFixedPoints(facet, range, high_degree_points);
  }

  /**
    * Check that all faces have at most two high-degree vertices: a facet with fewer than two high-degree
    * vertices can be perturbed by nudging the high-degree vertices, and pivoting the facet randomly
    * around these fixed points.
    */
  static bool isTiltCompatible(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);
    bool result = true;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->numHighDegreeVertices() > 2) {
        CGAL_SS3_TRANSF_TRACE("facet " << facet->getID() << " has too many high-degree vertices "
                                        << "(" << facet->numHighDegreeVertices() << ")");
        result = false;
        break;
      }
    }
    return result;
  }

  static void randTiltPlanes(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    double range = 1e-10;
    ConfigurationSPtr config = Configuration::getInstance();
    if (config->isLoaded()) {
      range = config->getDouble("main", "perturbation_epsilon");
    }

    // If we only nudged planes with fixed point constraints, we might not ensure generic position,
    // for example if two pairs of constraints are along the same line.
    //
    // @todo could restrict to only high-degree vertices in facets that have 2 high-degree vertices
    for (VertexSPtr vertex : polyhedron->vertices()) {
      Point3SPtr p = vertex->getPoint();
      std::array<double, 3> v_r = randVec(-range, range);

      double px = CGAL::to_double(p->x()) + v_r[0];
      double py = CGAL::to_double(p->y()) + v_r[1];
      double pz = CGAL::to_double(p->z()) + v_r[2];

      Point3SPtr p_t = KernelFactory::createPoint3(px, py, pz);
      vertex->setPoint(p_t);
    }

    for (FacetSPtr facet : polyhedron->facets()) {
        perturbPlaneCoefficientsHighDegrees(facet, range);
    }
  }

  static void randTiltPlanesv3(const PolyhedronSPtr& polyhedron)
  {
    CGAL_SS3_TRANSF_TRACE("Random Plane Tilt (v3)");
    CGAL_SS3_DEBUG_SPTR(polyhedron);

    double range = 1e-10;
    ConfigurationSPtr config = Configuration::getInstance();
    if (config->isLoaded()) {
      range = config->getDouble("main", "perturbation_epsilon");
    }

    CGAL_SS3_TRANSF_TRACE("  perturbation_epsilon = " << range);

    if (isTiltCompatible(polyhedron)) {
      CGAL_SS3_TRANSF_TRACE("Polyhedron can simply be tilted immediately");
      randTiltPlanes(polyhedron);
      CGAL_assertion_code(bool success =)
      Transformation::resetPoints(polyhedron);
      CGAL_assertion(success);
      return;
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/tilt-v3_input.obj", polyhedron,
                      false /*do_triangulate*/,
                      true /*convert_to_double*/);
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
      auto pcdt = constructFacetTriangulation<CDT2_Tag>(f);

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
      CGAL_SS3_TRANSF_TRACE_V(16, "Checking if V" << v->getID() << " (deg: " << v->degree() << ") is fixed");
      auto it = determining_facets.find(v);
      return (it != determining_facets.end() && it->second.size() == 3);
    };

    auto is_facet_fixed = [&](const FacetSPtr& f) -> bool
    {
      CGAL_SS3_TRANSF_TRACE_V(16, "Checking if F" << f->getID() << " (" << f->vertices().size() << " nv) is fixed");
      CGAL_SS3_TRANSF_TRACE_V(16, "  fixing_vertices size: " << fixing_vertices[f].size());
      CGAL_assertion(fixing_vertices[f].size() <= 3);
      return (f->isTriangle() && fixing_vertices[f].size() == 3) ||
              (!f->isTriangle() && fixing_vertices[f].size() == 2);
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

      CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->getID() << " has " << a_determined_n << " determined vertices");
      CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->getID() << " has " << b_determined_n << " determined vertices");

      if (a_determined_n != b_determined_n) {
        // Give priority to the one with the least determined vertices
        return a_determined_n < b_determined_n;
      }

      // same number of determined vertices, give priority to the facet with the most high-degree vertices
      unsigned int a_hdv_n = hdv_count(a);
      unsigned int b_hdv_n = hdv_count(b);

      CGAL_SS3_TRANSF_TRACE_V(64, "F" << a->getID() << " has " << a_hdv_n << " high-degree vertices");
      CGAL_SS3_TRANSF_TRACE_V(64, "F" << b->getID() << " has " << b_hdv_n << " high-degree vertices");

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
        CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << facet->getID());
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
          CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is determined by F" << facet->getID() << " (a)");
        }
      }
    }

    // This is the main list of facets that we will process
    std::list<FacetSPtr> facets_to_process;
    for (const FacetSPtr& facet : polyhedron->facets()) {
      if (facet->isTriangle() || !has_high_degree_vertices(facet)) {
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
          VertexSPtr tv = e->dst(f);
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

              CGAL_SS3_TRANSF_TRACE("Facet F" << fprime->getID() << " needs triangulating due to missing high-degree edge between V" << sv->getID() << " and V" << tv->getID());

              CGAL_SS3_TRANSF_TRACE("Triangulate F" << fprime->getID());
              CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

              Transformation::triangulateFacet(fprime, polyhedron);

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
    IO::OBJFile::save("results/tilt-v3_preprocessed.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    // Forward declarations for mutually recursive lambdas
    std::function<bool(FacetSPtr, VertexSPtr)> add_fixing_vertex;
    std::function<void(VertexSPtr)> determine_vertex;

    auto nudge_constrained_vertex = [&](const VertexSPtr& v)
    {
      CGAL_SS3_TRANSF_TRACE("  Nudging V" << v->getID() << " from " << *(v->getPoint()));

      std::vector<Plane3SPtr> constraining_planes;
      for (const FacetSPtr& df : determining_facets[v]) {
        if (is_facet_fixed(df)) {
          constraining_planes.push_back(df->getPlane());
          CGAL_SS3_TRANSF_TRACE("    F" << df->getID() << " constrains the nudge");
        }
      }

      CGAL_assertion(constraining_planes.size() <= 3);

      const size_t n_fixed = constraining_planes.size();
      if (n_fixed == 3) {
        Transformation::resetPoint(v, { constraining_planes[0], constraining_planes[1], constraining_planes[2] });
        CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " reset to " << *v->getPoint());
        return;
      }

      Point3SPtr p = v->getPoint();

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
        FT nv = CGAL::simplest_rational_in_interval<typename Traits::Exact_kernel::FT>(d1, d2);
        return nv;
      };

      FT x = nudge_to_simplest_rational_in_interval(p->x());
      FT y = nudge_to_simplest_rational_in_interval(p->y());
      FT z = nudge_to_simplest_rational_in_interval(p->z());
#else
      std::array<double, 3> v_r = randVec(-range/2.0, range/2.0);
      double x = CGAL::to_double(p->x()) + v_r[0];
      double y = CGAL::to_double(p->y()) + v_r[1];
      double z = CGAL::to_double(p->z()) + v_r[2];
#endif
      Point3SPtr p_nudged = KernelFactory::createPoint3(x, y, z);
      // CGAL_SS3_TRANSF_TRACE("base nudge: " << x << " " << y << " " << z);

      Point3SPtr p_new;

      if (n_fixed == 0) {
        p_new = p_nudged;
      } else if (n_fixed == 1) {
        Plane3SPtr plane = constraining_planes[0];
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        // something similar but a little more subtle:
        // 1. project the point onto the plane
        // 2. express the point as a linear combination of the plane's origin and basis: pp = o + l1 * b1 + l2 * b2
        // 3. nudge l1 and l2 to l1' and l2' with a random interval around l1 and l2, and
        //    simplest_rational_in_interval
        // 4. recompute the point as pp = o + l1' * b1 + l2' * b2
        Point3SPtr pp = KernelWrapper::projection(plane, p_nudged);
        const Point_3& o = plane->point();
        const Vector_3& b1 = plane->base1();
        const Vector_3& b2 = plane->base2();
        FT l1 = CGAL::scalar_product(*pp - o, b1);
        FT l2 = CGAL::scalar_product(*pp - o, b2);
        FT nl1 = nudge_to_simplest_rational_in_interval(l1);
        FT nl2 = nudge_to_simplest_rational_in_interval(l2);
        p_new = KernelFactory::createPoint3(o.x() + nl1 * b1.x() + nl2 * b2.x(),
                                            o.y() + nl1 * b1.y() + nl2 * b2.y(),
                                            o.z() + nl1 * b1.z() + nl2 * b2.z());
#else
        p_new = KernelWrapper::projection(plane, p_nudged);
#endif
      } else if (n_fixed == 2) {
        Plane3SPtr plane1 = constraining_planes[0];
        Plane3SPtr plane2 = constraining_planes[1];
        Line3SPtr line = KernelWrapper::intersection(plane1, plane2);
#ifdef CGAL_SS3_USE_SIMPLEST_RATIONAL_IN_INTERVAL
        // something similar but a little more subtle:
        // 1. project the point onto the line
        // 2. express the point as a linear combination of the line's origin and basis: pp = o + l * v
        // 3. nudge l to l' with a random interval around l, and simplest_rational_in_interval
        // 4. recompute the point as pp = o + l' * v
        Point3SPtr pp = KernelWrapper::projection(line, p_nudged);
        const Point_3& o = line->point();
        const Vector_3& d = line->to_vector();
        FT l = CGAL::scalar_product(*pp - o, d);
        FT nl = nudge_to_simplest_rational_in_interval(l);
        p_new = KernelFactory::createPoint3(o.x() + nl * d.x(),
                                            o.y() + nl * d.y(),
                                            o.z() + nl * d.z());
#else
        p_new = KernelWrapper::projection(line, p_nudged);
#endif
      }

      CGAL_SS3_TRANSF_TRACE("  Nudged V" << v->getID() << " to " << *p_new);

      v->setPoint(p_new);
    };

    // sometimes we could fix as a polygon, but we need triangualte for other reasons
    auto should_triangulate_facet = [&](const FacetSPtr& f) -> bool
    {
      if (f->isTriangle() || is_facet_fixed(f)) {
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
          std::size_t l = Size_shenanigans::length(*(v->getPoint()));
          if (l > max_length) {
            CGAL_SS3_TRANSF_TRACE("Vertex V" << v->getID() << " is too long");
            CGAL_SS3_TRANSF_TRACE(CGAL::exact(*(v->getPoint())) << " (l=" << l << ")");
            CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " should be triangulated");
            return true;
          }
        }

        // if the vertex will be determined by the fixation of this facet, check the facets length
        if (it->second.size() == 2) {
          for (const FacetSPtr& of : determining_facets[v]) {
            std::size_t l = Size_shenanigans::length(*(of->getPlane()));
            if (l > max_length) {
              CGAL_SS3_TRANSF_TRACE("Facet F" << of->getID() << " is too long");
              CGAL_SS3_TRANSF_TRACE(CGAL::exact(*(of->getPlane())) << " (l=" << l << ")");
              CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " should be triangulated");
              return true;
            }
          }
        }
      }

      return false;
    };

    auto is_facet_overconstrained = [&](const FacetSPtr& f) -> bool
    {
      if (f->isTriangle() || is_facet_fixed(f)) {
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
            } else if (determining_facets[v].size() == 2 && ft->hasVertex(v)) {
              ++constrain_n;
            }
          }

          if (constrain_n > 2) {
            CGAL_SS3_TRANSF_TRACE("F" << ft->getID() << " would be over constrained by fixing of F" << f->getID());
            return true;
          }
        }
      }

      return false;
    };

    auto triangulate_facet = [&](const FacetSPtr& facet_tt)
    {
      CGAL_SS3_TRANSF_TRACE("Triangulate F" << facet_tt->getID());

      CGAL_assertion(!is_facet_fixed(facet_tt));

      // the facet is not yet fixed, so no vertex can have it as determining facet
      CGAL_assertion_code(for (const VertexSPtr& v : facet_tt->vertices()) {)
      CGAL_assertion(determining_facets[v].size() <= 3);
      CGAL_assertion(determining_facets[v].count(facet_tt) == 0);
      CGAL_assertion_code(})

      CGAL_SS3_TRANSF_TRACE_CODE(++had_to_triangulate_n;)

      auto [local_vertices, new_facets] = Transformation::triangulateFacet(facet_tt, polyhedron);

      for (const VertexSPtr& v : local_vertices) {
        CGAL_SS3_TRANSF_TRACE("local vertex " << v->getID() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)");

        if (is_vertex_determined(v)) {
          CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is already determined, skipping");
          continue;
        }

        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr fptr = wf.lock()) {
            if (is_facet_fixed(fptr)) {
              CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << fptr->getID() << " (c)");
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
        CGAL_SS3_TRANSF_TRACE("spawned F" << nf->getID());

        for (const VertexSPtr& iv : nf->vertices()) {
          if (is_vertex_determined(iv)) {
            CGAL_SS3_TRANSF_TRACE("newborn F" << nf->getID() << " is constrained by V" << iv->getID());
            fixing_vertices[nf].insert(iv);
          }
        }
      }
    };

    add_fixing_vertex = [&](const FacetSPtr& f, const VertexSPtr& v) -> bool
    {
      CGAL_precondition(fixing_vertices[f].size() <= 3);

      CGAL_SS3_TRANSF_TRACE("  Fix F" << f->getID() << " with V" << v->getID());

      if (is_facet_fixed(f)) {
        CGAL_SS3_TRANSF_TRACE("  F" << f->getID() << " is already fixed");
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
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "F" << f->getID() << " is now fixed by");
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : fixing_vertices[f]))
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->getID() << " [measure=" << Size_shenanigans::length(*(fv->getPoint())) << "]");
      CGAL_SS3_TRANSF_TRACE(ss.str());

      if (f->isTriangle()) {
        CGAL_assertion(fixing_vertices[f].size() == 3); // just to be clear

        // for triangles, all vertices are determined, and there is nothing to nudge
        // (note that vertices were themselves nudged so the facet is nudged).
        f->initPlane();
        Transformation::normalizePlaneCoefficients(f);

#ifdef CGAL_SS3_DUMP_FILES
        dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_3.OFF", f);
#endif

        // Here we do not need to add the fixed facet to incident determined vertices
        // because all vertices are already fully determined
        return true;
      }

      std::vector<Point3SPtr> fixed_points;
      for (const VertexSPtr& fixed_v : fixing_vertices[f]) {
        fixed_points.push_back(fixed_v->getPoint());
      }

      perturbPlaneCoefficientsFixedPoints(f, range, fixed_points);

      CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " is now fixed at " << *(f->getPlane()) << " [measure=" << Size_shenanigans::length(*(f->getPlane())) << "]");

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_fixed_" + std::to_string(fixed_points.size()) + ".OFF", f);
#endif

      CGAL_SS3_TRANSF_TRACE("Newly fixed facet F" << f->getID() << " determines its high-degree incident vertices...");

      // Need to now tag the vertices of the facet
      for (const VertexSPtr& v : f->vertices()) {
        CGAL_SS3_TRANSF_TRACE("incident: " << v->getID() << " (deg=" << v->degree() << "; " << determining_facets[v].size() << " determining facets)");
        if (!is_vertex_determined(v)) {
          determining_facets[v].insert(f);
          CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << f->getID() << " (b)");
        }
      }

      return true;
    };

    determine_vertex = [&](const VertexSPtr& v)
    {
      CGAL_precondition(is_high_degree(v));
      CGAL_precondition(is_vertex_determined(v));

      CGAL_SS3_TRANSF_TRACE_CODE(auto it = determining_facets[v].begin();)
      CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is now fully determined by"
                            << " F" << (*it)->getID() << " [measure=" << Size_shenanigans::length(*((*it)->getPlane()))
                            << "] F" << (*std::next(it))->getID() << " [measure=" << Size_shenanigans::length(*((*std::next(it))->getPlane()))
                            << "] F" << (*std::next(it, 2))->getID() << " [measure=" << Size_shenanigans::length(*((*std::next(it, 2))->getPlane())) << "]");

      // set the nudged position for the vertex: a nudge constrained by already fixed incident facets
      nudge_constrained_vertex(v);

      CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " is now determined at " << *(v->getPoint()) << " [measure=" << Size_shenanigans::length(*(v->getPoint())) << "]");

      // compute the plane coefficients of any incident facet that becomes fixed
      // by this vertex becoming determined
      for (FacetWPtr wf : v->facets()) {
        if (FacetSPtr f = wf.lock()) {
          add_fixing_vertex(f, v);
        }
      }
    };

    CGAL_SS3_TRANSF_TRACE("== Main facet flood... ==");

    while (!facets_to_process.empty()) {
      facets_to_process.sort(facet_sorter); // @todo priority queue...
      FacetSPtr facet = facets_to_process.front();
      facets_to_process.pop_front();

      CGAL_SS3_TRANSF_TRACE("Pop F" << facet->getID());

      CGAL_assertion(!facet->isTriangle());
      CGAL_assertion(fixing_vertices[facet].size() <= 2);

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "  Fixing vertices:";)
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& fv : fixing_vertices[facet]))
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " V" << fv->getID();)
      CGAL_SS3_TRANSF_TRACE(ss.str());

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
          CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << facet->getID() << " (d)");
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
    CGAL_SS3_TRANSF_TRACE("== Deal with remaining facets with high degree vertices... ==");

    for (const FacetSPtr& f : polyhedron->facets()) {
      if (f->isTriangle() || is_facet_fixed(f)) {
        continue;
      }

      CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << f->getID() << " [remaining]");

      std::vector<Point3SPtr> fixed_points;
      for (const VertexSPtr& v : fixing_vertices[f]) {
        CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is a fixing vertex");
        fixed_points.push_back(v->getPoint());
      }

      perturbPlaneCoefficientsFixedPoints(f, range, fixed_points);

      // fixing the facet cannot determine a vertex because the facet has already been visited

      // add random vertices to mark the facet as fixed
      static int dummy_id = -1;
      while (!is_facet_fixed(f)) {
        VertexSPtr dummy_v = Vertex::create(KernelFactory::createPoint3(0,0,0));
        dummy_v->setID(dummy_id--);
        fixing_vertices[f].insert(dummy_v);
      }

#ifdef CGAL_SS3_DUMP_FILES
      dump_facet("results/nudged_face_" + std::to_string(nudged_face_id++) + "_remaining.OFF", f);
#endif

      CGAL_postcondition(is_facet_fixed(f));
    }

    // At this point, everything that is not a high-degree triangular facet should be fixed
    for (const FacetSPtr& f : polyhedron->facets()) {
      if (f->isTriangle() || !has_high_degree_vertices(f)) {
        continue;
      }
      CGAL_assertion(is_facet_fixed(f));
    }

    // Nudge vertices that can still be nudged, for randomness
    CGAL_SS3_TRANSF_TRACE("== Nudge undetermined high-degree vertices... ==");

    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (is_high_degree(v) && !is_vertex_determined(v)) {
        CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is high degree and not fully determined, nudge it");
        nudge_constrained_vertex(v);

        // determine the vertex
        // since we know only triangle facets are left, we don't need to cascade and check
        // if incident facets become fixed
        for (FacetWPtr wf : v->facets()) {
          if (FacetSPtr f = wf.lock()) {
            if (!is_facet_fixed(f)) {
              // the facet cannot be without high-degree vertices since v is high degree
              CGAL_assertion(f->isTriangle());
              fixing_vertices[f].insert(v);
            }
          }
        }

        // add dummy facets to mark the vertex as determined
        static int dummy_id = -1;
        while (!is_vertex_determined(v)) {
          FacetSPtr dummy_f = Facet::create();
          dummy_f->setID(dummy_id--);
          determining_facets[v].insert(dummy_f);
        }

        CGAL_postcondition(is_vertex_determined(v));
      }
    }

    // Now handle triangle faces with high degrees
    CGAL_SS3_TRANSF_TRACE("== Deal with remaining triangles... ==");

    for (const FacetSPtr& f : polyhedron->facets()) {
      if (!f->isTriangle() || !has_high_degree_vertices(f)) {
        continue;
      }

      CGAL_SS3_TRANSF_TRACE_CODE(std::stringstream ss;)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "Fix triangle F" << f->getID() << " [");
      CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : f->vertices()) {)
      CGAL_SS3_TRANSF_TRACE_CODE(ss << "V" << v->getID());
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " (" << v->degree() << ")");
      CGAL_SS3_TRANSF_TRACE_CODE(if (is_vertex_determined(v)) { ss << "*"; })
      CGAL_SS3_TRANSF_TRACE_CODE(ss << " "; } ss << "]";)
      CGAL_SS3_TRANSF_TRACE(ss.str());

      CGAL_SS3_TRANSF_TRACE("Nudge and fix F" << f->getID() << " [triangle]");

      std::vector<Point3SPtr> fixed_points;
      for (const VertexSPtr& v : fixing_vertices[f]) {
        CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is a fixing vertex");
        fixed_points.push_back(v->getPoint());
      }

      if (fixed_points.size() == 3) {
        f->initPlane();
        Transformation::normalizePlaneCoefficients(f);
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
          CGAL_SS3_TRANSF_TRACE("  V" << v->getID() << " is determined by F" << f->getID() << " (f)");
          // no need to cascade here, we know only triangles are left
        }
      }

      for (const VertexSPtr& v : f->vertices()) {
        fixing_vertices[f].insert(v);
      }

      CGAL_postcondition(is_facet_fixed(f));
    }

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/tilt_v3-pre_reset.obj", polyhedron,
                      false /*do_triangulate*/,
                      true /*convert_to_double*/);
#endif

    CGAL_SS3_TRANSF_TRACE("Reset the position of not-fully-constrained vertices...");

    // Recompute all points which were not fixed (degree 3 vertices)
    for (const VertexSPtr& v : polyhedron->vertices()) {
      if (!is_high_degree(v)) {
        Transformation::resetPoint(v);

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

    CGAL_SS3_TRANSF_TRACE("All facets processed");

#ifdef CGAL_SS3_DUMP_FILES
    IO::OBJFile::save("results/tilt_v3.obj", polyhedron,
                           false /*do_triangulate*/,
                           true /*convert_to_double*/);
    IO::OBJFile::save("results/tilt_v3-triangulated.obj", polyhedron,
                           true /*do_triangulate*/,
                           true /*convert_to_double*/);
#endif

    CGAL_assertion_code(for (const VertexSPtr& v : polyhedron->vertices()) {)
    CGAL_assertion(is_vertex_determined(v));
    CGAL_assertion_code(})

    CGAL_assertion_code(for (const FacetSPtr& f : polyhedron->facets()) {)
    CGAL_assertion(fixing_vertices[f].size() <= 3);
    CGAL_assertion_code(})

    CGAL_assertion_code(for (const FacetSPtr& facet : polyhedron->facets()) {)
    CGAL_assertion_code(for (const VertexSPtr& v : facet->vertices()) {)
    CGAL_assertion_code(std::cout << "check V" << v->getID() << " on F" << facet->getID() << std::endl;)
    CGAL_assertion(facet->getPlane()->has_on(*(v->getPoint())));
    CGAL_assertion_code(})
    CGAL_assertion_code(})

    CGAL_SS3_TRANSF_TRACE("Had to triangulate " << had_to_triangulate_n << " facets");

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " has depth " << CGAL::depth(*(v->getPoint())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " has depth " << CGAL::depth(*(f->getPlane())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const VertexSPtr& v : polyhedron->vertices()))
    CGAL_SS3_TRANSF_TRACE("V" << v->getID() << " has length " << Size_shenanigans::length(*(v->getPoint())));

    CGAL_SS3_TRANSF_TRACE_CODE(for (const FacetSPtr& f : polyhedron->facets()) )
    CGAL_SS3_TRANSF_TRACE("F" << f->getID() << " has length " << Size_shenanigans::length(*(f->getPlane())));
  }
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_POLYHEDRON_PERTURBATION_H */
