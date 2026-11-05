// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Qijia Huang

#ifndef CGAL_MEDIAL_SKELETON_OFFSET_FUNCTION_H
#define CGAL_MEDIAL_SKELETON_OFFSET_FUNCTION_H

#include <CGAL/license/Variational_medial_axis.h>

#include <CGAL/Variational_medial_axis_sampling.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <limits>
#include <utility>

namespace CGAL {


#ifndef DOXYGEN_RUNNING
template <class TriangleMesh, class GeomTraits = Default> class Medial_skeleton_offset_function
{
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using MSkeleton = Medial_skeleton<TriangleMesh, GeomTraits>;


public:
  Medial_skeleton_offset_function(const MSkeleton& skeleton)
      : skeleton_(skeleton)
      , spheres_(&skeleton_.vertices()) {
    const std::size_t n = skeleton_.number_of_vertices();
    radii_.reserve(n);
    for(const auto& sph : *spheres_)
      radii_.push_back(CGAL::approximate_sqrt(sph.squared_radius()));
  }


  // Evaluate the minimum distance from point p to the interpolated medial axis (cones and slabs)
  FT operator()(const Point_3& p) const {
    FT dmin = (std::numeric_limits<FT>::max)();

    // TODO: the filtering strategy doesn't work
    std::unordered_set<std::pair<std::size_t, std::size_t>, boost::hash<std::pair<std::size_t, std::size_t>>>
        visited_edges;

    for(const auto& tri : skeleton_.faces()) {
      const std::size_t a = tri[0], b = tri[1], c = tri[2];
      const Point_3& c1 = (*spheres_)[a].center();
      const Point_3& c2 = (*spheres_)[b].center();
      const Point_3& c3 = (*spheres_)[c].center();
      const FT& r1 = radii_[a];
      const FT& r2 = radii_[b];
      const FT& r3 = radii_[c];
      /*
      const FT& r_max = r1>r2 ? (r1 > r3 ? r1 : r3) : (r2 > r3 ? r2 : r3);
      if(dmin > eval_slab(c1, c2, c3, r_max, r_max, r_max,p)) {
         dmin = (std::min)(dmin, eval_slab(c1, c2, c3, r1, r2, r3, p));
      }
      */
      dmin = (std::min)(dmin, eval_slab(c1, c2, c3, r1, r2, r3, p));
      visited_edges.insert({a, b});
      visited_edges.insert({b, c});
      visited_edges.insert({c, a});
      if(dmin < 0)
        return dmin;
    }
    for(auto [i, j] : skeleton_.edges()) {
      if(visited_edges.find({i, j}) != visited_edges.end())
        continue;
      const Point_3& c1 = (*spheres_)[i].center();
      const Point_3& c2 = (*spheres_)[j].center();
      const FT& r1 = radii_[i];
      const FT& r2 = radii_[j];
      /*
      const FT r_max = r1 > r2 ? r1 : r2;
      if(dmin > eval_cone(c1, c2, r_max, r_max, p))
        dmin = (std::min)(dmin, eval_cone(c1, c2, r1, r2, p));
      */
      dmin = (std::min)(dmin, eval_cone(c1, c2, r1, r2, p));
      if(dmin < 0)
        return dmin;
    }

    return dmin;
  }

private:
  FT sphere_distance(const Point_3& p, const Point_3& c, const FT r) const {
    return CGAL::approximate_sqrt((p - c).squared_length()) - r;
  }

  FT cone_distance(const Point_3& p, const Point_3& c1, const Point_3& c2, const FT r1, const FT r2, FT t) const {
    const Point_3 c = Point_3((FT(1) - t) * c1.x() + t * c2.x(), (FT(1) - t) * c1.y() + t * c2.y(),
                              (FT(1) - t) * c1.z() + t * c2.z());
    const FT r = (FT(1) - t) * r1 + t * r2;
    return sphere_distance(p, c, r);
  }

  FT slab_distance(const Point_3& p,
                   const Point_3& c1,
                   const Point_3& c2,
                   const Point_3& c3,
                   const FT& r1,
                   const FT& r2,
                   const FT& r3,
                   FT t1,
                   FT t2) const {
    const Point_3 c = Point_3(t1 * c1.x() + t2 * c2.x() + (FT(1) - t1 - t2) * c3.x(),
                              t1 * c1.y() + t2 * c2.y() + (FT(1) - t1 - t2) * c3.y(),
                              t1 * c1.z() + t2 * c2.z() + (FT(1) - t1 - t2) * c3.z());
    const FT r = t1 * r1 + t2 * r2 + (FT(1) - t1 - t2) * r3;
    return sphere_distance(p, c, r);
  }

  // Evaluate the minimum of f(t) = ||c(t) - p|| - r(t) over the medial cone defined by two spheres (c1,r1) and (c2,r2),
  // with c(t) = (1-t) * c1 + t * c2, r(t) = (1-t) * r1 + t * r2, and t in [0,1].
  FT eval_cone(const Point_3& c1, const Point_3& c2, const FT r1, const FT r2, const Point_3& p) const {
    const Vector_3 c12 = c2 - c1;
    const Vector_3 c1p = p - c1;
    const FT sq_length = c12.squared_length();
    if(sq_length == FT(0))
      return (std::min)(sphere_distance(p, c1, r1), sphere_distance(p, c2, r2));
    const FT length = CGAL::approximate_sqrt(sq_length);
    const FT k = (r2 - r1) / length;
    if(k * k >= FT(1))
      return (std::min)(sphere_distance(p, c1, r1), sphere_distance(p, c2, r2));
    const FT s0 = CGAL::scalar_product(c1p, c12) / length;
    const FT h = CGAL::approximate_sqrt(CGAL::cross_product(c1p, c12).squared_length()) / length;
    const FT z = s0 + k * h / CGAL::approximate_sqrt(FT(1) - k * k);
    const FT t = std::clamp(z / length, FT(0), FT(1));
    return cone_distance(p, c1, c2, r1, r2, t);
  }

  // Evaluate the minimum of f(t1,t2) = ||c(t1,t2) - p|| - r(t1,t2) over the medial slab defined by three spheres
  // (c1,r1), (c2,r2) and (c3,r3), with c(t1,t2) = t1 * c1 + t2 * c2 + (1-t1-t2) * c3, r(t1,t2) defined likewise,
  // and (t1,t2) in the triangle t1 >= 0, t2 >= 0, t1 + t2 <= 1.
  FT eval_slab(const Point_3& c1,
               const Point_3& c2,
               const Point_3& c3,
               const FT& r1,
               const FT& r2,
               const FT& r3,
               const Point_3& p) const {
    const FT boundary =
        (std::min)({eval_cone(c1, c2, r1, r2, p), eval_cone(c1, c3, r1, r3, p), eval_cone(c2, c3, r2, r3, p)});

    const Vector_3 e1 = c1 - c3;
    const Vector_3 e2 = c2 - c3;
    const Vector_3 c3p = c3 - p;
    const Vector_3 normal = CGAL::cross_product(e1, e2);
    const FT sq_normal = normal.squared_length(); // determinant of the Gram matrix of (e1, e2)
    if(sq_normal == FT(0))
      return boundary;

    const FT a = CGAL::scalar_product(e1, e1);
    const FT b = CGAL::scalar_product(e2, e2);
    const FT c = CGAL::scalar_product(e1, e2);
    const FT d = CGAL::scalar_product(e1, c3p);
    const FT e = CGAL::scalar_product(e2, c3p);
    const FT r13 = r1 - r3;
    const FT r23 = r2 - r3;

    // m = alpha * e1 + beta * e2 with m * e1 = r13 and m * e2 = r23
    const FT alpha = (b * r13 - c * r23) / sq_normal;
    const FT beta = (a * r23 - c * r13) / sq_normal;
    const FT sq_normal_component = FT(1) - (alpha * r13 + beta * r23); // 1 - ||m||^2
    if(sq_normal_component <= FT(0))
      return boundary;

    // signed distance from p to the plane of the centers, and ||c - p|| at the stationary point
    const FT sigma = CGAL::scalar_product(c3p, normal) / CGAL::approximate_sqrt(sq_normal);
    const FT dist = CGAL::abs(sigma) / CGAL::approximate_sqrt(sq_normal_component);

    // e_i * (c - c3) = e_i * (p - c3 + dist * n)
    const FT rhs1 = dist * r13 - d;
    const FT rhs2 = dist * r23 - e;
    const FT t1 = (b * rhs1 - c * rhs2) / sq_normal;
    const FT t2 = (a * rhs2 - c * rhs1) / sq_normal;
    if(t1 >= FT(0) && t2 >= FT(0) && t1 + t2 <= FT(1))
      return (std::min)(boundary, slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2));
    return boundary;
  }

private:
  const MSkeleton& skeleton_;
  const std::vector<Sphere_3>* spheres_;

  std::vector<FT> radii_;
};

#endif // DoXYGEN_RUNNING


} // namespace CGAL
#endif // CGAL_MEDIAL_SKELETON_OFFSET_FUNCTION_H
