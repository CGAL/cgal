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
template <class TriangleMesh, class GeomTraits = Default>
class Medial_skeleton_offset_function
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
      , spheres_(&skeleton_.vertices())
      {
    const std::size_t n = skeleton_.number_of_vertices();
    radii_.reserve(n);
    for(const auto& sph : *spheres_)
      radii_.push_back(CGAL::approximate_sqrt(sph.squared_radius()));
  }


  // Evaluate the minimum distance from point p to the interpolated medial axis (cones and slabs)
  FT operator()(const Point_3& p) const {
    FT dmin = (std::numeric_limits<FT>::max)();

    //TODO: the filtering strategy doesn't work
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
      const FT& r_max = r1>r2 ? (r1 > r3 ? r1 : r3) : (r2 > r3 ? r2 : r3);
      /* if(dmin > eval_slab(c1, c2, c3, r_max, r_max, r_max,p)) {
         dmin = (std::min)(dmin, eval_slab(c1, c2, c3, r1, r2, r3, p));
       }*/
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
      const FT r_max = r1 > r2 ? r1 : r2;
      /*if(dmin > eval_cone(c1, c2, r_max, r_max, p))
        dmin = (std::min)(dmin, eval_cone(c1, c2, r1, r2, p));*/
      dmin = (std::min)(dmin, eval_cone(c1, c2, r1, r2, p));
      if(dmin < 0)
        return dmin;
    }

    return dmin;
  }

private:
  int solve_quadric(const FT a, const FT b, const FT c, FT& t1, FT& t2) const {
    if(a == FT(0)) {
      if(b == FT(0))
        return 0;
      t1 = -c / b;
      return 1;
    } else {
      FT delta = b * b - 4 * a * c;
      if(delta < FT(0))
        return 0;
      else if(delta == FT(0)) {
        t1 = -b / (FT(2) * a);
        return 1;
      } else {
        const FT sqrt_delta = CGAL::approximate_sqrt(delta);
        t1 = (-b - sqrt_delta) / (FT(2) * a);
        t2 = (-b + sqrt_delta) / (FT(2) * a);
        return 2;
      }
    }
  }

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

  // Evaluate the minimum distrance from point p to the medial cone defined by two spheres (c1,r1) and (c2,r2)
  FT eval_cone(const Point_3& c1, const Point_3& c2, const FT r1, const FT r2, const Point_3& p) const {

    FT dmin = (std::numeric_limits<FT>::max)();
    const Vector_3 c21 = c2 - c1;
    const Vector_3 c1p = p - c1;
    const FT r21 = r2 - r1;
    if(CGAL::abs(r21) < FT(1e-12)) {
      FT t = -(c21 * c1p) / c21.squared_length();
      t = std::clamp(t, FT(0), FT(1));
      return cone_distance(p, c1, c2, r1, r2, t);
    }
    const FT a = CGAL::scalar_product(c21, c21);
    const FT b = CGAL::scalar_product(c21, c1p);
    const FT c = CGAL::scalar_product(c1p, c1p);
    const FT A = a * (a - r21 * r21);
    const FT B = 2 * b * (r21 * r21 - a);
    const FT C = b * b - r21 * r21 * c;
    FT t1 = 0, t2 = 0, dist1, dist2;
    int root_nb = solve_quadric(A, B, C, t1, t2);
    if(root_nb != 0) {
      t1 = std::clamp(t1, FT(0), FT(1));
      dist1 = cone_distance(p, c1, c2, r1, r2, t1);
      dmin = (std::min)(dmin, dist1);
      if(root_nb != 1) {
        t2 = std::clamp(t2, FT(0), FT(1));
        dist2 = cone_distance(p, c1, c2, r1, r2, t2);
        dmin = (std::min)(dmin, dist2);
      }
    }
    dmin = (std::min)(dmin, cone_distance(p, c1, c2, r1, r2, FT(0)));
    dmin = (std::min)(dmin, cone_distance(p, c1, c2, r1, r2, FT(1)));
    return dmin;
  }

  // Evaluate the minimum distance from point p to the medial slab defined by three spheres (c1,r1), (c2,r2) and (c3,r3)
  FT eval_slab(const Point_3& c1,
               const Point_3& c2,
               const Point_3& c3,
               const FT& r1,
               const FT& r2,
               const FT& r3,
               const Point_3& p) const {
    constexpr FT EPS = FT(1e-12);
    auto is_zero = [&](FT x) { return CGAL::abs(x) <= EPS; };
    auto feasible = [&](FT u, FT v) { return u >= FT(0) && v >= FT(0) && u + v <= FT(1) && u <= FT(1) && v <= FT(1); };
    FT dmin = (std::numeric_limits<FT>::max)();

    // c(t1, t2) = t1 * c1 + t2 * c2+ (1-t1-t2) * c3 = (c1-c3) * t1 + (c2-c3) * t2 + c3
    // r(t1, t2) = t1 * r1 + t2 * r2 + (1-t1-t2) * r3 = (r1-r3) * t1 + (r2-r3) * t2 + r3
    // f(t1, t2) = ||c(t1, t2) - p || - r(t1, t2)
    // obj: argmin_(t1,t2) f(t1, t2)
    const Vector_3 c13 = c1 - c3;
    const Vector_3 c23 = c2 - c3;
    const Vector_3 c3p = c3 - p;
    const FT r13 = r1 - r3;
    const FT r23 = r2 - r3;
    // let x = c(t1, t2) - p = c13 * t1 + c23 * t2 + c3p
    //  f(t1, t2) = ||x|| - (r1-r3) * t1 + (r2-r3) * t2 + r3
    // df(t1,t2) / dt1 = (x*c13)/||x|| - r13 = 0 ===> x*c13 = ||x||* r13
    // df(t1,t2) / dt2 = (x*c23)/||x|| - r23 = 0 ===> x*c23 = ||x||* r23
    // x^2 = ||c13||^2 *t1^2 + ||c23||^2 *t2^2 + 2*(c13*c23)*t1*t2 + 2*(c13*c3p)*t1 + 2*(c23*c3p)*t2 + ||c3p||^2
    const FT a = CGAL::scalar_product(c13, c13);
    const FT b = CGAL::scalar_product(c23, c23);
    const FT c = CGAL::scalar_product(c13, c23);
    const FT d = CGAL::scalar_product(c13, c3p);
    const FT e = CGAL::scalar_product(c23, c3p);
    const FT f = CGAL::scalar_product(c3p, c3p);
    // x^2 = a*t1^2 + b*t2^2 + 2*c*t1*t2 + 2*d*t1 + 2*e*t2 + f
    FT t1 = 0, t2 = 0, dist1 = 0, dist2 = 0;

    // three spheres have the same radius
    if(is_zero(r13) && is_zero(r23)) {
      //     x*c13                                         = 0
      //===> ||c13||^2 * t1 + (c13 *c23)* t2 + (c13 * c3p) = 0
      //===> a*t1 + c*t2 + d                               = 0

      //     x*c23                                         = 0
      //===> ||c23||^2 * t2 + (c23 *c13)* t1 + (c23 * c3p) = 0
      //===> b*t2 + c*t1 + e                               = 0

      const FT denom = a * b - c * c;
      if(!is_zero(denom)) {
        t1 = (c * e - b * d) / denom;
        t2 = (c * d - a * e) / denom;
      }
    } else if(is_zero(r13) && !is_zero(r23)) {
      // x * c13 = 0 ===> a*t1 + c*t2 + d = 0 ===> t1 = -(c/a)*t2 - d/a
      const FT h = -c / a;
      const FT k = -d / a;

      //         df(t1,t2) / dt2 = 0
      //===>               x*c23 = ||x||* r23
      //===>         (x * c23)^2 = ||x||^2 * r23^2
      //===> (b*t2 + c*t1 + e)^2 = ||x||^2 * r23^2
      const FT A = (b + c * h) * (b + c * h) - r23 * r23 * (a * h * h + b + FT(2) * c * h);
      const FT B = (FT(2) * (b + c * h) * (c * k + e) -
                    r23 * r23 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
      const FT C = (c * k + e) * (c * k + e) - r23 * r23 * (a * k * k + FT(2) * d * k + f);
      FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
      solve_quadric(A, B, C, t2_1, t2_2);
      t1_1 = h * t2_1 + k;
      t1_2 = h * t2_2 + k;
      dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
      dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
      if(dist2 < dist1) {
        t1 = t1_2;
        t2 = t2_2;
      } else {
        t1 = t1_1;
        t2 = t2_1;
      }
    } else if(!is_zero(r13) && is_zero(r23)) {
      // x * c23 = 0 ===> b*t2 + c*t1 + e = 0 ===> t1 = -(b/c)*t2 - e/c
      const FT h = -b / c;
      const FT k = -e / c;

      //         df(t1,t2) / dt1 = 0
      //===>               x*c13 = ||x||* r13
      //===>         (x * c13)^2 = ||x||^2 * r13^2
      //===> (a*t1 + c*t2 + d)^2 = ||x||^2 * r13^2
      const FT A = (c + a * h) * (c + a * h) - r13 * r13 * (a * h * h + b + FT(2) * c * h);
      const FT B =
          (2 * (c + a * h) * (a * k + d) - r13 * r13 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
      const FT C = (a * k + d) * (a * k + d) - r13 * r13 * (a * k * k + FT(2) * d * k + f);
      FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
      solve_quadric(A, B, C, t2_1, t2_2);
      t1_1 = h * t2_1 + k;
      t1_2 = h * t2_2 + k;
      dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
      dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
      if(dist2 < dist1) {
        t1 = t1_2;
        t2 = t2_2;
      } else {
        t1 = t1_1;
        t2 = t2_1;
      }
    } else {
      // both spheres have different radii
      // x * c13 = ||x|| * r13
      // x * c23 = ||x|| * r23
      // ===> (x* c13)/r13   = (x * c23)/r23
      // ===> r23 * (x* c13) = r13 * (x * c23)
      // ===> (r23 * a- r13* c)*t1 + (r23*c-r13*b) *t2- (r23 * d - r13 * e) = 0
      const FT u = r23 * a - r13 * c;
      const FT v = r23 * c - r13 * b;
      const FT w = r23 * d - r13 * e;
      // ===> u*t1 + v*t2 + w = 0
      if(is_zero(u) && !is_zero(v)) {
        t2 = -w / v;
        const FT A = a * a - r13 * r13 * a;
        const FT B = FT(2) * a * (c * t2 + d) - r13 * r13 * (FT(2) * c * t2 + FT(2) * d);
        const FT C = (c * t2 + d) * (c * t2 + d) - r13 * r13 * (b * t2 * t2 + FT(2) * e * t2 + f);
        FT t1_1 = 0, t1_2 = 0;
        solve_quadric(A, B, C, t1_1, t1_2);
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2);
        if(dist2 < dist1) {
          t1 = t1_2;
        } else {
          t1 = t1_1;
        }
      } else if(u != 0 && v == 0) {
        t1 = -w / u;
        const FT A = b * b - r23 * r23 * b;
        const FT B = FT(2) * b * (c * t1 + e) - r23 * r23 * (FT(2) * c * t1 + FT(2) * e);
        const FT C = (c * t1 + e) * (c * t1 + e) - r23 * r23 * (a * t1 * t1 + FT(2) * d * t1 + f);
        FT t2_1 = 0, t2_2 = 0;
        solve_quadric(A, B, C, t2_1, t2_2);
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2_1);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2_2);
        if(dist2 < dist1) {
          t2 = t2_2;
        } else {
          t2 = t2_1;
        }
      } else {
        const FT h = -v / u;
        const FT k = -w / u;
        const FT A = (b + c * h) * (b + c * h) - r23 * r23 * (a * h * h + b + FT(2) * c * h);
        const FT B = (2 * (b + c * h) * (c * k + e) -
                      r23 * r23 * (FT(2) * a * h * k + FT(2) * c * k + FT(2) * d * h + FT(2) * e));
        const FT C = (c * k + e) * (c * k + e) - r23 * r23 * (a * k * k + FT(2) * d * k + f);
        FT t1_1 = 0, t1_2 = 0, t2_1 = 0, t2_2 = 0;
        solve_quadric(A, B, C, t2_1, t2_2);
        t1_1 = h * t2_1 + k;
        t1_2 = h * t2_2 + k;
        dist1 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_1, t2_1);
        dist2 = slab_distance(p, c1, c2, c3, r1, r2, r3, t1_2, t2_2);
        if(dist2 < dist1) {
          t1 = t1_2;
          t2 = t2_2;
        } else {
          t1 = t1_1;
          t2 = t2_1;
        }
      }
    }
    if(feasible(t1, t2)) {
      dmin = (std::min)(dmin, slab_distance(p, c1, c2, c3, r1, r2, r3, t1, t2));
      return dmin;
    }
    dmin = (std::min)(dmin, eval_cone(c1, c2, r1, r2, p));
    dmin = (std::min)(dmin, eval_cone(c1, c3, r1, r3, p));
    dmin = (std::min)(dmin, eval_cone(c2, c3, r2, r3, p));
    return dmin;
  }

private:
  const MSkeleton& skeleton_;
  const std::vector<Sphere_3>* spheres_;

  std::vector<FT> radii_;
};

#endif // DoXYGEN_RUNNING


} // namespace CGAL
#endif // CGAL_MEDIAL_SKELETON_OFFSET_FUNCTION_H
