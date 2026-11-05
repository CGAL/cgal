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
#include <CGAL/AABB_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>
#include <boost/functional/hash.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <algorithm>
#include <initializer_list>
#include <memory>
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
  using Iso_cuboid_3 = typename GT::Iso_cuboid_3;
  using MSkeleton = Medial_skeleton<TriangleMesh, GeomTraits>;

  // The primitives of the tree are the medial slabs (ids 0 to the number of faces - 1, in the order of `faces()`)
  // followed by the cones of the free edges (medial edges that do not belong to any slabs).
  // The datum of a primitive is its bounding box.
  using Box_map = typename Pointer_property_map<Iso_cuboid_3>::const_type;
  using Point_map = typename Pointer_property_map<Point_3>::const_type;
  using Primitive = AABB_primitive<std::size_t, Box_map, Point_map, Tag_true, Tag_false>;
  using AABB_traits = AABB_traits_3<GT, Primitive>;
  using Tree = AABB_tree<AABB_traits>;

  // The tree is not copyable, and the function object is copied by the meshing domain:
  // the copies share the tree and the data its property maps point to.
  struct Search_tree
  {
    std::vector<Iso_cuboid_3> boxes;
    std::vector<Point_3> reference_points;
    std::vector<FT> max_radii; // largest radius of the spheres of each primitive
    FT max_radius = FT(0);     // largest radius of all spheres
    Tree tree;
  };

  class Minimum_traversal_traits;


public:
  Medial_skeleton_offset_function(const MSkeleton& skeleton)
      : skeleton_(skeleton)
      , spheres_(&skeleton_.vertices()) {
    const std::size_t n = skeleton_.number_of_vertices();
    radii_.reserve(n);
    for(const auto& sph : *spheres_)
      radii_.push_back(CGAL::approximate_sqrt(sph.squared_radius()));

    // The cone of an edge of a face is contained in the slab of that face,
    // so only the edges that do not belong to any face are evaluated as cones.
    using Edge = std::pair<std::size_t, std::size_t>;
    auto make_edge = [](std::size_t i, std::size_t j) { return Edge((std::min)(i, j), (std::max)(i, j)); };
    std::unordered_set<Edge, boost::hash<Edge>> face_edges;
    for(const auto& tri : skeleton_.faces()) {
      face_edges.insert(make_edge(tri[0], tri[1]));
      face_edges.insert(make_edge(tri[1], tri[2]));
      face_edges.insert(make_edge(tri[2], tri[0]));
    }
    for(const auto& e : skeleton_.edges())
      if(face_edges.find(make_edge(e.first, e.second)) == face_edges.end())
        free_edges_.push_back(e);

    build_search_tree();
  }


  // Evaluate the minimum distance from point p to the interpolated medial axis (cones and slabs).
  // The primitives whose lower bound is not smaller than the current minimum are skipped,
  // which gives the same value as evaluating all primitives, up to rounding errors.
  FT operator()(const Point_3& p) const {
    if(search_tree_ == nullptr)
      return (std::numeric_limits<FT>::max)();
    Minimum_traversal_traits traits(*this);
    search_tree_->tree.traversal_with_priority(p, traits);
    return traits.minimum();
  }

private:
  void build_search_tree() {
    const std::size_t n = skeleton_.faces().size() + free_edges_.size();
    if(n == 0)
      return;

    auto data = std::make_shared<Search_tree>();
    data->boxes.reserve(n);
    data->reference_points.reserve(n);
    data->max_radii.reserve(n);
    // `Sphere_3::bbox()` rounds outward, so the box of a primitive contains all its spheres.
    auto add_primitive = [&](std::initializer_list<std::size_t> ids) {
      Bbox_3 box;
      FT max_radius = FT(0);
      for(std::size_t i : ids) {
        box += (*spheres_)[i].bbox();
        max_radius = (std::max)(max_radius, radii_[i]);
      }
      data->boxes.emplace_back(box);
      data->reference_points.push_back(CGAL::midpoint(data->boxes.back().min(), data->boxes.back().max()));
      data->max_radii.push_back(max_radius);
      data->max_radius = (std::max)(data->max_radius, max_radius);
    };
    for(const auto& tri : skeleton_.faces())
      add_primitive({tri[0], tri[1], tri[2]});
    for(const auto& [i, j] : free_edges_)
      add_primitive({i, j});

    data->tree.rebuild(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(n),
                       make_property_map(std::as_const(data->boxes)),
                       make_property_map(std::as_const(data->reference_points)));
    search_tree_ = std::move(data);
  }

  // Lower bound of the value at p of a primitive whose spheres are in `box` and have a radius at most `max_radius`:
  // the distance from p to the box if p is outside the box, and -max_radius otherwise.
  template <class Box>
  static FT lower_bound(const Point_3& p, const Box& box, const FT& max_radius) {
    FT sq_dist = FT(0);
    bool inside = true;
    for(int i = 0; i < 3; ++i) {
      const FT v = p.cartesian(i);
      if(v < FT(box.min_coord(i))) {
        sq_dist += CGAL::square(FT(box.min_coord(i)) - v);
        inside = false;
      } else if(v > FT(box.max_coord(i))) {
        sq_dist += CGAL::square(v - FT(box.max_coord(i)));
        inside = false;
      }
    }
    return inside ? -max_radius : CGAL::approximate_sqrt(sq_dist);
  }

  FT eval_primitive(std::size_t id, const Point_3& p) const {
    const std::size_t nf = skeleton_.faces().size();
    if(id < nf) {
      const auto& tri = skeleton_.faces()[id];
      return eval_slab((*spheres_)[tri[0]].center(), (*spheres_)[tri[1]].center(), (*spheres_)[tri[2]].center(),
                       radii_[tri[0]], radii_[tri[1]], radii_[tri[2]], p);
    }
    const auto& [i, j] = free_edges_[id - nf];
    return eval_cone((*spheres_)[i].center(), (*spheres_)[j].center(), radii_[i], radii_[j], p);
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

  // Traversal traits computing the minimum of the values of the primitives at the query point.
  // A node or a primitive is visited only if its lower bound is smaller than the current minimum,
  // and the child with the smaller lower bound is visited first.
  class Minimum_traversal_traits
  {
  public:
    using Priority = FT;

    explicit Minimum_traversal_traits(const Medial_skeleton_offset_function& f)
        : f_(f)
        , data_(*f.search_tree_)
        , minimum_((std::numeric_limits<FT>::max)()) {}

    bool go_further() const { return true; }

    void intersection(const Point_3& p, const Primitive& primitive) {
      const std::size_t id = primitive.id();
      if(lower_bound(p, data_.boxes[id], data_.max_radii[id]) < minimum_)
        minimum_ = (std::min)(minimum_, f_.eval_primitive(id, p));
    }

    template <class Node>
    std::pair<bool, Priority> do_intersect_with_priority(const Point_3& p, const Node& node) const {
      const FT lb = lower_bound(p, node.bbox(), data_.max_radius);
      return std::make_pair(lb < minimum_, -lb);
    }

    template <class Node>
    bool do_intersect(const Point_3& p, const Node& node) const {
      return do_intersect_with_priority(p, node).first;
    }

    FT minimum() const { return minimum_; }

  private:
    const Medial_skeleton_offset_function& f_;
    const Search_tree& data_;
    FT minimum_;
  };

private:
  const MSkeleton& skeleton_;
  const std::vector<Sphere_3>* spheres_;

  std::vector<FT> radii_;
  std::vector<std::pair<std::size_t, std::size_t>> free_edges_;
  std::shared_ptr<const Search_tree> search_tree_;
};

#endif // DoXYGEN_RUNNING


} // namespace CGAL
#endif // CGAL_MEDIAL_SKELETON_OFFSET_FUNCTION_H
