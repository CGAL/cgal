#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/internal/Octree_wrapper.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <cmath>
#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Vector = typename Kernel::Vector_3;
using Point = typename Kernel::Point_3;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

using Octree_wrapper = CGAL::Isosurfacing::internal::Octree_wrapper<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Octree_wrapper>;
using Gradients = CGAL::Isosurfacing::Gradient_function_3<Octree_wrapper>;
using Domain = CGAL::Isosurfacing::Isosurfacing_domain_3<Octree_wrapper, Data>;

struct Refine_one_eighth
{
  std::size_t min_depth_;
  std::size_t max_depth_;

  std::size_t octree_dim_;

  Refine_one_eighth(std::size_t min_depth,
                    std::size_t max_depth)
    : min_depth_(min_depth),
      max_depth_(max_depth)
  {
    octree_dim_ = std::size_t(1) << max_depth_;
  }

  Octree_wrapper::Uniform_coords uniform_coordinates(const Octree_wrapper::Octree::Node& node) const
  {
    auto coords = node.global_coordinates();
    const std::size_t depth_factor = std::size_t(1) << (max_depth_ - node.depth());
    for(int i=0; i < Octree_wrapper::Octree::Node::Dimension::value; ++i)
      coords[i] *= uint32_t(depth_factor);

    return coords;
  }

  bool operator()(const Octree_wrapper::Octree::Node& n) const
  {
    if(n.depth() < min_depth_)
      return true;

    if(n.depth() == max_depth_)
      return false;

    auto leaf_coords = uniform_coordinates(n);

    if(leaf_coords[0] >= octree_dim_ / 2)
      return false;

    if(leaf_coords[1] >= octree_dim_ / 2)
      return false;

    if(leaf_coords[2] >= octree_dim_ / 2)
      return false;

    return true;
  }
};

int main(int, char**)
{
  const CGAL::Bbox_3 bbox{-1., -1., -1.,  1., 1., 1.};
  Octree_wrapper octree_wrap { bbox };

  Refine_one_eighth split_predicate(3, 5);
  octree_wrap.refine(split_predicate);

  auto sphere_function = [](const Point& p) -> FT
  {
    return std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
  };

  auto sphere_gradient = [&](const Point& p) -> Vector
  {
    const Vector g = p - CGAL::ORIGIN;
    return g / std::sqrt(g.squared_length());
  };

  Data data(sphere_function, sphere_gradient);
  Domain domain(octree_wrap, data);

  Point_range points;
  Polygon_range polygons;

  CGAL::Isosurfacing::dual_contouring(domain, 0.8, points, polygons);

  CGAL::IO::write_polygon_soup("dual_contouring_octree.off", points, polygons);

  return EXIT_SUCCESS;
}
