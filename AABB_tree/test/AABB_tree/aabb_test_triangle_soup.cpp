#include <iostream>
#include <vector>
#include <array>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_indexed_triangle_primitive_3.h>

#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/IO/polygon_soup_io.h>

using K = CGAL::Simple_cartesian<double>;
using P = K::Point_3;
using L = K::Line_3;
using T = K::Triangle_3;
using R = K::Ray_3;

using PointRange = std::vector<P>;
using FaceRange = std::vector<std::array<std::size_t, 3> >;
using Primitive = CGAL::AABB_indexed_triangle_primitive_3<K, PointRange, FaceRange>;
using AABB_triangle_traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<AABB_triangle_traits>;
using Point_and_primitive_id = Tree::Point_and_primitive_id;

int main()
{
  P a(0.0, 0.0, 0.0);
  P b(0.0, 1.0, 0.0);
  P c(1.0, 0.0, 0.0);
  P d(1.0, 1.0, 0.0);
  P e(2.0, 0.0, 0.0);
  P f(2.0, 1.0, 0.0);

  PointRange points = { a, b, c, d, e, f };

  FaceRange triangles;
  triangles.push_back({ 0, 2, 1 });
  triangles.push_back({ 1, 2, 3 });
  triangles.push_back({ 3, 2, 4 });
  triangles.push_back({ 3, 4, 5 });

  // constructs AABB tree
  std::vector<std::size_t> indices(triangles.size(), 0);
  Tree tree(indices.begin(), indices.end(), points, triangles);

  // point sampling
  Point_and_primitive_id id;
  id = tree.closest_point_and_primitive(P(0.5, 0.4, 0));
  id = tree.closest_point_and_primitive(P(0.5, 0.6, 0));
  id = tree.closest_point_and_primitive(P(1.5, 0.4, 0));
  id = tree.closest_point_and_primitive(P(1.5, 0.6, 0));
  id = tree.closest_point_and_primitive(P(3.0, 0.5, 0));

  R ray(P(5.5, 0.5, 0), P(1.5, 0.4, 0));
  auto intersection = tree.first_intersection(ray);

  assert(intersection.has_value());
  assert(intersection->second == 3);

  std::vector<P> pts1, pts2;
  std::vector<std::array<std::size_t, 3> > trs1, trs2;
  if(!CGAL::IO::read_polygon_soup(CGAL::data_file_path("meshes/knot1.off"), pts1, trs1)){
    std::cout << "error reading knot1" << std::endl;
    exit(1);
  }
  if(!CGAL::IO::read_polygon_soup(CGAL::data_file_path("meshes/lion.off"), pts2, trs2)){
    std::cout << "error reading lion" << std::endl;
    exit(1);
  }

  Tree tree1(trs1.begin(), trs1.end(), pts1, trs1);
  Tree tree2(trs2.begin(), trs2.end(), pts2, trs2);
  tree1.build();
  tree2.build();

  std::vector< std::pair<std::size_t, std::size_t> > inter;
  CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(inter));
  assert(inter.size() == 1191);

  return EXIT_SUCCESS;
}
