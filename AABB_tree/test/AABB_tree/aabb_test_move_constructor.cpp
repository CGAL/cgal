
#include <iostream>
#include <list>
#include <utility>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/assertions.h>
#include <fstream>

template <int test_number> auto create_tree();
template <int test_number, typename T> void init_tree(T &tree) {}
template <int test_number, typename T> bool test_tree(T &tree);
template <int test_number> class TestUtils;

// test 0 is from "aabb_test_singleton_tree"
template <> struct TestUtils<0> {
  typedef CGAL::Simple_cartesian<double> K;
  typedef K::FT FT;
  typedef K::Point_3 Point;
  typedef K::Plane_3 Plane;
  typedef K::Segment_3 Segment;
  typedef K::Triangle_3 Triangle;
  typedef std::vector<Segment>::iterator Iterator;
  typedef CGAL::AABB_segment_primitive<K, Iterator> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  Point a = {1.0, 0.0, 0.0};
  Point b = {0.0, 1.0, 0.0};
  Point c = {0.0, 0.0, 1.0};
  Point d = {0.0, 0.0, 0.0};
  std::vector<Segment> segments = {Segment(Point(0, 0, 0), Point(2, 2, 2))};
};

template <> auto create_tree<0>() {
  using T = TestUtils<0>;
  T utils;
  return T::Tree(utils.segments.begin(), utils.segments.end());
}
using Test0Param = decltype(create_tree<0>());
template <> bool test_tree<0, Test0Param>(Test0Param &tree) {
  using T = TestUtils<0>;
  T utils;

  T::Plane plane_query(utils.a, utils.b, utils.d);
  T::Triangle triangle_query(utils.a, utils.b, utils.c);

  // Test calls to all functions
  CGAL::Emptyset_iterator devnull;
  tree.all_intersections(triangle_query, devnull);
  tree.all_intersected_primitives(triangle_query, devnull);
  assert(tree.any_intersected_primitive(triangle_query));
  assert(tree.any_intersection(triangle_query));
  const CGAL::Bbox_3 bbox = tree.bbox();
  assert(bbox == CGAL::Bbox_3(0, 0, 0, 2, 2, 2));
  tree.clear();
  tree.insert(utils.segments.begin(), utils.segments.end());
  tree.build();
  assert(tree.closest_point(T::Point(-0.1, -0.1, -0.1)) == T::Point(0, 0, 0));
  assert(tree.closest_point(T::Point(-0.1, -0.1, -0.1), T::Point(0, 0, 0)) ==
         T::Point(0, 0, 0));
  assert(tree.closest_point_and_primitive(T::Point(-0.1, -0.1, -0.1)).second ==
         utils.segments.begin());
  assert(tree.do_intersect(plane_query) == true);
  assert(tree.do_intersect(triangle_query) == true);
  assert(!tree.empty());
  assert(tree.size() == 1);
  tree.clear();
  assert(tree.size() == 0);
  tree.insert(utils.segments.begin(), utils.segments.end());
  assert(tree.size() == 1);
  assert(tree.number_of_intersected_primitives(plane_query) == 1);
  tree.rebuild(utils.segments.begin(), utils.segments.end());
  assert(tree.size() == 1);
  assert(tree.number_of_intersected_primitives(triangle_query) == 1);
  assert(tree.squared_distance(T::Point(0, 0, 0)) == 0);
  return true;
}

// test 1 is from "aabb_test_all_intersected_primitives"
template <> struct TestUtils<1> {
  typedef CGAL::Epick K;
  typedef K::FT FT;
  typedef K::Point_3 Point;
  typedef K::Vector_3 Vector;
  typedef K::Segment_3 Segment;
  typedef K::Ray_3 Ray;
  typedef CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick>> Mesh;
  typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh, CGAL::Default,
                                                      CGAL::Tag_false>
      S_Primitive;
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh, CGAL::Default,
                                                   CGAL::Tag_false>
      T_Primitive;
  typedef CGAL::AABB_traits<K, T_Primitive> T_Traits;
  typedef CGAL::AABB_traits<K, S_Primitive> S_Traits;
  typedef CGAL::AABB_tree<T_Traits> T_Tree;
  typedef CGAL::AABB_tree<S_Traits> S_Tree;
  typedef T_Tree::Primitive_id T_Primitive_id;
  typedef S_Tree::Primitive_id S_Primitive_id;
};

template <> auto create_tree<1>() {
  using T = TestUtils<1>;

  static CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick>> m1 = {};
  static CGAL::Surface_mesh<CGAL::Point_3<CGAL::Epick>> m2 = {};
  static bool mesh_loaded = false;
  if (!mesh_loaded) {
    std::ifstream in("data/cube.off");
    assert(in);
    in >> m1;
    in.close();
    in.open("data/tetrahedron.off");
    assert(in);
    in >> m2;
    in.close();
    mesh_loaded = true;
  }
  return std::make_pair(T::T_Tree{faces(m1).first, faces(m1).second, m1},
                        T::S_Tree{edges(m2).first, edges(m2).second, m2});
}
using Test1Param = decltype(create_tree<1>());
template <> void init_tree<1, Test1Param>(Test1Param &trees) {
  trees.first.build();
  trees.second.build();
}
template <> bool test_tree<1, Test1Param>(Test1Param &trees) {
  using T = TestUtils<1>;

  auto &cube_tree = trees.first;
  auto &tet_tree = trees.second;

  std::list<T::T_Tree::Primitive::Id> t_primitives;
  std::list<T::S_Tree::Primitive::Id> s_primitives;
  cube_tree.all_intersected_primitives(tet_tree,
                                       std::back_inserter(t_primitives));
  CGAL_assertion(t_primitives.size() == 6);
  tet_tree.all_intersected_primitives(cube_tree,
                                      std::back_inserter(s_primitives));
  CGAL_assertion(s_primitives.size() == 6);
  CGAL_assertion(tet_tree.do_intersect(cube_tree));
  CGAL_assertion(cube_tree.do_intersect(tet_tree));

  std::vector<T::T_Tree::Primitive::Id> all_primitives;
  cube_tree.all_intersected_primitives(tet_tree,
                                       std::back_inserter(all_primitives));
  bool found_f5 = false;
  for (auto prim : all_primitives) {
    if ((int)prim.first == 5)
      found_f5 = true;
  }
  CGAL_assertion(found_f5);
  CGAL_USE(found_f5);
  return true;
}

template <int test_number> bool run_test() {
  // create_tree should return prvalue for guaranteed copy elision
  auto tree_1 = create_tree<test_number>();
  init_tree<test_number>(tree_1);
  auto tree_2 = create_tree<test_number>();
  init_tree<test_number>(tree_2);
  auto tree_3 = create_tree<test_number>();
  init_tree<test_number>(tree_3);

  decltype(tree_1) tree_ctor{std::move(tree_2)};
  decltype(tree_1) tree_assig{};
  tree_assig = std::move(tree_3);

  bool normal = test_tree<test_number>(tree_1);
  bool move_ctor = test_tree<test_number>(tree_ctor);
  bool move_ass =
      test_tree<test_number>(tree_assig); // test move assignment operator

  if (!normal)
    std::cout << "Test " << test_number << "failed on the original tree\n";
  if (!move_ctor)
    std::cout << "Test " << test_number << "failed on move constructed tree\n";
  if (!move_ass)
    std::cout << "Test " << test_number << "failed on move assigned tree\n";
  return normal && move_ctor && move_ass;
}

int main() { return (run_test<0>() && run_test<1>()) ? 0 : 1; }
