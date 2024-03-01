
#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Orthtree_traits.h>
#include <CGAL/Orthtree/Split_predicates.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Simple_cartesian.h>

#include <iostream>
#include <cassert>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using FT = Kernel::FT;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;
using Octree_without_data = CGAL::Orthtree<CGAL::Orthtree_traits<Kernel, 3>>;

template<typename Tree>
int test(Tree &tree)
{
  assert (tree.is_leaf(tree.root())); // tree is not refined yet

  Tree copy1 (tree);
  assert (copy1.is_leaf(copy1.root())); // copy1 is thus not refined either
  assert (tree == copy1); // tree should be equal to copy1

  tree.refine(CGAL::Orthtrees::Maximum_depth(5));
  assert (!tree.is_leaf(tree.root())); // tree is now refined
  assert (copy1.is_leaf(copy1.root())); // copy1 should be unaffected and still unrefined
  assert (tree != copy1); // tree should be different from copy1

  Tree copy2 (tree);
  assert (!copy2.is_leaf(copy2.root())); // copy2 should be refined
  assert (tree == copy2); // tree should be equal to copy2

  Tree move (std::move(tree));
  assert (!move.is_leaf(move.root())); // move should be refined

  assert (tree.is_leaf(tree.root())); // tree should be back to init state (unrefined)
  assert (copy1.is_leaf(copy1.root())); // copy1 still unaffected and still unrefined
  assert (!copy2.is_leaf(copy2.root())); // copy2 unaffected by move and still refined
  assert (move == copy2); // move should be equal to copy2

  return EXIT_SUCCESS;
}

int main()
{
  std::size_t nb_pts = 100;
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++i)
    points.insert(*(generator++));

  Octree base(points, points.point_map());
  test(base);

  Octree_without_data base2({});
  test(base2);
}
