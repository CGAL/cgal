
#define CGAL_TRACE_STREAM std::cerr

#include <iostream>
#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>
#include <CGAL/point_generators_3.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using FT = Kernel::FT;
using Point_set = CGAL::Point_set_3<Point>;
using Octree = CGAL::Octree<Kernel, Point_set, typename Point_set::Point_map>;

int main(void)
{
  std::size_t nb_pts = 100;
  Point_set points;
  CGAL::Random_points_in_cube_3<Point> generator;
  points.reserve(nb_pts);
  for (std::size_t i = 0; i < nb_pts; ++i)
    points.insert(*(generator++));

  Octree base ({points, points.point_map()});
  assert (base.is_leaf(base.root())); // base is not refined yet

  Octree copy1 (base);
  assert (copy1.is_leaf(copy1.root())); // copy1 is thus not refined either
  assert (base == copy1); // base should be equal to copy1

  base.refine();
  assert (!base.is_leaf(base.root())); // base is now refined
  assert (copy1.is_leaf(copy1.root())); // copy1 should be unaffected and still unrefined
  assert (base != copy1); // base should be different from copy1

  Octree copy2 (base);
  assert (!copy2.is_leaf(copy2.root())); // copy2 should be refined
  assert (base == copy2); // base should be equal to copy2

  Octree move (std::move(base));
  assert (!move.is_leaf(move.root())); // move should be refined
  // fixme: my linter isn't happy about use-after-move
  assert (base.is_leaf(base.root())); // base should be back to init state (unrefined)
  assert (copy1.is_leaf(copy1.root())); // copy1 still unaffected and still unrefined
  assert (!copy2.is_leaf(copy2.root())); // copy2 unaffected by move and still refined
  assert (move == copy2); // move should be equal to copy2

  return EXIT_SUCCESS;
}
