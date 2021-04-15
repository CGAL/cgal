#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;

using Traits = CGAL::Search_traits_3<Kernel>;
using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Tree = Neighbor_search::Tree;
using Point_with_distance = Neighbor_search::Point_with_transformed_distance;

using Generator = CGAL::Random_points_in_sphere_3<Point_3>;

int main()
{
  const unsigned int N = 1000;
  const unsigned int k = 6;

  // Generate N points in a sphere
  std::vector<Point_3> points;
  points.reserve (N);
  Generator generator;
  for (unsigned int i = 0; i < N; ++ i)
    points.push_back (*(generator++));

  // Build tree in parallel
  Tree tree(points.begin(), points.end());
  tree.build<CGAL::Parallel_tag>();

  // Query tree in parallel
  std::vector<std::vector<Point_3> > neighbors (points.size());
  tbb::parallel_for (tbb::blocked_range<std::size_t> (0, points.size()),
                     [&](const tbb::blocked_range<std::size_t>& r)
                     {
                       for (std::size_t s = r.begin(); s != r.end(); ++ s)
                       {
                         // Neighbor search can be instantiated from
                         // several threads at the same time
                         Neighbor_search search (tree, points[s], k);
                         neighbors[s].reserve(k);

                         // neighbor search returns a set of pair of
                         // point and distance <Point_3,FT>, here we
                         // keep the points only
                         for (const Point_with_distance& pwd : search)
                           neighbors[s].push_back (pwd.first);
                       }
                     });

  return 0;
}
