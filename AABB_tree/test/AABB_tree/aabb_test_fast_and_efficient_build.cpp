#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Triangle_3 = Kernel::Triangle_3;
using Segment_3 = Kernel::Segment_3;

using Triangles = std::vector<Triangle_3>;
using Primitive_iterator = Triangles::const_iterator;
using Primitive = CGAL::AABB_triangle_primitive<Kernel, Primitive_iterator>;
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

static double height(int i, int j)
{
  return 0.05 * ((i * 7 + j * 11) % 9);
}

static Triangles make_triangles(const int nx, const int ny)
{
  std::vector<Point_3> points;
  points.reserve((nx + 1) * (ny + 1));

  for(int j = 0; j <= ny; ++j)
    for(int i = 0; i <= nx; ++i)
      points.emplace_back(i, j, height(i, j));

  auto point = [&](int i, int j) -> const Point_3&
  {
    return points[static_cast<std::size_t>(j) * (nx + 1) + i];
  };

  Triangles triangles;
  triangles.reserve(2 * nx * ny);

  for(int j = 0; j < ny; ++j)
  {
    for(int i = 0; i < nx; ++i)
    {
      const Point_3& p00 = point(i,     j);
      const Point_3& p10 = point(i + 1, j);
      const Point_3& p01 = point(i,     j + 1);
      const Point_3& p11 = point(i + 1, j + 1);

      triangles.emplace_back(p00, p10, p11);
      triangles.emplace_back(p00, p11, p01);
    }
  }

  return triangles;
}

static std::vector<Segment_3> make_queries(const int nx, const int ny)
{
  std::vector<Segment_3> queries;

  // Vertical segments through cell interiors and grid boundaries.
  for(int j = 0; j <= ny; ++j)
  {
    for(int i = 0; i <= nx; ++i)
    {
      const double x = (i == nx) ? i - 0.25 : i + 0.25;
      const double y = (j == ny) ? j - 0.25 : j + 0.25;
      queries.emplace_back(Point_3(x, y, -1.0), Point_3(x, y, 2.0));
    }
  }

  // Queries that should miss the triangle set.
  queries.emplace_back(Point_3(-2, -2, -1), Point_3(-2, -2, 2));
  queries.emplace_back(Point_3(nx + 2, ny + 2, -1),
                       Point_3(nx + 2, ny + 2, 2));

  // A few oblique segments crossing the input.
  queries.emplace_back(Point_3(-1, 0.5, 0.2),
                       Point_3(nx + 1, ny - 0.5, 0.2));
  queries.emplace_back(Point_3(0.5, -1, 0.3),
                       Point_3(nx - 0.5, ny + 1, 0.3));

  return queries;
}

static void same_query_results(const Tree& reference,
                               const Tree& candidate,
                               const std::vector<Segment_3>& queries,
                               const char* candidate_name)
{
  for(std::size_t i = 0; i < queries.size(); ++i)
  {
    std::vector<Primitive_iterator> reference_ids;
    std::vector<Primitive_iterator> candidate_ids;

    reference.all_intersected_primitives( queries[i], std::back_inserter(reference_ids));
    candidate.all_intersected_primitives( queries[i], std::back_inserter(candidate_ids));

    // Query traversal order is not part of the result: compare sorted IDs.
    std::sort(reference_ids.begin(), reference_ids.end());
    std::sort(candidate_ids.begin(), candidate_ids.end());

    assert(reference_ids == candidate_ids);
  }
}

int main()
{
  constexpr int nx = 30;
  constexpr int ny = 25;

  const Triangles triangles = make_triangles(nx, ny);
  const std::vector<Segment_3> queries = make_queries(nx, ny);

  Tree default_tree;
  Tree fast_tree;
  Tree efficient_tree;

  // All trees receive exactly the same primitive range.
  default_tree.insert(triangles.cbegin(), triangles.cend());
  fast_tree.insert(triangles.cbegin(), triangles.cend());
  efficient_tree.insert(triangles.cbegin(), triangles.cend());

  default_tree.build();
  fast_tree.fast_build();
  efficient_tree.efficient_build();

  same_query_results(default_tree, fast_tree, queries, "fast_build");
  same_query_results(default_tree, efficient_tree, queries, "efficient_build");

  return EXIT_SUCCESS;
}