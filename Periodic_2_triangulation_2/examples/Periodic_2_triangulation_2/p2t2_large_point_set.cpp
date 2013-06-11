#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_traits_2<K> GT;

typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> PDT;

typedef PDT::Point          Point;

int main()
{
  CGAL::Timer t;
  typedef CGAL::Creator_uniform_2<double, Point> Creator;
  CGAL::Random random(7);
  CGAL::Random_points_in_square_2<Point, Creator> in_square(.5, random);

  int n = 10000;
  std::vector<Point> pts;

  PDT PT1, PT2, PT3;

  // Generating n random points
  for (int i = 0 ; i < n ; i++)
    {
      Point p = *in_square;
      in_square++;
      pts.push_back(Point(p.x() + .5, p.y() + .5));
    }

  // Standard insertion
  t.start();
  for (int i = 0 ; i < n ; i++)
    {
      PT1.insert(pts[i]);
    }
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (Standard insertion)" << std::endl;
  t.reset();

  // Iterator range insertion using spatial sorting but no dummy points
  t.start();
  PT2.insert(pts.begin(), pts.end()); // third parameter defaults to false
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (with spatial sorting)" << std::endl;
  t.reset();

  // Iterator range insertion using spatial sorting and dummy point heuristic
  t.start();
  PT3.insert(pts.begin(), pts.end(), true);
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (Dummy point heuristic)" << std::endl;

  return 0;
}
