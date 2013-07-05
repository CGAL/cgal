// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#define CGAL_USE_DELAUNAY

#include "./types.h"
#include <CGAL/point_generators_2.h>


void test_insertion()
{
  // Create point sets
  typedef CGAL::Creator_uniform_2<double, Point>  Creator;
  CGAL::Random rnd(7);
  CGAL::Random_points_on_circle_2<Point, Creator> on_circle(0.5, rnd);

  // Center of the circle at (0.5, 0.5)
  std::vector<Point> pts_rnd1000;
  for (int i = 0 ; i < 1000 ; i++)
    {
      pts_rnd1000.push_back(*on_circle++ + Vector(0.5, 0.5));
    }

  Triangulation pt(pts_rnd1000.begin(), pts_rnd1000.end());
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.number_of_vertices() == 1000);
  pt.clear();
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.empty());
  CGAL_assertion(pt.number_of_vertices() == 0);

  pt.insert(pts_rnd1000.begin(), pts_rnd1000.end());
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.number_of_vertices() == 1000);
  pt.clear();
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.empty());
  CGAL_assertion(pt.number_of_vertices() == 0);

  // Center of the circle around the origin
  pts_rnd1000.clear();
  for (int i = 0 ; i < 1000 ; i++)
    {
      Point p = *on_circle++;
      pts_rnd1000.push_back(Point(p.x() < 0 ? 1 + p.x() : p.x(), p.y() < 0 ? 1 + p.y()  : p.y()));
    }

  pt.insert(pts_rnd1000.begin(), pts_rnd1000.end());
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.number_of_vertices() == 1000);
  pt.clear();
  CGAL_assertion(pt.is_valid());
  CGAL_assertion(pt.empty());
  CGAL_assertion(pt.number_of_vertices() == 0);
}

int main()
{
  test_insertion();

  return 0;
}
