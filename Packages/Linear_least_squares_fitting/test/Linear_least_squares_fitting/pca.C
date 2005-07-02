// Test program for the linear_least_square_fitting() functions.
// Pierre Alliez

#include <vector>
#include <cassert>
#include <stdlib.h>

#include <CGAL/Cartesian.h>

#include <CGAL/copy_n.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/point_generators_2.h>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Point_2           Point_2;
typedef K::Triangle_2        Triangle_2;
typedef K::Line_2            Line_2;
typedef K::Point_3           Point_3;

#include "pca_utils.h"


void test_2_point_set(const unsigned int nb_points,
                      const FT epsilon)
{
  std::cout << "2D: fit a line to a point set" << std::endl;

  // create random points nearby a segment
  std::vector<Point_2> points;
  Point_2 p = random_point_2();
  Point_2 q = random_point_2();
  std::cout << "  generate " << nb_points << 
       " 2D random points on a segment...";
  points_on_segment_2(p,q,nb_points,std::back_inserter(points));
  perturb_points_2(points.begin(),points.end(),epsilon);
  std::cout << "done" << std::endl;

  // fit a line
  std::cout << "  fit a 2D line...";
  Line_2 line;
  FT quality = linear_least_squares_fitting_2(points.begin(),points.end(),line);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  // dump to ps
  std::cout << "  dump to ps...";
  dump_ps("test_2d_points.ps",points,line);
  std::cout << "done" << std::endl;

  std::vector<Triangle_2> triangles;
}

/*
void test_2_triangle_set(const unsigned int nb_triangles,
                         const FT epsilon)
{
  std::cout << "2D: fit a line to a triangle set" << std::endl;

  // create 2D triangles
  std::vector<Triangle_2> triangles;

  std::cout << "  generate 2D triangles...";
  Point_2 a(0.2,0.2);
  Point_2 b(0.8,0.2);
  Point_2 c(0.8,0.4);
  Point_2 d(0.2,0.4);
  triangles.push_back(Triangle_2(a,b,c));
  triangles.push_back(Triangle_2(a,c,d));
  std::cout << "done" << std::endl;

  // fit a line
  std::cout << "  fit a 2D line...";
  Line_2 line;
  FT quality = linear_least_squares_fitting_2(triangles.begin(),triangles.end(),line);
  std::cout << "done (quality: " << quality << ")" << std::endl;

  // dump to ps
  std::cout << "  dump to ps...";
  dump_ps("test_2d_triangles.ps",triangles,line);
  std::cout << "done" << std::endl;
}
*/

int main()
{
  test_2_point_set(500,0.05);
  //test_2_triangle_set(100,0.01);
  return 0;
}
