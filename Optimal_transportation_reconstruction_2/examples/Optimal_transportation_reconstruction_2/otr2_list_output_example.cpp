#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Reconstruction_simplification_2.h>

#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef CGAL::Reconstruction_simplification_2<K>            Rs_2;

void load_xy_file(const std::string& filename, std::vector<Point>& points)
{
  std::ifstream ifs(filename.c_str());
  Point point;
  while (ifs >> point)
    points.push_back(point);

  ifs.close();
}

void list_output(Rs_2& rs2)
{
  std::cout << "(-------------List output---------- )" << std::endl;

  std::vector<Point> isolated_points;
  std::vector<Segment> segments;

  rs2.list_output(
    std::back_inserter(isolated_points), std::back_inserter(segments));

  std::vector<Point>::iterator pit;
  for (pit = isolated_points.begin(); pit != isolated_points.end(); pit++) 
    std::cout  << *pit << std::endl;

  std::vector<Segment>::iterator sit;
  for (sit = segments.begin(); sit != segments.end(); sit++) 
    std::cout << *sit << std::endl;
}

int main ()
{
  std::vector<Point> points;

  load_xy_file("data/stair-noise00.xy", points);

  Rs_2 rs2(points);
  rs2.run(100); // 100 steps
  list_output(rs2);

  return 0;
}
