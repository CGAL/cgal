#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_2                                          Point;
typedef K::Segment_2                                        Segment;

typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr_2;

void load_xy_file(const std::string& filename, std::vector<Point>& points)
{
  std::ifstream ifs(filename.c_str());
  Point point;
  while (ifs >> point)
    points.push_back(point);

  ifs.close();
}

void list_output(Otr_2& otr2)
{
  std::cout << "(-------------List output---------- )" << std::endl;

  std::vector<Point> isolated_points;
  std::vector<Segment> segments;

  otr2.list_output(
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

  Otr_2 otr2(points);
  otr2.run(100); // 100 steps
  list_output(otr2);

  return 0;
}
