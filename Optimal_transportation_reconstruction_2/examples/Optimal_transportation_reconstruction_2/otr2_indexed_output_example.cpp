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

typedef CGAL::Optimal_transportation_reconstruction_2<K>    Otr_2;

void load_xy_file(const std::string& filename, std::vector<Point>& points)
{
  std::ifstream ifs(filename.c_str());
  Point point;
  while (ifs >> point)
    points.push_back(point);

  ifs.close();
}

void indexed_output(Otr_2& otr2)
{
  std::cout << "(-------------Off output---------- )" << std::endl;

  std::vector<Point> points;
  std::vector<std::size_t> isolated_vertices;
  std::vector<std::pair<std::size_t,std::size_t> > edges;

  otr2.indexed_output(
      std::back_inserter(points),
      std::back_inserter(isolated_vertices),
      std::back_inserter(edges));

  std::cout << "OFF " << points.size() << " 0 " << edges.size()  << std::endl;

  // points
  std::vector<Point>::iterator pit;
  for (pit = points.begin(); pit != points.end(); pit++)
    std::cout << *pit << std::endl;

  // isolated vertices
  std::vector<std::size_t>::iterator vit;
  for (vit = isolated_vertices.begin(); vit != isolated_vertices.end(); vit++)
    std::cout << "1 "  << *vit << std::endl;

  // edges
  std::vector<std::pair<std::size_t, std::size_t> >::iterator eit;
  for (eit = edges.begin(); eit != edges.end(); eit++)
    std::cout << "2 "  << eit->first << " " << eit->second << std::endl;
}

int main ()
{
  std::vector<Point> points;

  load_xy_file("data/stair-noise00.xy", points);

  Otr_2 otr2(points);
  otr2.run(100); // 100 steps
  indexed_output(otr2);

  return 0;
}
