#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/STL.h>

#include <iostream>
#include <fstream>
#include <vector>

int main ()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_3 Point;
  typedef std::vector<std::size_t> Polygon;
  //make a tetrahedorn soup.
  std::vector<Point> ps(4);
  ps[0] = Point(0,0,0); ps[1] = Point(1,0,0);
  ps[2] = Point(0,1,0); ps[3] = Point(0,0,1);
  std::vector<Polygon> faces(4);

  faces[0].push_back(0);
  faces[0].push_back(2);
  faces[0].push_back(1);

  faces[1].push_back(0);
  faces[1].push_back(3);
  faces[1].push_back(2);

  faces[2].push_back(1);
  faces[2].push_back(2);
  faces[2].push_back(3);

  faces[3].push_back(0);
  faces[3].push_back(1);
  faces[3].push_back(3);

  std::ofstream os("tetra.stl");
  CGAL::write_STL(os, ps, faces);
  if(!os)
  {
    std::cerr<<"error during STL writing."<<std::endl;
    return 1;
  }
  os.close();
  ps.clear();
  faces.clear();
  std::ifstream is("tetra.stl");
  if(!CGAL::read_STL(is, ps, faces))
  {
    std::cerr<<"error during STL reading."<<std::endl;
    return 1;
  }
  if(ps.size() != 4 || faces.size() != 4)
  {
    std::cerr<<"error during STL file interpretation."<<std::endl;
    return 1;
  }
  return 0;
}
