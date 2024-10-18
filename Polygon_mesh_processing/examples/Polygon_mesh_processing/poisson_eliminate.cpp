#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/poisson_eliminate.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;


int main(int argc, char* argv[])
{
  std::string filename = std::filesystem::path(argv[1]).stem().string();
  Mesh sm;
  CGAL::IO::read_polygon_mesh(argv[1], sm);

  std::vector<Point_3> points;

  CGAL::Polygon_mesh_processing::poisson_eliminate(sm, std::back_inserter(points));

  std::string poisson_points = filename+"-poisson.xyz";
  std::ofstream out(poisson_points);
  out.precision(17);
  for(const Point_3& p : points){
    out << p << std::endl;
  }
  return 0;
}
