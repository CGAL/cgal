#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                          Point;
typedef CGAL::Surface_mesh<Point>                           Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef _WIN32
int wmain(int argc, wchar_t* argv[])
#else
int main(int argc, char* argv[])
#endif
{
  std::filesystem::path data = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/测试模型.off");

  Mesh mesh;
  PMP::IO::read_polygon_mesh(data, mesh);

  std::ifstream in(data);
  std::string line;
  while(std::getline(in, line))
  {
    std::cout << line << std::endl;
  }
  return 0;
}
