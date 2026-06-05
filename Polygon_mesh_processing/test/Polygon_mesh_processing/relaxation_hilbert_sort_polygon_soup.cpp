#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/hilbert_sort.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

#include <vector>
#include <array>
#include <string>
#include <filesystem>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Points = std::vector<Point_3>;
using Face = std::array<std::size_t,3>;
using Polygons = std::vector<Face>;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Timer = CGAL::Timer;


int main(int argc, char* argv[])
{
  bool perform_hilbert_sort = argc > 2;
  const std::string input = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/Nefer.off");
  const std::string output = std::filesystem::path(input).stem().string() + "_sorted.off";
  Points points;
  Polygons polygons;
  Surface_mesh sm;
  Timer timer;

  CGAL::IO::read_polygon_soup(input, points, polygons);

  if(perform_hilbert_sort){
    timer.start();
    CGAL::Polygon_mesh_processing::hilbert_sort_polygon_soup(points, polygons);
    timer.stop();
    std::cout << "Hilbert sort took " << timer.time() << " seconds." << std::endl;
    timer.reset();
  }

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);
  timer.start();
  CGAL::Polygon_mesh_processing::tangential_relaxation(sm, CGAL::parameters::number_of_iterations(40));
  timer.stop();
  std::cout << "Tangential relaxation took " << timer.time() << " seconds." << std::endl;

  CGAL::IO::write_polygon_mesh(output, sm);

  return EXIT_SUCCESS;
}

