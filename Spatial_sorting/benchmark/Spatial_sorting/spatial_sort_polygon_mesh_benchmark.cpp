// #define CGAL_PROFILE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/hilbert_sort.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
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
  const std::string input = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
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


  std::size_t total = 0;
  for(int i = 0; i < 100; ++i)
  {
    for(auto v : vertices(sm)){
      std::pair<int,int>  minmax = std::make_pair(std::numeric_limits<int>::max(), 0);
      for(auto h : halfedges_around_target(halfedge(v,sm), sm))
      {
        ++total;
        if(h.idx()<minmax.first) minmax.first = h.idx();
        if(h.idx()>minmax.second) minmax.second = h.idx();
      }
        CGAL_HISTOGRAM_PROFILER("spread", minmax.second - minmax.first);
    }
  }

  timer.stop();
  std::cout << "halfedges_around_target() " << timer.time() << " seconds." << std::endl;

  std::cout << total << std::endl;

  return EXIT_SUCCESS;
}


