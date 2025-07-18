#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  std::vector<Point> input_points;
  std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;
  if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }
  PMP::repair_polygon_soup(input_points, input_triangles);
  PMP::triangulate_polygons(input_points, input_triangles);

  CGAL::Real_timer t;
  t.start();
  PMP::autorefine_triangle_soup(input_points, input_triangles,
                                CGAL::parameters::concurrency_tag(CGAL::Parallel_if_available_tag())
                                                 /* .apply_iterative_snap_rounding(true) */);
  t.stop();
  std::cout << "#points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " in " << t.time() << " sec." << std::endl;
  //  CGAL::IO::write_polygon_soup("autorefined.off", input_points, input_triangles, CGAL::parameters::stream_precision(17));


  #if 1
  PMP::orient_polygon_soup(input_points, input_triangles);

  CGAL::Surface_mesh<Kernel::Point_3> mesh;
  PMP::polygon_soup_to_polygon_mesh(input_points, input_triangles, mesh);
  PMP::remove_almost_degenerate_faces(mesh);

  input_points.clear();
  input_triangles.clear();
  PMP::polygon_mesh_to_polygon_soup(mesh, input_points, input_triangles);
  #endif

  t.reset();
  t.start();
  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(input_points, input_triangles);\
  t.stop();

  std::cout << "call to CDT3 done in " << t.time() << " sec." << std::endl;

  return ccdt.triangulation().dimension() == 3 ? 0 : 1;
}
