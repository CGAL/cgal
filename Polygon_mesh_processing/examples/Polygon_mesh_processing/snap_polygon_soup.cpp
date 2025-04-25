#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPICK;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;
typedef CGAL::Simple_cartesian<double> Cartesian;

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename Kernel>
struct Example{
  void operator()(const std::string& filename, int grid_size){
    const std::string out_file = "rounded_soup.off";
    std::vector<typename Kernel::Point_3> input_points;
    std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;
    std::cout << "Snap rounding on " << filename << " with " << typeid(Kernel).name() << "\n";
    if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
    {
      std::cerr << "Cannot read " << filename << "\n";
      return;
    }
    std::cout << "#points = " << input_points.size() << " and #triangles = " << input_triangles.size() << std::endl;
    std::cout << "Is 2-manifold: " << PMP::is_polygon_soup_a_polygon_mesh(input_triangles) << std::endl;

    std::vector<std::pair<std::size_t, std::size_t>> pairs_of_intersecting_triangles;
    PMP::triangle_soup_self_intersections(input_points, input_triangles, std::back_inserter(pairs_of_intersecting_triangles));
    std::cout << "Nb of pairs of intersecting triangles: " << pairs_of_intersecting_triangles.size() << std::endl;

    PMP::repair_polygon_soup(input_points, input_triangles);
    PMP::triangulate_polygons(input_points, input_triangles);

    CGAL::Real_timer t;
    t.start();
    bool success=PMP::autorefine_triangle_soup(input_points, input_triangles, CGAL::parameters::apply_iterative_snap_rounding(true).erase_all_duplicates(false).concurrency_tag(CGAL::Parallel_if_available_tag()).snap_grid_size(grid_size).number_of_iterations(15));
    t.stop();
    std::cout << "\nOutput:" << std::endl;
    std::cout << "#points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " in " << t.time() << " sec." << std::endl;
    if(success)
      std::cout << "Does self-intersect: " << PMP::does_triangle_soup_self_intersect<CGAL::Parallel_if_available_tag>(input_points, input_triangles) << std::endl;
    else
      std::cout << "ROUNDING FAILED" << std::endl;

    std::vector<Cartesian::Point_3> output_points;
    for(auto &p: input_points)
      output_points.emplace_back(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
    CGAL::IO::write_polygon_soup(out_file, output_points, input_triangles, CGAL::parameters::stream_precision(17));
    std::cout << "Is 2-manifold: " << PMP::orient_polygon_soup(input_points, input_triangles) << "\n\n"  << std::endl;
  }
};

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  const int grid_size = argc <= 2 ? 23
                                  : std::stoi(std::string(argv[2]));

  const bool epeck = argc <= 3 ? false
                               : std::string(argv[3])=="EPECK";

  if(epeck)
    Example<EPECK>()(filename, grid_size);
  else
    Example<EPICK>()(filename, grid_size);
  return 0;
}
