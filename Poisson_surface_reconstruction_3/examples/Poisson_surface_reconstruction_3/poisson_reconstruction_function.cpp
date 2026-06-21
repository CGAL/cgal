#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/IO/read_points.h>

#include <vector>
#include <fstream>
#include <filesystem>
#include <string>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Pwn;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char** argv)
{
  std::vector<Pwn> points;
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");
  if(!CGAL::IO::read_points(filename, std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                             .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read input file!" << std::endl;
    return EXIT_FAILURE;
  }

  Polyhedron output_mesh;

  double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
    (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

  if (CGAL::poisson_surface_reconstruction_delaunay
      (points.begin(), points.end(),
       CGAL::First_of_pair_property_map<Pwn>(),
       CGAL::Second_of_pair_property_map<Pwn>(),
       output_mesh, 0.25 * average_spacing))
    {
       std::string fname = std::filesystem::path(filename).filename().stem().string()+"-20-30-0.375.off";

        std::ofstream out(fname.c_str());
        out.precision(17);
        out << output_mesh;
        std::cout << "Output written to " << fname << std::endl;
    }
  else
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
