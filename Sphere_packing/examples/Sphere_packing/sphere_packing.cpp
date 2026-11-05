#include <CGAL/Sphere_packing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Sphere_3 = CGAL::Sphere_3<Kernel>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, const char **argv) {
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/dino.off");

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<Sphere_3> spheres;
  spheres.clear();

  const float target_coverage = 0.98f;

  pack_spheres(mesh, std::back_inserter(spheres),
    CGAL::parameters::number_of_spheres(1000).target_coverage(target_coverage));

  std::cout << "Number of spheres: " << spheres.size() << std::endl;

  std::ofstream out("dino.spheres");
  double volume = 0;
  for (const auto& s : spheres) {
    double radius = sqrt(s.squared_radius());
    volume += (4.0 / 3.0) * CGAL_PI * radius * radius * radius;
    out << s.center() << " " << CGAL::sqrt(s.squared_radius()) << std::endl;
  }

  out.close();

  std::cout << "Total volume of spheres: " << volume << std::endl;
  std::cout << "Volume of the mesh: " << PMP::volume(mesh) << std::endl;
  std::cout << "Coverage: " << volume / PMP::volume(mesh) << std::endl;
  std::cout << "Coverage target: " << target_coverage << std::endl;
  return 0;
}
