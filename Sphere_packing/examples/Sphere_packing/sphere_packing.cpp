#include <CGAL/Sphere_packing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

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
    pack_spheres(mesh, std::back_inserter(spheres),
      CGAL::parameters::number_of_spheres(1000).target_coverage(0.98f));

  std::ofstream out("dino.spheres");
  for (const auto& s : spheres)
    out << s.center() << " " << CGAL::sqrt(s.squared_radius()) << std::endl;

  out.close();

  return 0;
}
