#include <CGAL/Sphere_packing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <cmath>
#include <iomanip>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Sphere_3 = CGAL::Sphere_3<Kernel>;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits_3<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, const char** argv) {
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

  CGAL::Side_of_triangle_mesh<Mesh, Kernel> inside(mesh);
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  // The sphere packing runs on GPU and thus does not use exact arithmetic and an exact inside check for each sphere will not work.
  // Thus, a tolerance is used to check if the sphere is inside the mesh. The tolerance is set to 1e-4 times the diagonal length of the bounding box of the mesh.
  double tol = std::sqrt(CGAL::to_double(tree.bbox().squared_diagonal_length())) * 1e-4;

  bool pass = true;
  for (const auto& s : spheres)
    if (CGAL::ON_BOUNDED_SIDE != inside(s.center())) {
      std::cerr << std::setprecision(17) << "Sphere center is not inside the mesh. c: (" << s.center() << ") r: " << CGAL::approximate_sqrt(s.squared_radius()) << std::endl;
      pass = false;
    }
    else {
      Point cp = tree.closest_point(s.center());
      FT d = CGAL::approximate_sqrt(CGAL::squared_distance(cp, s.center()));
      FT r = CGAL::approximate_sqrt(s.squared_radius());
      if (d - r < -tol) {
        std::cerr << std::setprecision(17) << "Sphere is not fully inside the mesh. c: (" << s.center() << ") r: " << r << " d: " << d << std::endl;
        pass = false;
      }
    }

  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
