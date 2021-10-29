#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/IO/OFF.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_3  = typename Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;

using TAG                     = CGAL::Sequential_tag;
using Surface_mesh            = CGAL::Surface_mesh<Point_3>;
using Polyhedron              = CGAL::Polyhedron_3<Kernel>;
using Affine_transformation_3 = CGAL::Aff_transformation_3<Kernel>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {

  const double error_bound = 1e-4;
  const std::string filepath = (argc > 1 ? argv[1] : CGAL::data_file_path("meshes/blobby.off"));

  // We create a tetrahedron, remesh it, and compute the distance.
  // The expected distance is error_bound.
  std::cout << std::endl << "* remeshing tetrahedron example:" << std::endl;

  Surface_mesh mesh1, mesh2;
  CGAL::make_tetrahedron(
    Point_3(0, 0, 0), Point_3(2, 0, 0),
    Point_3(1, 1, 1), Point_3(1, 0, 2), mesh1);
  mesh2 = mesh1;

  using edge_descriptor = typename boost::graph_traits<Surface_mesh>::edge_descriptor;
  Surface_mesh::Property_map<edge_descriptor, bool> is_constrained_map =
    mesh2.add_property_map<edge_descriptor, bool>("e:is_constrained", true).first;

  const double target_edge_length = 0.05;
  PMP::isotropic_remeshing(
    mesh2.faces(), target_edge_length, mesh2,
    PMP::parameters::edge_is_constrained_map(is_constrained_map));

  std::cout << "* one-sided bounded-error Hausdorff distance: " <<
    PMP::bounded_error_Hausdorff_distance<TAG>(mesh1, mesh2, error_bound) << std::endl;

  // We load a mesh, save it in two different containers, and
  // translate the second mesh by 1 unit. The expected distance is 1.
  std::cout << std::endl << "* moving mesh example:" << std::endl;

  Surface_mesh surface_mesh;
  CGAL::IO::read_OFF(filepath, surface_mesh);

  Polyhedron polyhedron;
  CGAL::IO::read_OFF(filepath, polyhedron);

  PMP::transform(Affine_transformation_3(CGAL::Translation(),
    Vector_3(FT(0), FT(0), FT(1))), polyhedron);

  std::cout << "* symmetric bounded-error Hausdorff distance: " <<
    PMP::bounded_error_symmetric_Hausdorff_distance<TAG>(surface_mesh, polyhedron, error_bound)
  << std::endl << std::endl;
  return EXIT_SUCCESS;
}
