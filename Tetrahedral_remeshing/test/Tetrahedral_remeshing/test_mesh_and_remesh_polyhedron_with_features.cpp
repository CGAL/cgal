#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tetrahedral_remeshing.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface_mesh> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  std::ifstream input(fname);
  Surface_mesh mesh;
  input >> mesh;
  if (input.fail()) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  if (!CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(mesh);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.025,
    facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
    cell_radius_edge_ratio = 3, cell_size = 0.05);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  auto t3 = CGAL::convert_to_triangulation_3(std::move(c3t3));

  double target_edge_length = 0.1;//coarsen the mesh
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    number_of_iterations(3).remesh_boundaries(false));

  assert(t3.is_valid());

  target_edge_length = 0.04;//re-refine the mesh
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    number_of_iterations(1).remesh_boundaries(true));

  assert(t3.is_valid());

  return EXIT_SUCCESS;
}
