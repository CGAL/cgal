#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>
#include <CGAL/IO/File_medit.h>

#include <string>
#include <unordered_set>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Surface_mesh = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface_mesh>;

#ifdef CGAL_CONCURRENT_MESH_3
using Concurrency_tag = CGAL::Parallel_if_available_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<
  Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// Triangulation for Remeshing
using Triangulation_3 = CGAL::Triangulation_3<Tr::Geom_traits,
                                              Tr::Triangulation_data_structure>;
using Vertex_handle = Triangulation_3::Vertex_handle;
using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 5;

  std::ifstream input(fname);
  Surface_mesh mesh;

  std::string filename(fname);
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
  const double size = 0.072;
  Mesh_criteria criteria(edge_size = size,
                         facet_angle = 25,
                         facet_size = size,
                         facet_distance = 0.1 * size,
                         cell_radius_edge_ratio = 2,
                         cell_size = size);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  Triangulation_3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3),
    edge_is_constrained_map(constraints_pmap));

  // Remeshing
  CGAL::tetrahedral_isotropic_remeshing(tr, size,
    number_of_iterations(nb_iter)
    .edge_is_constrained_map(constraints_pmap));

  std::ofstream out("out_remeshed.mesh");
  CGAL::IO::write_MEDIT(out, tr);
  out.close();

  return EXIT_SUCCESS;
}
