#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;

typedef CGAL::Sequential_tag Concurrency_tag;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,
                                   CGAL::Default,
                                   Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Mesh_3::Mesher_3<C3t3,
                               Mesh_criteria,
                               Mesh_domain>   Mesher;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  const std::string fname = CGAL::data_file_path("meshes/halfcube.off");
  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain

  std::vector<Polyhedron*> polyhedra;
  polyhedra.push_back(&polyhedron);

  Mesh_domain domain(polyhedra.begin(), polyhedra.end());
  domain.detect_features();

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle = 25,
                         facet_size = 0.17,
                         facet_distance = 0.017,
                         edge_size = 0.17);

  // Mesh generation
  namespace p = CGAL::parameters;
  C3t3 c3t3;
  CGAL::Mesh_3::internal::C3t3_initializer<C3t3,
    Mesh_domain,
    Mesh_criteria,
    true>()(c3t3, domain, criteria, true);

  Mesher mesher(c3t3, domain, criteria, CGAL::FACET_VERTICES_ON_SURFACE);
  mesher.initialize();
  mesher.display_number_of_bad_elements();
  while ( ! mesher.is_algorithm_done() ) mesher.one_step();
  assert(c3t3.triangulation().number_of_vertices() > 200);

  // Output
  mesher.display_number_of_bad_elements();

  return EXIT_SUCCESS;
}
