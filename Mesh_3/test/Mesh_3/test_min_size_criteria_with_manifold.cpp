#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron>;


// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("meshes/wedge.off");
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(params::edge_size(0.07).
                                 facet_angle(25).facet_size(0.07).facet_distance(0.007).
                                 facet_min_size(0.07).
                                 cell_radius_edge_ratio(3).cell_size(0.07));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::manifold().no_exude().no_perturb());

  // Output
  std::ofstream file("out-sm.vtu");
  CGAL::IO::output_to_vtu(file, c3t3, CGAL::IO::ASCII);

  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, polyhedron);
  if(!CGAL::Polygon_mesh_processing::does_self_intersect(polyhedron)) {
    std::cerr << "ERROR: No self-intersection detected!\n";
    return EXIT_FAILURE;
  }
  else
    std::cerr <<
R"(OK: Self-intersection detected as expected. The manifold criterion of Mesh_3
was hampered by the facet minimal size criterion.
)";
  return EXIT_SUCCESS;
}
