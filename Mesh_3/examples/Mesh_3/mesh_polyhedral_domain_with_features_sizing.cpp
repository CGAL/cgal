#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Sizing_field_with_aabb_tree.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron>;
using Features_sizing_field = CGAL::Sizing_field_with_aabb_tree<K, Mesh_domain>;


#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("meshes/fandisk.off");
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Features_sizing_field edges_sizing_field(0.07, domain);
  Mesh_criteria criteria(params::edge_size(edges_sizing_field).
                                 facet_distance(0.0072).
                                 cell_radius_edge_ratio(3));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_exude().no_perturb());

  // Output
  CGAL::dump_c3t3(c3t3, "out_sizing_field_with_aabb_tree");

  return EXIT_SUCCESS;
}
