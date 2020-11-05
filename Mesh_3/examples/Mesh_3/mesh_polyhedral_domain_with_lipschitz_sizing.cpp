#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>

//Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

// Domain
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

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

// Sizing field
typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain, Mesh_domain::AABB_tree> Lip_sizing;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/fandisk.off";
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if (input.fail()){
    std::cerr << "Error: Cannot read file " << fname << std::endl;
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

  // Create Lipschitz sizing field
  Lip_sizing lip_sizing(domain, &domain.aabb_tree());

  FT min_size = 0.02;
  lip_sizing.add_parameters_for_subdomain(1,       //subdomain id
                                          0.3,     //k
                                          min_size,//min_size
                                          0.5);    //max_size

  // Mesh criteria
  Mesh_criteria criteria(edge_size = min_size,
                         facet_angle = 25,
                         facet_size = min_size,
                         facet_distance = 0.005,
                         cell_radius_edge_ratio = 3,
                         cell_size = lip_sizing);
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return EXIT_SUCCESS;
}
