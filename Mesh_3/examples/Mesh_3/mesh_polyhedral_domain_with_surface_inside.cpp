#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Timer.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;


// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,CGAL::Sequential_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char*argv[])
{
  std::cout.precision(17);
  std::cerr.precision(17);
  const char* fname = (argc>1)?argv[1]:"data/horizons.off";
  std::ifstream input(fname);
  const char* fname2 = (argc>2)?argv[2]:"data/horizons-domain.off";
  std::ifstream input2(fname2);
  Polyhedron sm, smbounding;
  input >> sm;
  input2 >> smbounding;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  CGAL::Timer t;
  t.start();
  // Create domain
  Mesh_domain domain(sm, smbounding);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.025,
                         facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
                         cell_radius_edge_ratio = 3, cell_size = 0.05);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      no_perturb(), no_exude());

  std::cerr << t.time() << " sec." << std::endl;
  // Output
  dump_c3t3(c3t3, "out");
}
