#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/lloyd_optimize_mesh_3.h>
#include <CGAL/odt_optimize_mesh_3.h>
#include <CGAL/perturb_mesh_3.h>
#include <CGAL/exude_mesh_3.h>

#include <sstream>
#include <cstring>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int, char*[])
{
  // Collect options
  std::size_t nb_runs   = 2;
  unsigned int nb_lloyd = 2;
  unsigned int nb_odt   = 2;
  double perturb_bound  = 10.;
  double exude_bound    = 15.;

  // Domain
  std::cout << "\tSeed is\t 0" << std::endl;
  Mesh_domain domain("data/cube.off");
    //no random generator is given, so CGAL::Random(0) is used

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.2,
                         facet_angle = 25,
                         facet_size = 0.2,
                         facet_distance = 0.002,
                         cell_radius_edge_ratio = 3,
                         cell_size = 0.2);

  // iterate
  std::vector<std::string> output_c3t3;
  output_c3t3.reserve(5 * nb_runs);
  for(std::size_t i = 0; i < nb_runs; ++i)
  {
    std::cout << "------- Iteration " << (i+1) << " -------" << std::endl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        no_perturb(),
                                        no_exude());
    std::ostringstream oss;
    c3t3.output_to_medit(oss);
    output_c3t3.push_back(oss.str()); //[5*i]
    oss.clear();

    //LLOYD (1)
    CGAL::lloyd_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_lloyd);
    c3t3.output_to_medit(oss);
    output_c3t3.push_back(oss.str());//[i*5+1]
    oss.clear();

    //ODT (2)
    CGAL::odt_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_odt);
    c3t3.output_to_medit(oss);
    output_c3t3.push_back(oss.str());//[i*5+2]
    oss.clear();

    //PERTURB (3)
    CGAL::perturb_mesh_3(c3t3, domain, sliver_bound=perturb_bound);
    c3t3.output_to_medit(oss);
    output_c3t3.push_back(oss.str());//[i*5+3]
    oss.clear();

    //EXUDE (4)
    CGAL::exude_mesh_3(c3t3, sliver_bound=exude_bound);
    c3t3.output_to_medit(oss);
    output_c3t3.push_back(oss.str());//[i*5+4]
    oss.clear();

    if(i == 0)
      continue;
    //else check
    for(std::size_t j = 0; j < 5; ++j)
    {
      if(0 != output_c3t3[5*(i-1)+j].compare(output_c3t3[5*i+j]))
      {
        std::cerr << "Meshing operation " << j << " is not deterministic.\n";
        CGAL_assertion(false);
      }
    }
  }

  return 0;
}
