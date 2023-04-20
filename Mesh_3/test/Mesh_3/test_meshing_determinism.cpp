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
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <cstring>
#ifdef CGAL_LINKED_WITH_TBB
#define TBB_PREVIEW_GLOBAL_CONTROL 1
#  include <tbb/global_control.h>
#endif

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

template <typename Concurrency_tag>
void test()
{
  // Collect options
  std::size_t nb_runs   = 2;
  unsigned int nb_lloyd = 2;
  unsigned int nb_odt   = 2;
  double perturb_bound  = 10.;
  double exude_bound    = 15.;

  // Domain
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

  // Triangulation
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, K, Concurrency_tag>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<
    Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

  // Mesh Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // Domain
  std::cout << "\tSeed is\t 0" << std::endl;
  std::ifstream input(CGAL::data_file_path("meshes/cube.off"));
  Polyhedron polyhedron;
  input >> polyhedron;
  Mesh_domain domain(polyhedron);
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
  std::vector<std::string> output_surfaces;

  const std::size_t nb_operations = 5;

  output_c3t3.reserve(nb_operations * nb_runs);
  for(std::size_t i = 0; i < nb_runs; ++i)
  {
    std::cout << "------- Iteration " << (i+1) << " -------" << std::endl;
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                        no_perturb(),
                                        no_exude());
    std::ostringstream oss;
    CGAL::IO::write_MEDIT(oss, c3t3);
    output_c3t3.push_back(oss.str()); //[5*i]
    oss.clear();
    Polyhedron out_poly;
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_poly);
    oss << out_poly;
    output_surfaces.push_back(oss.str());//[5*i]
    out_poly.clear();
    oss.clear();

    //LLOYD (1)
    CGAL::lloyd_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_lloyd);
    CGAL::IO::write_MEDIT(oss, c3t3);
    output_c3t3.push_back(oss.str());//[i*5+1]
    oss.clear();
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_poly);
    oss << out_poly;
    output_surfaces.push_back(oss.str());//[i*5+1]
    out_poly.clear();
    oss.clear();

    //ODT (2)
    CGAL::odt_optimize_mesh_3(c3t3, domain, max_iteration_number = nb_odt);
    CGAL::IO::write_MEDIT(oss, c3t3);
    output_c3t3.push_back(oss.str());//[i*5+2]
    oss.clear();
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_poly);
    oss << out_poly;
    output_surfaces.push_back(oss.str());//[i*5+2]
    out_poly.clear();
    oss.clear();

    //PERTURB (3)
    CGAL::perturb_mesh_3(c3t3, domain, sliver_bound=perturb_bound);
    CGAL::IO::write_MEDIT(oss, c3t3);
    output_c3t3.push_back(oss.str());//[i*5+3]
    oss.clear();
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_poly);
    oss << out_poly;
    output_surfaces.push_back(oss.str());//[i*5+3]
    out_poly.clear();
    oss.clear();

    //EXUDE (4)
    CGAL::exude_mesh_3(c3t3, sliver_bound=exude_bound);
    CGAL::IO::write_MEDIT(oss, c3t3);
    output_c3t3.push_back(oss.str());//[i*5+4]
    oss.clear();
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, out_poly);
    oss << out_poly;
    output_surfaces.push_back(oss.str());//[i*5+4]
    out_poly.clear();
    oss.clear();

    if(i == 0)
      continue;
    //else check
    for(std::size_t j = 0; j < nb_operations; ++j)
    {
      std::size_t id1 = nb_operations * (i - 1) + j;
      std::size_t id2 = nb_operations * i + j;
      if(0 != output_c3t3[id1].compare(output_c3t3[id2]))
      {
        std::cerr << "Meshing operation " << j << " is not deterministic.\n";
        assert(false);
      }
      if (0 != output_surfaces[id1].compare(output_surfaces[id2]))
      {
        std::cerr << "Output surface after operation " << j << " is not deterministic.\n";
        assert(false);
      }
    }
  }
}

int main(int, char*[])
{
  std::cout << "Sequential test" << std::endl;
  test<CGAL::Sequential_tag>();

#ifdef CGAL_LINKED_WITH_TBB
  std::cout << "\n\nParallel with 1 thread test" << std::endl;
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
  test<CGAL::Parallel_tag>();
#endif
}
