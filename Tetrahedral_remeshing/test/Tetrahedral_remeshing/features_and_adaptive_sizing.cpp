#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/property_map.h>

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Adaptive_remeshing_sizing_field.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/Real_timer.h>

#include <string>
#include <unordered_set>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// To avoid verbose function and named parameters call
namespace p = CGAL::parameters;

template <typename Concurrency_tag>
int mesh(const std::string fname)
{
  // Triangulation for Meshing
  using Tr = typename CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type;
  using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>;

  // Criteria
  using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

  // Triangulation for Remeshing
  using Triangulation_3 = CGAL::Triangulation_3<typename Tr::Geom_traits, typename Tr::Triangulation_data_structure>;
  using Vertex_handle = typename Triangulation_3::Vertex_handle;
  using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
  using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
  using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if (input.fail() || !CGAL::is_triangle_mesh(polyhedron)) {
    std::cerr << "Error: Input invalid " << fname << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  // Create domain
  Mesh_domain domain(polyhedron);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(p::edge_size = 0.1,
                         p::facet_size = 0.1,
                         p::facet_distance = 0.005);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, p::no_perturb().no_exude());

  timer.stop();
  std::cout << "Meshing done (" << timer.time()
            << " seconds)" << std::endl;

  //CGAL::dump_c3t3(c3t3, "out_meshing");

  //Remeshing : extract triangulation
  timer.reset();
  timer.start();

  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  Triangulation_3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3),
      p::edge_is_constrained_map(constraints_pmap));

  //note we use the move semantic, with std::move(c3t3),
  //  to avoid a copy of the triangulation by the function
  //  `CGAL::convert_to_triangulation_3()`
  //  After the call to this function, c3t3 is an empty and valid C3t3.
  //It is possible to use :  CGAL::convert_to_triangulation_3(c3t3),
  //  Then the triangulation is copied and duplicated, and c3t3 remains as is.

  CGAL::tetrahedral_isotropic_remeshing(tr,
    CGAL::create_adaptive_remeshing_sizing_field(tr, p::edge_is_constrained_map(constraints_pmap)),
    p::number_of_iterations(3)
    .edge_is_constrained_map(constraints_pmap)
    .smooth_constrained_edges(true));

  timer.stop();
  std::cout << "Remeshing done (" << timer.time()
            << " seconds)" << std::endl;

  //std::ofstream out("out_remeshing.mesh");
  //CGAL::IO::write_MEDIT(out, tr);
  //out.close();
  return EXIT_SUCCESS;

}

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");

  bool success = true;
#ifdef CGAL_LINKED_WITH_TBB
  std::cout << "Running test with the CGAL::Parallel_if_available_tag tag." << std::endl;
  success = (EXIT_SUCCESS == mesh<CGAL::Parallel_if_available_tag>(fname));
  std::cout << std::endl;
#endif

  std::cout << "Running test with the CGAL::Sequential_tag tag." << std::endl;
  success = success && (EXIT_SUCCESS == mesh<CGAL::Sequential_tag>(fname));

  return success;
}
