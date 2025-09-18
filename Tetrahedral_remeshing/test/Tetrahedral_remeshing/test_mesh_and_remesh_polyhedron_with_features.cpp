#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/property_map.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <unordered_set>

#include <CGAL/Timer.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface_mesh> Mesh_domain;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

template <typename Concurrency_tag>
void mesh(const Surface_mesh& mesh)
{
  if(std::is_convertible<Concurrency_tag, CGAL::Parallel_tag>::value)
    std::cout << "* Concurrent Remeshing *" << std::endl;
  else
    std::cout << "* Sequential Remeshing *" << std::endl;
  CGAL::Real_timer timer;

  // Triangulation
  typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  // Create domain
  Mesh_domain domain(mesh);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.025,
    facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
    cell_radius_edge_ratio = 3, cell_size = 0.05);

  // Mesh generation
  timer.start();
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb().no_exude());
  timer.stop();
  std::cout << "\tmeshing time     : " << timer.time() << " seconds" << std::endl;

  // Property maps of constraints
  using Vertex_handle = typename Tr::Vertex_handle;
  using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
  using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
  using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

  using Corners_set = std::unordered_set<Vertex_handle, boost::hash<Vertex_handle>>;
  using Corners_pmap = CGAL::Boolean_property_map<Corners_set>;

  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  Corners_set corners;
  Corners_pmap corners_pmap(corners);

  auto t3 = CGAL::convert_to_triangulation_3(std::move(c3t3),
                                             edge_is_constrained_map(constraints_pmap)
                                            .vertex_is_constrained_map(corners_pmap));

  double target_edge_length = 0.2;//coarsen the mesh
  timer.reset();
  timer.start();
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
     number_of_iterations(3)
    .remesh_boundaries(false)
    .edge_is_constrained_map(constraints_pmap)
    .vertex_is_constrained_map(corners_pmap));
  timer.stop();
  std::cout << "\tremeshing time 1 : " << timer.time() << " seconds" << std::endl;

  assert(t3.is_valid());

  target_edge_length = 0.02;//re-refine the mesh
  timer.reset();
  timer.start();

  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
     number_of_iterations(1)
    .remesh_boundaries(true)
    .edge_is_constrained_map(constraints_pmap)
    .vertex_is_constrained_map(corners_pmap));

  timer.stop();
  std::cout << "\tremeshing time 2 : " << timer.time() << " seconds " << std::endl;

  assert(t3.is_valid());
}

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  std::ifstream input(fname);
  Surface_mesh m;
  input >> m;
  if(input.fail()) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(m)) {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Running test with the CGAL::Sequential_tag tag." << std::endl;
  mesh<CGAL::Sequential_tag>(m);

#if defined CGAL_ACTIVATE_CONCURRENT_MESHING && defined CGAL_LINKED_WITH_TBB
  std::cout << std::endl;
  std::cout << "Running test with the CGAL::Parallel_tag tag." << std::endl;
  mesh<CGAL::Parallel_tag>(m);
#endif

  return EXIT_SUCCESS;
}
