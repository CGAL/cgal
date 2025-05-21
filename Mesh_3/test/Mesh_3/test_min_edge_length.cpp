#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Sizing_field_with_aabb_tree.h>


template <typename K,
          typename Polyhedron,
          typename Concurrency_tag>
struct Tester
{
  // Domain
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;

  // Triangulation
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain, K, Concurrency_tag>::type Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  void operator()(const std::string fname, const std::string out_fname)
  {
    std::ifstream input(fname);
    using namespace CGAL::parameters;

    Polyhedron polyhedron;
    input >> polyhedron;
    if (input.fail()) {
      std::cerr << "Error: Cannot read file " << fname << std::endl;
      return;
    }

    if (!CGAL::is_triangle_mesh(polyhedron)) {
      std::cerr << "Input geometry is not triangulated." << std::endl;
      return;
    }

    // Create domain
    Mesh_domain domain(polyhedron);

    domain.detect_features(40);

    CGAL::Sizing_field_with_aabb_tree<K, Mesh_domain> sf(0.074, domain);

    // Mesh criteria
    Mesh_criteria criteria(edge_size = sf,
                           edge_min_size = 0.1,
                           facet_distance = 0.0074,
                           facet_angle = 25,
                           facet_size = 0.074,
                           cell_radius_edge_ratio = 3,
                           cell_size = 0.074);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

    // Output
    CGAL::dump_c3t3(c3t3, out_fname);
  }
};

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/dragknob.off");

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Surface_mesh<K::Point_3>                      Surface_mesh;

  std::cout << "Sequential test" << std::endl;
  Tester<K, Surface_mesh, CGAL::Sequential_tag>()(fname, "out-dragknob-sequential");
#ifdef CGAL_LINKED_WITH_TBB
  std::cout << "Parallel test" << std::endl;
  Tester<K, Surface_mesh, CGAL::Parallel_tag>()(fname, "out-dragknob-parallel");
#endif

  return EXIT_SUCCESS;
}
