#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/property_map.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_medit.h>
#include <string>
#include <unordered_set>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation for Meshing
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// Triangulation for Remeshing
typedef CGAL::Triangulation_3<typename Tr::Geom_traits,
  typename Tr::Triangulation_data_structure> Triangulation_3;
using Vertex_handle = Triangulation_3::Vertex_handle;

using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

struct Distance_from_corner_sizing_field
{
  typedef K::FT FT;
  typedef K::Point_3 Point_3;

  const Point_3 corner{ -1., -1., -1 }; //lower corner of the cube

  template<typename Index>
  FT operator()(const Point_3& p, const int dim, const Index&) const
  {
    if (dim < 3)
      return 0.;
    const FT d_to_origin = CGAL::approximate_sqrt(CGAL::squared_distance(p, corner));
    return 0.02 + 0.2 * d_to_origin;
  }
};

// To avoid verbose function and named parameters call
namespace p = CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if (input.fail() || !CGAL::is_triangle_mesh(polyhedron)) {
    std::cerr << "Error: Input invalid " << fname << std::endl;
    return EXIT_FAILURE;
  }

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
    Distance_from_corner_sizing_field(),
    p::number_of_iterations(3)
    .edge_is_constrained_map(constraints_pmap)
    .smooth_constrained_edges(true));

  std::ofstream out("out_remeshing.mesh");
  CGAL::IO::write_MEDIT(out, tr);
  out.close();

  return EXIT_SUCCESS;
}
