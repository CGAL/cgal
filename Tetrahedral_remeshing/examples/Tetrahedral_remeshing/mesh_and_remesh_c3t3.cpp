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
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Mesh_polyhedron_3<K>::type;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K>;

#ifdef CGAL_CONCURRENT_MESH_3
using Concurrency_tag = CGAL::Parallel_if_available_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation for Meshing
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<
  Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// Triangulation for Remeshing
using Triangulation_3 = CGAL::Triangulation_3<Tr::Geom_traits,
                             Tr::Triangulation_data_structure>;
using Vertex_handle = Triangulation_3::Vertex_handle;

using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

using Corners_set = std::unordered_set<Vertex_handle, boost::hash<Vertex_handle>>;
using Corners_pmap = CGAL::Boolean_property_map<Corners_set>;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
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
  Mesh_criteria criteria(edge_size = 0.025,
    facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
    cell_radius_edge_ratio = 3, cell_size = 0.05);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Property map of constraints
  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  Corners_set corners;
  Corners_pmap corners_pmap(corners);

  Triangulation_3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3),
                      edge_is_constrained_map(constraints_pmap).
                      vertex_is_constrained_map(corners_pmap));

  //note we use the move semantic, with std::move(c3t3),
  //  to avoid a copy of the triangulation by the function
  //  `CGAL::convert_to_triangulation_3()`
  //  After the call to this function, c3t3 is an empty and valid C3t3.
  //It is possible to use :  CGAL::convert_to_triangulation_3(c3t3),
  //  Then the triangulation is copied and duplicated, and c3t3 remains as is.

  const double target_edge_length = 0.1;//coarsen the mesh
  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
   number_of_iterations(5)
   .smooth_constrained_edges(true)
   .edge_is_constrained_map(constraints_pmap));

  std::ofstream out("out_remeshed.mesh");
  CGAL::IO::write_MEDIT(out, tr);
  out.close();

  return EXIT_SUCCESS;
}
