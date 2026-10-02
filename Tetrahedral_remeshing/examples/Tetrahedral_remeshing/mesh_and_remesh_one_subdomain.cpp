#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Mesh_3/Detect_features_in_image.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>

#include <unordered_set>
#include <iostream>


using K = CGAL::Exact_predicates_inexact_constructions_kernel;
/// [Domain definition]
using Image_domain = CGAL::Labeled_mesh_domain_3<K>;
using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Image_domain>;
/// [Domain definition]

#ifdef CGAL_CONCURRENT_MESH_3
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// Triangulation for Remeshing
using Triangulation_3 = CGAL::Triangulation_3<Tr::Geom_traits, Tr::Triangulation_data_structure>;
// Constraints for protecting 1D-features
using Vertex_handle = Triangulation_3::Vertex_handle;
using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;


template <typename Tr>
struct Cells_of_subdomain_pmap
{
  using Subdomain_index
    = typename Tr::Triangulation_data_structure::Cell::Subdomain_index;
  const Subdomain_index m_subdomain;

public:
  using key_type = typename Tr::Cell_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;

  friend value_type get(const Cells_of_subdomain_pmap& map, const key_type& c) {
    return (map.m_subdomain == c->subdomain_index());
  }
  friend void put(Cells_of_subdomain_pmap&, const key_type&, const value_type) {
    ; // nothing to do : subdomain indices are updated in remeshing
  }
};


// To avoid verbose function and named parameters call
namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc>1) ? argv[1] : CGAL::data_file_path("images/420.inr");

  /// [Loads image]
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  /// [Loads image]

  /// [Domain creation]
  Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image,
         params::features_detector = CGAL::Mesh_3::Detect_features_in_image());
  /// [Domain creation]

  const CGAL::Bbox_3 bbox = domain.bbox();
  double diag = CGAL::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                           CGAL::square(bbox.ymax() - bbox.ymin()) +
                           CGAL::square(bbox.zmax() - bbox.zmin()));
  double sizing_default = diag * 0.05;

  /// [Mesh criteria]
  /// Note that `edge_size` is needed with 1D-features
  Mesh_criteria criteria(params::edge_size = sizing_default,
    params::facet_size = sizing_default,
    params::facet_distance = sizing_default / 10
  );
  /// [Mesh criteria]

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::no_exude().no_perturb());
  /// [Meshing]

  // Output
  std::ofstream out0("out_meshed.mesh");
  CGAL::IO::write_MEDIT(out0, c3t3.triangulation());
  out0.close();


  /// [Collect edge constraints before remeshing]
  Constraints_set constraints;
  Constraints_pmap constraints_pmap(constraints);

  Triangulation_3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3),
                         params::edge_is_constrained_map(constraints_pmap));
  /// [Collect edge constraints before remeshing]

  /// [Remeshing]
  const double remeshing_target_size = sizing_default * 0.5;
  std::cout << "Remeshing with target edge length: " << remeshing_target_size << std::endl;
  CGAL::tetrahedral_isotropic_remeshing(tr, remeshing_target_size,
      params::cell_is_selected_map(Cells_of_subdomain_pmap<Triangulation_3>{4})
     .number_of_iterations(5)
     .edge_is_constrained_map(constraints_pmap)
     .remesh_boundaries(false));
  /// [Remeshing]

  // Output
  std::ofstream out("out_remeshed.mesh");
  CGAL::IO::write_MEDIT(out, tr);
  out.close();

  return EXIT_SUCCESS;
}
