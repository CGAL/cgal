#include <vector>
#include <iostream>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

/// [Domain definition]
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Image_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Image_domain> Mesh_domain;
/// [Domain definition]

#include <CGAL/Mesh_3/Detect_features_in_image.h>

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/420.inr");

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

  CGAL::Bbox_3 bbox = domain.bbox();
  double diag = CGAL::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                           CGAL::square(bbox.ymax() - bbox.ymin()) +
                           CGAL::square(bbox.zmax() - bbox.zmin()));
  double sizing_default = diag * 0.05;

  /// [Mesh criteria]
  /// Note that `edge_size` is needed with 1D-features
  Mesh_criteria criteria(params::edge_size = sizing_default,
    params::facet_angle = 30,
    params::facet_size = sizing_default,
    params::facet_distance = sizing_default / 10,
    params::facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
    params::cell_radius_edge_ratio = 0,
    params::cell_size = 0
  );
  /// [Mesh criteria]

  /// [Meshing]
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      params::no_exude(),
                                      params::no_perturb());
  /// [Meshing]

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
