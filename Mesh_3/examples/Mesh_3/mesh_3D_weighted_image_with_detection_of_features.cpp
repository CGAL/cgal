#define CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
#define CGAL_MESH_3_VERBOSE 1

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
#include <CGAL/Mesh_3/generate_label_weights.h>
#include <CGAL/Mesh_3/Detect_features_in_image.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Image_domain = CGAL::Labeled_mesh_domain_3<K>;
using Mesh_domain = CGAL::Mesh_domain_with_polyline_features_3<Image_domain>;
/// [Domain definition]

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
using namespace CGAL::parameters;

/// [Add 1D features]

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("images/liver.inr.gz");
  // Loads image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  /// [Domain creation]
  const float sigma = (std::max)(image.vx(), (std::max)(image.vy(), image.vz()));
  CGAL::Image_3 img_weights =
    CGAL::Mesh_3::generate_label_weights(image, sigma);

  Mesh_domain domain
    = Mesh_domain::create_labeled_image_mesh_domain(image,
                                                    weights = std::ref(img_weights),
                                                    relative_error_bound = 1e-6,
                                                    features_detector = CGAL::Mesh_3::Detect_features_in_image());
  /// [Domain creation]

  CGAL::Bbox_3 bbox = domain.bbox();
  double diag = CGAL::sqrt((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
                         + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin())
                         + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));
  double sizing_default = diag * 0.05;

  /// Note that `edge_size` is needed with 1D-features [Mesh criteria]
  Mesh_criteria criteria(edge_size = sizing_default,
      facet_angle = 30, facet_size = sizing_default, facet_distance = sizing_default / 10,
      //facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH,
      cell_radius_edge_ratio = 0, cell_size = 0
  );
  /// [Mesh criteria]

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_exude(),
                                      CGAL::parameters::no_perturb());

  // Output
  CGAL::dump_c3t3(c3t3, "out");

  return 0;
}
