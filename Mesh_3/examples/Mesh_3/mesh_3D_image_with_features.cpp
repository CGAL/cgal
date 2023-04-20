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

#include <CGAL/Mesh_3/Detect_features_on_image_bbox.h>

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

namespace params = CGAL::parameters;

// Read input features
#include "read_polylines.h"

int main(int argc, char* argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/420.inr");
  // Loads image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  /// Load 1D-features
  const std::string lines_fname = (argc > 2) ? argv[2] : CGAL::data_file_path("images/420.polylines.txt");
  std::vector<std::vector<K::Point_3> > features_inside;
  if (!read_polylines(lines_fname, features_inside)) // see file "read_polylines.h"
  {
    std::cerr << "Error: Cannot read file " << lines_fname << std::endl;
    return EXIT_FAILURE;
  }

  /// [Domain creation]
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image,
    params::features_detector = CGAL::Mesh_3::Detect_features_on_image_bbox(),
    params::input_features  = std::cref(features_inside));//use std::cref to avoid a copy
  /// [Domain creation]

  /// Note that `edge_size` is needed with 1D-features [Mesh criteria]
  Mesh_criteria criteria(params::edge_size = 6.,
    params::facet_angle = 30,
    params::facet_size = 6,
    params::facet_distance = 4,
    params::cell_radius_edge_ratio = 3,
    params::cell_size = 8);
  /// [Mesh criteria]

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return EXIT_SUCCESS;
}
