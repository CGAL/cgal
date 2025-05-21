
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <cstdlib>
#include <iostream>

#include <boost/container/flat_set.hpp>

typedef float Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

template <typename Set>
struct Image_to_multiple_iso_level_sets {
  const Set& set;
  Image_to_multiple_iso_level_sets(const Set& set) : set(set) {}

  int operator()(double v) const {
       return int(std::distance(set.begin(),
                                set.lower_bound(float(v))));
  }
};

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1)?argv[1]:CGAL::data_file_path("images/skull_2.9.inr");
  // Load image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  typedef boost::container::flat_set<float> Flat_set;
  Flat_set iso_values;
  if(argc < 2) {
    iso_values.insert(1.5f);
    iso_values.insert(2.9f);
    iso_values.insert(3.5f);
  } else {
    for(int i = 2; i < argc; ++i) {
      iso_values.insert(static_cast<float>(std::atof(argv[i])));
    }
  }

  // Domain
  Mesh_domain domain =
    Mesh_domain::create_gray_image_mesh_domain
    (params::image(image).
             image_values_to_subdomain_indices(
               Image_to_multiple_iso_level_sets<Flat_set>(iso_values)).
             value_outside(0.f)
     );

  // Mesh criteria
  Mesh_criteria criteria(params::facet_angle(30).facet_size(6).facet_distance(2).
                         cell_radius_edge_ratio(3).cell_size(8));

  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::IO::write_MEDIT(medit_file, c3t3);
  medit_file.close();

  return 0;
}
