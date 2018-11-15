#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_to_labeled_subdomains_function_wrapper.h>
#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_constant_domain_field_3.h>

#include <CGAL/number_type_config.h> // CGAL_PI

#include <cmath>
#include <fstream>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_to_labeled_subdomains_function_wrapper<Function, K> Function_wrapper;
typedef CGAL::Labeled_mesh_domain_3<K>                                     Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                       C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT schwarz_p(const Point& p)
{
  const FT x2 = std::cos( p.x() * 2 * CGAL_PI ),
           y2 = std::cos( p.y() * 2 * CGAL_PI ),
           z2 = std::cos( p.z() * 2 * CGAL_PI );

  return x2 + y2 + z2;
}

// This example uses a cell sizing field that depends on the subdomain index
typedef CGAL::Mesh_constant_domain_field_3<Periodic_mesh_domain::R,
                                           Periodic_mesh_domain::Index> Field;

int main(int argc, char** argv)
{
  int number_of_copies_in_output = (argc > 1) ? atoi(argv[1]) : 4; // can be 1, 2, 4, or 8

  Iso_cuboid canonical_cube(0, 0, 0, 1, 1, 1);

  Function_wrapper wrapper(schwarz_p);
  Periodic_mesh_domain domain(wrapper, canonical_cube);

  // Write the two different sizing fields for cells
  int volume_dimension = 3;
  Field size(2);
  size.set_size(0.1, volume_dimension, domain.index_from_subdomain_index(2)); // exterior
  size.set_size(0.03, volume_dimension, domain.index_from_subdomain_index(1)); // interior

  Periodic_mesh_criteria criteria(facet_angle = 30.,
                                  facet_size = 0.05,
                                  facet_distance = 0.025,
                                  cell_radius_edge_ratio = 2.,
                                  cell_size = size);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("output_implicit_with_subdomains.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file, c3t3, number_of_copies_in_output);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

