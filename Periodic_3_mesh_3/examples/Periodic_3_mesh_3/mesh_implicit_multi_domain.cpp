#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_function_wrapper.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/number_type_config.h> // CGAL_PI

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (*Function)(const Point&);
  // This wrapper is needed to make 'sphere_function' periodic.
typedef CGAL::Periodic_3_function_wrapper<Function, K>      Periodic_function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Periodic_function> Multi_domain_wrapper;
typedef CGAL::Labeled_mesh_domain_3<K>                      Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type  Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                        C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                           Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Implicit functions
FT sphere_function(const Point& p)
{
  return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5)) - 0.15;
}

FT schwarz_p(const Point& p)
{
  const FT x2 = std::cos( p.x() * 2 * CGAL_PI ),
           y2 = std::cos( p.y() * 2 * CGAL_PI ),
           z2 = std::cos( p.z() * 2 * CGAL_PI );

  return x2 + y2 + z2;
}

int main(int argc, char** argv)
{
  int domain_size = (argc > 1) ? atoi(argv[1]) : 1;
  int number_of_copies_in_output = (argc > 2) ? atoi(argv[2]) : 4; // can be 1, 2, 4, or 8

  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  std::vector<Periodic_function> funcs;
  funcs.push_back(Periodic_function(&schwarz_p, canonical_cube));
  funcs.push_back(Periodic_function(&sphere_function, canonical_cube));

  // The vector of vectors of sign is passed as a vector of strings (since a string
  // is a vector of chars)
  std::vector<std::string> vps;
  vps.push_back("--");
  vps.push_back("-+");

  Multi_domain_wrapper multi_domain_function(funcs, vps);
  Periodic_mesh_domain domain(multi_domain_function, canonical_cube);

  Periodic_mesh_criteria criteria(facet_angle = 30,
                                  facet_size = 0.04,
                                  facet_distance = 0.025,
                                  cell_radius_edge_ratio = 2.,
                                  cell_size = 0.04);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("output_multi_domain.mesh");
  CGAL::output_periodic_mesh_to_medit<C3t3>(medit_file, c3t3, number_of_copies_in_output,
                                            false /*do not associate different colors to each copy*/,
                                            false /*do not rebind*/, true /*show patches*/);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
