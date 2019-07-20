#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/number_type_config.h> // CGAL_PI

#include <cmath>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K>                      Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                       C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                           Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Implicit function
FT double_p(const Point& p)
{
  const FT cx = std::cos(2 * CGAL_PI * p.x()),
           cy = std::cos(2 * CGAL_PI * p.y()),
           cz = std::cos(2 * CGAL_PI * p.z());
  const FT c2x = std::cos(4 * CGAL_PI * p.x()),
           c2y = std::cos(4 * CGAL_PI * p.y()),
           c2z = std::cos(4 * CGAL_PI * p.z());

  return 0.5 * (cx*cy + cy*cz + cz*cx) + 0.2 * (c2x + c2y + c2z);
}

int main(int argc, char** argv)
{
  // 'int' because the 'double_p' function is periodic over the domain only if
  // the length of the side of the domain is an integer
  int domain_size = (argc > 1) ? atoi(argv[1]) : 1;
  int number_of_copies_in_output = (argc > 2) ? atoi(argv[2]) : 8; // can be 1, 2, 4, or 8

  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  // there is no need for periodicity... ?
  Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(double_p, canonical_cube);

  Periodic_mesh_criteria criteria(facet_angle = 30,
                                  facet_size = 0.05 * domain_size,
                                  facet_distance = 0.025 * domain_size,
                                  cell_radius_edge_ratio = 2.,
                                  cell_size = 0.05);

  // Mesh generation with optimizers
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria,
                                                 odt(convergence=0.03, freeze_bound=0.02, time_limit=30),
                                                 lloyd(max_iteration_number=10),
                                                 perturb(sliver_bound=10, time_limit=30),
                                                 exude(sliver_bound=10, time_limit=0));

  std::ofstream medit_file("output_implicit_shape_optimized.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file, c3t3);

  // Below, the mesh generation and the optimizations are done in several calls
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria,
                                                     no_odt(), no_lloyd(),
                                                     no_perturb(), no_exude());

  std::ofstream medit_file_bis("output_implicit_shape_non-optimized.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file_bis, c3t3_bis);

  // Now, call each optimizer with its global function
  CGAL::odt_optimize_periodic_3_mesh_3(c3t3_bis, domain, convergence=0.03, freeze_bound=0.02, time_limit=30);
  CGAL::lloyd_optimize_periodic_3_mesh_3(c3t3_bis, domain, max_iteration_number=10);
  CGAL::perturb_periodic_3_mesh_3(c3t3_bis, domain, sliver_bound=10, time_limit=30);
  CGAL::exude_periodic_3_mesh_3(c3t3_bis, sliver_bound=10, time_limit=0);

  std::ofstream medit_file_ter("output_implicit_shape_two_steps.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file_ter, c3t3_bis, number_of_copies_in_output);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

