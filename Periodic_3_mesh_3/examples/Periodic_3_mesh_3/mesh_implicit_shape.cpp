#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

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
namespace params = CGAL::parameters;

// Implicit function
FT schwarz_p(const Point& p)
{
  const FT x2 = std::cos( p.x() * 2 * CGAL_PI ),
           y2 = std::cos( p.y() * 2 * CGAL_PI ),
           z2 = std::cos( p.z() * 2 * CGAL_PI );
  return x2 + y2 + z2;
}

int main(int argc, char** argv)
{
  // 'atoi' because the 'schwarz_p' function is periodic over the domain
  // only if the length of the side of the domain is an integer.
  const int x_span = (argc > 1) ? atoi(argv[1]) : 1;
  const int y_span = (argc > 2) ? atoi(argv[2]) : x_span;
  const int z_span = (argc > 3) ? atoi(argv[3]) : x_span;
  const int min_span = (std::min)({x_span, y_span, z_span});

  const int number_of_copies_in_output = (argc > 4) ? atoi(argv[4]) : 4; // can be 1, 2, 4, or 8

  const Iso_cuboid canonical_cube(0, 0, 0, x_span, y_span, z_span);

  const Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(schwarz_p, canonical_cube);

  Periodic_mesh_criteria criteria(params::facet_angle(30)
                                         .facet_size(0.035 * min_span)
                                         .facet_distance(0.025 * min_span)
                                         .cell_radius_edge_ratio(2.)
                                         .cell_size(0.05 * min_span));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);

  std::cout << "Created mesh with " << c3t3.number_of_vertices_in_complex() << " vertices"
                        << " and " <<  c3t3.number_of_cells_in_complex() << " cells" << std::endl;

  std::ofstream medit_file("output_implicit_shape.mesh");
  CGAL::IO::output_periodic_mesh_to_medit(medit_file, c3t3, number_of_copies_in_output);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

