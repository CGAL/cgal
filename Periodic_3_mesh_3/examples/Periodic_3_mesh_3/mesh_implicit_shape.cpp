#define CGAL_MESH_3_VERBOSE
#define _DEBUG
#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
#define CGAL_MESH_3_DEBUG_FACET_CRITERIA
#define CGAL_MESH_3_DEBUG_CELL_CRITERIA

#include <CGAL/Mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Medit_IO.h>
#include <CGAL/Implicit_periodic_3_mesh_domain_3.h>
#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/optimize_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <CGAL/number_type_config.h> // CGAL_PI

#include <cmath>
#include <iostream>
#include <fstream>

namespace P3M3_IO = CGAL::Periodic_3_mesh_3::IO;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_periodic_3_mesh_domain_3<Function,K> Periodic_implicit_mesh_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Periodic_implicit_mesh_domain> Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                       C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Implicit function
FT schwarz_p(const Point& p) {
  const FT x2 = std::cos( p.x() * 2 * CGAL_PI ),
           y2 = std::cos( p.y() * 2 * CGAL_PI ),
           z2 = std::cos( p.z() * 2 * CGAL_PI );
  return x2 + y2 + z2;
}

int main(int argc, char** argv)
{
  // 'int' because the 'schwarz_p' function is periodic over the domain only if
  // the length of the side of the domain is an integer
  int domain_size = 1;

  if (argc > 1)
    domain_size = atoi(argv[1]);

  Periodic_mesh_domain domain(schwarz_p, CGAL::Iso_cuboid_3<K>(0, 0, 0, domain_size, domain_size, domain_size));

  Periodic_mesh_criteria criteria(facet_angle = 30,
                                  facet_size = 0.05 * domain_size,
                                  facet_distance = 0.025 * domain_size,
                                  cell_radius_edge_ratio = 2.,
                                  cell_size = 0.05);

  // Mesh generation
#if 1
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria,
//                                                 odt(),
                                                 no_odt(),
                                                 lloyd(),
//                                                 no_lloyd(),
//                                                 exude(),
                                                 no_exude(),
//                                                 perturb(sliver_bound=10, time_limit=10)
                                                 no_perturb()
                                                 );

  std::cout << "Single call done" << std::endl;
  std::ofstream medit_file("output_implicit_shape.mesh");
  P3M3_IO::write_complex_to_medit(medit_file, c3t3);
  CGAL_assertion(c3t3.is_valid());
  CGAL_assertion(c3t3.triangulation().is_valid());
#endif

//  for(C3t3::Cells_in_complex_iterator cicit = c3t3.cells_in_complex_begin();
//      cicit != c3t3.cells_in_complex_end(); ++cicit)
//  {
//    std::cout << "Cell: " << &(*cicit) << std::endl;
//    for(int i=0; i<4; ++i)
//      std::cout << c3t3.triangulation().point(cicit, i).point() << std::endl;
//  }

//  for(C3t3::Facets_in_complex_iterator ficit = c3t3.facets_in_complex_begin();
//      ficit != c3t3.facets_in_complex_end(); ++ficit)
//  {
//    for(int i=0; i<3; ++i)
//      std::cout << ficit->first->vertex((ficit->second + i + 1)%4)->point().point() << std::endl;
//  }

#if 0
  // Mesh generation and optimization in several call
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria,
                                                     no_perturb(), no_exude());

  std::ofstream medit_file_bis("output_implicit_shape_non-optimized.mesh");
  P3M3_IO::write_complex_to_medit(medit_file_bis, c3t3_bis);

  CGAL_assertion(c3t3_bis.is_valid());
  CGAL_assertion(c3t3_bis.triangulation().is_valid());

  std::cout << "Optimizing..." << std::endl;
  CGAL::lloyd_optimize_periodic_3_mesh_3(c3t3_bis, domain);
//  CGAL::perturb_periodic_3_mesh_3(c3t3_bis, domain, sliver_bound=10, time_limit=10);

  std::cout << "Two-step call done" << std::endl;
  std::ofstream medit_file_ter("output_implicit_shape_two_steps.mesh");
  P3M3_IO::write_complex_to_medit(medit_file_ter, c3t3_bis);
#endif

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}

