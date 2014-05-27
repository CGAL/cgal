#include <CGAL/Periodic_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_3/Implicit_to_labeled_function_wrapper.h>
#include <CGAL/Periodic_mesh_facet_criteria_3.h>
#include <CGAL/Periodic_mesh_cell_criteria_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Periodic_mesh_3/Implicit_to_labeled_subdomains_function_wrapper.h>

#include <CGAL/Periodic_implicit_mesh_domain_3.h>
#include <CGAL/make_periodic_mesh_3.h>
#include <CGAL/Mesh_3_periodic_triangulation_3.h>

#include <CGAL/Mesh_constant_domain_field_3.h>


// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;

// Function
//typedef FT (*Function)(const Point&);
//typedef CGAL::Mesh_3::Implicit_multi_domain_to_labeled_function_wrapper<Function, K> Wrapper;
typedef FT (Function)(const Point&);

// Domain
typedef CGAL::Mesh_3::Implicit_to_labeled_function_wrapper<Function, K> Wrapper;
typedef CGAL::Periodic_implicit_mesh_domain_3<Function, K, Wrapper>                 Periodic_mesh_domain;

// Triangulation
typedef CGAL::Mesh_3_periodic_triangulation_3_generator<Periodic_mesh_domain>::type Mesh_3_periodic_triangulation_3;
typedef Mesh_3_periodic_triangulation_3                                             Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                                 C3t3;

// Edge criteria
typedef CGAL::Mesh_edge_criteria_3<Tr> Edge_criteria;
// Facet criteria
typedef CGAL::Periodic_mesh_facet_criteria_3<Tr> Periodic_facet_criteria;
// Cell criteria
typedef CGAL::Periodic_mesh_cell_criteria_3<Tr> Periodic_cell_criteria;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr, Edge_criteria, Periodic_facet_criteria, Periodic_cell_criteria> Mesh_criteria;


// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

const FT PI = std::acos(-1.);

// Function
FT sphere_function_1 (const Point& p)
{ return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5))-0.2; }
FT sphere_function_2 (const Point& p)
{ return CGAL::squared_distance(p, Point(0.66, 0.5, 0.5))-0.2; }


typedef CGAL::Mesh_constant_domain_field_3<Periodic_mesh_domain::R, Periodic_mesh_domain::Index> Field;

int main(int argc, char** argv)
{
  int domain_size = 1;
  std::string index = "1";

  if (argc > 1)
  {
	  domain_size = atof(argv[1]);
	  index = argv[1];
  }

//  std::vector<Function> funcs;
//  funcs.push_back(&sphere_function_1);
//  funcs.push_back(&sphere_function_2);
//  Periodic_mesh_domain domain(Wrapper(funcs, "--"), CGAL::Bbox_3(0, 0, 0, domain_size, domain_size, domain_size));
  Periodic_mesh_domain domain(sphere_function_1, CGAL::Iso_cuboid_3<K>(0, 0, 0, domain_size, domain_size, domain_size));
  
  Mesh_criteria criteria(domain
//  , facet_angle=30, facet_size=0.05 * domain_size, facet_distance=0.025 * domain_size,
//                              cell_radius_edge_ratio=2, cell_size = 0.05
//                              );
    , facet_angle=30, facet_size=0.10 * domain_size, facet_distance=0.25 * domain_size,
                                cell_radius_edge_ratio=2, cell_size = 0.10
                                 );

  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_mesh_3<C3t3>(domain, criteria);
  
  // Output
  std::ofstream medit_file( (std::string("output_") + index + ".mesh").data() );
  
  write_complex_to_medit(medit_file, c3t3, 1);
  
  medit_file.close();
  
  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
