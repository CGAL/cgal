#include <CGAL/Periodic_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Periodic_mesh_facet_criteria_3.h>
#include <CGAL/Periodic_mesh_cell_criteria_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Periodic_mesh_3/Implicit_to_labeled_subdomains_function_wrapper.h>

#include <CGAL/Periodic_implicit_mesh_domain_3.h>
#include <CGAL/Periodic_labeled_mesh_domain_3.h>
#include <CGAL/make_periodic_mesh_3.h>
#include <CGAL/Mesh_3_periodic_triangulation_3.h>

#include <CGAL/Mesh_constant_domain_field_3.h>


// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;

// Function
typedef FT (*Function)(const Point&);
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function> Wrapper;

// Domain
typedef CGAL::Periodic_labeled_mesh_domain_3<Wrapper, K>                 Periodic_mesh_domain;

// Triangulation
typedef CGAL::Mesh_periodic_3_triangulation_3<Periodic_mesh_domain>::type Mesh_3_periodic_triangulation_3;
typedef Mesh_3_periodic_triangulation_3                                             Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                                 C3t3;

// Edge criteria
typedef CGAL::Mesh_edge_criteria_3<Tr> Edge_criteria;
// Facet criteria
typedef CGAL::Periodic_mesh_facet_criteria_3<Tr> Periodic_facet_criteria;
// Cell criteria
typedef CGAL::Periodic_mesh_cell_criteria_3<Tr> Periodic_cell_criteria;
// Criteria
typedef CGAL::Periodic_3_mesh_criteria_3<Tr, Edge_criteria, Periodic_facet_criteria, Periodic_cell_criteria> Mesh_criteria;


// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

const FT PI = 3.14159265358979;

// Function
FT sphere_function_1 (const Point& p)
{ return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5))-0.15; }
FT sphere_function_2 (const Point& p)
{ return CGAL::squared_distance(p, Point(0.66, 0.5, 0.5))-0.2; }
FT schwarz_p(const Point& p) {
  const FT x2=std::cos( p.x() * 2*PI ),
  y2=std::cos( p.y() * 2*PI ),
  z2=std::cos( p.z() * 2*PI );
  return x2 + y2 + z2;
}

int main(int argc, char** argv)
{
  int domain_size = 1;
  std::string index = ".";

  if (argc > 1)
  {
    index = argv[1];
    if (argc > 2)
      domain_size = atof(argv[2]);
  }

  std::vector<Function> funcs;
  funcs.push_back(&schwarz_p);
  funcs.push_back(&sphere_function_1);
  std::vector<std::string> vps;
  vps.push_back("--");
  vps.push_back("-+");
  Wrapper wrapper(funcs, vps);
  Periodic_mesh_domain domain(wrapper, CGAL::Iso_cuboid_3<K>(0, 0, 0, domain_size, domain_size, domain_size));
  
  Mesh_criteria criteria(domain
      , facet_angle=30, facet_size=0.05, facet_distance=0.025,
      cell_radius_edge_ratio=2, cell_size = 0.05
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
