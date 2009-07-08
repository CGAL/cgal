#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Implicit_mesh_domain_3<Function,K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;

// Function
FT sphere_function (const Point& p)
{
  const FT x2=p.x()*p.x();
  const FT y2=p.y()*p.y();
  const FT z2=p.z()*p.z();

  return x2+y2+z2-1;
}

int main()
{
  // Domain (Warning: Sphere_3 constructor uses square radius !)
  Mesh_domain domain(sphere_function, K::Sphere_3(CGAL::ORIGIN, 2.));

  // Criteria
  Facet_criteria facet_criteria(30, 0.1, 0.025); // angle, size, approximation
  Cell_criteria cell_criteria(3, 0.1); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);

  return 0;
}

