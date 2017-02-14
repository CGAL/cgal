#define MESH_3_VERBOSE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_3_periodic_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_periodic_facet_criteria_3.h>
#include <CGAL/Mesh_periodic_cell_criteria_3.h>

#include <CGAL/Periodic_implicit_mesh_domain_3.h>
#include <CGAL/make_periodic_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Periodic_implicit_mesh_domain_3<Function,K> Periodic_mesh_domain;

// Triangulation
typedef CGAL::Mesh_periodic_3_triangulation_3<Periodic_mesh_domain>::type Mesh_3_periodic_triangulation_3;
typedef Mesh_3_periodic_triangulation_3 Tr;

typedef Tr::Iso_cuboid Iso_cuboid;

typedef CGAL::Bbox_3 Bbox_3;

// Mesh complex
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Facet criteria
typedef CGAL::Mesh_periodic_facet_criteria_3<Tr> Periodic_facet_criteria;
// Cell criteria
typedef CGAL::Mesh_periodic_cell_criteria_3<Tr> Periodic_cell_criteria;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr, Periodic_facet_criteria, Periodic_cell_criteria> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// to be deleted
Mesh_3_periodic_triangulation_3 tr;
//

Point canonicalize_point(const Point& p) {
  const Point& result = Point( ( p.x() < 1. ? p.x() : p.x()-1. ),
    ( p.y() < 1. ? p.y() : p.y() - 1. ),
    ( p.z() < 1. ? p.z() : p.z() - 1. ) );
  assert( !(result.x() < 0) && (result.x() < 1) );
  assert( !(result.y() < 0) && (result.y() < 1) );
  assert( !(result.z() < 0) && (result.z() < 1) );
  return result;
}

FT sphere_function (const Point& p)
{ return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5))-0.2; }

const FT& PI = 3.14159265358979;

FT schwarz_p(const Point_3& point) {
  Point p = tr.canonicalize_point(point);
  
  const FT x2=std::cos( p.x() * 2*PI ), 
    y2=std::cos( p.y() * 2*PI ),
    z2=std::cos( p.z() * 2*PI ); 
  return x2 + y2 + z2;
}

int main()
{
  // Periodic mesh domain (Warning: Sphere_3 constructor uses squared radius !)
  Periodic_mesh_domain domain(schwarz_p, Bbox_3(0, 0, 0, 1, 1, 1));

  // Mesh criteria
  Periodic_facet_criteria facet_criteria(domain, 30, 0.05, 0.05);
  Periodic_cell_criteria cell_criteria(domain, 2, 0.05);

  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_periodic_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  //c3t3.output_to_medit(medit_file);
  write_complex_to_medit_8(medit_file, c3t3);

  return 0;
}
