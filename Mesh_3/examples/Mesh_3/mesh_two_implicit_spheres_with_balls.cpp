#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Domain
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Mesh_domain_with_polyline_features_3<
  CGAL::Implicit_mesh_domain_3<Function,K> >              Mesh_domain;

// Polyline
typedef std::vector<Point>        Polyline_3;
typedef std::list<Polyline_3>       Polylines;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere_function1 (const Point& p)
{ return CGAL::squared_distance(p, Point(CGAL::ORIGIN))-2; }

FT sphere_function2 (const Point& p)
{ return CGAL::squared_distance(p, Point(2, 0, 0))-2; }

FT sphere_function (const Point& p)
{
  if(sphere_function1(p) < 0 || sphere_function2(p) < 0)
    return -1;
  else 
    return 1;
}

#include <cmath>

int main()
{
  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(sphere_function, K::Sphere_3(Point(1, 0, 0), 6.));

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.15,
                         facet_angle = 25, facet_size = 0.15,
                         cell_radius_edge_ratio = 2, cell_size = 0.15);
  
  // Create edge that we want to preserve
  Polylines polylines (1);
  Polyline_3& polyline = polylines.front();
  
  for(int i = 0; i < 360; ++i)
  {
    Point p (1, std::cos(i*CGAL_PI/180), std::sin(i*CGAL_PI/180));
    polyline.push_back(p);
  }
  polyline.push_back(polyline.front()); // close the line

  // Insert edge in domain
  domain.add_features(polylines.begin(), polylines.end());
  
  // Mesh generation without feature preservation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_features());

  std::ofstream medit_file("out-no-protection.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();
  c3t3.clear();

  // Mesh generation with feature preservation
  c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
  
  // Output
  medit_file.open("out-with-protection.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  return 0;
}
