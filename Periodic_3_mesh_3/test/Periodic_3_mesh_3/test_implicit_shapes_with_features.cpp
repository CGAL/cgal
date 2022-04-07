#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_function_wrapper.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;

typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point;
typedef K::Iso_cuboid_3                                             Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);
typedef CGAL::Periodic_3_function_wrapper<Function, K>              Periodic_function;
typedef CGAL::Mesh_domain_with_polyline_features_3<
          CGAL::Labeled_mesh_domain_3<K> >                          Periodic_mesh_domain;

// Polyline
typedef std::vector<Point>                                          Polyline_3;
typedef std::list<Polyline_3>                                       Polylines;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type             Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
          Tr, Periodic_mesh_domain::Corner_index, Periodic_mesh_domain::Curve_index>  C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                                   Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

////////////////////////////////////////////////////////////////////////////////
static const FT r = 0.6;
FT sphere_function(const Point& p)
{
  return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5)) - r*r;
}

void fill_sphere_polylines(Polylines& polylines)
{
  Polyline_3 polyline;

  for(int i = 0; i < 360; ++i)
  {
    Point p(0., 0.5 + 0.33 * std::cos(i*CGAL_PI/180), 0.5 + 0.33 * std::sin(i*CGAL_PI/180));
    polyline.push_back(p);
  }
  polyline.push_back(polyline.front()); // close the line

  polylines.push_back(polyline);
}

void test_protected_sphere()
{
  std::cout << "test_protected_sphere..." << std::endl;

  int domain_size = 1;
  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(
      Periodic_function(sphere_function, canonical_cube), canonical_cube);

  // Mesh criteria
  Periodic_mesh_criteria criteria(edge_size = 0.02 * domain_size,
                                  facet_angle = 0.05 * domain_size,
                                  facet_size = 0.02 * domain_size,
                                  cell_radius_edge_ratio = 2, cell_size = 0.5);

  // Create edge that we want to preserve
  Polylines polylines;
  fill_sphere_polylines(polylines);

  // Insert edge in domain
  domain.add_features(polylines.begin(), polylines.end());

  // Mesh generation with feature preservation
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, manifold());

  // Output
  std::ofstream medit_file_bis("protected_sphere.mesh");
  CGAL::IO::output_periodic_mesh_to_medit(medit_file_bis, c3t3_bis);
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
static const FT mx = 0.2501;
static const FT Mx = 0.74;
static const FT my = 0.2501;
static const FT My = 0.74;

FT squary_cylinder_function(const Point& p)
{
  if(p.x() < mx || p.x() > Mx ) return 1;
  if(p.y() < my || p.y() > My ) return 1;
  return -1;
}

void fill_squary_cylinder_polylines(Polylines& polylines)
{
  Polyline_3 polyline;

  polyline.push_back(Point(mx, mx, 0.));
  polyline.push_back(Point(mx, mx, 1.));
  polylines.push_back(polyline);

  FT z = 2.;
  polyline.clear();
  polyline.push_back(Point(mx, my, z));
  polyline.push_back(Point(Mx, my, z));
  polyline.push_back(Point(Mx, My, z));
  polyline.push_back(Point(mx, My, z));
  polyline.push_back(Point(mx, my, z));
  polylines.push_back(polyline);
}

void test_protected_squary_cylinder()
{
  std::cout << "test_protected_squary_cylinder..." << std::endl;

  int domain_size = 1;
  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(
      Periodic_function(squary_cylinder_function, canonical_cube), canonical_cube);

  // Mesh criteria
  Periodic_mesh_criteria criteria(edge_size = 0.02 * domain_size,
                                  facet_angle = 0.05 * domain_size,
                                  facet_size = 0.02 * domain_size,
                                  cell_radius_edge_ratio = 2, cell_size = 0.5);

  // Create edge that we want to preserve
  Polylines polylines;
  fill_squary_cylinder_polylines(polylines);

  // Insert edge in domain
  domain.add_features(polylines.begin(), polylines.end());

  // Mesh generation with feature preservation
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, manifold());

  // Output
  std::ofstream medit_file_bis("squary_cylinder.mesh");
  CGAL::IO::output_periodic_mesh_to_medit(medit_file_bis, c3t3_bis);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
FT squary_cylinder_function_2(const Point& p)
{
  FT min = 0.3, max = 0.7;

  if(p.x() > min && p.x() < max) return 1;
  if(p.y() > min && p.y() < max) return 1;

  return -1;
}

void fill_squary_cylinder_polylines_2(Polylines& polylines)
{
  FT min = 0.3, max = 0.7;
  Polyline_3 polyline;

  // vertical
  polyline.push_back(Point(min, min, -1.));
  polyline.push_back(Point(min, min, 0.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(1 + max, 1 + max, 2.));
  polyline.push_back(Point(1 + max, 1 + max, 3.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(max, min, 4.));
  polyline.push_back(Point(max, min, 5.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(min, max, 6.));
  polyline.push_back(Point(min, max, 7.));
  polylines.push_back(polyline);

  // vertical, on border
  polyline.clear();
  polyline.push_back(Point(min, 0, 0.));
  polyline.push_back(Point(min, 0, 1.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(0, min, 0.));
  polyline.push_back(Point(0, min, 1.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(0, max, 0.));
  polyline.push_back(Point(0, max, 1.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(max, 0, 0.));
  polyline.push_back(Point(max, 0, 1.));
  polylines.push_back(polyline);

  // horizontal
  polyline.clear();
  polyline.push_back(Point(min, min, 0.));
  polyline.push_back(Point(min, -min, 0.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(min, -min, 1.));
  polyline.push_back(Point(-min, -min, 1.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(-min, -min, 2.));
  polyline.push_back(Point(-min, min, 2.));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(-min, min, 3.));
  polyline.push_back(Point(min, min, 3.));
  polylines.push_back(polyline);
}

void test_protected_squary_cylinder_2()
{
  std::cout << "test_protected_squary_cylinder_2..." << std::endl;

  int domain_size = 1;
  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(
      Periodic_function(squary_cylinder_function_2, canonical_cube), canonical_cube);

  // Mesh criteria
  Periodic_mesh_criteria criteria(edge_size = 0.02 * domain_size,
                                  facet_angle = 0.05 * domain_size,
                                  facet_size = 0.02 * domain_size,
                                  cell_radius_edge_ratio = 2,
                                  cell_size = 0.05 * domain_size);

  // Create edge that we want to preserve
  Polylines polylines;
  fill_squary_cylinder_polylines_2(polylines);

  // Insert edge in domain
  domain.add_features(polylines.begin(), polylines.end());

  // Mesh generation with feature preservation
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, manifold());

  // Output
  std::ofstream medit_file_bis("squary_cylinder_2.mesh");
  CGAL::IO::output_periodic_mesh_to_medit(medit_file_bis, c3t3_bis);
}
////////////////////////////////////////////////////////////////////////////////

int main()
{
  test_protected_sphere();
  test_protected_squary_cylinder();
  test_protected_squary_cylinder_2();

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
