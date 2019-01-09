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

// Domain
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point;
typedef K::Iso_cuboid_3                                             Iso_cuboid;

typedef FT (Function)(const Point&);
typedef CGAL::Periodic_3_function_wrapper<Function, K>              Periodic_function;
typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<K> > Periodic_mesh_domain;

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

// Implicit function
static const FT cx = 0.51, cy = 0.51, cz = 0.5;
static const FT scale = 0.9;
FT cone_function(const Point& p)
{
  const FT x = p.x(), y = p.y(), z = p.z();

  if (((x-cx)*(x-cx) + (y-cy)*(y-cy)) > scale * (z-cz)*(z-cz))
    return 1; // outside
  else
    return -1; // inside
}

// To obtain a good looking mesh at the base of the cone, we protect the base circle
void cone_polylines(Polylines& polylines)
{
  const FT z = 0.;
  const FT radius_at_z = CGAL::sqrt(scale * cz * cz);

  Polyline_3 polyline;
  for(int i=0; i<360; ++i)
  {
    polyline.push_back(Point(cx + radius_at_z * std::sin(i*CGAL_PI/180),
                             cy + radius_at_z * std::cos(i*CGAL_PI/180),
                             z));
  }
  polyline.push_back(polyline.front()); // close the line

  polylines.push_back(polyline);
}

int main(int argc, char** argv)
{
  int number_of_copies_in_output = (argc > 1) ? atoi(argv[1]) : 8; // can be 1, 2, 4, or 8

  // Domain
  const int domain_size = 1;
  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  Periodic_mesh_domain domain =
    Periodic_mesh_domain::create_implicit_mesh_domain(
      Periodic_function(cone_function, canonical_cube), canonical_cube);

  // Mesh criteria
  Periodic_mesh_criteria criteria(edge_size = 0.02 * domain_size,
                                  facet_angle = 0.05 * domain_size,
                                  facet_size = 0.02 * domain_size,
                                  cell_radius_edge_ratio = 2,
                                  cell_size = 0.5);

  // Create the features that we want to preserve
  Polylines polylines;
  cone_polylines(polylines);

  // Insert the features in the domain
  domain.add_features(polylines.begin(), polylines.end());

  // Insert a corner to make sure the apex of the cone is present in the mesh
  domain.add_corner(Point(0.51, 0.51, 0.5));

  // Mesh generation WITHOUT feature preservation (and no optimizers)
  C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, no_features(),
                                                 no_exude(), no_perturb());
  std::ofstream medit_file("output_implicit_shape_without_protection.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file, c3t3, number_of_copies_in_output);

  // Mesh generation WITH feature preservation (and no optimizers)
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria, features(),
                                                     no_exude(), no_perturb());
  std::ofstream medit_file_bis("output_implicit_shape_with_protection.mesh");
  CGAL::output_periodic_mesh_to_medit(medit_file_bis, c3t3_bis, number_of_copies_in_output);

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
