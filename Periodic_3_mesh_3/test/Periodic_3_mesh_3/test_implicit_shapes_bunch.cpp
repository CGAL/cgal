#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_3/IO/File_medit.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/Periodic_3_function_wrapper.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/array.h>
#include <CGAL/number_type_config.h> // CGAL_PI

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Point_3                                          Point;
typedef K::Iso_cuboid_3                                     Iso_cuboid;

// Domain
typedef FT (Function)(const Point&);
  // the wrapper is needed because not all the functions are triply periodic
typedef CGAL::Periodic_3_function_wrapper<Function, K>      Periodic_function;
typedef CGAL::Labeled_mesh_domain_3<K>                      Periodic_mesh_domain;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Periodic_mesh_domain>::type   Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>                         C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                           Periodic_mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Function
FT sphere(const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5))-0.2;
}

FT schwarz_p(const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT x2 = std::cos( p.x() * 2 * CGAL_PI ),
           y2 = std::cos( p.y() * 2 * CGAL_PI ),
           z2 = std::cos( p.z() * 2 * CGAL_PI );
  return x2 + y2 + z2;
}

// not triply periodic
FT scherk(const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  return std::exp( p.z() ) * std::cos( p.y() * 2 * CGAL_PI) - std::cos( p.x() * 2 * CGAL_PI );
}

// Triply Implicit Periodic Functions for meshing
FT schwarz_p_transl (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT x2 = std::cos(p.x() * 2 * CGAL_PI + CGAL_PI / 2.0),
           y2 = std::cos(p.y() * 2 * CGAL_PI + CGAL_PI / 2.0),
           z2 = std::cos(p.z() * 2 * CGAL_PI + CGAL_PI / 2.0);
  return x2 + y2 + z2;
}

FT gyroid (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  return cx * sy + cy * sz + cz * sx;
}

FT diamond (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz;
}

FT double_p (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT c2x = std::cos(2 * p.x() * 2 * CGAL_PI),
           c2y = std::cos(2 * p.y() * 2 * CGAL_PI),
           c2z = std::cos(2 * p.z() * 2 * CGAL_PI);
  return 0.5 * (cx * cy  + cy * cz + cz * cx ) + 0.2 * (c2x + c2y + c2z);
}

FT G_prime (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT c2x = std::cos(2 * p.x() * 2 * CGAL_PI),
           c2y = std::cos(2 * p.y() * 2 * CGAL_PI),
           c2z = std::cos(2 * p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  const FT s2x = std::sin(2 * p.x() * 2 * CGAL_PI),
           s2y = std::sin(2 * p.y() * 2 * CGAL_PI),
           s2z = std::sin(2 * p.z() * 2 * CGAL_PI);
  return 5 * (s2x * sz * cy + s2y * sx * cz + s2z * sy * cx)
           + 1 * (c2x * c2y + c2y * c2z + c2z * c2x);
}

FT lidinoid (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT c2x = std::cos(2 * p.x() * 2 * CGAL_PI),
           c2y = std::cos(2 * p.y() * 2 * CGAL_PI),
           c2z = std::cos(2 * p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  const FT s2x = std::sin(2 * p.x() * 2 * CGAL_PI),
           s2y = std::sin(2 * p.y() * 2 * CGAL_PI),
           s2z = std::sin(2 * p.z() * 2 * CGAL_PI);
  return 1 * (s2x * sz * cy + s2y * sx * cz + s2z * sy * cx)
           - 1 * (c2x * c2y + c2y * c2z + c2z * c2x) + 0.3;
}

FT D_prime (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT c2x = std::cos(2 * p.x() * 2 * CGAL_PI),
           c2y = std::cos(2 * p.y() * 2 * CGAL_PI),
           c2z = std::cos(2 * p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  return 1 * ( sx * sy * sz) + 1 * ( cx * cy * cz)
           - 1 * ( c2x * c2y + c2y * c2z + c2z * c2x) - 0.4;
}

FT split_p (const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT cx = std::cos(p.x() * 2 * CGAL_PI),
           cy = std::cos(p.y() * 2 * CGAL_PI),
           cz = std::cos(p.z() * 2 * CGAL_PI);
  const FT c2x = std::cos(2 * p.x() * 2 * CGAL_PI),
           c2y = std::cos(2 * p.y() * 2 * CGAL_PI),
           c2z = std::cos(2 * p.z() * 2 * CGAL_PI);
  const FT sx = std::sin(p.x() * 2 * CGAL_PI),
           sy = std::sin(p.y() * 2 * CGAL_PI),
           sz = std::sin(p.z() * 2 * CGAL_PI);
  const FT s2x = std::sin(2 * p.x() * 2 * CGAL_PI),
           s2y = std::sin(2 * p.y() * 2 * CGAL_PI),
           s2z = std::sin(2 * p.z() * 2 * CGAL_PI);
  return 1.1 * (s2x * sz * cy + s2y * sx * cz + s2z * sy * cx)
             - 0.2 * (c2x * c2y + c2y * c2z + c2z * c2x)
             - 0.4 * (cx + cy + cz);
}

FT cylinder(const Point& p)
{
  assert(p.x() >= 0 && p.y() >= 0 && p.z() >= 0 && p.x() < 1 && p.y() < 1 && p.z() < 1);
  const FT x2 = (p.y() - 0.5) * (p.y() - 0.5),
           y2 = (p.z() - 0.5) * (p.z() - 0.5);
  return x2 + y2 - 0.1;
}

typedef CGAL::Mesh_constant_domain_field_3<Periodic_mesh_domain::R,
                                           Periodic_mesh_domain::Index> Sizing_field;

int main()
{
  int domain_size = 1;
  Iso_cuboid canonical_cube(0, 0, 0, domain_size, domain_size, domain_size);

  // Array of the functions
  const int functions_count = 11;

  std::array<Periodic_function, functions_count> implicit_functions =
  {{
    Periodic_function(cylinder, canonical_cube),
    Periodic_function(D_prime, canonical_cube),
    Periodic_function(diamond, canonical_cube),
    Periodic_function(double_p, canonical_cube),
    Periodic_function(G_prime, canonical_cube),
    Periodic_function(gyroid, canonical_cube),
    Periodic_function(lidinoid, canonical_cube),
    Periodic_function(scherk, canonical_cube),
    Periodic_function(schwarz_p, canonical_cube),
    Periodic_function(sphere, canonical_cube),
    Periodic_function(split_p, canonical_cube)
  }};

  for(int i=0; i<functions_count; ++i)
  {
    // Periodic mesh domain
    Periodic_mesh_domain domain =
      Periodic_mesh_domain::create_implicit_mesh_domain(
        implicit_functions[i], canonical_cube);

    // Mesh criteria
    Periodic_mesh_criteria criteria(facet_angle = 30,
                                    facet_size = 0.05,
                                    facet_distance = 0.025,
                                    cell_radius_edge = 2,
                                    cell_size = 0.05);

    // Mesh generation
    std::cout << "Meshing implicit shape n°: " << i << std::endl;
    C3t3 c3t3 = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria);

    // File name
    std::stringstream index;
    index << i;
    std::string file_name = "output_implicit_shape_" + index.str() + ".mesh";

    // Output
    std::ofstream medit_file(file_name.c_str());
    CGAL::IO::output_periodic_mesh_to_medit(medit_file, c3t3);
  }

  std::cout << "EXIT SUCCESS" << std::endl;
  return 0;
}
