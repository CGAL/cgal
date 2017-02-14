#include <CGAL/Periodic_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Periodic_mesh_facet_criteria_3.h>
#include <CGAL/Periodic_mesh_cell_criteria_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Periodic_implicit_mesh_domain_3.h>
#include <CGAL/make_periodic_mesh_3.h>
#include <CGAL/Mesh_3_periodic_triangulation_3.h>

#include <CGAL/Mesh_constant_domain_field_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Periodic_implicit_mesh_domain_3<Function,K> Periodic_mesh_domain;

// Triangulation
typedef CGAL::Mesh_periodic_3_triangulation_3<Periodic_mesh_domain>::type Mesh_3_periodic_triangulation_3;
typedef Mesh_3_periodic_triangulation_3 Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

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

// Function
FT sphere(const Point& p)
{ return CGAL::squared_distance(p, Point(0.5, 0.5, 0.5))-0.2; }

const FT PI = std::acos(-1.);

FT schwarz_p(const Point& p) {
  const FT x2=std::cos( p.x() * 2*PI ), 
  y2=std::cos( p.y() * 2*PI ),
  z2=std::cos( p.z() * 2*PI ); 
  return x2 + y2 + z2;
}

// not triply periodic
FT scherk(const Point& p) {
  return std::exp( p.z() ) * std::cos( p.y() * 2 * PI) - std::cos( p.x() * 2 * PI );
}

// Triply Implicit Periodic Functions for meshing
FT schwarz_p_transl (const Point& p) {
  const FT x2 = std::cos(p.x() * 2 * PI + PI / 2.0), 
  y2 = std::cos(p.y() * 2 * PI + PI / 2.0),
  z2 = std::cos(p.z() * 2 * PI + PI / 2.0); 
  return x2 + y2 + z2;
}

FT gyroid (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI); 
  return cx * sy + cy * sz + cz * sx;
  
}

FT diamond (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI); 
  return sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz;
}

FT double_p (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI), 
  c2y = std::cos(2 * p.y() * 2 * PI),
  c2z = std::cos(2 * p.z() * 2 * PI); 
  return 0.5 * (cx * cy  + cy * cz + cz * cx ) + 0.2*(c2x + c2y + c2z);
}

FT G_prime (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI), 
  c2y = std::cos(2 * p.y() * 2 * PI),
  c2z = std::cos(2 * p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI), 
  s2y = std::sin(2 * p.y() * 2 * PI),
  s2z = std::sin(2 * p.z() * 2 * PI);  
  return 5 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx) + 
  1*(c2x * c2y  + c2y * c2z  + c2z * c2x);
}

FT lidinoid (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI), 
  c2y = std::cos(2 * p.y() * 2 * PI),
  c2z = std::cos(2 * p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI), 
  s2y = std::sin(2 * p.y() * 2 * PI),
  s2z = std::sin(2 * p.z() * 2 * PI);  
  return 1 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx) - 
  1 * (c2x * c2y  + c2y * c2z  + c2z * c2x) + 0.3;
}

FT D_prime (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI), 
  c2y = std::cos(2 * p.y() * 2 * PI),
  c2z = std::cos(2 * p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI);
  return 1 * ( sx * sy * sz) + 1 * ( cx * cy * cz) -
  1 * ( c2x * c2y  + c2y * c2z  + c2z * c2x) - 0.4;
}

FT split_p (const Point& p) {
  const FT cx = std::cos(p.x() * 2 * PI), 
  cy = std::cos(p.y() * 2 * PI),
  cz = std::cos(p.z() * 2 * PI);
  const FT c2x = std::cos(2 * p.x() * 2 * PI), 
  c2y = std::cos(2 * p.y() * 2 * PI),
  c2z = std::cos(2 * p.z() * 2 * PI); 
  const FT sx = std::sin(p.x() * 2 * PI), 
  sy = std::sin(p.y() * 2 * PI),
  sz = std::sin(p.z() * 2 * PI);
  const FT s2x = std::sin(2 * p.x() * 2 * PI), 
  s2y = std::sin(2 * p.y() * 2 * PI),
  s2z = std::sin(2 * p.z() * 2 * PI); 
  return  1.1 * (s2x * sz * cy  + s2y * sx * cz  + s2z * sy * cx) 
  - 0.2 * (c2x * c2y  + c2y * c2z  + c2z * c2x) 
  - 0.4 * (cx + cy + cz);
}

// simple implementation
Point canonicalize_point(const Point& p) {
  const FT& xmin = 0., xmax = 1; 
  const FT& ymin = 0., ymax = 1;
  const FT& zmin = 0., zmax = 1;
  
  const FT xsize = xmax - xmin; 
  const FT ysize = ymax - ymin; 
  const FT zsize = zmax - zmin; 
  
  FT x = p.x(), y = p.y(), z = p.z();
  if(p.x() < xmin) {
    x += xsize; 
  } else if( !(p.x() < xmax) ) {
    x -= xsize;
  }
  
  if(p.y() < ymin) {
    y += ysize; 
  } else if( !(p.y() < ymax) ) {
    y -= ysize;
  }
  
  if(p.z() < zmin) {
    z += zsize; 
  } else if( !(p.z() < zmax) ) {
    z -= zsize;
  }
  
  assert( !(x < xmin) && (x < xmax) );
  assert( !(y < ymin) && (y < ymax) );
  assert( !(z < zmin) && (z < zmax) );
  return Point(x, y, z);
}

FT cylinder(const Point& p) {
  Point p_ = canonicalize_point( p );
  const FT x2 = (p_.y() - 0.5)*(p_.y() - 0.5), 
  y2=(p_.z() - 0.5)*(p_.z() - 0.5); 
  return x2+y2-0.1;
}

typedef CGAL::Mesh_constant_domain_field_3<Periodic_mesh_domain::R,
Periodic_mesh_domain::Index> Sizing_field;

int main()
{
  const int functions_count = 11;
  // Array of the functions
  Function* implicit_function[functions_count] =
  {
    &sphere,
    &cylinder,
    &scherk,
    &schwarz_p,
    &gyroid,
    &diamond,
    &double_p,
    &G_prime,
    &lidinoid,
    &D_prime,
    &split_p
  };
    
  for( int i = 0; i < functions_count; i++ ) {
    
    try {
    
      // Periodic mesh domain (Warning: Sphere_3 constructor uses squared radius !)
      Periodic_mesh_domain domain(*implicit_function[i], CGAL::Iso_cuboid_3<K>(0, 0, 0, 1, 1, 1));
      
      /*
      double kidney_size = 3.;
      int volume_dimension = 3;
      Sizing_field size(8);
      size.set_size(kidney_size, volume_dimension, 
                    domain.index_from_subdomain_index(2));
      
      size.set_size(0.05, volume_dimension, 
                    domain.index_from_subdomain_index(1));
       */
       
      // Mesh criteria
      Mesh_criteria criteria(domain, facet_angle=30, facet_size=0.05, facet_distance=0.025,
                             cell_radius_edge=2, cell_size = 0.05/*size*/);
      
      
      // Mesh generation
      C3t3 c3t3 = CGAL::make_periodic_mesh_3<C3t3>(domain, criteria, no_odt(), no_perturb(), no_lloyd());
      
      // File name
      std::stringstream index;
      index << i;
      std::string file_name = std::string("out") + index.str() + ".mesh"; 
      
      // Output
      std::ofstream medit_file(file_name.c_str());
      
      write_complex_to_medit(medit_file, c3t3);
    
      std::cout << i << ": passed" << std::endl;
    } catch(...) {
      std::cout << i << ": failed" << std::endl;
    }
    
  }
  
  return 0;
}

