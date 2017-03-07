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
typedef K::Segment_3 Segment;
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

const FT& PI = 3.14159265358979;

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

FT sphere(const Point& p) {
  Point p_ = canonicalize_point( p );
  const FT x2 = (p_.x() - 0.5) * (p_.x() - 0.5), 
    y2 = (p_.y() - 0.5) * (p_.y() - 0.5), 
    z2 = (p_.z() - 0.5) * (p_.z() - 0.5);
  return x2 + y2 + z2 - 0.2;
}

FT cylinder(const Point& p) {
  Point p_ = canonicalize_point( p );
  const FT x2 = (p_.y() - 0.5)*(p_.y() - 0.5), 
    y2=(p_.z() - 0.5)*(p_.z() - 0.5); 
  return x2+y2-0.1;
}

// user function

struct Segments_function : public std::unary_function<Point, FT>
{
  typedef std::vector<Segment> Segments;

  Segments_function() : segments()
  {
    const Point pmid(0.5,0.5,0.5);
    nb_evals = 0;
    segments.push_back( Segment(Point(0,0.5,0), pmid) );
    segments.push_back( Segment(Point(1,0.5,0), pmid) );
    segments.push_back( Segment(Point(0,0.5,1), pmid) );
    segments.push_back( Segment(Point(1,0.5,1), pmid) );
    segments.push_back( Segment(Point(0.5,0,0), pmid) );
    segments.push_back( Segment(Point(0.5,1,0), pmid) );
    segments.push_back( Segment(Point(0.5,0,1), pmid) );
    segments.push_back( Segment(Point(0.5,1,1), pmid) );
  }

  FT operator()(const Point& p_)
  {
    Point p = canonicalize_point(p_);
    ++nb_evals;

    FT min_distance = 1000000;
    for (Segments::const_iterator si = segments.begin(); 
      si != segments.end(); ++si)
      min_distance = std::min(CGAL::squared_distance(p, *si), min_distance);

    return min_distance - 0.01; // Change the squared beam radius here
  }

  Segments segments;
  int nb_evals;
};

FT segments(const Point& p) {
  return Segments_function()(p);
}

int main() {

  const int functions_count = 12;
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
    &split_p,
    &segments
  };

  Iso_cuboid box(0, 0, 0, 1, 1, 1);

  // Mesh criteria
  Periodic_facet_criteria facet_criteria(box, 30, 0.03, 0.03);
  Periodic_cell_criteria cell_criteria(box, 2, 0.05);

  Mesh_criteria criteria(facet_criteria, cell_criteria);

  for( int i = 0; i < functions_count; i++ ) {
    
    try {
      // Periodic mesh domain (Warning: Sphere_3 constructor uses squared radius !)
      Periodic_mesh_domain domain(*implicit_function[i], box);
      
      // Mesh generation
      C3t3 c3t3 = CGAL::make_periodic_mesh_3<C3t3>(domain, criteria);
      
      // File name
      std::stringstream index;
      index << i;
      std::string file_name = std::string("out") + index.str() + ".mesh"; 
      
      // Output
      std::ofstream medit_file(file_name.c_str());
      //c3t3.output_to_medit(medit_file);
      write_complex_to_medit_8(medit_file, c3t3);
      
      std::cout << i << ": passed" << std::endl;
    
    } catch (...) {
      std::cout << i << ": failed" << std::endl;
    }
    
  }
}
