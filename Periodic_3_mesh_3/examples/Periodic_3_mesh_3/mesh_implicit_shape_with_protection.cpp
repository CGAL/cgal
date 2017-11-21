#define CGAL_MESH_3_VERBOSE
//#define CGAL_MESHES_DEBUG_REFINEMENT_POINTS
//#define CGAL_MESH_3_DEBUG_FACET_CRITERIA
//#define CGAL_MESH_3_DEBUG_CELL_CRITERIA
#define CGAL_MESH_3_PROFILING
//#define CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
#define CGAL_MESH_3_PROTECTION_DEBUG 1111

#include <CGAL/Mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Medit_IO.h>
#include <CGAL/Implicit_periodic_3_mesh_domain_3.h>
#include <CGAL/make_periodic_3_mesh_3.h>
#include <CGAL/Periodic_3_mesh_triangulation_3.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_domain_with_polyline_features_3.h>

#include <algorithm>
#include <cmath>

namespace P3M3_IO = CGAL::Periodic_3_mesh_3::IO;

// Kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;

// Domain
typedef K::FT                                                       FT;
typedef K::Point_3                                                  Point;
typedef FT (Function)(const Point&);
typedef CGAL::Mesh_domain_with_polyline_features_3<
          CGAL::Implicit_periodic_3_mesh_domain_3<Function,K> >     Mesh_domain;

// Polyline
typedef std::vector<Point>                                          Polyline_3;
typedef std::list<Polyline_3>                                       Polylines;

// Triangulation
typedef CGAL::Periodic_3_mesh_triangulation_3<Mesh_domain>::type    Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<
          Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>  C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                                   Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// A pythagoras triple
static const FT r = 0.6;

static const FT mx = 0.26;
static const FT Mx = 0.74;
static const FT my = 0.26;
static const FT My = 0.74;

// Function
static Tr dummy_tr; // default constructs a iso_cuboid (0,0,0,1,1,1)

Point canonicalize_point(const Point& p)
{
  Point canonical_p = dummy_tr.robust_canonicalize_point(p);

  CGAL_postcondition( !(canonical_p.x() < 0) && (canonical_p.x() < 1) );
  CGAL_postcondition( !(canonical_p.y() < 0) && (canonical_p.y() < 1) );
  CGAL_postcondition( !(canonical_p.z() < 0) && (canonical_p.z() < 1) );
  return canonical_p;
}

FT sphere_function(const Point& p)
{
  return CGAL::squared_distance(canonicalize_point(p), Point(0.5, 0.5, 0.5)) - r*r;
}

void fill_sphere_polylines(Polylines& polylines)
{
  Polyline_3 polyline;

  for(int i = 0; i < 360; ++i)
  {
    Point p(0., 0.5 + 0.33 * std::cos(i*CGAL_PI/180), 0.5 + 0.33 * std::sin(i*CGAL_PI/180));
    std::cout << "Adding " << p << " to polyline" << std::endl;
    polyline.push_back(p);
  }
  polyline.push_back(polyline.front()); // close the line

  polylines.push_back(polyline);
}

FT squary_cylinder_function(const Point& p)
{
  if(p.x() < mx || p.x() > Mx ) return 1;
  if(p.y() < my || p.y() > My ) return 1;
  return -1;
}

void fill_squary_cylinder_polylines(Polylines& polylines, const FT edge_length)
{
  Polyline_3 polyline;

  polyline.push_back(Point(mx, mx, 0.));
  polyline.push_back(Point(mx, mx, 1. - 0.5 * edge_length));
  polylines.push_back(polyline);

  polyline.clear();
  polyline.push_back(Point(mx, my, 0.));
  polyline.push_back(Point(Mx, my, 0.));
  polyline.push_back(Point(Mx, My, 0.));
  polyline.push_back(Point(mx, My, 0.));
  polyline.push_back(Point(mx, my, 0.));
  polylines.push_back(polyline);
}

int main()
{
  int domain_size = 1;

  // Domain (Warning: Sphere_3 constructor uses squared radius !)
  Mesh_domain domain(squary_cylinder_function,
                     CGAL::Iso_cuboid_3<K>(0, 0, 0, domain_size, domain_size, domain_size));

  // Mesh criteria
  FT edge_length = 5 * r * sin(CGAL_PI/180.);
  const FT max_size = 0.125 * FT(domain_size) - 1e-10;
  edge_length = (std::min)(edge_length, max_size);

  Mesh_criteria criteria(edge_size = edge_length,
                         facet_angle = 0.05 * domain_size,
                         facet_size = 0.02 * domain_size,
                         cell_radius_edge_ratio = 2, cell_size = 0.5);

  // Create edge that we want to preserve
  Polylines polylines;

  // fill_sphere_polylines(polylines);
  fill_squary_cylinder_polylines(polylines, edge_length);

  // Insert edge in domain
  domain.add_features(polylines.begin(), polylines.end());

  // Mesh generation with feature preservation
  C3t3 c3t3_bis = CGAL::make_periodic_3_mesh_3<C3t3>(domain, criteria
                                                     //,no_exude(), no_perturb()
                                                     );

  // Output
  std::ofstream medit_file_bis("towers.mesh");
  P3M3_IO::write_complex_to_medit(medit_file_bis, c3t3_bis);

  return 0;
}
