
//******************************************************************************
// File Description :
// Outputs to out.mesh a mesh of implicit domains.
//******************************************************************************



#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

// IO
#include <CGAL/IO/File_medit.h>

using namespace CGAL::parameters;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::FT FT;
typedef FT (*Function)(const Point&);
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
                                                        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;
typedef CGAL::Mesh_domain_with_polyline_features_3<Mesh_domain> Mesh_domain_with_features;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;


double cube_function_1 (const Point& p)
{
  if( p.x() >= 0 && p.x() <= 2 &&
      p.y() >= 0 && p.y() <= 2 &&
      p.z() >= 0 && p.z() <= 2 )
    return -1.;
  return 1.;
}

double cube_function_2 (const Point& p)
{
  if( p.x() >= 1 && p.x() <= 3 &&
      p.y() >= 1 && p.y() <= 3 &&
      p.z() >= 1 && p.z() <= 3 )
    return -1.;
  return 1.;
}

// Polyline
typedef std::vector<Point>  Polyline_3;
typedef std::list<Polyline_3> Polylines;

void create_polylines (Polylines& polylines)
{
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,1,1));
    polyline.push_back(Point(2,1,1));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,1,1));
    polyline.push_back(Point(2,2,1));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,2,1));
    polyline.push_back(Point(1,2,1));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,2,1));
    polyline.push_back(Point(1,1,1));
    polylines.push_back(polyline);
  }
//----------
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,1,2));
    polyline.push_back(Point(2,1,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,1,2));
    polyline.push_back(Point(2,2,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,2,2));
    polyline.push_back(Point(1,2,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,2,2));
    polyline.push_back(Point(1,1,2));
    polylines.push_back(polyline);
  }
  //----------
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,1,1));
    polyline.push_back(Point(1,1,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,1,1));
    polyline.push_back(Point(2,1,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(2,2,1));
    polyline.push_back(Point(2,2,2));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point(1,2,1));
    polyline.push_back(Point(1,2,2));
    polylines.push_back(polyline);
  }
}

int main()
{
  // Define functions
  Function f1 = cube_function_1;
  Function f2 = cube_function_2;

  Function_vector v;
  v.push_back(f1);
  v.push_back(f2);

  std::vector<std::string> vps;
  vps.push_back("--");

  // Domain (Warning: Sphere_3 constructor uses square radius !)
  Mesh_domain_with_features domain(Function_wrapper(v, vps), K::Sphere_3(CGAL::ORIGIN, 5.*5.));
  Polylines polylines;
  create_polylines(polylines);
  domain.add_features(polylines.begin(),polylines.end());

  // Set mesh criteria
  Mesh_criteria criteria(edge_size = 0.15,
      facet_angle = 30, facet_size = 0.2,
      cell_radius_edge_ratio = 2, cell_size = 0.4);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
  CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
  
  // Exudation
  CGAL::exude_mesh_3(c3t3,12);
  
  // Output
  std::ofstream medit_file("out_cubes_intersection_with_features.mesh");
  CGAL::output_to_medit(medit_file, c3t3);

  return 0;
}
