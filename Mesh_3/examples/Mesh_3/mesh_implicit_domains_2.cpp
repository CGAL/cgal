#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include "implicit_functions.h"

// IO
#include <CGAL/IO/File_medit.h>

using namespace CGAL::parameters;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
                                                        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;


int main()
{
  // Define functions
  Function f1(&torus_function);
  Function f2(&sphere_function<3>);

  Function_vector v;
  v.push_back(f1);
  v.push_back(f2);

  std::vector<std::string> vps;
  vps.push_back("+-");

  /// [Domain creation] (Warning: Sphere_3 constructor uses square radius !)
  namespace param = CGAL::parameters;
  Mesh_domain domain(param::function = Function_wrapper(v, vps),
                     param::bounding_object = K::Sphere_3(CGAL::ORIGIN,
                                                          5.*5.),
                     param::relative_error_bound = 1e-6);
  /// [Domain creation]

  // Set mesh criteria
  Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
  Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

  // Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
  CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);

  // Exudation
  CGAL::exude_mesh_3(c3t3,12);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::output_to_medit(medit_file, c3t3);

  return 0;
}
