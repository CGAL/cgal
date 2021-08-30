#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

typedef std::pair<Point_3, Vector_3> Pwn;
typedef CGAL::First_of_pair_property_map<Pwn> Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn> Vector_map;

typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson;

typedef CGAL::Labeled_mesh_domain_3<Kernel> Implicit_domain;
typedef CGAL::Mesh_triangulation_3<Implicit_domain, CGAL::Default,
                                   CGAL::Parallel_tag>::type Mesh_tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Mesh_tr> C3t3;
typedef CGAL::Mesh_criteria_3<Mesh_tr> Mesh_criteria;

int main(int, char**)
{
  std::vector<Pwn> points;

  std::ifstream stream("data/oni.pwn");
  if (!stream ||
      !CGAL::IO::read_XYZ
      (stream, std::back_inserter(points),
       CGAL::parameters::
       point_map(Point_map()).
       normal_map(Vector_map())))
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }

  Poisson poisson (points.begin(), points.end(), Point_map(), Vector_map());
  if (!poisson.compute_implicit_function())
  {
    std::cerr << "Error: cannot compute implicit function" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Bbox_3 bbox
    = CGAL::bbox_3 (boost::make_transform_iterator
                    (points.begin(),
                     CGAL::Property_map_to_unary_function<Point_map>()),
                    boost::make_transform_iterator
                    (points.end(),
                     CGAL::Property_map_to_unary_function<Point_map>()));

  Implicit_domain domain
    = Implicit_domain::create_implicit_mesh_domain
    (poisson, bbox);

  Mesh_criteria criteria (CGAL::parameters::facet_angle = 30,
                          CGAL::parameters::facet_size = 4,
                          CGAL::parameters::facet_distance = 0.1);
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3> (domain, criteria,
                                       CGAL::parameters::manifold());

  return EXIT_SUCCESS;
}
