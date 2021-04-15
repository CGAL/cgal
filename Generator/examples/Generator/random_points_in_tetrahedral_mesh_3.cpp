#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
using namespace CGAL;

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Polyhedron_3<K>                                     Polyhedron;
// Domain
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>             Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<
  Mesh_domain,CGAL::Default,Concurrency_tag>::type                Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>               C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>                                 Mesh_criteria;

// (Unweighted) point type
typedef Tr::Bare_point                                            Point;

int main()
{
 // Generated points are in that vector
  std::vector<Point> points;
  // Create input polyhedron
  Polyhedron polyhedron;
  polyhedron.make_tetrahedron(Point(-1,0,0), Point(0,1,0), Point(1,0,0), Point(0,0,-1));
  // Create domain
  Mesh_domain domain(polyhedron);
   using namespace CGAL::parameters;
  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
                         cell_radius_edge_ratio=3);
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Create the generator, input is the C3t3 c3t3
  Random_points_in_tetrahedral_mesh_3<C3t3> g(c3t3);
  // Get 100 random points in cdt
  std::copy_n( g, 100, std::back_inserter(points));

  // Check that we have really created 100 points.
  assert( points.size() == 100);

  // test if the generated points are close
  std::cout << points[0] << std::endl;

  return 0;
}



