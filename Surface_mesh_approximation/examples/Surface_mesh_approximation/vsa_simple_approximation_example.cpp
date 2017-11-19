#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // create polyhedral surface and read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;


  // The output will be an indexed triangle mesh
  std::vector<Kernel::Point_3> points;
  std::vector<std::vector<std::size_t> > triangles;
 
  // free function interface with named parameters
  CGAL::VSA::mesh_approximation(input,
                                std::back_inserter(points),
                                std::back_inserter(triangles),
                                CGAL::VSA::parameters::max_nb_proxies(200));

  std::cout << "#vertices: " << points.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
