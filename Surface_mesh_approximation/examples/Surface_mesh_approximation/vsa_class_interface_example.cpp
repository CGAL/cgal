#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/VSA_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;
typedef CGAL::VSA_approximation<Polyhedron, Vertex_point_map> Mesh_approximation;

// L21 error metric 
typedef Mesh_approximation::Error_metric Metric;
typedef Mesh_approximation::Proxy_fitting Proxy_fitting;

int main()
{
  // read input surface triangle mesh
  Polyhedron input;
  std::ifstream file("data/mask.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // create VSA algorithm instance
  Mesh_approximation approx(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));

  // set error and fitting functors
  Metric metric(input);
  Proxy_fitting proxy_fitting(input);
  approx.set_metric(metric, proxy_fitting);

  // seeding 100 random proxies
  approx.seeding(CGAL::Random, 100);
  
  // run 30 iterations 
  approx.run(30);

  // add 3 proxies to the one with the maximum fitting error
  // run 5 iterations between each addition
  approx.add_proxies_furthest(3, 5);

  // run 10 iterations
  approx.run(10);

  // teleport 2 proxies to tunnel out of local minima
  // run 5 iterations between each teleport
  approx.teleport_proxies(2, 5);
  approx.run(10);

  // meshing with default parameters
  approx.extract_mesh();

  return EXIT_SUCCESS;
}
