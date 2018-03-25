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
typedef Mesh_approximation::Error_metric L21_metric;

int main()
{
  // reads input surface triangle mesh
  Polyhedron input;
  std::ifstream file("data/mask.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // creates VSA algorithm instance
  Mesh_approximation approx(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));

  // sets error and fitting functors
  L21_metric metric(input,
    get(boost::vertex_point, const_cast<Polyhedron &>(input)));
  approx.set_metric(metric);

  // seeds 100 random proxies
  approx.seeding(CGAL::Random, 100);
  
  // runs 30 iterations 
  approx.run(30);

  // adds 3 proxies to the one with the maximum fitting error,
  // running 5 iterations between each addition
  approx.add_to_furthest_proxies(3, 5);

  // runs 10 iterations
  approx.run(10);

  // teleports 2 proxies to tunnel out of local minima, 
  // running 5 iterations between each teleport
  approx.teleport_proxies(2, 5);

  // runs 10 iterations
  approx.run(10);

  // generates output mesh with default parameters
  approx.extract_mesh(CGAL::Surface_mesh_approximation::parameters::all_default());

  return EXIT_SUCCESS;
}
