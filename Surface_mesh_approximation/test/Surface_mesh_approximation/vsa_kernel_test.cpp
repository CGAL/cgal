#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Simple_cartesian<double> Sckernel;

template <typename K, typename TM>
int test() {
  TM tm;
  std::ifstream input("./data/cube_meshed.off");
  if (!input || !(input >> tm) || num_vertices(tm) == 0) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Polyhedron_3<K> Polyhedron;
  Polyhedron out_mesh;
  CGAL::vsa_mesh_approximation(tm, out_mesh,
    CGAL::VSA::parameters::seeding_by_number(6).
      iterations(30).
      inner_iterations(5).
      chord_subdivide(0.5));

  return EXIT_SUCCESS;
}

/**
 * This file tests the VSA free function API with different configuration of kernel and surface mesh.
 */
int main()
{
  if (test<Epic, CGAL::Polyhedron_3<Epic> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test<Epic, CGAL::Surface_mesh<Epic::Point_3> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test<Sckernel, CGAL::Polyhedron_3<Sckernel> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  if (test<Sckernel, CGAL::Surface_mesh<Sckernel::Point_3> >() == EXIT_FAILURE)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
