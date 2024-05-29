#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

int main() {

  std::string filename("data/issue_8213.off");

  Surface_mesh surface_mesh;

  std::ifstream in(filename);
  in >> surface_mesh;

  const size_t target_number_of_faces = 3;
  SMS::Face_count_stop_predicate<Surface_mesh> stop(target_number_of_faces);

  std::cout << "Input mesh number of faces: " << surface_mesh.number_of_faces() << ", target number of faces: " << target_number_of_faces << std::endl;

  SMS::edge_collapse(
    surface_mesh,
    stop,
    CGAL::parameters::
       filter(SMS::Bounded_normal_change_filter<>())
      .get_cost(SMS::LindstromTurk_cost<Surface_mesh>())
      .get_placement(SMS::LindstromTurk_placement<Surface_mesh>())
    );

  std::cout.precision(17);
  std::cout <<  surface_mesh << std::endl;
  return EXIT_SUCCESS;
}
