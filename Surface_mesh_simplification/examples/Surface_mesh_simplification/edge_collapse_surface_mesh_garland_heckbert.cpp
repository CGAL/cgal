#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/Surface_mesh_simplification/GarlandHeckbert_edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_placement.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;


int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  std::ifstream is(argv[1]);
  is >> surface_mesh;
  if(!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.1);

  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh>::garland_heckbert_map_type map;
  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh> vis(map);

  int r = SMS::edge_collapse(surface_mesh, stop,
          CGAL::parameters::get_cost(SMS::GarlandHeckbert_cost <Surface_mesh>(map))
                           .get_placement(SMS::GarlandHeckbert_placement<Surface_mesh>(map))
                           .visitor(vis));

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
