#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;
struct Dummy_placement {

  template <typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile&) const
  {
    return boost::none;
  }

 template <typename Profile>
 boost::optional<typename Profile::Point> operator()(const Profile&, const boost::optional<typename Profile::Point>& op) const
  {
    return op;
  }


};
int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fold.off");
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number
  const std::size_t stop_n = (argc > 2) ? std::stoi(argv[2]) : num_halfedges(surface_mesh)/2 - 1;
  SMS::Count_stop_predicate<Surface_mesh> stop(stop_n);

  typedef SMS::LindstromTurk_placement<Surface_mesh> Placement;

  CGAL::Timer t;
  t.start();
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  std::cout << "Collapsing edges of mesh: " << filename << ", aiming for " << stop_n << " final edges..." << std::endl;
  SMS::Bounded_normal_change_filter<> filter;
  SMS::edge_collapse(surface_mesh, stop,
                     CGAL::parameters::get_cost(SMS::LindstromTurk_cost<Surface_mesh>())
                                      .filter(filter)
                                      .get_placement(Placement()));

  std::cout << t.time() << " sec" << std::endl;
  CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out.off", surface_mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
