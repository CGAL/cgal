#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                              Kernel;
typedef Kernel::Point_3                                             Point;

// Setup an enriched polyhedron type which stores an id() field in the items
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor        vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor      halfedge_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/cube.off";
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

  // The items in this polyhedron have an "id()" field
  // which the default index maps used in the algorithm
  // need to get the index of a vertex/edge.
  // However, the Polyhedron_3 class doesn't assign any value to
  // this id(), so we must do it here:
  int index = 0;
  for(halfedge_descriptor hd : halfedges(surface_mesh))
    hd->id() = index++;
  index = 0;

  for(vertex_descriptor vd : vertices(surface_mesh))
    vd->id() = index++;

  // In this example, the simplification stops when the number of undirected edges
  // drops below xx% of the initial count
  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  std::cout << "Collapsing edges of mesh: " << filename << ", aiming for " << 100 * ratio << "% of the input edges..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh, stop);

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n";

  std::ofstream os((argc > 3) ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
