#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Midpoint placement policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

//Placement wrapper
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;

typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;

// BGL property map which indicates whether an edge is marked as non-removable
struct Border_is_constrained_edge_map
{
  const Surface_mesh* sm_ptr;
  typedef edge_descriptor                                       key_type;
  typedef bool                                                  value_type;
  typedef value_type                                            reference;
  typedef boost::readable_property_map_tag                      category;

  Border_is_constrained_edge_map(const Surface_mesh& sm) : sm_ptr(&sm) {}

  friend bool get(Border_is_constrained_edge_map m, const key_type& edge) {
    return CGAL::is_border(edge, *m.sm_ptr);
  }
};

// Placement class
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                                   Border_is_constrained_edge_map > Placement;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/mesh_with_border.off";
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

  Surface_mesh::Property_map<halfedge_descriptor, std::pair<Point_3, Point_3> > constrained_halfedges;
  constrained_halfedges = surface_mesh.add_property_map<halfedge_descriptor,std::pair<Point_3, Point_3> >("h:vertices").first;

  std::size_t nb_border_edges=0;
  for(halfedge_descriptor hd : halfedges(surface_mesh))
  {
    if(CGAL::is_border(hd, surface_mesh))
    {
      constrained_halfedges[hd] = std::make_pair(surface_mesh.point(source(hd, surface_mesh)),
                                                 surface_mesh.point(target(hd, surface_mesh)));
      ++nb_border_edges;
    }
  }

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);
  Border_is_constrained_edge_map bem(surface_mesh);

  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  std::cout << "Collapsing as many edges of mesh: " << filename << " as possible..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::edge_is_constrained_map(bem)
                                              .get_placement(Placement(bem)));

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n";

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  os.precision(17);
  os << surface_mesh;

  // now check!
  for(halfedge_descriptor hd : halfedges(surface_mesh))
  {
    if(CGAL::is_border(hd,surface_mesh))
    {
      --nb_border_edges;
      if(constrained_halfedges[hd] != std::make_pair(surface_mesh.point(source(hd, surface_mesh)),
                                                     surface_mesh.point(target(hd, surface_mesh))))
      {
        std::cerr << "oops. send us a bug report\n";
      }
    }
  }
  assert(nb_border_edges==0);

  return EXIT_SUCCESS;
}
