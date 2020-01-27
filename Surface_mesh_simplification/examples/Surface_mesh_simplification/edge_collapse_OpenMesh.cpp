#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <iostream>
#include <fstream>

typedef OpenMesh::PolyMesh_ArrayKernelT</* default traits*/>     Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::edge_descriptor       edge_descriptor;

class Constrained_edge_map
{
public:
  typedef boost::read_write_property_map_tag    category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(Surface_mesh& sm)
    : sm_(sm)
  {
    sm_.add_property(constraint);
  }

  inline friend reference get(const Constrained_edge_map& em, key_type e)
  {
    bool b = em.sm_.property(em.constraint,em.sm_.edge_handle(e.idx()));
    return b;
  }

  inline friend void put(const Constrained_edge_map& em, key_type e, value_type b)
  {
    em.sm_.property(em.constraint,em.sm_.edge_handle(e.idx())) = b;
  }

private:
  Surface_mesh& sm_;
  OpenMesh::EPropHandleT<bool> constraint;
};

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  Constrained_edge_map constraints_map(surface_mesh);

  const char* filename = (argc > 1) ? argv[1] : "data/cube-meshed.off";
  OpenMesh::IO::read_mesh(surface_mesh, filename);

  if(!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // For the pupose of the example we mark 100 edges as constrained edges
  int count = 0;
  for(edge_descriptor e : edges(surface_mesh))
    put(constraints_map, e, (count++ < 100));

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  const std::size_t stop_n = (argc > 2) ? std::stoi(argv[2]) : 1000;
  SMS::Count_stop_predicate<Surface_mesh> stop(stop_n);

  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  std::cout << "Collapsing edges of mesh: " << filename << ", aiming for " << stop_n << " final edges..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh, stop,
                             CGAL::parameters::halfedge_index_map(get(CGAL::halfedge_index,surface_mesh))
                                              .vertex_point_map(get(boost::vertex_point, surface_mesh))
                                              .edge_is_constrained_map(constraints_map));

  surface_mesh.garbage_collection();
  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << num_edges(surface_mesh) << " final edges.\n";

  OpenMesh::IO::write_mesh(surface_mesh, "out.off");

  return EXIT_SUCCESS;
}
