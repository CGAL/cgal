#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <CGAL/Surface_mesh_simplification/GarlandHeckbert_edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_cost_stop_predicate.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  if(argc<3)
  {
    std::cerr << "Usage: " << argv[0] << " input.off minimal_quadric_error [out.off]\n";
    return EXIT_FAILURE;
  }

  Surface_mesh surface_mesh;

  std::ifstream is(argv[1]);
  assert(is);

  is >> surface_mesh;
  double threshold = atof(argv[2]);

  std::cout << num_vertices(surface_mesh) << " vertices and " << num_edges(surface_mesh) << " edges (input)" << std::endl;

  if(!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh>::garland_heckbert_state_type state;

  int r = SMS::edge_collapse(surface_mesh,
                             CGAL::Surface_mesh_simplification::GarlandHeckbert_cost_stop_predicate<double>(threshold),
                             CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, surface_mesh))
                                              .halfedge_index_map(get(CGAL::halfedge_external_index  ,surface_mesh))
                                              .get_cost(SMS::GarlandHeckbert_cost <Surface_mesh>(state))
                                              .get_placement(SMS::GarlandHeckbert_placement<Surface_mesh>(state))
                                              .visitor(SMS::GarlandHeckbert_edge_collapse_visitor_base<Surface_mesh>(state)));

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
