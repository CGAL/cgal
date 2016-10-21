#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv )
{
  if (argc<3)
  {
    std::cerr << "Usage: " << argv[0] << " input.off minimal_edge_length [out.off]\n";
    return 1;
  }

  Surface_mesh surface_mesh;

  std::ifstream is(argv[1]) ; is >> surface_mesh ;
  double threshold = atof(argv[2]);

  int r = SMS::edge_collapse
            (surface_mesh
             , CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double>(threshold)
             , CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh))
                               .halfedge_index_map  (get(CGAL::halfedge_external_index  ,surface_mesh))
                               .get_cost (SMS::Edge_length_cost <Surface_mesh>())
                               .get_placement(SMS::Midpoint_placement<Surface_mesh>())
            );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;

  std::ofstream os( argc > 3 ? argv[3] : "out.off" ) ; os << surface_mesh ;

  return 0 ;
}
