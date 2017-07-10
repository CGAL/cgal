#include <iostream>
#include <fstream>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, MyTraits>::type LCC;
namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv )
{
  if (argc<2 || argc>3)
  {
    std::cout<<"Usage: simplification_Linear_cell_complex inofffile [outofffile]"<<std::endl;
    return EXIT_FAILURE;
  }

  LCC lcc;
  CGAL::read_off(argv[1], lcc);

  lcc.display_characteristics(std::cout)<<", is_valid="<<CGAL::is_valid(lcc)<<std::endl;

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number (1000)
  SMS::Count_stop_predicate<LCC> stop(1000);
  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  int r = SMS::edge_collapse
    (lcc
     ,stop
     ,CGAL::parameters::halfedge_index_map(get(CGAL::halfedge_index, lcc))
          .vertex_index_map(get(boost::vertex_index, lcc))
          .get_cost(SMS::Edge_length_cost<LCC>())
          .get_placement(SMS::Midpoint_placement<LCC>())
     );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (lcc.number_of_darts()/2) << " final edges.\n" ;

  lcc.display_characteristics(std::cout)<<", is_valid="<<CGAL::is_valid(lcc)<<std::endl;

  CGAL::write_off((argc > 2 ? argv[2] : "out.off"), lcc);
  return EXIT_SUCCESS;
}
// EOF //
