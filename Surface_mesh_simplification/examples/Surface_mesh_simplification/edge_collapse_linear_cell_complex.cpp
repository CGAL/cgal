#include <iostream>
#include <fstream>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_incremental_builder_v2.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex.h>
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;

struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, Refs > Dart;
    typedef CGAL::Cell_attribute_with_point< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple<Vertex_attribute, void, Face_attribute> Attributes;
  };
};

typedef CGAL::Linear_cell_complex<2, 3, MyTraits, Myitem> LCC;
namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv )
{
  if (argc!=2)
  {
    std::cout<<"Usage: simplification_Linear_cell_complex inofffile [outofffile]"<<std::endl;
    return EXIT_FAILURE;
  }

  LCC lcc;
  std::ifstream is(argv[1]);
  CGAL::load_off_v2(lcc, is);

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
     ,CGAL::parameters::halfedge_index_map(get(CGAL::halfedge_external_index, lcc))
          .vertex_index_map(get(CGAL::vertex_external_index,lcc))
          .get_cost(SMS::Edge_length_cost<LCC>())
          .get_placement(SMS::Midpoint_placement<LCC>())
     );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (lcc.number_of_darts()/2) << " final edges.\n" ;

  lcc.display_characteristics(std::cout)<<", is_valid="<<CGAL::is_valid(lcc)<<std::endl;

  std::ofstream os(argc > 2 ? argv[2] : "out.off");
  CGAL::write_off(lcc, os);
  return EXIT_SUCCESS;
}
// EOF //
