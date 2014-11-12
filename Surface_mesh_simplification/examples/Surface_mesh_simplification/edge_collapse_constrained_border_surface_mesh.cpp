#include <iostream>
#include <fstream>
#include <map>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>


// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Midpoint placement policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

//Placement wrapper
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::edge_iterator edge_iterator;

namespace SMS = CGAL::Surface_mesh_simplification ;

//
// BGL property map which indicates whether an edge is marked as non-removable
//
struct Border_is_constrained_edge_map{
  const Surface_mesh* sm_ptr;
  typedef edge_descriptor key_type;
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  Border_is_constrained_edge_map()
  {}

  Border_is_constrained_edge_map(const Surface_mesh& sm)
    : sm_ptr(&sm)
  {}

  friend bool get(Border_is_constrained_edge_map m, const key_type& edge) {
    return  m.sm_ptr->is_border(edge);
  }
};

//
// Placement class
//
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                                   Border_is_constrained_edge_map > Placement;

int main( int argc, char** argv )
{
  Surface_mesh surface_mesh;

  if (argc!=2){
    std::cerr<< "Usage: " << argv[0] << " input.off\n";
    return 1;
  }

  std::ifstream is(argv[1]);
  if(!is){
    std::cerr<< "Filename provided is invalid\n";
    return 1;
  }

  is >> surface_mesh  ;

  
  Surface_mesh::Property_map<edge_descriptor,std::pair<Point_3, Point_3> > constrained_edges;

  constrained_edges = surface_mesh.add_property_map<edge_descriptor,std::pair<Point_3, Point_3>>("e:vertices").first;

  std::size_t nb_border_edges=0;
  BOOST_FOREACH(edge_descriptor ed, surface_mesh.edges()){
    if(surface_mesh.is_border(ed)){
      constrained_edges[ed] = std::make_pair(surface_mesh.point(source(ed,surface_mesh)),
                                               surface_mesh.point(target(ed,surface_mesh)));
      ++nb_border_edges;
    }
  }

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);

  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  int r = SMS::edge_collapse
            (surface_mesh
            ,stop
            ,CGAL::edge_is_constrained_map(Border_is_constrained_edge_map(surface_mesh))
                  .get_placement(Placement())
            );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << surface_mesh.number_of_edges() << " final edges.\n" ;

  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;


  return 0 ;
}
