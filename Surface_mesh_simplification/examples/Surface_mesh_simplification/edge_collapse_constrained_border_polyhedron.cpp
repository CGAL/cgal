#include <iostream>
#include <fstream>
#include <map>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

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
typedef CGAL::Polyhedron_3<Kernel> Surface;

namespace SMS = CGAL::Surface_mesh_simplification ;

//
// BGL property map which indicates whether an edge is marked as non-removable
//
struct Border_is_constrained_edge_map{
  typedef boost::graph_traits<Surface>::edge_descriptor key_type;
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  friend bool get(Border_is_constrained_edge_map, key_type edge) {
    return edge->is_border_edge();
  }
};

//
// Placement class
//
typedef SMS::Constrained_placement<SMS::Midpoint_placement<Surface>,
                                   Border_is_constrained_edge_map > Placement;

int main( int argc, char** argv )
{
  Surface surface;

  if (argc!=2){
    std::cerr<< "Usage: " << argv[0] << " input.off\n";
    return 1;
  }

  std::ifstream is(argv[1]);
  if(!is){
    std::cerr<< "Filename provided is invalid\n";
    return 1;
  }

  is >> surface ;

  // map used to check that constrained_edges and the points of its vertices
  // are preserved at the end of the simplification
  std::map<Surface::Halfedge_handle,std::pair<Point_3, Point_3> >constrained_edges;
  std::size_t nb_border_edges=0;

  for (Surface::Halfedge_iterator hit=surface.halfedges_begin(),
                                  hit_end=surface.halfedges_end(); hit!=hit_end;
                                  ++hit )
  {
    if ( hit->is_border() ){
      constrained_edges[hit]=std::make_pair( hit->opposite()->vertex()->point(),
                                             hit->vertex()->point() );
      ++nb_border_edges;
    }
  }

  // Contract the surface as much as possible
  SMS::Count_stop_predicate<Surface> stop(0);

  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface lack an "id()" field.
  int r = SMS::edge_collapse
            (surface
            ,stop
            ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,surface))
                  .edge_index_map  (boost::get(CGAL::edge_external_index  ,surface))
                  .edge_is_constrained_map(Border_is_constrained_edge_map())
                  .get_placement(Placement())
            );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;

  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;

  // now check!
  for (Surface::Halfedge_iterator hit=surface.halfedges_begin(),
                                  hit_end=surface.halfedges_end(); hit!=hit_end;
                                  ++hit )
  {
    if (hit->is_border()){
      --nb_border_edges;
      assert( constrained_edges[hit] ==
              std::make_pair( hit->opposite()->vertex()->point(),
                              hit->vertex()->point() ) );
    }
  }
  assert( nb_border_edges==0 );

  return 0 ;
}
