#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

namespace SMS = CGAL::Surface_mesh_simplification ;

//
// BGL property map which indicates whether an edge is border OR is marked as non-removable
//
class Constrains_map : public boost::put_get_helper<bool,Constrains_map>
{
public:

  typedef boost::readable_property_map_tag                    category;
  typedef bool                                                value_type;
  typedef bool                                                reference;
  typedef boost::graph_traits<Surface const>::edge_descriptor key_type;

  Constrains_map() : mConstrains(false) {}

  reference operator[](key_type const& e) const { return e->is_border() || is_constrained(e) ; }
  
  void set_is_constrained ( key_type const& e, bool is )  { mConstrains[e]=is; }
  
  bool is_constrained( key_type const& e ) const { return mConstrains.is_defined(e) ? mConstrains[e] : false ; }
  
private:
  
  CGAL::Unique_hash_map<key_type,bool> mConstrains ;
  
};

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface drops below the specified number (1000)
  SMS::Count_stop_predicate<Surface> stop(10);
     
  Constrains_map constrains_map ;
     
     
  // This example marks ALL edges as non-removable, but a real world application would mark only selected ones.   
  for( Surface::Halfedge_iterator eb = surface.halfedges_begin(), ee = surface.halfedges_end() ; eb != ee ; ++ eb ) 
    constrains_map.set_is_constrained(eb,true);
     
  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface lack an "id()" field.
  int r = SMS::edge_collapse
            (surface
            ,stop
            ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,surface)) 
                  .edge_index_map  (boost::get(CGAL::edge_external_index  ,surface)) 
                  .edge_is_border_map(constrains_map)
            );
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}
