#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

typedef Kernel::Point_3 Point ;

typedef Surface::Halfedge_handle Halfedge_handle ;


// === EXAMPLE SPECIFIC DETAILS ENDS HERE ===

namespace SMS = CGAL::Surface_mesh_simplification ;

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  Visitor vis ;

  int r = SMS::edge_collapse(surface
                             ,SMS::Count_ratio_stop_predicate<Surface>(0.10) 
                             ,CGAL::edge_index_map(boost::get(CGAL::edge_external_index,surface))
                             .visitor(&vis)
                             );

  
  std::cout << "\nEdges collected: " << vis.collected
            << "\nEdges proccessed: " << vis.processed
            << "\nEdges collapsed: " << vis.collapsed
            << std::endl
            << "\nEdges not collapsed due to topological constrians: " 
            << vis.non_collapsable
            << "\nEdge not collapsed due to cost computation constrians: " 
            << vis.cost_uncomputable 
            << "\nEdge not collapsed due to placement computation constrians: " 
            << vis.placement_uncomputable 
            << std::endl ; 
            
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===
  
  std::cout << "\nFinished...\n" << r << " edges removed.\n"  
            << (surface.size_of_halfedges()/2) << " final edges." ;
  
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
