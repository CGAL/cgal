//
// First a bunch of general support headers, used in all the examples.
//
#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Unique_hash_map.h>

//
// Then the package-specific headers:
//

// 1.The header defining how a concrete surface type is adapted to the TSM concept.
#include <CGAL/Surface_mesh_simplification/Polyhedron.h>

// 2.The header defining the algorithm.
#include <CGAL/Surface_mesh_simplification/Edge_collapse.h>

// 3.The header defining the stop condition policy (which has no default)
// 4.Any additional header optionally needed by a non-default policy.
// Each example includes these policy headers between these marks:

// === EXAMPLE SPECIFIC HEADERS BEGINS HERE ===
// === EXAMPLE SPECIFIC HEADERS ENDS HERE ===


//
// Instantiation of the surface type -- same in all these examples
//
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Surface; 

// === EXAMPLE SPECIFIC DETAILS BEGINS HERE ===
// === EXAMPLE SPECIFIC DETAILS ENDS HERE ===

// All the classes and functions used in these examples are found in the following
// nested namespace, but a namespace alias is used for brevity.
namespace TSMS = CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

//
// main function. All the examples share the following structure.
// The concrete details relevant to each example are placed 
// between the marks inside main():
//
int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  // === CONCRETE USAGE EXAMPLE BEGINS HERE ===
  
  // === CONCRETE USAGE EXAMPLE ENDS HERE ===

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

// EOF //
