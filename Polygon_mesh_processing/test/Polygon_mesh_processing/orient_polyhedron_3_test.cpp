#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/orient_polygon_mesh.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;


void test(const char* file_name)
{
  std::ifstream input(file_name);
  Polyhedron poly; // file should contain oriented polyhedron
  
  if ( !input || !(input >> poly) || poly.empty() )
  {
    std::cerr << "Error: can not read file: " << file_name;
    CGAL_assertion(false);
  }

  bool before = CGAL::Polygon_mesh_processing::is_outward_oriented(poly);
  CGAL_assertion(before);

  poly.inside_out();

  bool after = CGAL::Polygon_mesh_processing::is_outward_oriented(poly);
  CGAL_assertion(!after);

  std::cerr << file_name << " passed the test." << std::endl;
}

int main(int argc, char** argv) {
  //files.push_back("data/elephant.off");
  //files.push_back("data/camel.off");

  for(int i=1;i<argc;++i) {
    test(argv[i]);
  }
  std::cerr << "All done." << std::endl;
}
