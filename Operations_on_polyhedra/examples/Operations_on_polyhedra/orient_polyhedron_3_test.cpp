#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/orient_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;


void test(const char* file_name) {

  std::ifstream input(file_name);
  Polyhedron oriented_poly; // file should contain oriented poly
  if ( !input || !(input >> oriented_poly) || oriented_poly.empty() ){
    std::cerr << "Error: can not read file: " << file_name;
    assert(false);
  }

  assert(CGAL::is_oriented(oriented_poly));
  oriented_poly.inside_out();
  assert(!CGAL::is_oriented(oriented_poly));

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
