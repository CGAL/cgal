#include <fstream>
#include <iostream>

template<class Polyhedron>
bool read_to_polyhedron(const char* file_name, Polyhedron& mesh)
{
  std::ifstream input(file_name);

  if ( !input || !(input >> mesh) || mesh.empty() ){
    std::cerr << "Problem occured while reading off file";
    return false;
  }
  return true;
}
