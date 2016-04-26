#include <fstream>
#include <iostream>

template<class Polyhedron>
bool read_to_polyhedron(const char* file_name, Polyhedron& mesh)
{
  std::ifstream input(file_name);
    
  bool ok = true;
  if(!input)
   ok = false;
  else if(!(input>>mesh))
   ok = false;
  else if(mesh.empty())
   ok = false;
  if ( !ok ){
    std::cerr << "Problem occured while reading off file";
    return false;
  }
  return true;
}
