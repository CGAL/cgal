#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include "draw_facewidth.h"

using LCC_3            =CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using CST              =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
using Dart_const_handle=LCC_3::Dart_const_handle;

int main(int argc, char* argv[])
{
  std::cout<<"Program facewidth_on_unweighted_map started."<<std::endl;
  std::string filename(argc==1?"data/double-torus.off":argv[1]);
  bool draw=(argc<3?false:std::string(argv[2])=="-draw");

  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  LCC_3 lcc;
  CGAL::load_off(lcc, inp);
  std::cout<<"File '"<<filename<<"' loaded. Finding the facewidth..."<<std::endl;

  CST cst(lcc, true);
  std::vector<Dart_const_handle> cycle=cst.compute_face_width(true);

  if (cycle.size()==0)
  { std::cout<<"  Cannot find such cycle."<<std::endl; }
  else
  {
    std::cout<<"  Number of faces: "<<cycle.size()<<std::endl;

    if (draw) { draw_facewidth(lcc, cycle); }
  }

  return EXIT_SUCCESS;
}
