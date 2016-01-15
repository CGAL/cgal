#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/internal/Polyhedron_plane_clipping_3.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char** argv)
{
  const char* fname = argc>1?argv[1]:"data/elephant.off";
  std::ifstream input(fname);

  Polyhedron P;
  input >> P;
  std::cout << "Polyhedron with " << P.size_of_vertices() << " vertices\n";
// test clipping by a plane
{
  Polyhedron* P_clipped=
    CGAL::corefinement::clip_polyhedron(P, Kernel::Plane_3(1,1,0,0));
  std::ofstream output("clipped1.off");
  output << *P_clipped;
  delete P_clipped;
}
// test clipping with the opposite plane
{
  Polyhedron* P_clipped=
    CGAL::corefinement::clip_polyhedron(P, Kernel::Plane_3(-1,-1,0,0));
  std::ofstream output("clipped2.off");
  output << *P_clipped;
  delete P_clipped;
}
// test clipping with a plane with no intersection having the elephant
// on its positive side, result should be empty
{
  Polyhedron* P_clipped=
    CGAL::corefinement::clip_polyhedron(P, Kernel::Plane_3(1,0,0,1));
  if (!P_clipped->empty()){
    std::cerr << "Error: Polyhedron should be empty!\n";
    return 1;
  }
  delete P_clipped;
}
// test clipping with a plane with no intersection having the elephant
// on its negative side, result should be empty
{
  Polyhedron* P_clipped=
    CGAL::corefinement::clip_polyhedron(P, Kernel::Plane_3(-1,0,0,-1));
  if (P_clipped->size_of_vertices()!=P.size_of_vertices()){
    std::cerr << "Error: Polyhedron should be full!\n";
    return 1;
  }
  delete P_clipped;
}

// clip with a plane but do not close the polyhedron
{
  CGAL::corefinement::inplace_clip_open_polyhedron(P, Kernel::Plane_3(1,1,0,0));
  std::ofstream output("open_clipped.off");
  output << P;
}

  return 0;
}
