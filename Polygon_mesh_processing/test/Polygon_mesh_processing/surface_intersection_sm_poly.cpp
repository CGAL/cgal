#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/iterator.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Timer.h>

namespace PMP=CGAL::Polygon_mesh_processing;

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

template <class TriangleMesh>
void run(const char* filename1, const char* filename2, const char* msg)
{
  TriangleMesh mesh1;
  std::ifstream input(filename1);
  if ( !input || !(input >> mesh1) ) {
    std::cerr << filename1 << " is not a valid off file.\n";
    exit(1);
  }
  input.close();

  input.open(filename2);
  TriangleMesh mesh2;
  if ( !input || !(input >> mesh2) ) {
    std::cerr << filename2 << " is not a valid off file.\n";
    exit(1);
  }

  CGAL::Timer time;
  time.start();
  PMP::surface_intersection(mesh1, mesh2, CGAL::Emptyset_iterator());
  time.stop();

  std::cout << "Runtime for " << msg << " " << time.time() << "\n";
}

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/XXXXXXX.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/XXXXXXX.off";

  run<Mesh>(filename1,filename2,"Surface_mesh");
  run<Polyhedron>(filename1,filename2,"Polyhedron_3");

  return 0;
}
