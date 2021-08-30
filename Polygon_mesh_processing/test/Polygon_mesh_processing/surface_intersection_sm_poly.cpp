#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/iterator.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, MyTraits>::type LCC;

template <class TriangleMesh>
void run(const char* filename1, const char* filename2, const char* msg)
{
  TriangleMesh mesh1;
  if ( !PMP::IO::read_polygon_mesh(filename1, mesh1) ) {
    std::cerr << filename1 << " is not a valid off file.\n";
    exit(1);
  }

  TriangleMesh mesh2;
  if ( !PMP::IO::read_polygon_mesh(filename2, mesh2) ) {
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
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";

  run<Mesh>(filename1,filename2,"Surface_mesh");
  run<Polyhedron>(filename1,filename2,"Polyhedron_3");
  run<LCC>(filename1,filename2,"Linear_cell_complex");

  return 0;
}
