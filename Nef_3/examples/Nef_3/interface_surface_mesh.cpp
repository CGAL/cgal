#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Timer.h>

#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3                                        Point;
typedef CGAL::Surface_mesh<Point>                         Mesh;
typedef CGAL::Nef_polyhedron_3<K>                         Nef_polyhedron;


int main(int argc, char* argv[])
{
  const std::string filenameP = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/knot.off");
  const std::string filenameQ = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/knot2.off");
  const std::string filenameR = (argc > 3) ? argv[3] : "union.off";

  Mesh meshP, meshQ, result;
  CGAL::IO::read_polygon_mesh(filenameP, meshP);
  CGAL::IO::read_polygon_mesh(filenameQ, meshQ);

  CGAL::Timer timer;
  timer.start();
  Nef_polyhedron NP(meshP, get(CGAL::halfedge_index, meshP), get(CGAL::face_index, meshP));
  Nef_polyhedron NQ(meshQ, get(CGAL::halfedge_index, meshQ), get(CGAL::face_index, meshQ));


  NP += NQ;
  timer.stop();
  std::cout << timer.time() << " sec" << std::endl;
  if(NP.is_simple()) {
    CGAL::convert_nef_polyhedron_to_polygon_mesh(NP, result);
    CGAL::IO::write_polygon_mesh(filenameR, result);
  }
  else
    std::cerr << "result is not a 2-manifold." << std::endl;

    return 0;
}

