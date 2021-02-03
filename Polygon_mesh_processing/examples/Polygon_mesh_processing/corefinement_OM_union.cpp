#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef OpenMesh::PolyMesh_ArrayKernelT< >                    Mesh;


namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";

  Mesh mesh1, mesh2;

  OpenMesh::IO::read_mesh(mesh1, filename1);
  OpenMesh::IO::read_mesh(mesh2, filename2);

  Mesh out;
  bool valid_union = PMP::corefine_and_compute_union(mesh1,mesh2, out);

  if (valid_union)
  {
    std::cout << "Union was successfully computed\n";
    OpenMesh::IO::write_mesh(out, "union.off");
    return 0;
  }
  std::cout << "Union could not be computed\n";
  return 1;
}
