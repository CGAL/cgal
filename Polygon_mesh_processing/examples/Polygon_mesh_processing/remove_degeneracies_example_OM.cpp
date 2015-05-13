#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef OpenMesh::PolyMesh_ArrayKernelT< > Mesh;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/degtri_sliding.off";

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);


  std::size_t nb
    = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh,
                                                             CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                                                             geom_traits(K()));

  std::cout << "There were " << nb << " degenerate faces in this mesh" << std::endl;
  mesh.garbage_collection();
  OpenMesh::IO::write_mesh(mesh, "repaired.off");
  
  return 0;
}
