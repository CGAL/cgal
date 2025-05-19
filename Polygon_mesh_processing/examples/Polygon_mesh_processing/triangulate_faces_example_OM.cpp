#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef OpenMesh::PolyMesh_ArrayKernelT< >                  Mesh;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube_quad.off");
  const char* outfilename = (argc > 2) ? argv[2] : "cube_tri.off";

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);

  CGAL::Polygon_mesh_processing::triangulate_faces(mesh,
     CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                                                   geom_traits(Kernel()));

  mesh.garbage_collection();
  OpenMesh::IO::write_mesh(mesh,outfilename);

  return 0;
}
