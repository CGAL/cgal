#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef OpenMesh::PolyMesh_ArrayKernelT< > Mesh;
int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/quads_to_stitch.off");

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);

  std::cout << "Before stitching : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << num_vertices(mesh) << std::endl;
  std::cout << "\t Number of halfedges :\t" << num_halfedges(mesh) << std::endl;
  std::cout << "\t Number of facets    :\t" << num_faces(mesh) << std::endl;

  CGAL::Polygon_mesh_processing::stitch_borders(mesh, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));

  mesh.garbage_collection();
  std::cout << "Stitching done : " << std::endl;
  std::cout << "\t Number of vertices  :\t" << num_vertices(mesh) << std::endl;
  std::cout << "\t Number of halfedges :\t" << num_halfedges(mesh) << std::endl;
  std::cout << "\t Number of facets    :\t" << num_faces(mesh) << std::endl;

  return 0;
}
