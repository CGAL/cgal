#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

typedef OpenMesh::PolyMesh_ArrayKernelT< > Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;



int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";

  Mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);

  CGAL::OM_pmap<Mesh, face_descriptor, Vector> fnormals(mesh);
  CGAL::OM_pmap<Mesh, vertex_descriptor, Vector> vnormals(mesh);
  mesh.add_property(fnormals.handle());
  mesh.add_property(vnormals.handle());

  Vector v(0, 0, 0);
  for(vertex_descriptor vd : vertices(mesh)){
      put(vnormals, vd, v);
    }
  for(face_descriptor fd : faces(mesh)){
      put(fnormals, fd, v);
    }
  CGAL::Polygon_mesh_processing::compute_normals
    (mesh,
     vnormals,
     fnormals,
     CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
     geom_traits(K()));

  std::cout << "Face normals :" << std::endl;
  for(face_descriptor fd : faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
  std::cout << "Normals at vertices :" << std::endl;
  for(vertex_descriptor vd : vertices(mesh)){
    std::cout << vnormals[vd] << std::endl;
  }
  return 0;
}
