#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main(int argc, char** argv )
{
  Mesh mesh;

  const std::string filename = (argc>1)?argv[1]:CGAL::data_file_path("meshes/in.off");
  const char* outname= (argc>2)?argv[2]:"out.om";
  CGAL::IO::read_polygon_mesh(filename, mesh);

  mesh.request_vertex_status();

  int i = 0;
  for(auto v : vertices(mesh)){
        mesh.status(v).set_selected((i%2) == 0);
        ++i;
    }


  OpenMesh::IO::write_mesh(mesh, outname, OpenMesh::IO::Options::Status);

  Mesh mesh2;
  OpenMesh::IO::Options options  = OpenMesh::IO::Options::Status;
  bool read = OpenMesh::IO::read_mesh(mesh2, outname, options);
  std::cout << num_vertices(mesh2) << std::endl;
  assert(read);
  if(options.vertex_has_status()){
    for(auto v : vertices(mesh2)){
        std::cout << std::boolalpha <<  mesh2.status(v).selected() << std::endl;
    }
  }else{
    std::cout << "no vertex status" << std::endl;
  }

  return 0;
}
