#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main(int argc, char** argv )
{
  Mesh mesh;

  std::vector<vertex_descriptor> V;
  const char* filename = (argc>1)?argv[1]:"in.off";
  const char* outname= (argc>2)?argv[2]:"out.off";
  CGAL::IO::read_polygon_mesh(filename, mesh);

  for(vertex_descriptor vd : vertices(mesh)){
    for(halfedge_descriptor hd : CGAL::halfedges_around_target(vd,mesh)){
      if(! CGAL::is_border(edge(hd,mesh),mesh)){
        CGAL::Euler::flip_edge(hd,mesh);
        CGAL::IO::write_polygon_mesh(outname, mesh);
        return 0;
      }
    }
  }
  return 0;
}
