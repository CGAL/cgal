#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>

#include <unordered_map>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> Mesh;


typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

int main()
{
  {
    std::unordered_map<vertex_descriptor, int> vmap;
    std::unordered_map<halfedge_descriptor, int> hmap;
    std::unordered_map<edge_descriptor, int> emap;
    std::unordered_map<face_descriptor, int> fmap;
  }
  {
    boost::unordered_map<vertex_descriptor, int> vmap;
    boost::unordered_map<halfedge_descriptor, int> hmap;
    boost::unordered_map<edge_descriptor, int> emap;
    boost::unordered_map<face_descriptor, int> fmap;
  }
  return 0;
}
