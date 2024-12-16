#include <CGAL/basic.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

typedef OpenMesh::PolyMesh_ArrayKernelT<> Om;
typedef boost::graph_traits< Om >         Traits;
typedef Traits::edge_descriptor           edge_descriptor;
typedef Traits::halfedge_descriptor       halfedge_descriptor;
typedef Traits::vertex_descriptor         vertex_descriptor;
typedef Traits::face_descriptor           face_descriptor;
//typedef Kernel::Point_3 Point_3;

int main()
{
  Om om;
  for (Om::EdgeHandle ed : edges(om)) {
          std::cout << "edge" << std::endl;
  }
  return 0;
}
