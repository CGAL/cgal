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

int main()
{
  Om om;
  for (Om::EdgeHandle ed : edges(om))
  {
    CGAL_USE(ed);
  }
  for (edge_descriptor ed : edges(om))
  {
    CGAL_USE(ed);
  }
  for (halfedge_descriptor hd : edges(om))
  {
    CGAL_USE(hd);
  }
  for (face_descriptor fd : faces(om))
  {
    CGAL_USE(fd);
  }
  for (vertex_descriptor vd : vertices(om))
  {
    CGAL_USE(vd);
  }
  return 0;
}
