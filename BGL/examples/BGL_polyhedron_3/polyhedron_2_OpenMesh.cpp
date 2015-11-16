
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <OpenMesh/Core/IO/MeshIO.hh>

#if 1
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#else
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#endif

#include <CGAL/boost/graph/convert_surface_mesh.h>

#include <boost/unordered_map.hpp>

#include <iostream>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Source;

#if 1
typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Target;
#else
typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> Target;
#endif
typedef boost::graph_traits<Source>::vertex_descriptor sm_vertex_descriptor;
typedef boost::graph_traits<Target>::vertex_descriptor tm_vertex_descriptor;

typedef boost::graph_traits<Source>::halfedge_descriptor sm_halfedge_descriptor;
typedef boost::graph_traits<Target>::halfedge_descriptor tm_halfedge_descriptor;

namespace OpenMesh {

inline  std::size_t hash_value(const VertexHandle&  i)
  {
    return i.idx();
  }

inline  std::size_t hash_value(const HalfedgeHandle&  i)
  {
    return i.idx();
  }

inline  std::size_t hash_value(const FaceHandle&  i)
  {
    return i.idx();
  }

}

int main(int argc, char* argv[])
{
  Source S;
  Target T;
  std::ifstream in((argc>1)?argv[1]:"cube.off");
  in >> S;

  {
    boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor> v2v;
    boost::unordered_map<sm_halfedge_descriptor, tm_halfedge_descriptor> h2h;
    
    convert_surface_mesh(S,T,v2v,h2h);
    OpenMesh::IO::write_mesh(T, "om.off");
  }
  S.clear();
  {
    boost::unordered_map<tm_vertex_descriptor, sm_vertex_descriptor> v2v;
    boost::unordered_map<tm_halfedge_descriptor, sm_halfedge_descriptor> h2h;
    
    convert_surface_mesh(T,S,v2v,h2h);
    std::ofstream out("reverse.off");
    out << S << std::endl;
  }

  
  return 0;
}
