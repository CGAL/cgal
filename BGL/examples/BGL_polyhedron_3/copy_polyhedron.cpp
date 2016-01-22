#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#if defined(CGAL_USE_OPENMESH)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

namespace OpenMesh { // auxiliary functions so OpenMesh Handles can be hashed
inline std::size_t hash_value(const VertexHandle&  i) { return i.idx(); }
inline std::size_t hash_value(const HalfedgeHandle&  i) { return i.idx(); }
inline std::size_t hash_value(const FaceHandle&  i) { return i.idx(); }
}

#endif

#include <CGAL/boost/graph/copy_face_graph.h>

#include <boost/unordered_map.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef CGAL::Polyhedron_3<Kernel>                       Source;
typedef boost::graph_traits<Source>::vertex_descriptor   sm_vertex_descriptor;
typedef boost::graph_traits<Source>::halfedge_descriptor sm_halfedge_descriptor;

typedef CGAL::Surface_mesh<Point>                         Target1;

int main(int argc, char* argv[])
{
  Source S;

  std::ifstream in((argc>1)?argv[1]:"cube.off");
  in >> S;

  Target1 T1;
  {
    typedef boost::graph_traits<Target1>::vertex_descriptor   tm_vertex_descriptor;
    typedef boost::graph_traits<Target1>::halfedge_descriptor tm_halfedge_descriptor;

    boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor> v2v;
    boost::unordered_map<sm_halfedge_descriptor, tm_halfedge_descriptor> h2h;
    
    CGAL::copy_face_graph(S,T1,v2v,h2h);
    std::ofstream out("sm.off");
    out << T1;
  }

#if defined(CGAL_USE_OPENMESH)
  typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Target2;
  Target2 T2;
  {
    typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/>    Target2;

    typedef boost::graph_traits<Target2>::vertex_descriptor   tm_vertex_descriptor;
    typedef boost::graph_traits<Target2>::halfedge_descriptor tm_halfedge_descriptor;

    boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor> v2v;
    boost::unordered_map<sm_halfedge_descriptor, tm_halfedge_descriptor> h2h;
    
    CGAL::copy_face_graph(S,T2,v2v,h2h);

    OpenMesh::IO::write_mesh(T2, "om.off");
  }
#endif
  S.clear();
  {
    typedef boost::graph_traits<Target1>::vertex_descriptor   source_vertex_descriptor;
    typedef boost::graph_traits<Target1>::halfedge_descriptor source_halfedge_descriptor;

    typedef boost::graph_traits<Source>::vertex_descriptor   tm_vertex_descriptor;
    typedef boost::graph_traits<Source>::halfedge_descriptor tm_halfedge_descriptor;

    boost::unordered_map<source_vertex_descriptor, tm_vertex_descriptor> v2v;
    boost::unordered_map<source_halfedge_descriptor, tm_halfedge_descriptor> h2h;
    
    CGAL::copy_face_graph(T1,S,v2v,h2h);
    std::ofstream out("reverse.off");
    out << T1;
  }
  return 0;
}
