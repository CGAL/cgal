#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#if defined(CGAL_USE_OPENMESH)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#endif

#include <CGAL/boost/graph/copy_face_graph.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/unordered_map.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Polyhedron_3<Kernel>                       Source;
typedef boost::graph_traits<Source>::vertex_descriptor   sm_vertex_descriptor;
typedef boost::graph_traits<Source>::halfedge_descriptor sm_halfedge_descriptor;
typedef boost::graph_traits<Source>::face_descriptor     sm_face_descriptor;

typedef CGAL::Exact_predicates_exact_constructions_kernel Other_kernel;
typedef Other_kernel::Point_3                             Point;

int main(int argc, char* argv[])
{
  Source S;

  std::ifstream in((argc>1)?argv[1]:"cube.off");
  in >> S;

  // Note that the vertex_point property of the Source and Target1
  // come from different kernels.
  typedef CGAL::Surface_mesh<Point> Target1;
  Target1 T1;
  {
    CGAL::copy_face_graph(S, T1);
    std::ofstream out("sm.off");
    out << T1;
  }

#if defined(CGAL_USE_OPENMESH)
  typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Target2;
  Target2 T2;
  {
    typedef boost::graph_traits<Target2>::vertex_descriptor   tm_vertex_descriptor;
    typedef boost::graph_traits<Target2>::halfedge_descriptor tm_halfedge_descriptor;
    typedef boost::graph_traits<Target2>::face_descriptor     tm_face_descriptor;

    // Use an unordered_map to keep track of elements.
    boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor>     v2v;
    boost::unordered_map<sm_halfedge_descriptor, tm_halfedge_descriptor> h2h;
    boost::unordered_map<sm_face_descriptor, tm_face_descriptor>         f2f;
    
    CGAL::copy_face_graph(S, T2, std::inserter(v2v, v2v.end()),
                                 std::inserter(h2h, h2h.end()),
                                 std::inserter(f2f, f2f.end()));
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
    
    CGAL::copy_face_graph(T1, S, std::inserter(v2v, v2v.end()), std::inserter(h2h, h2h.end()));
    std::ofstream out("reverse.off");
    out << S;
  }
  return 0;
}
