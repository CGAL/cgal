
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/boost/graph/convert_surface_mesh.h>

#include <boost/unordered_map.hpp>

#include <iostream>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Source;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Target;

typedef boost::graph_traits<Source>::vertex_descriptor sm_vertex_descriptor;
typedef boost::graph_traits<Target>::vertex_descriptor tm_vertex_descriptor;

typedef boost::graph_traits<Source>::halfedge_descriptor sm_halfedge_descriptor;
typedef boost::graph_traits<Target>::halfedge_descriptor tm_halfedge_descriptor;

boost::unordered_map<sm_vertex_descriptor, tm_vertex_descriptor> v2v;
boost::unordered_map<sm_halfedge_descriptor, tm_halfedge_descriptor> h2h;

int main(int, char* argv[])
{
  Source S;
  Target T;
  std::ifstream in(argv[1]);
  in >> S;

  convert_surface_mesh(S,T,v2v,h2h);

  OpenMesh::IO::write_mesh(T, "om.off");

  return 0;
}
