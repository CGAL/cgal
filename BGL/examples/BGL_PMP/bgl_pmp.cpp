
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <pmp/SurfaceMesh.h>

#include <CGAL/boost/graph/graph_traits_SurfaceMesh.h>
#include <CGAL/boost/graph/properties_SurfaceMesh.h>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

typedef boost::graph_traits<pmp::SurfaceMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::face_descriptor face_descriptor;

int main(int argc, char* argv[])
{
  pmp::SurfaceMesh sm;
  vertex_descriptor v0, v1, v2, v3;
  v0 = sm.add_vertex(pmp::Point(0,0,0));
  v1 = sm.add_vertex(pmp::Point(1,0,0));
  v2 = sm.add_vertex(pmp::Point(0,2,0));

  sm.add_triangle(v0,v1,v2);

  std::list<edge_descriptor> mst;

  boost::kruskal_minimum_spanning_tree(sm,
                                       std::back_inserter(mst));

 for(edge_descriptor e : mst){
    vertex_descriptor s = source(e,sm);
    vertex_descriptor t = target(e,sm);
    std::cout << sm.position(s) << " -- " << sm.position(t) << std::endl;
 }
  return 0;
}
