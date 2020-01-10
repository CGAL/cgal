
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <pmp/SurfaceMesh.h>

namespace std {

template<>
struct iterator_traits<pmp::SurfaceMesh::VertexIterator> {
  typedef int difference_type;
  typedef pmp::Vertex value_type;
  typedef const pmp::Vertex& reference;
  typedef pmp::Vertex* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
  
template<>
struct iterator_traits<pmp::SurfaceMesh::EdgeIterator> {
  typedef int difference_type;
  typedef pmp::Edge value_type;
  typedef const pmp::Edge& reference;
  typedef pmp::Edge* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
template<>
struct iterator_traits<pmp::SurfaceMesh::HalfedgeIterator> {
  typedef int difference_type;
  typedef pmp::Halfedge value_type;
  typedef const pmp::Halfedge& reference;
  typedef pmp::Halfedge* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
template<>
struct iterator_traits<pmp::SurfaceMesh::FaceIterator> {
  typedef int difference_type;
  typedef pmp::Face value_type;
  typedef const pmp::Face& reference;
  typedef pmp::Face* pointer;
  typedef std::random_access_iterator_tag iterator_category;
};
}

#include <CGAL/boost/graph/graph_traits_SurfaceMesh.h>
#include <CGAL/boost/graph/properties_SurfaceMesh.h>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

typedef boost::graph_traits<pmp::SurfaceMesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<pmp::SurfaceMesh>::halfedge_descriptor halfedge_descriptor;

int main()
{
  typedef boost::property_map<pmp::SurfaceMesh, boost::vertex_index_t>::type VIM;
  typedef boost::property_map<pmp::SurfaceMesh, boost::edge_index_t>::type EIM;
  typedef boost::property_map<pmp::SurfaceMesh, CGAL::vertex_point_t>::type VPM;
  pmp::SurfaceMesh sm;
  pmp::Vertex v0, v1, v2, v3;
  v0 = sm.add_vertex(pmp::Point(0,0,0));
  v1 = sm.add_vertex(pmp::Point(1,0,0));
  v2 = sm.add_vertex(pmp::Point(0,2,0));

  sm.add_triangle(v0,v1,v2);

  std::list<edge_descriptor> mst;

  boost::kruskal_minimum_spanning_tree(sm, 
                                       std::back_inserter(mst));

  edge_descriptor e = *(edges(sm).begin());
  std::cout << e << std::endl;
  halfedge_descriptor h = halfedge(e,sm);
  std::cout << h << std::endl;
  edge_descriptor e2 = edge(h,sm);
  std::cout << e2 << std::endl;
  halfedge_descriptor h2 = halfedge(e2,sm);
  std::cout << h2 << std::endl;

  h = opposite(h,sm);
  e = edge(h,sm);
  std::cout << e << std::endl;
  h = halfedge(e,sm);
  std::cout << h << std::endl;
  
  
  VIM vim = get(boost::vertex_index, sm);
  VPM vpm = get(CGAL::vertex_point, sm);
  for(auto v : vertices(sm)){
    std::cout << degree(v,sm) << " " << get(vim,v) << " " << get(vpm,v) << std::endl;
  }

  EIM eim = get(boost::edge_index,sm);
  
  for(auto e : edges(sm)){
    std::cout << e  << " " <<  get(eim,e) << std::endl;
  }

  typedef boost::property_map<pmp::SurfaceMesh, CGAL::dynamic_vertex_property_t<int>>::type V_index_map;
  V_index_map dvim;
  dvim = get(CGAL::dynamic_vertex_property_t<int>(), sm);

  return 0;
}
