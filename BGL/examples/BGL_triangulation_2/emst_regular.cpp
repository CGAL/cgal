#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>

#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <fstream>
#include <iostream>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::FT                                                   FT;
typedef K::Weighted_point_2                                     Weighted_point;

typedef CGAL::Regular_triangulation_2<K>                        Triangulation;

typedef boost::graph_traits<Triangulation>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator     vertex_iterator;
typedef boost::graph_traits<Triangulation>::edge_descriptor     edge_descriptor;

template <typename T>
struct Compute_edge_weight
{
  const T& t_;

  Compute_edge_weight(const T& t) : t_(t) { }

  typedef typename boost::graph_traits<T>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<T>::edge_descriptor   edge_descriptor;

  FT operator()(const edge_descriptor ed) const {
    vertex_descriptor svd = source(ed, t_);
    vertex_descriptor tvd = target(ed, t_);
    typename T::Vertex_handle sv = svd;
    typename T::Vertex_handle tv = tvd;
    return CGAL::power_product(sv->point(), tv->point());
  }
};

// The BGL makes use of indices associated to the vertices
// We use a std::map to store the index
typedef std::map<vertex_descriptor,int> VertexIndexMap;

// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap> VertexIdPropertyMap;

int main(int argc,char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/weighted_points.xyw";
  std::ifstream input(filename);
  Triangulation tr;

  Weighted_point wp;
  while(input >> wp)
    tr.insert(wp);

  // Note that with the input "data/weighted_points.xyw", there is one hidden vertex
  std::cout << "number of hidden vertices: " << tr.number_of_hidden_vertices() << std::endl;

  // Associate indices to the vertices
  VertexIndexMap vertex_id_map;
  VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
  int index = 0;

  for(vertex_descriptor vd : vertices(tr))
    vertex_id_map[vd]= index++;

  // We use a custom edge length property map that computes the power distance
  // between the extremities of an edge of the regular triangulation
  typedef Compute_edge_weight<Triangulation> Edge_weight_functor;

  // In the function call you can see a named parameter: vertex_index_map
  std::list<edge_descriptor> mst;
  boost::kruskal_minimum_spanning_tree(tr, std::back_inserter(mst),
                                       vertex_index_map(vertex_index_pmap)
                                       .weight_map(boost::make_function_property_map<
                                          edge_descriptor, FT, Edge_weight_functor>(Edge_weight_functor(tr))));

  std::cout << "The edges of the Euclidean mimimum spanning tree:" << std::endl;
  for(edge_descriptor ed : mst)
  {
    vertex_descriptor svd = source(ed, tr);
    vertex_descriptor tvd = target(ed, tr);
    Triangulation::Vertex_handle sv = svd;
    Triangulation::Vertex_handle tv = tvd;
    std::cout << "[ " << sv->point() << "  |  " << tv->point() << " ] " << std::endl;
  }

  return EXIT_SUCCESS;
}
