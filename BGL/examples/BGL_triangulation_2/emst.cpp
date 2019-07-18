#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <fstream>
#include <iostream>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef K::Point_2                                                  Point;

typedef CGAL::Delaunay_triangulation_2<K>                           Triangulation;

typedef boost::graph_traits<Triangulation>::vertex_descriptor       vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator         vertex_iterator;
typedef boost::graph_traits<Triangulation>::edge_descriptor         edge_descriptor;

// The BGL makes use of indices associated to the vertices
// We use a std::map to store the index
typedef std::map<vertex_descriptor,int>                             VertexIndexMap;

// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap>             VertexIdPropertyMap;

int main(int argc,char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/points.xy";
  std::ifstream input(filename);
  Triangulation tr;

  Point p;
  while(input >> p)
    tr.insert(p);

  // Associate indices to the vertices
  VertexIndexMap vertex_id_map;
  VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
  int index = 0;

  for(vertex_descriptor vd : vertices(tr))
    vertex_id_map[vd] = index++;

  // We use the default edge weight which is the squared length of the edge
  // This property map is defined in graph_traits_Triangulation_2.h

  // In the function call you can see a named parameter: vertex_index_map
  std::list<edge_descriptor> mst;
  boost::kruskal_minimum_spanning_tree(tr, std::back_inserter(mst),
                                       vertex_index_map(vertex_index_pmap));

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
