#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel           K;
typedef K::Point_2                                                    Point;

typedef CGAL::Triangulation_2<K>                                      Triangulation;

typedef boost::graph_traits<Triangulation>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator           vertex_iterator;

typedef std::map<vertex_descriptor, int>                              VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;

int main(int argc,char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/points.xy";
  std::ifstream input(filename);
  Triangulation tr;

  Point p ;
  while(input >> p)
    tr.insert(p);

  vertex_iterator vit, ve;

  // Associate indices to the vertices
  VertexIndexMap vertex_id_map;
  VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
  int index = 0;

  for(vertex_descriptor vd : vertices(tr))
    vertex_id_map[vd] = index++;

  // Dijkstra's shortest path needs property maps for the predecessor and distance
  // We first declare a vector
  std::vector<vertex_descriptor> predecessor(num_vertices(tr));

  // and then turn it into a property map
  boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap>
      predecessor_pmap(predecessor.begin(), vertex_index_pmap);

  std::vector<double> distance(num_vertices(tr));
  boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap>
    distance_pmap(distance.begin(), vertex_index_pmap);

  // start at an arbitrary vertex
  vertex_descriptor source = *vertices(tr).first;
  std::cout << "\nStart dijkstra_shortest_paths at " << source->point() <<"\n";

  boost::dijkstra_shortest_paths(tr, source,
                                 distance_map(distance_pmap)
                                 .predecessor_map(predecessor_pmap)
                                 .vertex_index_map(vertex_index_pmap));

  for(vertex_descriptor vd : vertices(tr))
  {
    std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "] ";
    std::cout << " has distance = " << boost::get(distance_pmap,vd)
              << " and predecessor ";
    vd = boost::get(predecessor_pmap,vd);
    std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "]\n ";
  }

  return EXIT_SUCCESS;
}
