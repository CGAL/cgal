#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>

#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef K::Point_2                                                      Point;

typedef CGAL::Triangulation_vertex_base_with_id_2<K>                    Tvb;
typedef CGAL::Triangulation_face_base_2<K>                              Tfb;
typedef CGAL::Triangulation_data_structure_2<Tvb, Tfb>                  Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                          Triangulation;

typedef boost::graph_traits<Triangulation>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Triangulation>::vertex_iterator             vertex_iterator;

typedef boost::property_map<Triangulation, boost::vertex_index_t>::type VertexIdPropertyMap;

int main(int argc,char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/points.xy";
  std::ifstream input(filename);
  Triangulation tr;

  Point p;
  while(input >> p)
    tr.insert(p);

  // associate indices to the vertices
  int index = 0;
  for(vertex_descriptor vd : vertices(tr))
    vd->id()= index++;

  VertexIdPropertyMap vertex_index_pmap = get(boost::vertex_index, tr);

  // Dijkstra's shortest path needs property maps for the predecessor and distance
  std::vector<vertex_descriptor> predecessor(num_vertices(tr));
  boost::iterator_property_map<std::vector<vertex_descriptor>::iterator, VertexIdPropertyMap>
    predecessor_pmap(predecessor.begin(), vertex_index_pmap);

  std::vector<double> distance(num_vertices(tr));
  boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap>
    distance_pmap(distance.begin(), vertex_index_pmap);

  vertex_descriptor source = *vertices(tr).first;
  std::cout << "\nStart dijkstra_shortest_paths at " << source->point() << std::endl;

  boost::dijkstra_shortest_paths(tr, source, distance_map(distance_pmap)
                                            .predecessor_map(predecessor_pmap));

  for(vertex_descriptor vd : vertices(tr))
  {
    std::cout << vd->point() << " [" << vd->id() << "] ";
    std::cout << " has distance = "  << get(distance_pmap,vd) << " and predecessor ";

    vd = get(predecessor_pmap,vd);
    std::cout << vd->point() << " [" << vd->id() << "]\n";
  }

  return EXIT_SUCCESS;
}
