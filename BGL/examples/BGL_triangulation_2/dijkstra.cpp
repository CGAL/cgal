#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/filtered_graph.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

typedef CGAL::Triangulation_2<K> Triangulation;

// As we want to run Dijskra's shortest path algorithm we only
// consider finite vertices and edges.

template <typename T>
struct Is_finite {

  const T* t_;

  Is_finite()
    : t_(NULL)
  {}

  Is_finite(const T& t)
    : t_(&t)
  { }

  template <typename VertexOrEdge>
  bool operator()(const VertexOrEdge& voe) const {
    return ! t_->is_infinite(voe);
  }
};

typedef Is_finite<Triangulation> Filter;
typedef boost::filtered_graph<Triangulation,Filter,Filter> Finite_triangulation;
typedef boost::graph_traits<Finite_triangulation>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Finite_triangulation>::vertex_iterator vertex_iterator;

typedef std::map<vertex_descriptor,int> VertexIndexMap;
VertexIndexMap vertex_id_map;

typedef boost::associative_property_map<VertexIndexMap> VertexIdPropertyMap;
VertexIdPropertyMap vertex_index_pmap(vertex_id_map);

int
main(int,char*[])
{
  Triangulation t;
  Filter is_finite(t);
  Finite_triangulation ft(t, is_finite, is_finite);

  t.insert(Point(0,0));
  t.insert(Point(1,0));
  t.insert(Point(0.2,0.2));
  t.insert(Point(0,1));
  t.insert(Point(0,2));

  vertex_iterator vit, ve;
  // Associate indices to the vertices
  int index = 0;
  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
    vertex_descriptor  vd = *vit;
    vertex_id_map[vd]= index++;
    }

  // Dijkstra's shortest path needs property maps for the predecessor and distance
  // We first declare a vector
  std::vector<vertex_descriptor> predecessor(boost::num_vertices(ft));
  // and then turn it into a property map
  boost::iterator_property_map<std::vector<vertex_descriptor>::iterator,
                               VertexIdPropertyMap>
    predecessor_pmap(predecessor.begin(), vertex_index_pmap);

  std::vector<double> distance(boost::num_vertices(ft));
  boost::iterator_property_map<std::vector<double>::iterator,
                               VertexIdPropertyMap>
    distance_pmap(distance.begin(), vertex_index_pmap);

  // start at an arbitrary vertex
  vertex_descriptor source = *boost::vertices(ft).first;
  std::cout << "\nStart dijkstra_shortest_paths at " << source->point() <<"\n";

  boost::dijkstra_shortest_paths(ft, source,
				 distance_map(distance_pmap)
				 .predecessor_map(predecessor_pmap)
				 .vertex_index_map(vertex_index_pmap));

  for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
    vertex_descriptor vd = *vit;
    std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "] ";
    std::cout << " has distance = "  << boost::get(distance_pmap,vd)
	      << " and predecessor ";
    vd =  boost::get(predecessor_pmap,vd);
    std::cout << vd->point() << " [" <<  vertex_id_map[vd] << "]\n ";
  }

  return 0;
}
