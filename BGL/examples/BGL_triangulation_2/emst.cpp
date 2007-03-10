#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

typedef CGAL::Delaunay_triangulation_2<K> Triangulation;

// As we only consider finite vertices and edges
// we need the following filter

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
typedef boost::graph_traits<Finite_triangulation>::edge_descriptor edge_descriptor;

// The BGL makes use of indices associated to the vertices
// We use a std::map to store the index
typedef std::map<vertex_descriptor,int> VertexIndexMap;
VertexIndexMap vertex_id_map;

// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<VertexIndexMap> VertexIdPropertyMap;
VertexIdPropertyMap vertex_index_pmap(vertex_id_map);

int
main(int,char*[])
{
  Triangulation t;
  Filter is_finite(t);
  Finite_triangulation ft(t, is_finite, is_finite);

  Point p ;
  while(std::cin >> p){
    t.insert(p);
  }

  vertex_iterator vit, ve;
  // Associate indices to the vertices
  int index = 0;
  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
    vertex_descriptor  vd = *vit;
    vertex_id_map[vd]= index++;
    }


  // We use the default edge weight which is the squared length of the edge
  // This property map is defined in graph_traits_Triangulation_2.h

  // In the function call you can see a named parameter: vertex_index_map
   std::list<edge_descriptor> mst;
   boost::kruskal_minimum_spanning_tree(t,
					std::back_inserter(mst),
					vertex_index_map(vertex_index_pmap));


   std::cout << "The edges of the Euclidean mimimum spanning tree:" << std::endl;

   for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it){
     edge_descriptor ed = *it;
     vertex_descriptor svd = boost::source(ed,t);
     vertex_descriptor tvd = boost::target(ed,t);
     Triangulation::Vertex_handle sv = svd;
     Triangulation::Vertex_handle tv = tvd;
     std::cout << "[ " << sv->point() << "  |  " << tv->point() << " ] " << std::endl;
   }

   return 0;
}
