#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>

#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>

#include <CGAL/internal/boost/function_property_map.hpp>

#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               FT;
typedef K::Weighted_point_2                                 Weighted_point;

typedef CGAL::Regular_triangulation_2<K> Triangulation;

// As we only consider finite vertices and edges
// we need the following filter

template <typename T>
struct Is_finite {

  const T* t_;

  Is_finite() : t_(NULL) { }
  Is_finite(const T& t) : t_(&t) { }

  template <typename VertexOrEdge>
  bool operator()(const VertexOrEdge& voe) const {
    return ! t_->is_infinite(voe);
  }
};

template <typename T>
struct Compute_edge_weight
{
  const T& t_;

  Compute_edge_weight(const T& t) : t_(t) { }

  typedef typename boost::graph_traits<T>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<T>::edge_descriptor edge_descriptor;

  FT operator()(const edge_descriptor ed) const {
    vertex_descriptor svd = source(ed, t_);
    vertex_descriptor tvd = target(ed, t_);
    typename T::Vertex_handle sv = svd;
    typename T::Vertex_handle tv = tvd;
    return CGAL::power_product(sv->point(), tv->point());
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
main(int argc,char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/weighted_points.xyw";
  std::ifstream input(filename);
  Triangulation t;
  Filter is_finite(t);
  Finite_triangulation ft(t, is_finite, is_finite);

  Weighted_point wp ;
  while(input >> wp){
    t.insert(wp);
  }

  // Note that with the input "data/weighted_points.xyw", there is one hidden vertex
  std::cout << "number of hidden vertices: " << t.number_of_hidden_vertices() << std::endl;

  vertex_iterator vit, ve;
  // Associate indices to the vertices
  int index = 0;
  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vit,ve)=boost::vertices(ft); vit!=ve; ++vit ){
    vertex_descriptor  vd = *vit;
    vertex_id_map[vd]= index++;
  }

  // We use a custom edge length property map that computes the power distance
  // between the extremities of an edge of the regular triangulation
  typedef Compute_edge_weight<Triangulation> Edge_weight_functor;

  // In the function call you can see a named parameter: vertex_index_map
   std::list<edge_descriptor> mst;
   boost::kruskal_minimum_spanning_tree(
            ft, std::back_inserter(mst),
            vertex_index_map(vertex_index_pmap).
            weight_map(CGAL::internal::boost_::make_function_property_map<
              edge_descriptor, FT, Edge_weight_functor>(Edge_weight_functor(t))));

   std::cout << "The edges of the Euclidean mimimum spanning tree:" << std::endl;
   for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it){
     edge_descriptor ed = *it;
     vertex_descriptor svd = source(ed,t);
     vertex_descriptor tvd = target(ed,t);
     Triangulation::Vertex_handle sv = svd;
     Triangulation::Vertex_handle tv = tvd;
     std::cout << "[ " << sv->point() << "  |  " << tv->point() << " ] " << std::endl;
   }

   return 0;
}
