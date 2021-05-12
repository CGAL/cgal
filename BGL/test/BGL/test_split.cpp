
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/graph/adjacency_list.hpp>

#include <map>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;


typedef boost::adjacency_list<boost::vecS, boost::setS, boost::undirectedS, Point_2 > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;

typedef std::map<Point_2, vertex_descriptor> Point_vertex_map;

typedef std::vector<Point_2> Polyline_2;


// inserts a polyline into a graph
void insert(const std::vector<Point_2>& poly, Graph& graph, Point_vertex_map& pvmap)
{
    vertex_descriptor u, v;
    for (std::size_t i = 0; i < poly.size(); i++) {
        // check if the point is not yet in the graph
        if (pvmap.find(poly[i]) == pvmap.end()) {
            v = add_vertex(graph);
            pvmap[poly[i]] = v;
        }
        else {
            v = pvmap[poly[i]];
        }
        graph[v] = poly[i];  // associate the point to the vertex
        if (i != 0) {
            add_edge(u, v, graph);
        }
        u = v;
    }
}

template <typename Graph>
struct Polyline_visitor
{
  std::list<Polyline_2>& polylines;
  const Graph& points_pmap;

  Polyline_visitor(std::list<Polyline_2>& lines,
                   const Graph& points_property_map)
    : polylines(lines),
      points_pmap(points_property_map)
  {}

  void start_new_polyline()
  {
    Polyline_2 V;
    polylines.push_back(V);
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    Polyline_2& polyline = polylines.back();
    polyline.push_back(points_pmap[vd]);
  }

  void end_polyline()
  {}

};


int main()
{
  Polyline_2 polyA = { Point_2(0,0), Point_2(1,0), Point_2(2,0), Point_2(3,0), Point_2(4,0)};
  Polyline_2 polyB = { Point_2(1,-1), Point_2(1,0), Point_2(2,0), Point_2(2,1), Point_2(2,2) };

  Graph graph;
  Point_vertex_map pvmap;

  insert(polyA, graph, pvmap);
  insert(polyB, graph, pvmap);

  std::list<Polyline_2> polylines;
  Polyline_visitor<Graph> polyline_visitor(polylines, graph);

  CGAL::split_graph_into_polylines( graph,
                                    polyline_visitor);


  for(std::list<Polyline_2>::iterator it = polylines.begin(); it!= polylines.end(); ++it){
     Polyline_2& poly = *it;
     std::size_t n;
     if(poly.front() == poly.back()){
       std::cout << "POLYGON" << std::endl;
       n = poly.size() -1;
     }else{
       std::cout << "POLYLINE" << std::endl;
       n = poly.size();
     }
     for(std::size_t j=0; j < n; j++){
       std::cout << poly[j] << std::endl;
     }
     std::cout << std::endl;
   }

  return 0;

}
