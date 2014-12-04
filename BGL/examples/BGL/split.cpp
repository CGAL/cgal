
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;

typedef boost::adjacency_list < boost::listS,
                                boost::vecS, 
                                boost::undirectedS,
                                Point_2 > G;

typedef boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

typedef std::vector<Point_2> Polyline_2;

struct Is_terminal
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor vd , const Graph& g )
  {
    return false; // degree(vd,g) != 2; is a bad test in case of parallel edges
  }
};


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
};


int main()
{
  G g;

  std::list<Polyline_2> polylines;  
  Polyline_visitor<G> polyline_visitor(polylines, g);
  std::map<Point_2, vertex_descriptor> p2vd;

  int n;
  std::cin >> n; // number of segments


  Point_2 p, q;
  vertex_descriptor  vdp, vdq; 
  for(int i=0; i < n; i++){
    std::cin >> p >> q;
   
    if(p2vd.find(p) == p2vd.end()){
      vdp = add_vertex(g);
      g[vdp] = p;
      p2vd[p] = vdp;
    } else {
      vdp = p2vd[p];
    }
    if(p2vd.find(q) == p2vd.end()){
      vdq = add_vertex(g);
      g[vdq] = q;
      p2vd[q] = vdq;
    } else {
      vdq = p2vd[q];
    }
    boost::add_edge(vdp, vdq, g);
  }

   CGAL::split_graph_into_polylines( g,
                                     polyline_visitor,
                                     Is_terminal() );
   std::cout.precision(17);

   
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
