#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

typedef CGAL::Polyhedron_3<K> G;


struct Is_border {
  const G& g;
  Is_border(const G& g)
    : g(g)
  {}

 template <typename Edge>
  bool operator()(const Edge& e) const {
   return is_border(e,g);
  }
};

typedef boost::filtered_graph<G,Is_border> FG;
typedef boost::graph_traits<FG>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

typedef std::vector<Point_3> Polyline_3;

struct Is_terminal
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor vd , const Graph& g )
  {
    return false;
  }
};


struct Polyline_visitor
{
  std::list<Polyline_3>& polylines;

  Polyline_visitor(std::list<Polyline_3>& lines)
    : polylines(lines)
  {}

  void start_new_polyline()
  {
    Polyline_3 V;
    polylines.push_back(V);
  }

  void add_node(boost::graph_traits<G>::vertex_descriptor vd)
  {
    Polyline_3& polyline = polylines.back();
    polyline.push_back(vd->point());
  }
};


int main()
{
  G g;

  std::cin >> g;

  Is_border ib(g);
  FG fg(g,ib);

  std::list<Polyline_3> polylines;
  Polyline_visitor polyline_visitor(polylines);

   CGAL::split_graph_into_polylines( fg,
                                     polyline_visitor,
                                     Is_terminal() );

   std::cout.precision(17);

   for(std::list<Polyline_3>::iterator it = polylines.begin(); it!= polylines.end(); ++it){
     Polyline_3& poly = *it;
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
