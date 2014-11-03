
#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <map>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;


typedef boost::adjacency_list < boost::listS,
                                boost::vecS, 
                                boost::undirectedS,
                                Point_2 > G;

typedef boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

typedef std::vector<Point_2> Polyline;

struct Is_terminal
{
  template <typename VertexDescriptor, typename Graph>
  bool operator ()(VertexDescriptor vd , const Graph& g )
  {
    return degree(vd,g) != 2;
  }
};


template <typename Graph> 
struct Polyline_visitor
{
  std::list<Polyline>& polylines;
  Graph& points_pmap;

  Polyline_visitor(std::list<Polyline>& lines,
                   Graph& points_property_map)
    : polylines(lines),
      points_pmap(points_property_map)
  {}

  void start_new_polyline()
  {
    Polyline V;
    polylines.push_back(V);
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    Polyline& polyline = polylines.back();
    polyline.push_back(points_pmap[vd]);
  }
};


struct LessSegment {
  bool operator()(const Segment_2& a, const Segment_2& b) const
  {
    if(a.source() < b.source()){
      return true;
    }
    if(a.source() > b.source()){
      return false;
    }
    if(a.target() < b.target()){
      return true;
    }
    if(a.target() > b.target()){
      return false;
    }
    return false;
  }
};

int main()
{
  G g;
  std::list<Polyline> polylines;  
  Polyline_visitor<G> polyline_visitor(polylines, g);
  std::map<Point_2, vertex_descriptor> p2vd;

  int n;
  std::cin >> n; // number of segments

  std::set<Segment_2,LessSegment> segments;

  Point_2 p, q;
  vertex_descriptor  vdp, vdq; 
  for(int i=0; i < n; i++){
    std::cin >> p >> q;
    if(p == q){
      //std::cerr << "ignore degenerate segment"<< std::endl;
      continue;
    }
    Segment_2 s = (p < q)? Segment_2(p,q) : Segment_2(q,p);
    if(segments.find(s) != segments.end()){
      //std::cerr << "ignore duplicate segment"<< std::endl;
      continue;
    } else {
      segments.insert(s);
    }

    if(p2vd.find(p) == p2vd.end()){
      vdp =  add_vertex(g);
      g[vdp] = p;
      p2vd[p] = vdp;
    } else {
      vdp = p2vd[p];
    }
    if(p2vd.find(q) == p2vd.end()){
      vdq =  add_vertex(g);
      g[vdq] = q;
      p2vd[q] = vdq;
    } else {
      vdq = p2vd[q];
    }
    boost::add_edge(vdp, vdq, g);
  }

  { 
    typedef boost::graph_traits<G>::vertex_iterator vertex_iterator; 
    vertex_iterator b,e;
    boost::tie(b,e) = vertices(g);
    for(; b!= e; ++b){
      std::cerr << degree(*b,g) << std::endl;
    }
  }
  std::cerr <<  "A  " << std::endl;
   CGAL::split_graph_into_polylines( g,
                                      polyline_visitor,
                                     CGAL::IsTerminalDefault() );
   std::cout.precision(17);

   
   for(std::list<Polyline>::iterator it = polylines.begin(); it!= polylines.end(); ++it){
     Polyline& poly = *it;
     std::cout << poly.size() ;
     for(int j=0; j < poly.size(); j++){
       std::cout << "  " << poly[j];
     }
     std::cout << std::endl;
   }
   
  return 0;
}
