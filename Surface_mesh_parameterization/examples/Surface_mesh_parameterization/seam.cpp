#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;


template <typename PolygonMesh>
vertex_descriptor closest_vertex(const Point& p, const PolygonMesh& mesh)
{  
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type PPmap;
  PPmap ppmap = get(CGAL::vertex_point, mesh);
  vertex_descriptor res;
  double dist = (std::numeric_limits<double>::max)();
  //std::cerr <<dist << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    //std::cerr << "vd = " << vd << std::endl;
    Point vp = get(ppmap,vd);
    //std::cerr << "vp = " << vp << std::endl;
    double d = CGAL::squared_distance(p,vp);
    //std::cerr << "d = " << d << std::endl;
    if(d < dist){
      res = vd;
      dist = d;
      //std::cerr <<res << std::endl;
      //std::cerr <<dist << std::endl;
    } 
  }
  return res;
}


int main(int, char* argv[]) 
{
  Mesh mesh;
  std::ifstream in(argv[1]);
  in >> mesh;

  typedef boost::property_map<Mesh, CGAL::vertex_point_t>::const_type PPmap;
  PPmap ppmap = get(CGAL::vertex_point, mesh);

  Point sp, tp;
  std::ifstream pin(argv[2]);
  pin >> sp >>  tp;

  vertex_descriptor sv, tv;
  sv = closest_vertex(sp, mesh);
  tv = closest_vertex(tp, mesh);
  
  std::cout << "2\n" << get(ppmap,sv)  << std::endl;
  std::cout <<get(ppmap,tv)  << std::endl;

  Mesh::Property_map<vertex_descriptor,vertex_descriptor> predecessor;
  predecessor = mesh.add_property_map<vertex_descriptor,vertex_descriptor>("v:predecessor").first;

  boost::dijkstra_shortest_paths(mesh, sv, predecessor_map(predecessor));

  std::vector<Point> points;
  do {
    points.push_back(get(ppmap,tv));
    tv = predecessor[tv];
  } while(predecessor[tv]!=tv);

  points.push_back(get(ppmap,sv));
  std::cout << points.size() << std::endl;
  BOOST_FOREACH(Point p, points){
    std::cout << p << std::endl;
  }
  return 0;
}
