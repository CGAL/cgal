#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Surface_mesh<Point> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator             vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;

typedef boost::property_map<Polyhedron,CGAL::vertex_point_t>::type PPmap;


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor                  vertex_desc;
typedef boost::graph_traits<Graph>::vertex_iterator                    vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator                      edge_iter;

typedef std::map<vertex_desc, std::vector<vertex_descriptor> >         Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>            Correspondence_PMap;

typedef std::map<vertex_desc, Point>                                   GraphPointMap;
typedef boost::associative_property_map<GraphPointMap>                 GraphPointPMap;


// This example extracts a medially centered skeleton from a given mesh.
int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.is_empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return 1;
  }

  PPmap ppmap = get(CGAL::vertex_point, mesh);
  Graph g;
  GraphPointMap points_map;
  GraphPointPMap points(points_map);

  Correspondence_map corr_map;
  Correspondence_PMap corr(corr_map);


  CGAL::extract_mean_curvature_flow_skeleton(mesh, g, points, corr);
#if 0
  vertex_iterator vb, ve;

  std::cout << "vertices: " << boost::num_vertices(g) << "\n";
  std::cout << "edges: " << boost::num_edges(g) << "\n";

  // Output all the edges.
  edge_iter ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    Point s = points[source(*ei, g)];
    Point t = points[target(*ei, g)];
    std::cout << s << " " << t << "\n";
  }


  // Output skeletal points and the corresponding surface points
  vertex_iter gvb, gve;
  for (boost::tie(gvb, gve) = boost::vertices(g); gvb != gve; ++gvb)
  {
    vertex_desc i = *gvb;
    Point skel = points[i];
    std::cout << skel << ": ";

    for (size_t j = 0; j < corr[i].size(); ++j)
    {
      Point surf = ppmap[corr[i][j]];
      std::cout << surf << " ";
    }
    std::cout << "\n";
  }
#endif
  return 0;
}

