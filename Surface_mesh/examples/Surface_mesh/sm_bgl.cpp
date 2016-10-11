#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

// workaround a bug in Boost-1.54
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int /* argc */, char* argv[]) 
{
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;
  Mesh::Property_map<vertex_descriptor,vertex_descriptor> predecessor;
  predecessor = sm.add_property_map<vertex_descriptor,vertex_descriptor>("v:predecessor").first;

  boost::prim_minimum_spanning_tree(sm, predecessor, boost::root_vertex(*vertices(sm).first));

  std::cout << "#VRML V2.0 utf8\n"
    "DirectionalLight {\n"
    "direction 0 -1 0\n"
    "}\n"
    "Shape {\n"
    "  appearance Appearance {\n"
    "    material Material { emissiveColor 1 0 0}}\n"
    "    geometry\n"
    "    IndexedLineSet {\n"
    "      coord Coordinate {\n"
    "        point [ \n";

  BOOST_FOREACH(vertex_descriptor vd, vertices(sm)){
    std::cout <<  "        " << sm.point(vd) << "\n";
  }
  
  std::cout << "        ]\n"
    "     }\n"
    "      coordIndex [\n"; 
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm)){
    if(predecessor[vd]!=vd){
      std::cout << "      " << std::size_t(vd) << ", " << std::size_t(predecessor[vd]) <<  ", -1\n";
    }
  }
  
  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";

  sm.remove_property_map(predecessor);
  return 0;
}
