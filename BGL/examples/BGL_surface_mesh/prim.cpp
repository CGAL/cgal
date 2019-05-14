#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int, char* argv[]) 
{
  Mesh P;
  //std::cin >> P;
  std::ifstream in(argv[1]);
  in >> P;
  Mesh::Property_map<vertex_descriptor,vertex_descriptor> predecessor;
  predecessor = P.add_property_map<vertex_descriptor,vertex_descriptor>("v:predecessor").first;

  boost::prim_minimum_spanning_tree(P, predecessor, boost::root_vertex(*vertices(P).first));

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

  BOOST_FOREACH(vertex_descriptor vd, vertices(P)){
    std::cout <<  "        " << P.point(vd) << "\n";
  }
  
  std::cout << "        ]\n"
    "     }\n"
    "      coordIndex [\n"; 
  BOOST_FOREACH(vertex_descriptor vd, vertices(P)){
    if(predecessor[vd]!=vd){
      std::cout << "      " << std::size_t(vd) << ", " << std::size_t(predecessor[vd]) <<  ", -1\n";
    }
  }
  
  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";

  P.remove_property_map(predecessor);
  return 0;
}
