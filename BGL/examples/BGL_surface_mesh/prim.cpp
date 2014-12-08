#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>
#include <list>

#include <boost/graph/prim_minimum_spanning_tree.hpp>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;

void
prim(Mesh& P)
{
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

  vertex_iterator vb,ve;
  for(boost::tie(vb, ve) = vertices(P); vb!=ve; ++vb){
    std::cout <<  "        " << P.point(*vb) << "\n";
  }

  std::cout << "        ]\n"
               "     }\n"
    "      coordIndex [\n"; 

  for(boost::tie(vb, ve) = vertices(P); vb!=ve; ++vb){
    if(predecessor[*vb]!=*vb){
      int s = *vb;
      int t = predecessor[*vb];
      std::cout << "      " << s << ", " << t <<  ", -1\n";
    }
  }

  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";
}


int main() {

  Mesh P;
  std::cin >> P;

  prim(P);

  return 0;
}
