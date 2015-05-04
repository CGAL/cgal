#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>
#include <list>

#include <boost/graph/kruskal_min_spanning_tree.hpp>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;

void
kruskal(const Mesh& sm)
{
   // We use the default edge weight which is the squared length of the edge

  std::list<edge_descriptor> mst;

  boost::kruskal_minimum_spanning_tree(sm, 
                                       std::back_inserter(mst));

  std::cout << "#VRML V2.0 utf8\n"
    "Shape {\n"
    "  appearance Appearance {\n"
    "    material Material { emissiveColor 1 0 0}}\n"
    "    geometry\n"
    "    IndexedLineSet {\n"
    "      coord Coordinate {\n"
    "        point [ \n";

  vertex_iterator vb,ve;
  for(boost::tie(vb, ve) = vertices(sm); vb!=ve; ++vb){
    std::cout <<  "        " << sm.point(*vb) << "\n";
  }

  std::cout << "        ]\n"
               "     }\n"
    "      coordIndex [\n";

  for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it)
  {
    edge_descriptor e = *it ;
    vertex_descriptor s = source(e,sm);
    vertex_descriptor t = target(e,sm);
    std::cout << "      " << s << ", " << t <<  ", -1\n";
  }

  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";
}


int main(int,char** argv) {

  Mesh sm;
  std::ifstream input(argv[1]);
  input >> sm;

  kruskal(sm);

  return 0;
}
