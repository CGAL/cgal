#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <list>

#include <boost/graph/kruskal_min_spanning_tree.hpp>


typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;

// The BGL makes heavy use of indices associated to the vertices
// We use a std::map to store the index

typedef std::map<vertex_descriptor,int> Vertex_index_map;
Vertex_index_map vertex_index_map;

// A std::map is not a property map, because it is not lightweight
typedef boost::associative_property_map<Vertex_index_map> Vertex_index_pmap;
Vertex_index_pmap vertex_index_pmap(vertex_index_map);

void
kruskal(const Polyhedron& P)
{
  // associate indices to the vertices
  vertex_iterator vb, ve;
  int index = 0;

  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vb and ve
  for(boost::tie(vb, ve)=vertices(P); vb!=ve; ++vb){
    vertex_index_pmap[*vb]= index++;
  }


  // We use the default edge weight which is the length of the edge
  // This property map is defined in graph_traits_Polyhedron_3.h

  // In the function call you can see a named parameter: vertex_index_map
  std::list<edge_descriptor> mst;

  boost::kruskal_minimum_spanning_tree(P,
                                       std::back_inserter(mst),
                                       boost::vertex_index_map(vertex_index_pmap));

  std::cout << "#VRML V2.0 utf8\n"
    "Shape {\n"
    "  appearance Appearance {\n"
    "    material Material { emissiveColor 1 0 0}}\n"
    "    geometry\n"
    "    IndexedLineSet {\n"
    "      coord Coordinate {\n"
    "        point [ \n";

  for(boost::tie(vb, ve) = vertices(P); vb!=ve; ++vb){
    std::cout <<  "        " << (*vb)->point() << "\n";
  }

  std::cout << "        ]\n"
               "     }\n"
    "      coordIndex [\n";

  for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it)
  {
    edge_descriptor e = *it ;
    vertex_descriptor s = source(e,P);
    vertex_descriptor t = target(e,P);
    std::cout << "      " << vertex_index_pmap[s] << ", " << vertex_index_pmap[t] <<  ", -1\n";
  }

  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";
}


int main() {

  Polyhedron P;
  Point a(1,0,0);
  Point b(0,1,0);
  Point c(0,0,1);
  Point d(0,0,0);

  P.make_tetrahedron(a,b,c,d);

  kruskal(P);

  return 0;
}
