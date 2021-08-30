#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <iostream>
#include <list>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef Kernel::Point_3                             Point;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;

typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, LCC_traits>::type LCC;

typedef boost::graph_traits<LCC>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<LCC>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<LCC>::edge_descriptor   edge_descriptor;



void kruskal(const LCC& lcc)
{
  // We use the default edge weight which is the length of the edge
  // This property map is defined in graph_traits_Linear_cell_complex_for_combinatorial_map.h

  // This function call requires a vertex_index_map named parameter which
  // when  ommitted defaults to "get(vertex_index,graph)".
  // That default works here because the vertex type has an "id()" method
  // field which is used by the vertex_index internal property.
  std::list<edge_descriptor> mst;
  boost::kruskal_minimum_spanning_tree(lcc,std::back_inserter(mst));

  std::cout << "#VRML V2.0 utf8\n"
    "Shape {\n"
    "appearance Appearance {\n"
    "material Material { emissiveColor 1 0 0}}\n"
    "geometry\n"
    "IndexedLineSet {\n"
    "coord Coordinate {\n"
    "point [ \n";

  vertex_iterator vb, ve;
  for(boost::tie(vb,ve) = vertices(lcc); vb!=ve; ++vb){
    std::cout << (*vb)->point() << "\n";
  }

  std::cout << "]\n"
    "}\n"
    "coordIndex [\n";

  for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it){
    std::cout << source(*it,lcc)->id()
              << ", " << target(*it,lcc)->id() <<  ", -1\n";
  }

  std::cout << "]\n"
    "}#IndexedLineSet\n"
    "}# Shape\n";
}


int main()
{
  LCC lcc;

  Point a(1,0,0);
  Point b(0,1,0);
  Point c(0,0,1);
  Point d(0,0,0);

  lcc.make_tetrahedron(a,b,c,d);
  kruskal(lcc);

  return 0;
}
