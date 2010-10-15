#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <iostream>
#include <list>

#include <boost/graph/kruskal_min_spanning_tree.hpp>

typedef CGAL::Cartesian<double>                                      Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor   edge_descriptor;



void
kruskal( const Polyhedron& P)
{

  // We use the default edge weight which is the squared length of the edge
  // This property map is defined in graph_traits_Polyhedron_3.h
  
  // This function call requires a vertex_index_map named parameter which
  // when  ommitted defaults to "get(vertex_index,graph)".
  // That default works here because the vertex type supports the "id()"
  // field which is used by the vertex_index internal property.
  std::list<edge_descriptor> mst;
  boost::kruskal_minimum_spanning_tree(P,std::back_inserter(mst));
  
  std::cout << "#VRML V2.0 utf8\n"
    "Shape {\n"
    "appearance Appearance {\n"
    "material Material { emissiveColor 1 0 0}}\n"
    "geometry\n"
    "IndexedLineSet {\n"
    "coord Coordinate {\n"
    "point [ \n";
  
  vertex_iterator vb, ve;
  for(boost::tie(vb,ve) = boost::vertices(P); vb!=ve; ++vb){
    std::cout << (*vb)->point() << "\n";
  }
  
  std::cout << "]\n"
    "}\n"
    "coordIndex [\n";
  
  for(std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it){
    std::cout << boost::source(*it,P)->id()
	      << ", " << boost::target(*it,P)->id() <<  ", -1\n";
  }
  
  std::cout << "]\n"
    "}#IndexedLineSet\n"
    "}# Shape\n";
}


int main() {
  
  Polyhedron P;
  
  Point a(1,0,0);
  Point b(0,1,0);
  Point c(0,0,1);
  Point d(0,0,0);
  
  P.make_tetrahedron(a,b,c,d);

  // associate indices to the vertices using the "id()" field of the vertex.
  vertex_iterator vb, ve;
  int index = 0;
  
  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vb,ve)=boost::vertices(P); vb!=ve; ++vb ){
    vertex_descriptor  vd = *vb;
    vd->id() = index++;
  }
  
  kruskal(P);
  return 0;
}
