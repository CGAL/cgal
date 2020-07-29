#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <boost/graph/breadth_first_search.hpp>

#include <fstream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;


int main(int argc, char** argv) {

  Polyhedron P;
  std::ifstream in((argc>1)?argv[1]:"cube.off");
  in >> P ;

  // associate indices to the vertices using the "id()" field of the vertex.
  vertex_iterator vb, ve;
  int index = 0;

  // boost::tie assigns the first and second element of the std::pair
  // returned by boost::vertices to the variables vit and ve
  for(boost::tie(vb,ve)=vertices(P); vb!=ve; ++vb ){
    vertex_descriptor  vd = *vb;
    vd->id() = index++;
  }

  // This is the vector where the distance gets written to
  std::vector<int> distance(P.size_of_vertices());


  // Here we start at an arbitrary vertex
  // Any other vertex could be the starting point
  boost::tie(vb,ve)=vertices(P);
  vertex_descriptor  vd = *vb;

  std::cout << "We compute distances to " << vd->point() << std::endl;


  // bfs = breadth first search explores the graph
  // Just as the distance_recorder there is a way to record the predecessor of a vertex
  boost::breadth_first_search(P,
                              vd,
                              visitor(boost::make_bfs_visitor(boost::record_distances(make_iterator_property_map(distance.begin(),
                                                                                                                 get(boost::vertex_index, P)),
                                                                                      boost::on_tree_edge()))));



  // Traverse all vertices and show at what distance they are
  for(boost::tie(vb,ve)=vertices(P); vb!=ve; ++vb ){
    vd = *vb;
    std::cout <<  vd->point() << "  is " << distance[vd->id()] << " hops away" << std::endl;
  }

  return 0;
}
