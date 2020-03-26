#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <boost/graph/breadth_first_search.hpp>

#include <fstream>

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef Kernel::Point_3                             Point;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;

typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, 3, LCC_traits>::type LCC;

typedef boost::graph_traits<LCC>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<LCC>::vertex_iterator   vertex_iterator;

int main(int argc, char** argv)
{
  LCC lcc;
  CGAL::read_off((argc>1)?argv[1]:"cube.off", lcc);

  // This is the vector where the distance gets written to
  std::vector<int> distance(lcc.vertex_attributes().size());

  // Here we start at an arbitrary vertex
  // Any other vertex could be the starting point
  vertex_iterator vb, ve;
  boost::tie(vb,ve)=vertices(lcc);
  vertex_descriptor  vd = *vb;

  std::cout << "We compute distances to " << vd->point() << std::endl;

  // bfs = breadth first search explores the graph
  // Just as the distance_recorder there is a way to record the predecessor of a vertex
  boost::breadth_first_search(lcc,
                              vd,
                              visitor(boost::make_bfs_visitor
                                      (boost::record_distances
                                       (make_iterator_property_map
                                        (distance.begin(),
                                         get(boost::vertex_index, lcc)),
                                        boost::on_tree_edge()))));

  // Traverse all vertices and show at what distance they are
  for(boost::tie(vb,ve)=vertices(lcc); vb!=ve; ++vb)
  {
    vd = *vb;
    std::cout<<vd->point()<<"  is "<<distance[vd->id()]<<" hops away."<<std::endl;
  }

  return 0;
}
