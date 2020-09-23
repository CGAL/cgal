#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/iterator.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;

typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, LCC_traits>::type LCC;

typedef boost::graph_traits<LCC>                   GraphTraits;
typedef GraphTraits::vertex_descriptor             vertex_descriptor;
typedef CGAL::Halfedge_around_target_iterator<LCC> halfedge_around_target_iterator;

template <typename OutputIterator>
OutputIterator
adjacent_vertices_V1(const LCC& g,
                     vertex_descriptor vd,
                     OutputIterator out)
{
  typename GraphTraits::halfedge_descriptor hb = halfedge(vd,g), done(hb);
  do
  {
    *out++ = source(hb,g);
    hb = opposite(next(hb,g),g);
  }
  while(hb!= done);
  return out;
}


template <typename OutputIterator>
OutputIterator adjacent_vertices_V2(const LCC& g,
                                    vertex_descriptor vd,
                                    OutputIterator out)
{
  halfedge_around_target_iterator hi, he;
  for(boost::tie(hi, he) = halfedges_around_target(halfedge(vd,g),g); hi != he; ++hi)
  {
    *out++ = source(*hi,g);
  }
  return out;
}


int main(int argc, char** argv)
{
  LCC lcc;
  CGAL::read_off((argc>1)?argv[1]:"cube.off", lcc);

  GraphTraits::vertex_iterator vi = vertices(lcc).first;
  std::list<vertex_descriptor> V;
  adjacent_vertices_V1(lcc, *vi, std::back_inserter(V));
  ++vi;
  adjacent_vertices_V2(lcc, *vi, std::back_inserter(V));
  std::cerr << "done\n";
  return 0;
}
