#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <fstream>

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> LCC_traits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, LCC_traits>::type LCC;

typedef boost::graph_traits<LCC>                   GraphTraits;
typedef GraphTraits::vertex_descriptor             vertex_descriptor;
typedef GraphTraits::halfedge_descriptor           halfedge_descriptor;
typedef CGAL::Halfedge_around_target_iterator<LCC> halfedge_around_target_iterator;


template <typename G>
struct Source {
  const G* g;

  Source()
    : g(NULL)
  {}

  Source(const G& g)
    : g(&g)
  {}

  typedef typename boost::graph_traits<G>::vertex_descriptor result_type;
  typedef typename boost::graph_traits<G>::halfedge_descriptor argument_type;

  result_type operator()(argument_type h) const
  {
    return source(h, *g);
  }
};

int main(int argc, char** argv)
{
  LCC lcc;
  CGAL::read_off((argc>1)?argv[1]:"cube.off", lcc);
  GraphTraits::vertex_descriptor vd = *(vertices(lcc).first);

  typedef boost::transform_iterator<Source<LCC>,halfedge_around_target_iterator> adjacent_vertex_iterator;

  halfedge_around_target_iterator hb,he;
  boost::tie(hb,he) = halfedges_around_target(halfedge(vd,lcc),lcc);
  adjacent_vertex_iterator avib, avie;
  avib = boost::make_transform_iterator(hb, Source<LCC>(lcc));
  avie = boost::make_transform_iterator(he, Source<LCC>(lcc));

  std::list<vertex_descriptor> V;
  std::copy(avib,avie, std::back_inserter(V));
  return 0;
}
