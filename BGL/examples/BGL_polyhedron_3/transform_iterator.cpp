#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <fstream>

#include <boost/iterator/transform_iterator.hpp>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef CGAL::Halfedge_around_target_iterator<Polyhedron> halfedge_around_target_iterator;


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

int main(int, char** argv)
{ 
  std::ifstream in(argv[1]);
  Polyhedron P;
  in >> P;
  GraphTraits::vertex_descriptor vd = *(vertices(P).first);

  typedef boost::transform_iterator<Source<Polyhedron>,halfedge_around_target_iterator> adjacent_vertex_iterator; 

  halfedge_around_target_iterator hb,he;
  boost::tie(hb,he) = halfedges_around_target(halfedge(vd,P),P);
  adjacent_vertex_iterator avib, avie;
  avib = boost::make_transform_iterator(hb, Source<Polyhedron>(P));
  avie = boost::make_transform_iterator(he, Source<Polyhedron>(P));
  
  std::list<vertex_descriptor> V;
  std::copy(avib,avie, std::back_inserter(V));
  return 0;
}
