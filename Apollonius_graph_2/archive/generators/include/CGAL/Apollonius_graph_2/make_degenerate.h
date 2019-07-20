#ifndef CGAL_APOLLONIUS_GRAPH_2_MAKE_DEGENERATE_H
#define CGAL_APOLLONIUS_GRAPH_2_MAKE_DEGENERATE_H 1

#include <CGAL/basic.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>

namespace CGAL {

template<class InputIterator, class OutputIterator, class Traits>
OutputIterator
make_degenerate(InputIterator first,
		InputIterator beyond,
		OutputIterator oit,
		const Traits& tr = Traits())
{
  typedef CGAL::Apollonius_graph_hierarchy_2<Traits> Apollonius_graph;

  typedef typename Apollonius_graph::Site_2 Site_2;

  typedef typename Apollonius_graph::Finite_faces_iterator
    Finite_faces_iterator;
  
  Apollonius_graph ag;
  Site_2 site;

  for (InputIterator it = first; it != beyond; ++it) {
    ag.insert(*it);
  }

  for (Finite_faces_iterator f = ag.finite_faces_begin();
       f != ag.finite_faces_end(); ++f) {
    *oit++ = ag.dual(f);
  }
  return oit;
}

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_MAKE_DEGENERATE_H
