#ifndef CGAL_NEF_NARY_UNION_BY_SUMMUP_H
#define CGAL_NEF_NARY_UNION_BY_SUMMUP_H

#include <CGAL/Nef_3/Nary_union_strategy.h>
#include <list>

CGAL_BEGIN_NAMESPACE

template<class Polyhedron>
class Nary_union_by_summup {

  Polyhedron sum;

 public:
  void add_polyhedron(const Polyhedron& P) {
    sum += P;
  }

  Polyhedron get_union() {
    return sum;
  }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_NARY_UNION_BY_SUMMUP_H
