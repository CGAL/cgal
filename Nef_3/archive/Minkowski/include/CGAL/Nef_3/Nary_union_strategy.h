#ifndef CGAL_NEF_NARY_UNION_STRATEGY_H
#define CGAL_NEF_NARY_UNION_STRATEGY_H

namespace CGAL {

template<class Polyhedron>
class Nary_union_strategy {
 public:
  void add_polyhedron(const Polyhedron& P) = 0;
  Polyhedron get_union() = 0;
};

} //namespace CGAL
#endif // CGAL_NEF_NARY_UNION_STRATEGY_H
