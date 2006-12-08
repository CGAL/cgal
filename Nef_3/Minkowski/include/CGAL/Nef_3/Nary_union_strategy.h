#ifndef CGAL_NEF_NARY_UNION_STRATEGY_H
#define CGAL_NEF_NARY_UNION_STRATEGY_H

CGAL_BEGIN_NAMESPACE

template<class Polyhedron>
class Nary_union_strategy {
 public:
  void add_polyhedron(const Polyhedron& P) = 0;
  Polyhedron get_union() = 0;
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_NARY_UNION_STRATEGY_H
