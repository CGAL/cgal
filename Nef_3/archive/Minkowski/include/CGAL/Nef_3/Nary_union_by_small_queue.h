#ifndef CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H
#define CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H

#ifdef CGAL_NEF_FAST_UNION
#include <CGAL/Nef_3/Binary_union.h>
#endif
#include <list>

namespace CGAL {

template<class Polyhedron>
class Nary_union_by_small_queue {

  int inserted;
  std::list<Polyhedron> queue;
  typedef typename std::list<Polyhedron>::iterator pit;
  Polyhedron empty;

 public:
  Nary_union_by_small_queue() : inserted(0) {}

  void unite() {
    pit i1(queue.begin()), i2(i1);
    ++i2;

#ifdef CGAL_NEF_FAST_UNION
    CGAL::Binary_union<Polyhedron> bu(*i1, *i2);
    Polyhedron tmp(empty);
    tmp.delegate(bu, false, false);
#else
    Polyhedron tmp(*i1 + *i2);
#endif

    queue.pop_front();
    queue.pop_front();
    queue.push_front(tmp);
  }

  void add_polyhedron(const Polyhedron& P) {
    queue.push_front(P);
    ++inserted;
    std::cerr << "inserted " << inserted << std::endl;
    for(int i=2;(inserted%i) == 0; i*=2) {
      unite();
    }
  }

  Polyhedron get_union() {

    while(queue.size() > 1)
      unite();
    inserted = 0;
    return queue.front();
  }
};

} //namespace CGAL
#endif // CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H
