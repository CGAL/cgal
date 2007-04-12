#ifndef CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H
#define CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H

#ifdef CGAL_NEF_FAST_UNION
#include <CGAL/Nef_3/Binary_union.h>
#endif
#include <list>

CGAL_BEGIN_NAMESPACE

template<class Polyhedron>
class Nary_union_by_small_queue {

  int inserted;
  std::list<Polyhedron> queue;
  typedef typename std::list<Polyhedron>::iterator pit;

 public:
  Nary_union_by_small_queue() : inserted(0) {}
  
  void unite() {
    pit i1(queue.begin()), i2(i1);
    ++i2;

#ifdef CGAL_NEF_FAST_UNION
    CGAL::Binary_union<Polyhedron> bu(*i1, *i2);
    Polyhedron tmp;
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
    for(int i=2;(inserted%i) == 0; i*=2) {
      unite();
    }
  }

  Polyhedron get_union() {

    while(queue.size() > 1)
      unite();
    return queue.front();
  }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_NARY_UNION_BY_SMALL_QUEUE_H
