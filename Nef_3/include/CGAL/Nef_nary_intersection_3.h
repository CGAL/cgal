#ifndef CGAL_NEF_NARY_INTERSECTION_H
#define CGAL_NEF_NARY_INTERSECTION_H

#include <list>

CGAL_BEGIN_NAMESPACE

template<class Polyhedron>
class Nary_intersection {

  int inserted;
  std::list<Polyhedron> queue;
  typedef typename std::list<Polyhedron>::iterator pit;
  Polyhedron empty;

 public:
  Nary_intersection() : inserted(0) {}
  
  void intersect() {
    pit i1(queue.begin()), i2(i1);
    ++i2;

    Polyhedron tmp(*i1 * *i2);

    queue.pop_front();
    queue.pop_front();
    queue.push_front(tmp);
  }

  void add_polyhedron(const Polyhedron& P) {
    queue.push_front(P);
    ++inserted;
    for(int i=2;(inserted%i) == 0; i*=2) {
      intersect();
    }
  }

  Polyhedron get_intersection() {

    while(queue.size() > 1)
      intersect();
    inserted = 0;
    return queue.front();
  }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF_NARY_INTERSECTION_H
