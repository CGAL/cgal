#ifndef CGAL_NEF_NARY_UNION_BY_QUEUE_H
#define CGAL_NEF_NARY_UNION_BY_QUEUE_H

#include <CGAL/Nef_3/Nary_union_strategy.h>
#include <list>

namespace CGAL {

template<class Polyhedron>
class Nary_union_by_queue {

  std::list<Polyhedron> queue;
  typedef typename std::list<Polyhedron>::iterator pit;

 public:
  void add_polyhedron(const Polyhedron& P) {
    queue.push_back(P);
    if(queue.size()%100 == 0)
        std::cerr << "queue size: " << queue.size() << std::endl;
  }

  Polyhedron get_union() {

    while(queue.size() > 1) {
      if(queue.size()%100 == 0 || queue.size()<30)
        std::cerr << queue.size() << " polyhedra in the queue" << std::endl;

      pit i1(queue.begin()), i2(i1);
      ++i2;
      Polyhedron tmp(*i1 + *i2);

      queue.pop_front();
      queue.pop_front();
      queue.push_back(tmp);
    }
    return queue.front();
  }
};

} //namespace CGAL
#endif // CGAL_NEF_NARY_UNION_BY_QUEUE_H
