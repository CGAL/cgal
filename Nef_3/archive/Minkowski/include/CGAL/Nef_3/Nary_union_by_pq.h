#ifndef CGAL_NEF_NARY_UNION_BY_PQ_H
#define CGAL_NEF_NARY_UNION_BY_PQ_H

#include <CGAL/Nef_3/Nary_union_strategy.h>
#include <list>

namespace CGAL {

template<class Polyhedron>
class Nary_union_by_pq {

  typedef typename std::multimap<typename Polyhedron::Size_type, Polyhedron> PQ;
  typedef typename PQ::iterator PQ_iterator;
  PQ pq;

 public:
  void add_polyhedron(const Polyhedron& P) {
    pq.insert(make_pair(P.number_of_vertices(),P));
#ifndef NDEBUG
    if(pq.size()%100 == 0 || pq.size()<30)
        std::cerr << "pq size: " << pq.size() << std::endl;
#endif
  }

  Polyhedron get_union() {

    PQ_iterator i1, i2;
    while(pq.size() > 1) {
      i1 = i2 = pq.begin();
      ++i2;

#ifndef NDEBUG
      if(pq.size()%100 == 0 || pq.size()<30)
        std::cerr << pq.size() << " polyhedra in the priority queue "
                  << i1->first << "," << i2->first << std::endl;
 #endif

      Polyhedron N1(i1->second);
      Polyhedron N2(i2->second);
      Polyhedron Ntmp(N1 + N2);

      pq.erase(i1);
      pq.erase(i2);
      pq.insert(make_pair(Ntmp.number_of_vertices(),Ntmp));
    }

    return pq.begin()->second;
  }
};

} //namespace CGAL
#endif // CGAL_NEF_NARY_UNION_BY_PQ_H
