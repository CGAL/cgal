#ifndef CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_COMBINED_H
#define CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_COMBINED_H

#include <CGAL/Nef_S2/Gaussian_map.h>
#include <CGAL/Nef_S2/gaussian_map_to_nef_3.h>
#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
#include <CGAL/Nef_3/Nary_union_by_summup.h>
#elif defined CGAL_NEF3_NARY_UNION_VIA_PQ
#include <CGAL/Nef_3/Nary_union_by_pq.h>
#elif defined CGAL_NEF3_NARY_UNION_BY_QUEUE
#include <CGAL/Nef_3/Nary_union_by_queue.h>
#elif defined CGAL_NEF3_NARY_UNION_USING_DU
#include <CGAL/Nef_3/Nary_union_using_distinct_uniter.h>
#else
#include <CGAL/Nef_3/Nary_union_by_small_queue.h>
#endif

namespace CGAL {

#ifdef CGAL_NEF_NARY_UNION_L1_SORT
template<typename Point_3>
struct L1_sort {
  bool operator()(Point_3 p0, Point_3 p1) {
    return
      p0.x() + p0.y() + p0.z() < p1.x() + p1.y() + p1.z();
  }
};
#define NARY_SORT L1_sort<Point_3>
#elif defined CGAL_NEF_NARY_UNION_L2_SORT
template<typename Point_3>
struct L2_sort {
  bool operator()(Point_3 p0, Point_3 p1) {
    return
      p0.x()*p0.x() + p0.y()*p0.y() + p0.z()*p0.z() <
      p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z();
  }
};
#define NARY_SORT L2_sort<Point_3>
#endif

template<typename Nef_polyhedron>
Nef_polyhedron
bipartite_nary_union_sorted_combined(Nef_polyhedron& N0,
                                     Nef_polyhedron& N1) {

  typedef typename Nef_polyhedron::Kernel Kernel;
  typedef typename Nef_polyhedron::Point_3 Point_3;
  typedef typename Nef_polyhedron::Volume_const_iterator  Volume_const_iterator;
  typedef typename Nef_polyhedron::Volume_const_handle  Volume_const_handle;
  typedef typename Nef_polyhedron::SFace_const_handle  SFace_const_handle;

  typedef CGAL::Gaussian_map<Kernel> Gaussian_map;

#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
 CGAL::Nary_union_by_summup<Nef_polyhedron>
   nary_union;
#elif defined CGAL_NEF3_NARY_UNION_VIA_PQ
  CGAL::Nary_union_by_pq<Nef_polyhedron>
    nary_union;
#elif defined CGAL_NEF3_NARY_UNION_USING_DU
  CGAL::Nary_union_using_distinct_uniter<Nef_polyhedron>
    nary_union;
#elif defined CGAL_NEF3_NARY_UNION_BY_QUEUE
  CGAL::Nary_union_by_queue<Nef_polyhedron>
    nary_union;
#else
  CGAL::Nary_union_by_small_queue<Nef_polyhedron>
    nary_union;
#endif

  CGAL::Timer t1, t2, t3, t4;

  typedef typename std::list<Gaussian_map> GM_list;
  typedef typename GM_list::const_iterator GM_iterator;
  typedef typename std::pair<GM_iterator, GM_iterator> GM_pair;
#ifdef NARY_SORT
  typedef typename std::multimap<Point_3, GM_pair, NARY_SORT> PQ;
#else
  typedef typename std::multimap<Point_3, GM_pair> PQ;
#endif
  typedef typename PQ::iterator PQ_iterator;

  GM_list GM0;
  int shells = N0.number_of_volumes();
  Volume_const_iterator c0;
  for(c0 = ++N0.volumes_begin();
      c0 != N0.volumes_end(); ++c0) {
    if(!c0->mark()) continue;
    std::cerr << "noch " << --shells << " shells" << std::endl;
    t2.start();
    GM0.push_back(Gaussian_map(N0, c0));
    t2.stop();
  }

  std::cerr << "GM0 done" << std::endl;

  GM_list GM1;
  shells = N1.number_of_volumes();
  Volume_const_iterator c1;
  for(c1 = ++N1.volumes_begin();
      c1 != N1.volumes_end(); ++c1) {
    if(!c1->mark()) continue;
    std::cerr << "noch " << --shells << " shells" << std::endl;
    t2.start();
    GM1.push_back(Gaussian_map(N1, c1));
    t2.stop();
  }

  std::cerr << "GM1 done" << std::endl;

  PQ pq;
  GM_iterator gi0, gi1;
  c0 = N0.volumes_begin();
  CGAL_assertion(!c0->mark());
  for(gi0 = GM0.begin(); gi0 != GM0.end(); ++gi0) {
    do { ++c0;} while(!c0->mark());
    c1 = N1.volumes_begin();
    CGAL_assertion(!c1->mark());
    for(gi1 = GM1.begin(); gi1 != GM1.end(); ++gi1) {
      do { ++c1;} while(!c1->mark());
      Point_3 p0(SFace_const_handle(c0->shells_begin())->center_vertex()->point());
      Point_3 p1(SFace_const_handle(c1->shells_begin())->center_vertex()->point());
      pq.insert(std::make_pair(p0 + (p1-CGAL::ORIGIN),
                               std::make_pair(gi0, gi1)));

    }
  }

  Nef_polyhedron empty;
  while(pq.size() > 0) {
    std::cerr << pq.size() << " Gaussian maps in priority queue" << std::endl;

    PQ_iterator pqi(pq.begin());

    t1.start();
    Gaussian_map GcG;
    GcG.minkowski_sum(*pqi->second.first,
                      *pqi->second.second);
    t1.stop();
    pq.erase(pqi);
    Nef_polyhedron Ntmp(empty);
    t3.start();
    CGAL::gaussian_map_to_nef_3<Nef_polyhedron> Converter(GcG);
    Ntmp.delegate(Converter, true);
    t3.stop();
    CGAL_assertion(Ntmp.is_valid());
    t4.start();
    nary_union.add_polyhedron(Ntmp);
    t4.stop();
  }

  t4.start();
  Nef_polyhedron result = nary_union.get_union();
  t4.stop();

  std::cerr << "c2GM : " << t2.time() << std::endl;
  std::cerr << "Convex Minkowski sums : " << t1.time() << std::endl;
  std::cerr << "GM2c : " << t3.time() << std::endl;
  std::cerr << "Union of subpolyhedra: " << t4.time() << std::endl;

  return result;
}

} //namespace CGAL
#endif // CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_COMBINED_H
