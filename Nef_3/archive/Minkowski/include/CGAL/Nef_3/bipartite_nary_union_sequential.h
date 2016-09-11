#ifndef CGAL_NEF_BIPARTITE_NARY_UNION_SEQUENTIAL_H
#define CGAL_NEF_BIPARTITE_NARY_UNION_SEQUENTIAL_H

#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_nef_3.h>

#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
#include <CGAL/Nef_3/Nary_union_by_summup.h>
#elif defined CGAL_NEF3_NARY_UNION_VIA_PQ
#include <CGAL/Nef_3/Nary_union_by_pq.h>
#elif defined CGAL_NEF3_NARY_UNION_BY_QUEUE
#include <CGAL/Nef_3/Nary_union_by_queue.h>
#elif defined CGAL_NEF3_NARY_UNION_BY_SHUFFLED_QUEUE
#include <CGAL/Nef_3/Nary_union_by_shuffled_queue.h>
#elif defined CGAL_NEF3_NARY_UNION_USING_DU
#include <CGAL/Nef_3/Nary_union_using_distinct_uniter.h>
#else
#include <CGAL/Nef_3/Nary_union_by_small_queue.h>
#endif

namespace CGAL {

template<typename Nef_polyhedron>
Nef_polyhedron 
bipartite_nary_union_sequential(Nef_polyhedron& N0,  
				 Nef_polyhedron& N1) {

  typedef typename Nef_polyhedron::Kernel Kernel;
  typedef typename Nef_polyhedron::Volume_const_iterator  Volume_const_iterator;
  typedef typename Nef_polyhedron::Volume_const_handle  Volume_const_handle;

  typedef CGAL::Gausian_map<Kernel> Gausian_map;

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

  CGAL::Timer t1, t2;

  t1.start();
  std::list<Gausian_map> GM;
  int shells = N0.number_of_volumes();
  Volume_const_iterator c0 = ++N0.volumes_begin();
  for(;c0!=N0.volumes_end();++c0) {
    std::cerr << "noch " << --shells << " shells" << std::endl;
    if(c0->mark() == false) continue;
    GM.push_back(Gausian_map(N0, c0));
  }

  shells = N1.number_of_volumes();
  Volume_const_iterator c1 = ++N1.volumes_begin();
  for(;c1!=N1.volumes_end();++c1) {
    --shells;
    if(c1->mark() == false) continue;
    Gausian_map G(N1, c1);
    
    int maps = GM.size()+1;
    typename std::list<Gausian_map>::const_iterator gi;
    for(gi = GM.begin(); gi != GM.end(); ++gi) {
      std::cerr << "noch "
		<< shells << ", " 
		<< --maps << " convex sums" << std::endl;
      Gausian_map GcG;
      GcG.minkowski_sum(G,*gi);
      Nef_polyhedron Ntmp;
      CGAL::gausian_map_to_nef_3<Nef_polyhedron> Converter(GcG);
      Ntmp.delegate(Converter, true);
      CGAL_assertion(Ntmp.is_valid());
      t1.stop();
      t2.start();
      nary_union.add_polyhedron(Ntmp);
      t2.stop();
      t1.start();
    }
  }
  t1.stop();

  t2.start();
  Nef_polyhedron result = nary_union.get_union();
  t2.stop();

  std::cerr << "Convex Minkowski sums : " << t1.time() << std::endl;
  std::cerr << "Union of subpolyhedra: " << t2.time() << std::endl;

  return result;
}

} //namespace CGAL
#endif // CGAL_NEF_BIPARTITE_NARY_UNION_SEQUENTIAL_H
