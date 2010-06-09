#ifndef CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_SEPARATELY_H
#define CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_SEPARATELY_H

#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_nef_3.h>

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

#include <list>

namespace CGAL {

template<class Const_decorator>
struct Sort_volumes_by_smallest_vertex {
  typedef typename Const_decorator::Volume_const_handle Volume_handle;
  typedef typename Const_decorator::SFace_const_handle SFace_handle;
  Sort_volumes_by_smallest_vertex() {}
  bool operator()(Volume_handle c0, Volume_handle c1) {
    return
      SFace_handle(c0->shells_begin())->center_vertex()->point() <
      SFace_handle(c1->shells_begin())->center_vertex()->point();
  }
};

template<typename Nef_polyhedron>
Nef_polyhedron
bipartite_nary_union_sorted_separately(Nef_polyhedron& N0,  
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

  CGAL::Timer t1, t2, t3;

  t1.start();
  std::list<Volume_const_handle> volumes;
  for(Volume_const_iterator c = ++N0.volumes_begin();
      c!=N0.volumes_end();++c) {  
    if(c->mark())
      volumes.push_back(c);
  }
  volumes.sort(Sort_volumes_by_smallest_vertex<Nef_polyhedron>());
  t1.stop();

  t2.start();
  std::list<Gausian_map> GM;
  int shells = volumes.size();
  typename std::list<Volume_const_handle>::iterator 
    lci(volumes.begin());
  for(;lci!=volumes.end();++lci) {
    std::cerr << "noch " << --shells << " sorted shells" << std::endl;
    GM.push_back(Gausian_map(N0, *lci));
  }
  t2.stop();

  t1.start();
  volumes.clear();
  for(Volume_const_iterator c = ++N1.volumes_begin();
      c!=N1.volumes_end();++c) {  
    if(c->mark())
      volumes.push_back(c);
  }
  volumes.sort(Sort_volumes_by_smallest_vertex<Nef_polyhedron>());
  t1.stop();

  t2.start();
  shells = volumes.size()+1;
  for(lci=volumes.begin();lci!=volumes.end();++lci) {
    --shells;
    Gausian_map G(N1, *lci);
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
      t2.stop();
      CGAL_assertion(Ntmp.is_valid());

      t3.start();
      nary_union.add_polyhedron(Ntmp);
      t3.stop();
      t2.start();
    }
  }
  t2.stop();

  t3.start();
  Nef_polyhedron result = nary_union.get_union();
  t3.stop();

  std::cerr << "Sorting volumes : " << t1.time() << std::endl;
  std::cerr << "Convex Minkowski sums : " << t2.time() << std::endl;
  std::cerr << "Union of subpolyhedra: " << t3.time() << std::endl;

  return result;
}

} //namespace CGAL
#endif // CGAL_NEF_BIPARTITE_NARY_UNION_SORTED_SEPARATELY_H
