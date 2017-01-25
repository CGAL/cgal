// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ 
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_MS3_BIPARTITE_NARY_UNION_SORTED_COMBINED_H
#define CGAL_MS3_BIPARTITE_NARY_UNION_SORTED_COMBINED_H

#include <CGAL/license/Minkowski_sum_3.h>


#include <CGAL/Minkowski_sum_3/Gaussian_map.h>
#include <CGAL/Minkowski_sum_3/Gaussian_map_to_nef_3.h>
#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
#include <CGAL/Nef_3/Nary_union_by_summup.h>
#elif defined CGAL_NEF3_NARY_UNION_VIA_PQ
#include <CGAL/Nef_3/Nary_union_by_pq.h>
#elif defined CGAL_NEF3_NARY_UNION_BY_QUEUE
#include <CGAL/Nef_3/Nary_union_by_queue.h>
#elif defined CGAL_NEF3_NARY_UNION_USING_DU
#include <CGAL/Nef_3/Nary_union_using_distinct_uniter.h>
#elif defined CGAL_NEF3_NARY_UNION_BY_SMALL_QUEUE
#include <CGAL/Nef_3/Nary_union_by_small_queue.h>
#else
#include <CGAL/Nef_nary_union_3.h>
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

  typedef typename Nef_polyhedron::SM_const_decorator 
    SM_const_decorator;
  typedef typename Nef_polyhedron::Kernel Kernel;
  typedef typename Nef_polyhedron::Point_3 Point_3;
  typedef typename Nef_polyhedron::Volume_const_iterator  Volume_const_iterator;
  typedef typename Nef_polyhedron::SFace_const_handle  SFace_const_handle;
  typedef typename Nef_polyhedron::SHalfedge_const_handle  SHalfedge_const_handle;
  typedef typename Nef_polyhedron::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Nef_polyhedron::Halfedge_const_iterator  Halfedge_const_iterator;
  typedef typename Nef_polyhedron::Halffacet_const_iterator  Halffacet_const_iterator;

  typedef CGAL::Gaussian_map<Kernel, Nef_polyhedron> Gaussian_map;

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
#elif defined CGAL_NEF3_NARY_UNION_BY_SMALL_QUEUE
  CGAL::Nary_union_by_small_queue<Nef_polyhedron>
    nary_union;
#else
  CGAL::Nef_nary_union_3<Nef_polyhedron> nary_union;
#endif

  typedef std::pair<Gaussian_map, Point_3> GMapPoint;
  typedef typename std::list<GMapPoint> GM_list;
  typedef typename GM_list::const_iterator GM_iterator;
  typedef typename std::pair<GM_iterator, GM_iterator> GM_pair;
#ifdef NARY_SORT
  typedef typename std::multimap<Point_3, GM_pair, NARY_SORT> PQ;
#else
  typedef typename std::multimap<Point_3, GM_pair> PQ;
#endif
  typedef typename PQ::iterator PQ_iterator;

  GM_list GM0;
  Volume_const_iterator c0;
  for(c0 = ++N0.volumes_begin(); 
      c0 != N0.volumes_end(); ++c0) {
    if(!c0->mark()) continue;
    Point_3 p(SFace_const_handle(c0->shells_begin())->center_vertex()->point());
    GM0.push_back(std::make_pair(Gaussian_map(N0, c0), p));
  }

  GM_list GM1;
  Volume_const_iterator c1;
  for(c1 = ++N1.volumes_begin(); 
      c1 != N1.volumes_end(); ++c1) {
    if(!c1->mark()) continue;
    Point_3 p(SFace_const_handle(c1->shells_begin())->center_vertex()->point());
    GM1.push_back(std::make_pair(Gaussian_map(N1, c1), p));
  }

  CGAL_assertion_msg(!GM0.empty() || !GM1.empty(), 
		     "one operand must be full-dimensional");

  Vertex_const_iterator vci;
  for(vci = N0.vertices_begin();
      vci != N0.vertices_end(); ++vci)
    if(vci->number_of_sfaces() == 1 &&
       vci->number_of_svertices() == 0 &&
       vci->mark() &&
       !vci->sfaces_begin()->mark())
      GM0.push_back(std::make_pair(Gaussian_map(vci), vci->point()));
  for(vci = N1.vertices_begin();
      vci != N1.vertices_end(); ++vci)
    if(vci->number_of_sfaces() == 1 &&
       vci->number_of_svertices() == 0 &&
       vci->mark() &&
       !vci->sfaces_begin()->mark())
      GM1.push_back(std::make_pair(Gaussian_map(vci), vci->point()));
  
  Halfedge_const_iterator eci;
  for(eci = N0.halfedges_begin();
      eci != N0.halfedges_end(); ++eci) {
    if(eci->is_twin()) continue;
    SM_const_decorator SD(&*eci->source());
    if(!SD.is_isolated(eci)) continue;
    if(eci->source()->point() < eci->twin()->source()->point())
      GM0.push_back(std::make_pair(Gaussian_map(eci), eci->source()->point()));
    else
      GM0.push_back(std::make_pair(Gaussian_map(eci), eci->twin()->source()->point()));
  }
  
  for(eci = N1.halfedges_begin();
      eci != N1.halfedges_end(); ++eci) {
    if(eci->is_twin()) continue;
    SM_const_decorator SD(&*eci->source());
    if(!SD.is_isolated(eci)) continue;
    if(eci->source()->point() < eci->twin()->source()->point())
      GM1.push_back(std::make_pair(Gaussian_map(eci), eci->source()->point()));
    else
      GM1.push_back(std::make_pair(Gaussian_map(eci), eci->twin()->source()->point()));
  }

  Halffacet_const_iterator fci;
  for(fci = N0.halffacets_begin();
      fci != N0.halffacets_end(); ++fci) {
    if(fci->is_twin()) continue;  
    if( fci->incident_volume() != fci->twin()->incident_volume() &&
        ( fci->incident_volume()->mark() || fci->twin()->incident_volume()->mark() )) 
	{
		continue;
	}
    SHalfedge_const_handle se(fci->facet_cycles_begin());
    GM0.push_back(std::make_pair(Gaussian_map(fci), 
				 se->source()->source()->point()));    
  }
  for(fci = N1.halffacets_begin();
      fci != N1.halffacets_end(); ++fci) {
    if(fci->is_twin()) continue;  
    if( fci->incident_volume() != fci->twin()->incident_volume() &&
        ( fci->incident_volume()->mark() || fci->twin()->incident_volume()->mark() )) 
	{
		continue;
	}
    SHalfedge_const_handle se(fci->facet_cycles_begin());
    GM1.push_back(std::make_pair(Gaussian_map(fci), 
				 se->source()->source()->point()));    
  }

  PQ pq;
  GM_iterator gi0, gi1;
  for(gi0 = GM0.begin(); gi0 != GM0.end(); ++gi0) {
    for(gi1 = GM1.begin(); gi1 != GM1.end(); ++gi1) {
      pq.insert(std::make_pair(gi0->second+(CGAL::ORIGIN-gi1->second),
			       std::make_pair(gi0, gi1)));

    }
  }

  CGAL_assertion_msg(!GM0.empty() && !GM1.empty(), 
		     "one operand must be full-dimensional");

  Nef_polyhedron empty;
  while(pq.size() > 0) {
    PQ_iterator pqi(pq.begin());

    Gaussian_map GcG;
    GcG.minkowski_sum(pqi->second.first->first, 
		      pqi->second.second->first);
    pq.erase(pqi);
    Nef_polyhedron Ntmp(empty);
    CGAL::Gaussian_map_to_nef_3<Nef_polyhedron> Convertor(GcG);
    Ntmp.delegate(Convertor, true);
    CGAL_assertion(Ntmp.is_valid());
    nary_union.add_polyhedron(Ntmp);
    delete GcG.sphere_map();
  }

  // clean up the spherical_mapS
  for(GM_iterator it = GM0.begin(); it != GM0.end(); ++it) {
    delete it->first.sphere_map();
  }

  for(GM_iterator it = GM1.begin(); it != GM1.end(); ++it) {
    delete it->first.sphere_map();
  }
  
  return nary_union.get_union();
}

} //namespace CGAL
#endif // CGAL_MS3_BIPARTITE_NARY_UNION_SORTED_COMBINED_H
