// Copyright(c) 2022  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_MARK_DOMAINS_H
#define CGAL_MARK_DOMAINS_H

#include <CGAL/license/Triangulation_2.h>

#include <list>
#include <deque>

#include <CGAL/Unique_hash_map.h>

namespace CGAL {


  template <typename CDT, typename InDomainPmap>
void
mark_domains(CDT& ct,
             Unique_hash_map<typename CDT::Face_handle,int>& nesting_level,
             typename CDT::Face_handle start,
             int index,
             std::list<typename CDT::Edge>& border,
             InDomainPmap ipm)
{
  typedef typename CDT::Face_handle Face_handle;
  typedef typename CDT::Edge Edge;

  if(nesting_level[start] != -1){
    return;
  }
  std::list<Face_handle> queue;
  queue.push_back(start);

  while(! queue.empty()){
    Face_handle fh = queue.front();
    queue.pop_front();
    if(nesting_level[fh] == -1){
      nesting_level[fh] = index;
      if(index %2 == 1){
        put(ipm, fh, true);
      }
      for(int i = 0; i < 3; i++){
        Edge e(fh,i);
        Face_handle n = fh->neighbor(i);
        if(nesting_level[n] == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
template <typename CDT, typename InDomainPmap>
void
mark_domains(CDT& cdt, InDomainPmap ipm)
{
  typedef typename CDT::Face_handle Face_handle;
  typedef typename CDT::Edge Edge;

  Unique_hash_map<Face_handle,int> nesting_level(-1, cdt.number_of_faces());

  for(Face_handle f : cdt.all_face_handles()){
    put(ipm, f, false);
  }

  std::list<Edge> border;
  mark_domains(cdt, nesting_level, cdt.infinite_face(), 0, border, ipm);
  while(! border.empty()){
    Edge e = border.front();
    border.pop_front();
    Face_handle n = e.first->neighbor(e.second);
    if(nesting_level[n] == -1){
      mark_domains(cdt, nesting_level, n, nesting_level[e.first]+1, border, ipm);
    }
  }
}

} // namespace CGAL

#endif // CGAL_MARK_DOMAINS_H
