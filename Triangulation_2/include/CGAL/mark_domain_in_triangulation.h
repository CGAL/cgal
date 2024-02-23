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

#ifndef CGAL_MARK_DOMAIN_IN_TRIANGULATION_H
#define CGAL_MARK_DOMAIN_IN_TRIANGULATION_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/assertions.h>
#include <CGAL/Triangulation_2/internal/In_domain.h>
#include <CGAL/Unique_hash_map.h>

#include <boost/property_map/property_map.hpp>

#include <deque>
#include <list>

namespace CGAL {

namespace internal {

template <typename CT, typename InDomainPmap>
void
mark_domain_in_triangulation(CT& ct,
                             Unique_hash_map<typename CT::Face_handle,int>& nesting_level,
                             typename CT::Face_handle start,
                             int index,
                             std::list<typename CT::Edge>& border,
                             InDomainPmap ipm)
{
  typedef typename CT::Face_handle Face_handle;
  typedef typename CT::Edge Edge;

  CGAL_precondition(ct.dimension() == 2);

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

} // namespace internal


//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
template <typename CT, typename InDomainPmap>
void
mark_domain_in_triangulation(CT& cdt, InDomainPmap ipm)
{
  typedef typename CT::Face_handle Face_handle;
  typedef typename CT::Edge Edge;

  Unique_hash_map<Face_handle,int> nesting_level(-1, cdt.number_of_faces());

  for(Face_handle f : cdt.all_face_handles()){
    put(ipm, f, false);
  }

  std::list<Edge> border;
  internal::mark_domain_in_triangulation(cdt, nesting_level, cdt.infinite_face(), 0, border, ipm);
  while(! border.empty()){
    Edge e = border.front();
    border.pop_front();
    Face_handle n = e.first->neighbor(e.second);
    if(nesting_level[n] == -1){
      internal::mark_domain_in_triangulation(cdt, nesting_level, n, nesting_level[e.first]+1, border, ipm);
    }
  }
}


template <typename CT>
void
mark_domain_in_triangulation(CT& cdt)
{
  internal::In_domain<CT> in_domain;
  mark_domain_in_triangulation(cdt, in_domain);
}

} // namespace CGAL

#endif // CGAL_MARK_DOMAIN_IN_TRIANGULATION_H
