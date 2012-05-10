// Copyright (c) 2012 Geometry Factory. All rights reserved.
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
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_POLYLINE_SIMPLIFICATION_2_MARK_VERTICES_UNREMOVBLE_H
#define CGAL_POLYLINE_SIMPLIFICATION_2_MARK_VERTICES_UNREMOVBLE_H

namespace CGAL {

template < class Tr >
class Polyline_constrained_triangulation_2;

namespace Polyline_simplification_2 {

// Fix the leftmost, rightmost, topmost and bottommost vertex
template <typename Tr>
void
mark_vertices_unremovable(CGAL::Polyline_constrained_triangulation_2<Tr>& pct,
                          typename CGAL::Polyline_constrained_triangulation_2<Tr>::Constraint_id cid)
{
  typedef typename CGAL::Polyline_constrained_triangulation_2<Tr> PCT;
  typedef typename PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
  typename PCT::Vertex_handle l,r,b,t;
  l = r = b = t = *(pct.vertices_in_constraint_begin(cid));
  for(Vertices_in_constraint_iterator it = pct.vertices_in_constraint_begin(cid);
        it != pct.vertices_in_constraint_end(cid);
        it++){ 
    if(it->point.x() < l->point().x()) l = *it; 
    if(it->point.x() > r->point().x()) r = *it; 
    if(it->point.y() < b->point().y()) b = *it; 
    if(it->point.y() > t->point().y()) t = *it; 
    }
    l->fixed = r->fixed = t->fixed = b->fixed = true;
}


// Fix the leftmost, rightmost, topmost and bottommost vertex
template <typename Tr>
void
mark_vertices_unremovable(CGAL::Polyline_constrained_triangulation_2<Tr>& pct)
{
  typedef typename CGAL::Polyline_constrained_triangulation_2<Tr> PCT;
  typedef typename PCT::Constraint_iterator Constraint_iterator;
  typedef typename PCT::Constraint_id Constraint_id;

  Constraint_iterator cit = pct.constraints_begin(), e = pct.constraints_end();
  for(; cit!=e; ++cit){
    Constraint_id cid = *cit;
    mark_vertices_unremovable(pct,cid);
  }
}

} // namespace polyline_simplification_2
} // namespace CGAL 

#endif

