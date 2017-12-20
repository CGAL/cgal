// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ophir Setter    <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_GPS_BFS_SCANNER_H
#define CGAL_GPS_BFS_SCANNER_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <queue>
#include <stack>

namespace CGAL {

template <class Arrangement_, class Visitor_>
class Gps_bfs_scanner
{
  typedef Arrangement_     Arrangement;

  typedef typename Arrangement::Inner_ccb_iterator    Inner_ccb_iterator;
  typedef typename Arrangement::Ccb_halfedge_circulator 
                                                  Ccb_halfedge_circulator;
  typedef typename Arrangement::Face_iterator    Face_iterator;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  
  typedef Visitor_         Visitor;

protected:
  Visitor*                     m_visitor;
  std::queue<Inner_ccb_iterator>   m_holes;
  std::stack<Ccb_halfedge_circulator>  m_ccb_stack;

public:

  Gps_bfs_scanner(Visitor& v): m_visitor(&v)
  {}

  void scan(Arrangement& arr)
  {
    Face_iterator ubf;
    for (ubf = arr.faces_begin(); ubf != arr.faces_end(); ++ubf)
    {
      if (ubf->number_of_outer_ccbs() != 0)
        continue;
      if (ubf->visited() == true)
        continue;

      ubf->set_visited(true);
      push_to_queue_holes_of_face(ubf);
      
      while(!m_holes.empty())
      {
        Inner_ccb_iterator hole = m_holes.front();
        m_holes.pop();
        scan(*hole);
      }
    }
  }

  void scan(Ccb_halfedge_circulator ccb)
  {
    _scan(ccb);
    while(!m_ccb_stack.empty())
    {
      Ccb_halfedge_circulator curr_ccb = m_ccb_stack.top();
      m_ccb_stack.pop();
      _scan(curr_ccb);
    }

  }
  void _scan(Ccb_halfedge_circulator ccb)
  {
    Ccb_halfedge_circulator ccb_circ = ccb;
    Ccb_halfedge_circulator ccb_end  = ccb;
	  Face_iterator new_f;
    do
    {
      Halfedge_iterator he  = ccb_circ;
      new_f = he->twin()->face();
      if(!new_f->visited())
      {
        push_to_queue_holes_of_face(he->twin()->face());
        new_f->set_visited(true);
        m_visitor->discovered_face(he->face(), new_f, he); 
        
        //scan(he->twin());
        m_ccb_stack.push(he->twin());
      }
      ++ccb_circ;
    }
    while(ccb_circ != ccb_end);
  }

  
  void push_to_queue_holes_of_face(Face_iterator f)
  {
    for(Inner_ccb_iterator hit = f->inner_ccbs_begin(); 
        hit!= f->inner_ccbs_end(); ++hit)
    {
      m_holes.push(hit);
    }
  }
};

} //namespace CGAL

#endif
