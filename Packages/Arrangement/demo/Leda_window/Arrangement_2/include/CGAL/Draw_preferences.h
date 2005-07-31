// Copyright (c) 2003  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef DRAW_PREFERENCES_H
#define DRAW_PREFERENCES_H

#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Pm_drawer.h>
#include <CGAL/IO/draw_pm.h>

CGAL_BEGIN_NAMESPACE

template< class Arr_2, 
          class Ccb_halfedge_circulator, 
          class Holes_iterator >
class My_Arr_drawer : public Pm_drawer< Arr_2, Window_stream > {
private:
  typedef Pm_drawer<Arr_2, Window_stream> Base;
public:
  My_Arr_drawer(Window_stream & W):
    Pm_drawer<Arr_2, Window_stream>(W), m_W(W) {}
  
  void draw_face(typename Base::Face_handle f) {
    if (f->does_outer_ccb_exist()) {
      Ccb_halfedge_circulator cc = f->outer_ccb();
      do {
	m_W << cc->curve();
      } while (++cc != f->outer_ccb());  
    }

    Holes_iterator hit = f->holes_begin(), eit = f->holes_end();
    for (;hit!=eit; ++hit) {
      Ccb_halfedge_circulator cc=*hit; 
      do {
	m_W << cc->curve();
      } while (++cc != *hit);  
    }      
  }

  void draw_vertices(typename Base::Vertex_const_iterator Vertices_begin, 
		     typename Base::Vertex_const_iterator Vertices_end) {
    m_W << GREEN;
    Base::draw_vertices(Vertices_begin, Vertices_end);
  }
  
  void draw_halfedges(typename Base::Halfedge_const_iterator Halfedges_begin, 
		      typename Base::Halfedge_const_iterator Halfedges_end) {
    m_W << BLUE;
    Base base(window());
    base.draw_halfedges(Halfedges_begin, Halfedges_end);
  }

  Window_stream & m_W;
};
 
CGAL_END_NAMESPACE
#endif
