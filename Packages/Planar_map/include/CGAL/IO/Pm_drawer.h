// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// Author(s)     : Eti Ezra <estere@post.tau.ac.il>

#ifndef CGAL_IO_PM_DRAWER_H
#define CGAL_IO_PM_DRAWER_H

#ifndef CGAL_PROTECT_IOSTREAM
#include <iostream>
#define CGAL_PROTECT_IOSTREAM
#endif

#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif

CGAL_BEGIN_NAMESPACE


template <class PM_, class Window_>
class Pm_drawer {
public:

  typedef   PM_                                  PM;
  typedef typename PM::Vertex_iterator           Vertex_iterator;
  typedef typename PM::Halfedge_iterator         Halfedge_iterator;
  typedef typename PM::Face_iterator             Face_iterator;
  typedef typename PM::Vertex_const_iterator     Vertex_const_iterator;
  typedef typename PM::Halfedge_const_iterator   Halfedge_const_iterator;
  typedef typename PM::Face_const_iterator       Face_const_iterator;

  typedef typename PM::Vertex_handle             Vertex_handle;
  typedef typename PM::Halfedge_handle           Halfedge_handle;
  typedef typename PM::Face_handle               Face_handle;
  typedef typename PM::Vertex_const_handle       Vertex_const_handle;
  typedef typename PM::Halfedge_const_handle     Halfedge_const_handle;  
  typedef typename PM::Face_const_handle         Face_const_handle;

  typedef Window_   Window;

  Pm_drawer   (Window& w)  { m_window = &w; }
  
  Window&  window()    { return *m_window; }
   
  
  //void draw_vertex_attributes(const Point& p) {}
  
  void draw_vertex(Vertex_handle v) {
    window() << v->point();
  }
  
  void draw_vertex(Vertex_const_handle v) {
    window() << v->point();
  }
  //void draw_halfedge_attributes(const Curve& cv) {}
  
  void draw_halfedge(Halfedge_handle h) {
    window() << h->curve();
  }
  
  void draw_halfedge(Halfedge_const_handle h) {
    window() << h->curve();
  }

  void draw_face(Face_handle f) {}
  
  void draw_face(Face_const_handle f) {}
  
  void draw_vertices(Vertex_iterator Vertices_begin,
                     Vertex_iterator Vertices_end)
  {
    for (Vertex_iterator v_iter = Vertices_begin; v_iter !=  Vertices_end;
         v_iter++)
      draw_vertex(v_iter);
  }

  void draw_vertices(Vertex_const_iterator Vertices_begin,
                     Vertex_const_iterator Vertices_end)
  {
    for (Vertex_const_iterator v_iter = Vertices_begin; v_iter != Vertices_end;
         v_iter++)
      draw_vertex(v_iter);
  }
   
  void draw_halfedges(Halfedge_iterator Halfedges_begin,
                      Halfedge_iterator Halfedges_end)
  {
    for (Halfedge_iterator h_iter = Halfedges_begin; h_iter != Halfedges_end;
         h_iter++)
      draw_halfedge(h_iter);
  }

  void draw_halfedges(Halfedge_const_iterator Halfedges_begin,
                      Halfedge_const_iterator Halfedges_end)
  {
    for (Halfedge_const_iterator h_iter = Halfedges_begin;
         h_iter != Halfedges_end; h_iter++)
      draw_halfedge(h_iter);
  }
  
  void draw_faces(Face_iterator Faces_begin, Face_iterator Faces_end) {
    for (Face_iterator f_iter = Faces_begin; f_iter != Faces_end; f_iter++)
      draw_face(f_iter);
  }

  void draw_faces(Face_const_iterator Faces_begin,
                  Face_const_iterator Faces_end)
  {
    for (Face_const_iterator f_iter = Faces_begin; f_iter != Faces_end;
         f_iter++)
      draw_face(f_iter);
  }

protected:
  Window         *m_window;
};

CGAL_END_NAMESPACE

#endif  // CGAL_IO_PM_DRAWER_H 
