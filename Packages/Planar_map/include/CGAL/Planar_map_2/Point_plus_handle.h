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

#ifndef CGAL_POINT_PLUS_HANDLE_H
#define CGAL_POINT_PLUS_HANDLE_H


#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

template <class traits, class vertexHandle>
class Point_plus_handle;

// Point_plus_rep:
// Point_plus_rep holds a Point plus a vertex handle of the vertex in the 
// subdivision that will hold that point.
// The reason we need the vertex handle information is to update the 
// subdivision by the time the sweep line progresses without makeing any 
// point location query. This class holds the representation, and the next 
// will hold the Handle to Point_plus.
template <class traits, class vertexHandle>
class Point_plus_rep {
public:
  typedef traits                   Traits;
  typedef typename Traits::Point   Point;
  typedef vertexHandle             Vertex_handle;
  
  Point_plus_rep() {}
  
  Point_plus_rep(const Point & p) : p_(p), v_(Vertex_handle(NULL)) {}
  
  Point_plus_rep(const Point & p, Vertex_handle v) : p_(p), v_(v) {}
  
  ~Point_plus_rep() {}
  
protected:
  friend class Point_plus_handle<Traits, Vertex_handle>;    
  
  Point p_;
  Vertex_handle v_;
};

// Point_plus:
// The handle to Point_plus.
template <class traits, class vertexHandle>
class Point_plus_handle :
  public Handle_for<Point_plus_rep<traits, vertexHandle> > 
{
  typedef Handle_for<Point_plus_rep<traits, vertexHandle> > 
                                                Handle_for_Point_plus_rep;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_3
  using Handle_for_Point_plus_rep::ptr;
#endif

public:
  typedef traits                                Traits;
  typedef typename Traits::Point                Point;
  typedef vertexHandle                          Vertex_handle;
  typedef Point_plus_rep<Traits, Vertex_handle> Point_plus_rep_pm;
  
  Point_plus_handle() : Handle_for_Point_plus_rep() {}
  
  Point_plus_handle(const Point & p) : 
    Handle_for_Point_plus_rep(Point_plus_rep_pm(p)) 
  {  
  }
  
  Point_plus_handle(const Point & p, Vertex_handle v) : 
    Handle_for_Point_plus_rep(Point_plus_rep_pm(p, v)) 
  { 
  }
  
  Point_plus_handle(const Point_plus_handle & p_plus) : 
    Handle_for_Point_plus_rep(p_plus) {}
  
  ~Point_plus_handle() {}
  
  Point_plus_handle & operator=(const Point_plus_handle & p_plus) {
    Handle_for_Point_plus_rep::operator=(p_plus);
    return *this;
  }
  
  bool operator==(const Point_plus_handle & p_plus) const
  { return ptr()->p_ == p_plus.point(); }
  
  bool operator!=(const Point_plus_handle & p_plus) const
  { return !(operator==(p_plus)); }
  
  void set_point(const Point & p) { ptr()->p_ = p; }
  
  void set_vertex (Vertex_handle v) { ptr()->v_ = v; }
  
  const Point & point() const { return ptr()->p_; }
  
  Vertex_handle vertex() const { return ptr()->v_; }
};

CGAL_END_NAMESPACE

#endif
