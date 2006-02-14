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
// Author(s)     : Eti Ezra

#ifndef CGAL_POINT_HANDLE_PLUS_H
#define CGAL_POINT_HANDLE_PLUS_H


#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

template <class PM_>
class Point_handle_plus;

// Point_rep_plus:
// Point_rep_plus holds a Point plus a vertex handle of the vertex in the 
// subdivision that will hold that point.
// The reason we need the vertex handle information is to update the 
// subdivision by the time the sweep line progresses without makeing any 
// point location query. This class holds the representation, and the next 
// will hold the Handle to Point_plus.
template <class PM_>
class Point_rep_plus /*: public Ref_counted*/ {
public:
  typedef PM_                          PM;
  typedef typename PM::Traits          Traits;
  typedef typename Traits::Point       Point;
  typedef typename PM::Vertex_handle   Vertex_handle;
  
  Point_rep_plus() {}
  
  Point_rep_plus(const Point& p) : p_(p), v_(Vertex_handle(NULL)) {}
  
  Point_rep_plus(const Point& p, Vertex_handle v) : p_(p), v_(v) {}
  
  ~Point_rep_plus() {}
  
protected:
  friend class Point_handle_plus<PM_>;    
  
  Point p_;
  Vertex_handle v_;
};

// Point_plus:
// The handle to Point_plus.
template <class PM_>
class Point_handle_plus : public  Handle_for<Point_rep_plus<PM_> > {
  typedef Handle_for<Point_rep_plus<PM_> > Handle_for_Point_rep_plus;
public:
  typedef PM_                         PM;
  typedef typename PM::Traits         Traits;
  typedef typename Traits::Point      Point;
  typedef typename PM::Vertex_handle  Vertex_handle;
  typedef Point_rep_plus<PM>          Point_rep_plus_pm;
  
  Point_handle_plus() : Handle() {}
  
  Point_handle_plus(const Point& p) : 
    Handle_for_Point_rep_plus(Point_rep_plus_pm(p)) 
  {  
    //PTR = new Point_rep_plus_pm(p); 
  }
  
  Point_handle_plus(const Point& p, Vertex_handle v) : 
    Handle_for_Point_rep_plus(Point_rep_plus_pm(p,v)) 
  { 
    //PTR = new Point_rep_plus_pm(p,v); 
  }
  
  Point_handle_plus(const Point_handle_plus& p_plus) : 
    Handle_for_Point_rep_plus(p_plus) {}
  
  ~Point_handle_plus() {}
  
  Point_handle_plus& operator=(const Point_handle_plus &p_plus) {
    //Point::operator=(p);
    //ptr()->p_ = p_plus.point();
    //ptr()->v_ = p_plus.vertex();
    Handle_for_Point_rep_plus::operator=(p_plus);
    return *this;
  }
  
  bool operator==(const Point_handle_plus &p_plus) const
  { return ptr()->p_ == p_plus.point(); }
  
  bool operator!=(const Point_handle_plus &p_plus) const
  { return !(operator==(p_plus)); }
  
  void set_point(const Point& p) { ptr()->p_ = p; }
  
  void set_vertex (Vertex_handle v) { ptr()->v_ = v; }
  
  const Point& point() const { return ptr()->p_; }
  
  Vertex_handle vertex() const { return ptr()->v_; } 
  
private:
  //Point_rep_plus_pm* ptr() const { return (Point_rep_plus_pm*) PTR; }
};


CGAL_END_NAMESPACE

#endif
