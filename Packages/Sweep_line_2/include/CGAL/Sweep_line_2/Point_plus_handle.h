// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-65 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Sweep_line_2/Point_plus_handle.h
// package       : arr (1.87)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_POINT_PLUS_HANDLE_H
#define CGAL_POINT_PLUS_HANDLE_H


#ifndef CGAL_HANDLE_FOR_H
#include <CGAL/Handle_for.h>
#endif

#ifndef CGAL_ASSERTIONS_H
#include <CGAL/assertions.h>
#endif


CGAL_BEGIN_NAMESPACE

template <class PM_>
class Point_plus_handle;

// Point_plus_rep:
// Point_plus_rep holds a Point plus a vertex handle of the vertex in the 
// subdivision that will hold that point.
// The reason we need the vertex handle information is to update the 
// subdivision by the time the sweep line progresses without makeing any 
// point location query. This class holds the representation, and the next 
// will hold the Handle to Point_plus.
template <class PM_>
class Point_plus_rep /*: public Ref_counted*/ {
public:
  typedef PM_                          PM;
  typedef typename PM::Traits          Traits;
  typedef typename Traits::Point       Point;
  typedef typename PM::Vertex_handle   Vertex_handle;
  
  Point_plus_rep() {}
  
  Point_plus_rep(const Point& p) : p_(p), v_(Vertex_handle(NULL)) {}
  
  Point_plus_rep(const Point& p, Vertex_handle v) : p_(p), v_(v) {}
  
  ~Point_plus_rep() {}
  
protected:
  friend class Point_plus_handle<PM_>;    
  
  Point p_;
  Vertex_handle v_;
};

// Point_plus:
// The handle to Point_plus.
template <class PM_>
class Point_plus_handle : public  Handle_for<Point_plus_rep<PM_> > {
  typedef Handle_for<Point_plus_rep<PM_> > Handle_for_Point_plus_rep;
public:
  typedef PM_                         PM;
  typedef typename PM::Traits         Traits;
  typedef typename Traits::Point      Point;
  typedef typename PM::Vertex_handle  Vertex_handle;
  typedef Point_plus_rep<PM>          Point_plus_rep_pm;
  
  Point_plus_handle() : Handle() {}
  
  Point_plus_handle(const Point& p) : 
    Handle_for_Point_plus_rep(Point_plus_rep_pm(p)) 
  {  
    //PTR = new Point_plus_rep_pm(p); 
  }
  
  Point_plus_handle(const Point& p, Vertex_handle v) : 
    Handle_for_Point_plus_rep(Point_plus_rep_pm(p,v)) 
  { 
    //PTR = new Point_plus_rep_pm(p,v); 
  }
  
  Point_plus_handle(const Point_plus_handle& p_plus) : 
    Handle_for_Point_plus_rep(p_plus) {}
  
  ~Point_plus_handle() {}
  
  Point_plus_handle& operator=(const Point_plus_handle &p_plus) {
    //Point::operator=(p);
    //ptr()->p_ = p_plus.point();
    //ptr()->v_ = p_plus.vertex();
    Handle_for_Point_plus_rep::operator=(p_plus);
    return *this;
  }
  
  bool operator==(const Point_plus_handle &p_plus) const
  { return ptr()->p_ == p_plus.point(); }
  
  bool operator!=(const Point_plus_handle &p_plus) const
  { return !(operator==(p_plus)); }
  
  void set_point(const Point& p) { ptr()->p_ = p; }
  
  void set_vertex (Vertex_handle v) { ptr()->v_ = v; }
  
  const Point& point() const { return ptr()->p_; }
  
  Vertex_handle vertex() const { return ptr()->v_; } 
  
private:
  //Point_plus_rep_pm* ptr() const { return (Point_plus_rep_pm*) PTR; }
};


CGAL_END_NAMESPACE

#endif





