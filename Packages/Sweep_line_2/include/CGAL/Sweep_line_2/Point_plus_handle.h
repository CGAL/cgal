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
