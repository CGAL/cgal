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
// release       : $CGAL_Revision: CGAL-2.4-I-65 $
// release_date  : $CGAL_Date: 2002/03/19 $
//
// file          : include/CGAL/Pm_insert_utils.h
// package       : Planar_map (5.87)
// maintainer    : Tali Zvi <talizvi@post.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Tali Zvi          <talizvi@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_INSERT_UTILS_H
#define CGAL_PM_INSERT_UTILS_H

#include <CGAL/Handle_for.h>

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

template <class SweepLineTraits_2, 
          class Point_plus_>
class Pm_curve_node : public Handle_for< typename SweepLineTraits_2::X_curve > 
{
  typedef SweepLineTraits_2         Traits;
  typedef typename Traits::X_curve  X_curve; 
  typedef Handle_for<X_curve>       Handle_for_curve;
  typedef Pm_curve_node             Self;

  typedef Point_plus_               Point_plus;

public:

  Pm_curve_node(const X_curve& cv, const Point_plus& p, Traits* traits_) : 
    Handle_for_curve(cv), m_point(p), m_traits(traits_)
    {}

  ~Pm_curve_node() {}

  bool operator==(const Pm_curve_node & cv_node) const {
    return *ptr() == *(cv_node.ptr());
  }
  
  bool operator!=(const Pm_curve_node & cv_node) const {
    return !operator==(cv_node);
  }

  const X_curve & get_curve() const { return  *ptr(); }
  
  Point_plus & get_point() { return m_point; } 

  Self & operator=(const Self & cv_node) {
    Handle_for_curve::operator=(cv_node);
     return *this;
  }  

private:
  Point_plus m_point;
  Traits * m_traits;
};


/*!
    A class describing an intersection point.
    It contains the point itself and a list of
    all curves going through that point, ordered by
    compare_at_x and compare_at_x_right of the traits.
  */
template <class SweepLineTraits_2, 
          class Point_plus_> 
class Pm_point_node
{

public:
  typedef  Pm_point_node                             Self;

  typedef  SweepLineTraits_2                         Traits;
  typedef  typename Traits::X_curve                  X_curve; 
  typedef  typename Traits::Point                    Point;

  typedef  Point_plus_                               Point_plus;

  typedef  Pm_curve_node<SweepLineTraits_2,Point_plus_>
                                                     Curve_node_;
  typedef  std::vector<Curve_node_>                  Curve_node_container;
  typedef  typename Curve_node_container::iterator   Curve_node_iterator;
  typedef  typename Curve_node_container::const_iterator  
                                                     Curve_node_const_iterator;

  Pm_point_node(Traits *traits_) : m_traits(traits_) {}
    
  Pm_point_node(const Curve_node_& cv, Traits *traits_) : 
    m_intersect_p(cv.get_point()), m_traits(traits_) {
    m_curves.push_back(cv);  
  }
    
  Pm_point_node(const Curve_node_& cv, 
                const Point_plus & ref_point, 
                Traits * traits_) : 
    m_intersect_p(ref_point), m_traits(traits_) {
    m_curves.push_back(cv);  
  }

  void add_curve(Curve_node_ & cv) {
    m_curves.push_back(cv);
  }

  Point_plus & get_point() { return m_intersect_p; }
  const Point_plus & get_point() const { return m_intersect_p; }

  Curve_node_iterator curves_begin() { return m_curves.begin(); }
  Curve_node_iterator curves_end() { return m_curves.end(); }
    
  Curve_node_const_iterator curves_begin() const { return m_curves.begin(); }
  Curve_node_const_iterator curves_end() const { return m_curves.end(); }

 protected:
  Point_plus m_intersect_p;
  Curve_node_container m_curves;   
  Traits *m_traits;

};

template <class Point, class SweepLineTraits_2>
class Pm_less_point_xy 
{
 public:
  typedef SweepLineTraits_2  Traits;
    
  Pm_less_point_xy(Traits *traits_) : m_traits(traits_) {}
    
  inline  bool operator()(const Point & p1, const Point & p2) const 
  {
    return (m_traits->compare_xy(p1, p2) == SMALLER);
  }
private:
  Traits * m_traits;
};

CGAL_END_NAMESPACE

#endif
