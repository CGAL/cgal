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
// Author(s)     : Tali Zvi          <talizvi@post.tau.ac.il>
#ifndef CGAL_PM_INSERT_UTILS_H
#define CGAL_PM_INSERT_UTILS_H

#include <CGAL/Handle_for.h>
#include <CGAL/Planar_map_2/Point_plus_handle.h>

CGAL_BEGIN_NAMESPACE

template <class SweepLineTraits_2, 
          class Point_plus_>
class Pm_curve_node : public Handle_for< typename SweepLineTraits_2::X_curve > 
{
  typedef SweepLineTraits_2         Traits;
  typedef typename Traits::X_curve  X_curve; 
  typedef Handle_for<X_curve>       Handle_for_curve;
  typedef Pm_curve_node             Self;

  typedef Point_plus_               Point_plus;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Handle_for_curve::ptr;
#endif
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
