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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>,
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

#include <CGAL/assertions.h>
#include <CGAL/Sweep_line_2_old/Sweep_line_traits.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

template <class SweepLineTraits_2>
class Point_less_functor 
{
public:
  typedef SweepLineTraits_2           Traits;
  typedef typename Traits::Point_2    Point_2;
  
  Point_less_functor(Traits * t) : m_traits(t)
  {}
  
  bool operator()(const Point_2& p1,const Point_2& p2) const  
  { 
    return (m_traits->compare_xy(p1,p2) == SMALLER);
  }

private:

  /*! a pointer to a traits object */
  Traits * m_traits;
};




template <class SweepLineTraits_2, class Subcurve> 
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  

  Status_line_curve_less_functor(Traits *t) : m_traits(t) 
  {}

  Comparison_result operator()(const Subcurve * c1, const Subcurve * c2) const 
  {
    return m_traits->curve_compare_y_at_x (c1->get_left_end(),
                                           c2->get_last_curve());
  }

  Comparison_result operator()(const Point_2& pt, const Subcurve * c2) const 
  {
    return m_traits->curve_compare_y_at_x (pt,c2->get_last_curve());
  }

     

private:

  /*! a pointer to a traits object */
  Traits * m_traits;
};

CGAL_END_NAMESPACE

#endif
