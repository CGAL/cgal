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
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>

#ifndef CGAL_SWEEP_LINE_FUNCTORS_H
#define CGAL_SWEEP_LINE_FUNCTORS_H

#include <CGAL/assertions.h>

CGAL_BEGIN_NAMESPACE

template <class Point, class SweepLineTraits_2>
class Point_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  
  Point_less_functor(Traits * traits) : m_traits(traits){}
  
  bool operator()(const Point* p1, const Point* p2) const 
  { 
    return (m_traits->compare_xy(*p1,*p2) == SMALLER);
  }

         bool operator()(const Point& p1,const Point& p2)const  
  { 
    return (m_traits->compare_xy(p1,p2) == SMALLER);
  }

private:
  /*! a pointer to a traits class */
  Traits * m_traits;
};




template <class SweepLineTraits_2, class Subcurve> 
class Status_line_curve_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef bool (Status_line_curve_less_functor::*func)
    (const Subcurve*, const Subcurve*) const;

  struct Compare_param {
    Compare_param(Traits *t) : m_compare_func(1), m_traits(t)  {}

    int m_compare_func;
    Traits * m_traits;
  };

  Status_line_curve_less_functor(Compare_param *p ,Point_2 *pt) : m_compare_param(p),
                                                                  m_sweep_pos(pt)
                                                                  //,counter(0)
  {
    m_compare[0] = &Status_line_curve_less_functor::compare_at;
    m_compare[1] = &Status_line_curve_less_functor::compare_right;
  }

  bool operator()(const Subcurve * c1, const Subcurve * c2) const {
    return (this->*m_compare[m_compare_param->m_compare_func])(c1, c2);
  }

  bool compare_at(const Subcurve * c1, const Subcurve * c2)  const 
  {
   /* const Point_2 * p = &(c2->get_last_point());
    if ( m_compare_param->m_traits->compare_x(c1->get_last_point(),  
                                              c2->get_last_point()) == LARGER )
      p = &(c1->get_last_point());*/ //Baruch
    //now 'p' holds the right most point

    Comparison_result r = 
      m_compare_param->m_traits->curves_compare_y_at_x(c1->get_curve(), 
                                       c2->get_curve(), 
                                       *m_sweep_pos);
    if ( r == SMALLER) {
      return true;
    } 
    return false;
  }

  bool compare_right(const Subcurve * c1, const Subcurve * c2)  const 
  {

    //counter++;
    const X_monotone_curve_2 &cv1 = c1->get_curve();
    const X_monotone_curve_2 &cv2 = c2->get_curve();
    Traits *t = m_compare_param->m_traits;

    

      
    /*
    if ( t->curve_is_vertical(cv1) )
    {
      if (t->point_in_x_range(cv2, c1->get_source()) &&
          t->curve_compare_y_at_x(c1->get_right_end(), cv2) == SMALLER )
      {
        return true;
      }
      return false;
    }
    if ( t->curve_is_vertical(cv2))
    {
      if (t->point_in_x_range(cv1, c2->get_source()) &&
          t->curve_compare_y_at_x(c2->get_left_end(), cv1) == LARGER)
      {
        return true;
      }
      return false;
    }*/

    // non of the curves are vertical...
   /* const Point_2 *p = &(c2->get_last_point());
    if ( m_compare_param->m_traits->compare_x(c1->get_last_point(),  
                                              c2->get_last_point()) == LARGER )
      p = &c1->get_last_point();*/

     
    /*Comparison_result r =  t->curves_compare_y_at_x (c1->get_curve(), 
                                                     c2->get_curve(), 
                                                     *p);*/

    Comparison_result r =  t->curves_compare_y_at_x (c1->get_curve(), 
                                                     c2->get_curve(), 
                                                     *m_sweep_pos);

    if (r == EQUAL)
    {
      if(t->curve_is_vertical(cv1))
      {
        if(t->curve_is_vertical(cv2))
          return (c1 < c2);  // both cv1, cvs are vertical
        
        //cv1 is vertical but cv2 is not vertical
        r = t->curve_compare_y_at_x(*m_sweep_pos, cv2);
        if(r == EQUAL)
          return false;
        if(r == LARGER)
          return false;
        else // if we've reached here , r == SMALLER
          return true;
      }
      else 
        if(t->curve_is_vertical(cv2))
        {
          //cv2 is vertical but cv1 is not vertical
          r = t->curve_compare_y_at_x(*m_sweep_pos, cv1);
          if(r == EQUAL)
            return true;
          if(r == LARGER)
            return true;
          else // if we've reached here , r == SMALLER
            return false; 
        }

      // cv1 and cv2 are both not vertical
      r = t->curves_compare_y_at_x_right(c1->get_curve(), 
                                         c2->get_curve(), 
                                         *m_sweep_pos);
    }
    if ( r == SMALLER) {
      return true;
    } 
    if ( r == LARGER ) {
      return false;
    }

    // r = EQUAL (e.g overlap curves)
    return (c1 < c2);
    //return ( c1->getId() < c2->getId() );
  }

  Compare_param * m_compare_param;

private:
  func m_compare[2];

  const Point_2 *m_sweep_pos;

  public:
 // mutable unsigned int counter;
};

CGAL_END_NAMESPACE

#endif
