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

#ifndef CGAL_SWEEP_LINE_EVENT_H
#define CGAL_SWEEP_LINE_EVENT_H


#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>
#include <list>
#include <set>

 
CGAL_BEGIN_NAMESPACE

/*! @class Sweep_line_event
 *
 * A class associated with an event in a sweep line algorithm.
 * An intersection point in the sweep line algorithm is refered to as an event.
 * This class contains the information that is associated with any given 
 * event point. This information contains the following:
 * - the actual point 
 * - a list of curves that pass through the event point and defined to 
 *   the left of the event point.
 * - a list of curves that pass through the event point and defined to 
 *   the right of the event point.
 *
 * The class mostly exists to store information and does not have any 
 * significant functionality otherwise.
 * 
 */



template<class SweepLineTraits_2, class CurveWrap>
class Sweep_line_event
{
public:
  typedef SweepLineTraits_2                                 Traits;
  typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;
  typedef typename Traits::Point_2                          Point_2;

  typedef CurveWrap                                         SubCurve;
  typedef std::list<SubCurve*>                              SubcurveContainer; 
  typedef typename SubcurveContainer::iterator              SubCurveIter;

  typedef Status_line_curve_less_functor<Traits, SubCurve>  StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess>          StatusLine;
  typedef typename StatusLine::iterator                     StatusLineIter;

  typedef std::pair<bool, SubCurveIter>                     Pair;



  Sweep_line_event(){}

  
  void init(const Point_2 &point)
  {
    m_point = point;
  }


  /*! Destructor */
  ~Sweep_line_event() 
  {}


  

  
  void add_curve_to_left(SubCurve *curve)
  {
    // look for the curve, and if exists, nothing to do
    if(std::find(m_leftCurves.begin(), m_leftCurves.end(), curve) ==
      m_leftCurves.end())
      m_leftCurves.push_back(curve);
  }



    void push_back_curve_to_left(SubCurve *curve)
    {
      m_leftCurves.push_back(curve);
    }



  std::pair<bool, SubCurveIter>  add_curve_to_right(SubCurve *curve) 
  {
    if (m_rightCurves.empty()) {
      m_rightCurves.push_back(curve);
      return Pair(false, m_rightCurves.begin());
    }


    SubCurveIter iter = m_rightCurves.begin();
    while ( iter != m_rightCurves.end() ) 
    {
      if ( (*iter) ==  curve)
        return Pair(false, m_rightCurves.end());
      ++iter;
    }

    iter = m_rightCurves.begin();
    Comparison_result res;
    while ((res = traits()->curves_compare_y_at_x_right(curve->get_last_curve(),
                                                      (*iter)->get_last_curve(), 
                                                      m_point)) == LARGER)
    {
      ++iter;
      if ( iter == m_rightCurves.end())
      {
              m_rightCurves.insert(iter, curve);
              return Pair(false, --iter);
      }
    }
    
    if ( res == EQUAL ) //overlap !!
    {
      if((*iter)->is_parent(curve))
        return Pair(false,m_rightCurves.end());
      return Pair(true, iter);
     }
     m_rightCurves.insert(iter, curve);
     return Pair(false,--iter);
  }
  



  void remove_curve_from_left(SubCurve* curve)
  {
    for(SubCurveIter iter = m_leftCurves.begin();
        iter!= m_leftCurves.end();
        ++iter)
    {
      if((*iter)== curve)
      {
        m_leftCurves.erase(iter);
        return;
      }
    }
  }




  void replace_right_curve(SubCurve* sc1, SubCurve* sc2)
  {
    for(SubCurveIter iter = m_rightCurves.begin();
        iter!= m_rightCurves.end();
        ++iter)
    {
      if(*iter == sc1)
      {
        *iter = sc2;
        return;
      }
    }
  }


  /*! Returns an iterator to the first curve to the left of the event */
  SubCurveIter left_curves_begin() {
    return m_leftCurves.begin();
  }

  /*! Returns an iterator to the one past the last curve to the left 
      of the event */
  SubCurveIter left_curves_end() {
    return m_leftCurves.end();
  }

  /*! Returns an iterator to the first curve to the right of the event */
  SubCurveIter right_curves_begin() {
    return m_rightCurves.begin();
  }

  /*! Returns an iterator to the one past the last curve to the right 
      of the event */
  SubCurveIter right_curves_end() {
    return m_rightCurves.end();
  }

  /*! Returns the number of intersecting curves that are defined
      to the right of the event point. */
  int get_num_right_curves() {
    return m_rightCurves.size();
  }

  /*! Returns the number of intersecting curves that are defined
      to the left of the event point. */
  int get_num_left_curves() {
    return m_leftCurves.size();
  }

  /*! Returns true if at least one intersecting curve is defined to 
      the left of the point. */
  bool has_left_curves() {
    return !m_leftCurves.empty();
  }

  /*! Returns the actual point of the event */
  const Point_2 &get_point() const {
    return m_point;
  }

 
  

#ifndef NDEBUG
  void Print() ;
#endif
 

  protected:


  /*! The point of the event */
  Point_2 m_point;

  /*! A list of curves on the left side of the event (or traverse at event)*/
  SubcurveContainer m_leftCurves;

  /*! A list of curves on the right side of the event(or traverse at event),
      sorted by their y value to the right of the point */
  SubcurveContainer m_rightCurves;

 
  protected:


  Traits* traits() const
  {
    return Sweep_line_traits<Traits>::get_traits();
  }

  
#ifndef NDEBUG
public:
  int id;
#endif

};





#ifndef NDEBUG
template<class SweepLineTraits_2, class CurveWrap>
void 
Sweep_line_event<SweepLineTraits_2, CurveWrap>::
Print() 
{
  std::cout << "\tEvent id: " << id << "\n" ;
  std::cout << "\t" << m_point << "\n" ;
  std::cout << "\tLeft curves: \n" ;
  for ( SubCurveIter iter = m_leftCurves.begin() ;
        iter != m_leftCurves.end() ; ++iter )
  {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
  std::cout << "\tRight curves: \n" ;
  for ( SubCurveIter iter1 = m_rightCurves.begin() ;
        iter1 != m_rightCurves.end() ; ++iter1 )
  {
    std::cout << "\t";
    (*iter1)->Print();
    std::cout << "\n";
  }
 
  std::cout << std::endl;
}


 
#endif // NDEBUG

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_EVENT_H
