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
 * and some more data that is used to help with the algorithm.
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
  typedef std::list<SubCurve *> SubcurveContainer; 
  typedef typename SubcurveContainer::iterator              SubCurveIter;

  typedef Status_line_curve_less_functor<Traits, SubCurve>  StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess>          StatusLine;
  typedef typename StatusLine::iterator                     StatusLineIter;

  typedef std::pair<bool, SubCurveIter>                     Pair;



  Sweep_line_event(){}

  /*! Constructor */
  Sweep_line_event(const Point_2 &point) :
    m_point(point),
    m_isInitialized(false),
    m_isInternalIntersectionPoint(false)
  {}


  void init(const Point_2 &point)
  {
    m_point = point;
    m_isInitialized = false;
    m_isInternalIntersectionPoint = false;
  }


  /*! Destructor */
  ~Sweep_line_event() 
  {}


  

  SubCurveIter add_curve_to_left(SubCurve *curve)
  {
    //insert the curve at the right place...
    if (m_leftCurves.empty())
    {
      m_leftCurves.push_back(curve);
      if(!traits()->curve_is_vertical(curve->get_curve()))
      {
        m_isInitialized = true;
        m_rightmostPointToLeft = curve->get_left_end();
      }
      return m_leftCurves.end();
    }

    //update_rightmost_point(ref) if curve is not vertical
    if(!traits()->curve_is_vertical(curve->get_curve()))
    {
      if(! m_isInitialized)
      {
        m_rightmostPointToLeft = curve->get_last_point();
        m_isInitialized = true;
      }
      else
        if ( traits()->compare_x(curve->get_last_point() , m_rightmostPointToLeft) == LARGER )
          m_rightmostPointToLeft = curve->get_last_point();  
    


      SubCurveIter iter = m_leftCurves.begin();
      const X_monotone_curve_2 &cv = curve->get_curve();
    
      // look for the curve, and if exists, erase it.
      while ( iter != m_leftCurves.end() ) 
      {
        if ( (*iter) ==  curve)
        {
          return ++iter;
        }
        ++iter;
      }
    
      // insert the curve so that the list remains sorted...
      Comparison_result res = SMALLER;
      iter = m_leftCurves.begin();
      while(iter != m_leftCurves.end() &&
            traits()->curve_is_vertical((*iter)->get_curve()))
      {
        ++iter;
      }


      while ( iter != m_leftCurves.end() )
      {
        res = traits()->curves_compare_y_at_x(cv, 
                                             (*iter)->get_curve(),
                                              m_rightmostPointToLeft);
        if (res == EQUAL)
        {
          res = traits()->curves_compare_y_at_x_right(cv, 
                                                      (*iter)->get_curve(),
                                                      m_rightmostPointToLeft);
        }
        if ( res != LARGER )
          break;
        ++iter;
      }
    
      while ( iter != m_leftCurves.end() &&
              res == EQUAL )
    
      {
        if((*iter)->is_parent(curve))
          return ++iter; 
        ++iter;
        if ( iter == m_leftCurves.end())
          break;

        res = traits()->curves_compare_y_at_x(cv, 
                                             (*iter)->get_curve(),
                                              m_rightmostPointToLeft);
      }
    
      // insert the curve. If the curve is already in the list, it is not added
      m_leftCurves.insert(iter, curve);
      return iter;
      }
      else // the curve is vertical  
      {
        SubCurveIter iter = m_leftCurves.begin();

         // look for the curve, and if exists, erase it.
        while ( iter != m_leftCurves.end() ) 
        {
          if ( (*iter) ==  curve)
          {
                  m_leftCurves.erase(iter);
                  break;
          }
          ++iter;
        }
        iter = m_leftCurves.begin();
        while ( iter != m_leftCurves.end() &&
                traits()->curve_is_vertical((*iter)->get_curve()))
        {
          if((*iter)->is_parent(curve))
            return ++iter;
          ++iter;
          if ( iter == m_leftCurves.end())
            break;
        }
        m_leftCurves.insert(iter, curve);
        return iter;
      }
    }

    void add_curve_to_left(SubCurveIter pos, SubCurve *curve)
    {
      m_leftCurves.insert(pos, curve);
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
    while ((res = traits()->curves_compare_y_at_x_right(curve->get_curve(),
                                                      (*iter)->get_curve(), 
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

 
  /*! Mark the event as an intersection point at an interior of a curve.
   */
  void mark_internal_intersection_point() {
    m_isInternalIntersectionPoint = true;
  }

  /*!
    @return returns true if the event is an intersection point at the 
    interior of at least one of the curves passing throuogh the event 
    point.
   */
  bool is_internal_intersection_point() const {
    return m_isInternalIntersectionPoint;
  }

 

#ifndef NDEBUG
  void Print();
#endif
 

  protected:


  /*! The point of the event */
  Point_2 m_point;

  /*! A list of curves on the left side of the event, sorted by their y value
      to the left of the point */
  SubcurveContainer m_leftCurves;

  /*! A list of curves on the right side of the event, sorted by their y value
      to the right of the point */
  SubcurveContainer m_rightCurves;

  /*! The rightmost curve end point that is to the left of the event
      point. This point is used as a reference point when curves are compared
      to the left of the event point. 
  */
  Point_2 m_rightmostPointToLeft;

  /*! A boolean indicating that m_rightmostPointToLeft is initialized */
  bool m_isInitialized;

  
  /*! a flag that inidcates whether the event is an "interior" intersection 
      point, or just an end point of all curves passing through it.
  */
  bool m_isInternalIntersectionPoint;

 


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
