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

#include <list>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>

 
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
  typedef SweepLineTraits_2                                Traits;
  typedef typename Traits::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Traits::Point_2                         Point_2;

  typedef CurveWrap                                        SubCurve;
  template<typename SC>
  struct SC_container { typedef std::list<SC> other; };

  typedef std::list<SubCurve*>                             SubcurveContainer; 
  typedef typename SubcurveContainer::iterator             SubCurveIter;
  typedef typename SubcurveContainer::reverse_iterator     SubCurveRevIter;
 
  typedef std::pair<bool, SubCurveIter>                    Pair;


  /*! The type of the event */
  typedef enum 
  {
    DEFAULT = 0,
    
    LEFT_END = 1, // a curve's left-end is on the event point
    
    RIGHT_END = 2, // a curve's right-end is on the event point
    
    ACTION = 4,   // action point 
    
    QUERY = 8,    //query point
    
    INTERSECTION = 16,     // two curves intersects at their interior 
    
    WEAK_INTERSECTION = 32, // when a curve's end-point is on the interior
                           //of another curve (also may indicate overlap)
    OVERLAP = 64 // end-point of an overlap subcurve

  }Attribute;


  Sweep_line_event(){}

  
  void init(const Point_2 &point, Attribute type)
  {
    m_point = point;
    m_type = type;
  }


  /*! Destructor */
  ~Sweep_line_event() 
  {}


  

  
  void add_curve_to_left(SubCurve *curve)
  {
    // look for the curve, and if exists, nothing to do
    for(SubCurveIter iter = m_leftCurves.begin();
        iter != m_leftCurves.end();
        ++iter)
    {
      if((curve == *iter) || (*iter)->is_inner_node(curve))// || (curve)->is_leaf(*iter) || (curve)->has_same_leaf(*iter))
      {
        //*iter = curve;
        return;
      }
      if((curve)->is_inner_node(*iter))
      {
        *iter = curve;
        return;
      }
    }
    m_leftCurves.push_back(curve);
  }



    void push_back_curve_to_left(SubCurve *curve)
    {
      m_leftCurves.push_back(curve);
    }



  std::pair<bool, SubCurveIter>  add_curve_to_right(SubCurve *curve, Traits* tr) 
  {
    if (m_rightCurves.empty()) {
      m_rightCurves.push_back(curve);
      return Pair(false, m_rightCurves.begin());
    }


    SubCurveIter iter = m_rightCurves.begin();
    
    Comparison_result res;
    while ((res = tr->compare_y_at_x_right_2_object()
                    (curve->get_last_curve(),
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
  

  // add two Subcurves to the right of the event.
  //precondition: no right curves, the order of sc1, sc2 is correct
  std::pair<bool, SubCurveIter> add_pair_curves_to_right(SubCurve *sc1,
                                                         SubCurve *sc2)
  {
    m_rightCurves.push_back(sc1);
    m_rightCurves.push_back(sc2);

    SubCurveIter iter = m_rightCurves.end(); --iter;
    return (std::make_pair(false,iter));
  }



  void remove_curve_from_left(SubCurve* curve)
  {
    for(SubCurveIter iter = m_leftCurves.begin();
        iter!= m_leftCurves.end();
        ++iter)
    {
      /*if((*iter)== curve || curve->is_parent(*iter) || (*iter)->is_leaf(curve))
      {
        m_leftCurves.erase(iter);
        return;
      }*/
      if(curve->has_common_leaf(*iter))
      {
         m_leftCurves.erase(iter);
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

  /*! Returns a reverse_iterator to the first curve of the reversed list
      of the right curves of the event */
  SubCurveRevIter right_curves_rbegin()
  {
    return m_rightCurves.rbegin();
  }

  /*! Returns a reverse_iterator to the past-end curve of the reversed list
      of the right curves of the event */
  SubCurveRevIter right_curves_rend()
  {
    return m_rightCurves.rend();
  }

  /*! Returns a reverse_iterator to the first curve of the reversed list
      of the left curves of the event */
  SubCurveRevIter left_curves_rbegin()
  {
    return m_leftCurves.rbegin();
  }

  /*! Returns a reverse_iterator to the past-end curve of the reversed list
      of the left curves of the event */
  SubCurveRevIter left_curves_rend()
  {
    return m_leftCurves.rend();
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
  bool has_left_curves() const{
    return !m_leftCurves.empty();
  }

  /*! Returns true if at least one intersecting curve is defined to 
      the right of the point. */
  bool has_right_curves() const{
    return !m_rightCurves.empty();
  }

  /*! Returns the actual point of the event */
  const Point_2 &get_point() const {
    return m_point;
  }

  /*! Returns the actual point of the event (non-const) */
  Point_2& get_point()
  {
    return m_point;
  }

  /*! change the point of the event. */
  void set_point(const Point_2& pt)
  {
    m_point = pt;
  }

  bool is_left_end() const
  {
    return ((m_type & LEFT_END) != 0);
  }

  bool is_right_end() const
  {
    return ((m_type & RIGHT_END) != 0);
  }

  bool is_intersection() const
  {
    return ((m_type & INTERSECTION ) != 0);
  }

  bool is_action() const
  {
    return ((m_type & ACTION ) != 0);
  }

  bool is_query() const
  {
    return ((m_type & QUERY ) != 0);
  }

  bool is_weak_intersection() const
  {
    return((m_type & WEAK_INTERSECTION) != 0);
  }

  bool is_overlap() const
  {
    return ((m_type & OVERLAP ) != 0);
  }

  void set_left_end()
  {
    m_type |= LEFT_END;
  }

  void set_right_end()
  {
    m_type |= RIGHT_END;
  }

  void set_intersection()
  {
    m_type |= INTERSECTION;
  }

  void set_action()
  {
    m_type |= ACTION;
  }

  void set_query()
  {
    m_type |= QUERY;
  }

  void set_weak_intersection()
  {
    m_type |= WEAK_INTERSECTION;
  }

  void set_overlap()
  {
    m_type |= OVERLAP;
  }

  void set_attribute(Attribute type)
  {
    m_type |= type;
  }

  



  template <class InputIterator>
  void replace_left_curves(InputIterator begin, InputIterator end)
  {
    SubCurveIter left_iter = m_leftCurves.begin();
    for(InputIterator itr = begin; itr != end; ++itr , ++left_iter)
    {
      *left_iter = static_cast<SubCurve*>(*itr);
    }
    m_leftCurves.erase(left_iter, m_leftCurves.end());
  }


  bool is_right_curve_bigger(SubCurve* c1, SubCurve* c2)
  {
    for(SubCurveIter iter = m_rightCurves.begin();
        iter != m_rightCurves.end();
        ++iter)
    {
      if(*iter == c1 ||
         static_cast<SubCurve*>((*iter)->get_orig_subcurve1()) == c1 ||
         static_cast<SubCurve*>((*iter)->get_orig_subcurve2()) == c1)
        return false;
      if(*iter == c2 ||
         static_cast<SubCurve*>((*iter)->get_orig_subcurve1()) == c2 ||
         static_cast<SubCurve*>((*iter)->get_orig_subcurve2()) == c2)
        return true;
    }
    return true;
  }


 
  

#ifdef VERBOSE
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

  /*! */
  char m_type;

};





#ifdef VERBOSE
template<class SweepLineTraits_2, class CurveWrap>
void 
Sweep_line_event<SweepLineTraits_2, CurveWrap>::
Print() 
{
  std::cout << "\tEvent info: "  << "\n" ;
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
