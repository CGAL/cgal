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
#ifndef CGAL_SWEEP_LINE_EVENT_H
#define CGAL_SWEEP_LINE_EVENT_H

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <list>
#include <set>
#include<functional>


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
 * - a list of vertical curves that pass through the event
 * - a list of points that are intersection points on the vertical curves
 * and some more data that is used to help with the algorithm.
 *
 * The class mostly exists to store information and does not have any 
 * significant functionality otherwise.
 * 
 * TODO - implement this class with a set to hold the left and right curves
 * TODO - implement the vertical points array as a set
 */

template<class SweepLineTraits_2, class CurveWrap>
class Sweep_line_event
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;

  //typedef Sweep_line_subcurve<Traits> SubCurve;
  typedef CurveWrap SubCurve;
  typedef std::list<SubCurve *> SubcurveContainer; //TODO - change it to SET (faster?)
  typedef typename SubcurveContainer::iterator SubCurveIter;

  typedef Status_line_curve_less_functor<Traits, SubCurve> StatusLineCurveLess;
  typedef std::set<SubCurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;

	typedef  Point_less_functor<Point_2 ,SweepLineTraits_2 > PointLess;
	typedef std::set<Point_2 , PointLess> VerticalXPointSet;  
  typedef typename VerticalXPointSet::iterator VerticalXPointSetIter; 

  typedef std::list<SubCurve *> VerticalCurveList;
  typedef typename VerticalCurveList::iterator VerticalCurveListIter;

  /*! Constructor */
  Sweep_line_event(const Point_2 &point, Traits *traits) :
    m_point(point), m_traits(traits), m_isInitialized(false),
    m_isInternalIntersectionPoint(false), m_containsOverlap(false)
	

  { 
	  m_verticalCurveXPoints = new VerticalXPointSet(PointLess(m_traits));
    m_leftCurves = new SubcurveContainer();
    m_rightCurves = new SubcurveContainer();
  }

  /*! Destructor. Deletes the lists of curves, without deleting the 
      curves themselves. 
  */
  virtual ~Sweep_line_event() {
    delete m_leftCurves;
    delete m_rightCurves;
  }


  /*! Adds a new curve to the event. The curve is added only to the list/s
   *  in which it is defined (left or/and right).
   *  If the curve is vertical, it is added to the list of vertical curves.
   *
   *  Precondition: The event point has to be either the source or the 
   *  target of the curve.
   *  @param curve  a pointer to the curve.
   */
  void add_curve(SubCurve *scurve)
  {
    const X_monotone_curve_2 &curve = scurve->get_curve();
    const Point_2 &source = m_traits->curve_source(curve);
    const Point_2 &target = m_traits->curve_target(curve);
    
    if ( m_traits->curve_is_vertical(curve) ) 
    {
      m_verticalCurves.push_back(scurve);

    }
    else 
    {
      const Point_2 *rel = &(source);
      if ( m_traits->point_equal(m_point, source) )
	      rel = &(target);
      
      if ( m_traits->compare_x(m_point, *rel) == LARGER )
      {
	      add_curve_to_left(scurve, m_rightmostPointToLeft, true);
      } 
      else 
      {
      	add_curve_to_right(scurve);
      }
    }
  }

  /*! Adds a new curve that is defined to the left of the event point.
   *  The insertion is performed so that the curves remain sorted by their
   *  Y values to the left of the event.                             <br>
   *  If the curve is already in the list of curves, it is removed and 
   *  re-inserted. This way the curves remain sorted.
   *
   *  @param curve  a pointer to the curve.
   *  @pram ref a reference point to perform the compare by
   *  @param isInitStage true when thie method is called at the 
   *  initialization stage (in which case some extra tests are performed).
   *
   * TODO - check to see in which cases the curve is re-inserted with 
   *        a different ordering. Probably in case of conics.
   */
  void add_curve_to_left(SubCurve *curve, const Point_2 &ref, 
			 bool isInitStage=false) 
  {
    if ( isInitStage )
    {
      if ( !m_isInitialized ) 
      {
	if ( curve->is_source_left_to_target()) {
	  m_rightmostPointToLeft = curve->get_source();
	}
	else{
	  m_rightmostPointToLeft = curve->get_target();
	}
	m_isInitialized = true;

      } else {
	update_rightmost_point(curve);
      }

    } else if ( !curve->is_end_point(m_point) ) {
      m_isInternalIntersectionPoint = true;
    }

    // now insert the curve at the right place...

    if (m_leftCurves->empty()) {
      m_leftCurves->push_back(curve);
      return;
    }

    SubCurveIter iter = m_leftCurves->begin();
    const X_monotone_curve_2 &cv = curve->get_curve();
    
    // look for the curve, and if exists, erase it.
    while ( iter != m_leftCurves->end() ) {
      if ( (*iter)->getId() ==  curve->getId()) {
	m_leftCurves->erase(iter);
	break;
      }
      ++iter;
    }
    
    // insert the curve so that the list remains sorted...
    Comparison_result res = SMALLER;
    iter = m_leftCurves->begin();

    while ( iter != m_leftCurves->end() )
    {
      if ( m_traits->point_in_x_range((*iter)->get_curve(), ref))
      {
	const Point_2 &ref_point = largest_point(curve->get_last_point(), 
						 (*iter)->get_last_point());
        res = m_traits->curves_compare_y_at_x (cv, (*iter)->get_curve(), 
					       ref_point);
	if (res == EQUAL) {
	  res = m_traits->curves_compare_y_at_x_right(cv, (*iter)->get_curve(), 
						      ref_point);
	}
      }
      else
      {
	const Point_2 &ref_point = largest_point(curve->get_last_point(),
						 (*iter)->get_last_point());
        res = m_traits->curves_compare_y_at_x (cv, (*iter)->get_curve(), 
					       ref_point);
	if (res == EQUAL)
	  res = m_traits->curves_compare_y_at_x_right(cv, (*iter)->get_curve(), 
						      ref_point);
      }

      if ( res != LARGER )
        break;
      ++iter;
    }
    
    while ( iter != m_leftCurves->end() &&
	    res == EQUAL &&
	    curve->getId() > (*iter)->getId() )
    {
      m_containsOverlap = true;
      ++iter;
      if ( iter == m_leftCurves->end())
	break;

      const Point_2 &ref_point = largest_point(curve->get_last_point(), 
					       (*iter)->get_last_point());
      res = m_traits->curves_compare_y_at_x (cv, (*iter)->get_curve(), 
					     ref_point);
      if (res == EQUAL)
	res = m_traits->curves_compare_y_at_x_right(cv, (*iter)->get_curve(), 
						    ref_point);
    }
    
    // insert the curve. If the curve is already in the list, it is not added
    m_leftCurves->insert(iter, curve);
  }


  /*! Adds a new curve that is defined to the right of the event point.
   *  The insertion is performed so that the curves remain sorted by their Y 
   *  values to the right of the event.
   *  @param curve  a pointer to the curve.
   */
  void add_curve_to_right(SubCurve *curve) 
  {
    if ( !curve->is_end_point(m_point) )
      m_isInternalIntersectionPoint = true;

    if (m_rightCurves->empty()) {
      m_rightCurves->push_back(curve);
      return;
    }


    SubCurveIter iter = m_rightCurves->begin();
    Comparison_result res;
    while (((res = m_traits->curves_compare_y_at_x (curve->get_curve(),
						(*iter)->get_curve(), 
						 m_point)) == LARGER) ||
	   (res == EQUAL &&
	    (res = m_traits->curves_compare_y_at_x_right(curve->get_curve(),
						      (*iter)->get_curve(), 
						      m_point)) == LARGER))
    {
      ++iter;
      if ( iter == m_rightCurves->end()) {
	m_rightCurves->insert(iter, curve);
	return;
      }
    }
    
    while ( res == EQUAL && curve->getId() > (*iter)->getId() )
    {
      m_containsOverlap = true;
      ++iter;
      if ( iter == m_rightCurves->end() ) {
	m_rightCurves->insert(iter, curve);
	return;
      }

      res = m_traits->curves_compare_y_at_x (curve->get_curve(),
					  (*iter)->get_curve(), 
					  m_point);
      if (res == EQUAL)
	res = m_traits->curves_compare_y_at_x_right(curve->get_curve(),
						 (*iter)->get_curve(), 
						 m_point);
    }
    
    // insert the curve only if it is not already in...
    if ( (*iter)->getId() !=  curve->getId()) {
      m_rightCurves->insert(iter, curve);
    }
  }
  

  /*! Returns an iterator to the first curve to the left of the event */
  SubCurveIter left_curves_begin() {
    return m_leftCurves->begin();
  }

  /*! Returns an iterator to the one past the last curve to the left 
      of the event */
  SubCurveIter left_curves_end() {
    return m_leftCurves->end();
  }

  /*! Returns an iterator to the first curve to the right of the event */
  SubCurveIter right_curves_begin() {
    return m_rightCurves->begin();
  }

  /*! Returns an iterator to the one past the last curve to the right 
      of the event */
  SubCurveIter right_curves_end() {
    return m_rightCurves->end();
  }

  /*! Returns the number of intersecting curves that are defined
      to the right of the event point. */
  int get_num_right_curves() {
    return m_rightCurves->size();
  }

  /*! Returns the number of intersecting curves that are defined
      to the left of the event point. */
  int get_num_left_curves() {
    return m_leftCurves->size();
  }

  /*! Returns true if at least one intersecting curve is defined to 
      the left of the point. */
  bool has_left_curves() {
    return !m_leftCurves->empty();
  }

  /*! Returns the actual point of the event */
  const Point_2 &get_point() {
    return m_point;
  }

  /*! 
    @return returns true if at least one of the curves passign 
    through the event is vertical.
  */
  bool does_contain_vertical_curve() const {
    return !m_verticalCurves.empty();
  }


  /*!returns the list of vertical curves passing through the event point.
    @return a reference to the list of curves.
  */
  VerticalCurveList &get_vertical_curves() {
    return m_verticalCurves;
  }


  /*! Insert a new intersection point on any of the vertical curves.
   *  The list of points is sorted by their y values.              <br>
   *  If the requireSort flag is true, the appripriate place in the list 
   *  is searched for. If not, the point is assumed to have the largest y 
   *  value, and is inserted at the end of the list.               <br>
   *  If the pioint already exists, the point is not inserted again.
   *  @param p a reference to the point
   *  @param requireSort false if the point is to be added at the end
   *  of the list.
   *  
   */
  void add_vertical_curve_x_point(const Point_2 &p, bool requireSort=false) 
  {
    m_verticalCurveXPoints->insert(p);
  }
	  /*
    if ( m_verticalCurveXPoints.empty() ) 
    {
      m_verticalCurveXPoints.push_back(p); 
      return;
    }

    if ( !requireSort ) 
    {
      if (!m_traits->point_equal(p, m_verticalCurveXPoints.back())) {
	m_verticalCurveXPoints.push_back(p);
      }
    } else
    {
      VerticalXPointSetIter iter = m_verticalCurveXPoints.begin();
      while ( iter != m_verticalCurveXPoints.end() )
      {
	if ( m_traits->compare_xy(*iter, p) == SMALLER )
	  ++iter; 
	else
	  break;
      }
      if ( iter == m_verticalCurveXPoints.end() )
	m_verticalCurveXPoints.push_back(p);
      else if (!m_traits->point_equal(p, *iter)) {
	m_verticalCurveXPoints.insert(iter, p);
      }
    }
  }
  */

  /*! 
   *  Returns a referece to the list of intersection points on the 
   * vertical curves passign through the event. If no vertical curves 
   * pass through the event or no intersection curves exist, the list 
   * will be empty.
   * @return a reference to the list of points.
   */
  VerticalXPointSet &get_vertical_x_point_list() {
    return *m_verticalCurveXPoints;
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

  /*! 
    @return true if the any two curves in the event overlap, false otherwise.
  */
  bool does_contain_overlap() const {
    return m_containsOverlap;
  }

#ifndef NDEBUG
  void Print();
  void PrintVerticalXPoints();
#endif
 
protected:

  /*! Whenever a new curve is added to the event at the initialization 
   * stage, the right most end point to the left of the event point is 
   * updated.
   * Precondition: the event is either the source or destination of the curve.
   * @param curve a pointer to a new curve added to the event.
   */
  void update_rightmost_point(SubCurve *curve)
  {
    if ( curve->is_source_left_to_target())
    {
      if ( curve->is_target(m_point) )
	if ( m_traits->compare_x(curve->get_source(), 
				 m_rightmostPointToLeft) == LARGER )
	  m_rightmostPointToLeft = curve->get_source();
    } else
    {
      if ( curve->is_source(m_point) )
	if ( m_traits->compare_x(curve->get_target(), 
				 m_rightmostPointToLeft) == LARGER )
	  m_rightmostPointToLeft = curve->get_target();
    }
  }




  


  /*! The point of the event */
  Point_2 m_point;

  /*! A pointer to a traits class */
  Traits *m_traits;

  /*! A list of curves on the left side of the event, sorted by their y value
      to the left of the point */
  SubcurveContainer *m_leftCurves;

  /*! A list of curves on the right side of the event, sorted by their y value
      to the right of the point */
  SubcurveContainer *m_rightCurves;

  /*! The rightmost curve end point that is to the left of the event
      point. This point is used as a reference point when curves are compared
      to the left of the event point. 
  */
  Point_2 m_rightmostPointToLeft;

  /*! An indication whether this event has been initialized. The event is
      initialized after the first curve has been added to the left of the 
      event. 
  */
  bool m_isInitialized;

  /*! a list of vertical curves going through this event */
  VerticalCurveList m_verticalCurves; 

  /*! a list of intersection points on the vertical curves */
  VerticalXPointSet* m_verticalCurveXPoints;

  /*! a flag that inidcates whether the event is an "interior" intersection 
      point, or just an end point of all curves passing through it.
  */
  bool m_isInternalIntersectionPoint;

  /*! true if any two curves passing through the event overlap. */
  bool m_containsOverlap;

  const Point_2 &largest_point(const Point_2 &p1, const Point_2 &p2)
  {
    if ( m_traits->compare_x(p1, p2) == LARGER )
      return p1;
    return p2;
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
  for ( SubCurveIter iter = m_leftCurves->begin() ;
	iter != m_leftCurves->end() ; ++iter )
  {
    std::cout << "\t";
    (*iter)->Print();
    std::cout << "\n";
  }
  std::cout << std::endl;
  std::cout << "\tRight curves: \n" ;
  for ( SubCurveIter iter1 = m_rightCurves->begin() ;
	iter1 != m_rightCurves->end() ; ++iter1 )
  {
    std::cout << "\t";
    (*iter1)->Print();
    std::cout << "\n";
  }
  std::cout <<"\tVertical curves: \n" ;
  for( VerticalCurveListIter iter2 = m_verticalCurves.begin() ; 
       iter2 != m_verticalCurves.end() ; 
       ++iter2)
  {
    std::cout<<"\t";
    (*iter2)->Print();
    std::cout<<"\n";
  }
  std::cout << std::endl;
}

template<class SweepLineTraits_2, class CurveWrap>
void 
Sweep_line_event<SweepLineTraits_2, CurveWrap>::
PrintVerticalXPoints()
{
  std::cout << "Vertical intersection points for " << m_point << ":\n";
  typename std::list<Point_2>::iterator iter = m_verticalCurveXPoints->begin();
  while ( iter != m_verticalCurveXPoints->end() )
  {
    std::cout << "\t" << *iter << "\n";
    ++iter;
  }
}
 
#endif // NDEBUG

CGAL_END_NAMESPACE

#endif // CGAL_SWEEP_LINE_EVENT_H
