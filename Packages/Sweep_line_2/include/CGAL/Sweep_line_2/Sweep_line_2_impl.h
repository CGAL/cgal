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
#ifndef CGAL_SWEEP_LINE_2_IMPL_H
#define CGAL_SWEEP_LINE_2_IMPL_H

#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/assertions.h>
#include <map>
#include <set>

#ifndef VERBOSE

#define SL_DEBUG(a)
#define PRINT_INSERT(a)
#define PRINT_ERASE(a)
#define PRINT_NEW_EVENT(p, e) 
#define DBG(a)
#define STORE_RESULT(a)

#else

#define SL_DEBUG(a) {a}
#define PRINT_INSERT(a) { std::cout << "+++ inserting "; \
                          (a)->Print(); \
                          std::cout << "    currentPos = " << m_currentPos \
                                    << "\n"; \
                          }
#define PRINT_ERASE(a)  { std::cout << "--- erasing " ; \
                          (a)->Print(); }
#define PRINT_NEW_EVENT(p, e) \
{ std::cout << "%%% a new event was created at " << (p) << std::endl; \
  (e)->Print(); }
#define DBG(a) { std::cout << a << std::endl; }
#define STORE_RESULT(a) {a}
#endif

#include <list>

CGAL_BEGIN_NAMESPACE

/*!
  Sweep_line_2_impl is a class that implements the sweep line algorithm
  based on the algorithm of Bentley and Ottmann.
  It extends the algorithm to support not only segments but polylines and 
  general curves as well.
  The curves are defined by the traits class that is one of the template 
  arguments.

  The algorithm is also extended to support the following degenerate cases:
  - non x-monotone curves
  - vertical segments
  - multiple (more then two) segments intersecting at one point
  - curves beginning and ending on other curves.
  - overlapping curves

  There are two main functionalities supported by this algorithm:
  1. calculate the non intersecting curves that are a product of 
     intersecting a set of input curves
  2. calculate all the intersection points between the curves specified

  An extension of this algorithm is used to produce a planar map containing
  the input curves efficiently. \sa Pmwx_aggregate_insert.

  General flow:
  After the initialization stage, the events are handled from left to right.

  For each event

    First pass - handles special cases in which curves start or end 
                 at the interior of another curve
    Handle vertical curve (bottom) - special caring if the event 
                 is a bottom end of a vertical curve
    Handle vertical overlapping curves - special caring if there are
                 overlapping vertical curves passing through the event.
    Handle left curves - iterate over the curves that intersect 
                 at the event point and defined to the left of the 
                 event. This is mostly where the output (sub curves 
		 or intersection points) is produced.
    Handle vertical curve (top) - special caring if the event id 
                 the upper end of a vertical curve. In case of vertical 
		 curves, this is where the output is produced.
    Handle right curves - iterate over the curves that intersect 
                 the event point and defined to the right of the 
		 event point. This is where new intersection points 
		 are calculated.
  End

  Convensions through out the code:
  In order to make the code as readable as possible, some convensions were 
  made in regards to variable naming:

    xp - is the intersection point between two curves
    slIter - an iterator to the status line, always points to a curve.
    out - always refer to the output iterator

*/

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
class Sweep_line_2_impl
{
public:
  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  typedef SweepEvent Event;
  typedef Point_less_functor<Point_2, Traits> PointLess;
  typedef std::map<const Point_2 *, Event*, PointLess> EventQueue; 
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;
  typedef std::vector<Event*> EventPtrContainer;
  typedef typename EventPtrContainer::iterator EventPtrContainerIter;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef CurveWrap Subcurve;
  typedef std::list<Subcurve*> SubCurveList;
  typedef typename SubCurveList::iterator SubCurveListIter;

  typedef Status_line_curve_less_functor<Traits, Subcurve> StatusLineCurveLess;
  typedef typename StatusLineCurveLess::Compare_param CompareParams;
  typedef std::set<Subcurve*, StatusLineCurveLess> StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;
  typedef typename Event::VerticalCurveList VerticalCurveList;
  typedef typename Event::VerticalCurveListIter VerticalCurveListIter;
  typedef typename Event::VerticalXPointSet VerticalXPointList;
  typedef typename Event::VerticalXPointSetIter VerticalXPointListIter;

  typedef std::list<Event *> EventList;
  typedef typename EventList::iterator EventListIter;

  typedef std::list<X_monotone_curve_2> CurveList;
  typedef typename CurveList::iterator CurveListIter;

  class  SweepLineGetSubCurves {};
  class  SweepLineGetPoints {};
  class  SweepLineGetInterCurveList {};
  class  SweepLinePlanarmap {};

  Sweep_line_2_impl()  : m_traits(new Traits()), m_traitsOwner(true),
    m_includeEndPoints(true), m_found_intersection(false), is_first_point(true)
  {
  }
  Sweep_line_2_impl(Traits *t) : m_traits(t), m_traitsOwner(false),
    m_includeEndPoints(true), m_found_intersection(false), is_first_point(true)
  {
  }

  virtual ~Sweep_line_2_impl();

  /*!
   *  Given a container of curves, this function returns a list of curves
   *  that are created by intersecting the input curves.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param overlapping indicates whether overlapping curves should be 
   *                   reported once or multiple times. If false, the 
   *                   overlapping curves are reported once only.
   */
  template <class OutpoutIterator>
  void  get_subcurves(CurveInputIterator begin, CurveInputIterator end, 
		      OutpoutIterator subcurves, bool overlapping = false)
  { 
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_overlapping = overlapping;
    sweep(subcurves, SweepLineGetSubCurves());
  }

  /*!
   *  Given a range of curves, this function returns a list of points 
   *  that are the intersection points of the curves.
   *  The intersections are calculated using the sweep algorithm.
   *  \param curves_begin the input iterator that points to the first curve 
   *                      in the range.
   *  \param curves_end the input past-the-end iterator of the range.
   *  \param subcurves an iterator to the first curve in the range
   *                   of curves created by intersecting the input curves.
   *  \param endpoints if true, the end points of the curves are reported
   *                   as intersection points. Defaults to true.
   *  \param overlapping indicates whether there are overlapping curves
   *                     in the input range. Defaults to false.
   */
  template <class OutpoutIterator>
  void  get_intersection_points(CurveInputIterator begin, 
                                CurveInputIterator end, 
                                OutpoutIterator points,
                                bool includeEndPoints = true)
  { 
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_includeEndPoints = includeEndPoints;
    m_found_intersection = false;
    sweep(points, SweepLineGetPoints());
  }

 /*!
  *  Given a range of curves, this function returns an iterator 
  *  to the beginning of a range that contains the list of curves 
  *  for each intersection point between any two curves in the 
  *  specified range.
  *  The intersections are calculated using the sweep algorithm.
  *  \param begin the input iterator that points to the first curve 
  *                      in the range.
  *  \param end the input past-the-end iterator of the range.
  *  \param intersecting_curves an iterator to the output
  *  \param endpoints if true, the end points of the curves are reported
  *                   as intersection points. Defaults to true.
  */
  template <class OutputIterator>
  void  get_intersecting_curves(CurveInputIterator begin, 
				CurveInputIterator end, 
				OutputIterator intersecting_curves,
				bool endpoints = true)
  { 
    // TODO - implement...
    SweepLineGetInterCurveList tag;
    (void) tag;
  }

  bool do_curves_intersect(CurveInputIterator begin, 
                            CurveInputIterator end)
  {
    init(begin, end);
    SL_DEBUG(
      PrintSubCurves();
      PrintEventQueue();
    )
    m_includeEndPoints = false;
    std::vector<Point_2> dummy;
    m_stop_at_first_int = true;
    m_found_intersection = false;
    sweep(std::back_inserter(dummy), SweepLineGetPoints(), true);
    return m_found_intersection;
  }


protected:


  void init(CurveInputIterator begin, CurveInputIterator end);
  void init_curve(X_monotone_curve_2 &curve);

  /*! The main loop to calculate intersections among the curves
   *  Looping over the events in the queue, for each event we first
   *  handle the curves that are tothe left of the event point (i.e., 
   *  curves that we are done with), and then we look at the curves 
   *  to the right of the point, which means we attept to find intersections
   *  between them and their neighbours on the sweep line.
   */

  template <class OutpoutIterator, class Op>
  void sweep(OutpoutIterator out, Op tag, bool stop_at_first_int=false)
  {
    EventQueueIter eventIter = m_queue->begin();
    m_prevPos = *eventIter->first;
    m_sweepLinePos = *(eventIter->first);

    while ( eventIter != m_queue->end() )
    {
      const Point_2 *p = (eventIter->first);
      if ( m_traits->compare_x(m_sweepLinePos, *p) == SMALLER ) {
        m_prevPos = m_sweepLinePos;
	m_verticals.clear();
	m_verticalSubCurves.clear();
      }
      m_sweepLinePos = *p;
      m_currentPos = *p;

      p = (eventIter->first);
      m_currentEvent = eventIter->second;
      SL_DEBUG(std::cout << "------------- " << *p << " --------------"
	       << std::endl;
	       PrintStatusLine();
	       m_currentEvent->Print();
      )
      
      if ( m_traits->compare_x(*(eventIter->first), m_sweepLinePos) != EQUAL) {
	SL_DEBUG(std::cout << "====================== clearing miniq " 
		 << eventIter->first  << " "
		 << m_prevPos << "\n";)
	m_miniq.clear();
      }
      m_miniq.push_back(m_currentEvent);
      

      handle_vertical_curve_bottom(tag);
      handle_vertical_overlap_curves();
      handle_left_curves(out, tag);
      
      m_queue->erase(eventIter);
      
      handle_vertical_curve_top(out, tag);
      handle_right_curves(out, tag);
      eventIter = m_queue->begin();
    }

    if ( stop_at_first_int && m_found_intersection )
      return;
  }

  void handle_vertical_curve_bottom(SweepLineGetSubCurves &tag);
  void handle_vertical_curve_bottom(SweepLineGetPoints &tag);
  void handle_vertical_overlap_curves();



  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void handle_left_curves(OutpoutIterator out, 
			SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)

    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    m_currentPos = m_prevPos;
    const Point_2 &eventPoint = m_currentEvent->get_point();

    m_use_hint_for_erase = false;
    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    // ** fix here
    {
      Subcurve *leftCurve = *leftCurveIter; 
      const X_monotone_curve_2 &cv = leftCurve->get_curve();
      const Point_2 &lastPoint = leftCurve->get_last_point();

      if ( leftCurve->is_source(eventPoint))
      {      
        if ( !leftCurve->is_target(lastPoint) )
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          add_curve_to_output(a, leftCurve, out);
        } else {
          add_curve_to_output(cv, leftCurve, out);
        }
      } else if ( leftCurve->is_target(eventPoint))
      {
        if ( !leftCurve->is_source(lastPoint))
        {
          X_monotone_curve_2 a,b;
          m_traits->curve_split(cv, a, b, lastPoint);
          add_curve_to_output(b, leftCurve, out);
        } else {
          add_curve_to_output(cv, leftCurve, out);
        }

      } else { 
        X_monotone_curve_2 a,b;
        if ( leftCurve->is_source(lastPoint)) {
          m_traits->curve_split(cv, a, b, eventPoint);
          add_curve_to_output(a, leftCurve, out);
        } else if ( leftCurve->is_target(lastPoint)) {
          m_traits->curve_split(cv, b, a, eventPoint);
          add_curve_to_output(a, leftCurve, out);
        } else {
          const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
          if ( leftCurve->is_source_left_to_target() ) {
            m_traits->curve_split(lastCurve, a, b, eventPoint);
            add_curve_to_output(a, leftCurve, out);
          } else {
	    m_traits->curve_split(lastCurve, b, a, eventPoint);
	    add_curve_to_output(a, leftCurve, out);
          }
        }
        leftCurve->set_last_point(eventPoint);
        leftCurve->set_last_curve(b); 
      
      }
      
      // remove curve from the status line (also checks intersection 
      // between the neighbouring curves)
      remove_curve_from_status_line(leftCurve);
      m_use_hint_for_erase = true;

      m_currentPos = m_prevPos;
      ++leftCurveIter;
    }
    SL_DEBUG(std::cout << "Handling left curve END" << std::endl;)
      
  }

  /*!
   *  Handle a vertical curve when the event being processed is the top end 
   *  of the curve. In this situation, the event contains a list of
   *  intersection points on the vertical curve. We go through this list and
   *  outpt the subcurves induced by these intersection points.
   *  If the curve is not vertical, returns without doing anything.
   * 
   *  @param out an iterator to the output
   *  @param tag a tag that indicates the version of the method
   */
  template <class OutpoutIterator>
  void handle_vertical_curve_top(OutpoutIterator out, 
				 SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout<<"handle_vertical_curve_top... (" 
	              << m_currentEvent->get_point() << ")\n";)
    if ( !m_currentEvent->does_contain_vertical_curve() ) {
      SL_DEBUG(std::cout<<"exiting\n ";)
      return;
    }
    SL_DEBUG(std::cout<<"\n ";)

    VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
    VerticalCurveListIter vciter = vcurves.begin();

    while ( vciter !=vcurves.end() )
    {

      Subcurve *vcurve = *vciter;
      const Point_2 &topPoint = m_currentEvent->get_point();
      // if this is the bottom point, nothing to do here
      if ( vcurve->is_bottom_end(topPoint)) {
	SL_DEBUG(std::cout<<"this is the bottom. skipping.\n";)
	++vciter;
	continue;
      }

      SL_DEBUG(std::cout<<"handling top point of vertical curve\n";)


      // the following while loop comes to handle  | 
      // in the case where a new curve begins at   |------
      // a vertical curve                          |

      // find the "position" of the curve of the status line
      StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
      
      if ( slIter != m_statusLine->end() ) 
      {
	SL_DEBUG(std::cout<<"starting at curve \n";)
	SL_DEBUG((*slIter)->Print();)

	while (slIter != m_statusLine->end() &&
	       m_traits->point_in_x_range((*slIter)->get_curve(), topPoint) &&
	       m_traits->curve_compare_y_at_x(
		 topPoint, 
		 (*slIter)->get_curve()) == LARGER &&
	       m_traits->point_in_x_range((*slIter)->get_curve(), 
					  vcurve->get_bottom_end()) &&
	       m_traits->curve_compare_y_at_x(
		 vcurve->get_bottom_end(), 
		 (*slIter)->get_curve()) == SMALLER)
	{
	  SL_DEBUG(std::cout<<"checking \n";)
	  SL_DEBUG((*slIter)->Print();) 
	  if ( m_traits->compare_x(
		 (*slIter)->get_left_end(), topPoint) == EQUAL)
          {
	    m_currentEvent->add_vertical_curve_x_point(
	      (*slIter)->get_left_end(), true);
	  }
	  ++slIter;
	}   
      }

      // now we go over the list of intersection points on the vertical
      // curve in at the event and process them...
      SL_DEBUG(std::cout<<"handling the splitting now\n";)
      VerticalXPointList &pointList = 
	                    m_currentEvent->get_vertical_x_point_list();
      if ( pointList.empty() )
      {
	add_vertical_curve_to_output(out, vcurve->get_curve());
	++vciter;
	continue;
      }
    
      X_monotone_curve_2 a, b, c;
      a = vcurve->get_curve();
      SL_DEBUG(std::cout << "there are " << pointList.size() << " points\n";)
      SL_DEBUG(m_currentEvent->PrintVerticalXPoints();)
      for ( VerticalXPointListIter i = pointList.begin() ;
	    i != pointList.end(); ++i )
      {
	SL_DEBUG(std::cout<< "splitting: " << a << " at " << *i ;)
	if ( !vcurve->is_point_in_range(*i) )
	{
	  SL_DEBUG(std::cout << " not !\n";)
	  continue;
	}
	SL_DEBUG(std::cout << " yes! \n";)
	m_traits->curve_split(a, b, c, *i);
	if ( vcurve->is_source_left_to_target()) {
	  add_vertical_curve_to_output(out, b);
	  a = c;
	} else {
	  add_vertical_curve_to_output(out, c);
	  a = b;
	}
      }
      if ( vcurve->is_source_left_to_target() ) {
	add_vertical_curve_to_output(out, c);
      }
      else {
	add_vertical_curve_to_output(out, b);
      }
      ++vciter;
    }
  }


  /*! Loop over the curves to the right of the status line and handle them:
   * - if we are at the beginning of the curve, we insert it to the status 
   *   line, then we look if it intersects any of its neighbours.
   * - if we are at an intersection point between two curves, we add them
   *   to the status line and attempt to intersect them with their neighbours. 
   * - We also check to see if the two intersect again to the right of the 
   *   point.
   */
  template <class OutpoutIterator>
  void handle_right_curves(OutpoutIterator out, SweepLineGetSubCurves &tag)
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;)
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";)
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 )
      return;
      
    m_currentPos = m_sweepLinePos;
    if ( numRightCurves == 1 )
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;)
	  
      SL_DEBUG(
	    Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
	    PRINT_INSERT(tmp1);
	    )
	  
      StatusLineIter slIter = m_statusLine->insert(
	                        m_status_line_insert_hint, 
				*(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); //xxx
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;

      SL_DEBUG(PrintStatusLine();)

      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
	return;

      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;

      SubCurveList mylist;
      if ( slIter != m_statusLine->begin() )
      {
	--prev;
	
	StatusLineIter tmp = prev;
	mylist.push_back(*prev);
	while ( tmp != m_statusLine->begin() ) 
	{
	  --tmp;
	  if ( do_curves_overlap(*prev, *tmp) )
	    mylist.push_back(*tmp);
	  else 
	    break;
	}
      }
      
      if ( next != m_statusLine->end() )
      { 
	
	StatusLineIter tmp = next;
	mylist.push_back(*next);
	++tmp;
	while ( tmp != m_statusLine->end() ) 
	{
	  if ( do_curves_overlap(*next, *tmp) )
	  {
	    mylist.push_back(*tmp);
	    ++tmp;
	  }
	  else
	    break;
	}
      }
      intersect_curve_group(*(m_currentEvent->right_curves_begin()), 
			    mylist, out);
      
    } else
    {
      
      /* this block takes care of 
      //
      //           /
      //          /
      //       --------
      //          \
      //           \
      */
      int numLeftCurves = m_currentEvent->get_num_left_curves();
      if ( numLeftCurves == 0 ) {
	
	SL_DEBUG(std::cout << " - handling special case " << std::endl;)
	StatusLineIter slIter;
						     
	EventCurveIter currentOne = m_currentEvent->right_curves_begin();
	while ( currentOne != m_currentEvent->right_curves_end() ) {
	  slIter = m_statusLine->lower_bound(*currentOne);
	  if ( slIter != m_statusLine->end() ) {
	    Subcurve *c = *slIter;
	    if ( CurveStartsAtCurve(*currentOne, c)) {
	      m_currentEvent->add_curve_to_left(c, m_sweepLinePos);
	      m_currentEvent->add_curve_to_right(c);
	      X_monotone_curve_2 a,b;
	      if ( c->is_source_left_to_target() ) {
		m_traits->curve_split(c->get_last_curve(), a, b, 
				      m_currentEvent->get_point());
	      } else {
		m_traits->curve_split(c->get_last_curve(), b, a, 
				      m_currentEvent->get_point());
	      }
	      c->set_last_point(m_currentEvent->get_point());
	      c->set_last_curve(b);
	      
	      add_curve_to_output(a, c, out);
	      break;
	    }
	  }
	  currentOne++;
	}
      }
      // end block ...
      
      
      
      SubCurveList mylist;
      SubCurveList prevlist;
      SubCurveList currentlist;
      
      SL_DEBUG(std::cout << " - intersection point " << std::endl;);
      EventCurveIter firstOne = m_currentEvent->right_curves_begin();
      EventCurveIter lastOne = m_currentEvent->right_curves_end(); --lastOne;
      EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();
      
      PRINT_INSERT(*firstOne);

      StatusLineIter slIter = m_statusLine->insert(m_status_line_insert_hint, 
						   *firstOne);
      (*firstOne)->set_hint(slIter);
      
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine->begin() )
      { 
	StatusLineIter prev = slIter; --prev;
	
	// find all curves that are overlapping with the prev curve
	StatusLineIter tmp = prev;
	prevlist.push_back(*prev);
	while ( tmp != m_statusLine->begin() ) 
	{
	  --tmp;
	  if ( do_curves_overlap(*prev, *tmp))
	    prevlist.push_back(*tmp);
	  else
	    break;
	}
	
	intersect_curve_group(*slIter, prevlist, out);
      }
      currentlist.push_back(*firstOne);
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      
      while ( currentOne != rightCurveEnd )
      {
	m_currentPos = m_sweepLinePos;
	PRINT_INSERT(*currentOne);
	++slIter;
	slIter = m_statusLine->insert(slIter, *currentOne);
	(*currentOne)->set_hint(slIter);
	
	SL_DEBUG(PrintStatusLine(););
	if ( do_curves_overlap(*currentOne, *prevOne))
	{
	  intersect_curve_group(*currentOne, currentlist, out);
	  currentlist.push_back(*currentOne);
	} else {
	  prevlist = currentlist;
	  currentlist.clear();
	  currentlist.push_back(*currentOne);
	}
	
	intersect_curve_group(*currentOne, prevlist, out);
	prevOne = currentOne;
	++currentOne;
      }
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      lastOne = currentOne; --lastOne;
      m_currentPos = m_sweepLinePos;
      //PRINT_INSERT(*lastOne);
      
      SL_DEBUG(PrintStatusLine();)
	StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() ) {
	intersect_curve_group(*next, currentlist, out, true);
	StatusLineIter tmp = next; ++tmp;
	while ( tmp != m_statusLine->end() ) 
	{
	  if ( do_curves_overlap(*next, *tmp))
	  {
	    intersect_curve_group(*tmp, currentlist, out, true);
	    ++tmp;
	  }
	  else
	    break;
	}
      }
    }
  }
  

  template <class OutpoutIterator>
  void handle_right_curves(OutpoutIterator out, SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;);
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";);
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if ( numRightCurves == 0 )
      return;
    
    m_currentPos = m_sweepLinePos;
    if ( numRightCurves == 1 )
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);

      SL_DEBUG(
	Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
	PRINT_INSERT(tmp1);
	);

      StatusLineIter slIter = m_statusLine->insert(
                                      m_status_line_insert_hint, 
				      *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); 
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      SL_DEBUG(PrintStatusLine(););
      
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine->size() == 1 )
	return;
      
      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;

      SubCurveList mylist;
      if ( slIter != m_statusLine->begin() )
      {
	--prev;
	
	if ( CurveStartsAtCurve(*slIter, *prev) && !m_includeEndPoints) {
	  SL_DEBUG(std::cout << "Reporting point (4): " 
                             << (*slIter)->get_left_end() 
                 	     << "\n";);
	  //*out = (*slIter)->get_left_end(); ++out;
	  add_point_to_output((*slIter)->get_left_end(), out);
	  m_found_intersection = true;
	}
	StatusLineIter tmp = prev;
	mylist.push_back(*prev);
	while ( tmp != m_statusLine->begin() ) 
	{
	  --tmp;
	  if ( do_curves_overlap(*prev, *tmp) )
	    mylist.push_back(*tmp);
	  else 
	    break;
	}
      }
      
      if ( next != m_statusLine->end() )
      { 
	if ( CurveStartsAtCurve(*slIter, *next)  && !m_includeEndPoints) {
	  SL_DEBUG(std::cout << "Reporting point (5): " 
		   << (*slIter)->get_left_end() 
		   << "\n";);
	  //*out = (*slIter)->get_left_end(); ++out;
	  add_point_to_output((*slIter)->get_left_end(), out);
	  m_found_intersection = true;
	}
	StatusLineIter tmp = next;
	mylist.push_back(*next);
	++tmp;
	while ( tmp != m_statusLine->end() ) 
	{
	  if ( do_curves_overlap(*next, *tmp) )
	  {
	    mylist.push_back(*tmp);
	    ++tmp;
	  }
	  else
	    break;
	}
      }
      intersect_curve_group(*(m_currentEvent->right_curves_begin()), mylist);
      
    } else
    {
      /* this block takes care of 
      //
      //           /
      //          /
      //       --------
      //          \
      //           \
      */
      int numLeftCurves = m_currentEvent->get_num_left_curves();
      if ( numLeftCurves == 0 ) {
	
	SL_DEBUG(std::cout << " - handling special case " << std::endl;);
	StatusLineIter slIter;
	EventCurveIter currentOne = m_currentEvent->right_curves_begin();
	while ( currentOne != m_currentEvent->right_curves_end() ) {
	  slIter = m_statusLine->lower_bound(*currentOne);
	  if ( slIter != m_statusLine->end() ) {
	    Subcurve *c = *slIter;
	    if ( CurveStartsAtCurve(*currentOne, c) && !m_includeEndPoints) {
	      SL_DEBUG(std::cout << "Reporting point (8): "
                                 << (*currentOne)->get_left_end()
		                 << "\n";);
	      //*out = (*currentOne)->get_left_end(); ++out;
	      add_point_to_output((*currentOne)->get_left_end(), out);
	      m_found_intersection = true;
	      break;
	    }
	  }
	  currentOne++;
	}
      }
      // end block ..
      
      
      SubCurveList mylist;
      SubCurveList prevlist;
      SubCurveList currentlist;
      
      
      SL_DEBUG(std::cout << " - intersection point " << std::endl;);
      EventCurveIter firstOne = m_currentEvent->right_curves_begin();
      EventCurveIter lastOne = m_currentEvent->right_curves_end(); --lastOne;
      EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();
      
      PRINT_INSERT(*firstOne);
      
      StatusLineIter slIter = m_statusLine->insert(m_status_line_insert_hint, 
						   *firstOne);
      (*firstOne)->set_hint(slIter);
      
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine->begin() )
      { 
	StatusLineIter prev = slIter; --prev;
	
	if ( CurveStartsAtCurve(*slIter, *prev) && !m_includeEndPoints) {
	  SL_DEBUG(std::cout << "Reporting point (6): " 
		   << (*slIter)->get_left_end() 
		   << "\n";);
	  //*out = (*slIter)->get_left_end(); ++out;
	  add_point_to_output((*slIter)->get_left_end(), out);
	  m_found_intersection = true;
	}
	
	// find all curves that are overlapping with the prev curve
	StatusLineIter tmp = prev;
	prevlist.push_back(*prev);
	while ( tmp != m_statusLine->begin() ) 
	{
	  --tmp;
	  if ( do_curves_overlap(*prev, *tmp))
	    prevlist.push_back(*tmp);
	  else
	    break;
	}
	
	intersect_curve_group(*slIter, prevlist);
      }
      currentlist.push_back(*firstOne);
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      
      while ( currentOne != rightCurveEnd )
      {
	m_currentPos = m_sweepLinePos;
	PRINT_INSERT(*currentOne);
	++slIter;
	slIter = m_statusLine->insert(slIter, *currentOne);
	(*currentOne)->set_hint(slIter);
	
	SL_DEBUG(PrintStatusLine(););
	if ( do_curves_overlap(*currentOne, *prevOne))
	{
	  intersect_curve_group(*currentOne, currentlist);
	  currentlist.push_back(*currentOne);
	} else {
	  prevlist = currentlist;
	  currentlist.clear();
	  currentlist.push_back(*currentOne);
	}
	
	intersect_curve_group(*currentOne, prevlist);
	prevOne = currentOne;
	++currentOne;
      }
      m_status_line_insert_hint = slIter; ++m_status_line_insert_hint;
      
      lastOne = currentOne; --lastOne;
      m_currentPos = m_sweepLinePos;
      PRINT_INSERT(*lastOne);
      
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine->end() ) {
	
	if ( CurveStartsAtCurve(*slIter, *next)  && !m_includeEndPoints) {
	  SL_DEBUG(std::cout << "Reporting point (7): " 
                             << (*slIter)->get_left_end() 
  		             << "\n";);
	  //*out = (*slIter)->get_left_end(); ++out;
	  add_point_to_output((*slIter)->get_left_end(), out);
	  m_found_intersection = true;
	}
	
	intersect_curve_group(*next, currentlist);
	StatusLineIter tmp = next; ++tmp;
	while ( tmp != m_statusLine->end() ) 
	{
	  if ( do_curves_overlap(*next, *tmp))
	  {
	    intersect_curve_group(*tmp, currentlist);
	    ++tmp;
	  }
	  else
	    break;
	}
      }
    }
  }
  


  // utility methods 
  bool intersect(Subcurve *c1, Subcurve *c2);
  void intersect_curve_group(Subcurve *c1, SubCurveList &mylist);


  // reverse = false ==> we check if the curve starts at the list of curves
  // reverse = true ==> we check if any of the curves in the list start at 
  // the single curve
  template <class OutpoutIterator>
  void intersect_curve_group(Subcurve *c1, SubCurveList &mylist,
			   OutpoutIterator out, bool reverse=false)
  {
    m_tmpOut.clear();
    SL_DEBUG(std::cout << "intersect_curve_group (with out)\n";)
    SL_DEBUG(std::cout << "intersecting with " << mylist.size()
             << " curves\n";)
    SubCurveListIter i = mylist.begin();
    while ( i != mylist.end())
    {
      bool flag;
      if ( reverse ) {
	flag = CurveStartsAtCurve(*i, c1);
	if ( flag && (c1->get_last_point() != m_currentEvent->get_point()) ) {
	  SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";)
	  m_currentEvent->add_curve_to_right(c1);
	  m_currentEvent->add_curve_to_left(c1, m_prevPos);
	  X_monotone_curve_2 a,b;
	  SL_DEBUG(std::cout << "splitting " << (c1)->get_last_curve() 
                             << " at " 
		             << m_currentEvent->get_point() << "\n";)
	  if ( (c1)->is_source_left_to_target() ) 
	    m_traits->curve_split((c1)->get_last_curve(), a, b, 
				  m_currentEvent->get_point());
	  else
	    m_traits->curve_split((c1)->get_last_curve(), b, a, 
				  m_currentEvent->get_point());
	  (c1)->set_last_point(m_currentEvent->get_point());
	  (c1)->set_last_curve(b); 
	  (c1)->set_last_subcurve(a); 
	  m_tmpOut.push_back(c1);
	}
      }
      else
      {
	flag = CurveStartsAtCurve(c1, *i);
	if ( flag && ((*i)->get_last_point() != m_currentEvent->get_point())) {
	  
	  SL_DEBUG(std::cout << "CurveStartsAtCurve 3 \n";)
	  m_currentEvent->add_curve_to_right(*i);
	  m_currentEvent->add_curve_to_left(*i, m_prevPos);
	  X_monotone_curve_2 a,b;
	  SL_DEBUG(std::cout << "splitting " << (*i)->get_last_curve() 
                             << " at " 
		             << m_currentEvent->get_point() << "\n";)
	  if ( (*i)->is_source_left_to_target() ) 
	    m_traits->curve_split((*i)->get_last_curve(), a, b, 
				  m_currentEvent->get_point());
	  else
	    m_traits->curve_split((*i)->get_last_curve(), b, a, 
				  m_currentEvent->get_point());
	  (*i)->set_last_point(m_currentEvent->get_point());
	  (*i)->set_last_curve(b); 
	  (*i)->set_last_subcurve(a); 
	  m_tmpOut.push_back(*i);
	  
	}
      }
      
      intersect(c1, *i);
      ++i;
    }
    
    for ( SubCurveListIter iter = m_tmpOut.begin() ; iter != m_tmpOut.end();
	  ++iter) {
      add_curve_to_output((*iter)->get_last_subcurve(), *iter, out);
    }
  }

  void remove_curve_from_status_line(Subcurve *leftCurve);

  bool is_internal_x_point(const Point_2 &p);
  bool handle_vertical_curve_x_at_end(Subcurve *vcurve, Subcurve *curve, 
				      Event *topEndEvent,
				      SweepLineGetSubCurves tag);
  bool handle_vertical_curve_x_at_end(Subcurve *vcurve, Subcurve *curve, 
				      Event *topEndEvent, 
				      SweepLineGetPoints tag);


  bool do_curves_overlap(Subcurve *c1, Subcurve *c2);
  bool similar_curves(const X_monotone_curve_2 &a, 
		      const X_monotone_curve_2 &b);
  bool vertical_subcurve_exists(const X_monotone_curve_2 &a);


  /*! Adds a new curve to the output list. If the overlapping flag is false,
      each unique curve is reported only once.

      @param out an iterator to the output container
      @param cv the curve to be added
  */
  template <class OutpoutIterator>
  void add_curve_to_output(const X_monotone_curve_2 &cv, Subcurve *curve, 
                           OutpoutIterator out)
  {
    static Subcurve *prevCurve = 0;
    static X_monotone_curve_2 prevXCv;

    if ( m_overlapping ) {
      *out = cv;
      ++out;
    } else {
      if ( prevCurve && similar_curves(cv, prevXCv)) {
	SL_DEBUG(std::cout << " curve already reported... " << std::endl;)
	return;
      }
      prevCurve = curve;
      prevXCv = cv;
      *out = cv;
      ++out;
    }
  }
  template <class OutpoutIterator>
  void add_point_to_output(const Point_2 p, OutpoutIterator out) {

    SL_DEBUG(std::cout << "Reporting point " << p;)
    if ( is_first_point ) {
      is_first_point = false;
      *out = p;
      m_lastReportedPoint = p;
      ++out;
      SL_DEBUG(std::cout << " YES - first time\n";)
    } else {
      if ( m_lastReportedPoint != p ) {
	*out = p;
	++out;
	m_lastReportedPoint = p;
	SL_DEBUG(std::cout << " YES\n";)
      } else {
	SL_DEBUG(std::cout << " NO\n";)
      }
    }
  }


  template <class OutpoutIterator>
  void add_vertical_curve_to_output(OutpoutIterator out, 
				const X_monotone_curve_2 &cv)
  {
    if ( m_overlapping ) {
      *out = cv;
      ++out;
    } else {
      if ( vertical_subcurve_exists(cv)) {
	SL_DEBUG(std::cout << " curve already reported... " << std::endl;)
	return;
      }
      m_verticalSubCurves.push_back(cv);
      *out = cv;
      ++out;
    }
  }

  /*! 
   * Returns true if the point is in the interior of the curve.
   * Returns false if the point is outside the range of the curve or
   * if the point is either the source or the target of the curve.
   * @return true if the point is int he interior of the curve.
   */
  bool is_point_in_curve_interior(const X_monotone_curve_2 &c, 
				  const Point_2 &p)
  {
    if (! m_traits->point_in_x_range(c,p) || 
	m_traits->curve_compare_y_at_x(p, c) != EQUAL)
      return false;
    if ( is_end_point(p) )
      return false;
    return true;
  }


  //
  // a second implementation for the get_intersection_points functionality
  //

  /*!
   *  Handle a vertical curve when the event being processed is the top end 
   *  of the curve. In this situation, the event contains a list of
   *  intersection points on the vertical curve. We go through this list and
   *  output the intersection points.
   *  If the curve is not vertical, returns without doing anything.
   * 
   *  @param out an iterator to the output
   *  @param tag a tag that indicates the version of the method
   */
  template <class OutpoutIterator>
  void handle_vertical_curve_top(OutpoutIterator out, 
			      SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout<<"handle_vertical_curve_top... ";)
    if ( !m_currentEvent->does_contain_vertical_curve() )
    {
      SL_DEBUG(std::cout<<"exiting\n ";)
      return;
    }
    SL_DEBUG(std::cout<<"\n ";)

    VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
    VerticalCurveListIter vciter = vcurves.begin();

    while ( vciter != vcurves.end() )
    {
      Subcurve *vcurve = *vciter; 
      const Point_2 &topPoint = m_currentEvent->get_point();
      if ( vcurve->is_bottom_end(topPoint)) {
	SL_DEBUG(std::cout<<"this is the bottom. skipping.\n";)
	++vciter;
	continue;
      }

      SL_DEBUG(std::cout<<"handling top point of vertical curve\n";)
      StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
      if ( slIter != m_statusLine->end() )
      {
	SL_DEBUG(std::cout<<"starting at curve \n";)
	SL_DEBUG((*slIter)->Print();)
	  
        const Point_2 &bottomPoint = vcurve->get_bottom_end();

	while (slIter != m_statusLine->end() &&
	       m_traits->point_in_x_range((*slIter)->get_curve(),
                                          topPoint) &&
	       m_traits->curve_compare_y_at_x(topPoint,
                                              (*slIter)->get_curve()) ==
               LARGER &&
	       m_traits->point_in_x_range((*slIter)->get_curve(), 
                                          bottomPoint) && 
	       m_traits->curve_compare_y_at_x(bottomPoint,
                                              (*slIter)->get_curve()) ==
               SMALLER)
	{
	  SL_DEBUG(std::cout<<"checking \n";)
	  SL_DEBUG((*slIter)->Print();) 
	  if ( m_traits->compare_x((*slIter)->get_left_end(),topPoint) 
	                                                       == EQUAL)
	  {
	    m_currentEvent->add_vertical_curve_x_point(
	                                    (*slIter)->get_left_end());
	    // if this point was not an event, we need to report the point
	    // test40/42
	    if ( !m_includeEndPoints && 
		 !is_internal_x_point((*slIter)->get_left_end())) {
	      SL_DEBUG(std::cout << "Reporting point (1): " 
                                 << (*slIter)->get_left_end() << "\n";)
	      //*out = (*slIter)->get_left_end(); ++out;
	      add_point_to_output((*slIter)->get_left_end(), out);
	      m_found_intersection = true;
	    }
	  }
	  ++slIter;
	}   
      }
      ++vciter;
    }
  }
  /*! For each left-curve, if it is the "last" subcurve, i.e., the 
   * event point is the right-edge of the original curve, the 
   * last sub curve is created and added to the result. Otherwise
   * the curve is added as is to the result.
   */
  template <class OutpoutIterator>
  void handle_left_curves(OutpoutIterator out, SweepLineGetPoints &tag)
  {
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;)
    SL_DEBUG(m_currentEvent->Print();)
    const Point_2 &eventPoint = m_currentEvent->get_point();
    if ( !m_currentEvent->has_left_curves())
    {
      if (m_includeEndPoints || 
	  m_currentEvent->is_internal_intersection_point())
      {
	SL_DEBUG(std::cout << "Reporting point (2): " << eventPoint << "\n";)
	  //*out = eventPoint; ++out;    
	add_point_to_output(eventPoint, out);
	m_found_intersection = true;
      }
      return;
    }

    // delete the curve from the status line
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    m_currentPos = m_prevPos;
    m_use_hint_for_erase = false;
    while ( leftCurveIter != m_currentEvent->left_curves_end() )
    {
      // before deleting check new neighbors that will become after deletion
      remove_curve_from_status_line(*leftCurveIter);
      m_use_hint_for_erase = true;
      PRINT_ERASE((*leftCurveIter));
      m_currentPos = m_prevPos;
      Subcurve *leftCurve = *leftCurveIter; 
      leftCurve->set_last_point(eventPoint);
      ++leftCurveIter;
    }

    if ( m_includeEndPoints || 
	 m_currentEvent->is_internal_intersection_point() )
    {	
      SL_DEBUG(std::cout << "Reporting point (3): " << eventPoint << "\n";)
      //*out = eventPoint; ++out;
      add_point_to_output(eventPoint, out);
      m_found_intersection = true;
    }
  }

  
#ifndef NDEBUG
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
  void PrintVerticals();
#endif

protected:
  /*! a pointer to a traits object */
  Traits *m_traits;

  /*! an indication to whether the traits should be deleted in the distructor
   */
  bool m_traitsOwner;

  /*! if false, overlapping subcurves are reported only one. 
    Otherwise, they are reported as many times as they appeard. */
  bool m_overlapping;

  /*! if true, end points are reported as intersection points */
  bool m_includeEndPoints;

  /*! this is true when get_intersection_points() is called */
  bool m_stop_at_first_int;

  /*! used to hold all event, processed or not, so that they can 
`     all be deleted at the end */
  EventPtrContainer m_events; 

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves, as created on the fly */
  SubCurveList m_subCurves;

  /*! The status line */
  StatusLine *m_statusLine;

  /*! a struct that holds params associated with the curve compare functor */
  CompareParams *m_comp_param;

  /*! A reference point that is used for comapring curves. It is used
      when inserting/erasing curves from the status line. */
  Point_2 m_currentPos;

  /*! A reference point that is used for comapring curves */
  Point_2 m_prevPos;

  /*! The current position (in X) of the sweep line */
  Point_2 m_sweepLinePos;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_monotone_curve_2> m_xcurves;

  /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! a queue that holds all the events that have the same x coordinate as 
      the status line. */
  EventList m_miniq;

  /*! a list of vertical curves at the x coordinate of the current event 
      point.*/
  SubCurveList m_verticals;
  CurveList m_verticalSubCurves;

  /*! when an intersection point is found this is turned to true */
  bool m_found_intersection;

  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;

  /*! A an indication wether the hint can be used to erase from the status
    line */
  bool m_use_hint_for_erase;
  SubCurveList m_tmpOut;
  Point_2 m_lastReportedPoint;
  bool is_first_point;

  /*! a counter the is used to assign unique ids to the curves. */
  int m_curveId;

#ifndef NDEBUG
  int m_eventId;
#endif

  bool CurveStartsAtCurve(Subcurve *one, Subcurve *two)
  {
    SL_DEBUG(std::cout << "CurveStartsAtCurve: \n";)
    SL_DEBUG(std::cout << one->get_curve() << "\n" << two->get_curve() 
                       << "\n";)

    if ( m_traits->point_equal(one->get_left_end(),two->get_left_end()))
      return false;

    if ( !m_traits->point_equal(one->get_left_end(),
				m_currentEvent->get_point()) )
      return false;
    if ( m_traits->curve_compare_y_at_x(one->get_left_end(), 
					two->get_curve()) == EQUAL )
      return true;
    return false;
  }

};

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
~Sweep_line_2_impl() 
{
  if ( m_traitsOwner ) delete m_traits;

  for ( SubCurveListIter sci = m_subCurves.begin() ; 
	sci != m_subCurves.end() ; ++sci)
  {
    delete *sci;
  }

  for ( EventPtrContainerIter ei = m_events.begin();
	ei != m_events.end() ; ++ei)
  {
    delete *ei;
  }
  delete m_queue;
  delete m_statusLine;
  delete m_comp_param;
}



/*! initializes the data structures to work with:
 *  - x-monotonize the inf\put curves
 *  - for each end point of each curve create an event
 *  - initialize the event queue
 *  -
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
init(CurveInputIterator begin, CurveInputIterator end)
{
  PointLess pred(m_traits);
  m_queue = new EventQueue(pred);
  m_comp_param = new CompareParams(m_traits);
  StatusLineCurveLess slcurveless(m_comp_param);
  m_statusLine = new StatusLine(slcurveless);
  m_status_line_insert_hint = m_statusLine->begin();

#ifndef NDEBUG
  m_eventId = 0;
#endif
  m_curveId = 0;

  int count = 0;
  CurveInputIterator iter;
  for ( iter = begin ; iter != end ; ++iter)
  {
    std::list<X_monotone_curve_2> xcurves;
    m_traits->curve_make_x_monotone(*iter, std::back_inserter(xcurves));
    SL_DEBUG(
             std::cout << "curve " << *iter << " was split into " 
             << xcurves.size() << " curves." << std::endl;
             )

    for (typename std::list<X_monotone_curve_2>::iterator i = xcurves.begin();
         i != xcurves.end() ; ++i )
    {
      m_xcurves.push_back(*i);
      init_curve(m_xcurves[count]);
      count++;
    }
  }
}



/*! Given an x-monotone curve, create events for each end (if 
 *  one doesn't exist already). 
 *  For each curve create a Subcurve instance.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
init_curve(X_monotone_curve_2 &curve)
{
  const Point_2 &source = m_traits->curve_source(curve);
  const Point_2 &target = m_traits->curve_target(curve);
  Event *e = 0;
  
  Subcurve *subCv = new Subcurve(m_curveId++, curve, &m_currentPos, m_traits);
  m_subCurves.push_back(subCv);
  
  // handle the source point
  EventQueueIter eventIter = m_queue->find(&source); 
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << source << " already exists\n";)
    e = eventIter->second;
  } else  {
    e = new Event(source, m_traits); 
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(&(e->get_point()), e));
  }
  e->add_curve(subCv);
  PRINT_NEW_EVENT(source, e);
    
  // handle the target point
  eventIter = m_queue->find(&target);
  if ( eventIter != m_queue->end() ) {
    SL_DEBUG(std::cout << "event " << target << " already exists\n";)
    e = eventIter->second;
  } else  {
    e = new Event(target, m_traits); 
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    m_events.push_back(e);
    m_queue->insert(EventQueueValueType(&(e->get_point()), e));
  }
  e->add_curve(subCv);
  PRINT_NEW_EVENT(target, e);
}




/*!
 * Handles the degenerate case of vertical curves. Most of the cases
 * that occur with vertical curves are handled by this method and 
 * handle_vertical_curve_top method.
 * When the current event is the bottom end of a vertical curve, we look
 * for intersection points between the vertical curve and any curve
 * in the status line that in the y-range that is defined by the bottom 
 * and top ends of the vertical curve. When those are found, we create
 * new events, unless ones already exist, in which case we update the events.
 * 
 * @param tag a tag that indicates the version of this method
 * \sa handle_vertical_curve_top
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
handle_vertical_curve_bottom(SweepLineGetSubCurves &tag)
{
  SL_DEBUG(std::cout<<"\nhandle_vertical_curve_bottom... ("
                    << m_currentEvent->get_point() << ")\n";)
  if ( !m_currentEvent->does_contain_vertical_curve() )
  {
    SL_DEBUG(std::cout<<" - not vertical - exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
  VerticalCurveListIter vciter = vcurves.begin();
  const Point_2 &currentPoint = m_currentEvent->get_point();

  SL_DEBUG(std::cout << vcurves.size() << " vertical curves in event\n";)
  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    SL_DEBUG(std::cout << "working on " << vcurve->get_curve() << "\n";)
    if ( vcurve->is_top_end(currentPoint))
    {
      vciter++;
      continue;
    }
    
    SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";)
    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
    if ( slIter == m_statusLine->end() ) {
      SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";)
      vciter++;
      continue;
    }    
    
    SL_DEBUG(std::cout<<"starting at curve \n";)
    SL_DEBUG((*slIter)->Print();)
    const Point_2 &topEnd = vcurve->get_top_end();
    EventQueueIter topEndEventIter = m_queue->find(&topEnd);
    CGAL_assertion(topEndEventIter!=m_queue->end());
    Event *topEndEvent = topEndEventIter->second;

    bool lastEventCreatedHere = false;
    Event *prevEvent = 0;

    while (slIter != m_statusLine->end() &&
	   (! m_traits->point_in_x_range((*slIter)->get_curve(), 
					    topEnd) ||
	    m_traits->curve_compare_y_at_x(
	           topEnd, (*slIter)->get_curve()) != SMALLER) &&
	   (! m_traits->point_in_x_range((*slIter)->get_curve(), 
					 currentPoint) ||
	    m_traits->curve_compare_y_at_x(currentPoint,
                                           (*slIter)->get_curve()) != LARGER))
    {
      SL_DEBUG(std::cout<<"intersecting with \n";)
      SL_DEBUG((*slIter)->Print();) 
	
      if ( handle_vertical_curve_x_at_end(vcurve, *slIter, topEndEvent, tag)) {
	++slIter;
	continue;
      }
      
      // handle a curve that goes through the interior of the vertical curve
      const X_monotone_curve_2 &cv1 = vcurve->get_curve();
      const X_monotone_curve_2 &cv2 = (*slIter)->get_curve();

      Object res = m_traits->nearest_intersection_to_right(cv1, cv2,
                                                           currentPoint);
      CGAL_assertion(!res.is_empty());
      Point_2 xp;
      if (!CGAL::assign(xp, res))
        CGAL_assertion(0);
      
      EventQueueIter eqi = m_queue->find(&xp);
      Event *e = 0;
      if ( eqi == m_queue->end() )
      {
	e = new Event(xp, m_traits); 
#ifndef NDEBUG
	e->id = m_eventId++;
#endif
	m_events.push_back(e);
	
	e->add_curve_to_left(*slIter, m_sweepLinePos);
	e->add_curve_to_right(*slIter);
	PRINT_NEW_EVENT(xp, e);
	m_queue->insert(EventQueueValueType(&(e->get_point()), e));

	lastEventCreatedHere = true;

      } else {
	e = eqi->second;
	
	// the only time we need to update the event is when the event
	// is created here (which also includes overlapping curves)
	if ( e == prevEvent ) {
	  if ( lastEventCreatedHere )
	  {
	    if ( !(*slIter)->is_left_end(xp) ) 
	      e->add_curve_to_left(*slIter, m_sweepLinePos);
	    if ( !(*slIter)->is_right_end(xp) ) 
	      e->add_curve_to_right(*slIter);
	  } 
	}
	else
	  lastEventCreatedHere = false;

	SL_DEBUG(std::cout << "Updating event \n";)
	SL_DEBUG(e->Print();)
      }
      
      topEndEvent->add_vertical_curve_x_point(xp);
      ++slIter;
      prevEvent = e;
    }    
    vciter++;
  }

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}

/*!
 * Handles overlapping vertical curves. 
 * If the current event point does not contain vertical curves, nothing is
 * done here.
 * Fo the current event point, we go through the list of vertical curves 
 * defined in the same x coordinate (m_verticals). For each curve, we check 
 * if the event point is in the interior of the vertical curve. If so, 
 * the event is set to be an intersection point (between the two 
 * vertical curves). 
 * While going through the vertical curves, if we reach a curve that the 
 * event point is above the curve, we remove the curve from the list.
 * 
 * Finally, we go thorugh the vertical curves of the event. If the event 
 * point is the bottom end of a vertical curve, we add the vertical curve 
 * to the list of vertical curves (m_verticals).
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
handle_vertical_overlap_curves()
{
  SL_DEBUG(std::cout<<"\nhandle_vertical_overlap_curves... (" 
                    << m_currentEvent->get_point() << ")";)

  if ( !m_currentEvent->does_contain_vertical_curve() ) {
    SL_DEBUG(std::cout << "no vertical - exiting\n";)
    return;
  }
  SL_DEBUG(std::cout << "\n";)
  SL_DEBUG(PrintVerticals();)

  const Point_2 &point = m_currentEvent->get_point();
  SubCurveListIter iter = m_verticals.begin();
  while ( iter != m_verticals.end() )
  {
    Subcurve *curve = *iter;
   
    if (m_traits->point_in_x_range(curve->get_curve(), point) &&
	m_traits->curve_compare_y_at_x(point, curve->get_curve()) == LARGER)
    {
      iter = m_verticals.erase(iter);

    } else if (!curve->is_end_point(point)) {

      EventQueueIter eventIter = m_queue->find(&(curve->get_top_end()));
      CGAL_assertion(eventIter!=m_queue->end());
      (eventIter->second)->add_vertical_curve_x_point(point, true);
      m_currentEvent->mark_internal_intersection_point();
      ++iter;
    } else {
      ++iter;
    }
  }

  VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
  VerticalCurveListIter vciter = vcurves.begin();
  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    if ( vcurve->is_bottom_end(point) ) {
      m_verticals.push_back(vcurve);
    }
    ++vciter;
  }
}




/*!
 * Perform intersection between the specified curve and all curves in the 
 * given group of curves.
 */ 
template <class CurveInputIterator, class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
intersect_curve_group(Subcurve *c1, SubCurveList &mylist)
{
  m_tmpOut.clear();
  SL_DEBUG(std::cout << "intersecting with " << mylist.size() << " curves\n";)
  SubCurveListIter i = mylist.begin();
  while ( i != mylist.end())
  {

    intersect(c1, *i);
    ++i;
  }
}



/*!
 * When a curve is removed from the status line for good, its top and
 * bottom neighbors become neighbors. This method finds these cases and
 * looks for the itnersection point, if one exists.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return an iterator to the position where the curve will be removed from.
 */
template <class CurveInputIterator, class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
remove_curve_from_status_line(Subcurve *leftCurve)
{
  SL_DEBUG(std::cout << "remove_curve_from_status_line\n";)
  SL_DEBUG(PrintStatusLine();)
  SL_DEBUG(leftCurve->Print();)

  StatusLineIter sliter;
  sliter = leftCurve->get_hint(); 

  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint;
  if ( !leftCurve->is_end_point(m_currentEvent->get_point())) {
    m_statusLine->erase(sliter);
    SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
    return;
  }

  m_currentPos = m_prevPos;
  CGAL_assertion(sliter!=m_statusLine->end());
  StatusLineIter end = m_statusLine->end(); --end;
  if ( sliter != m_statusLine->begin() && sliter != end ) 
  {
    SubCurveList mylist;
    StatusLineIter prev = sliter; --prev;
    
    // collect all curves that overlap with *prev
    StatusLineIter tmp = prev;
    mylist.push_back(*prev);
    while ( tmp != m_statusLine->begin() ) 
    {
      --tmp;
      if ( do_curves_overlap(*prev, *tmp))
	mylist.push_back(*tmp);
      else
	break;
    }
    
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with the the *prev curve and all overlaps
    tmp = next;
    intersect_curve_group(*tmp, mylist);
    // if there are curves that overlap with the *next curve, intersect
    // them with the *prev curve and all overlaps
    ++tmp;
    while ( tmp != m_statusLine->end() ) 
    {
      if ( do_curves_overlap(*next, *tmp))
      {
	intersect_curve_group(*tmp, mylist);
	++tmp;
      }
      else 
	break;
    }
  }
  m_statusLine->erase(sliter);
  SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
} 

/*! 
 * Finds intersection between two curves. 
 * If the two curves intersect, create a new event (or use the event that 
 * already exits in the intersection point) and insert the curves to the
 * event.
 * @param curve1 a pointer to the first curve
 * @param curve2 a pointer to the second curve
 * @return true if the two curves overlap.
*/
template <class CurveInputIterator, class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
intersect(Subcurve *c1, Subcurve *c2)
{
  SL_DEBUG(std::cout << "Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  SL_DEBUG(std::cout << "\t";)
  SL_DEBUG(c2->Print();)
  SL_DEBUG(std::cout << "\n";)
  SL_DEBUG(std::cout << "relative to " << m_currentEvent->get_point() << "\n";)

  if ( c1->getId() == c2->getId() ) {
    SL_DEBUG(std::cout << "same curve, returning....\n";)
    return false;
  }

  Subcurve *scv1 = c1;
  Subcurve *scv2 = c2;
  const X_monotone_curve_2 &cv1 = scv1->get_curve();
  const X_monotone_curve_2 &cv2 = scv2->get_curve();

  bool isOverlap = false;

  Object res =
    m_traits->nearest_intersection_to_right(cv1, cv2, 
                                            m_currentEvent->get_point());
  if (!res.is_empty())
  {
    Point_2 xp;
    if (!CGAL::assign(xp, res)) {
      X_monotone_curve_2 cv;
      if (CGAL::assign(cv, res)) {
        xp = m_traits->curve_source(cv);
        Point_2 xp1 = m_traits->curve_target(cv);
        if ( m_traits->compare_x(xp1, xp) == LARGER )
          xp = xp1;
        SL_DEBUG(std::cout << "overlap detected\n";)
          isOverlap = true;
      }
    }

    SL_DEBUG(
      std::cout << " a new event is created between:\n\t";
      scv1->Print();
      std::cout << "\t";
      scv2->Print();
      std::cout << "\trelative to ("
                << m_sweepLinePos << ")\n\t at (" 
                << xp << ")" << std::endl;
    )

    // check to see if an event at this point already exists...
    EventQueueIter eqi = m_queue->find(&xp);
    Event *e = 0;
    if ( eqi == m_queue->end() )
    {
      e = new Event(xp, m_traits); 
#ifndef NDEBUG
      e->id = m_eventId++;
#endif
      m_events.push_back(e);
      
      e->add_curve_to_left(c1, m_sweepLinePos);
      e->add_curve_to_left(c2, m_sweepLinePos);
      
      e->add_curve_to_right(c1);
      e->add_curve_to_right(c2);
      
      PRINT_NEW_EVENT(xp, e);
      m_queue->insert(EventQueueValueType(&(e->get_point()), e));
      return isOverlap;
    } else 
    {
      SL_DEBUG(std::cout << "event already exists,updating.. (" << xp <<")\n";)
      e = eqi->second;
      e->add_curve_to_left(c1, m_sweepLinePos);
      if ( !scv1->is_end_point(xp)) { 
	e->add_curve_to_right(c1);
      }
      e->add_curve_to_left(c2, m_sweepLinePos);
      if ( !scv2->is_end_point(xp) ) {
	e->add_curve_to_right(c2);
      }
      SL_DEBUG(e->Print();)
    }
    return isOverlap;
  } 
  SL_DEBUG(std::cout << "not found 2\n";)
  return isOverlap;
}



template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
is_internal_x_point(const Point_2 &p)
{
  EventListIter itt = m_miniq.begin();
  while ( itt != m_miniq.end() )
  {
    if ( m_traits->point_equal(p, (*itt)->get_point())) 
    {
      if ((*itt)->is_internal_intersection_point()) {
	return true;
      }
      (*itt)->mark_internal_intersection_point(); 
                                     // this is to handle cases: |/ .
      return false;                  // (test 50/51)             |\ .
    } 
    ++itt;
  }
  CGAL_assertion(0);
  return false;
}



/*!
 * Handles the case in which a curve ont he status line passes through
 * one of the end points of the vertical curve.
 *
 * @param vcurve the vertical curve we are dealing with
 * @param curve a cerve that intersects with the vertical curve
 * @param topEndEvent the event attached to the top end of the vertical curve
 * @param tag 
 * @return returns true if the curve passed through one of the ends of the 
 *              vertical curve. Returns false otherwise.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
handle_vertical_curve_x_at_end(Subcurve *vcurve, Subcurve *curve, 
			  Event *topEndEvent, SweepLineGetSubCurves tag)
{
  const Point_2 &topEnd = vcurve->get_top_end();
  // handle a curve that goes through the top point of the vertical curve
  if (m_traits->point_in_x_range(curve->get_curve(), topEnd) &&
      m_traits->curve_compare_y_at_x(topEnd, curve->get_curve()) == EQUAL)
  {
    if ( !curve->is_left_end(topEnd)) {
      topEndEvent->add_curve_to_left(curve, m_prevPos);
    }
    if ( ! curve->is_right_end(topEnd)) {
      topEndEvent->add_curve_to_right(curve);
    }
    return true;
  } 
  
  // handle a curve that goes through the bottom point of the vertical curve
  const Point_2 &currentPoint = m_currentEvent->get_point();
  if (m_traits->point_in_x_range((curve)->get_curve(), currentPoint) &&
      m_traits->curve_compare_y_at_x(currentPoint, (curve)->get_curve()) ==
      EQUAL)
  {
    if ( !(curve)->is_left_end(currentPoint)) {
      m_currentEvent->add_curve_to_left(curve, m_prevPos);
    }
    if ( ! (curve)->is_right_end(currentPoint)) {
      m_currentEvent->add_curve_to_right(curve);
    }
    return true;;
  }
  return false;
}



template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
do_curves_overlap(Subcurve *c1, Subcurve *c2)
{
  SL_DEBUG(std::cout << "do_curves_overlap " << m_sweepLinePos << "\n" 
	             << "\t" << c1->get_curve() << "\n"
	             << "\t" << c2->get_curve() << "\n";)

  const Point_2 *p = &(c2->get_last_point());
  if ( m_traits->compare_x(c1->get_last_point(), 
			   c2->get_last_point()) == LARGER )
    p = &(c1->get_last_point());

  if ((m_traits->curves_compare_y_at_x(c1->get_curve(),
				       c2->get_curve(),
				       *p) != EQUAL))
    return false;

  if ( m_traits->curves_overlap(c1->get_curve(),c2->get_curve()) )
    return true;

  return false;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
similar_curves(const X_monotone_curve_2 &a, const X_monotone_curve_2 &b)
{
  if ( m_traits->curve_equal(a, b))
    return true;
  return false;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
vertical_subcurve_exists(const X_monotone_curve_2 &a)
{
  for (typename std::list<X_monotone_curve_2>::iterator iter =
         m_verticalSubCurves.begin() ;
       iter != m_verticalSubCurves.end() ; ++iter)
  {
    if (similar_curves(*iter, a)) 
      return true;
  }
  return false;
}



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//                 point implementation of the methods                    //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


/*!
 * Handles the degenerate case of vertical curves. Most of the cases
 * that occur with vertical curves are handled by this method and the
 * handle_vertical_curve_top method.
 *
 * When the current event is the bottom end of a vertical curve, we look
 * for intersection points between the vertical curve and any curve
 * in the status line that in the y-range that is defined by the bottom 
 * and top ends of the vertical curve. When those are found, we create
 * new events, unless ones already exist, in which case we update the events.
 * 
 * @param tag a tag that indicates the version of this method
 * \sa handle_vertical_curve_top
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
handle_vertical_curve_bottom(SweepLineGetPoints &tag)
{
  SL_DEBUG(std::cout<<"handle_vertical_curve_bottom... ";)
  if ( !m_currentEvent->does_contain_vertical_curve() )
  {
    SL_DEBUG(std::cout<<"exiting\n ";)
    return;
  }
  SL_DEBUG(std::cout<<"\n ";)

  VerticalCurveList &vcurves = m_currentEvent->get_vertical_curves();
  VerticalCurveListIter vciter = vcurves.begin();
  const Point_2 &currentPoint = m_currentEvent->get_point();

  while ( vciter != vcurves.end() )
  {
    Subcurve *vcurve = *vciter;
    if ( vcurve->is_top_end(currentPoint))
    {
      ++vciter;
      continue;
    }

    SL_DEBUG(std::cout<<"handling bottom point of vertical curve\n";);
    StatusLineIter slIter = m_statusLine->lower_bound(vcurve);
    if ( slIter == m_statusLine->end() ) {
      SL_DEBUG(std::cout<<"no curves intersecting. exiting\n";);
      ++vciter;
      continue;
    }    

    SL_DEBUG(std::cout<<"starting at curve \n";);
    SL_DEBUG((*slIter)->Print(););

    const Point_2 &topEnd = vcurve->get_top_end();
    EventQueueIter topEndEventIter = m_queue->find(&topEnd);
    CGAL_assertion(topEndEventIter!=m_queue->end());
    Event *topEndEvent = topEndEventIter->second;

    while (slIter != m_statusLine->end() &&
	   (! m_traits->point_in_x_range((*slIter)->get_curve(), 
					    topEnd) ||
	    m_traits->curve_compare_y_at_x(topEnd, (*slIter)->get_curve()) !=
            SMALLER) &&
	   (! m_traits->point_in_x_range((*slIter)->get_curve(), 
					    currentPoint) ||
	    m_traits->curve_compare_y_at_x(currentPoint,
                                           (*slIter)->get_curve()) != LARGER))
    {
      SL_DEBUG(std::cout<<"intersecting with \n";)
      SL_DEBUG((*slIter)->Print();) 
	
      if ( handle_vertical_curve_x_at_end(vcurve, *slIter, topEndEvent, tag))
      {
	++slIter;
	continue;
      }

      Object res = 
	m_traits->nearest_intersection_to_right(vcurve->get_curve(), 
						(*slIter)->get_curve(), 
						currentPoint);
      CGAL_assertion(!res.is_empty());
      Point_2 xp;
      if (!CGAL::assign(xp, res))
        CGAL_assertion(0);

      EventQueueIter eqi = m_queue->find(&xp); 
      Event *e = 0;
      if ( eqi == m_queue->end() )
      {
	e = new Event(xp, m_traits); 
#ifndef NDEBUG
	e->id = m_eventId++;
#endif
	m_events.push_back(e);
      
	e->add_curve_to_left(*slIter, m_sweepLinePos);
	e->add_curve_to_right(*slIter);

	PRINT_NEW_EVENT(xp, e);
	m_queue->insert(EventQueueValueType(&(e->get_point()), e)); 
      } else {
	e = eqi->second;
	e->mark_internal_intersection_point();
	SL_DEBUG(std::cout << "Updating event \n";)
	SL_DEBUG(e->Print();)
	e->add_curve(vcurve); // test41
	e->add_curve_to_left(*slIter, (*slIter)->get_left_end());
	if ( m_traits->compare_x((*slIter)->get_right_end(), 
				 m_currentEvent->get_point()) == LARGER )
	  e->add_curve_to_right(*slIter);
      }
      
      topEndEvent->add_vertical_curve_x_point(xp);
      ++slIter;
    }    
    ++vciter;
  }

  SL_DEBUG(std::cout<<"Done Handling vertical\n";)
}



/*!
 * Handles the case in which a curve ont he status line passes through
 * one of the end points of the vertical curve.
 *
 * @param vcurve the vertical curve we are dealing with
 * @param curve a cerve that intersects with the vertical curve
 * @param topEndEvent the event attached to the top end of the vertical curve
 * @param tag 
 * @return returns true if the curve passed through one of the ends of the 
 *              vertical curve. Returns false otherwise.
 */
template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline bool
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
handle_vertical_curve_x_at_end(Subcurve *vcurve, Subcurve *curve, 
			  Event *topEndEvent, SweepLineGetPoints tag)
{
  const Point_2 &topEnd = vcurve->get_top_end();
  // handle a curve that goes through the top point of the vertical curve
  if (m_traits->point_in_x_range((curve)->get_curve(), topEnd) &&
      m_traits->curve_compare_y_at_x(topEnd, (curve)->get_curve()) == EQUAL)
  {
    if ( !curve->is_end_point(topEnd)) {
      topEndEvent->mark_internal_intersection_point();
    }
    return true;
  } 

  // handle a curve that goes through the bottom point of the vertical curve
  if (m_traits->point_in_x_range((curve)->get_curve(),
				    m_currentEvent->get_point()) &&
      m_traits->curve_compare_y_at_x(m_currentEvent->get_point(),
                                     (curve)->get_curve()) == EQUAL)
  {
    if ( !curve->is_end_point(m_currentEvent->get_point())) {
      m_currentEvent->mark_internal_intersection_point();
    }
    return true;
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//                         DEBUG UTILITIES                                //
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

#ifndef NDEBUG

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintEventQueue()
{
  SL_DEBUG(std::cout << std::endl << "Event queue: " << std::endl;)
  EventQueueIter iter = m_queue->begin();
  while ( iter != m_queue->end() )
  {
    SL_DEBUG(std::cout << "Point (" << iter->first << ")" << std::endl;)
    Event *e = iter->second;
    e->Print();
    ++iter;
  }
  SL_DEBUG(std::cout << "--------------------------------" << std::endl;)
}

template <class CurveInputIterator,  class SweepLineTraits_2,
          class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintSubCurves()
{
  SL_DEBUG(std::cout << std::endl << "Sub curves: " << std::endl;)
  SubCurveListIter iter = m_subCurves.begin();
  while ( iter != m_subCurves.end() )
  {
    (*iter)->Print();
    ++iter;
  }
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintStatusLine()
{
  if ( m_statusLine->size() == 0) {
    std::cout << std::endl << "Status line: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Status line: (" 
	    << m_currentPos << ")" << std::endl;
  StatusLineIter iter = m_statusLine->begin();
  while ( iter != m_statusLine->end() )
  {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Status line - end" << std::endl;
}

template <class CurveInputIterator,  class SweepLineTraits_2,
         class SweepEvent, class CurveWrap>
inline void 
Sweep_line_2_impl<CurveInputIterator,SweepLineTraits_2,SweepEvent,CurveWrap>::
PrintVerticals()
{
  if ( m_verticals.size() == 0) {
    std::cout << std::endl << "Verticals: empty" << std::endl;
    return;
  }
  std::cout << std::endl << "Verticals: " << m_verticals.size() << " (" 
	    << m_currentEvent->get_point() << ")" << std::endl;
  SubCurveListIter iter = m_verticals.begin();
  while ( iter != m_verticals.end() )
  {
    (*iter)->Print();
    ++iter;
  }
  std::cout << "Verticals - end" << std::endl;
}

#endif

CGAL_END_NAMESPACE

#endif
