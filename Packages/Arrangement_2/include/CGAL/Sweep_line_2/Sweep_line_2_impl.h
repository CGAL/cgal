// Copyright (c) 2005  Tel-Aviv University (Israel).
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

#ifndef CGAL_SWEEP_LINE_2_IMPL_H
#define CGAL_SWEEP_LINE_2_IMPL_H

#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <CGAL/assertions.h>
#include <CGAL/memory.h>
#include <CGAL/Object.h>
#include <CGAL/Sweep_line_2/Sweep_line_functors.h>
#include <CGAL/Sweep_line_2/Sweep_line_traits.h>
#include <CGAL/Sweep_line_2/Sweep_line_event.h>
#include <CGAL/Sweep_line_2/Sweep_line_subcurve.h>
#include <CGAL/Sweep_line_2/Sweep_line_rb_tree.h>



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
                          std::cout << "    currentPos = "  \
                                    << m_currentEvent->get_point() \
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

  An extension of this algorithm is used to produce a planar map containing
  the input curves efficiently. \sa Arr_aggregate_insert.

  General flow:
  After the initialization stage, the events are handled from left to right.

  For each event

    First pass - handles special cases in which curves start or end 
                 at the interior of another curve
    Handle left curves - iterate over the curves that intersect 
                 at the event point and defined to the left of the 
                 event. 
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

*/
template < class SweepLineTraits_2,
           class SweepEvent,
           class CurveWrap,
           class SweepVisitor,
           typename Allocator = CGAL_ALLOCATOR(int) >
class Sweep_line_2_impl
{
public:

  typedef SweepLineTraits_2 Traits;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

  typedef SweepEvent Event;
  typedef Point_less_functor<Traits> PointLess;
  typedef std::map< Point_2 , Event*, PointLess> EventQueue; 
  typedef typename EventQueue::iterator EventQueueIter;
  typedef typename EventQueue::value_type EventQueueValueType;

  typedef typename Event::SubCurveIter EventCurveIter;

  typedef CurveWrap Subcurve;
  typedef std::list<Subcurve*> SubCurveList;
  typedef typename SubCurveList::iterator SubCurveListIter;
   
  typedef Status_line_curve_less_functor<Traits, Subcurve> StatusLineCurveLess;
  typedef Red_black_tree<Subcurve*, StatusLineCurveLess, 
			 CGAL_ALLOCATOR(int)>              StatusLine;
  typedef typename StatusLine::iterator StatusLineIter;
 
  typedef typename Allocator::template rebind<Event>      EventAlloc_rebind;
  typedef typename EventAlloc_rebind::other               EventAlloc;

  typedef typename Allocator::template rebind<Subcurve>   SubcurveAlloc_rebind;
  typedef typename SubcurveAlloc_rebind::other            SubCurveAlloc;

  typedef Curves_pair<Subcurve>                           CurvesPair;
  typedef Curves_pair_less_functor<Subcurve>              CurvesPairLess;
  typedef std::set<CurvesPair,CurvesPairLess>             CurvesTable;

  typedef random_access_input_iterator<std::vector<Object> > vector_inserter;

  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2_impl (SweepVisitor* visitor) :
      //m_traits(new Traits()),
      m_sweep_line_traits(&m_traits),
      //m_traitsOwner(true),
      m_statusLineCurveLess(&m_traits),
      m_queue(new EventQueue(PointLess(&m_traits))),
      m_statusLine(m_statusLineCurveLess),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
#ifndef NDEBUG
      m_eventId(0),
#endif
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }


  /*!
   * Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2_impl(const Traits *traits, SweepVisitor* visitor) :
      m_traits(*traits),
      m_sweep_line_traits(&m_traits),
      //m_traitsOwner(false),
      m_statusLineCurveLess(&m_traits),
      m_queue(new EventQueue(PointLess(&m_traits))),
      m_statusLine(m_statusLineCurveLess),
      m_xcurves(0),
      m_status_line_insert_hint(m_statusLine.begin()),
      m_num_of_subCurves(0),
#ifndef NDEBUG
      m_eventId(0),
#endif
      m_visitor(visitor)
  {
    m_visitor->attach(this);
  }

  /*! Destrcutor. */
  virtual ~Sweep_line_2_impl();

  /*!
   * Initialize the sweep-line with a range of curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of CurveInputIterator is the traits' Curve_2.
   */
  template<class CurveInputIterator>
  void init(CurveInputIterator curves_begin, CurveInputIterator curves_end)
  {
    // Subdivide each input curve into x-monotone subcurves.
    for(CurveInputIterator iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits.make_x_monotone_2_object() (*iter,
					   std::back_inserter(m_xcurves));
    }

    // Initialize the x-monotone curves and give them indices.
    unsigned int   index = 0;

    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);

    for (index = 0; index < m_num_of_subCurves ; ++index)
    {
      _init_curve (m_xcurves[index], index);
    }

    return;
  }

  /*!
   * Initialize the sweep-line with a range of x-monotone curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for the range.
   * \pre The value-type of XCurveInputIterator is the traits'
   *      X_monotone_curve_2.
   */
  template<class XCurveInputIterator>
  void init_x_curves (XCurveInputIterator curves_begin,
		      XCurveInputIterator curves_end)
  {
    // Initialize the x-monotone curves and give them indices.
    unsigned int         index = 0;
    XCurveInputIterator  iter;

    m_num_of_subCurves = std::distance(curves_begin,curves_end);
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);

    for (iter = curves_begin; iter != curves_end; ++iter, ++index)
    {
      _init_curve (*iter,index);
    }

    return;
  }

  /*!
   * Initialize the sweep-line with a range of curves and a range of
   * x-monotone curves.
   * \param curves_begin An iterator for the first curve in the range.
   * \param curves_end A past-the-end iterator for this range.
   * \param xcurves_begin  An iterator for the first x-monotone curve in the
   *                      range.
   * \param xcurves_end A past-the-end iterator for this range.
   * \pre The value-type of CurveInputIterator is the traits' Curve_2, and the
   *      value-type of XCurveInputIterator is the traits' X_monotone_curve_2.
   */
  template<class CurveInputIterator, class XCurveInputIterator>
  void init(CurveInputIterator curves_begin, CurveInputIterator curves_end,
            XCurveInputIterator xcurves_begin, XCurveInputIterator xcurves_end)
  {
    // Subdivide each input curve into x-monotone subcurves.
    CurveInputIterator iter;

    for (iter = curves_begin; iter != curves_end; ++iter)
    {
      m_traits.make_x_monotone_2_object() (*iter,
					   std::back_inserter(m_xcurves));
    }

    // Add the input x-monotone curves.
    std::copy(xcurves_begin, xcurves_end, std::back_inserter(m_xcurves));

    // Initialize the x-monotone curves and give them indices.
    unsigned int    index;

    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);
    
    for(index = 0; index < m_num_of_subCurves ; ++index)
    {
      _init_curve (m_xcurves[index],index);
    }
  }

  /*!
   * Initialize the sweep-line with a range of event points and a range of
   * x-monotone curves.
   * \param xcurves_begin  An iterator for the first x-monotone curve in the
   *                      range.
   * \param xcurves_end A past-the-end iterator for this range.
   * \param points_begin An iterator for the first point in the range.
   * \param points_end A past-the-end iterator for this range.
   * \pre The value-type of XCurveInputIterator is the traits' 
   *      X_monotone_curve_2, and the value-type of PointInputIterator is the 
   *      traits' Point_2.
   */
  template<class CurveInputIterator, class PointInputIterator>
  void init (CurveInputIterator curves_begin, CurveInputIterator curves_end,
	     PointInputIterator points_begin, PointInputIterator points_end,
	     bool)
  {
    // Initialize the x-monotone curves and give them indices.
    unsigned int   index;

    std::copy(curves_begin, curves_end, std::back_inserter(m_xcurves));
    m_num_of_subCurves = m_xcurves.size();
    m_subCurves = m_subCurveAlloc.allocate(m_num_of_subCurves);

    for (index = 0; index < m_num_of_subCurves ; ++index)
    {
      _init_curve (m_xcurves[index], index);
    }

    // Initialize the events associated with the given points.
    PointInputIterator  iter;

    for (iter = points_begin; iter != points_end; ++iter)
    {
      _init_point (*iter);
    }

    return;
  }

  /*! Preform the main sweep-line loop. */
  void sweep()
  {
    // Looping over the events in the queue.
    EventQueueIter eventIter = m_queue->begin();

    while (eventIter != m_queue->end())
    {
      // Get the next event from the queue.
      m_currentEvent = eventIter->second;

      SL_DEBUG(std::cout << "------------- " 
                         << m_currentEvent->get_point() 
                         << " --------------"
                         << std::endl;
               PrintStatusLine();
               m_currentEvent->Print(););
      
      // Handle the subcurves that are to the left of the event point (i.e., 
      // subcurves that we are done with).
      _handle_left_curves();

      // Handle the subcurves to the right of the event point, reorder them
      // and test for intersections between them and their immediate neighbors
      // on the status line.
      _handle_right_curves();

      // Inform the visitor about the event.
      if (m_visitor->after_handle_event(m_currentEvent,
					m_status_line_insert_hint,
					m_is_event_on_above))
      {
	// It is possible o deallocate the event:
        deallocate_event(m_currentEvent);
      }

      // We are done with the current event - remove it from the queue.
      m_queue->erase(eventIter);
      eventIter = m_queue->begin();
    }

    return;
  }

  /*! Get an iterator for the first subcurve in the status line. */
  StatusLineIter  StatusLine_begin()
  {
    return m_statusLine.begin();
  }

  /*! Get a past-the-end iterator for the subcurves in the status line. */
  StatusLineIter  StatusLine_end()
  {
    return m_statusLine.end();
  }

protected:

  /*!
   * Initialize an event associated with a point.
   * \param p The given point.
   */
  void _init_point(const Point_2& pt)
  {
    // Try inserting the event to the queue and act according to the result:
    Event                                 *e;
    const std::pair<EventQueueIter, bool>& insertion_res =
      m_queue->insert(EventQueueValueType(pt,0));

    if (insertion_res.second == true)
    {
      // We have a new event:
      e = allocate_event(pt);
      m_visitor -> init_event(e);
     
    #ifndef NDEBUG
      e->id = m_eventId++;
    #endif
     
      (insertion_res.first)->second = e;
    }
    else
    {
      // The event already exsits:
      e = (insertion_res.first)->second;
      m_visitor -> init_event(e);
    }

    return;
  }
  
  /*!
   * Initialize an event associated with an x-monotone curve.
   * \param curve The given x-monotone curve.
   * \param indec Its unique index.
   */
  void _init_curve(X_monotone_curve_2 &curve,unsigned int index);
  
  /*! Handle the subcurve to the left of the current event point. */
  void _handle_left_curves()
  { 
    SL_DEBUG(std::cout << "Handling left curve" << std::endl;);

    // For each left-curve, if it is the "last" subcurve, i.e., the 
    // event point is the right-edge of the original curve, the 
    // last sub curve is created and added to the result. Otherwise
    // the curve is added as is to the result.
    if(! m_currentEvent->has_left_curves())
    {
      /* this block takes care of
      //
      //           /
      //          /
      //       --------
      //          \
      //           \
      */
      SL_DEBUG(std::cout << " - handling special case " << std::endl;)
                                                     
      const std::pair<StatusLineIter, bool>& res =
        m_statusLine.lower_bound (m_currentEvent->get_point(), 
				  m_statusLineCurveLess);
      m_status_line_insert_hint = res.first;

      if(res.second)
      {
	// The current event point starts at the interior of a curve at the 
	// y-structure (can also indicates overlap).
        if(!m_currentEvent -> has_right_curves())
        {
          m_is_event_on_above = true;
          m_visitor->before_handle_event(m_currentEvent);
          return;
        }
        m_is_event_on_above = false;
        Subcurve* sc = *m_status_line_insert_hint;
        X_monotone_curve_2 last_curve = sc->get_last_curve();
        m_currentEvent->add_curve_to_left(sc); 
        bool is_overlap = _add_curve_to_right(m_currentEvent, sc);
        X_monotone_curve_2 a,b;
       
        m_traits.split_2_object() (last_curve,
				   m_currentEvent->get_point(), 
				   a, b);
        ++m_status_line_insert_hint; 
        m_statusLine.remove_at(res.first);
        if(!is_overlap)
        {
          sc->set_last_curve(b);
        }
        m_visitor->before_handle_event(m_currentEvent);
        m_visitor->add_subcurve(a, sc);
        return;
      }
      m_is_event_on_above = false;
      m_visitor->before_handle_event(m_currentEvent);
      return;
    }
        
    m_is_event_on_above = false; 

    SL_DEBUG(std::cout<<"left curves before sorting: "<<"\n";);
    SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });
    _sort_left_curves();
    m_visitor->before_handle_event(m_currentEvent);

    SL_DEBUG(std::cout<<"left curves after sorting: "<<"\n";);
    SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });

     // Check if the curve should be removed for good.
    bool remove_for_good = false; 

    EventCurveIter left_iter = m_currentEvent->left_curves_begin();
    while(left_iter != m_currentEvent->left_curves_end())
    {
      Subcurve *leftCurve = *left_iter; 
    
      if((Event*)leftCurve->get_right_event() == m_currentEvent)
      {  
        remove_for_good = true;
        m_visitor->add_subcurve(leftCurve->get_last_curve(), leftCurve);

        if(leftCurve->get_orig_subcurve1() != NULL)
        {
          leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);
          leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);

          // clip the two subcurves according to end_overlap point
          Subcurve* res;
          if(( res = (Subcurve*) (leftCurve->get_orig_subcurve1() -> 
               clip(m_currentEvent))) != NULL)
          {
            _add_curve_to_right(m_currentEvent, res);

          }

          if(( res = (Subcurve*)(leftCurve->get_orig_subcurve2() -> 
               clip(m_currentEvent))) != NULL)
          {
            _add_curve_to_right(m_currentEvent, res);
          }
        }
      }
      else
      { 
        X_monotone_curve_2 a, b;
        const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
       
        m_traits.split_2_object()(lastCurve,
				  m_currentEvent->get_point(), 
				  a, b);
        m_visitor->add_subcurve(a, leftCurve);

        leftCurve->set_last_curve(b);
      }
       ++left_iter;

      // remove curve from the status line (also checks intersection 
      // between the neighbouring curves,only if the curve is removed for good)
      _remove_curve_from_status_line(leftCurve, remove_for_good);    
    }
    SL_DEBUG(std::cout << "Handling left curve END" << std::endl;);
      
    return;
  }

  /*!
   * Sort the left subcurves of an event point according to their order in
   * their status line (no geometric comprasions are needed).
   */
  void _sort_left_curves()
  {
    CGAL_assertion(m_currentEvent->has_left_curves());
    EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();
    while ( leftCurveIter != m_currentEvent->left_curves_end() )  
    {
      if( (*leftCurveIter)->get_overlap_subcurve() != NULL &&
             m_traits.compare_xy_2_object()(
             (*leftCurveIter)->get_overlap_subcurve()->get_left_end(),
             m_currentEvent->get_point()) == SMALLER)
      {
        Subcurve *leftCurve = (Subcurve*)((*leftCurveIter)->getSubcurve());
        m_currentEvent->replace_right_curve((*leftCurveIter), leftCurve);
        *leftCurveIter = leftCurve;
      }
      ++leftCurveIter;
    }

    Subcurve *curve = *(m_currentEvent->left_curves_begin());
    StatusLineIter slIter = curve->get_hint();
    CGAL_assertion(*slIter == curve);
   
    
    for (++slIter; slIter != m_statusLine.end(); ++slIter)
    {
      if(std::find(m_currentEvent->left_curves_begin(),
                   m_currentEvent->left_curves_end(),
                   *slIter) ==
                   m_currentEvent->left_curves_end())
         break;
    }
    StatusLineIter end (slIter);

    slIter = curve->get_hint();
    if(slIter == m_statusLine.begin())
    {
      m_currentEvent->replace_left_curves(slIter,end);
      return;
    }
    --slIter;
    for(;slIter != m_statusLine.begin(); --slIter)
    {
      if( std::find(m_currentEvent->left_curves_begin(),
                    m_currentEvent->left_curves_end(),
                    *slIter) == m_currentEvent->left_curves_end())
      {
        m_currentEvent->replace_left_curves(++slIter,end);
        return;
      }
    }
    if(std::find(m_currentEvent->left_curves_begin(),
                  m_currentEvent->left_curves_end(),
                  *slIter) == m_currentEvent->left_curves_end())
    {
      m_currentEvent->replace_left_curves(++slIter,end);;
    }
    else
    {
      m_currentEvent->replace_left_curves(slIter,end);
    }
  }

  /*! Handle the subcurve to the left of the current event point. */
  void _handle_right_curves()
  {
    SL_DEBUG(std::cout << "Handling right curves (" ;)
    SL_DEBUG(std::cout << m_currentEvent->get_point() << ")\n";)
    
    if(! m_currentEvent->has_right_curves())
      return;


    // Loop over the curves to the right of the status line and handle them:
    // - If we are at the beginning of the curve, we insert it to the status 
    //   line, then we look if it intersects any of its neighbors.
    // - If we are at an intersection point between two curves, we add them
    //   to the status line and attempt to intersect them with their neighbors.
    // - We also check to see if the two intersect again to the right of the 
    //   point.
    int numRightCurves = m_currentEvent->get_num_right_curves();
    if(numRightCurves == 1)
    {
      SL_DEBUG(std::cout << " - beginning of curve " << std::endl;);
      SL_DEBUG(
        Subcurve *tmp1 = *(m_currentEvent->right_curves_begin());
        PRINT_INSERT(tmp1);
               );

      StatusLineIter slIter = 
        m_statusLine.insert_predecessor(m_status_line_insert_hint, 
                                *(m_currentEvent->right_curves_begin()));
      
      (*(m_currentEvent->right_curves_begin()))->set_hint(slIter); 
     
      SL_DEBUG(PrintStatusLine(););
     
      // if this is the only curve on the status line, nothing else to do
      if ( m_statusLine.size() == 1 )
        return;

      StatusLineIter prev = slIter;
      StatusLineIter next = slIter;
      ++next;
      if ( slIter != m_statusLine.begin() )
      {
        --prev;
        _intersect(*prev, *(m_currentEvent->right_curves_begin()));
      }
      if ( next != m_statusLine.end() )
      { 
        _intersect(*(m_currentEvent->right_curves_begin()), *next);
      } 
    }
    else  // numRightCurves > 1
    {
      EventCurveIter firstOne = m_currentEvent->right_curves_begin();
      EventCurveIter lastOne = m_currentEvent->right_curves_end(); --lastOne;
      EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();
      PRINT_INSERT(*firstOne);

      StatusLineIter slIter = m_statusLine.insert_predecessor
	(m_status_line_insert_hint, *firstOne);
      ((Subcurve*)(*firstOne))->set_hint(slIter);
        
      SL_DEBUG(PrintStatusLine(););
      if ( slIter != m_statusLine.begin() )
      { 
        //  get the previous curve in the y-str
        StatusLineIter prev = slIter; --prev;
        _intersect(*prev, *slIter);
      }
      
      EventCurveIter currentOne = firstOne; ++currentOne;
      EventCurveIter prevOne = firstOne;
      while ( currentOne != rightCurveEnd )
      {
        PRINT_INSERT(*currentOne);
        slIter = m_statusLine.insert_predecessor
	  (m_status_line_insert_hint, *currentOne);
        ((Subcurve*)(*currentOne))->set_hint(slIter);
          
        SL_DEBUG(PrintStatusLine(););
      
        _intersect(*prevOne, *currentOne);
        prevOne = currentOne;
        ++currentOne;
      }        
      lastOne = currentOne; --lastOne;
        
      SL_DEBUG(PrintStatusLine(););
      StatusLineIter next = slIter; ++next;
      if ( next != m_statusLine.end() )
        _intersect( *prevOne,*next);
    }
  }

  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  bool _add_curve_to_right (Event* event, Subcurve* curve)
  {
    std::pair<bool, SubCurveListIter> pair_res = 
                                              event->add_curve_to_right(curve);

    if (! pair_res.first)
      // No overlap occurs:
      return (false);

    // An overlap occurs:
    // TODO: take care of polylines in which overlap can happen anywhere
    SL_DEBUG(std::cout<<"Overlap detected at right insertion...\n";);
    SubCurveListIter iter = pair_res.second;
      
    X_monotone_curve_2 overlap_cv;
    
    vector_inserter vi(m_x_objects);
    vector_inserter vi_end(m_x_objects);
    vi_end = m_traits.intersect_2_object()(curve->get_last_curve(),
					   (*iter)->get_last_curve(),
					   vi);
    
    CGAL::assign(overlap_cv, *vi);
    
    // Alocate a new Subcure for the overlap
    Subcurve *overlap_sc = m_subCurveAlloc.allocate(1);
    m_subCurveAlloc.construct(overlap_sc,m_masterSubcurve);
    overlap_sc->init(overlap_cv);
    m_overlap_subCurves.push_back(overlap_sc);
      
    // Get the right end of overlap_cv
    Point_2 end_overlap = m_traits.construct_max_vertex_2_object()(overlap_cv);
     
    //find the event assiciated with end_overlap point (right end point)
    EventQueueIter q_iter = m_queue->find( end_overlap );
    CGAL_assertion(q_iter!=m_queue->end());

    // set the members left event and right event of overlap_sc
    overlap_sc->set_left_event(event);
    overlap_sc->set_right_event((*q_iter).second);

    m_visitor->init_subcurve(overlap_sc);

    // Remove curve, *iter from the left curves of end_overlap event
    ((*q_iter).second)->remove_curve_from_left(curve);
    ((*q_iter).second)->remove_curve_from_left(*iter);

    // Add overlap_sc to the left curves
    ((*q_iter).second)->add_curve_to_left(overlap_sc);

    curve  -> set_overlap_subcurve(overlap_sc);
    (*iter)-> set_overlap_subcurve(overlap_sc);

    overlap_sc -> set_orig_subcurve1(*iter);
    overlap_sc -> set_orig_subcurve2(curve);  

    // Replace current sub-curve (*iter) with the new sub-curve
    (*iter) = overlap_sc;

    // Inidicate that an overlap has occured:
    return (true);
  }

  /*! Compute intersections between the two given curves. */ 
  void _intersect(Subcurve *c1, Subcurve *c2, bool after_remove = false);

  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line (Subcurve *leftCurve,
				       bool remove_for_good);

#ifndef NDEBUG
  void PrintEventQueue();
  void PrintSubCurves();
  void PrintStatusLine();
#endif

protected:

  /*! a  traits object */
  Traits m_traits;

  /*! an object that holds a static trait object 
   *  to be used by events and subcurves
   */
  Sweep_line_traits<Traits> m_sweep_line_traits;

  /*! Y-str comprasion functor */
  StatusLineCurveLess m_statusLineCurveLess;

  /*! the queue of events (intersection points) to handle */
  EventQueue *m_queue;

  /*! The subcurves array */
  Subcurve *m_subCurves;

  /*! The status line (Y-str) */
  StatusLine m_statusLine;

  /*! if non x-monotone are specified, this hold the x-monotone 
    curves created when splitting them into x-monotone curves. */
  std::vector<X_monotone_curve_2> m_xcurves;

  /*! a pointer to the current event */
  Event *m_currentEvent;

  /*! An iterator of the  status line that is used as a hint for inserts. */
  StatusLineIter m_status_line_insert_hint;

  /*! indicates if current event is on the interior of existing curve*/
  bool m_is_event_on_above;

  /*! An allocator for the events objects */
  EventAlloc m_eventAlloc;

  /*! An allocator for the Subcurve objects */
  SubCurveAlloc m_subCurveAlloc;

  /*! a master Event (created once by the constructor) for the allocator's 
   *  usgae. */
  Event m_masterEvent;

  /*! a master Subcurve (created once by the constructor) for the allocator's 
   *  usgae. */
  Subcurve m_masterSubcurve;

  /*! The num of subcurves  */
  unsigned int m_num_of_subCurves;

  /*! contains all of the new sub-curve creaed by overlap */
  SubCurveList m_overlap_subCurves;

#ifndef NDEBUG
  int m_eventId;
#endif

  /*! a pointer to the visitor object which will be notidifed during sweep */
  SweepVisitor* m_visitor;

  /*! a lookup table of pairs of Subcurves that ahve been intersected */
  //TODO: replace set with hash 
  CurvesTable m_curves_table;

  /*! a vector holds the intersection objects */
  std::vector<Object> m_x_objects;

 

  Event* allocate_event(const Point_2& pt)
  {
    Event *e =  m_eventAlloc.allocate(1); 
    m_eventAlloc.construct(e, m_masterEvent);
    e->init(pt);
    return e;
  }

  public:
  
  void deallocate_event(Event* event)
  {
    m_eventAlloc.destroy(event);
    m_eventAlloc.deallocate(event,1);
  }

};

template <class SweepLineTraits_2, class SweepEvent, class CurveWrap, 
	  class SweepVisitor, typename Allocator>
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
		  Allocator>::
~Sweep_line_2_impl() 
{  
  for(unsigned int i=0 ; i < m_num_of_subCurves; ++i)
    m_subCurveAlloc.destroy(m_subCurves+i);

  if(m_num_of_subCurves) //if its zero, nothing to deallocate
    m_subCurveAlloc.deallocate(m_subCurves,m_num_of_subCurves);

  for(SubCurveListIter itr = m_overlap_subCurves.begin();
      itr != m_overlap_subCurves.end();
      ++itr)
  {
    m_subCurveAlloc.destroy(*itr);
    m_subCurveAlloc.deallocate(*itr, 1);
  }
  delete m_queue;
}

/*! Given an x-monotone curve, create events for each end (if 
 *  one doesn't exist already). 
 *  For each curve create a Subcurve instance.
 */
template < class SweepLineTraits_2,
          class SweepEvent, class CurveWrap,class SweepVisitor,
          typename Allocator>
inline void 
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
_init_curve(X_monotone_curve_2 &curve,unsigned int j)
{ 
  Event *e ;
 
  m_subCurveAlloc.construct(m_subCurves+j,m_masterSubcurve);
  (m_subCurves+j)->init(curve);
 
  const Point_2 &left_end = (m_subCurves+j)->get_left_end();
  const Point_2 &right_end = (m_subCurves+j)->get_right_end();
  
  // Handle the right endpoint of the curve.
  const std::pair<EventQueueIter, bool>& insertion_res =
    m_queue->insert(EventQueueValueType(right_end,0));
  
  if (insertion_res.second == true)
  {
    // Create a new event:
    e = allocate_event(right_end);
    
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    
    (insertion_res.first)->second = e;
  }
  else
  {
    // The event already exists:
    e = (insertion_res.first)->second;
  }
  
  e->add_curve_to_left(m_subCurves+j);
  PRINT_NEW_EVENT(right_end, e);

  (m_subCurves+j)->set_right_event(e);

  // Handle the left endpoint of the curve.
  const std::pair<EventQueueIter, bool>& insertion_res2 =
    m_queue->insert(EventQueueValueType(left_end,0));
  
  if(insertion_res2.second == true)
  {
    // Create a new event:
    e = allocate_event(left_end);
#ifndef NDEBUG
    e->id = m_eventId++;
#endif
    
    (insertion_res2.first)->second = e;  
  }
  else
  {
    // The event already exists:
    e = (insertion_res2.first)->second;
  }
  _add_curve_to_right(e, m_subCurves+j);
  PRINT_NEW_EVENT(left_end, e);
  
  (m_subCurves+j)->set_left_event(e);
  
  m_visitor->init_subcurve(m_subCurves+j);
  
  return;
}

/*!
 * When a curve is removed from the status line for good, its top and
 * bottom neighbors become neighbors. This method finds these cases and
 * looks for the intersection point, if one exists.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return an iterator to the position where the curve will be removed from.
 */
template <class SweepLineTraits_2,
          class SweepEvent, class CurveWrap, class SweepVisitor,
          typename Allocator>
inline void
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
_remove_curve_from_status_line(Subcurve *leftCurve, bool remove_for_good)
                              
{
  SL_DEBUG(std::cout << "remove_curve_from_status_line\n";);
  SL_DEBUG(PrintStatusLine(););
  SL_DEBUG(leftCurve->Print(););

  StatusLineIter sliter = leftCurve->get_hint(); 
  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint;

  if(! remove_for_good)
  {
    m_statusLine.remove_at(sliter);
    SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
    return;
  }

  CGAL_assertion(sliter!=m_statusLine.end());
  if (sliter != m_statusLine.begin() && sliter != m_statusLine.maximum()) 
  {
    StatusLineIter prev = sliter; --prev;
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with  *prev 
    _intersect(*prev, *next, true);
  }
  m_statusLine.remove_at(sliter);
  SL_DEBUG(std::cout << "remove_curve_from_status_line Done\n";)
} 

/*! 
 * Finds intersection between two curves. 
 * If the two curves intersect, create a new event (or use the event that 
 * already exits in the intersection point) and insert the curves to the
 * event.
 * @param curve1 a pointer to the first curve
 * @param curve2 a pointer to the second curve
 * @param after_remove a boolean indicates if the intersection inovoked 
 *        by removing some curve from the Y-str     
 * @return true if the two curves overlap.
*/

template < class SweepLineTraits_2,
           class SweepEvent, class CurveWrap,class SweepVisitor,
           typename Allocator >
inline void 
Sweep_line_2_impl<SweepLineTraits_2,SweepEvent,CurveWrap,SweepVisitor,
                  Allocator>::
_intersect(Subcurve *c1, Subcurve *c2, bool after_remove)
{
  SL_DEBUG(std::cout << "Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  SL_DEBUG(std::cout << "\t";)
  SL_DEBUG(c2->Print();)
  SL_DEBUG(std::cout << "\n";)

  CGAL_assertion(c1 != c2);
  
  // look up for (c1,c2) in the table
  CurvesPair curves_pair(c1,c2);
  if(m_curves_table.find(curves_pair) != m_curves_table.end())
    return;  //the curves have already been checked for intersection

  //insert curves_pair to the table to avoid future checkings for intersection
  m_curves_table.insert(curves_pair);


  vector_inserter vi (m_x_objects) ;
  vector_inserter vi_end (m_x_objects);
  vi_end = m_traits.intersect_2_object()(c1->get_last_curve(),
                                         c2->get_last_curve(),
                                         vi);
 
  if(vi == vi_end) 
  {
    SL_DEBUG(std::cout<<"no intersection...\n";);
    return; // no intersection at all
  }
  
  // BZBZ
  // the two subCurves may start at the same point, in that case we will
  // ignore the first intersection point (if we got to that stage, they cannot 
  // be overlap )
  if((SweepEvent*)c1->get_left_event() == m_currentEvent &&
     (SweepEvent*)c2->get_left_event() == m_currentEvent)
  {
     ++vi;
  }


  if(after_remove)
  {
    for( ; vi != vi_end ; ++vi)
    {
      std::pair<Point_2,unsigned int> xp_point;
      CGAL::assign(xp_point, *vi); // TODO: add an assertion
      if(m_traits.compare_xy_2_object()(m_currentEvent->get_point(),
                                      xp_point.first) ==  SMALLER)
        break;
    }
  }

  for( ; vi != vi_end ; ++vi)
  {
    Object                          xp_object = *vi;
    std::pair<Point_2,unsigned int> xp_point;
    Point_2                         xp;
    unsigned int                    multiplicity = 0;

    if(!CGAL::assign(xp_point, xp_object))
    {
      X_monotone_curve_2 cv;
      // TODO : add some assetions to make sure its X_monotone_curve_2
      CGAL::assign(cv, xp_object);
      xp = m_traits.construct_max_vertex_2_object()(cv);
    }
    else
    {
      xp = xp_point.first;
      multiplicity = xp_point.second;
      SL_DEBUG(std::cout<<"found an intersection point: " << xp << "\n";);
    }

  
    // insert the event and check if an event at this point already exists
   
    const std::pair<EventQueueIter,bool>& insert_res = 
      (m_queue->insert(EventQueueValueType(xp,0)));

    Event *e ;
    if(insert_res.second)    
    {                                   
      // a new event is creatd , which inidicates 
      // that the intersection point cannot be one 
      //of the end-points of two curves
      e = allocate_event(xp);
      
#ifndef NDEBUG
      e->id = m_eventId++;
#endif
      
      e->push_back_curve_to_left(c1);
      e->push_back_curve_to_left(c2);
      
      //multiplicity = 0;  //TODO : remove !!!
      // Act according to the multiplicity:
      if (multiplicity == 0)
      {
	// The multiplicity of the intersection point is unkown or undefined:
        _add_curve_to_right(e, c1);
        _add_curve_to_right(e, c2);
      }
      else
      {
        if((multiplicity % 2) == 1)
	{
	  // The mutiplicity of the intersection point is odd: Swap their
	  // order to the right of this point.
    std::swap(c1,c2);
          e->add_pair_curves_to_right(c1,c2);
        }
	else
        {
	  // The mutiplicity of the intersection point is even, so they
	  // maintain their order to the right of this point.
          CGAL_assertion((multiplicity % 2) == 0);
          e->add_pair_curves_to_right(c1,c2);
        }
      }
	
      PRINT_NEW_EVENT(xp, e);
      (insert_res.first)->second = e;
    } 
    else   // the event already exists
    {
      SL_DEBUG(std::cout << "event already exists,updating.. (" << xp <<")\n";)
      e = (insert_res.first)->second;
    
      e->add_curve_to_left(c1);

      if ( !c1->is_end_point(e))
      { 
        _add_curve_to_right(e, c1);
      }

      if(e->get_num_left_curves() == 1)
        e->push_back_curve_to_left(c2);
      else
        e->add_curve_to_left(c2); 
      
      if ( !c2->is_end_point(e) ) 
      {
        _add_curve_to_right(e, c2);
      }
      SL_DEBUG(e->Print();)
    }
  } 
}

//DEBUG UTILITIES
#ifndef NDEBUG
  #include<CGAL/Sweep_line_2/Sweep_line_2_impl_debug.h>
#endif

CGAL_END_NAMESPACE

#endif
