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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 (based on old version by Tali Zvi)

#ifndef CGAL_SWEEP_LINE_2_H
#define CGAL_SWEEP_LINE_2_H


#include <list>
#include <CGAL/Object.h>
#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Arrangement_2/Open_hash.h>


CGAL_BEGIN_NAMESPACE

/*!
  Sweep_line_2 is a class that implements the sweep line algorithm
  based on the algorithm of Bentley and Ottmann.
  It extends the algorithm to support not only segments but general curves
  as well and isolated points.
  The curves are defined by the traits class that is one of the template 
  arguments.

  The algorithm is also extended to support the following degenerate cases:
  - non x-monotone curves
  - vertical segments
  - multiple (more then two) segments intersecting at one point
  - curves beginning and ending on other curves.
  - overlapping curves

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
template < class Traits_,
           class SweepVisitor,
           class CurveWrap = Sweep_line_subcurve<Traits_>,
           class SweepEvent = Sweep_line_event<Traits_, CurveWrap>,
           typename Allocator = CGAL_ALLOCATOR(int) >
class Sweep_line_2 : public Basic_sweep_line_2<Traits_,
                                               SweepVisitor,
                                               CurveWrap,
                                               SweepEvent,
                                               Allocator>                                               
{
public:

  typedef Traits_                                 Traits;
  typedef typename Traits::Point_2                Point_2;
  typedef typename Traits::X_monotone_curve_2     X_monotone_curve_2;

  typedef Basic_sweep_line_2<Traits,
                             SweepVisitor,
                             CurveWrap,
                             SweepEvent,
                             Allocator>           Base;


  typedef SweepEvent                              Event;
  typedef typename Base::EventQueueIter           EventQueueIter;
  typedef typename Event::SubCurveIter            EventCurveIter;

  typedef typename Base::Base_event               Base_event;
  typedef typename Base_event::Attribute          Attribute;

  typedef typename Base::Base_subcurve            Base_subcurve;
  
  typedef CurveWrap                               Subcurve;

  typedef std::list<Subcurve*>                    SubCurveList;
  typedef typename SubCurveList::iterator         SubCurveListIter; 

  typedef typename Base::StatusLineIter           StatusLineIter;

  
  typedef Curves_pair<Subcurve>                   CurvesPair;
  typedef Curves_pair_hash_functor<Subcurve>      CurvesPairHasher;
  typedef Curves_pair_equal_functor<Subcurve>     CurvesPairEqual;
  typedef Open_hash<CurvesPair,
                    CurvesPairHasher,
                    CurvesPairEqual>              CurvesPairSet;

  typedef random_access_input_iterator<std::vector<Object> >
                                                  vector_inserter;



  /*!
   * Constructor.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2 (SweepVisitor* visitor) : Base(visitor),
                                         m_curves_pair_set(0)
  {}


  /*!
   * Constructor.
   * \param traits A pointer to a sweep-line traits object.
   * \param visitor A pointer to a sweep-line visitor object.
   */
  Sweep_line_2(Traits *traits, SweepVisitor* visitor) :
      Base(traits, visitor),
      m_curves_pair_set(0)
  {}

  /*! Destrcutor. */
  virtual ~Sweep_line_2()
  {}

  

  protected:

  /*! Init the data structures for the sweep algoeithm */
  void _init_structures()
  {
    Base::_init_structures();

     // Resize the hash to be O(2*n), where n is the number of input curves.
    m_curves_pair_set.resize (2 * m_num_of_subCurves);
  }


  /*! Compete the sweep (compete data strcures) */
  void _complete_sweep()
  {
    Base::_complete_sweep();

    // We can clean the set of curve pairs.
    m_curves_pair_set.clear();

    for(SubCurveListIter itr = m_overlap_subCurves.begin();
        itr != m_overlap_subCurves.end();
        ++itr)
    {
      m_subCurveAlloc.destroy(*itr);
      m_subCurveAlloc.deallocate(*itr, 1);
    }

    m_overlap_subCurves.clear();
  }



  /*! Handle the subcurve to the left of the current event point. */
  void _handle_left_curves()
  { 
    PRINT("Handling left curve" << std::endl;);

   

    m_is_event_on_above = false;

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
      PRINT(" - handling special case " << std::endl;);
             
      m_statusLineCurveLess.set_is_equal(false);
      m_status_line_insert_hint =
        m_statusLine.lower_bound (m_currentEvent->get_point(), 
				                          m_statusLineCurveLess);
      StatusLineIter temp = m_status_line_insert_hint;

      m_is_event_on_above = m_statusLineCurveLess.is_equal();

      if(m_is_event_on_above)
      {
	      // The current event point starts at the interior of a curve at the 
	      // y-structure (can also indicates overlap).
        if(!m_currentEvent -> has_right_curves())
        {
          // event of isolated point
          if(m_currentEvent -> is_query())
          {
            m_is_event_on_above = true;
            m_visitor->before_handle_event(m_currentEvent);
            return;
          }
          CGAL_assertion(m_currentEvent -> is_action());
          m_currentEvent->set_weak_intersection();
        }  

        Subcurve *sc = static_cast<Subcurve*>(*m_status_line_insert_hint);
        const X_monotone_curve_2&  last_curve = sc->get_last_curve();
        m_currentEvent->set_weak_intersection();
        m_visitor->update_event(m_currentEvent, sc);
        m_currentEvent->add_curve_to_left(sc);
 
        bool       is_overlap = _add_curve_to_right(m_currentEvent, sc);

        m_traits->split_2_object() (last_curve,
				   m_currentEvent->get_point(), 
				   sub_cv1, sub_cv2);

        ++m_status_line_insert_hint; 
        
        if(is_overlap)
        {
          m_statusLine.erase (temp);
          m_visitor->before_handle_event (m_currentEvent);
          m_visitor->add_subcurve (sub_cv1, sc);
          return;
        }
      }
      else // no left curves for sure
      {
        m_visitor->before_handle_event(m_currentEvent);
        return;
      }
    }
    
        

    PRINT("left curves before sorting: "<<"\n";);
    SL_DEBUG(if (m_currentEvent->left_curves_begin() != 
		 m_currentEvent->left_curves_end() )
             {
               m_currentEvent->Print();
             });
    _fix_overlap_subcurves();
    _sort_left_curves();
    m_visitor->before_handle_event(m_currentEvent);

    PRINT("left curves after sorting: "<<"\n";);
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
        //if(leftCurve->get_orig_subcurve1() != NULL)
        //  leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);//XXX
        //if(leftCurve->get_orig_subcurve2() != NULL)
        //  leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);//XXX
      }
      else
      {
        const X_monotone_curve_2 &lastCurve = leftCurve->get_last_curve();
       
        m_traits->split_2_object()(lastCurve,
                                   m_currentEvent->get_point(),
                                   sub_cv1,
                                   sub_cv2);
        m_visitor->add_subcurve(sub_cv1, leftCurve);
        leftCurve->set_last_curve(sub_cv2);
      }
      ++left_iter;

      //remove curve from the status line (also checks intersection 
      //between the neighbouring curves,only if the curve is removed for good)
      _remove_curve_from_status_line(leftCurve, remove_for_good);    
    }
    PRINT( "Handling left curve END" << std::endl;);
      
    return;
  }

  

  

  /*! Handle the subcurve to the left of the current event point. */
  void _handle_right_curves()
  {
    PRINT("Handling right curves (" ;)
    PRINT(m_currentEvent->get_point() << ")\n";)
    
    if(! m_currentEvent->has_right_curves())
      return;


    // Loop over the curves to the right of the status line and handle them:
    // - If we are at the beginning of the curve, we insert it to the status 
    //   line, then we look if it intersects any of its neighbors.
    // - If we are at an intersection point between two curves, we add them
    //   to the status line and attempt to intersect them with their neighbors
    // - We also check to see if the two intersect again to the right of the 
    //   point.

    EventCurveIter currentOne = m_currentEvent->right_curves_begin();
    EventCurveIter rightCurveEnd = m_currentEvent->right_curves_end();

    PRINT_INSERT(*currentOne);

    StatusLineIter slIter = 
      m_statusLine.insert_before(m_status_line_insert_hint, *currentOne);
    ((Subcurve*)(*currentOne))->set_hint(slIter);
  
    SL_DEBUG(PrintStatusLine(););
    if ( slIter != m_statusLine.begin() )
    { 
      //  get the previous curve in the y-str
      StatusLineIter prev = slIter; --prev;
      _intersect(static_cast<Subcurve*>(*prev),
                 static_cast<Subcurve*>(*slIter));
    }
    
    
    EventCurveIter prevOne = currentOne;
    ++currentOne;
    while ( currentOne != rightCurveEnd )
    {
      PRINT_INSERT(*currentOne);
      slIter = m_statusLine.insert_before
	(m_status_line_insert_hint, *currentOne);
      ((Subcurve*)(*currentOne))->set_hint(slIter);
        
      SL_DEBUG(PrintStatusLine(););
    
      //BZBZ
      if(reinterpret_cast<Event*>((*currentOne)->get_left_event()) ==
          m_currentEvent ||
          reinterpret_cast<Event*>((*prevOne)->get_left_event()) ==
          m_currentEvent ) 
          _intersect(*prevOne, *currentOne);
      prevOne = currentOne;
      ++currentOne;
    }        
      
    SL_DEBUG(PrintStatusLine(););

    //the next Subcurve at the Y-str 
    ++slIter;
    if ( slIter != m_statusLine.end() )
      _intersect( static_cast<Subcurve*>(*prevOne),
                  static_cast<Subcurve*>(*slIter));
  }


  /*!
   * Add a subcurve to the right of an event point.
   * \param event The event point.
   * \param curve The subcurve to add.
   * \return (true) if an overlap occured; (false) otherwise.
   */
  bool _add_curve_to_right (Event* event, Subcurve* curve,
                            bool overlap_exist = false)
  {
    std::pair<bool, EventCurveIter> pair_res = 
      event->add_curve_to_right(curve, m_traits);

    if (! pair_res.first)
      // No overlap occurs:
      return (false);

    _handle_overlap(event, curve, pair_res.second, overlap_exist);

    // Inidicate that an overlap has occured:
    return (true);
  }


  /*! Fix overlap Subcurves before handling current event */
  void _fix_overlap_subcurves();

  /* Handle overlap at right insertion to event */
  void _handle_overlap(Event* event, Subcurve* curve, EventCurveIter iter, bool overlap_exist);
  
  /*! Compute intersections between the two given curves. */ 
  void _intersect(Subcurve *c1, Subcurve *c2, bool after_remove = false);

  /*! Remove a curve from the status line. */
  void _remove_curve_from_status_line (Subcurve *leftCurve,
				                               bool remove_for_good);

  void _create_intersection_point(const Point_2& xp,
                                  unsigned int multiplicity,
                                  Subcurve* &c1,
                                  Subcurve* &c2,
                                  bool is_overlap = false);


protected:

  /*! contains all of the new sub-curve creaed by overlap */
  SubCurveList m_overlap_subCurves;

  /*! a lookup table of pairs of Subcurves that have been intersected */
  CurvesPairSet m_curves_pair_set;

  /*! Auxiliary vector to hold the intersection objects */
  std::vector<Object> m_x_objects;

  /*! Auxiliary varibales (for splitting curves). */
  X_monotone_curve_2  sub_cv1;
  X_monotone_curve_2  sub_cv2;

 

  template<class SweepCurve>
  Point_2 get_left_end(SweepCurve* sc) const
  {
    return m_traits->construct_min_vertex_2_object()(sc->get_last_curve());
  }
};




/*!
 * When a curve is removed from the status line for good, its top and
 * bottom neighbors become neighbors. This method finds these cases and
 * looks for the intersection point, if one exists.
 *
 * @param leftCurve a pointer to the curve that is about to be deleted
 * @return an iterator to the position where the curve will be removed from.
 */
template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
inline void Sweep_line_2<Traits_,
                         SweepVisitor,
                         CurveWrap,
                         SweepEvent,
                         Allocator>::
_remove_curve_from_status_line(Subcurve *leftCurve, bool remove_for_good)
                              
{
  PRINT("remove_curve_from_status_line\n";);
  SL_DEBUG(PrintStatusLine(););
  SL_DEBUG(leftCurve->Print(););

  StatusLineIter sliter = leftCurve->get_hint(); 
  m_status_line_insert_hint = sliter; ++m_status_line_insert_hint; 

  if(! remove_for_good)
  {
    m_statusLine.erase(sliter);
    PRINT("remove_curve_from_status_line Done\n";)
    return;
  }

  CGAL_assertion(sliter!=m_statusLine.end());
  StatusLineIter lastOne = m_statusLine.end();
  --lastOne;

  if (sliter != m_statusLine.begin() && sliter != lastOne) 
  {
    StatusLineIter prev = sliter; --prev;
    StatusLineIter next = sliter; ++next;
    
    // intersect *next with  *prev 
    _intersect(static_cast<Subcurve*>(*prev),
               static_cast<Subcurve*>(*next),
               true);
  }
  m_statusLine.erase(sliter);
  PRINT("remove_curve_from_status_line Done\n";)
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

template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
void Sweep_line_2<Traits_,
                  SweepVisitor,
                  CurveWrap,
                  SweepEvent,
                  Allocator>::
 _intersect(Subcurve *c1, Subcurve *c2, bool after_remove)
{
  PRINT("Looking for intersection between:\n\t";)
  SL_DEBUG(c1->Print();)
  PRINT("\t";)
  SL_DEBUG(c2->Print();)
  PRINT("\n";)

  CGAL_assertion(c1 != c2);
  
  // special treatment for overlap subcurves
  if(c1->get_orig_subcurve1() != NULL || c2->get_orig_subcurve1() != NULL)
  {
    std::list<Base_subcurve*> c1_origs;
    std::list<Base_subcurve*> c2_origs;

    c1->get_all_inner_noes(std::back_inserter(c1_origs));
    c2->get_all_inner_noes(std::back_inserter(c2_origs));
    typedef typename std::list<Base_subcurve*>::iterator  Itr;
    for(Itr iter1 = c1_origs.begin(); iter1 != c1_origs.end(); ++iter1)
    {
      for(Itr iter2 = c2_origs.begin(); iter2 != c2_origs.end(); ++iter2)
      {
        CurvesPair curves_pair(static_cast<Subcurve*>(*iter1),
                               static_cast<Subcurve*>(*iter2));
        if(m_curves_pair_set.find(curves_pair) !=
           m_curves_pair_set.end())
          return;
      }
    }
  }

  // look up for (c1,c2) in the table and insert if doesnt exist
  CurvesPair curves_pair(c1,c2);
  if(! (m_curves_pair_set.insert(curves_pair)).second )
    return;  //the curves have already been checked for intersection

  float load_factor = static_cast<float>(m_curves_pair_set.size()) /
                        m_curves_pair_set.bucket_count();
  // after lot of benchemarks, keeping load_factor<=6 is optimal
  if(load_factor > 6.0f)
    m_curves_pair_set.resize(m_curves_pair_set.size() * 6);

 
  vector_inserter vi (m_x_objects) ;
  vector_inserter vi_end (m_x_objects);
  vi_end = m_traits->intersect_2_object()(c1->get_last_curve(),
                                         c2->get_last_curve(),
                                         vi);
 
  if(vi == vi_end) 
  {
    PRINT("no intersection...\n";);
    return; // no intersection at all
  }
  
  // BZBZ
  //  the two subCurves may start at the same point,in that case we will
  // ignore the first intersection point (if we got to that stage, they cannot 
  // be overlap )
  if((SweepEvent*)c1->get_left_event() == m_currentEvent &&
     (SweepEvent*)c2->get_left_event() == m_currentEvent)
  {
     ++vi;
  }

  //BZBZ
  // if the two subcurves have a common right-event, 
  // we can ignore last intersection (re-computing the intersection point
  // can crash the sweep later with inexact number types

  if (reinterpret_cast<SweepEvent*>(c1->get_right_event()) ==
      reinterpret_cast<SweepEvent*>(c2->get_right_event()))
  {
    --vi_end; 
  }  

  const std::pair<Point_2,unsigned int>  *xp_point;

  if (after_remove)
  {
    for( ; vi != vi_end ; ++vi)
    {
      xp_point = object_cast<std::pair<Point_2,unsigned int> > (&(*vi));
      CGAL_assertion (xp_point != NULL);

      if(m_traits->compare_xy_2_object() (m_currentEvent->get_point(),
					  xp_point->first) ==  SMALLER)
        break;
    }
  }

  for( ; vi != vi_end ; ++vi)
  {
    const X_monotone_curve_2 *icv;
    Point_2                   xp;
    unsigned int              multiplicity = 0;

    xp_point = object_cast<std::pair<Point_2,unsigned int> > (&(*vi));
    if (xp_point != NULL)
    {
      xp = xp_point->first;
      multiplicity = xp_point->second;
      PRINT("found an intersection point: " << xp << "\n";);
      _create_intersection_point(xp, multiplicity, c1, c2);
    }
    else
    {
      icv = object_cast<X_monotone_curve_2> (&(*vi));
      CGAL_assertion (icv != NULL);

      Point_2 left_xp = m_traits->construct_min_vertex_2_object()(*icv);
      xp = m_traits->construct_max_vertex_2_object()(*icv);
      
      sub_cv1 = *icv;
      _create_intersection_point(xp, 0 , c1 , c2);
      _create_intersection_point(left_xp, 0 , c1 ,c2, true);
    } 
  }
}



template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
void Sweep_line_2<Traits_,
                  SweepVisitor,
                  CurveWrap,
                  SweepEvent,
                  Allocator>::
_create_intersection_point(const Point_2& xp,
                           unsigned int multiplicity,
                           Subcurve* &c1,
                           Subcurve* &c2,
                           bool is_overlap)
{
   // insert the event and check if an event at this point already exists.   
    const std::pair<Event*, bool>& pair_res = 
      push_event(xp, Base_event::DEFAULT);
    
    Event *e = pair_res.first;
    if(pair_res.second)    
    {                                   
      // a new event is creatd , which inidicates 
      // that the intersection point cannot be one 
      //of the end-points of two curves

      e->set_intersection();
      
      m_visitor ->update_event(e, c1, c2, true);
      e->push_back_curve_to_left(c1);
      e->push_back_curve_to_left(c2);
      
      // Act according to the multiplicity:
      if (multiplicity == 0)
      {
        // The multiplicity of the intersection point is unkown or undefined:
        _add_curve_to_right(e, c1, is_overlap);
        _add_curve_to_right(e, c2, is_overlap);
        if(! is_overlap)
        {
          if(e->is_right_curve_bigger(c1, c2))
            std::swap(c1, c2);
        }
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
    } 
    else   // the event already exists
    {
      PRINT("event already exists,updating.. (" << xp <<")\n";);
      if( e == m_currentEvent)  //BZBZ
      {
        //it can happen when c1 starts at the interior of c2 (or the opposite)
        return;
      }

      e->add_curve_to_left(c1);
      e->add_curve_to_left(c2); 

      if ( !c1->is_end_point(e) && !c2->is_end_point(e))
      {
        _add_curve_to_right(e, c1, is_overlap);
        _add_curve_to_right(e, c2, is_overlap);
        e->set_intersection();
        m_visitor ->update_event(e, c1, c2);
      }
      else
      {
        if(!c1->is_end_point(e) && c2->is_end_point(e))
        {
          _add_curve_to_right(e, c1, is_overlap);
          e->set_weak_intersection();
          m_visitor ->update_event(e, c1);
        }
        else 
          if(c1->is_end_point(e) && !c2->is_end_point(e))
          {
            _add_curve_to_right(e, c2, is_overlap);
            e->set_weak_intersection();
            m_visitor ->update_event(e, c2);
          }
      }
     if(! is_overlap)
     {
       if(e->is_right_curve_bigger(c1, c2))
         std::swap(c1, c2);
     }
   
      SL_DEBUG(e->Print();)
    }
}





template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
void Sweep_line_2<Traits_,
                  SweepVisitor,
                  CurveWrap,
                  SweepEvent,
                  Allocator>::
_fix_overlap_subcurves()
{
  CGAL_assertion(m_currentEvent->has_left_curves());
  EventCurveIter leftCurveIter = m_currentEvent->left_curves_begin();

  //special treatment for Subcuves that store overlaps
  while ( leftCurveIter != m_currentEvent->left_curves_end() )  
  {
    Subcurve *leftCurve = *leftCurveIter;
    if( (*leftCurveIter)->get_overlap_subcurve() != NULL &&
            m_traits->compare_xy_2_object()(
            get_left_end((*leftCurveIter)->get_overlap_subcurve()),
            m_currentEvent->get_point()) == SMALLER)
    {
      leftCurve = (Subcurve*)((*leftCurveIter)->getSubcurve());
      while(leftCurve->get_overlap_subcurve() != NULL)
      {
        leftCurve = (Subcurve*)(leftCurve->getSubcurve());
      }

      m_currentEvent->replace_right_curve((*leftCurveIter), leftCurve);
      *leftCurveIter = leftCurve;
    }
    if((Event*)leftCurve->get_right_event() == m_currentEvent)
    {
      if(leftCurve->get_orig_subcurve1() != NULL)
      {
        leftCurve->get_orig_subcurve1()->set_overlap_subcurve(NULL);
        leftCurve->get_orig_subcurve2()->set_overlap_subcurve(NULL);

        std::list<Base_subcurve*> list_of_sc;
        leftCurve->get_all_leaves(std::back_inserter(list_of_sc));
        for(typename std::list<Base_subcurve*>::iterator itr = list_of_sc.begin();
            itr != list_of_sc.end();
            ++itr)
        {
          Base_subcurve *sc = *itr;
          if((void*)sc->get_right_event() != (void*)m_currentEvent && 
              m_traits->compare_xy_2_object()
                (m_traits->construct_min_vertex_2_object()(sc->get_last_curve()),
                 m_currentEvent->get_point()) == SMALLER)
          {
            m_traits->split_2_object() (sc->get_last_curve(),
                                				m_currentEvent->get_point(),
                                        sub_cv1, sub_cv2);
            sc->set_last_curve(sub_cv2);
            sc->set_overlap_subcurve(NULL);
            _add_curve_to_right(m_currentEvent, (Subcurve*)sc);
           
            m_currentEvent->set_weak_intersection();
            m_visitor ->update_event(m_currentEvent,(Subcurve*) sc);
          }
        }
      }
    }     
    ++leftCurveIter;
  }
}


template <class Traits_,
          class SweepVisitor,
          class CurveWrap,
          class SweepEvent,
          typename Allocator>
void Sweep_line_2<Traits_,
                  SweepVisitor,
                  CurveWrap,
                  SweepEvent,
                  Allocator>::
_handle_overlap(Event* event,
                Subcurve* curve,
                EventCurveIter iter,
                bool overlap_exist)
{
    // An overlap occurs:
    // TODO: take care of polylines in which overlap can happen anywhere
    PRINT("Overlap detected at right insertion...\n";);
    //EventCurveIter iter = pair_res.second;
       
    X_monotone_curve_2 overlap_cv;
    if(overlap_exist)
      overlap_cv = sub_cv1;
    else
    {
      std::vector<Object>  obj_vec; 
      vector_inserter vit(obj_vec);
      vector_inserter vit_end = vit;
      vit_end = m_traits->intersect_2_object()(curve->get_last_curve(),
                                              (*iter)->get_last_curve(),
                                              vit);
    
      //BZBZ 06/07/05
      if(obj_vec.empty())
        return;

      overlap_cv = object_cast<X_monotone_curve_2> (obj_vec.front());
    }
    
    // Get the right end of overlap_cv
    Point_2 end_overlap =
      m_traits->construct_max_vertex_2_object()(overlap_cv);

    //find the event assiciated with end_overlap point (right end point)
    EventQueueIter q_iter = m_queue->find( end_overlap );
    CGAL_assertion(q_iter!=m_queue->end());

    Event* right_end = (*q_iter).second;

    // Alocate a new Subcure for the overlap
    Subcurve *overlap_sc = m_subCurveAlloc.allocate(1);
    m_subCurveAlloc.construct(overlap_sc,m_masterSubcurve);
    overlap_sc->init(overlap_cv, event, right_end );
    m_overlap_subCurves.push_back(overlap_sc);

    PRINT(curve<<" + " <<*iter<<" => " <<overlap_sc<<"\n");
    // Set the two events' attribute to intersection
    event -> set_overlap();
    right_end -> set_overlap();

    // Remove curve, *iter from the left curves of end_overlap event
    right_end->remove_curve_from_left(curve);
    right_end->remove_curve_from_left(*iter);

    // Add overlap_sc to the left curves
    right_end->add_curve_to_left(overlap_sc);

    curve  -> set_overlap_subcurve(overlap_sc);
    (*iter)-> set_overlap_subcurve(overlap_sc);

    overlap_sc -> set_orig_subcurve1(*iter);
    overlap_sc -> set_orig_subcurve2(curve);  

    // Replace current sub-curve (*iter) with the new sub-curve
    (*iter) = overlap_sc;
  }




CGAL_END_NAMESPACE

#endif
