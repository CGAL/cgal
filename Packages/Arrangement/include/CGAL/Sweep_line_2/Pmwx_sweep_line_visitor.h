#ifndef PMWX_SWEEP_LINE_VISITOR_H
#define PMWX_SWEEP_LINE_VISITOR_H

#include <CGAL/Sweep_line_2/Pmwx_sweep_line_event.h>
#include <CGAL/Sweep_line_2/Pmwx_sweep_line_curve.h>
#include <CGAL/Sweep_line_2/Pmwx_insert_info.h>

CGAL_BEGIN_NAMESPACE

template <class Traits, class Arr, class Arr_notif>
class Pmwx_sweep_line_visitor
{
  typedef typename Arr::Halfedge_handle                            Halfedge_handle;
  typedef Pmwx_sweep_line_visitor<Traits,
                                  Arr,
                                  Arr_notif>                       Self;
  typedef Pmwx_sweep_line_curve<Traits,Self, Halfedge_handle>      Subcurve;
  typedef Pmwx_sweep_line_event<Traits, Subcurve, Self>            Event;
  typedef typename Traits::X_monotone_curve_2                      X_monotone_curve_2;
  typedef typename Traits::Point_2                                 Point_2;

  
  typedef Sweep_line_2_impl<Traits,
                            Event,
                            Subcurve,
                            Self,
                            CGAL_ALLOCATOR(int)>                   Sweep_line;

  typedef Pmwx_insert_info<Halfedge_handle>                        PmwxInsertInfo;
  typedef std::list<Subcurve *>                                    SubcurveContainer;
  typedef typename SubcurveContainer::iterator                     SubCurveIter;
public:

  Pmwx_sweep_line_visitor(Arr *arr, Arr_notif *notif, Traits* traits):
      m_arr(arr),
      m_notif(notif),
      m_traits(traits)
  {}

  void attach(Sweep_line *sl)
  {
    m_sweep_line = sl;
  }

 

  void before_handle_event(Event* event)
  {
    m_currentEvent = event;
  }

  void after_handle_event(Event* event)
  {
    for(SubCurveIter itr = event->left_curves_begin();
        itr != event->left_curves_end();
        ++itr)
    {
      (*itr)->set_last_event(event);
    }

    if(event->get_num_right_curves() == 0)
      m_sweep_line->deallocate_event(event);
    else
    {
      for(SubCurveIter itr = event->right_curves_begin();
        itr != event->right_curves_end();
        ++itr)
      {
        (*itr)->set_last_event(event);
      }
    }
  }

  void add_subcurve(const X_monotone_curve_2& cv,Subcurve* sc)
  {
    Event *lastEvent = (sc)->get_last_event();
    PmwxInsertInfo *insertInfo = lastEvent->get_insert_info();
    Halfedge_handle res; 
    PmwxInsertInfo *currentInfo = m_currentEvent -> get_insert_info();
    Halfedge_handle hhandle = currentInfo->get_halfedge_handle();

    int jump = lastEvent->get_halfedge_jump_count(sc);
    Point_2 p1 ( lastEvent->get_point() );
    Point_2 p2 ( m_currentEvent->get_point() );
    // if the previous event on the curve is not in the planar map yet
    if ( insertInfo->get_halfedge_handle() == Halfedge_handle(NULL) ) 
    {
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) )
      {
        if ( !m_traits->point_equal( hhandle->target()->point(), p2 ))
          hhandle = hhandle->twin();
        res = m_arr->non_intersecting_insert_from_vertex(cv, hhandle, m_notif);
      }
      else
      {
        // if this is the first left curve being inserted
        res = m_arr->insert_in_face_interior(cv, m_arr->unbounded_face(), m_notif);
      }
    } 
    else 
      // the previous event on the curve is already in the planar map. 
      // Let's use it.
    {
      Halfedge_handle prev = insertInfo->get_halfedge_handle();
      if ( !m_traits->point_equal( prev->target()->point(), p1 ))
        prev = prev->twin();
      // skip to the right halfedge
      for ( int i = 0 ; i < jump ; i++ )
        prev = (prev->next_halfedge())->twin();
      
      // we have a handle from the previous insert
      if ( hhandle != Halfedge_handle(NULL) ) 
      {
        if ( !m_traits->point_equal( hhandle->target()->point(), p2 ))
          hhandle = hhandle->twin();

        CGAL_assertion(prev->face() == hhandle->face());
       
         
        res = m_arr->non_intersecting_insert_at_vertices(cv, prev, hhandle, 
          m_notif);
      }
      else
      {
        res = m_arr->non_intersecting_insert_from_vertex(cv, prev, m_notif);                          
      }
    }
    if ( lastEvent->get_num_left_curves() == 0 &&  
      lastEvent->is_curve_largest((Subcurve*)sc) )
    {
      insertInfo->set_halfedge_handle(res->twin());
    }
    insertInfo = m_currentEvent->get_insert_info();
    insertInfo->set_halfedge_handle(res);

    if(lastEvent->get_insert_info()->dec_right_curves_counter() == 0)
    {
      m_sweep_line->deallocate_event(lastEvent);
    }
  }


  void init_subcurve(Subcurve* sc)
  {
    sc -> set_last_event((Event*)(sc->get_left_event()));
  }


   protected:
     
     
  Arr         *m_arr;
  Arr_notif   *m_notif;
  Traits      *m_traits;
  Sweep_line  *m_sweep_line;
  Event       *m_currentEvent;
};



CGAL_END_NAMESPACE

#endif




