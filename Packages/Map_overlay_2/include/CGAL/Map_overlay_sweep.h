// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Map_overlay_sweep.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_MAP_OVERLAY_SWEEP_H
#define CGAL_MAP_OVERLAY_SWEEP_H

#include <vector>
#include <list>

#ifndef CGAL_IN_PLACE_LIST_H
#include <CGAL/In_place_list.h>
#endif

#ifndef CGAL_HANDLE_H
#include <CGAL/Handle.h>
#endif

#ifndef CGAL_ASSERTIONS_H
#include <CGAL/assertions.h>
#endif

#ifndef CGAL_MAP_OVERLAY_BASE_H
#include <CGAL/Map_overlay_base.h>
#endif

#ifndef CGAL_SWEEP_CURVES_BASE_H
#include <CGAL/Sweep_line_2/Sweep_curves_base_2.h>
#endif

#ifndef CGAL_POINT_PLUS_HANDLE_H
#include <CGAL/Sweep_line_2/Point_plus_handle.h>
#endif


CGAL_BEGIN_NAMESPACE

template <class PM_>
class X_curve_plus_id_handle;

template <class PM_>
class Point_plus_handle;


struct Map_overlay_sweep_utils { 
  // X_curve_plus_id_handle:
  // holds a curve and its id number. 
  // The addition of id number to a curve was made due to overlapping 
  // (in which some binary predicates return EQUAL, 
  // while we are interseted in sharp order relation.
  template <class PM_> 
  class X_curve_plus_id_handle : public PM_::Traits::X_curve
  {
  public:
    typedef PM_                           PM;
    typedef typename PM::Traits           Traits;
    //typedef SweepLineTraits_2             Traits;
    typedef typename Traits::X_curve      curve; 
    typedef typename Traits::Point        Point;
    
    //typedef typename Arrangement::Halfedge_handle Halfedge_handle;
    typedef typename PM::Halfedge_const_handle Halfedge_const_handle;
    
    X_curve_plus_id_handle() : 
      curve(), parent(), first_map_(true), flipped_(false) {};
    
    X_curve_plus_id_handle(const curve &cv, Halfedge_const_handle h, 
                           bool first_map, bool flipped, unsigned int id = 0) :
      curve(cv), parent(h), first_map_(first_map), flipped_(flipped), id_(id){}
    
    //Arr_X_curve_plus(Halfedge_handle h, bool b) : curve(h->curve()), parent(h), first_map(b) {}
    
    X_curve_plus_id_handle(Halfedge_const_handle h, 
                           bool  first_map, bool flipped, 
                           unsigned int id = 0) : 
      curve(h->curve()), parent(h), 
      first_map_(first_map), flipped_(flipped), id_(id) {}
    
    // used when no Halfedge_handle is supplied.
    X_curve_plus_id_handle(const curve &cv, bool first_map, 
                           bool flipped, unsigned int id = 0) : 
      curve(cv), parent(),  first_map_(first_map), 
      flipped_(flipped), id_(id) {};
    
    X_curve_plus_id_handle(const X_curve_plus_id_handle &cv) : 
      curve(cv), parent(cv.parent), first_map_(cv.first_map()), 
      flipped_(cv.flipped()), id_(cv.id()) {}
    
    ~X_curve_plus_id_handle(){}
    
    X_curve_plus_id_handle& operator=(const X_curve_plus_id_handle &cv)
    {
      curve::operator=(cv);
      parent=cv.get_parent();
      first_map_ = cv.first_map();
      flipped_= cv.flipped();
      id_ = cv.id();
      return *this;
    }
  
    bool operator==(const X_curve_plus_id_handle &cv) const
    {
      //return curve::operator==(cv);
      Traits traits;
      
      return (id_ == cv.id() && traits.curve_is_same(*this, cv));
    }
    
    void  set_id(unsigned int i) { id_ = i; } 
    
    Halfedge_const_handle get_parent() const
    {
      return parent;
    }
  
    bool first_map() const { return first_map_; }
    
    bool flipped() const { return flipped_; }

    unsigned int id() const { return id_; }

  protected:
    Halfedge_const_handle parent;
    bool  first_map_;
    bool  flipped_;
    unsigned int id_;
  };
};

template <class PM_, 
          class Map_overlay_change_notification_, 
          class Curve_iterator_ = std::list<Map_overlay_sweep_utils::X_curve_plus_id_handle<PM_> >::iterator >
class Map_overlay_sweep : public Map_overlay_base<PM_, Map_overlay_change_notification_>,
                          public Sweep_curves_base_2<Curve_iterator_,
                                                     typename PM_::Traits,
                                                     Point_plus_handle<PM_>, 
                        Map_overlay_sweep_utils::X_curve_plus_id_handle<PM_> >
{
public:
  typedef PM_                                    PM;
  typedef Map_overlay_change_notification_       Map_overlay_change_notification;
  typedef Curve_iterator_                        Curve_iterator;
  typedef Point_plus_handle<PM>                  Point_plus;
  typedef Map_overlay_sweep_utils::X_curve_plus_id_handle<PM> X_curve_plus;

  typedef typename PM::Vertex_handle             Vertex_handle;
  typedef typename PM::Halfedge_handle           Halfedge_handle;
  typedef typename PM::Face_handle               Face_handle;
  typedef typename PM::Vertex_const_handle       Vertex_const_handle;
  typedef typename PM::Halfedge_const_handle     Halfedge_const_handle;
  typedef typename PM::Face_const_handle         Face_const_handle;
  
  typedef typename  PM::Vertex                    Vertex;
  typedef typename  PM::Halfedge                  Halfedge;
  typedef typename  PM::Face                      Face;
  
  typedef typename  PM::Vertex_iterator           Vertex_iterator; 
  typedef typename  PM::Vertex_const_iterator     Vertex_const_iterator; 
  typedef typename  PM::Halfedge_iterator         Halfedge_iterator; 
  typedef typename  PM::Halfedge_const_iterator   Halfedge_const_iterator; 
  typedef typename  PM::Face_iterator             Face_iterator;
  typedef typename  PM::Face_const_iterator       Face_const_iterator;
  typedef typename  PM::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
 
  typedef typename  PM::Locate_type               Locate_type;
  
  //typedef SweepLineTraits_2                                Traits;
  typedef typename  PM::Traits                    Traits;
  typedef typename  Traits::X_curve               X_curve;
  typedef typename  Traits::Point                 Point;

  typedef Sweep_curves_base_2<Curve_iterator, Traits, Point_plus, X_curve_plus>
                                                                         Base;
  
  typedef typename Base::Curve_node                  Curve_node;
  typedef typename Base::Intersection_point_node     Intersection_point_node;
  typedef typename Base::Points_iterator             Points_iterator;  
  typedef typename Base::Points_const_iterator       Points_const_iterator;
  typedef typename Base::Curve_node_iterator         Curve_node_iterator;
  typedef typename Base::Curve_node_const_iterator   Curve_node_const_iterator;
  typedef typename Base::X_curve_list                X_curve_list;
  typedef typename Base::X_curve_list_iterator       X_curve_list_iterator;
  typedef typename Base::Vertices_points_plus        Vertices_points_plus;
  typedef typename Base::Event_queue                 Event_queue;
  typedef typename Base::Status_line                 Status_line;
  
  typedef typename Vertices_points_plus::value_type  
                                             Vertices_points_plus_value_type; 
  typedef typename Event_queue::value_type           Event_queue_value_type;
  typedef typename Status_line::value_type           Status_line_value_type;
  
  typedef typename Vertices_points_plus::iterator   
                                               Vertices_points_plus_iterator;  
  typedef typename Event_queue::iterator            Event_queue_iterator;  
  typedef typename Status_line::iterator            Status_line_iterator; 
  typedef typename std::list<Curve_node>::iterator  list_Curve_node_iterator;
  
  typedef std::list<X_curve_plus>                   X_curve_plus_list;
  typedef X_curve_plus_list::iterator               X_curve_plus_list_iterator;
  
  //typedef Curve_node::Points_iterator Points_iterator; typedef
  //Curve_node::Points_const_iterator Points_const_iterator; typedef
  //Intersection_point_node::Curve_node_iterator Curve_node_iterator;
  //typedef Intersection_point_node::Curve_node_const_iterator
  //Curve_node_const_iterator;
  //      typedef std::pair<Curve_node, X_curve_plus> Curve_pair;
  //      typedef std::map<Point, Point_plus, less_xy<Point> >
  //      Vertices_points_plus; typedef std::multimap<Point,
  //      Intersection_point_node, less_xy<Point> > Event_queue;
  //      typedef std::multimap<Curve_node, X_curve_plus,
  //      less_yx<Curve_node> > Status_line;
  
  void  map_overlay(const PM &a1, 
                    const PM &a2, 
                    Map_overlay_change_notification *pm_change_notf, 
                    PM &result)
  {
    
    std::list<X_curve_plus>   curves;
    unsigned int id = 0;
    
    // updating curevs to contain all updated curves plus.
    Halfedge_const_iterator h_iter;
    for (h_iter = a1.halfedges_begin(); 
         h_iter != a1.halfedges_end(); ++h_iter, ++h_iter, ++id){
      
      bool b;
      X_curve cv(h_iter->curve());
      if ((b = is_right(traits->curve_source(h_iter->curve()), 
                        traits->curve_target(h_iter->curve())) ))
        cv = traits->curve_flip(h_iter->curve());

      Halfedge_const_handle  h = h_iter;
      
      
      curves.push_back(X_curve_plus(cv, h, true, b, id));
                       //is_right(traits.curve_source(cv), 
                       //                      traits.curve_target(cv))) );
    }
    
    for(h_iter = a2.halfedges_begin(); 
        h_iter != a2.halfedges_end(); ++h_iter, ++h_iter, ++id){
      
      bool b;
      X_curve cv(h_iter->curve());
      if ((b = is_right(traits->curve_source(h_iter->curve()), 
                        traits->curve_target(h_iter->curve())) ))
        cv = traits->curve_flip(h_iter->curve());

      Halfedge_const_handle  h = h_iter;

      curves.push_back(X_curve_plus(cv, h, false, b, id));
                       //is_right(traits.curve_source(cv), 
                       //                     traits.curve_target(cv))) );
    }

    sweep_curves_to_planar_map(curves.begin(), 
                               curves.end(), 
                               pm_change_notf, 
                               result);

    pm_change_notf->update_all_faces(result);
  }
  
private:

  void  sweep_curves_to_planar_map(Curve_iterator curves_begin, 
                                   Curve_iterator curves_end, 
                                   Map_overlay_change_notification *pm_change_notf,
                                   PM &result)
  {
    //Traits                traits;
    typename Base::Less_xy  event_queue_pred(traits);
    Event_queue           event_queue(event_queue_pred);
    typename Base::Less_yx  status_pred(traits);
    Status_line           status(status_pred);
    //Event_queue           event_queue;
    //Status_line           status;

    //int c_sweep_t;
    //c_sweep_t = clock();
    
#ifdef  CGAL_SWEEP_LINE_DEBUG
    unsigned int n = 0;
    for (Curve_iterator cv_iter = curves_begin; 
         cv_iter !=  curves_end; ++cv_iter, ++n);
    cout<<"number of edges on input "<< n <<std::endl;
#endif
    
    /*
      // now adding to the x-monotone container all the curves 
      // in the original subdivision.
      for (X_curve_list_iterator  cv_iter = subdivision_curves.begin(); 
      cv_iter !=  subdivision_curves.end(); ++cv_iter)
      x_monotone_curves.push_back(*cv_iter);*/
    
    typename Base::Less_xy  pred(traits);
    Vertices_points_plus  input_vertices(pred);
    for (Curve_iterator cv_iter = curves_begin; 
         cv_iter != curves_end; ++cv_iter){
      if (input_vertices.find(traits->curve_source(*cv_iter)) == 
          input_vertices.end())  
        input_vertices.insert( Vertices_points_plus_value_type
                               (traits->curve_source(*cv_iter), 
                                Point_plus(traits->curve_source(*cv_iter))) );
      if (input_vertices.find(traits->curve_target(*cv_iter)) == 
          input_vertices.end())  
        input_vertices.insert( Vertices_points_plus_value_type
                               (traits->curve_target(*cv_iter), 
                                Point_plus(traits->curve_target(*cv_iter))) );
    }
    // end of input_vertices construction.
    
    // now creating the Curve_node handles and the event queue.
   
    for(Curve_iterator cv_iter = curves_begin; 
        cv_iter != curves_end; ++cv_iter){
      
      //X_curve cv(*cv_iter);
      //if (is_right(traits.curve_source(*cv_iter), 
      //             traits.curve_target(*cv_iter)) )
      //  cv = traits.curve_flip(*cv_iter);
      
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<< *cv_iter <<std::endl;
#endif
   
      Vertices_points_plus_iterator curr_point_plus = 
        input_vertices.find( traits->curve_source(*cv_iter) );
      //assert(traits.curve_source(cv) ==  curr_point_plus->second.point());
      
      Event_queue_iterator  edge_point = 
        event_queue.find( traits->curve_source(*cv_iter) );
      // defining one cv_node for both source and target event points.  
      //X_curve_plus  cv_plus(cv, id);  // to satisfy BCC.
      Curve_node  cv_node = Curve_node(*cv_iter, 
                                       curr_point_plus->second,traits ); 
      Intersection_point_node  source_point_node = 
        Intersection_point_node(cv_node, curr_point_plus->second,traits );
       
      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != source_point_node.get_point())
        event_queue.insert(Event_queue_value_type
			   (traits->curve_source(*cv_iter), 
			    source_point_node));
      else
        edge_point->second.merge(source_point_node);
      
      
      edge_point = event_queue.find( traits->curve_target(*cv_iter) );
      curr_point_plus = input_vertices.find( traits->curve_target(*cv_iter) );
      //assert(traits.curve_target(cv) ==  curr_point_plus->second.point());

      Intersection_point_node  target_point_node = 
        Intersection_point_node(cv_node, curr_point_plus->second, traits );

      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != target_point_node.get_point())
        event_queue.insert(Event_queue_value_type(traits->curve_target(*cv_iter), 
						  target_point_node));
      else
        edge_point->second.merge(target_point_node);
    }

    // now starting the sweeping.
    unsigned int queue_size = 0;
    bool         event_terminated = true;
    //bool         event_overlap_terminated = true;
    while ( !(event_queue.empty()) ){
      ++queue_size;
      // fetch the next event.
      Event_queue_iterator  event = event_queue.begin();  
      
      const Point&              event_point = event->first;
      Intersection_point_node&  point_node = event->second;
      //bool                     event_terminated = true;
      
#ifdef  CGAL_SWEEP_LINE_DEBUG     
      cout<<"* * * event point is "<<event_point<<
        " and point node is "<<point_node.get_point().point()<<std::endl;
      CGAL_assertion(event_point == point_node.get_point().point());
#endif
      
      event_terminated = true; // reinitializing event_terminated to true only after the updating of the subdivision.
     
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, status, 
                                           event_point, point_node);
      
      // handling overlapping curves. 
      // On each overlapping group, we remove iteratively each curve and check for new events after the removing.
      // when finish, we reinsert to the status all the overlappting removed curves.

      handle_overlapping_curves(event_queue, status, event_point, point_node);
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"Printing status line "<<std::endl;
      print_status(status);   
#endif
      
      for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
           cv_iter != point_node.curves_end(); ++cv_iter){
        if (event_point != traits->curve_source(cv_iter->get_curve()) && 
            event_point == cv_iter->get_rightmost_point().point())
          cv_iter->erase_rightmost_point();
      }

      // now, updating the planar map (or arrangement) according the curves 
      // enemating from the currnet event point.
      update_subdivision(point_node, pm_change_notf, result);

      // updating all the new intersection nodes of the curves 
      // participating within the event.
      for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
           cv_iter != point_node.curves_end(); ++cv_iter){
        if (event_point != cv_iter->get_rightmost_point().point())
          cv_iter->push_event_point(point_node.get_point());
      }

      //if (event_terminated)
      event_queue.erase(event);
    }
    
#ifdef  CGAL_SWEEP_LINE_DEBUG  
    std::cout<<"the number of events was "<<queue_size<<std::endl;
#endif

    CGAL_expensive_postcondition_code(is_valid(status)); 
    
    //c_sweep_t = clock() - c_sweep_t;
    //std::cout<<"The time required by sweep proccess: "<< (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;
  }
  
//    void  sweep_curves_to_planar_map(Curve_iterator curves_begin, 
//                                     Curve_iterator curves_end, 
//                                     Map_overlay_change_notification *pm_change_notf, 
//                                     PM &result)
//    {
//      Traits                traits;
//      Event_queue           event_queue;
//      Status_line           status;
//      Vertices_points_plus  input_vertices;
    
//      #ifdef  CGAL_SWEEP_LINE_DEBUG
//      unsigned int n = 0;
//      for (Curve_iterator cv_iter = curves_begin; 
//      cv_iter !=  curves_end; ++cv_iter, ++n);
//      cout<<"number of edges on input "<< n <<std::endl;
//  #endif

//  /*
//      // first handling the case of which results is not empty: since we are sweeping the curves we have to take all 
//      //the curves of result and 'paste' it to the input curves, then we have to clear result.
    
//      X_curve_plus_list  subdivision_curves;
//      for (Halfedge_iterator h_iter = result.halfedges_begin(); 
//      h_iter != result.halfedges_end(); h_iter++, h_iter++)
//      subdivision_curves.push_back(h_iter->curve());
    
//      result.clear();
    
//      // Now, creating all the point_plus handle: for any pair of overlapping points from the input we ensure we have only one handle. - not having such a structure as input_vertices caused a bug.
//      Vertices_points_plus  input_vertices;
//      for (X_curve_list_iterator  cv_iter = subdivision_curves.begin(); 
//      cv_iter != subdivision_curves.end(); cv_iter++){
    
//      if (input_vertices.find(traits.curve_source(*cv_iter)) == 
//      input_vertices.end())  
//      input_vertices.insert(Vertices_points_plus_value_type(traits.curve_source(*cv_iter), 
//      Point_plus(traits.curve_source(*cv_iter))) );
//      if (input_vertices.find(traits.curve_target(*cv_iter)) == 
//      input_vertices.end())  
//      input_vertices.insert(Vertices_points_plus_value_type
//      (traits.curve_target(*cv_iter), 
//      Point_plus(traits.curve_target(*cv_iter))) );
//      }*/
    
//      // splitting all curves to x-monotone curves.
//      X_curve_plus_list  x_monotone_curves;
//      for (Curve_iterator cv_iter = curves_begin; 
//           cv_iter != curves_end; ++cv_iter){
      
//        /*if (!traits.is_x_monotone(*cv_iter)) {
//          X_curve_list x_monotone_subcurves;
//          traits.make_x_monotone(*cv_iter, x_monotone_subcurves);
        
//          #ifdef  CGAL_SWEEP_LINE_DEBUG
//          std::cout<<"printing x-monotone parts"<<std::endl;
//          #endif
//          for(X_curve_list_iterator iter = x_monotone_subcurves.begin(); 
//          iter != x_monotone_subcurves.end(); iter++){
//          #ifdef  CGAL_SWEEP_LINE_DEBUG
//          std::cout<<*iter<<endl;
//          #endif
//          x_monotone_curves.push_back(X_curve_plus(*iter, cv_iter->parent(), cv_iter->first_map(), cv_iter->id()));  
//          }
//          }
//          else*/
      
//        x_monotone_curves.push_back(*cv_iter);
//      }
 
//      /*
//      // now adding to the x-monotone container all the curves 
//      // in the original subdivision.
//      for (X_curve_plus_list_iterator  cv_iter = subdivision_curves.begin(); 
//      cv_iter != subdivision_curves.end(); cv_iter++)
//      x_monotone_curves.push_back(*cv_iter);*/
    
//      // Now, creating all the point_plus handle: for any pair of overlapping points 
//    // from the input we ensure we have only one handle. - not having such 
//    // a structure as input_vertices caused a bug.
//      for (X_curve_plus_list_iterator cv_iter = x_monotone_curves.begin(); 
//           cv_iter != x_monotone_curves.end(); ++cv_iter){
//        if (input_vertices.find(traits.curve_source(*cv_iter)) == 
//            input_vertices.end())  
//          input_vertices.insert( Vertices_points_plus_value_type
//                                 (traits.curve_source(*cv_iter), 
//                                  Point_plus(traits.curve_source(*cv_iter))) );
//        if (input_vertices.find(traits.curve_target(*cv_iter)) == 
//            input_vertices.end())  
//          input_vertices.insert( Vertices_points_plus_value_type
//                                 (traits.curve_target(*cv_iter), 
//                                  Point_plus(traits.curve_target(*cv_iter))) );
//      }
//      // end of input_vertices construction.
    
//      // now creating the Curve_node handles and the event queue.
//      unsigned int id = 0;
//      for(X_curve_plus_list_iterator cv_iter = x_monotone_curves.begin(); 
//          cv_iter != x_monotone_curves.end(); ++cv_iter, ++id){
      
//        X_curve cv(*cv_iter);
//        //Halfedge_const_handle parent = cv_iter->get_parent();
      
//        if (is_right(traits.curve_source(*cv_iter), 
//                     traits.curve_target(*cv_iter)) )
//          cv = traits.curve_flip(*cv_iter);
      
//  #ifdef  CGAL_SWEEP_LINE_DEBUG
//        cout<<cv<<std::endl;
//  #endif
   
//        Vertices_points_plus_iterator curr_point_plus = 
//          input_vertices.find( traits.curve_source(cv) );
//        //assert(traits.curve_source(cv) ==  curr_point_plus->second.point());
      
//        Event_queue_iterator  edge_point = 
//          event_queue.find( traits.curve_source(cv) );
//        // defining one cv_node for both source and target event points.  
//        //X_curve_plus  cv_plus(cv, id);  // to satisfy BCC.
//        Curve_node  cv_node = Curve_node(X_curve_plus(cv, cv_iter->get_parent(), 
//                                                      cv_iter->first_map(), 
//                                                      cv_iter->flipped(), 
//                                                      id), 
//                                         curr_point_plus->second); 
      
//        Intersection_point_node  source_point_node = 
//          Intersection_point_node(cv_node, curr_point_plus->second );
       
//        if (edge_point == event_queue.end() || 
//            edge_point->second.get_point() != source_point_node.get_point())
//          event_queue.insert(Event_queue_value_type
//  			   (traits.curve_source(cv), 
//  			    source_point_node));
//        else
//          edge_point->second.merge(source_point_node);
      
      
//        edge_point = event_queue.find( traits.curve_target(cv) );
//        curr_point_plus = input_vertices.find( traits.curve_target(cv) );
//        //assert(traits.curve_target(cv) ==  curr_point_plus->second.point());

//        Intersection_point_node  target_point_node = 
//          Intersection_point_node(cv_node, curr_point_plus->second );

//        if (edge_point == event_queue.end() || 
//            edge_point->second.get_point() != target_point_node.get_point())
//          event_queue.insert(Event_queue_value_type(traits.curve_target(cv), 
//  						  target_point_node));
//        else
//          edge_point->second.merge(target_point_node);
//      }

//      int c_sweep_t;
//      c_sweep_t = clock();

//      // now starting the sweeping.
//      unsigned int queue_size = 0;
//      bool         event_terminated = true;
//      //bool         event_overlap_terminated = true;
//      while ( !(event_queue.empty()) ){
//        queue_size++;
//        // fetch the next event.
//        Event_queue_iterator  event = event_queue.begin();  
      
//        const Point&              event_point = event->first;
//        Intersection_point_node&  point_node = event->second;
//        //bool                     event_terminated = true;
      
//  #ifdef  CGAL_SWEEP_LINE_DEBUG     
//        cout<<"* * * event point is "<<event_point<<
//          " and point node is "<<point_node.get_point().point()<<std::endl;
//        CGAL_assertion(event_point == point_node.get_point().point());
//  #endif
      
//        event_terminated = true; // reinitializing event_terminated to true only after the updating of the subdivision.
     
//        // now continue with the sweep line.
//        event_terminated = handle_one_event (event_queue, status, 
//                                             event_point, point_node);
        
//        // handling overlapping curves. 
//        // On each overlapping group, we remove iteratively each curve and 
//        // check for new events after the removing.
//        // when finish, we reinsert to the status all the overlappting removed curves.

//        handle_overlapping_curves(event_queue, status, event_point, point_node);
      
//        if (!event_terminated){
//          handle_one_event (event_queue, status, event_point, point_node);
//        }
      
//  #ifdef  CGAL_SWEEP_LINE_DEBUG  
//        cout<<"Printing status line "<<std::endl;
//        print_status(status);   
//  #endif
      
//        for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
//             cv_iter != point_node.curves_end(); ++cv_iter){
//          if (event_point != traits.curve_source(cv_iter->get_curve()) && 
//              event_point == cv_iter->get_rightmost_point().point())
//            cv_iter->erase_rightmost_point();
//        }

//        // now, updating the planar map (or arrangement) according the curves 
//        // enemating from the currnet event point.
      
//        update_subdivision(point_node, pm_change_notf, result);

//        // updating all the new intersection nodes of the curves 
//        // participating within the event.
//        for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
//             cv_iter != point_node.curves_end(); ++cv_iter){
//          if (event_point != cv_iter->get_rightmost_point().point())
//            cv_iter->push_event_point(point_node.get_point());
//        }

//        //if (event_terminated)
//        event_queue.erase(event);
//      }
    
//  #ifdef  CGAL_SWEEP_LINE_DEBUG  
//      std::cout<<"the number of events was "<<queue_size<<std::endl;
//  #endif

//      CGAL_expensive_postcondition_code(is_valid(status)); 
    
//      c_sweep_t = clock() - c_sweep_t;
//      std::cout<<"The time required by sweep proccess: "<< (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;
//    }

  
  /*  bool  check_status_neighbors_intersections(Event_queue &event_queue, Status_line& status,  Status_line::iterator lower_neighbor, Point& point)
  {

#ifdef  CGAL_SWEEP_LINE_DEBUG  
    cout<<"Printing status line "<<std::endl;
    print_status(status);   
#endif

    Traits traits;
    //for (Status_line::iterator status_iter = status.begin(); status_iter != status.end(); status_iter++)
    Curve_node cv1 = lower_neighbor->first;
    
    Status_line::iterator next_neighbor = ++lower_neighbor;
    if (next_neighbor == status.end())
      return false;
    
    lower_neighbor--;
    
    Curve_node  cv2 = next_neighbor->first;
    
#ifdef  CGAL_SWEEP_LINE_DEBUG  
    cout<<"cv1 and cv2 are "<<cv1.get_curve()<<" "<<cv2.get_curve()<<std::endl;
#endif

    // in each node - checking intersections between two adjacent curves and updating the event queue if needed.  
    Point   xp1, xp2; 
    Point ref_point(lower_neighbor->first.get_rightmost_point().point());
    if (is_left(next_neighbor->first.get_rightmost_point().point(), lower_neighbor->first.get_rightmost_point().point()))
      ref_point = next_neighbor->first.get_rightmost_point().point();
    
    if (traits.nearest_intersection_to_right (cv1.get_curve(), cv2.get_curve(), ref_point, xp1, xp2) ){  
      
      // we choose the leftmost point because the intersection point may be one of the points on lower_neighbor or next_neighbor, specially it can be a tangent point, and the function nearest_intersection_to_right may return false (it's third parameter is the ittersection point itself!).

#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"rightmost pointon status, xp1 and xp2 are "<<lower_neighbor->first.get_rightmost_point().point()<<" "<<xp1<<" "<<xp2<<std::endl;
#endif

      // have to handle overlapping.
      if (xp1 == lower_neighbor->first.get_rightmost_point().point())
        xp1 = xp2;
      
      // for debugging.
      //if (traits.curve_get_point_status(cv1.get_curve(), xp1) != Traits::ON_CURVE)
      //  cout<<"The point "<<xp1<<" is not on the curve "<<cv1.get_curve()<<std::endl;
      //if (traits.curve_get_point_status(cv2.get_curve(), xp1) != Traits::ON_CURVE)
      //  cout<<"The point "<<xp1<<" is not on the curve "<<cv2.get_curve()<<std::endl;
      // end debugging.
      
      // if cv1 and cv2 have common edge point - we do not consider it as an intersection point. 
      if ( !( (xp1 == traits.curve_source(cv1.get_curve()) || xp1 == traits.curve_target(cv1.get_curve()) )
              && 
              (xp1 == traits.curve_source(cv2.get_curve()) || xp1 == traits.curve_target(cv2.get_curve()) )) ){
        
        //status_iter->first.push_event_point(xp1);
        //next_iter->first.push_event_point(xp1);
        
        Event_queue::iterator  xp_event = event_queue.find(xp1);
        bool xp_cv1_in_queue = false,  xp_cv2_in_queue = false;
        if (xp_event == event_queue.end())
          event_queue.insert(Event_queue::value_type(xp1, Intersection_point_node(cv1, cv2, xp1)));
        else{
          // have to check the event is a new event. (we might calculated with point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (traits.curve_is_same(cv_iter->get_curve(), cv1.get_curve()) ) // fix later : change it to compare the curve nodes it self!
              xp_cv1_in_queue = true;
            if (traits.curve_is_same(cv_iter->get_curve(), cv2.get_curve()) )
              xp_cv2_in_queue = true;
          }
          
          if (!xp_cv1_in_queue && !xp_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, cv2, Point_plus(xp1)));          
          else if (!xp_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, Point_plus(xp1)));
          else if (!xp_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, Point_plus(xp1)));
        }
        
        // updating the curve_node holding the right point of cv1 (so when we get to this event point 
        //we can find easily the curve node and remove it from status.
        
        point = xp1;
        return (!xp_cv1_in_queue || !xp_cv2_in_queue);
      }
    }
    return false;
    }*/

  void  update_subdivision(Intersection_point_node& point_node, 
                           Map_overlay_change_notification *pm_change_notf, 
                           PM &pm)
  {
    
#ifdef  CGAL_SWEEP_LINE_DEBUG
    cout<<"--------- updating map with point node"<<
      point_node.get_point().point() <<std::endl;
    for (Curve_node_iterator cv_iter1= point_node.curves_begin(); 
         cv_iter1 != point_node.curves_end(); cv_iter1++){
      cv_iter1++;
      Curve_node_iterator cv_iter2 = cv_iter1;
      cv_iter1--;
      for ( ; cv_iter2 != point_node.curves_end(); cv_iter2++){   
        if (traits->curves_overlap(cv_iter1->get_curve(), 
                                   cv_iter2->get_curve()))
          cout<<"update_subdivision "<<cv_iter1->get_curve()<<
            " and "<< cv_iter2->get_curve() <<" are overlapping"<<endl;
      }
    }
#endif    

    X_curve prev_sub_cv;
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); ++cv_iter){
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"now handling "<<cv_iter->get_curve()<<endl;
#endif
      
      bool  overlap=false;
      
      if (is_left(cv_iter->get_rightmost_point().point(), 
                  point_node.get_point().point())) { 
        // means we have a new sub curve to insert.
        Halfedge_handle h;
        
        // first splitting the curve in the points 
        // cv_iter->get_rightmost_point().point() 
        // and point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), right_cv = cv;

        //cout<<"cv is "<<cv<<endl;
        if (traits->curve_source(cv) != 
            cv_iter->get_rightmost_point().point() &&
            traits->curve_target(cv) != 
            cv_iter->get_rightmost_point().point() ) {
          traits->curve_split(cv, sub_cv, right_cv, 
                             cv_iter->get_rightmost_point().point());
          
          cv = right_cv;
        }

        if (traits->curve_source(cv) != point_node.get_point().point() &&
            traits->curve_target(cv) != point_node.get_point().point())
          traits->curve_split(cv, sub_cv, right_cv, 
                             point_node.get_point().point());
        else
          sub_cv = right_cv;
        
        // for debugging!
        //assert(sub_cv == X_curve(cv_iter->get_rightmost_point().point(), 
        // point_node.get_point().point()));
        //cout<<"sub curve calculates is "<<endl;
        //cout<<sub_cv<<endl;
        //cout<<"And the right curve is supposed to be"<<endl;
        //cout<<X_curve (cv_iter->get_rightmost_point().point(), 
        // point_node.get_point().point())<<endl;

        if (cv_iter != point_node.curves_begin()){
          if (traits->curves_overlap(sub_cv, prev_sub_cv)){
            //cout<<sub_cv<<" and "<< prev_sub_cv<<" are overlapping"<<endl;
            overlap=true;
          }
        }

#ifdef  CGAL_SWEEP_LINE_DEBUG
        cout<<"inserting "<<sub_cv<<endl;
#endif
        if (cv_iter->get_curve().flipped()){
          //cout<<"update subdivision: curve is flipped "<<sub_cv<<endl;
          sub_cv = traits->curve_flip(sub_cv);
          //cout<<sub_cv<<endl;
        }
        
        if (pm_change_notf)
          pm_change_notf->set_curve_attributes(sub_cv,  
                                               cv_iter->get_curve().get_parent(),
                                               cv_iter->get_curve().first_map());
        
        if (overlap){
          // special case of overlapping:
          // We do not insert the overlapped curve. 
          // However, we have to call add_edge of the notifier in order to update attributes 
          // of the current halfedge.
          if (pm_change_notf){
            h = find_halfedge(sub_cv, pm);
            pm_change_notf->add_edge(sub_cv, h, true, true);
          }
        }
        else {
          prev_sub_cv = sub_cv;
          
          if (cv_iter->get_rightmost_point().vertex() != Vertex_handle(NULL)){
            //assert(cv_iter->get_rightmost_point().point() == 
            // cv_iter->get_rightmost_point().vertex()->point());
            
            if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
              //assert(point_node.get_point().point() == 
              // point_node.get_point().vertex()->point());
              
              if (cv_iter->get_curve().flipped()) // opposite orientation.
                h = pm.insert_at_vertices(sub_cv, 
                                          point_node.get_point().vertex(), 
                                          cv_iter->get_rightmost_point().vertex(), 
                                          pm_change_notf);
              else
                h = pm.insert_at_vertices(sub_cv, 
                                          cv_iter->get_rightmost_point().vertex(), 
                                          point_node.get_point().vertex(), 
                                          pm_change_notf);
            }
            else {
              if (cv_iter->get_curve().flipped())
                h = pm.insert_from_vertex (sub_cv, 
                                           cv_iter->get_rightmost_point().vertex(), 
                                           false, 
                                           pm_change_notf);
              else
                h = pm.insert_from_vertex (sub_cv, 
                                           cv_iter->get_rightmost_point().vertex(), 
                                           true, 
                                           pm_change_notf);
            }
          }
          else if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
            //assert(point_node.get_point().point() 
            // == point_node.get_point().vertex()->point());
            if (cv_iter->get_curve().flipped())
              h = pm.insert_from_vertex (sub_cv, 
                                         point_node.get_point().vertex(), 
                                         true, 
                                         pm_change_notf);
            else
              h = pm.insert_from_vertex (sub_cv, 
                                         point_node.get_point().vertex(), 
                                         false, 
                                         pm_change_notf);
          }
          else{
            h = pm.insert_in_face_interior (sub_cv, 
                                            pm.unbounded_face(), 
                                            pm_change_notf);
            
            // the point is that if the curve has no source to start the insertion from, it has to be inserted to the unbounded face, because all the curves to the right of it have not inserted yet, and in that stage of the sweep line, the curve is on the unbounded face - later on it will be updated automatically by the Planar map (Arrangement) insert functions.
          }
        }
        //assert(h->source()->point() == cv_iter->get_rightmost_point().point() || 
        // (h->target()->point() == cv_iter->get_rightmost_point().point()));
        
        // now update the vertex handle of each point.
        //if (cv_iter->get_rightmost_point().vertex() == Vertex_handle(NULL))
        if (!overlap || pm_change_notf){
          if (h->source()->point() == 
              cv_iter->get_rightmost_point().point())
            //cv_iter->set_vertex_of_rightmost_point(h->source());
            cv_iter->get_rightmost_point().set_vertex(h->source());
          else if (h->target()->point() == 
                   cv_iter->get_rightmost_point().point())
              //cv_iter->set_vertex_of_rightmost_point(h->target());
            cv_iter->get_rightmost_point().set_vertex(h->target());
        }
        
        //assert(h->source()->point() == point_node.get_point().point() || (h->target()->point() == point_node.get_point().point()));
        
        //if (point_node.get_point().vertex() == Vertex_handle(NULL))
        if (!overlap || pm_change_notf){
          if (h->source()->point() == 
              point_node.get_point().point())
            point_node.get_point().set_vertex(h->source());
          else if (h->target()->point() == 
                   point_node.get_point().point())
            point_node.get_point().set_vertex(h->target());
        }
      }  
      // else - no new sub curve is inserted to the subdivision.
    }
  }
  
  Halfedge_handle find_halfedge(const X_curve& cv, PM& pm)
  {
    //cout<<"In find_halfedge"<<endl;
    
    Locate_type lt;
    Halfedge_handle h = pm.locate(traits->curve_source(cv),lt);

    //cout<<"cv="<<cv<<endl;
    //cout<<"h->curve()="<<h->curve()<<endl;
    
    //if (h->curve() == cv || h->curve() == traits.curve_flip(cv))
    //  return h;

    if (traits->curve_is_same(h->curve(),cv) || 
        traits->curve_is_same(h->curve(),traits->curve_flip(cv)) )
      return h;
     
    Vertex_handle v;
    if (h->source()->point() == traits->curve_source(cv))
      v = h->source();
    else
      v = h->target();
    
    typename PM::Halfedge_around_vertex_circulator 
      circ = v->incident_halfedges();

    do {
      //cout<<"circ->curve()="<<circ->curve()<<endl;
      if (traits->curve_is_same(circ->curve(),cv) || 
          traits->curve_is_same(circ->curve(),traits->curve_flip(cv)))
        return Halfedge_handle(circ);
      
      //if (circ->curve() == cv || circ->curve() == traits.curve_flip(cv))
      //  return Halfedge_handle(circ);
      
    } while (++circ != v->incident_halfedges());

    CGAL_assertion(0);
    return Halfedge_handle(0);
  }
  
  /*  void  update_subdivision(Intersection_point_node& point_node, 
                           Map_overlay_change_notification *pm_change_notf, 
                           PM &pm)
  {
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); cv_iter != point_node.curves_end(); cv_iter++){
      
      if (is_left(cv_iter->get_rightmost_point().point(), point_node.get_point().point())) { // means we have a new sub curve to insert.
        Halfedge_handle h;
        
        pm_change_notf->set_curve_attributes(X_curve(cv_iter->get_rightmost_point().point(), point_node.get_point().point()),  
                                               cv_iter->get_curve().get_parent(), 
                                               cv_iter->get_curve().get_parent()->twin(), cv_iter->get_curve().first_map());
        
        if (cv_iter->get_rightmost_point().vertex() != Vertex_handle(NULL)){
          //assert(cv_iter->get_rightmost_point().point() == cv_iter->get_rightmost_point().vertex()->point());
          
          if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
            //assert(point_node.get_point().point() == point_node.get_point().vertex()->point());
            
            h = pm.insert_at_vertices(X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()), 
                                       cv_iter->get_rightmost_point().vertex(), 
                                       point_node.get_point().vertex(), 
                                       pm_change_notf);
          }
          else
            h = pm.insert_from_vertex (X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()), 
                                        cv_iter->get_rightmost_point().vertex(), 
                                        true, pm_change_notf);
        }
        else if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
          //assert(point_node.get_point().point() == point_node.get_point().vertex()->point());
          
          h = pm.insert_from_vertex (X_curve( point_node.get_point().point(), cv_iter->get_rightmost_point().point()), 
                                      point_node.get_point().vertex(), 
                                      true, pm_change_notf);
        }
        else{
          h = pm.insert_in_face_interior (X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()),
                                           pm.unbounded_face(), 
                                           pm_change_notf);
          
        // the point is that if the curve has no source to start the insertion from, it has to be inserted to the unbounded face, because all the curves to the right of it have not inserted yet, and in that stage of the sweep line, the curve is on the unbounded face - later on it will be updated automatically by the Planar map (Arrangement) insert functions.
        
        }
        
        //assert(h->source()->point() == cv_iter->get_rightmost_point().point() || (h->target()->point() == cv_iter->get_rightmost_point().point()));
        
        // now update the vertex handle of each point.
        //if (cv_iter->get_rightmost_point().vertex() == Vertex_handle(NULL))
        {
          if (h->source()->point() == cv_iter->get_rightmost_point().point())
            cv_iter->set_vertex_of_rightmost_point(h->source());
          else if (h->target()->point() == cv_iter->get_rightmost_point().point())
            cv_iter->set_vertex_of_rightmost_point(h->target());
        }
        
        //assert(h->source()->point() == point_node.get_point().point() || (h->target()->point() == point_node.get_point().point()));
        
        //if (point_node.get_point().vertex() == Vertex_handle(NULL))
        {
          if (h->source()->point() == point_node.get_point().point())
            point_node.set_vertex(h->source());
          else if (h->target()->point() == point_node.get_point().point())
            point_node.set_vertex(h->target());
        }
      }  
      // else - no new sub curve is inserted to the subdivision.
    }
    }*/
  
  void  update_subdivision_debug(Intersection_point_node& point_node, PM &arr)
  {
    cout<<"--------- updating map"<<std::endl;
    
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); cv_iter != point_node.curves_end(); cv_iter++){
      
      if (is_left(cv_iter->get_rightmost_point().point(), point_node.get_point().point())) { // means we have a new sub curve to insert.
        Halfedge_handle h;
        
        cout<<"handling curves "<<cv_iter->get_curve()<<"with right most point so far "<<cv_iter->get_rightmost_point().point()<<std::endl;
        
        if (cv_iter->get_rightmost_point().vertex() != Vertex_handle(NULL)){
          assert(cv_iter->get_rightmost_point().point() == cv_iter->get_rightmost_point().vertex()->point());

          std::cout<<"Point of cv_iter is on map"<< cv_iter->get_rightmost_point().point() <<std::endl;
          
          if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
            std::cout<<"And Point of point_node is on map"<< point_node.get_point().point() <<std::endl;
            
            assert(point_node.get_point().point() == point_node.get_point().vertex()->point());
            
            h = arr.insert_at_vertices(X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()), 
                                       cv_iter->get_rightmost_point().vertex(), 
                                       point_node.get_point().vertex());
          }
          else
            h = arr.insert_from_vertex (X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()), 
                                        cv_iter->get_rightmost_point().vertex(), true);
        }
        else if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
          
          std::cout<<"Only Point of point_node is on map"<< point_node.get_point().point() <<std::endl;
          
          assert(point_node.get_point().point() == point_node.get_point().vertex()->point());
          
          h = arr.insert_from_vertex (X_curve( point_node.get_point().point(), cv_iter->get_rightmost_point().point()), 
                                      point_node.get_point().vertex(), true);
        }
        else{
          
          // debugging.
          Locate_type  lt;
          Halfedge_handle  tmp_h = arr.locate(cv_iter->get_rightmost_point().point(), lt);
          std::cout<<"curve right most point locate type is " <<lt<<", Point is "<<cv_iter->get_rightmost_point().point()<<std::endl;

          tmp_h = arr.locate(point_node.get_point().point(), lt);
          std::cout<<"point intersection node locate type is " <<lt<<", Point is "<<point_node.get_point().point()<<std::endl;
          
          //assert(lt == Arrangement::UNBOUNDED_FACE);
          // end debugging.
          
          h = arr.insert_in_face_interior (X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point()),
                                           arr.unbounded_face());
        // the point is that if the curve has no source to start the insertion from, it has to be inserted to the unbounded face, because all the curves to the right of it have not inserted yet, and in that stage of the sweep line, the curve is on the unbounded face - later on it will be updated automatically by the Planar map (Arrangement) insert functions.
        
        }
        assert(h->source()->point() == cv_iter->get_rightmost_point().point() || (h->target()->point() == cv_iter->get_rightmost_point().point()));
        
        // now update the vertex handle of each point.
        //if (cv_iter->get_rightmost_point().vertex() == Vertex_handle(NULL))
        //for (Curve_node_iterator curr_cv_iter = point_node.curves_begin(); curr_cv_iter != point_node.curves_end(); curr_cv_iter++)
        {
          if (h->source()->point() == cv_iter->get_rightmost_point().point())
            cv_iter->set_vertex_of_rightmost_point(h->source());
          else if (h->target()->point() == cv_iter->get_rightmost_point().point())
            cv_iter->set_vertex_of_rightmost_point(h->target());
          else
            assert(0);
        }
        
        assert(h->source()->point() == point_node.get_point().point() || (h->target()->point() == point_node.get_point().point()));
        
        //if (point_node.get_point().vertex() == Vertex_handle(NULL))
        {
          if (h->source()->point() == point_node.get_point().point())
            point_node.set_vertex(h->source());
          else if (h->target()->point() == point_node.get_point().point())
            point_node.set_vertex(h->target());
          else
            assert(0);
        }
      }  
      // else - no new sub curve is inserted to the subdivision.
    }

    // checking that the handles were correctly updated.
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); cv_iter != point_node.curves_end(); cv_iter++){
      if (is_left(cv_iter->get_rightmost_point().point(), point_node.get_point().point()))  // means we have a new sub curve to insert.
        assert (cv_iter->get_rightmost_point().vertex() != Vertex_handle(NULL));
    }

    cout<<"--------- finish updating map"<<std::endl;
  }
  
  /*
  void  build_overlay_subdivision(const std::list<Curve_node>& curves, Map_overlay_change_notification *pmwx_change_notf, Arrangement &arr)
  {

#ifdef  CGAL_SWEEP_LINE_DEBUG 
    cout<<"Inserting "<<curves.size()<<" Curve nodes"<<std::endl;
#endif

    for (std::list<Curve_node>::const_iterator cv_iter = curves.begin(); cv_iter != curves.end(); cv_iter++){

#ifdef  CGAL_SWEEP_LINE_DEBUG 
      cout<<"The Curve is "<<cv_iter->get_curve()<<std::endl;
#endif
      
      X_curve cv = cv_iter->get_curve(), left_cv, right_cv =  cv_iter->get_curve();
      for (Curve_node::Points_const_iterator points_iter = cv_iter->points_begin(); points_iter != cv_iter->points_end(); points_iter++){
        // make surve the splitting is not at the edge points.
        if (points_iter == cv_iter->points_begin())  
          continue;
        
        ++points_iter;
        if (points_iter == cv_iter->points_end())
          break;
        points_iter--;
        
#ifdef  CGAL_SWEEP_LINE_DEBUG 
        cout<<"Right sub curve is "<<right_cv<<"the splitting point is"<< *points_iter<<std::endl;
#endif
        // debugging.
        if (traits.curve_get_point_status(cv, *points_iter) != Traits::ON_CURVE){
          cout<<"The point "<<*points_iter<<" is not on the curve "<<cv<<std::endl;
          break;
        }

        traits.curve_split(cv, left_cv, right_cv, *points_iter);
#ifdef  CGAL_SWEEP_LINE_DEBUG 
        cout<<"Sub curve is "<<left_cv<<std::endl;
#endif
        //cout<<"Sub curve is "<<*points_iter<<" "<< *next_point<<std::endl;
        //X_curve  left_cv(*points_iter, *next_point);
        
        pmwx_change_notf->set_curve_attributes(left_cv, 
                                               cv_iter->get_curve().get_parent(), 
                                               cv_iter->get_curve().get_parent()->twin(), cv_iter->get_curve().is_first_map());
        arr.insert(left_cv, pmwx_change_notf);

        cv = right_cv;
      }

      pmwx_change_notf->set_curve_attributes(right_cv, 
                                             cv_iter->get_curve().get_parent(), 
                                             cv_iter->get_curve().get_parent()->twin(), cv_iter->get_curve().is_first_map());
      arr.insert(right_cv, pmwx_change_notf);
    }
  }
  
  void  build_pm_overlay(const std::list<Curve_node>& curves, Arrangement &arr)
  {
#ifdef  CGAL_SWEEP_LINE_DEBUG 
    cout<<"Inserting "<<curves.size()<<" Curve nodes"<<std::endl;
#endif
    
    srand(1);  // measuring the time for constructing planar map: we would like to have consistent results. 
    
    for (std::list<Curve_node>::const_iterator cv_iter = curves.begin(); cv_iter != curves.end(); cv_iter++){

#ifdef  CGAL_SWEEP_LINE_DEBUG 
      cout<<"The Curve is "<<cv_iter->get_curve()<<std::endl;
#endif 

      X_curve  cv = cv_iter->get_curve(), left_cv, right_cv =  cv_iter->get_curve();
      Halfedge_handle  h;
      for (Curve_node::Points_const_iterator points_iter = cv_iter->points_begin(); points_iter != cv_iter->points_end(); points_iter++){
        // make surve the splitting is not at the edge points.
        if (points_iter == cv_iter->points_begin())  
          continue;
       
        
#ifdef  CGAL_SWEEP_LINE_DEBUG 
        cout<<"Right sub curve is "<<right_cv<<"the splitting point is"<< *points_iter<<std::endl;
#endif
        bool last = false;
        ++points_iter;
        if (points_iter == cv_iter->points_end())
          last = true;
        points_iter--;
        
        if (last)
          left_cv = right_cv; // if the point is the last one, it's the target point and we can't make split.
        else
          traits.curve_split(cv, left_cv, right_cv, *points_iter);
        
          #ifdef  CGAL_SWEEP_LINE_DEBUG 
        cout<<"Sub curve is "<<left_cv<<std::endl;
#endif
        
        // now inserting left_cv to the map. We would like to make the insertion efficiently, means we know the source vertex in which it has to be inserted if this is not the first curve that is inserted. Hence we use the functions insert_at_vertices or insert_from_vertex.
        bool  second_intersection_point = false;
        points_iter--;
        if (points_iter == cv_iter->points_begin())
          second_intersection_point = true;
        points_iter++;
        
        // h always holds in its source the left point and in its target the right point.
        if (second_intersection_point)
          h = arr.insert(left_cv);
        else{
          Locate_type lt;
          
          Halfedge_handle located_h = arr.locate(*points_iter, lt);
          if (lt == Arrangement::VERTEX){
            Vertex_handle v;
            
            if (located_h->source()->point() == *points_iter)
              v = located_h->source();
            else
              v =  located_h->target();
            
            h = arr.insert_at_vertices (left_cv, h->target(), v);
          }
          else
            h = arr.insert_from_vertex (left_cv, h->target(), true);
        }
        
        cv = right_cv;
      }
      //arr.insert(right_cv);
    }
    }*/
  
/*  //--------------------------------------------------------- debuging functions ------------------------------------------
  bool  is_valid(const Status_line& status){
    std::list<Curve_node> curve_nodes;
    std::list<X_curve>    subcurves;
    
    for (Status_line::const_iterator iter = status.begin(); iter != status.end(); iter++)
      curve_nodes.push_back(iter->first);

    get_subcurves(curve_nodes, subcurves);
    
    return (!do_intersect_subcurves(subcurves));
  }

  // Output function.
  void  write_pm_overlay(const std::list<Curve_node>& curves)
  {
#ifdef  CGAL_SWEEP_LINE_DEBUG 
    cout<<"Inserting "<<curves.size()<<" Curve nodes"<<std::endl;
#endif
    
    std::ofstream f_curves("subcurves_of_sweep.txt"); 
    unsigned int num_of_subcurves = 0;

    for (std::list<Curve_node>::const_iterator cv_iter = curves.begin(); cv_iter != curves.end(); cv_iter++){

      X_curve cv = cv_iter->get_curve(), left_cv, right_cv =  cv_iter->get_curve();
      for (Curve_node::Points_const_iterator points_iter = cv_iter->points_begin(); points_iter != cv_iter->points_end(); points_iter++){
        // make surve the splitting is not at the edge points.
        
        if (points_iter == cv_iter->points_begin())  
          continue;
        
        num_of_subcurves++;
        
        ++points_iter;
        if (points_iter == cv_iter->points_end()){
          f_curves << right_cv.source().xcoordD() <<" "<< right_cv.source().ycoordD() <<"  "<< 
            right_cv.target().xcoordD() <<" "<< right_cv.target().ycoordD() <<std::endl;
          break;
        }
        points_iter--;

        traits.curve_split(cv, left_cv, right_cv, *points_iter);
        
        f_curves << left_cv.source().xcoordD() <<" "<< left_cv.source().ycoordD() <<"  "<< 
          left_cv.target().xcoordD() <<" "<< left_cv.target().ycoordD() <<std::endl;
        
        cv = right_cv;
      }
    }

    std::cout<<"Total number of subcurves is "<<num_of_subcurves<<std::endl;
  }
  
  // graphic function.
  void  draw_subcurves(const std::list<Curve_node>& curves)
  { 
    double    max_x = 10, min_x = -10, min_y = -10;    
    
    // debuging.
    //cout<<"Inserting "<<curves.size()<<" Curve nodes"<<std::endl;
    
    CGAL::Window_stream W(700, 700, "subcurves of overlay");
    W.init(min_x-1, max_x+1, min_y-1);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.button("finish",2);
    W.display();
    
    for (;;) {
      double  x,y;
      
      int b = W.read_mouse(x,y);
      if (b == 2) 
        break;

      for (std::list<Curve_node>::const_iterator cv_iter = curves.begin(); cv_iter != curves.end(); cv_iter++){
        X_curve cv = cv_iter->get_curve(), left_cv, right_cv =  cv_iter->get_curve();
        for (Curve_node::Points_const_iterator points_iter = cv_iter->points_begin(); points_iter != cv_iter->points_end(); points_iter++){
          
          W<<CGAL::BLUE;
          W << *points_iter;
          // make surve the splitting is not at the edge points.
          if (points_iter == cv_iter->points_begin())  
            continue;
          
          ++points_iter;
          if (points_iter == cv_iter->points_end())
            break;
          points_iter--;
          
          //cout<<"Right sub curve is "<<right_cv<<"the splitting point is"<< *points_iter<<std::endl;
          
          // debugging.
          //if (traits.curve_get_point_status(cv, *points_iter) != Traits::ON_CURVE){
          //  cout<<"The point "<<*points_iter<<" is not on the curve "<<cv<<std::endl;
          //  break;
          // }
          
          traits.curve_split(cv, left_cv, right_cv, *points_iter);
        
          //cout<<"Sub curve is "<<left_cv<<std::endl;
          
          W << CGAL::RED;
          W << left_cv;
          
          cv = right_cv;
        }
        
        W << CGAL::RED;
        W << right_cv;
      }
    }
  }
  
  void  get_subcurves(const std::list<Curve_node>& curves,  std::list<X_curve>& subcurves)
  {
    for (std::list<Curve_node>::const_iterator cv_iter = curves.begin(); cv_iter != curves.end(); cv_iter++){
      X_curve cv = cv_iter->get_curve(), left_cv, right_cv =  cv_iter->get_curve();
      for (Curve_node::Points_const_iterator points_iter = cv_iter->points_begin(); points_iter != cv_iter->points_end(); points_iter++){
        // make surve the splitting is not at the edge points.
        if (points_iter == cv_iter->points_begin())  
          continue;
        
        if (*points_iter == rightmost(traits.curve_source(cv_iter->get_curve()), traits.curve_target(cv_iter->get_curve())) )
          left_cv = right_cv;
        else
          traits.curve_split(cv, left_cv, right_cv, points_iter->point());
        
        subcurves.push_back(left_cv);
        
        cv = right_cv;
      }
      //subcurves.push_back(right_cv);
    }
  }
  
  bool  do_intersect_subcurves(const std::list<X_curve>& subcurves)
  {
    for (std::list<X_curve>::const_iterator scv_iter1 = subcurves.begin(); scv_iter1 != subcurves.end(); scv_iter1++){
      for (std::list<X_curve>::const_iterator  scv_iter2 = scv_iter1; scv_iter2 != subcurves.end(); scv_iter2++){
        if (scv_iter2 ==  scv_iter1)
          continue;
        
        Point  ref_point1, ref_point2, xp1, xp2;

        ref_point1 = leftmost(traits.curve_source(*scv_iter1), traits.curve_target(*scv_iter1));
        ref_point2 = leftmost(traits.curve_source(*scv_iter2), traits.curve_target(*scv_iter2));
        if (ref_point1 != ref_point2)
          ref_point1 = leftmost(ref_point1, ref_point2);
        
        if (traits.nearest_intersection_to_right(*scv_iter1, *scv_iter2, ref_point1, xp1, xp2)){
          
          if (xp1 ==  ref_point1)
            xp1 = xp2;
          
          if ( !( (xp1 == traits.curve_source(*scv_iter1) || xp1 == traits.curve_target(*scv_iter1) ) && 
                  (xp1 == traits.curve_source(*scv_iter2) || xp1 == traits.curve_target(*scv_iter2) )) ){

#ifdef  CGAL_SWEEP_LINE_DEBUG   
            cout<<"The curves "<<*scv_iter1<<" "<<*scv_iter2<<" are intersected in the point "<<xp1<<"\n";
#endif
            return true;
          }
        }
      }
    }
    return false;
  }
    
  bool is_left(const Point &p1, const Point &p2) const 
  { return (CGAL::compare_lexicographically_xy(p1, p2) == SMALLER); }
  bool is_right(const Point &p1, const Point &p2) const 
  { return (CGAL::compare_lexicographically_xy(p1, p2) == LARGER); }
  bool is_lower(const Point &p1, const Point &p2) const 
  { return (CGAL::compare_lexicographically_yx(p1, p2) == SMALLER); }
  bool is_higher(const Point &p1, const Point &p2) const 
  { return (CGAL::compare_lexicographically_yx(p1, p2) == LARGER); }

  bool is_same_x(const Point &p1, const Point &p2) const 
  { return (traits.compare_x(p1, p2) == EQUAL); }
  bool is_same_y(const Point &p1, const Point &p2) const 
  { return (traits.compare_y(p1, p2) == EQUAL); }
  bool is_same(const Point &p1, const Point &p2) const
  {
    return (traits.compare_x(p1, p2) == EQUAL) &&
      (traits.compare_y(p1, p2) == EQUAL);
  }
  const Point& leftmost(const Point &p1, const Point &p2) const
  {
    Comparison_result rx = CGAL::compare_lexicographically_xy(p1, p2);
    if (rx == SMALLER)
      return p1;
    else if (rx == LARGER)
      return p2;
   else
     assert(0);
  }

  const Point& rightmost(const Point &p1, const Point &p2) const
  { 
    Comparison_result rx = CGAL::compare_lexicographically_xy(p1, p2);
    if (rx == SMALLER)
      return p2;
    else if (rx == LARGER)
      return p1;
    else
      assert(0);
  }

  //-------------------------------- debuging function.
  void  print_status(const Status_line& status){

    cout<<"Curves on status are\n"; 
    for (Status_line::const_iterator status_iter = status.begin(); status_iter != status.end(); status_iter++){
      X_curve_plus cv1 = status_iter->second;
      
      cout<<cv1<<" ";
      print_points_on_curves(status_iter->first);
    }
    cout<<"\n";
  }

  void  print_points_on_curves(const Curve_node& cv){
    cout<<"list of points is "<<std::endl;

    for (Points_const_iterator  p_iter = cv.points_begin(); p_iter != cv.points_end(); p_iter++)
      cout<<p_iter->point()<<" ";
    
    cout<<std::endl;
}
Traits                traits; */

};

CGAL_END_NAMESPACE

#endif
