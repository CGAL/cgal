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
// release       : $CGAL_Revision: CGAL-2.3-I-81 $
// release_date  : $CGAL_Date: 2001/07/10 $
//
// file          : include/CGAL/Sweep_line_2/Sweep_curves_to_planar_map_2.h
// package       : Arrangement (2.07)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_SWEEP_CURVES_TO_PLANAR_MAP_2_H
#define CGAL_SWEEP_CURVES_TO_PLANAR_MAP_2_H

#include <vector>
#include <list>

#include <CGAL/In_place_list.h>
#include <CGAL/Handle.h>
#include <CGAL/assertions.h>
#include <CGAL/Sweep_line_2/Sweep_curves_base_2.h>
#include <CGAL/Sweep_line_2/Point_plus_handle.h>

CGAL_BEGIN_NAMESPACE

template <class PM_>
class Point_plus_handle;

struct Sweep_curves_to_planar_map_utils { 
  // X_curve_plus_id:
  // holds a curve and its id number. 
  // The addition of id number to a curve was made due to overlapping 
  // (in which some binary predicates return EQUAL, 
  // while we are interseted in sharp order relation. 
  template <class SweepLineTraits_2> 
  class X_curve_plus_id : public SweepLineTraits_2::X_curve
  {
  public:
    typedef SweepLineTraits_2             Traits;
    typedef typename Traits::X_curve      curve; 
    typedef typename Traits::Point        Point;
    
    X_curve_plus_id() : curve() {};
    
    X_curve_plus_id(const curve &cv, unsigned int i) : curve(cv) , _id(i) {}
    
    X_curve_plus_id(const X_curve_plus_id &cv) : curve(cv) , _id(cv.id()) {}
    
    ~X_curve_plus_id(){}
    
    X_curve_plus_id& operator=(const X_curve_plus_id &cv)
    {
      curve::operator=(cv);
      _id = cv.id();
      return *this;
    }
    
  bool operator==(const X_curve_plus_id &cv) const
    {
      Traits traits;
      
      return (_id == cv.id() && traits.curve_is_same(*this, cv));
      
      //return curve::operator==(cv);
    }
    
    void  set_id(unsigned int i) { _id = i; }
    
    unsigned int id() const { return _id; }
    
  protected:
    unsigned int _id;
  };
};

template <class CurveInputIterator, class SweepLineTraits_2, 
          class PM_, class Change_notification_>
class Sweep_curves_to_planar_map_2  : 
  public Sweep_curves_base_2<CurveInputIterator, 
                             SweepLineTraits_2, 
                             Point_plus_handle<PM_>, 
  Sweep_curves_to_planar_map_utils::X_curve_plus_id<SweepLineTraits_2> >
{
public:
  typedef CurveInputIterator                     Curve_iterator;
  typedef PM_                                    PM;
  typedef Change_notification_                   Change_notification;
  typedef Point_plus_handle<PM>                  Point_plus;
  typedef Sweep_curves_to_planar_map_utils::X_curve_plus_id<SweepLineTraits_2>
                                                 X_curve_plus;
  
  typedef typename PM::Vertex_handle              Vertex_handle;
  typedef typename PM::Halfedge_handle            Halfedge_handle;
  typedef typename PM::Vertex_const_handle        Vertex_const_handle;
  typedef typename PM::Halfedge_const_handle      Halfedge_const_handle;
  typedef typename  PM::Vertex_iterator           Vertex_iterator; 
  typedef typename  PM::Vertex_const_iterator     Vertex_const_iterator; 
  typedef typename  PM::Halfedge_iterator         Halfedge_iterator; 
  typedef typename  PM::Halfedge_const_iterator   Halfedge_const_iterator; 
 
  typedef SweepLineTraits_2                                Traits;
  typedef typename  Traits::X_curve                        X_curve;
  typedef typename  Traits::Point                          Point;

  typedef Sweep_curves_base_2<Curve_iterator, Traits, 
    Point_plus, X_curve_plus>                        Base;
  
  //typedef typename Base::X_curve_plus                X_curve_plus;
  
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

  Sweep_curves_to_planar_map_2() : 
    Base() {}
  
  Sweep_curves_to_planar_map_2(Traits *traits_) : 
    Base(traits_) {} 
  
  ~Sweep_curves_to_planar_map_2() {}
  
  void  sweep_curves_to_planar_map(Curve_iterator curves_begin, 
                                   Curve_iterator curves_end, 
                                   PM &result,
                                   Change_notification* change_notification)
                                   
  {
    // Remove out of loop scope to compile undef MSVC:
    Curve_iterator cv_iter;
    X_curve_list_iterator xcv_iter;
    
    typename Base::Less_xy  event_queue_pred(traits);
    Event_queue event_queue(event_queue_pred);
    typename Base::Less_yx  status_pred(traits);
    Status_line status(status_pred);
    
    //Event_queue           event_queue;
    //Status_line           status;
    
#ifdef  CGAL_SWEEP_LINE_DEBUG
    unsigned int n = 0;
    for (cv_iter = curves_begin; cv_iter !=  curves_end; ++cv_iter, ++n);
    cout<<"number of edges on input "<< n <<std::endl;
#endif

    // first handling the case of which results is not empty: since we
    // are sweeping the curves we have to take all
    //the curves of result and 'paste' it to the input curves, 
    // then we have to clear result.
    X_curve_list  subdivision_curves;
    for (Halfedge_iterator h_iter = result.halfedges_begin(); 
         h_iter != result.halfedges_end(); ++h_iter, ++h_iter)
      subdivision_curves.push_back(h_iter->curve());
    
    result.clear();
    
    // Now, creating all the point_plus handle: for any pair of
    // overlapping points from the input we ensure we have only one
    // handle. - not having such a structure as input_vertices caused
    // a bug.
    typename Base::Less_xy  pred(traits);
    Vertices_points_plus  input_vertices(pred);
    for (xcv_iter = subdivision_curves.begin(); 
         xcv_iter !=  subdivision_curves.end(); ++xcv_iter){

      if (input_vertices.find(traits->curve_source(*xcv_iter)) == 
          input_vertices.end())  
        input_vertices.insert(Vertices_points_plus_value_type(
                               traits->curve_source(*xcv_iter), 
                               Point_plus(traits->curve_source(*xcv_iter))) );
      if (input_vertices.find(traits->curve_target(*xcv_iter)) == 
          input_vertices.end())  
        input_vertices.insert(Vertices_points_plus_value_type
                              (traits->curve_target(*xcv_iter), 
                               Point_plus(traits->curve_target(*xcv_iter))) );
    }
    
    // splitting all curves to x-monotone curves.
    X_curve_list x_monotone_curves;
    for (cv_iter = curves_begin; cv_iter != curves_end; ++cv_iter){
      if (!traits->is_x_monotone(*cv_iter)) {
        X_curve_list x_monotone_subcurves;
        traits->make_x_monotone(*cv_iter, x_monotone_subcurves);

#ifdef  CGAL_SWEEP_LINE_DEBUG
        std::cout<<"printing x-monotone parts"<<std::endl;
#endif
        for(X_curve_list_iterator iter = x_monotone_subcurves.begin(); 
            iter != x_monotone_subcurves.end(); ++iter){
#ifdef  CGAL_SWEEP_LINE_DEBUG
          std::cout<<*iter<<endl;
#endif
          x_monotone_curves.push_back(*iter);  
        }
      }
      else
        x_monotone_curves.push_back(*cv_iter);
    }
    
    // now adding to the x-monotone container all the curves 
    // in the original subdivision.
    for (xcv_iter = subdivision_curves.begin(); 
         xcv_iter != subdivision_curves.end(); ++xcv_iter)
      x_monotone_curves.push_back(*cv_iter);
    
    for (xcv_iter = x_monotone_curves.begin(); 
         xcv_iter != x_monotone_curves.end(); ++xcv_iter){
      if (input_vertices.find(traits->curve_source(*xcv_iter)) == 
          input_vertices.end())  
        input_vertices.insert( Vertices_points_plus_value_type
                               (traits->curve_source(*xcv_iter), 
                                Point_plus(traits->curve_source(*xcv_iter))) );
      if (input_vertices.find(traits->curve_target(*xcv_iter)) == 
          input_vertices.end())  
        input_vertices.insert( Vertices_points_plus_value_type
                               (traits->curve_target(*xcv_iter), 
                                Point_plus(traits->curve_target(*xcv_iter))) );
    }
    // end of input_vertices construction.
    
    // now creating the Curve_node handles and the event queue.
    unsigned int id = 0;
    for(xcv_iter = x_monotone_curves.begin(); 
        xcv_iter != x_monotone_curves.end(); ++xcv_iter, ++id){
      
      X_curve cv(*xcv_iter);
      if (is_right(traits->curve_source(*xcv_iter), 
                   traits->curve_target(*xcv_iter)) )
        cv = traits->curve_flip(*xcv_iter);
      
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<cv<<std::endl;
#endif
   
      Vertices_points_plus_iterator curr_point_plus = 
        input_vertices.find( traits->curve_source(cv) );
#ifdef  CGAL_SWEEP_LINE_DEBUG
      CGAL_assertion(curr_point_plus != input_vertices.end());
      CGAL_assertion(traits->curve_source(cv) ==
                     curr_point_plus->second.point());
#endif

      Event_queue_iterator  edge_point = 
        event_queue.find( traits->curve_source(cv) );
      // defining one cv_node for both source and target event points.  
      //X_curve_plus  cv_plus(cv, id);  // to satisfy BCC.
      Curve_node  cv_node = Curve_node(X_curve_plus(cv, id), 
                                       curr_point_plus->second, traits); 
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<< "bogi" <<std::endl;
#endif

      Intersection_point_node  source_point_node = 
        Intersection_point_node(cv_node, curr_point_plus->second, traits );
      
      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != source_point_node.get_point())
        event_queue.insert(Event_queue_value_type
                           (traits->curve_source(cv), 
                            source_point_node));
      else
        edge_point->second.merge(source_point_node);
      
      edge_point = event_queue.find( traits->curve_target(cv) );
      curr_point_plus = input_vertices.find( traits->curve_target(cv) );
      //assert(traits->curve_target(cv) ==  curr_point_plus->second.point());

      Intersection_point_node  target_point_node = 
        Intersection_point_node(cv_node, curr_point_plus->second, traits );

      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != target_point_node.get_point())
        event_queue.insert(Event_queue_value_type(traits->curve_target(cv), 
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
      
      event_terminated = true; 
      // reinitializing event_terminated to true only after the
      // updating of the subdivision.
     
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, status, 
                                           event_point, point_node);
      
      handle_overlapping_curves(event_queue,status,event_point,point_node);
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"Printing status line "<<std::endl;
      print_status(status);   
#endif

      Curve_node_iterator cvn_iter;
      for (cvn_iter = point_node.curves_begin(); 
           cvn_iter != point_node.curves_end(); ++cvn_iter){
        if (event_point != traits->curve_source(cvn_iter->get_curve()) && 
            event_point == cvn_iter->get_rightmost_point().point())
          cvn_iter->erase_rightmost_point();
      }

      // now, updating the planar map (or arrangement) according the curves 
      // enemating from the currnet event point.
      update_subdivision(point_node, change_notification, result);

      // updating all the new intersection nodes of the curves 
      // participating within the event.
      for (cvn_iter = point_node.curves_begin(); 
           cvn_iter != point_node.curves_end(); ++cvn_iter){
        if (event_point != cvn_iter->get_rightmost_point().point())
          cvn_iter->push_event_point(point_node.get_point());
      }

      //if (event_terminated)
      event_queue.erase(event);
    }
    
#ifdef  CGAL_SWEEP_LINE_DEBUG  
    std::cout<<"the number of events was "<<queue_size<<std::endl;
#endif

    CGAL_expensive_postcondition_code(is_valid(status)); 
  }
  
private:
  void  update_subdivision(Intersection_point_node& point_node, 
                           Change_notification *pm_change_notf,
                           PM &pm)
  {

#ifdef  CGAL_SWEEP_LINE_DEBUG
    cout<<"--------- updating map with point node"<<
      point_node.get_point().point() <<std::endl;
    for (Curve_node_iterator cv_iter1= point_node.curves_begin(); 
         cv_iter1 != point_node.curves_end(); ++cv_iter1){
      ++cv_iter1;
      Curve_node_iterator cv_iter2 = cv_iter1;
      cv_iter1--;
      for ( ; cv_iter2 != point_node.curves_end(); ++cv_iter2){   
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
        X_curve  cv = cv_iter->get_curve();
        X_curve  sub_cv = cv_iter->get_curve(), right_cv = cv;

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

        if (cv_iter != point_node.curves_begin()){
          if (traits->curves_overlap(sub_cv, prev_sub_cv)){
            //cout<<sub_cv<<" and "<< prev_sub_cv<<" are overlapping"<<endl;
            overlap=true;
          }
        }
        
        if (overlap){
          // special case of overlapping: We do not insert the
          // overlapped curve.  
          // However, we have to call add_edge of
          // the notifier in order to update attributes of the current
          // halfedge.
          if (pm_change_notf){
            h = find_halfedge(sub_cv, pm);
            pm_change_notf->add_edge(sub_cv, h, true, true);
          }
        }
        else {

#ifdef  CGAL_SWEEP_LINE_DEBUG
          cout<<"inserting "<<sub_cv<<endl;
#endif

          prev_sub_cv = sub_cv;
          if (cv_iter->get_rightmost_point().vertex() != Vertex_handle(NULL)){
            //assert(cv_iter->get_rightmost_point().point() == 
            // cv_iter->get_rightmost_point().vertex()->point());
            
            if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
              //assert(point_node.get_point().point() == 
              // point_node.get_point().vertex()->point());
              
              h = pm.insert_at_vertices(
                                    sub_cv, 
                                    cv_iter->get_rightmost_point().vertex(),
                                    point_node.get_point().vertex(), 
                                    pm_change_notf);
            }
            else
              h = pm.insert_from_vertex(
                                     sub_cv, 
                                     cv_iter->get_rightmost_point().vertex(), 
                                     true,
                                     pm_change_notf);
          }
          else if (point_node.get_point().vertex() != Vertex_handle(NULL)) {
            //assert(point_node.get_point().point() 
            // == point_node.get_point().vertex()->point());
            
            h = pm.insert_from_vertex(sub_cv, 
                                      point_node.get_point().vertex(), 
                                      false, 
                                      pm_change_notf);
          }
          else{
            h = pm.insert_in_face_interior(sub_cv,  
                                           pm.unbounded_face(),
                                           pm_change_notf);
            
            // the point is that if the curve has no source to start
            // the insertion from, it has to be inserted to the
            // unbounded face, because all the curves to the right of
            // it have not inserted yet, and in that stage of the
            // sweep line, the curve is on the unbounded face - later
            // on it will be updated automatically by the Planar map
            // (Arrangement) insert functions.
            
          }
        }
        
        //assert(h->source()->point() ==
        //cv_iter->get_rightmost_point().point() ||
        //(h->target()->point() ==
        //cv_iter->get_rightmost_point().point()));
        
        // now update the vertex handle of each point.
        //if (cv_iter->get_rightmost_point().vertex() == Vertex_handle(NULL))
        if (!overlap || pm_change_notf) {
          if (h->source()->point() == cv_iter->get_rightmost_point().point())
            //cv_iter->set_vertex_of_rightmost_point(h->source());
            cv_iter->get_rightmost_point().set_vertex(h->source());
          else if (h->target()->point() == 
                   cv_iter->get_rightmost_point().point())
            //cv_iter->set_vertex_of_rightmost_point(h->target());
            cv_iter->get_rightmost_point().set_vertex(h->target());
        }
        
        //assert(h->source()->point() ==
        //point_node.get_point().point() || (h->target()->point() ==
        //point_node.get_point().point()));
        
        //if (point_node.get_point().vertex() == Vertex_handle(NULL))
        if (!overlap || pm_change_notf) {
          if (h->source()->point() == point_node.get_point().point())
            point_node.get_point().set_vertex(h->source());
          else if (h->target()->point() == point_node.get_point().point())
            point_node.get_point().set_vertex(h->target());
        }
      }  
      // else - no new sub curve is inserted to the subdivision.
    }
  }

  Halfedge_handle find_halfedge(const X_curve& cv, PM& pm)
  { 
    typename PM::Locate_type lt;
    Halfedge_handle h = pm.locate(traits->curve_source(cv),lt);
    
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
      if (traits->curve_is_same(circ->curve(),cv) || 
          traits->curve_is_same(circ->curve(),traits->curve_flip(cv)))
        return Halfedge_handle(circ);
      
    } while (++circ != v->incident_halfedges());

    CGAL_assertion(0);
    return Halfedge_handle(0);
  }
  
};

CGAL_END_NAMESPACE

#endif
