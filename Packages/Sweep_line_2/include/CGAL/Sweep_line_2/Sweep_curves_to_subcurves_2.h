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
// release       : $CGAL_Revision: CGAL-2.3-I-44 $
// release_date  : $CGAL_Date: 2001/03/09 $
//
// file          : include/CGAL/Sweep_line_2/Sweep_curves_to_subcurves_2.h
// package       : arr (1.87)
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

#ifndef CGAL_SWEEP_CURVES_TO_SUBCURVES_2_H
#define CGAL_SWEEP_CURVES_TO_SUBCURVES_2_H

#include <vector>
#include <list>

#include <CGAL/In_place_list.h>
#include <CGAL/memory.h>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>
#include <CGAL/Sweep_line_2/Sweep_curves_base_2.h>

CGAL_BEGIN_NAMESPACE

/*! 
  Holds a curve and its id number. 
  The addition of id number to a curve was made due to overlapping 
  (in which some binary predicates return EQUAL, 
  while we are interseted in sharp order relation.
*/
template <class SweepLineTraits_2> 
class X_curve_plus_id : public SweepLineTraits_2::X_curve
{
public:
  typedef SweepLineTraits_2            Traits;
  typedef typename Traits::X_curve     X_Curve; 
  typedef typename Traits::Point       Point;
  
  X_curve_plus_id() : X_Curve() {};
  X_curve_plus_id(const X_Curve &cv, unsigned int i) : X_Curve(cv) , _id(i) {}
  X_curve_plus_id(const X_curve_plus_id &cv) : X_Curve(cv) , _id(cv.id()) {}
  ~X_curve_plus_id(){}
  
  X_curve_plus_id& operator=(const X_curve_plus_id &cv) {
    X_Curve::operator=(cv);
    _id = cv.id();
    return *this;
  }
  
  bool operator==(const X_curve_plus_id &cv) const {
    Traits traits;
    return (_id == cv.id() && traits.curve_is_same(*this, cv));
  }
  
  void  set_id(unsigned int i) { _id = i; }
  unsigned int id() const { return _id; }
    
 protected:
  unsigned int _id;
};

template <class SweepLineTraits_2>
class Point_handle;

/*!
  A wrapper class to Point from the specified traits class.
  This is the equivalent of the Point_plus_rep class defined
  for the planar map version.
*/
template <class SweepLineTraits_2>
class Point_rep : public Ref_counted
{
public:
  typedef SweepLineTraits_2             Traits;
  typedef typename Traits::Point        Point;
    
  Point_rep() {}
  Point_rep(const Point& p) : p_(p) {}
  ~Point_rep() {}
  
protected:
  
  friend class Point_handle<SweepLineTraits_2>;    

  Point p_;
};

/*! A handle class for the Point_rep class.
 * Handle_for is used instead of Handle, using the CGAL_allocator.
 */
template <class SweepLineTraits_2>
class Point_handle : public Handle_for<Point_rep<SweepLineTraits_2> > {

  typedef Handle_for<Point_rep<SweepLineTraits_2> > Handle_for_Point_rep;
  
public:
  typedef SweepLineTraits_2                Traits;
  typedef typename Traits::Point           Point;
  typedef Point_rep<Traits>                Point_rep_traits;
    
  Point_handle() : Handle_for_Point_rep() {
#ifdef CGAL_MAKE_PROFILING
    std::cout<<"allocating handle for DEFAULT Point_rep_traits"<<endl;
#endif
  }
  
  Point_handle(const Point& p) : Handle_for_Point_rep(Point_rep_traits(p)) {  
#ifdef CGAL_MAKE_PROFILING
    std::cout<<"allocating handle for Point_rep_traits(p)" << p << endl;
#endif
  }
  
  Point_handle(const Point_handle& p) : Handle_for_Point_rep(p) {}
    
  ~Point_handle() {}
  
  Point_handle& operator=(const Point_handle &p) {
    Handle_for_Point_rep::operator=(p);
    return *this;
  }
  
  bool operator==(const Point_handle &p) const {
    return ptr()->p_ == p.point(); 
  }
  
  bool operator!=(const Point_handle &p) const { 
    return  !(operator==(p));
  }
  
  void set_point(const Point& p) { ptr()->p_ = p; }
  const Point& point() const { return ptr()->p_; }
    
};

////////////////////////////////////////////////////////////////////////////
//      SWEEP_CURVES_TO_SUBCURVES_2

/*! 
 * A utility class that supplies several methods that perform opearation
 * related to intersection between curves. The intersection calculations
 * are all implemented using the Sweep Line algorithm.
 * No pre condition applies to the curves.
 *
 * The following operations are supported:
 * - report the intersection points
 * - produce all sub curves created by intersecting the curves
 * - report the intersection points and a list of curves for each point
 * - query whether curves intersect.
 */
template <class CurveInputIterator,  class SweepLineTraits_2>
class Sweep_curves_to_subcurves_2 : 
  public Sweep_curves_base_2<CurveInputIterator, 
                             SweepLineTraits_2, 
                             Point_handle<SweepLineTraits_2>, 
                             X_curve_plus_id<SweepLineTraits_2> > 
{
    
public:
  typedef CurveInputIterator                         Curve_iterator;
  typedef SweepLineTraits_2                          Traits;
  typedef Point_handle<Traits>                       Point_handle_traits;
  typedef X_curve_plus_id<Traits>                    X_curve_plus;
  
  typedef typename  Traits::X_curve                  X_curve;
  typedef typename  Traits::Point                    Point;

  typedef Sweep_curves_base_2<Curve_iterator, Traits, 
    Point_handle_traits, X_curve_plus>               Base;

  typedef typename Base::Point_plus                  Point_plus;
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
 
  typedef typename Event_queue::value_type           Event_queue_value_type;
  typedef typename Status_line::value_type           Status_line_value_type;
  
  
  Sweep_curves_to_subcurves_2() : 
    Base() {}
  
  Sweep_curves_to_subcurves_2(Traits *traits_) : 
    Base(traits_) {} 
  
  ~Sweep_curves_to_subcurves_2() {}
  
  template <class OutpoutIterator>
  void  sweep_curves_to_subcurves(Curve_iterator curves_begin, 
                                  Curve_iterator curves_end, 
                                  OutpoutIterator subcurves, 
                                  bool overlapping = false);
  
  template <class OutpoutIterator>
  void  sweep_curves_to_points(Curve_iterator curves_begin, 
                               Curve_iterator curves_end, 
                               OutpoutIterator points,
                               bool endpoints = true,
                               bool overlapping = false);
  
  bool  sweep_do_curves_intersect(Curve_iterator curves_begin, 
                                  Curve_iterator curves_end);
  
  template <class OutputIterator>
  void  sweep_to_report_intersecting_curves(Curve_iterator curves_begin, 
                                            Curve_iterator curves_end, 
                                            OutputIterator intersecting_curves,
                                            bool endpoints = true);
  
private:

  void Init(Curve_iterator curves_begin, 
            Curve_iterator curves_end,
            Event_queue &event_queue,
            Status_line &status);
    
  // Checking whether the point of point_node is an edge point.
  bool  is_endpoint(const Intersection_point_node& point_node)
  {
    const Point& p = point_node.get_point().point();
    
    for (typename Intersection_point_node::Curve_node_const_iterator 
           cv_iter = point_node.curves_begin();
         cv_iter != point_node.curves_end(); ++cv_iter){
      if (traits->curve_source(cv_iter->get_curve()) == p ||
          traits->curve_target(cv_iter->get_curve()) == p )
        return true;
    }

    return false;
  }

  template <class OutpoutIterator>
  void  update_subcurves(Intersection_point_node& point_node, 
                         OutpoutIterator subcurves, 
                         bool overlapping) 
  {
    if (overlapping)
      update_subcurves_with_overlappings(point_node, subcurves);
    else
      update_subcurves_without_overlappings(point_node, subcurves);
  }
  
  template <class OutpoutIterator>
  void  update_subcurves_with_overlappings(Intersection_point_node& point_node,
                                           OutpoutIterator subcurves);
  template <class OutpoutIterator> 
  void  update_subcurves_without_overlappings(
                    Intersection_point_node& point_node, 
                    OutpoutIterator subcurves);
};


////////////////////////////////////////////////////////////////////////
////            METHODS




/*!
 *  Given a range of curves, this function returns a list of curves
 *  that are created by intersecting the input curves.
 *  The intersections are calculated using the sweep algorithm.
 *  \param curves_begin the input iterator that points to the first curve 
 *                      in the range.
 *  \param curves_end the input past-the-end iterator of the range.
 *  \param subcurves an iterator to the first curve in the range
 *                   of curves created by intersecting the input curves.
 *  \param overlapping indicates whether there are overlapping curves
 *                     in the input range. Defaults to false.
 */
template <class CurveInputIterator, class SweepLineTraits_2>
template <class OutpoutIterator>
inline void 
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
sweep_curves_to_subcurves(CurveInputIterator curves_begin, 
			  CurveInputIterator curves_end, 
			  OutpoutIterator subcurves, 
			  bool overlapping)
{ 
    typename Base::Less_xy         event_queue_pred(traits);
    Event_queue                    event_queue(event_queue_pred);
    typename Base::Less_yx         status_pred(traits);
    Status_line                    status(status_pred);
    
    Init(curves_begin, curves_end, event_queue, status);
    
    // now starting the sweeping.
    bool event_terminated = true;
    while ( !(event_queue.empty()) ){
      // fetch the next event.
      typename Event_queue::iterator  event = event_queue.begin();

      const Point&              event_point = event->first;
      Intersection_point_node&  point_node = event->second;
      
      CGAL_SL_DEBUG(
        cout<<"* * * event point is "<<event_point<<
          " and point node is "<<point_node.get_point().point()<<std::endl;
        CGAL_assertion(event_point == point_node.get_point().point());
      )
      
      event_terminated = true; // reinitializing event_terminated 
      // to true only after the updating of the subdivision.
     
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, 
                                           status, 
                                           event_point, 
                                           point_node);
      
      handle_overlapping_curves(event_queue,status,event_point,point_node);
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
      CGAL_SL_DEBUG(
        cout<<"Printing status line "<<std::endl;
        print_status(status);
      )

      Curve_node_iterator ncv_iter;
      for (ncv_iter = point_node.curves_begin(); 
           ncv_iter != point_node.curves_end(); ++ncv_iter){
        if (event_point != traits->curve_source(ncv_iter->get_curve()) && 
            event_point == ncv_iter->get_rightmost_point().point())
          ncv_iter->erase_rightmost_point();
      }
      // now, updating the planar map (or arrangement) according the 
      // curves enemating from the currnet event point.
      update_subcurves(point_node, subcurves, overlapping);

      // updating all the new intersection nodes of the curves 
      // participating within the event.
      for (ncv_iter = point_node.curves_begin(); 
           ncv_iter != point_node.curves_end(); ++ncv_iter){
        if (event_point != ncv_iter->get_rightmost_point().point())
          ncv_iter->push_event_point(point_node.get_point());
      }

      //if (event_terminated)
      event_queue.erase(event);
    }

    CGAL_expensive_postcondition_code(is_valid(status)); 
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
template <class CurveInputIterator, class SweepLineTraits_2>
template <class OutpoutIterator>
inline void
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
sweep_curves_to_points(Curve_iterator curves_begin, 
		       Curve_iterator curves_end, 
		       OutpoutIterator points,
		       bool endpoints,
		       bool overlapping)
{ 
    typename Base::Less_xy  event_queue_pred(traits);
    Event_queue           event_queue(event_queue_pred);
    typename Base::Less_yx  status_pred(traits);
    Status_line           status(status_pred);

    Init(curves_begin, curves_end, event_queue, status);
    
    // now starting the sweeping.
    bool         event_terminated = true;
    while ( !(event_queue.empty()) ){
      typename Event_queue::iterator  event = event_queue.begin();

      const Point&              event_point = event->first;
      Intersection_point_node&  point_node = event->second;
      //bool                     event_terminated = true;
      
      CGAL_SL_DEBUG(
        cout<<"* * * event point is "<<event_point<<
          " and point node is "<<point_node.get_point().point()<<std::endl;
        CGAL_assertion(event_point == point_node.get_point().point());
      )
      
      event_terminated = true; // reinitializing event_terminated 
      // to true only after the updating of the subdivision.
      
      //int handle_one_event_t = clock();
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, 
                                           status, 
                                           event_point, 
                                           point_node);
      
      handle_overlapping_curves(event_queue,status,event_point,point_node);
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
      CGAL_SL_DEBUG(
        cout<<"Printing status line "<<std::endl;
        print_status(status);
      )

      Curve_node_iterator cvn_iter;
      for (cvn_iter = point_node.curves_begin(); 
           cvn_iter != point_node.curves_end(); ++cvn_iter){
        if (event_point != traits->curve_source(cvn_iter->get_curve()) && 
            event_point == cvn_iter->get_rightmost_point().point())
          cvn_iter->erase_rightmost_point();
      }

      // now, updating the points containter according the 
      // curves enemating from the currnet event point.
      
      if (endpoints || !is_endpoint(point_node))
        *points = point_node.get_point().point();
      ++points;
      
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
    
    CGAL_expensive_postcondition_code(is_valid(status)); 
  }


/*!
 *  Given a range of curves, this function returns true if any two 
 *  curves in the range intersects. Returns false otherwise.
 *  The intersections are calculated using the sweep algorithm.
 *  \param curves_begin the input iterator that points to the first curve 
 *                      in the range.
 *  \param curves_end the input past-the-end iterator of the range.
 */
template <class CurveInputIterator, class SweepLineTraits_2>
inline bool
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
sweep_do_curves_intersect(Curve_iterator curves_begin, 
                          Curve_iterator curves_end)
{ 
  typename Base::Less_xy         pred1(traits);
  Event_queue                    event_queue(pred1);
  typename Base::Less_yx         pred2(traits);
  Status_line                    status(pred2);
  
  Init(curves_begin, curves_end, event_queue, status);
    
  // now starting the sweeping.
  bool         event_terminated = true;
  while ( !(event_queue.empty()) ){
    typename Event_queue::iterator  event = event_queue.begin();

    const Point&              event_point = event->first;
    Intersection_point_node&  point_node = event->second;
      
    event_terminated = true; // reinitializing event_terminated 
    // to true only after the updating of the subdivision.
     
    // now continue with the sweep line.
    event_terminated = handle_one_event (event_queue, status, 
					 event_point, point_node);
      
    handle_overlapping_curves(event_queue,status,event_point,point_node);
      
    if (!event_terminated) {
      handle_one_event (event_queue, status, event_point, point_node);
    }
      
    if (intersection_exist_)
      return intersection_exist_;

    CGAL_SL_DEBUG(  
		  cout<<"Printing status line "<<std::endl;
		  print_status(status);)   

      event_queue.erase(event);
  }

  CGAL_expensive_postcondition_code(is_valid(status)); 

  return intersection_exist_;
}



/*!
 *  Given a range of curves, this function returns an iterator 
 *  to the beginning of a range that contains the list of curves 
 *  for each intersection point between any two curves in the 
 *  specified range.
 *  The intersections are calculated using the sweep algorithm.
 *  \param curves_begin the input iterator that points to the first curve 
 *                      in the range.
 *  \param curves_end the input past-the-end iterator of the range.
 *  \param intersecting_curves an iterator to the output
 *  \param endpoints if true, the end points of the curves are reported
 *                   as intersection points. Defaults to true.
 */
template <class CurveInputIterator, class SweepLineTraits_2>
template <class OutputIterator>
inline void
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
sweep_to_report_intersecting_curves(Curve_iterator curves_begin, 
                                    Curve_iterator curves_end, 
                                    OutputIterator intersecting_curves,
                                    bool endpoints)
  { 
    typename Base::Less_xy  event_queue_pred(traits);
    Event_queue           event_queue(event_queue_pred);
    typename Base::Less_yx  status_pred(traits);
    Status_line           status(status_pred);
    
    Init(curves_begin, curves_end, event_queue, status);
    
    // now starting the sweeping.
    bool         event_terminated = true;
    while ( !(event_queue.empty()) ){
      typename Event_queue::iterator  event = event_queue.begin();

      const Point&              event_point = event->first;
      Intersection_point_node&  point_node = event->second;
      
      CGAL_SL_DEBUG(
        cout<<"* * * event point is "<<event_point<<
          " and point node is "<<point_node.get_point().point()<<std::endl;
        CGAL_assertion(event_point == point_node.get_point().point());
      )
      
      event_terminated = true; // reinitializing event_terminated 
      // to true only after the updating of the subdivision.
      
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, 
                                           status, 
                                           event_point, 
                                           point_node);
      
      handle_overlapping_curves(event_queue,status,event_point,point_node);
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
      if (intersection_exist_ || 
          (endpoints && is_endpoint(point_node)))
        {
          std::pair<Point, std::list<X_curve> > intersection_point_node;
          intersection_point_node.first = point_node.get_point().point();
          
          Curve_node_iterator cv_iter;
          for (cv_iter = point_node.curves_begin(); 
               cv_iter != point_node.curves_end(); ++cv_iter)
            {
              intersection_point_node.second.push_back(cv_iter->get_curve());
            }
          
          *intersecting_curves = intersection_point_node;
          ++intersecting_curves;
        }

      CGAL_SL_DEBUG(
        cout<<"Printing status line "<<std::endl;
        print_status(status);
      )

      Curve_node_iterator cvn_iter;
      for (cvn_iter = point_node.curves_begin(); 
           cvn_iter != point_node.curves_end(); ++cvn_iter){
        if (event_point != traits->curve_source(cvn_iter->get_curve()) && 
            event_point == cvn_iter->get_rightmost_point().point())
          cvn_iter->erase_rightmost_point();
      }
      
      // updating all the new intersection nodes of the curves 
      // participating within the event.
      for (cvn_iter = point_node.curves_begin(); 
           cvn_iter != point_node.curves_end(); ++cvn_iter){
        if (event_point != cvn_iter->get_rightmost_point().point())
          cvn_iter->push_event_point(point_node.get_point());
      }
      
      event_queue.erase(event);
    }
    
    CGAL_expensive_postcondition_code(is_valid(status));
  }

/*!
  Initializes the data structures used in the sweep algorithm, including:
  - splitting all input curves so that they are x-monotone
  - initializing the event queue with the curves end points
*/
template <class CurveInputIterator, class SweepLineTraits_2>
inline void
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
Init(Curve_iterator curves_begin, 
     Curve_iterator curves_end,
     Event_queue &event_queue,
     Status_line &status)
{
  // Remove out of loop scope to compile undef MSVC:
  Curve_iterator cv_iter;
  X_curve_list_iterator xcv_iter;
    
    
  CGAL_SL_DEBUG(
    unsigned int n = 0;
    for (cv_iter = curves_begin; cv_iter !=  curves_end; ++cv_iter, ++n);
    cout<<"number of edges on input "<< n <<std::endl;)
        
    
  // splitting all curves to x-monotone curves.
  X_curve_list x_monotone_curves;
  for (cv_iter = curves_begin; cv_iter != curves_end; ++cv_iter){
    if (!traits->is_x_monotone(*cv_iter)) {
      X_curve_list x_monotone_subcurves;
      traits->make_x_monotone(*cv_iter, x_monotone_subcurves);
      
      CGAL_SL_DEBUG(std::cout<<"printing x-monotone parts"<<std::endl;)
        for (X_curve_list_iterator iter = x_monotone_subcurves.begin(); 
	     iter != x_monotone_subcurves.end(); ++iter)
	  {
	    CGAL_SL_DEBUG(std::cout<<*iter<<endl;)
	      x_monotone_curves.push_back(*iter);  
	  }
    }
    else
      x_monotone_curves.push_back(*cv_iter);
  }
  
  // now creating the Curve_node handles and the event queue.
  unsigned int id = 0;
  for(xcv_iter = x_monotone_curves.begin(); 
      xcv_iter != x_monotone_curves.end(); ++xcv_iter, ++id){
    
    X_curve cv(*xcv_iter);
    if (is_right(traits->curve_source(*xcv_iter), 
		 traits->curve_target(*xcv_iter)) )
      cv = traits->curve_flip(*xcv_iter);
    
    CGAL_SL_DEBUG(cout<<cv<<std::endl;)
      
    typename Event_queue::iterator  edge_point = 
    event_queue.find( traits->curve_source(cv) );
    // defining one cv_node for both source and target event points. 
    Curve_node cv_node = Curve_node(X_curve_plus(cv, id),  // ?? diff
				    traits->curve_source(cv), traits ); 
    
    Intersection_point_node  source_point_node = 
      Intersection_point_node(cv_node, traits->curve_source(cv), traits );
    
    if (edge_point == event_queue.end() || 
	edge_point->second.get_point() != source_point_node.get_point())
      event_queue.insert(Event_queue_value_type(traits->curve_source(cv), 
						source_point_node));
    else
      edge_point->second.merge(source_point_node);
    
    
    edge_point = event_queue.find( traits->curve_target(cv) );
    
    Intersection_point_node  target_point_node = 
      Intersection_point_node(cv_node, traits->curve_target(cv), traits );
    
    if (edge_point == event_queue.end() || 
	edge_point->second.get_point() != target_point_node.get_point())
      event_queue.insert(Event_queue_value_type(traits->curve_target(cv), 
						target_point_node));
    else
      edge_point->second.merge(target_point_node);
  }
}



/*!
  Updates the curves the hae overlapping.
*/
template <class CurveInputIterator, class SweepLineTraits_2>
template <class OutpoutIterator>
inline void
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
update_subcurves_with_overlappings(Intersection_point_node& point_node,
				   OutpoutIterator subcurves) 
{
  CGAL_SL_DEBUG(
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
  )

  X_curve prev_sub_cv;
  for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
       cv_iter != point_node.curves_end(); ++cv_iter){
    CGAL_SL_DEBUG(cout<<"now handling "<<cv_iter->get_curve()<<endl;)
      
      if (is_left(cv_iter->get_rightmost_point().point(), 
                  point_node.get_point().point())) { 
        // means we have a new sub curve to insert.
        
        // first splitting the curve in the points 
        // cv_iter->get_rightmost_point().point() 
        // and point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), 
          right_cv = cv;
	
        if (traits->curve_source(cv) != 
            cv_iter->get_rightmost_point().point() 
            &&
            traits->curve_target(cv) != 
            cv_iter->get_rightmost_point().point()) 
          {
            traits->curve_split(cv, sub_cv, right_cv, 
				cv_iter->get_rightmost_point().point());
	    
            cv = right_cv;
          }
        
        if (traits->curve_source(cv) != point_node.get_point().point() &&
            traits->curve_target(cv) != point_node.get_point().point())
          traits->curve_split(cv, sub_cv, 
			      right_cv, point_node.get_point().point());
        else
          sub_cv = right_cv;
        
        CGAL_SL_DEBUG(cout<<"inserting "<<sub_cv<<endl;)
	  
	*subcurves = sub_cv;
        ++subcurves;
      }
    // else - no new sub curve is inserted to the subdivision.
  }
}

template <class CurveInputIterator, class SweepLineTraits_2>
template <class OutpoutIterator> 
inline void
Sweep_curves_to_subcurves_2<CurveInputIterator, SweepLineTraits_2>::
update_subcurves_without_overlappings(
                    Intersection_point_node& point_node, 
                    OutpoutIterator subcurves)
  {
    CGAL_SL_DEBUG(cout<<"--------- updating map with point node"<<
      point_node.get_point().point() <<std::endl;
    for (Curve_node_iterator cv_iter1= point_node.curves_begin(); 
         cv_iter1 != point_node.curves_end(); ++cv_iter1){
      ++cv_iter1;
      Curve_node_iterator cv_iter2 = cv_iter1;
      --cv_iter1;
      for ( ; cv_iter2 != point_node.curves_end(); ++cv_iter2){   
        if (traits->curves_overlap(cv_iter1->get_curve(), 
                                  cv_iter2->get_curve()))
          cout<<"update_subdivision "<<cv_iter1->get_curve()<<
            " and "<< cv_iter2->get_curve() <<" are overlapping"<<endl;
      }
    })

    X_curve prev_sub_cv;
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); ++cv_iter){
      CGAL_SL_DEBUG(cout<<"now handling "<<cv_iter->get_curve()<<endl;)

      if (is_left(cv_iter->get_rightmost_point().point(), 
                  point_node.get_point().point())) { 
        // means we have a new sub curve to insert.
        
        // first splitting the curve in the points 
        // cv_iter->get_rightmost_point().point() and 
        // point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), 
          right_cv = cv;
        
        if (traits->curve_source(cv) != 
            cv_iter->get_rightmost_point().point() 
            &&
            traits->curve_target(cv) != 
            cv_iter->get_rightmost_point().point() ) 
          {
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
            continue;
          }
        }

        CGAL_SL_DEBUG(cout<<"inserting "<<sub_cv<<endl;)

        prev_sub_cv = sub_cv;

        *subcurves = sub_cv;
        ++subcurves;
      }
      // else - no new sub curve is inserted to the subdivision.
    }
}

CGAL_END_NAMESPACE

#endif
