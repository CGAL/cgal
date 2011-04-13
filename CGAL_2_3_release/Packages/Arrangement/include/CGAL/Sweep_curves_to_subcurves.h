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
// file          : include/CGAL/Sweep_curves_to_subcurves.h
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

#ifndef CGAL_SWEEP_CURVES_TO_SUBCURVES_H
#define CGAL_SWEEP_CURVES_TO_SUBCURVES_H

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

#ifndef CGAL_SWEEP_CURVES_BASE_H
#include <CGAL/Sweep_curves_base.h>
#endif

//#include <CGAL/IO/leda_window.h>  //used for visualization -

CGAL_BEGIN_NAMESPACE

// Point_rep:
// Point_rep holds only a Point. This class holds the representation, 
// and the next will hold the Handle to Point.
template <class Traits_>
class Point_handle;

template <class Traits_>
class Point_rep : public Rep {
public:
  typedef Traits_                       Traits;
  typedef typename Traits::Point        Point;
  

  Point_rep() {}
  
  Point_rep(const Point& p) : p_(p) {}
  
  ~Point_rep() {}
  
protected:
  friend class Point_handle<Traits>;    
  
  Point p_;
};

// Point_handle:
// The handle to Point.
template <class Traits_>
class Point_handle : public Handle {
public:
  typedef Traits_                          Traits;
  typedef typename Traits::Point           Point;
  typedef Point_rep<Traits>                Point_rep_traits;

  Point_handle() : Handle() {}
  
  Point_handle(const Point& p) {  PTR = new Point_rep_traits(p); }
  
  Point_handle(const Point_handle& p) : Handle(p) {}
  
  ~Point_handle() {}
  
  Point_handle& operator=(const Point_handle &p) {
    //Point::operator=(p);
      ptr()->p_ = p.point();
      return *this;
  }
  
  bool operator==(const Point_handle &p) const
  { return ptr()->p_ == p.point(); }

    bool operator!=(const Point_handle &p) const
    { return ptr()->p_ != p.point(); }

    void set_point(const Point& p) { ptr()->p_ = p; }
    
    const Point& point() const { return ptr()->p_; }
    
  private:
    Point_rep_traits* ptr() const { return (Point_rep_traits*) PTR; }
  };

template <class Curve_iterator_,  class Traits_, class Container_>
class Sweep_curves_to_subcurves : 
  public Sweep_curves_base<Curve_iterator_, Traits_, Point_handle<Traits_> > 
{
  /*
  // Arr_X_curve_plus:
  // holds a curve and its id number. 
  // The addition of id number to a curve was made due to overlapping 
  // (in which some binary predicates return EQUAL, 
  // while we are interseted in sharp order relation. 
  class Arr_X_curve_plus: public Traits_::X_curve
  {
  public:
    typedef Traits_                       Traits;
    typedef typename Traits::X_curve      curve; 
    typedef typename Traits::Point        Point;

    Arr_X_curve_plus() : curve() {};

    Arr_X_curve_plus(const curve &cv, unsigned int i) : curve(cv) , _id(i) {}
    
    Arr_X_curve_plus(const Arr_X_curve_plus &cv) : curve(cv) , _id(cv.id()) {}
    
    ~Arr_X_curve_plus(){}
    
    Arr_X_curve_plus& operator=(const Arr_X_curve_plus &cv)
    {
      curve::operator=(cv);
      _id = cv.id();
      return *this;
    }
    
    bool operator==(const Arr_X_curve_plus &cv) const
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
  
  // Point_plus_rep:
  // Point_plus_rep holds a Point plus a vertex handle of the vertex 
  // in the subdivision that will hold that point.
  // The reason we need the vertex handle information is to update 
  // the subdivision by the time the sweep line progresses without makeing any 
  // point location query. This class holds the representation, 
  // and the next will hold the Handle to Point_plus.
  class Point_plus_rep : public Rep {
  public:
    typedef Traits_                       Traits;
    typedef typename Traits::Point        Point;
    //typedef typename PM::Vertex_handle   Vertex_handle;

    Point_plus_rep() {}

    Point_plus_rep(const Point& p) : p_(p) {}
    
    ~Point_plus_rep() {}
    
 protected:
    friend class Point_plus;    
    
    Point p_;
    //Vertex_handle v_;
  };
  
  // Point_plus:
  // The handle to Point_plus.
  class Point_plus : public Handle {
  public:
    typedef Traits_                          Traits;
    typedef typename Traits::Point           Point;
    //typedef typename PM::Vertex_handle  Vertex_handle;

    Point_plus() : Handle() {}
    
    Point_plus(const Point& p) {  PTR = new Point_plus_rep(p); }

    Point_plus(const Point_plus& p_plus) : Handle(p_plus) {}

    ~Point_plus() {}
    
    Point_plus& operator=(const Point_plus &p_plus) {
      //Point::operator=(p);
      ptr()->p_ = p_plus.point();
      return *this;
    }
    
    bool operator==(const Point_plus &p_plus) const
    { return ptr()->p_ == p_plus.point(); }

    bool operator!=(const Point_plus &p_plus) const
    { return ptr()->p_ != p_plus.point(); }

    void set_point(const Point& p) { ptr()->p_ = p; }
    
    const Point& point() const { return ptr()->p_; }
    
  private:
    Point_plus_rep* ptr() const { return (Point_plus_rep*) PTR; }
  };

  // Curve_node_rep:
  // Curve_node_rep holds a curve participating in the sweep proccess. Curve_node_rep holds a curve and a container of points.
  // The points container refers all the intersedction points (including edge points) calculated so far by the sweep line.
  // These points are ordered from left to right on the curve, which means they are sorted in a way we get immidiately all 
  // the disjoint subcurves reduce by the curve.
  class Curve_node_rep : public Rep {
  public:
    typedef Traits_                              Traits;
    typedef typename Traits::X_curve             X_curve; 
    typedef typename Traits::Point               Point;
    typedef  Arr_X_curve_plus                    X_curve_plus;
    typedef std::list<Point_plus>                Points_container;

    Curve_node_rep(const X_curve_plus& cv) : cv_(cv) {
      Traits  traits;
      points.push_back ( Point_plus(leftmost(traits.curve_source(cv), traits.curve_target(cv))) );
    }

    Curve_node_rep(const X_curve_plus& cv, const Point& p) : cv_(cv) {
      points.push_back(Point_plus(p));
    }

    Curve_node_rep(const X_curve_plus& cv, const Point_plus& p) : cv_(cv) {
      points.push_back(p);
    }
    
    //Curve_node_rep(const Self& cv_node_rep) : 
    // cv_(cv_node_rep.get_curve()), points(cv_node_rep.points_begin(), cv_node_rep.points_end()) {}

    ~Curve_node_rep() {}

  protected:
    friend class Curve_node;

    X_curve_plus          cv_;  // hold the left most intersecting point.
    Points_container      points;

  private:
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
  };
  
  // Curve_node:
  // The handle to curve node.
  class Curve_node : public Handle {
  public:
    typedef Traits_                                         Traits;
    typedef typename Traits::X_curve                        X_curve; 
    typedef typename Traits::Point                          Point;
    typedef  Arr_X_curve_plus                               X_curve_plus;
    typedef  Curve_node                                     Self;

    typedef Curve_node_rep::Points_container                Points_container;
    typedef Points_container::iterator                      Points_iterator;
    typedef Points_container::const_iterator                Points_const_iterator;

    Curve_node() : Handle() {}
    
    Curve_node(const X_curve_plus& cv) {
      PTR = new Curve_node_rep(cv);
    }
    
    // Curve_node(const X_curve_plus& cv, Vertex_handle v) {
    //  PTR = new Curve_node_rep(cv, v);
    // }

    Curve_node(const X_curve_plus& cv, const Point& p) {
      PTR = new Curve_node_rep(cv, p);
    }

    Curve_node(const X_curve_plus& cv, const Point_plus& p) {
      PTR = new Curve_node_rep(cv, p);
    }
    
    Curve_node(const Curve_node& cv_node) : Handle(cv_node) {}

    ~Curve_node() {}

    void  push_event_point(const Point_plus &event_point) {
      ptr()->points.push_back(event_point);
    }
    
    void  erase_rightmost_point() { 
      Points_iterator  iter = ptr()->points.end();
      iter--;
      ptr()->points.erase(iter); 
    } 

    const X_curve_plus& get_curve() const { 
      return  ptr()->cv_; 
    }
    
    const Point_plus& get_rightmost_point() const { 
      assert ( ptr()->points.size() > 0);
      return *( ptr()->points.rbegin());
    }

    Points_iterator  points_begin() { return  ptr()->points.begin(); }
    Points_iterator  points_end() { return  ptr()->points.end(); }
    
    Points_const_iterator  points_begin() const { return  ptr()->points.begin(); }
    Points_const_iterator  points_end() const { return  ptr()->points.end(); }

    // fix this operator!
    Self& operator=(const Self &cv_node)
    {
      ptr()->cv_ = cv_node.get_curve();
      ptr()->points.assign(cv_node.points_begin(), cv_node.points_end());
      
      return *this;
    }
    
  private:
    Curve_node_rep* ptr() const {return (Curve_node_rep*)PTR;}
  };
  
  // Intersection_point_node:
  // Intersection_point_node holds an event point. This class contains a Point_plus refering the event point, 
  // and a container holding all the Curve_nodes  anamating from that point. the curves are order counter clock wise around the event point.
  // the ordering is done trivially in the constructors, and specifically in the merge function which merges two Intersection_point_node sharing 
  // the same event point.
  class Intersection_point_node {
  public:
    typedef  Traits_                                        Traits;
    typedef  typename Traits::X_curve                        X_curve; 
    typedef  typename Traits::Point                          Point;
    typedef  Arr_X_curve_plus                                X_curve_plus;
    typedef  Intersection_point_node                         Self;

    typedef  std::list<Curve_node>                           Curve_node_container;
    typedef  Curve_node_container::iterator                  Curve_node_iterator;
    typedef  Curve_node_container::const_iterator            Curve_node_const_iterator;

    Intersection_point_node() {}
    
    Intersection_point_node(const Curve_node& cv) : intersect_p(cv.get_rightmost_point()) {
      //intersect_p = cv.get_rightmost_point();
      curves.push_back(cv);  
    }
    
    Intersection_point_node(const Curve_node& cv, const Point_plus& ref_point) : intersect_p(ref_point) {
      curves.push_back(cv);  
    }
    
    Intersection_point_node(const Curve_node& cv1, const Curve_node& cv2, const Point_plus& ref_point) : intersect_p(ref_point) {
      
      Traits traits;
 
      Comparison_result result = traits.curve_compare_at_x_right (cv1.get_curve(), cv2.get_curve(), ref_point.point());

      if (result == SMALLER){
        curves.push_back(cv1);
        curves.push_back(cv2);
      }
      else if (result == LARGER){
        curves.push_back(cv2);
        curves.push_back(cv1);
      }
      else { //equal: if one of the curve is to left of ref_point - this curve will be pushed first.
        if (rightmost(traits.curve_source(cv1.get_curve()), traits.curve_target(cv1.get_curve())) == ref_point.point()){
          curves.push_back(cv1);
          curves.push_back(cv2);
        }
        else { //include the case : if (rightmost(traits.curve_source(cv2.get_curve()), traits.curve_target(cv2.get_curve())) == ref_point)
          curves.push_back(cv2);
          curves.push_back(cv1);
        } 
      }
    } 
    
    // getting a node as an input, and merge it with our Intersection point node. This function orders all the curves enemating from intersect_p in counter clock wise order, particularly, the order when comparing the curves to the right of the intersection point is monotonicaly increasing.
    void merge(const Self& point_node){
      assert (intersect_p == point_node.get_point());
      
      assert (curves.size() > 0 && point_node.curves_begin() != point_node.curves_end());
      
      Traits  traits;
      Curve_node_container    merged_left_curves, merged_right_curves;
      //Curve_node_container    merged_curves;
      
      Curve_node_iterator        cv_iter = curves.begin();
      Curve_node_const_iterator  point_node_cv_iter = point_node.curves_begin(); 
      
      for ( ; cv_iter != curves.end() &&  point_node_cv_iter != point_node.curves_end(); ) {
        // both curves are defined to the left to the intersection point.
        if (is_left(traits.curve_source(cv_iter->get_curve()), intersect_p.point()) &&
            is_left(traits.curve_source(point_node_cv_iter->get_curve()), intersect_p.point()) ) {
          // first handle with overlappings.
          if (traits.curves_overlap(cv_iter->get_curve(), point_node_cv_iter->get_curve())){
            
            Point p1, p2;
            
            nearest_intersection_to_left (cv_iter->get_curve(), point_node_cv_iter->get_curve(), 
                                          intersect_p.point(), 
                                          //rightmost(traits.curve_source(cv_iter->get_curve()), 
                                          //        traits.curve_source(point_node_cv_iter->get_curve()) ), 
                                          p1, p2);
            
            if (p1 == intersect_p.point()){
#ifdef  CGAL_SWEEP_LINE_DEBUG
              cout<<cv_iter->get_curve()<<" and "<<point_node_cv_iter->get_curve()<<" overlap"<<endl;
#endif
              merged_left_curves.push_back(*cv_iter);
              cv_iter++;
              merged_left_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
              continue;
            }
          }
          
          Comparison_result result = traits.curve_compare_at_x_left (cv_iter->get_curve(), point_node_cv_iter->get_curve(), intersect_p.point());
          if (result == LARGER){
            merged_left_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          else if (result == SMALLER){
            merged_left_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          
          else { // at least one of the curves is vertical.
            // Now we shall insert the non vertical curve before the vertical one.
            if ( !traits.curve_is_vertical(cv_iter->get_curve()) ){
              merged_left_curves.push_back(*cv_iter);
              cv_iter++;
            }
            else{  // *cv_iter is vertical , *point_node_cv_iter may or may not be vertical, if it is not - has to come first, else the order is trivial.
              merged_left_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
            }
            //if both curves are vertical, they overlap and hence this case is already has taken care.
          }
        }
        
        // if the curves are defined only to the right of the intersection point.
        else if (traits.curve_source(cv_iter->get_curve()) == intersect_p.point() && 
                 traits.curve_source(point_node_cv_iter->get_curve()) ==  intersect_p.point() ) {
          // first handle with overlappings.
          if (traits.curves_overlap(cv_iter->get_curve(), point_node_cv_iter->get_curve())){
            Point p1, p2;
            
            traits.nearest_intersection_to_right (cv_iter->get_curve(), point_node_cv_iter->get_curve(), 
                                                  intersect_p.point(), p1, p2);
            
            if (p2 == intersect_p.point()){
#ifdef  CGAL_SWEEP_LINE_DEBUG
              cout<<cv_iter->get_curve()<<" and "<<point_node_cv_iter->get_curve()<<" overlap"<<endl;
#endif
              merged_right_curves.push_back(*cv_iter);
              cv_iter++;
              merged_right_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
              continue;
            }
          }
          
          Comparison_result result = traits.curve_compare_at_x_right (cv_iter->get_curve(), point_node_cv_iter->get_curve(), intersect_p.point());
          if (result == SMALLER){
            merged_right_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          else if (result == LARGER){
            merged_right_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          
          else { //equal. We get here if one of the curves is vertical.
            
            // Now we shall insert the non vertical curve before the vertical one.
            if ( !traits.curve_is_vertical(cv_iter->get_curve()) ){
              merged_right_curves.push_back(*cv_iter);
              cv_iter++;
            }
            else{  // *cv_iter is vertical , *point_node_cv_iter may or may not be vertical, if it is not - has to come first, else the order is trivial.
              merged_right_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
            }
            //if both curves are vertical, they overlap and hence this case is already has taken care.
          }
        }
        
        else{
          // Checking whether each curves starts at intersect_p - it means that lexicographically it's not defined to the left of intersect_p.
          if (is_left(traits.curve_source(cv_iter->get_curve()), intersect_p.point()) ){
            merged_left_curves.push_back(*cv_iter);
            cv_iter++;
          }
          else  if (traits.curve_source(cv_iter->get_curve()) == intersect_p.point() ){
            merged_right_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          if (is_left(traits.curve_source(point_node_cv_iter->get_curve()), intersect_p.point() )){
            merged_left_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          else  if (traits.curve_source(point_node_cv_iter->get_curve()) == intersect_p.point() ){
            merged_right_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
        }
      }
      
      for (; cv_iter != curves.end(); cv_iter++){
        if (traits.curve_target(cv_iter->get_curve()) == intersect_p.point())
          merged_left_curves.push_back(*cv_iter);
        else if (is_right(traits.curve_target(cv_iter->get_curve()), intersect_p.point()) )
          merged_right_curves.push_back(*cv_iter);
      }
      
      for (; point_node_cv_iter != point_node.curves_end(); point_node_cv_iter++){
        if (traits.curve_target(point_node_cv_iter->get_curve()) == intersect_p.point())
          merged_left_curves.push_back(*point_node_cv_iter);
        else  if (is_right(traits.curve_target(point_node_cv_iter->get_curve()), intersect_p.point()) )
          merged_right_curves.push_back(*point_node_cv_iter);
      }
      
      // now, copying the two merged vector to curves.
      curves.clear();
      copy(merged_left_curves.begin(), merged_left_curves.end(), back_inserter(curves));
      copy(merged_right_curves.begin(), merged_right_curves.end(), back_inserter(curves));
      
      //curves.insert(curves.begin(),merged_curves.begin(), merged_curves.end());
    }
    
    const Point_plus& get_point() const { return intersect_p; }
    
    //unsigned int number_of_curves() { return cv_vec.size(); }

    Curve_node_iterator  curves_begin() { return curves.begin(); }
    Curve_node_iterator  curves_end() { return curves.end(); }
    
    Curve_node_const_iterator  curves_begin() const { return curves.begin(); }
    Curve_node_const_iterator  curves_end() const { return curves.end(); }

  private:
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
    
    bool is_right(const Point &p1, const Point &p2) const 
      { return (CGAL::compare_lexicographically_xy(p1, p2) == LARGER); }

    bool is_left(const Point &p1, const Point &p2) const 
    { return (CGAL::compare_lexicographically_xy(p1, p2) == SMALLER); }

    bool nearest_intersection_to_left(const X_curve& cv1,
                                      const X_curve& cv2,
                                      const Point& pt,
                                      Point& p1,
                                      Point& p2) const 
    {
      Traits traits;
      
      Point rpt = traits.point_reflect_in_x_and_y( pt);
      X_curve rcv1 = traits.curve_reflect_in_x_and_y( cv1);
      X_curve rcv2 = traits.curve_reflect_in_x_and_y( cv2);
      
      Point rp1, rp2;
      bool result = traits.nearest_intersection_to_right(rcv1, rcv2, rpt, rp1, rp2);
    
      p1 = traits.point_reflect_in_x_and_y( rp1);
      p2 = traits.point_reflect_in_x_and_y( rp2);
      
      return result;
    }
    
    Point_plus                 intersect_p;  // hold the left most intersecting point.
    Curve_node_container       curves;   
  };
  
  template <class Point>
  class less_xy {
  public:
    inline  bool operator()(const Point& p, const Point& q) { 
      return CGAL::compare_lexicographically_xy(p,q) == CGAL::SMALLER ;
    }
  };


  // A predicate ordering two Curve nodes.
  template <class Curve_node>
  class less_yx {
  public:
    typedef Traits_                                         Traits;
    typedef typename Traits::X_curve                        X_curve; 
    typedef typename Traits::Point                          Point;
    typedef typename Traits::Curve_point_status             Curve_point_status;

    inline  bool operator()(const Curve_node& cv1, const Curve_node& cv2) { 
      Traits traits;
      Comparison_result result;
      
      // first check if one of the curves is vertical and the other is not.  (the traits return EQUAL in this case).
      if (traits.curve_is_vertical(cv2.get_curve()))
        if (traits.curve_is_vertical(cv1.get_curve())) // equal - the curves are overlapping.
          //cout<<"overlapping"<<std::endl;          
          return ( cv1.get_curve().id() < cv2.get_curve().id() );
      
      //if (traits.curve_is_vertical(cv1.get_curve()))
      // if (!traits.curve_is_vertical(cv2.get_curve()))
      //  return false;
              
      //if ( cv1.get_rightmost_point() != cv2.get_rightmost_point() ) //also solving the case that one of the curves end at that point.
      //result = traits.curve_compare_at_x (cv1.get_curve(), cv2.get_curve(), 
      //                                    rightmost(cv1.get_rightmost_point(), cv2.get_rightmost_point()) );
      //else
      
      if ( cv1.get_rightmost_point() != cv2.get_rightmost_point() ){     
        result = curve_node_compare_at_x(cv1, cv2,  rightmost(cv1.get_rightmost_point().point(),cv2.get_rightmost_point().point()) );
      } 
      else{ // both curves enamting the same point.
        if ( traits.curve_is_vertical(cv1.get_curve()) )  
          // if one of the curves enamating from the point is vertical then it will have larger value on the status.
          return false;
        else if ( traits.curve_is_vertical(cv2.get_curve()) )
          return true;
          
        result = traits.curve_compare_at_x_right (cv1.get_curve(), cv2.get_curve(), cv2.get_rightmost_point().point() );
        if (result == EQUAL)
          result = traits.curve_compare_at_x (cv1.get_curve(), cv2.get_curve(), cv2.get_rightmost_point().point() );
        
      }
      
      if (result == SMALLER){
        return true;
      }
      else if (result == LARGER){
        return false;
      }
      
      else { // equal - the curves are overlapping.
        
        return ( cv1.get_curve().id() < cv2.get_curve().id() );  
        // comparing the pointers of the parent halfedges.
      }
    }  
    
  private:
    Comparison_result curve_node_compare_at_x(const Curve_node &cv1, const Curve_node &cv2, const Point &q) {
      Traits traits;
      Comparison_result result;
      X_curve first_cv = cv1.get_curve(), second_cv = cv2.get_curve();
        
      // taking care the edge case of vertical segments: if one is vertical, then the order is set according the 'rightmost' point (so far) of the vertical curve comparing the point that a vertical line hits the other curve. 
      if ( traits.curve_is_vertical(cv1.get_curve()) ){
          
        // if the curve is vertical and the point is its target we refer that curve as a point.
        if (traits.curve_target(cv1.get_curve()) == cv1.get_rightmost_point().point()){
          Curve_point_status p_status = traits.curve_get_point_status(cv2.get_curve(), traits.curve_target(cv1.get_curve()) );
          if (p_status == Traits::UNDER_CURVE || p_status == Traits::ON_CURVE) // if the target is on cv2 its means that it tangent to cv2 from below.
            return  SMALLER;
          else if (p_status == Traits::ABOVE_CURVE)
            return  LARGER;
          else 
            return  EQUAL;
        }
          
        X_curve tmp_cv;
        if (traits.curve_source(cv1.get_curve()) != cv1.get_rightmost_point().point())
          traits.curve_split(cv1.get_curve(), tmp_cv, first_cv, cv1.get_rightmost_point().point());
      }
      
      if ( traits.curve_is_vertical(cv2.get_curve()) ){
          
        // if the curve is vertical and the point is its target we refer that curve as a point.
        if (traits.curve_target(cv2.get_curve()) == cv2.get_rightmost_point().point()){
          Curve_point_status p_status = traits.curve_get_point_status(cv1.get_curve(), traits.curve_target(cv2.get_curve()) );
          if (p_status == Traits::UNDER_CURVE ||  p_status == Traits::ON_CURVE) // if the target is on cv1 its means that it tangent to cv1 from below.
            return  LARGER;
          else if (p_status == Traits::ABOVE_CURVE)
            return  SMALLER;
          else  
            return EQUAL;
        }
        
        X_curve tmp_cv;
        if (traits.curve_source(cv2.get_curve()) != cv2.get_rightmost_point().point())
          traits.curve_split(cv2.get_curve(), tmp_cv, second_cv, cv2.get_rightmost_point().point());
      }    
      
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"first cv and second cv are "<<first_cv<<" "<<second_cv<<std::endl; 
#endif
      
      result = traits.curve_compare_at_x_right (first_cv, second_cv, 
                                                rightmost(cv1.get_rightmost_point().point(), cv2.get_rightmost_point().point()));
      if (result == EQUAL)
        result = traits.curve_compare_at_x (first_cv, second_cv,
                                            rightmost(cv1.get_rightmost_point().point(),cv2.get_rightmost_point().point()));
      
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"first cv "<<first_cv<<" is ";
      if (result == SMALLER)
        cout<<" below ";
      else if (result == LARGER)
        cout<<" above ";
      else
        cout<<" equal ";
      cout<<" second cv are "<<second_cv<<std::endl;
      
      cout<<"rightmost point of both is "<<rightmost(cv1.get_rightmost_point().point(), cv2.get_rightmost_point().point())<<std::endl;
#endif    
      
      // if one of the curves is vertical - EQUAL means that the other curve is between the source and target point of the vertical curve, and for our definition - it means that the verical curve comes first.
      if (result == EQUAL && traits.curve_is_vertical(first_cv)){
       
        // if first_cv is vertical and its source tangent to second_cv - it means that first_cv is above second_cv.
        Curve_point_status p_status = traits.curve_get_point_status(second_cv, traits.curve_source(first_cv));
        if (p_status == Traits::ON_CURVE)
          result = LARGER;

        else
          result = SMALLER;
      }
      else if (result == EQUAL && traits.curve_is_vertical(second_cv)){

        // if second_cv is vertical and its source tangent to first_cv - it means that second_cv is above first_cv.
        Curve_point_status p_status = traits.curve_get_point_status(first_cv, traits.curve_source(second_cv));
        if (p_status == Traits::ON_CURVE)
          result = SMALLER;
        
        else
          result = LARGER;
      }
      
      return result;
    } 
    
    const Point& rightmost(const Point &p1, const Point &p2) const {
      Comparison_result rx = CGAL::compare_lexicographically_xy(p1, p2);
      if (rx == SMALLER)
        return p2;
      else if (rx == LARGER)
        return p1;
      else
      assert(0);
      }
      
      };*/

public:

  typedef Container_                                    Container;  
  typedef Curve_iterator_                               Curve_iterator;
  typedef Traits_                                       Traits;
  typedef Point_handle<Traits>                          Point_handle_traits;
  
  typedef typename  Traits::X_curve                        X_curve;
  typedef typename  Traits::Point                          Point;

  typedef Sweep_curves_base<Curve_iterator, Traits, Point_handle_traits>  Base;
  
  typedef typename Base::X_curve_plus                X_curve_plus;
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
  

  void  sweep_curves_to_subcurves(Curve_iterator curves_begin, 
                                  Curve_iterator curves_end, 
                                  Container &subcurves, 
                                  bool overlapping = false)
  { 
    Traits                traits;
    Event_queue           event_queue;
    Status_line           status;
    
    //int c_sweep_t;
    //c_sweep_t = clock();
    
#ifdef  CGAL_SWEEP_LINE_DEBUG
    unsigned int n = 0;
    for (Curve_iterator cv_iter = curves_begin; 
         cv_iter !=  curves_end; cv_iter++, n++);
    cout<<"number of edges on input "<< n <<std::endl;
#endif
    
    // splitting all curves to x-monotone curves.
    X_curve_list  x_monotone_curves;
    for(Curve_iterator cv_iter = curves_begin; 
        cv_iter != curves_end; cv_iter++){
      if (!traits.is_x_monotone(*cv_iter)) {
        X_curve_list x_monotone_subcurves;
        traits.make_x_monotone(*cv_iter, x_monotone_subcurves);

#ifdef  CGAL_SWEEP_LINE_DEBUG
        std::cout<<"printing x-monotone parts"<<std::endl;
#endif
        for(X_curve_list_iterator iter = x_monotone_subcurves.begin(); 
            iter != x_monotone_subcurves.end(); iter++){
#ifdef  CGAL_SWEEP_LINE_DEBUG
          std::cout<<*iter<<endl;
#endif
          x_monotone_curves.push_back(*iter);  
        }
      }
      else
        x_monotone_curves.push_back(*cv_iter);
    }
       
    // now creating the Curve_node handles and the event queue.
    unsigned int id = 0;
    for(X_curve_list_iterator cv_iter = x_monotone_curves.begin(); 
        cv_iter != x_monotone_curves.end(); cv_iter++, id++){
      
      X_curve cv(*cv_iter);
      if (is_right(traits.curve_source(*cv_iter), 
                   traits.curve_target(*cv_iter)) )
        cv = traits.curve_flip(*cv_iter);
      
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<cv<<std::endl;
#endif
      
      typename Event_queue::iterator  edge_point = 
        event_queue.find( traits.curve_source(cv) );
      // defining one cv_node for both source and target event points. 
      Curve_node  cv_node = Curve_node(X_curve_plus(cv, id), 
                                       traits.curve_source(cv)); 
      
      Intersection_point_node  source_point_node = 
        Intersection_point_node(cv_node, traits.curve_source(cv) );
       
      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != source_point_node.get_point())
        event_queue.insert(Event_queue_value_type(traits.curve_source(cv), 
						  source_point_node));
      else
        edge_point->second.merge(source_point_node);
      
      
      edge_point = event_queue.find( traits.curve_target(cv) );

      Intersection_point_node  target_point_node = 
        Intersection_point_node(cv_node, traits.curve_target(cv) );

      if (edge_point == event_queue.end() || 
          edge_point->second.get_point() != target_point_node.get_point())
        event_queue.insert(Event_queue_value_type(traits.curve_target(cv), 
						  target_point_node));
      else
        edge_point->second.merge(target_point_node);
    }
    
     // now starting the sweeping.
    unsigned int queue_size = 0;
    bool         event_terminated = true;
    bool         event_overlap_terminated = true;
    while ( !(event_queue.empty()) ){
      queue_size++;
      // fetch the next event.
      typename Event_queue::iterator  event = event_queue.begin();

      const Point&              event_point = event->first;
      Intersection_point_node&  point_node = event->second;
      //bool                     event_terminated = true;
      
#ifdef  CGAL_SWEEP_LINE_DEBUG     
      cout<<"* * * event point is "<<event_point<<
        " and point node is "<<point_node.get_point().point()<<std::endl;
      CGAL_assertion(event_point == point_node.get_point().point());
#endif
      
      event_terminated = true; // reinitializing event_terminated 
      // to true only after the updating of the subdivision.
     
      // now continue with the sweep line.
      event_terminated = handle_one_event (event_queue, 
                                           status, 
                                           event_point, 
                                           point_node);
        
      // handling overlapping curves. 
      // On each overlapping group, we remove iteratively each 
      // curve and check for new events after the removing.
      // when finish, we reinsert to the status all the overlappting 
      // removed curves.
      for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
           cv_iter != point_node.curves_end(); cv_iter++){
        typename Status_line::iterator curr_cv_node = status.find(*cv_iter);
        if (curr_cv_node != status.end()){
          if (curr_cv_node != status.begin()){
            std::list<Curve_node>  overlapping_curves;
            
            typename Status_line::iterator lower =  --curr_cv_node;
            curr_cv_node++;
            for ( ;curr_cv_node != status.begin() && 
                    traits.curves_overlap(lower->first.get_curve(), 
                                          curr_cv_node->first.get_curve()); 
                  lower--){
              
              Point p1, p2;
              traits.nearest_intersection_to_right (
                                   curr_cv_node->first.get_curve(), 
                                   lower->first.get_curve(), 
                                   event_point, p1, p2);
              
              if (is_left(event_point, p1) || is_right(event_point, p2)) 
                // means that the overlapping does not intersects the 
                // status line.
                break;
              
              overlapping_curves.push_back(curr_cv_node->first);
              status.erase(curr_cv_node);
              
              curr_cv_node = lower;
              Curve_node  cv_node = curr_cv_node->first;
              
              // now taking care of the events created by the overlapping 
              // curve.
              
              Point xp;
              if (check_status_neighbors_intersections(event_queue, 
                                                       status, 
                                                       curr_cv_node, 
                                                       event_point, xp))
                if (xp == event_point) 
                  // Edge case of tangency in the event point, if it is this 
                  // event will be taked cared again.
                  event_overlap_terminated = false;
              
              if (curr_cv_node != status.begin()){
                curr_cv_node--;
                if (check_status_neighbors_intersections(event_queue, 
                                                         status, 
                                                         curr_cv_node, 
                                                         event_point, xp))
                  if (xp == event_point)  
                    // Edge case of tangency in the event point, if it is - this event will be taked cared again.
                    event_overlap_terminated = false;
                curr_cv_node++;
              }
              
              if (curr_cv_node != status.end() && 
                  CGAL::compare_lexicographically_xy (event_point, 
                       traits.curve_target(curr_cv_node->first.get_curve()) ) 
                  == CGAL::EQUAL){
                
                typename Status_line::iterator prev_cv_node;
                bool first = true;
                // hold the (lower) neighbor element of the current.
                if (curr_cv_node != status.begin()){
                  prev_cv_node = --curr_cv_node;
                  curr_cv_node++;
                  first = false;
                }
#ifdef  CGAL_SWEEP_LINE_DEBUG  
                cout<<"the event point is the right point of the curve - the curve leaves the status"<<std::endl;
#endif
                status.erase(curr_cv_node);
          
                CGAL_expensive_postcondition_code(is_valid(status));
                
                if (!first){
                  //cout<<"checking neighbors after deletion\n";
                  if (check_status_neighbors_intersections(event_queue, 
                                                           status, 
                                                           prev_cv_node, 
                                                           event_point, xp))
                    if (xp == event_point) 
                      // Edge case of tangency in the event point, if it is 
                      // this event will be taked cared again.
                      event_overlap_terminated = false;
                }
              } 
            }
            // reinsert to the status line all the overlapping removed curves.
            for (typename std::list<Curve_node>::iterator  
                   ovlp_iter = overlapping_curves.begin(); 
                 ovlp_iter != overlapping_curves.end();  ovlp_iter++)
              status.insert(Status_line_value_type(*ovlp_iter, ovlp_iter->get_curve()));
          }
        }
        
        curr_cv_node = status.find(*cv_iter);
        if (curr_cv_node != status.end()){
          std::list<Curve_node>  overlapping_curves;
          
          typename Status_line::iterator upper =  ++curr_cv_node;
          curr_cv_node--;
        
          
          for ( ;upper != status.end() && curr_cv_node != status.end() && 
                  traits.curves_overlap(curr_cv_node->first.get_curve(), 
                                        upper->first.get_curve()); upper++){

            Point p1, p2;
            traits.nearest_intersection_to_right (
                                 curr_cv_node->first.get_curve(), 
                                 upper->first.get_curve(), 
                                 event_point , p1, p2);
            
            // means that the overlapping does not intersects the status line
            if (is_left(event_point, p1) || is_right(event_point, p2))
              break;
            
            overlapping_curves.push_back(curr_cv_node->first);
            status.erase(curr_cv_node);
            
            curr_cv_node = upper;
            Curve_node  cv_node = curr_cv_node->first;
            
            // now taking care of the events created by the overlapping curve.
            Point xp;
            if (check_status_neighbors_intersections(event_queue, status, 
                                                     curr_cv_node, 
                                                     event_point, xp))
              if (xp == event_point) 
                // Edge case of tangency in the event point, if it is this 
                // event will be taked cared again.
                event_overlap_terminated = false;
            
            if (curr_cv_node != status.begin()){
              curr_cv_node--;
              if (check_status_neighbors_intersections(event_queue, 
                                                       status, 
                                                       curr_cv_node, 
                                                       event_point, xp))
                if (xp == event_point)  
                  // Edge case of tangency in the event point, if it is - 
                  // this event will be taked cared again.
                  event_overlap_terminated = false;
             
              curr_cv_node++;
            }
            
            if (curr_cv_node != status.end() && 
                CGAL::compare_lexicographically_xy (event_point, 
                       traits.curve_target(curr_cv_node->first.get_curve()) ) 
                == CGAL::EQUAL){
              
              typename Status_line::iterator prev_cv_node;
              bool first = true;
              // hold the (lower) neighbor element of the current.
              if (curr_cv_node != status.begin()){
                prev_cv_node = --curr_cv_node;
                curr_cv_node++;
                first = false;
              }
#ifdef  CGAL_SWEEP_LINE_DEBUG  
              cout<<"the event point is the right point of the curve - the curve leaves the status"<<std::endl;
#endif
              status.erase(curr_cv_node);
              
              CGAL_expensive_postcondition_code(is_valid(status));
              if (!first){
                //cout<<"checking neighbors after deletion\n";
                if (check_status_neighbors_intersections(event_queue, 
                                                         status, 
                                                         prev_cv_node, 
                                                         event_point, xp))
                  if (xp == event_point) 
                    // Edge case of tangency in the event point, if it is this event will be taked cared again.
                    event_overlap_terminated = false;
              }
            } 
          }
          // reinsert to the status line all the overlapping removed curves.
          for (typename std::list<Curve_node>::iterator  ovlp_iter = 
                 overlapping_curves.begin(); 
               ovlp_iter != overlapping_curves.end();  ovlp_iter++)
            status.insert(Status_line_value_type(*ovlp_iter, ovlp_iter->get_curve()));
        }
      }
      
      if (!event_terminated){
        handle_one_event (event_queue, status, event_point, point_node);
      }
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"Printing status line "<<std::endl;
      print_status(status);   
#endif
      
      for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
           cv_iter != point_node.curves_end(); cv_iter++){
        if (event_point != traits.curve_source(cv_iter->get_curve()) && 
            event_point == cv_iter->get_rightmost_point().point())
          cv_iter->erase_rightmost_point();
      }

      // now, updating the planar map (or arrangement) according the 
      // curves enemating from the currnet event point.
      update_subcurves(point_node, subcurves, overlapping);

      // updating all the new intersection nodes of the curves 
      // participating within the event.
      for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
           cv_iter != point_node.curves_end(); cv_iter++){
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
    //std::cout<<"The time required by sweep proccess: "
    // << (double) c_sweep_t / (double) CLOCKS_PER_SEC<<std::endl;
  }  
   
private:

  /*bool  handle_one_event (Event_queue &event_queue, Status_line& status, const Point& event_point, Intersection_point_node& point_node)
    {
    bool event_terminated = true;
    // now continue with the sweep line.
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); cv_iter != point_node.curves_end(); cv_iter++){
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"curve at event is "<<cv_iter->get_curve()<<std::endl;
      print_points_on_curves(*cv_iter);
#endif
      
      Status_line::iterator curr_cv_node = status.find(*cv_iter);
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      if (curr_cv_node != status.end()){
        cout<<"found curve "<<curr_cv_node->second<<std::endl;
        cout<<"with right most point "<<curr_cv_node->first.get_rightmost_point().point()<<std::endl;
      }
#endif
      // the event point is not the right point of the curve - insert curve to status if it's not already there.  
      if (curr_cv_node != status.end()){
          
        // remove cv_iter from status and reinsert it with the new point.
        //cout<<"found curve on status, erasing it from status"<<std::endl;
        status.erase(curr_cv_node);
      }
      
      CGAL_expensive_postcondition_code(is_valid(status));
      
      // for dubugging!
      //assert(traits.curve_get_point_status(cv_iter->get_curve(), event_point) == Traits::ON_CURVE);
      
      if (event_point != cv_iter->get_rightmost_point().point())
        cv_iter->push_event_point(point_node.get_point());
      
      // inserting the new event point to curve node, and then to the status line.
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"inserting the new event point to curve node, and then to the status line"<<std::endl;
#endif
      Status_line::iterator new_cv_node = status.insert(Status_line::value_type(*cv_iter, cv_iter->get_curve()));
      CGAL_expensive_postcondition_code(is_valid(status));
      
      Point xp;
      if (check_status_neighbors_intersections(event_queue, status, new_cv_node, event_point, xp))
        if (xp == event_point) 
          // Edge case of tangency in the event point, if it is this event will be taked cared again.
          event_terminated = false;
      
      if (new_cv_node != status.begin()){
        new_cv_node--;
        if (check_status_neighbors_intersections(event_queue, status, new_cv_node, event_point, xp))
          if (xp == event_point)  
            // Edge case of tangency in the event point, if it is - this event will be taked cared again.
            event_terminated = false;
        
        new_cv_node++;
        }
        
        if (new_cv_node != status.end() && 
        CGAL::compare_lexicographically_xy (event_point, traits.curve_target(cv_iter->get_curve()) ) == CGAL::EQUAL){
        
        //Status_line::iterator curr_cv_node = status.find(*cv_iter);
        
        // inserting the curve with its intersection points to the reults list.
        //Curve_node final_cv(new_cv_node->first);
        //final_cv.push_event_point(event_point);
        //disjoint_interior_curves.push_back(final_cv);
        
        Status_line::iterator prev_cv_node;
        bool first = true;
        // hold the (lower) neighbor element of the current.
        if (new_cv_node != status.begin()){
          prev_cv_node = --new_cv_node;
          new_cv_node++;
          first = false;
        }
#ifdef  CGAL_SWEEP_LINE_DEBUG  
        cout<<"the event point is the right point of the curve - the curve leaves the status"<<std::endl;
#endif
        status.erase(new_cv_node);
        
        CGAL_expensive_postcondition_code(is_valid(status));
        if (!first){
          //cout<<"checking neighbors after deletion\n";
          if (check_status_neighbors_intersections(event_queue, status, prev_cv_node, event_point, xp))
            if (xp == event_point) 
              // Edge case of tangency in the event point, if it is this event will be taked cared again.
              event_terminated = false;
        }
      }
    }
    return event_terminated;
  }
    
  bool  point_is_on_curve_interior(const Point& p, const X_curve& cv) {
  return (p != traits.curve_source(cv) || p != traits.curve_target(cv) );
  }
  
  // check whether two curves are tangent on a given point.
  bool  curves_tangent_at_point(const X_curve& cv1, const X_curve& cv2, const Point& point)
  {
  Traits traits;
  
  return (traits.curve_get_point_status(cv1, point) == Traits::ON_CURVE && 
  traits.curve_get_point_status(cv2, point) == Traits::ON_CURVE && 
  (point_is_on_curve_interior(point, cv1) || point_is_on_curve_interior(point, cv2) ) );
  
  }
  
  bool  check_status_neighbors_intersections(Event_queue &event_queue, Status_line& status,  
  Status_line::iterator lower_neighbor, 
                                             const Point &event_point, Point& point)
  {
    
#ifdef  CGAL_SWEEP_LINE_DEBUG  
    cout<<"Printing status line "<<std::endl;
    print_status(status);   
#endif

    Traits traits;
    const Curve_node& cv1 = lower_neighbor->first;
    
    Status_line::iterator next_neighbor = ++lower_neighbor;
    if (next_neighbor == status.end())
      return false;
    
    lower_neighbor--;
    
    const Curve_node&  cv2 = next_neighbor->first;
    
#ifdef  CGAL_SWEEP_LINE_DEBUG  
    cout<<"cv1 and cv2 are "<<cv1.get_curve()<<" "<<cv2.get_curve()<<std::endl;
#endif

    // in each node - checking intersections between two adjacent curves and updating the event queue if needed.  
    Point xp1 = event_point, xp2, xp3;
    
    //Point ref_point(lower_neighbor->first.get_rightmost_point().point());
    //if (is_left(next_neighbor->first.get_rightmost_point().point(), lower_neighbor->first.get_rightmost_point().point()))
    //  ref_point = next_neighbor->first.get_rightmost_point().point();
    
    //ref_point = event_point;  //I added this to avoid calcuating points that are left to the event point. When dealing with overlapping this bug may apear.
    
    bool  curves_tangent =  curves_tangent_at_point(cv1.get_curve(), cv2.get_curve(), event_point);
    bool  curves_intersect = traits.nearest_intersection_to_right(cv1.get_curve(), cv2.get_curve(), event_point, xp2, xp3);
    
    if (curves_tangent || curves_intersect){  
      
      // first arrange the intersection points xp1, xp2, xp3 if on of the conditoins (tangency or intersection) does not hold.
      if (!curves_tangent){
        xp1 = xp2;
        xp2 = xp3;
      }
      
      if (!curves_intersect)
        xp2 = xp3 = xp1;
      
      //if (is_left(xp1, event_point))
      // cout<<"interection point " <<xp1<<" is left to event point "<<point<<endl;
      
      // we choose the leftmost point because the intersection point may be one of the points on lower_neighbor or next_neighbor, specially it can be a tangent point, and the function nearest_intersection_to_right may return false (it's third parameter is the intersection point itself!).
      
#ifdef  CGAL_SWEEP_LINE_DEBUG  
      cout<<"rightmost point on status, xp1 and xp2 are "<<lower_neighbor->first.get_rightmost_point().point()<<" "<<xp1<<" "<<xp2<<std::endl;
#endif

      // have to handle overlapping.
      //if (xp1 == lower_neighbor->first.get_rightmost_point().point())
      //  xp1 = xp2;
      
      // for debugging.
      //if (traits.curve_get_point_status(cv1.get_curve(), xp1) != Traits::ON_CURVE)
      //  cout<<"The point "<<xp1<<" is not on the curve "<<cv1.get_curve()<<std::endl;
      //if (traits.curve_get_point_status(cv2.get_curve(), xp1) != Traits::ON_CURVE)
      //  cout<<"The point "<<xp1<<" is not on the curve "<<cv2.get_curve()<<std::endl;
      // end debugging.
      
      // checking if the curves are already participate within the intersection node, if at least one is not - we have a new event here or a new curve particitating on it and the function will return true.
     
       // handling overlapping.
      if ( (xp2 != xp3) && 
           !( (xp3 == traits.curve_source(cv1.get_curve()) || xp3 == traits.curve_target(cv1.get_curve()) ) && 
              (xp3 == traits.curve_source(cv2.get_curve()) || xp3 == traits.curve_target(cv2.get_curve()) )) ){
        
        //status_iter->first.push_event_point(xp3);
        //next_iter->first.push_event_point(xp3);
        
        Event_queue::iterator  xp_event = event_queue.find(xp3);
        bool xp3_cv1_in_queue = false,  xp3_cv2_in_queue = false;
        if (xp_event == event_queue.end())
          xp_event = event_queue.insert(Event_queue::value_type(xp3, Intersection_point_node(cv1, cv2, xp3)));
        else{
          // have to check whether the event is a new event. (we might calculated this point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (cv_iter->get_curve() == cv1.get_curve()) // traits.curve_is_same(cv_iter->get_curve(), cv1.get_curve()) && cv_iter->get_curve().id() == cv1.get_curve().id()) 
              {
                xp3_cv1_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
                std::cout<<cv1.get_curve()<<" was found "<<std::endl;
#endif
              }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              //traits.curve_is_same(cv_iter->get_curve(), cv2.get_curve()) && cv_iter->get_curve().id() == cv2.get_curve().id()){
              xp3_cv2_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
              std::cout<<cv2.get_curve()<<" was found "<<std::endl;
#endif
            }
          }
          
          if (!xp3_cv1_in_queue && !xp3_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, cv2, Point_plus(xp3)));          
          else if (!xp3_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, Point_plus(xp3)));
          else if (!xp3_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, Point_plus(xp3)));
        }
        //  return (!xp1_cv1_in_queue || !xp1_cv2_in_queue);
      }

      // handling overlapping.
      if ( (xp1 != xp2) && 
           !( (xp2 == traits.curve_source(cv1.get_curve()) || xp2 == traits.curve_target(cv1.get_curve()) ) && 
              (xp2 == traits.curve_source(cv2.get_curve()) || xp2 == traits.curve_target(cv2.get_curve()) )) ){
        
        //status_iter->first.push_event_point(xp2);
        //next_iter->first.push_event_point(xp2);
        
        Event_queue::iterator  xp_event = event_queue.find(xp2);
        bool xp2_cv1_in_queue = false,  xp2_cv2_in_queue = false;
        if (xp_event == event_queue.end())
          xp_event = event_queue.insert(Event_queue::value_type(xp2, Intersection_point_node(cv1, cv2, xp2)));
        else{
          // have to check whether the event is a new event. (we might calculated this point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (cv_iter->get_curve() == cv1.get_curve()) // traits.curve_is_same(cv_iter->get_curve(), cv1.get_curve()) && cv_iter->get_curve().id() == cv1.get_curve().id()) 
              {
                xp2_cv1_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
                std::cout<<cv1.get_curve()<<" was found "<<std::endl;
#endif
              }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              //traits.curve_is_same(cv_iter->get_curve(), cv2.get_curve()) && cv_iter->get_curve().id() == cv2.get_curve().id()){
              xp2_cv2_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
              std::cout<<cv2.get_curve()<<" was found "<<std::endl;
#endif
            }
          }
          
          if (!xp2_cv1_in_queue && !xp2_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, cv2, Point_plus(xp2)));          
          else if (!xp2_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, Point_plus(xp2)));
          else if (!xp2_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, Point_plus(xp2)));
        }
        //  return (!xp1_cv1_in_queue || !xp1_cv2_in_queue);
      }
      
      bool xp1_cv1_in_queue = false,  xp1_cv2_in_queue = false; 
      // if cv1 and cv2 have common edge point - we do not consider it as an intersection point. 
      if ( !( (xp1 == traits.curve_source(cv1.get_curve()) || xp1 == traits.curve_target(cv1.get_curve()) ) && 
              (xp1 == traits.curve_source(cv2.get_curve()) || xp1 == traits.curve_target(cv2.get_curve()) ) ) ){
        
        //status_iter->first.push_event_point(xp1);
        //next_iter->first.push_event_point(xp1);
        
        Event_queue::iterator  xp_event = event_queue.find(xp1);
        //bool xp1_cv1_in_queue = false,  xp1_cv2_in_queue = false;
        if (xp_event == event_queue.end())
          xp_event = event_queue.insert(Event_queue::value_type(xp1, Intersection_point_node(cv1, cv2, xp1)));
        else{
          // have to check whether the event is a new event. (we might calculated this point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (cv_iter->get_curve() == cv1.get_curve()) // traits.curve_is_same(cv_iter->get_curve(), cv1.get_curve()) && cv_iter->get_curve().id() == cv1.get_curve().id()) 
              {
                xp1_cv1_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
                std::cout<<cv1.get_curve()<<" was found "<<std::endl;
#endif
              }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              //traits.curve_is_same(cv_iter->get_curve(), cv2.get_curve()) && cv_iter->get_curve().id() == cv2.get_curve().id()){
              xp1_cv2_in_queue = true;
#ifdef  CGAL_SWEEP_LINE_DEBUG  
              std::cout<<cv2.get_curve()<<" was found "<<std::endl;
#endif
            }
          }
          
          if (!xp1_cv1_in_queue && !xp1_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, cv2, Point_plus(xp1)));
          else if (!xp1_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, Point_plus(xp1)));
          else if (!xp1_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, Point_plus(xp1)));
        }
        point = xp1;
        return (!xp1_cv1_in_queue || !xp1_cv2_in_queue);
        }
        }
        return false;
        }*/
  
  template <class Container>
  void  update_subcurves(Intersection_point_node& point_node, 
                         Container &subcurves, 
                         bool overlapping) 
  {
    if (overlapping)
      update_subcurves_with_overlappings(point_node, subcurves);
    else
      update_subcurves_without_overlappings(point_node, subcurves);
    
        /*#ifdef  CGAL_SWEEP_LINE_DEBUG
    cout<<"--------- updating map with point node"<<point_node.get_point().point() <<std::endl;
    for (Curve_node_iterator cv_iter1= point_node.curves_begin(); cv_iter1 != point_node.curves_end(); cv_iter1++){
      cv_iter1++;
      Curve_node_iterator cv_iter2 = cv_iter1;
      cv_iter1--;
      for ( ; cv_iter2 != point_node.curves_end(); cv_iter2++){   
        if (traits.curves_overlap(cv_iter1->get_curve(), cv_iter2->get_curve()))
          cout<<"update_subdivision "<<cv_iter1->get_curve()<<" and "<< cv_iter2->get_curve() <<" are overlapping"<<endl;
      }
    }
#endif    

    X_curve prev_sub_cv;
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); cv_iter != point_node.curves_end(); cv_iter++){
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"now handling "<<cv_iter->get_curve()<<endl;
#endif
      
      if (is_left(cv_iter->get_rightmost_point().point(), point_node.get_point().point())) { // means we have a new sub curve to insert.
        
        // first splitting the curve in the points cv_iter->get_rightmost_point().point() and point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), right_cv = cv;

        //cout<<"cv is "<<cv<<endl;
        if (traits.curve_source(cv) != cv_iter->get_rightmost_point().point() &&
            traits.curve_target(cv) != cv_iter->get_rightmost_point().point() ) {
          traits.curve_split(cv, sub_cv, right_cv, cv_iter->get_rightmost_point().point());
          
          cv = right_cv;
        }
        
        //cout<<"after first splitting cv is "<<cv<<endl;
        //cout<<"source and target of cv are"<<traits.curve_source(cv)<<" "<<traits.curve_target(cv)<<endl;
        if (traits.curve_source(cv) != point_node.get_point().point() &&
            traits.curve_target(cv) != point_node.get_point().point())
          traits.curve_split(cv, sub_cv, right_cv, point_node.get_point().point());
        else
          sub_cv = right_cv;
        
#ifdef  CGAL_SWEEP_LINE_DEBUG
        cout<<"inserting "<<sub_cv<<endl;
#endif

        //prev_sub_cv = sub_cv;

        subcurves.push_back(sub_cv);
        
        // else - no new sub curve is inserted to the subdivision.
      }
      }*/
  }
  
  template <class Container>
  void  update_subcurves_with_overlappings(Intersection_point_node& point_node, 
                                           Container &subcurves) 
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
        if (traits.curves_overlap(cv_iter1->get_curve(), 
                                  cv_iter2->get_curve()))
          cout<<"update_subdivision "<<cv_iter1->get_curve()<<
            " and "<< cv_iter2->get_curve() <<" are overlapping"<<endl;
      }
    }
#endif    

    X_curve prev_sub_cv;
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); cv_iter++){
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"now handling "<<cv_iter->get_curve()<<endl;
#endif
      
      if (is_left(cv_iter->get_rightmost_point().point(), 
                  point_node.get_point().point())) { 
        // means we have a new sub curve to insert.
        
        // first splitting the curve in the points 
        // cv_iter->get_rightmost_point().point() 
        // and point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), 
          right_cv = cv;

        //cout<<"cv is "<<cv<<endl;
        if (traits.curve_source(cv) != cv_iter->get_rightmost_point().point() 
            &&
            traits.curve_target(cv) != cv_iter->get_rightmost_point().point()) 
          {
            traits.curve_split(cv, sub_cv, right_cv, 
                               cv_iter->get_rightmost_point().point());
          
            cv = right_cv;
          }
        
        //cout<<"after first splitting cv is "<<cv<<endl;
        //cout<<"source and target of cv are"<<traits.curve_source(cv)
        // <<" "<<traits.curve_target(cv)<<endl;
        
        if (traits.curve_source(cv) != point_node.get_point().point() &&
            traits.curve_target(cv) != point_node.get_point().point())
          traits.curve_split(cv, sub_cv, 
                             right_cv, point_node.get_point().point());
        else
          sub_cv = right_cv;

        /*  if (cv_iter != point_node.curves_begin()){
            if (traits.curves_overlap(sub_cv, prev_sub_cv)){
            //cout<<sub_cv<<" and "<< prev_sub_cv<<" are overlapping"<<endl;
            continue;
            }
            }*/
        
#ifdef  CGAL_SWEEP_LINE_DEBUG
        cout<<"inserting "<<sub_cv<<endl;
#endif

        //prev_sub_cv = sub_cv;

        subcurves.push_back(sub_cv);
      }
      // else - no new sub curve is inserted to the subdivision.
    }
  }
  
  template <class Container> 
  void  update_subcurves_without_overlappings(
                          Intersection_point_node& point_node, 
                          Container &subcurves)
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
        if (traits.curves_overlap(cv_iter1->get_curve(), 
                                  cv_iter2->get_curve()))
          cout<<"update_subdivision "<<cv_iter1->get_curve()<<
            " and "<< cv_iter2->get_curve() <<" are overlapping"<<endl;
      }
    }
#endif    

    X_curve prev_sub_cv;
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); cv_iter++){
#ifdef  CGAL_SWEEP_LINE_DEBUG
      cout<<"now handling "<<cv_iter->get_curve()<<endl;
#endif

      if (is_left(cv_iter->get_rightmost_point().point(), 
                  point_node.get_point().point())) { 
        // means we have a new sub curve to insert.
        
        // first splitting the curve in the points 
        // cv_iter->get_rightmost_point().point() and 
        // point_node.get_point().point().
        X_curve cv = cv_iter->get_curve(), sub_cv = cv_iter->get_curve(), 
          right_cv = cv;
        
        //cout<<"cv is "<<cv<<endl;
        if (traits.curve_source(cv) != cv_iter->get_rightmost_point().point() 
            &&
            traits.curve_target(cv) != cv_iter->get_rightmost_point().point() ) 
          {
            traits.curve_split(cv, sub_cv, right_cv, 
                               cv_iter->get_rightmost_point().point());
          
            cv = right_cv;
          }
        
        //cout<<"after first splitting cv is "<<cv<<endl;
        //cout<<"source and target of cv are"<<traits.curve_source(cv)<<" "<<traits.curve_target(cv)<<endl;
        if (traits.curve_source(cv) != point_node.get_point().point() &&
            traits.curve_target(cv) != point_node.get_point().point())
          traits.curve_split(cv, sub_cv, right_cv, 
                             point_node.get_point().point());
        else
          sub_cv = right_cv;
        
        // for debugging!
        //assert(sub_cv == X_curve(cv_iter->get_rightmost_point().point(), point_node.get_point().point()));
        //cout<<"sub curve calculates is "<<endl;
        //cout<<sub_cv<<endl;
        //cout<<"And the right curve is supposed to be"<<endl;
        //cout<<X_curve( cv_iter->get_rightmost_point().point(), point_node.get_point().point())<<endl;

        if (cv_iter != point_node.curves_begin()){
          if (traits.curves_overlap(sub_cv, prev_sub_cv)){
            //cout<<sub_cv<<" and "<< prev_sub_cv<<" are overlapping"<<endl;
            continue;
          }
        }

#ifdef  CGAL_SWEEP_LINE_DEBUG
        cout<<"inserting "<<sub_cv<<endl;
#endif

        prev_sub_cv = sub_cv;

        subcurves.push_back(sub_cv);
      }
      // else - no new sub curve is inserted to the subdivision.
    }
  }
  
  //--------------------------------------------------------- debuging functions ------------------------------------------
/*  bool  is_valid(const Status_line& status){
    std::list<Curve_node> curve_nodes;
    std::list<X_curve>    subcurves;
    
    for (Status_line::const_iterator iter = status.begin(); iter != status.end(); iter++)
      curve_nodes.push_back(iter->first);

    get_subcurves(curve_nodes, subcurves);
    
    return (!do_intersect_subcurves(subcurves));
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
  
  Traits                traits;*/
};

CGAL_END_NAMESPACE

#endif
