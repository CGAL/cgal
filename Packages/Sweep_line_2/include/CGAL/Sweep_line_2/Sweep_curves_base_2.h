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
// file          : include/CGAL/Sweep_line_2/Sweep_curves_base_2.h
// package       : arr (1.87)
// maintainer    : Eyal Flato <flato@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra <estere@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================

#ifndef CGAL_SWEEP_CURVES_BASE_2_H
#define CGAL_SWEEP_CURVES_BASE_2_H

#include <vector>
#include <list>
#include <map>
#include <set>

#include <CGAL/In_place_list.h>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

#ifdef CGAL_SWEEP_LINE_DEBUG
#define CGAL_SL_DEBUG(x) x
#else
#define CGAL_SL_DEBUG(x) 
#endif

//#include <CGAL/IO/leda_window.h>  //used for visualization -

CGAL_BEGIN_NAMESPACE

template <class CurveInputIterator, class SweepLineTraits_2, 
  class Point_plus_, class X_curve_plus_>
class Sweep_curves_base_2 
{
public:
  typedef typename  SweepLineTraits_2::X_curve    X_curve;
protected:
  typedef Sweep_curves_base_2<CurveInputIterator,
    SweepLineTraits_2,Point_plus_,X_curve_plus_>  Self;

  /*!
  Curve_node_rep holds a curve participating in the sweep proccess. 
  It contains a curve and a container of points.
  The points container refers to all the intersedction points 
  (including edge points) calculated so far by the sweep line.
  These points are ordered from left to right on the curve, 
  which means they are sorted in a way we get immidiately all 
  the disjoint subcurves reduce by the curve.
  */

  class Curve_node_rep;
  class Curve_node;
  class Intersection_point_node;

  // SUNPRO requires these friends.
  friend class Curve_node_rep;
  friend class Curve_node;
  friend class Intersection_point_node;
  
  class Curve_node_rep { 
  public:
    typedef SweepLineTraits_2                   Traits;
    typedef Point_plus_                         Point_plus; 
    typedef X_curve_plus_                       X_curve_plus;
    
    typedef typename Traits::Point              Point;
    typedef std::vector<Point_plus>             Points_container;

    Curve_node_rep(const X_curve_plus& cv, const Point& p, 
		   Traits* traits_) : 
      cv_(cv), traits(traits_) {
      points.push_back( Point_plus(p) );
    }

    Curve_node_rep(const X_curve_plus& cv, 
                   const Point_plus& p, Traits* traits_) : 
      cv_(cv), traits(traits_)  {
      points.push_back(p);
    }

    ~Curve_node_rep() {}

    // here we only compare addresses of represenattion.
    bool operator==(const Curve_node_rep& rep) const 
    {
      // first at all comparing pointers to accelarate this operation.
      if (this == &rep)
        return true;
      if (&cv_ == &rep.cv_ && &points == &(rep.points))
        return true;

      return false;
    }

    bool operator!=(const Curve_node_rep& rep) const 
    {
      return !operator==(rep);
    }
    
  protected:
    friend class Curve_node;

    X_curve_plus         cv_;  // hold the left-most intersecting point.
    Points_container     points;
    Traits               *traits;

  private:
    const Point& leftmost(const Point &p1, const Point &p2) const
    {
      Comparison_result res = traits->compare_xy(p1,p2);

      return ((res == SMALLER) ? p1 : p2);
    }
  };
  
  /*!
    A handle to a curve node. This is just a wrapper class with no members.
  */
  // The handle to curve node.
  class Curve_node : public Handle_for<Curve_node_rep> {
    typedef Handle_for<Curve_node_rep>           Handle_for_Curve_node_rep;
  public:
    typedef  Point_plus_                         Point_plus;
    typedef  X_curve_plus_                       X_curve_plus;
    typedef  SweepLineTraits_2                   Traits;
    typedef  typename Traits::X_curve            X_curve; 
    typedef  typename Traits::Point              Point;
    
    typedef typename Curve_node_rep::Points_container   Points_container;
    typedef typename Points_container::iterator         Points_iterator;
    typedef typename Points_container::const_iterator  
	    					Points_const_iterator;

    Curve_node(Traits *traits_) : 
      Handle_for_Curve_node_rep(traits_) {
#ifdef CGAL_MAKE_PROFILING
      std::cout<<"alloc handle for DEFAULT Curve_node_rep"<<endl;
#endif
    }
    
    Curve_node(const X_curve_plus& cv, Traits *traits_) : 
      Handle_for_Curve_node_rep(Curve_node_rep(cv, traits_)) {
#ifdef CGAL_MAKE_PROFILING
      std::cout<<"alloc handle for Curve_node_rep(cv)" << cv << endl;
#endif
    }
    
    Curve_node(const X_curve_plus& cv, const Point& p, Traits *traits_) : 
      Handle_for_Curve_node_rep(Curve_node_rep(cv,p, traits_)) {
#ifdef CGAL_MAKE_PROFILING
      std::cout << "alloc handle for Curve_node_rep(cv,p)" << cv << " "
                << p << endl;
#endif
    }

    Curve_node(const X_curve_plus& cv, const Point_plus& p, 
	       Traits *traits_) : 
      Handle_for_Curve_node_rep(Curve_node_rep(cv,p,traits_)) {
#ifdef CGAL_MAKE_PROFILING
      std::cout<<"allocating handle for Curve_node_rep(cv,p_plus)" 
               << cv <<" "<< p.point() << endl;
#endif
    }
    
    Curve_node(const Curve_node& cv_node) : 
      Handle_for_Curve_node_rep(cv_node) {
#ifdef CGAL_MAKE_PROFILING
      std::cout<<"allocating handle for cv_node" << endl;
#endif
    }

    ~Curve_node() {}

    bool operator==(const Curve_node& cv_node) const 
    {
      return *ptr() == *(cv_node.ptr());
    }

    bool operator!=(const Curve_node& cv_node) const 
    {
      return !operator==(cv_node);
    }

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
    
    Point_plus& get_rightmost_point() { 
      CGAL_assertion ( ptr()->points.size() > 0);
      return *(ptr()->points.rbegin());
    }
    
    const Point_plus& get_rightmost_point() const { 
      CGAL_assertion ( ptr()->points.size() > 0);
      return *(ptr()->points.rbegin());
    }

    Points_iterator  points_begin() { 
      return  ptr()->points.begin(); 
    }
    
    Points_iterator  points_end() { 
      return  ptr()->points.end(); 
    }
    
    Points_const_iterator  points_begin() const { 
      return  ptr()->points.begin(); 
    }
    
    Points_const_iterator  points_end() const { 
      return  ptr()->points.end(); 
    }
    
    Self& operator=(const Self &cv_node)
    {
      Handle_for_Curve_node_rep::operator=(cv_node);
      
      return *this;
    }
    
  };

  /*!
    A class describing an intersection point.
    It contains the point itself and a list of
    all curves going through that point, ordered by
    compare_at_x and compare_at_x_right of the traits.
  */
  class Intersection_point_node {
    typedef Curve_node                                  Curve_node_;
  public:
    typedef  SweepLineTraits_2                          Traits;
    typedef  Point_plus_                                Point_plus;
    typedef  X_curve_plus_                              X_curve_plus;
    typedef  typename Traits::X_curve_2                 X_curve_2; 
    typedef  typename Traits::Point_2                   Point_2;
    typedef  Intersection_point_node                    Self;
    typedef  std::vector<Curve_node_>                   Curve_node_container;
    typedef  typename Curve_node_container::iterator    Curve_node_iterator;
    typedef  typename Curve_node_container::const_iterator  
                                                Curve_node_const_iterator;
    typedef typename Traits::Has_left_category          Has_left_category;

      // Obsolete
    typedef Point_2                                     Point;
    typedef X_curve_2                                   X_curve;

    Intersection_point_node(Traits *traits_) : traits(traits_) {}
    
    Intersection_point_node(const Curve_node_& cv, 
                            Traits *traits_) : 
      intersect_p(cv.get_rightmost_point()), traits(traits_) {
      curves.push_back(cv);  
    }
    
    Intersection_point_node(const Curve_node_& cv, 
                            const Point_plus& ref_point, 
                            Traits *traits_) : 
      intersect_p(ref_point), traits(traits_) {
      curves.push_back(cv);  
    }
    
    Intersection_point_node(const Curve_node_& cv1, 
                            const Curve_node_& cv2, 
                            const Point_plus& ref_point, 
                            Traits *traits_) : 
      intersect_p(ref_point), traits(traits_) {
      
      Comparison_result result = _curve_compare_at_x_right(cv1.get_curve(), 
							   cv2.get_curve(), 
							   ref_point.point());

      if (result == SMALLER){
        curves.push_back(cv1);
        curves.push_back(cv2);
      }
      else if (result == LARGER){
        curves.push_back(cv2);
        curves.push_back(cv1);
      }
      else { //equal: if one of the curve is to left of ref_point - 
        // this curve will be pushed first.
        if (traits->point_is_same(
              rightmost(traits->curve_source(cv1.get_curve()), 
                        traits->curve_target(cv1.get_curve())),
              ref_point.point())) {
          curves.push_back(cv1);
          curves.push_back(cv2);
        }
        else { //include the case : 
          // if (rightmost(traits.curve_source(cv2.get_curve()), 
          // traits.curve_target(cv2.get_curve())) == ref_point)
          curves.push_back(cv2);
          curves.push_back(cv1);
        } 
      }
    } 
    
    /*!
      Merges the curves of the "this" point and the input node.
      The point in the input point node has to be the same as
      the point in this instance.
  
      This function orders all the curves enemating from 
      intersect_p in counter clock wise order, particularly, the order 
      when comparing the curves to the right of the intersection point 
      is monotonicaly increasing.
    */
    void merge(const Self& point_node)
    {
      CGAL_assertion (intersect_p == point_node.get_point());
      
      CGAL_assertion (curves.size() > 0 && 
                      point_node.curves_begin() != point_node.curves_end());
      
      Curve_node_container    merged_left_curves, merged_right_curves;
      
      Curve_node_iterator       cv_iter=curves.begin();
      Curve_node_const_iterator point_node_cv_iter=point_node.curves_begin(); 
      
      for ( ; cv_iter != curves.end() &&  
              point_node_cv_iter != point_node.curves_end(); ) {
        // both curves are defined to the left to the intersection point.
        if (is_left(traits->curve_source(cv_iter->get_curve()), 
                    intersect_p.point()) &&
            is_left(traits->curve_source(point_node_cv_iter->get_curve()), 
                    intersect_p.point()) ) {
          // first handle with overlappings.
          if (traits->curves_overlap(cv_iter->get_curve(), 
                                    point_node_cv_iter->get_curve())){
            
            Point p1, p2;
            
            nearest_intersection_to_left (cv_iter->get_curve(), 
                                          point_node_cv_iter->get_curve(), 
                                          intersect_p.point(), 
                                          p1, p2);
            
            if (traits->point_is_same(p1, intersect_p.point())) {
              merged_left_curves.push_back(*cv_iter);
              cv_iter++;
              merged_left_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
              continue;
            }
          }
          
          Comparison_result result = _curve_compare_at_x_left
	                             (cv_iter->get_curve(), 
				      point_node_cv_iter->get_curve(), 
				      intersect_p.point());

          if (result == LARGER){
            merged_left_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          else if (result == SMALLER){
            merged_left_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          
          else { // at least one of the curves is vertical.
            // Now we shall insert the non vertical curve before the 
            // vertical one.
            if ( !traits->curve_is_vertical(cv_iter->get_curve()) ){
              merged_left_curves.push_back(*cv_iter);
              cv_iter++;
            }
            else{  // *cv_iter is vertical , *point_node_cv_iter may or may 
              // not be vertical, if it is not - has to come first, else the 
              // order is trivial.
              merged_left_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
            }
            // if both curves are vertical, they overlap and hence this case 
            // is already has taken care.
          }
        }
        
        // if the curves are defined only to the right of the intersection 
        // point.
        else if (traits->point_is_same(
                   traits->curve_source(cv_iter->get_curve()),
                   intersect_p.point()) && 
                 traits->point_is_same(
                   traits->curve_source(point_node_cv_iter->get_curve()),  
                   intersect_p.point())) {
          // first handle with overlappings.
          if (traits->curves_overlap(cv_iter->get_curve(), 
                                     point_node_cv_iter->get_curve())){
            Point p1, p2;
            
            traits->
              nearest_intersection_to_right(cv_iter->get_curve(), 
                                            point_node_cv_iter->get_curve(), 
                                            intersect_p.point(), p1, p2);
            
            if (traits->point_is_same(p2, intersect_p.point())) {
              merged_right_curves.push_back(*cv_iter);
              cv_iter++;
              merged_right_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
              continue;
            }
          }
          
          Comparison_result result = _curve_compare_at_x_right(
                                              cv_iter->get_curve(), 
                                              point_node_cv_iter->get_curve(), 
                                              intersect_p.point());

          if (result == SMALLER){
            merged_right_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          else if (result == LARGER){
            merged_right_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          
          else { //equal. We get here if one of the curves is vertical.
            
            // Now we shall insert the non vertical curve before the 
            // vertical one.
            if ( !traits->curve_is_vertical(cv_iter->get_curve()) ){
              merged_right_curves.push_back(*cv_iter);
              cv_iter++;
            }
            else{  // *cv_iter is vertical , *point_node_cv_iter may or may 
              // not be vertical, if it is not - has to come first, else the 
              // order is trivial.
              merged_right_curves.push_back(*point_node_cv_iter);
              point_node_cv_iter++;
            }
            // if both curves are vertical, they overlap and hence this case 
            // is already has taken care.
          }
        }
        
        else{
          // Checking whether each curves starts at intersect_p - 
          // it means that lexicographically it's not defined to the left of 
          // intersect_p.
          if (is_left(traits->curve_source(cv_iter->get_curve()), 
                      intersect_p.point()) ){
	    merged_left_curves.push_back(*cv_iter);
            cv_iter++;
          }
          else if (traits->point_is_same(
                     traits->curve_source(cv_iter->get_curve()),
                     intersect_p.point())) {
            merged_right_curves.push_back(*cv_iter);
            cv_iter++;
          }
          
          if (is_left(traits->curve_source(point_node_cv_iter->get_curve()), 
                      intersect_p.point() )){
            merged_left_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
          else if (traits->point_is_same(
                     traits->curve_source(point_node_cv_iter->get_curve()),
                     intersect_p.point())) {
            merged_right_curves.push_back(*point_node_cv_iter);
            point_node_cv_iter++;
          }
        }
      }
      
      for (; cv_iter != curves.end(); cv_iter++){
        if (traits->point_is_same(traits->curve_target(cv_iter->get_curve()),
                                  intersect_p.point()))
          merged_left_curves.push_back(*cv_iter);
        else if (is_right(traits->curve_target(cv_iter->get_curve()), 
                          intersect_p.point()) )
          merged_right_curves.push_back(*cv_iter);
      }
      
      for (; point_node_cv_iter != point_node.curves_end(); 
           point_node_cv_iter++){
        if (traits->point_is_same(
              traits->curve_target(point_node_cv_iter->get_curve()),
              intersect_p.point()))
          merged_left_curves.push_back(*point_node_cv_iter);
        else  if (is_right(traits->curve_target(
                                         point_node_cv_iter->get_curve()), 
                           intersect_p.point()) )
          merged_right_curves.push_back(*point_node_cv_iter);
      }
      
      // now, copying the two merged vector to curves.
      curves.clear();
      std::copy(merged_left_curves.begin(), 
           merged_left_curves.end(), 
           std::back_inserter(curves));
      std::copy(merged_right_curves.begin(), 
           merged_right_curves.end(), 
           std::back_inserter(curves));
    }

    
    Point_plus& get_point() { return intersect_p; }
    const Point_plus& get_point() const { return intersect_p; }

    Curve_node_iterator  curves_begin() { return curves.begin(); }
    Curve_node_iterator  curves_end() { return curves.end(); }
    
    Curve_node_const_iterator  curves_begin() const { return curves.begin(); }
    Curve_node_const_iterator  curves_end() const { return curves.end(); }

  protected:
    const Point& rightmost(const Point &p1, const Point &p2) const
    {
      Comparison_result res = traits->compare_xy(p1,p2);

      return ((res == SMALLER) ? p2 : p1);
    }
    
    bool is_right(const Point &p1, const Point &p2) const 
    { return (Compare_lexicographically_xy(p1, p2) == LARGER); }
    
    bool is_left(const Point &p1, const Point &p2) const 
    { return (Compare_lexicographically_xy(p1, p2) == SMALLER); }
    
    bool nearest_intersection_to_left(const X_curve_2 & cv1,
                                      const X_curve_2 & cv2,
                                      const Point_2 & pt,
                                      Point_2 & p1, Point_2 & p2) const 
    {
      return nearest_intersection_to_left_imp(cv1, cv2, pt, p1, p2,
                                              Has_left_category());
    }

    bool nearest_intersection_to_left_imp(const X_curve_2 & cv1,
                                          const X_curve_2 & cv2,
                                          const Point_2 & pt,
                                          Point_2 & p1, Point_2 & p2,
                                          Tag_true) const
    {
      return traits->nearest_intersection_to_left(cv1, cv2, pt, p1, p2);
    }
    
    bool nearest_intersection_to_left_imp(const X_curve_2 & cv1,
                                          const X_curve_2 & cv2,
                                          const Point_2 & pt,
                                          Point_2 & p1, Point_2 & p2,
                                          Tag_false) const 
    { 
      Point_2 rpt = traits->point_reflect_in_x_and_y( pt);
      X_curve_2 rcv1 = traits->curve_reflect_in_x_and_y( cv1);
      X_curve_2 rcv2 = traits->curve_reflect_in_x_and_y( cv2);
      
      Point_2 rp1, rp2;
      bool result = traits->nearest_intersection_to_right(rcv1, rcv2, rpt, 
                                                          rp1, rp2);
    
      p1 = traits->point_reflect_in_x_and_y( rp1);
      p2 = traits->point_reflect_in_x_and_y( rp2);
      
      return result;
    }
      
    Comparison_result Compare_lexicographically_xy(const Point& p1, 
                                                   const Point& p2) const
    {
      return (traits->compare_xy(p1,p2));
    }

    // Imitate the previous behaviour of the two functions: that is, return
    // EQUAL in undefined cases.
    Comparison_result _curve_compare_at_x_right (const X_curve& cv1,
						 const X_curve& cv2,
						 const Point& p) const
    {
      // In case the point is not in the x-range of both curves or that
      // one of the curves is not defined at p's right, return EQUAL.
      Comparison_result res1 = traits->compare_x(traits->curve_source(cv1),
						 traits->curve_target(cv1));
      Comparison_result res2 = traits->compare_x(traits->curve_source(cv2),
						 traits->curve_target(cv2));
      const Point_2& left1 = (res1 == SMALLER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point_2& right1 = (res1 == LARGER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point_2& left2 = (res2 == SMALLER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);
      const Point_2& right2 = (res2 == LARGER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);

      if (traits->compare_x (right1, p) != LARGER)
	return (EQUAL);
      else if (traits->compare_x (left1, p) == LARGER)
	return (EQUAL);

      if (traits->compare_x (right2, p) != LARGER)
	return (EQUAL);
      else if (traits->compare_x (left2, p) == LARGER)
	return (EQUAL);

      // Compare using the traits function:
      Comparison_result res = traits->curve_compare_at_x (cv1, cv2, p);
      if (res != EQUAL)
	return (res);

      return (traits->curve_compare_at_x_right (cv1, cv2, p));
    }

    Comparison_result _curve_compare_at_x_left (const X_curve& cv1,
						const X_curve& cv2,
						const Point& p) const
    {
      // In case the point is not in the x-range of both curves or that
      // one of the curves is not defined at p's left, return EQUAL.
      Comparison_result res1 = traits->compare_x(traits->curve_source(cv1),
						 traits->curve_target(cv1));
      Comparison_result res2 = traits->compare_x(traits->curve_source(cv2),
						 traits->curve_target(cv2));
      const Point_2& left1 = (res1 == SMALLER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point_2& right1 = (res1 == LARGER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point_2& left2 = (res2 == SMALLER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);
      const Point_2& right2 = (res2 == LARGER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);

      if (traits->compare_x (left1, p) != SMALLER)
	return (EQUAL);
      else if (traits->compare_x (right1, p) == SMALLER)
	return (EQUAL);

      if (traits->compare_x (left2, p) != SMALLER)
	return (EQUAL);
      else if (traits->compare_x (right2, p) == SMALLER)
	return (EQUAL);

      // Compare using the traits function:
      Comparison_result res = traits->curve_compare_at_x (cv1, cv2, p);
      if (res != EQUAL)
	return (res);

      return (traits->curve_compare_at_x_left (cv1, cv2, p));
    }

    // hold the left most intersecting point.
    Point_plus               intersect_p;
    Curve_node_container     curves;   
    Traits                  *traits;
  };
  
  template <class Point>
  class less_point_xy {
  public:
    typedef SweepLineTraits_2                   Traits;
    
    less_point_xy(Traits *traits_) : traits(traits_) {}
    
    inline  bool operator()(const Point& p1, const Point& p2) const
    { 
      return (traits->compare_xy(p1,p2) == SMALLER);
    }

  private:
    Traits *traits;
  };


  // A predicate ordering two Curve nodes.
  template <class _Curve_node>
  class less_curve_xy {
  public:
    typedef SweepLineTraits_2                       Traits;
    typedef typename Traits::X_curve                X_curve; 
    typedef typename Traits::Point                  Point;
    
    less_curve_xy(Traits *traits_) : traits(traits_) {}
    
    /*!
      Returns true if the first curve is "less" than the second curve.
      If both curves are "the same", it returns the one with a smaller id.
    */
    inline  bool operator()(const _Curve_node& cv1, 
                            const _Curve_node& cv2) const
    {
      Comparison_result result;
  
      const Point& ref_point = rightmost(cv1.get_rightmost_point().point(),
                                         cv2.get_rightmost_point().point());
  
      // if both curves are vartical, we look at the ids
      if (traits->curve_is_vertical(cv2.get_curve()) &&
        (traits->curve_is_vertical(cv1.get_curve())))
        return ( cv1.get_curve().id() < cv2.get_curve().id() );
  
      // if the rightmost points are different, we compare the curves relative
      // to the righmost of the two
      if ( cv1.get_rightmost_point() != cv2.get_rightmost_point() ) {
        result = curve_node_compare_at_x(cv1, cv2, ref_point);
      } else
      { 
        // both curves enamting the same point.
        if ( traits->curve_is_vertical(cv1.get_curve()) )  
          // if one of the curves enamating from the point is vertical
          // then it will have larger value on the status.
          return false;
        else if ( traits->curve_is_vertical(cv2.get_curve()) )
          return true;
    
        result = _curve_compare_at_x_right (cv1.get_curve(), 
					    cv2.get_curve(), 
					    ref_point);
        if (result == EQUAL)
          result = traits->curve_compare_at_x (cv1.get_curve(), 
                                               cv2.get_curve(), 
                                               ref_point);
      }
  
      if (result == SMALLER){
        return true;
      }
      else if (result == LARGER){
        return false;
      }
      else { // equal - the curves are overlapping.
        return ( cv1.get_curve().id() < cv2.get_curve().id() );  
      }
    }
      
  private:
    Comparison_result curve_node_compare_at_x(const _Curve_node &cv1, 
                                              const _Curve_node &cv2, 
                                              const Point &q) const
    {
      Comparison_result result;
      bool update_first_cv = false, update_second_cv = false;
      X_curve_plus first_cv , second_cv;  // making an optimization.
        
      // taking care the edge case of vertical segments: if one is
      // vertical, then the order is set according the 'rightmost'
      // point (so far) of the vertical curve comparing the point that
      // a vertical line hits the other curve.
      if ( traits->curve_is_vertical(cv1.get_curve()) ) {
        
        // if the curve is vertical and the point is its target we
        // refer that curve as a point.
        if (traits->point_is_same(traits->curve_target(cv1.get_curve()), 
                                  cv1.get_rightmost_point().point()))
	{
	  if (! traits->curve_is_in_x_range
	      (cv2.get_curve(), 
	       traits->curve_target(cv1.get_curve())))
	    return (EQUAL);
	  
	  Comparison_result cres = traits->curve_get_point_status
	    (cv2.get_curve(), 
	     traits->curve_target(cv1.get_curve()));

          if (cres == LARGER || cres == EQUAL)
            // if the target is on cv2 its means that it tangent to
            // cv2 from below.
            return  SMALLER;
          else
            return  LARGER;
        }
          
        X_curve tmp_cv;
        if (!traits->point_is_same(traits->curve_source(cv1.get_curve()), q)) {
          traits->curve_split(cv1.get_curve(), tmp_cv, first_cv, q);

          update_first_cv = true;
        }
      }
      
      if ( traits->curve_is_vertical(cv2.get_curve()) ){
        
        // if the curve is vertical and the point is its target we
        // refer that curve as a point.
        if (traits->point_is_same(traits->curve_target(cv2.get_curve()),
                                  cv2.get_rightmost_point().point()))
	{
	  if (! traits->curve_is_in_x_range
	      (cv1.get_curve(), 
	       traits->curve_target(cv2.get_curve())))
	    return (EQUAL);

          Comparison_result cres = traits->curve_get_point_status
	    (cv1.get_curve(), 
	     traits->curve_target(cv2.get_curve()));

          // if the target is on cv1 its means that it tangent to cv1 
          // from below.
	  if (cres == LARGER || cres == EQUAL)
            return  LARGER;
          else
            return  SMALLER;
        }
        
        X_curve tmp_cv;
        if (!traits->point_is_same(traits->curve_source(cv2.get_curve()), q)) {
          traits->curve_split(cv2.get_curve(), tmp_cv, second_cv, q);

          update_second_cv = true;
        }
      }    
      
      // making this four cases in order to make an optimization for not 
      // copying to first_cv and second_cv the original curves is not needed.
      if (!update_first_cv && !update_second_cv) 
        result =
          _curve_compare_at_x_right
	      (cv1.get_curve(), cv2.get_curve(),  
	       // making an optimization attemp:
	       rightmost(cv1.get_rightmost_point().point(),
			 cv2.get_rightmost_point().point()));
      
      else if (update_first_cv && !update_second_cv)
        result =
          _curve_compare_at_x_right
             (first_cv, cv2.get_curve(), 
	      rightmost(cv1.get_rightmost_point().point(),
			cv2.get_rightmost_point().point()));
      else if (!update_first_cv && update_second_cv)
        result =
          _curve_compare_at_x_right
	      (cv1.get_curve(), second_cv, 
	       rightmost(cv1.get_rightmost_point().point(),
			 cv2.get_rightmost_point().point()));
      else               // update_first_cv && update_second_cv is true.
        result =
          _curve_compare_at_x_right
	      (first_cv, second_cv, 
	       rightmost(cv1.get_rightmost_point().point(),
			 cv2.get_rightmost_point().point()));
      
      if (result == EQUAL){
        if (!update_first_cv && !update_second_cv)
          result =
            traits->
              curve_compare_at_x(cv1.get_curve(), cv2.get_curve(),
                                 rightmost(cv1.get_rightmost_point().point(), 
                                           cv2.get_rightmost_point().point()));
        else if (update_first_cv && !update_second_cv)
          result =
            traits->
              curve_compare_at_x(first_cv, cv2.get_curve(), 
                                 rightmost(cv1.get_rightmost_point().point(), 
                                           cv2.get_rightmost_point().point()));
        else if (!update_first_cv && update_second_cv)
          result =
            traits->
              curve_compare_at_x(cv1.get_curve(),second_cv, 
                                 rightmost(cv1.get_rightmost_point().point(), 
                                           cv2.get_rightmost_point().point()));
        else           // update_first_cv && update_second_cv is true.
          result =
            traits->
              curve_compare_at_x(first_cv, second_cv,
                                 rightmost(cv1.get_rightmost_point().point(), 
                                           cv2.get_rightmost_point().point()));
      }
    
      
      // if one of the curves is vertical - EQUAL means that the other curve 
      // is between the source and target point of the vertical curve, and 
      // for our definition - it means that the verical curve comes first.
      if ( result == EQUAL && 
         traits->curve_is_vertical(update_first_cv ? first_cv: 
                                   cv1.get_curve()) )
      {
	if (! traits->curve_is_in_x_range
	    (update_second_cv ? second_cv : cv2.get_curve(),
	     traits->curve_source(update_first_cv ? 
				  first_cv : cv1.get_curve())))
	  result = SMALLER;

        // if first_cv is vertical and its source tangent to second_cv - 
        // it means that first_cv is above second_cv.
	else if (traits->curve_get_point_status
		 (update_second_cv ? second_cv : cv2.get_curve(),
		  traits->curve_source(update_first_cv ? 
				       first_cv : cv1.get_curve())) == EQUAL)
          result = LARGER;
        else
          result = SMALLER;
      }
      else if (result == EQUAL && 
               traits->curve_is_vertical(update_second_cv? second_cv: 
                                         cv2.get_curve()))
      {
	if (! traits->curve_is_in_x_range
	    (update_first_cv ? first_cv : cv1.get_curve(), 
	     traits->curve_source(update_second_cv ? 
				  second_cv : cv2.get_curve())))
	  result = LARGER;

        // if second_cv is vertical and its source tangent to first_cv - 
        // it means that second_cv is above first_cv.
	else if (traits->curve_get_point_status
		 (update_first_cv ? first_cv: cv1.get_curve(), 
		  traits->curve_source(update_second_cv ? 
				       second_cv : cv2.get_curve())) == EQUAL)
          result = SMALLER;
        else
          result = LARGER;
      }
      return result;
    }
      
    const Point& rightmost(const Point &p1, const Point &p2) const
    { 
      Comparison_result res = traits->compare_xy(p1,p2);
      
      return ((res == SMALLER) ? p2 : p1);
    }

    // Imitate the previous behaviour of the function: that is, return
    // EQUAL in undefined cases.
    Comparison_result _curve_compare_at_x_right (const X_curve& cv1,
						 const X_curve& cv2,
						 const Point& p) const
    {
      // In case the point is not in the x-range of both curves or that
      // one of the curves is not defined at p's right, return EQUAL.
      Comparison_result res1 = traits->compare_x(traits->curve_source(cv1),
						 traits->curve_target(cv1));
      Comparison_result res2 = traits->compare_x(traits->curve_source(cv2),
						 traits->curve_target(cv2));
      const Point& left1 = (res1 == SMALLER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point& right1 = (res1 == LARGER) ? traits->curve_source(cv1) :
	traits->curve_target(cv1);
      const Point& left2 = (res2 == SMALLER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);
      const Point& right2 = (res2 == LARGER) ? traits->curve_source(cv2) :
	traits->curve_target(cv2);

      if (traits->compare_x (right1, p) != LARGER)
	return (EQUAL);
      else if (traits->compare_x (left1, p) == LARGER)
	return (EQUAL);

      if (traits->compare_x (right2, p) != LARGER)
	return (EQUAL);
      else if (traits->compare_x (left2, p) == LARGER)
	return (EQUAL);

      // Compare using the traits function:
      Comparison_result res = traits->curve_compare_at_x (cv1, cv2, p);
      if (res != EQUAL)
	return (res);

      return (traits->curve_compare_at_x_right (cv1, cv2, p));
    }
      
    Traits  *traits;
  };

public:
  typedef CurveInputIterator                    Curve_iterator;
  typedef SweepLineTraits_2                     Traits;
  typedef Point_plus_                           Point_plus;
  typedef typename  Traits::Point               Point;
  
  typedef X_curve_plus_                         X_curve_plus;
  typedef typename Curve_node::Points_iterator  Points_iterator;  
  typedef typename Curve_node::Points_const_iterator     
	                                        Points_const_iterator;
  typedef typename Intersection_point_node::Curve_node_iterator     
                                                Curve_node_iterator;
  typedef typename Intersection_point_node::Curve_node_const_iterator  
                                                Curve_node_const_iterator;

protected:
  typedef  std::pair<const Point,Point_plus>    pair_Point_Point_plus;
  typedef  std::pair<const Point,Intersection_point_node>  
                                        pair_Point_Intersection_point_node; 
  typedef  std::pair<const Curve_node,X_curve_plus>   
                                                pair_Curve_node_X_curve_plus;

  typedef less_point_xy<Point>                  Less_xy;
  typedef less_curve_xy<Curve_node>             Less_yx;
  
public:
  typedef std::map<Point, Point_plus, Less_xy>  Vertices_points_plus;

  typedef std::map<Point,Intersection_point_node, Less_xy>
                                                Event_queue;

  typedef std::set<Curve_node, Less_yx>         Status_line; 
  
  typedef std::list<X_curve>                    X_curve_list;
  typedef typename X_curve_list::iterator       X_curve_list_iterator;

  typedef typename Event_queue::value_type      Event_queue_value_type;
  typedef typename Status_line::value_type      Status_line_value_type;

  typedef typename Event_queue::iterator        Event_queue_iterator;
  typedef typename Status_line::iterator        Status_line_iterator;
  
  typedef typename std::list<Curve_node>::iterator
                                                list_Curve_node_iterator;
  
  Sweep_curves_base_2() : 
    traits(new Traits), use_delete_traits(true), intersection_exist_(false) {}
  
  Sweep_curves_base_2(Traits *traits_) : 
    traits(traits_), use_delete_traits(false), intersection_exist_(false) {}
  
  ~Sweep_curves_base_2()
  {
    if (use_delete_traits)
      delete traits;
  }

protected:
  
  void reset() { intersection_exist_ = false; }

  // change it back to the original one.
  bool  handle_one_event (Event_queue & event_queue, 
                          Status_line & status, 
                          const Point & event_point, 
                          Intersection_point_node & point_node)
  {
    bool event_terminated = true;
    bool point_node_first_on_status = false, point_node_last_on_status = false;
    
    // a container to hold the point_node curve nodes ordered as they
    // should on status.
    typedef std::set<Curve_node,Less_yx>         Local_status_line; 
    typedef typename Local_status_line::iterator Local_status_line_iterator; 
    Less_yx            pred(traits);
    Local_status_line local_status(pred);   // reserve the size of the
                                            // point_node enamating curves.
    
    
    Status_line_iterator hint = status.end();
    Curve_node_iterator cv_iter = point_node.curves_begin();
    if (cv_iter != point_node.curves_begin())
      --cv_iter;
    for ( ; hint == status.end() && cv_iter != point_node.curves_end(); 
          ++cv_iter)
      hint = std::find(status.begin(), status.end(), *cv_iter);
    
    if (hint == status.begin())
      point_node_first_on_status = true;
    
    if (hint != status.end()){
      ++hint;  // getting the first curve node on status above the
               // highest *cv_iter of curve_node.
      if (hint == status.end())
        point_node_last_on_status = true;
    }
    
    //bool hint_in_point_node = false;
    while (hint != status.end()){
      for (cv_iter = point_node.curves_begin() ; 
           cv_iter != point_node.curves_end(); ++cv_iter)
        if (*hint == *cv_iter) {
	  break;
	}
    
      if (cv_iter != point_node.curves_end())
        ++hint;
      else
        break;
    }
    
    Local_status_line_iterator local_status_hint = local_status.end();
    for (cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); ++cv_iter)
    {
      
      // trying to use general find in order to force it to use
      // operator == of Curve node. We only care here if *cv_iter is
      // on the status, we don't care where it should be put on status
      // according less_curve_xy predicate.  Since we work with handles only
      // address comparing will be done.
      Status_line_iterator curr_cv_node = std::find(status.begin(), 
                                                    status.end(), 
                                                    *cv_iter);
      
      // the event point is not the right point of the curve - 
      // insert curve to status if it's not already there.  
      if (curr_cv_node != status.end()){
        // remove cv_iter from status and reinsert it with the new point.
        status.erase(curr_cv_node);
      }
      
      CGAL_expensive_postcondition_code(is_valid(status));
      
      // reinserting curve only if it's still on status.
      if (!traits->point_is_same(event_point,
                                 cv_iter->get_rightmost_point().point()))
        cv_iter->push_event_point(point_node.get_point());
      
      if (local_status_hint != local_status.end()){
          ++local_status_hint;
      }
      
      if (Compare_lexicographically_xy(
                   event_point, 
                   traits->curve_target(cv_iter->get_curve()) ) == 
          CGAL::SMALLER){
        local_status_hint = local_status.insert(local_status_hint, *cv_iter);
        
        //if (local_status_hint == local_status.end()) cout << "ERROR
        //- returned value of insert to local status is
        //local_status.end()"<<endl;
      }
    }
    
    if (local_status.empty()) {
      if (status.empty())
        return event_terminated;
      
      if (!(point_node_first_on_status || point_node_last_on_status) 
          && hint != status.begin())  
        // if the neighbor of the curves enamting from point node is
        // not the last one or the first one on status
        if (check_status_neighbors_intersections(event_queue, 
                                                 status, 
                                                 --hint, 
                                                 event_point)) 
          event_terminated = false;
        
      
      return event_terminated;
    }
    
    Status_line_iterator new_cv_node_begin = status.insert(
                                               hint, 
                                               *(local_status.begin()) );
    
    hint = new_cv_node_begin;
    Local_status_line_iterator local_status_iter = local_status.begin();

    for (++local_status_iter; 
         local_status_iter != local_status.end(); ++local_status_iter){
      
       if (hint != status.end()){
         ++hint;
       }
       hint = status.insert(hint, *local_status_iter);
    }

    Status_line_iterator new_cv_node_end = new_cv_node_begin;
    std::advance(new_cv_node_end, local_status.size());
    Status_line_iterator new_cv_node_it;
    for (new_cv_node_it = new_cv_node_begin; 
         new_cv_node_it != new_cv_node_end; ++new_cv_node_it){

      if (check_status_neighbors_intersections(event_queue, status, 
                                               new_cv_node_it, event_point))
        // Edge case of tangency in the event point, 
        // if it is this event will be taked cared again.
        event_terminated = false;
      
      
      // only if it's the first one, we will also check the preseding
      // curve node adjancy.
      if (new_cv_node_it ==  new_cv_node_begin){
        if (new_cv_node_it != status.begin()){
          --new_cv_node_it;
          
          if (check_status_neighbors_intersections(event_queue, 
                                                   status, 
                                                   new_cv_node_it, 
                                                   event_point))
            // Edge case of tangency in the event point, if it is - 
            // this event will be taked cared again.
            event_terminated = false;
          ++new_cv_node_it;
        }
      }
    }
    
    return event_terminated;
  }
 
  
// Handling overlapping curves.  
// On each overlapping group, we remove
// iteratively each curve and check for new events after the removing.
// when finish, we reinsert to the status all the overlapping removed
// curves.
  void  handle_overlapping_curves(Event_queue& event_queue,
                                  Status_line& status,  
                                  const Point &event_point,
                                  Intersection_point_node& point_node)
  {
    // bool event_overlap_terminated = true;
    
    for (Curve_node_iterator cv_iter = point_node.curves_begin(); 
         cv_iter != point_node.curves_end(); ++cv_iter) {

      Status_line_iterator curr_cv_node = std::find(status.begin(), 
                                                    status.end(), 
                                                    *cv_iter);
      if (curr_cv_node != status.end()){
        if (curr_cv_node != status.begin()){
          std::list<Curve_node>  overlapping_curves;
          
          Status_line_iterator lower =  --curr_cv_node;
          ++curr_cv_node;
          for ( ;curr_cv_node != status.begin() && 
                  traits->curves_overlap(lower->get_curve(), 
                                         curr_cv_node->get_curve()); --lower){
            
            Point p1, p2;
            traits->nearest_intersection_to_right (curr_cv_node->get_curve(), 
                                                  lower->get_curve(), 
                                                  event_point, p1, p2);
            
            if (is_left(event_point, p1) || is_right(event_point, p2)) 
              // means that the overlapping does not intersects the 
              // status line.
              break;
            
            overlapping_curves.push_back(*curr_cv_node);
            status.erase(curr_cv_node);
            
            curr_cv_node = lower;
            Curve_node  cv_node = *curr_cv_node;
            
            // now taking care of the events created by the overlapping curve
            
            Point xp;
            if (check_status_neighbors_intersections(event_queue, 
                                                     status, 
                                                     curr_cv_node, 
                                                     event_point))
              // Edge case of tangency in the event point, 
              // if it is this event will be taked cared again.
              // event_overlap_terminated = false;
            
            if (curr_cv_node != status.begin()){
              --curr_cv_node;
              if (check_status_neighbors_intersections(event_queue, 
                                                       status, 
                                                       curr_cv_node, 
                                                       event_point))
                // Edge case of tangency in the event point, 
                // if it is - this event will be taked cared again.
                // event_overlap_terminated = false;
              ++curr_cv_node;
            }
            
            if (curr_cv_node != status.end() && 
                Compare_lexicographically_xy (
                         event_point, 
                         traits->curve_target(curr_cv_node->get_curve()) ) ==
                CGAL::EQUAL){
              
              Status_line_iterator prev_cv_node;
              bool first = true;
              // hold the (lower) neighbor element of the current.
              if (curr_cv_node != status.begin()){
                prev_cv_node = --curr_cv_node;
                ++curr_cv_node;
                first = false;
              }
              status.erase(curr_cv_node);
              
              CGAL_expensive_postcondition_code(is_valid(status));
              
              if (!first){
                if (check_status_neighbors_intersections(event_queue, 
                                                         status, 
                                                         prev_cv_node, 
                                                         event_point))
                  // Edge case of tangency in the event point, if it
                  // is this event will be taked cared again.
                  // event_overlap_terminated = false;
                    ;
              }
            } 
          }
          // reinsert to the status line all the overlapping removed curves.
          for (list_Curve_node_iterator ovlp_iter=overlapping_curves.begin(); 
               ovlp_iter != overlapping_curves.end(); ++ovlp_iter)
            
            status.insert(Status_line_value_type(*ovlp_iter));
        }
      }
      
      curr_cv_node = std::find(status.begin(), status.end(), *cv_iter);
      if (curr_cv_node != status.end()){
        std::list<Curve_node>  overlapping_curves;
        
        Status_line_iterator upper =  ++curr_cv_node;
        --curr_cv_node;
        
        
        for ( ;upper != status.end() && curr_cv_node != status.end() && 
                traits->curves_overlap(curr_cv_node->get_curve(), 
                                      upper->get_curve()); ++upper){
          
          Point p1, p2;
          traits->nearest_intersection_to_right (curr_cv_node->get_curve(), 
                                                upper->get_curve(), 
                                                event_point , p1, p2);
          
          // means that the overlapping does not intersects the status line.
          if (is_left(event_point, p1) || is_right(event_point, p2)) 
            break;
          
          overlapping_curves.push_back(*curr_cv_node);
          status.erase(curr_cv_node);
          
          curr_cv_node = upper;
          Curve_node  cv_node = *curr_cv_node;
          
          // now taking care of the events created by the overlapping curve.
          Point xp;
          if (check_status_neighbors_intersections(event_queue, status, 
                                                   curr_cv_node, 
                                                   event_point))
            // Edge case of tangency in the event point, 
            // if it is this event will be taked cared again.
            // event_overlap_terminated = false;
          
          if (curr_cv_node != status.begin()){
            --curr_cv_node;
            if (check_status_neighbors_intersections(event_queue, 
                                                     status, 
                                                     curr_cv_node, 
                                                     event_point))
              // Edge case of tangency in the event point, 
              // if it is - this event will be taked cared again.
              // event_overlap_terminated = false;
            
            ++curr_cv_node;
          }
          
          if (curr_cv_node != status.end() && 
              Compare_lexicographically_xy (
                       event_point, 
                       traits->curve_target(curr_cv_node->get_curve()) ) == 
              CGAL::EQUAL){
            
            Status_line_iterator prev_cv_node;
            bool first = true;
            // hold the (lower) neighbor element of the current.
            if (curr_cv_node != status.begin()){
              prev_cv_node = --curr_cv_node;
              ++curr_cv_node;
              first = false;
            }
            status.erase(curr_cv_node);
            
            CGAL_expensive_postcondition_code(is_valid(status));
            if (!first){
              //cout<<"checking neighbors after deletion\n";
              if (check_status_neighbors_intersections(event_queue, 
                                                       status, 
                                                       prev_cv_node, 
                                                       event_point))
                // Edge case of tangency in the event point, 
                // if it is this event will be taked cared again.
                // event_overlap_terminated = false;
                  ;
            }
          } 
        }
        // reinsert to the status line all the overlapping removed curves.
        for (list_Curve_node_iterator  ovlp_iter = 
               overlapping_curves.begin(); 
             ovlp_iter != overlapping_curves.end();  ++ovlp_iter)
          status.insert(Status_line_value_type(*ovlp_iter));
      }
    }
  }
  
  
  bool  point_is_on_curve_interior(const Point & p, const X_curve & cv) {
    return (!traits->point_is_same(p, traits->curve_source(cv)) ||
            !traits->point_is_same(p, traits->curve_target(cv)));
  }
  
  // check whether two curves are tangent on a given point.
  bool  curves_tangent_at_point(const X_curve& cv1, 
                                const X_curve& cv2, 
                                const Point& point)
  { 
    return (traits->curve_is_in_x_range(cv1, point) &&
	    traits->curve_get_point_status(cv1, point) == EQUAL && 
	    traits->curve_is_in_x_range(cv2, point) &&
            traits->curve_get_point_status(cv2, point) == EQUAL && 
            (point_is_on_curve_interior(point, cv1) || 
             point_is_on_curve_interior(point, cv2) ) );
    
  }

  bool check_status_neighbors_intersections(Event_queue& event_queue, 
                                            Status_line& status,  
                                            Status_line_iterator lower_neighbor
                                            ,const Point &event_point)
  { 
    const Curve_node& cv1 = *lower_neighbor;
    
    Status_line_iterator next_neighbor = ++lower_neighbor;
    if (next_neighbor == status.end())
      return false;
    
    --lower_neighbor;
    
    const Curve_node&  cv2 = *next_neighbor;
    
    // in each node - checking intersections between two 
    // adjacent curves and updating the event queue if needed.  
    const Point &xp1 = event_point;
    Point xp2, xp3;
    
    bool  curves_tangent = curves_tangent_at_point(cv1.get_curve(), 
                                                   cv2.get_curve(), 
                                                   event_point);
    bool curves_intersect = 
      traits->nearest_intersection_to_right(cv1.get_curve(), cv2.get_curve(),
                                            event_point, xp2, xp3);
      
    if (curves_tangent || curves_intersect) {   
      
      // checking if the curves are already participate within the
      // intersection node, if at least one is not - we have a new
      // event here or a new curve particitating on it and the
      // function will return true.
     
      // handling overlapping.
      if (curves_intersect && (!traits->point_is_same(xp2, xp3))) {
        
        Event_queue_iterator  xp_event = event_queue.find(xp3);
        bool xp3_cv1_in_queue = false,  xp3_cv2_in_queue = false;
        if (xp_event == event_queue.end())
          xp_event =
            event_queue.insert(Event_queue_value_type
                               (xp3, 
                                Intersection_point_node(cv1, cv2, xp3, traits)
                                )).first;
        else{
          // have to check whether the event is a new event. 
          // (we might calculated this point before).
          for ( Curve_node_iterator cv_iter = 
                  xp_event->second.curves_begin(); 
                cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (cv_iter->get_curve() == cv1.get_curve()) {
	      xp3_cv1_in_queue = true;
	    }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              xp3_cv2_in_queue = true;
            }
          }
          
          if (!xp3_cv1_in_queue && !xp3_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, 
                                                           cv2, 
                                                           Point_plus(xp3), 
                                                           traits)); 
      
          else if (!xp3_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, 
                                                           Point_plus(xp3), 
                                                           traits));
          else if (!xp3_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, 
                                                           Point_plus(xp3), 
                                                           traits));
        }
      }

      // handling intersection.
      bool xp2_cv1_in_queue = false,  xp2_cv2_in_queue = false;
      if ( curves_intersect )
      {
        Event_queue_iterator  xp_event = event_queue.find(xp2);
        if (xp_event == event_queue.end())
          xp_event = event_queue.insert(
	                 Event_queue_value_type(xp2, 
					 Intersection_point_node(cv1, 
								cv2, 
								xp2, 
								traits))).first;
        else
	{
          // have to check whether the event is a new event. (we might
          // calculated this point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); 
                cv_iter != xp_event->second.curves_end(); cv_iter++){
            if (cv_iter->get_curve() == cv1.get_curve()) {
	      xp2_cv1_in_queue = true;
	    }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              xp2_cv2_in_queue = true;
            }
          }
          
          if (!xp2_cv1_in_queue && !xp2_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, cv2, 
                                                           Point_plus(xp2), 
                                                           traits));  
        
          else if (!xp2_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, 
                                                           Point_plus(xp2), 
                                                           traits));
          else if (!xp2_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, 
                                                           Point_plus(xp2), 
                                                           traits));
        }
      }
      
      bool xp1_cv1_in_queue = false,  xp1_cv2_in_queue = false; 
      // if cv1 and cv2 have common edge point - we do not consider it 
      // as an intersection point. 
      if ( curves_tangent ) {
        Event_queue_iterator  xp_event = event_queue.find(xp1);
        if (xp_event == event_queue.end())
          xp_event=event_queue.insert(Event_queue_value_type(
                                       xp1, 
                                       Intersection_point_node(cv1, 
                                                               cv2, 
                                                               xp1, 
                                                               traits))).first;
        else{
          // have to check whether the event is a new event. 
          // (we might calculated this point before).
          for ( Curve_node_iterator cv_iter = xp_event->second.curves_begin(); 
                cv_iter != xp_event->second.curves_end(); ++cv_iter){
            if (cv_iter->get_curve() == cv1.get_curve())  {
	      xp1_cv1_in_queue = true;
	    }
            if (cv_iter->get_curve() ==  cv2.get_curve()) {
              xp1_cv2_in_queue = true;
            }
          }
          
          if (!xp1_cv1_in_queue && !xp1_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, 
                                                           cv2, 
                                                           Point_plus(xp1), 
                                                           traits));
          else if (!xp1_cv1_in_queue)
            xp_event->second.merge(Intersection_point_node(cv1, 
                                                           Point_plus(xp1), 
                                                           traits));
          else if (!xp1_cv2_in_queue)
            xp_event->second.merge(Intersection_point_node(cv2, 
                                                           Point_plus(xp1), 
                                                           traits));
        }
      }  
    
      if (curves_intersect || curves_tangent)
        intersection_exist_ = true;
      
      if (curves_intersect && !curves_tangent)
        return ((!xp2_cv1_in_queue || !xp2_cv2_in_queue) &&
                traits->point_is_same(xp2, event_point));
      else if (curves_tangent)
        return ((!xp1_cv1_in_queue || !xp1_cv2_in_queue) &&
                traits->point_is_same(xp1, event_point));
    }
    
    return false;
  }

  // defining a lexicographically compare function here in order to
  // make a usage in traits.
  Comparison_result Compare_lexicographically_xy(const Point & p1, 
                                                 const Point & p2) const
  {
    return (traits->compare_xy(p1,p2));
  }

#if defined SWEEP_DEBUG_MODE
  //------------------------------------------------------- debuging functions
  //-----------------------------------------
  bool  is_valid(Status_line& status) {

    std::list<Curve_node> curve_nodes;
    std::list<X_curve>    subcurves;
    
    for (typename Status_line::const_iterator iter = status.begin(); 
         iter != status.end(); iter++)
      curve_nodes.push_back(*iter);
    
    get_subcurves(curve_nodes, subcurves);
    
    return (!do_intersect_subcurves(subcurves));
  }
  
  void  get_subcurves(const std::list<Curve_node >& curves,  
                      std::list<X_curve>& subcurves)
  {
    for (std::list<Curve_node>::const_iterator cv_iter = curves.begin();
	 cv_iter != curves.end(); cv_iter++){
      X_curve cv = cv_iter->get_curve(), left_cv, 
        right_cv = cv_iter->get_curve();
      for (typename Curve_node::Points_const_iterator points_iter = 
             cv_iter->points_begin(); 
           points_iter != cv_iter->points_end(); points_iter++){
        // make surve the splitting is not at the edge points.
        if (points_iter == cv_iter->points_begin())  
          continue;
        
        if (traits->
            point_is_same(*points_iter,
                          rightmost(
                            traits->curve_source(cv_iter->get_curve()),
                            traits->curve_target(cv_iter->get_curve()))))
          left_cv = right_cv;
        else
          traits->curve_split(cv, left_cv, right_cv, points_iter->point());
        
        subcurves.push_back(left_cv);
        
        cv = right_cv;
      }
    }
  }
  
  bool  do_intersect_subcurves(const std::list<X_curve>& subcurves)
  {
    for (std::list<X_curve>::const_iterator scv_iter1 = subcurves.begin(); 
	 scv_iter1 != subcurves.end(); ++scv_iter1){
      for (std::list<X_curve>::const_iterator  scv_iter2 = scv_iter1; 
           scv_iter2 != subcurves.end(); ++scv_iter2){
        if (scv_iter2 ==  scv_iter1)
          continue;
        
        Point  ref_point1, ref_point2, xp1, xp2;

        ref_point1 = leftmost(traits->curve_source(*scv_iter1), 
                              traits->curve_target(*scv_iter1));
        ref_point2 = leftmost(traits->curve_source(*scv_iter2), 
                              traits->curve_target(*scv_iter2));
        if (!traits->point_is_same(ref_point1, ref_point2))
          ref_point1 = leftmost(ref_point1, ref_point2);

        if (traits->nearest_intersection_to_right(*scv_iter1, *scv_iter2, 
                                                 ref_point1, xp1, xp2)) {
          if (traits->point_is_same(xp1, ref_point1))
            xp1 = xp2;
          
          if (!((traits->point_is_same(xp1,
                                       traits->curve_source(*scv_iter1)) ||
                 traits->point_is_same(xp1,
                                       traits->curve_target(*scv_iter1))) && 
                (traits->point_is_same(xp1,
                                       traits->curve_source(*scv_iter2)) ||
                 traits->point_is_same(xp1,
                                       traits->curve_target(*scv_iter2)))))
          {
            return true;
          }
        }
      }
    }
    return false;
  }
#endif // SWEEP_DEBUG_MODE
    
  bool is_left(const Point &p1, const Point &p2) const 
  { return (Compare_lexicographically_xy(p1, p2) == SMALLER); }
  bool is_right(const Point &p1, const Point &p2) const 
  { return (Compare_lexicographically_xy(p1, p2) == LARGER); }
  bool is_lower(const Point &p1, const Point &p2) const 
  { return (Compare_lexicographically_yx(p1, p2) == SMALLER); }
  bool is_higher(const Point &p1, const Point &p2) const 
  { return (Compare_lexicographically_yx(p1, p2) == LARGER); }

  const Point& leftmost(const Point &p1, const Point &p2) const
  {
    if ( is_left(p1, p2) )
      return p1;
    return p2;
  }
    
  const Point& rightmost(const Point &p1, const Point &p2) const
  { 
    if ( is_right(p1, p2) )
      return p1;
    else
      return p2;
  }

#if defined SWEEP_DEBUG_MODE
  //-------------------------------- debuging functions.
  template <class Status>
  void  print_status(const Status &status);
  template <class _Curve_node>
  void  print_points_on_curves(const _Curve_node& cv);
#endif

  Traits   *traits;
  bool      use_delete_traits;
  bool      intersection_exist_;
};



////////////////////////////////////////////////////////////////////////
////            METHODS

#if defined SWEEP_DEBUG_MODE
// ---- debug functions -----
template <class CurveInputIterator, class SweepLineTraits_2, 
          class Point_plus_, class X_curve_plus_>
template <class Status>
void  
Sweep_curves_base_2<CurveInputIterator, SweepLineTraits_2,
                    Point_plus_, X_curve_plus_>::
print_status(const Status &status)
{

    cout<<"Curves on status are\n"; 
    for (typename Status_line::const_iterator 
           status_iter = status.begin(); status_iter != status.end(); 
         status_iter++){
      X_curve_plus cv1 = status_iter->get_curve();
      
      cout<<cv1<<" ";
      print_points_on_curves(*status_iter);
    }
    cout<<"\n";
  }


template <class CurveInputIterator, class SweepLineTraits_2, 
          class Point_plus_, class X_curve_plus_>
template <class _Curve_node>
void  
Sweep_curves_base_2<CurveInputIterator, SweepLineTraits_2,
                    Point_plus_, X_curve_plus_>::
print_points_on_curves(const _Curve_node& cv)
{
    cout<<"list of points is "<<std::endl;

    for (Points_const_iterator  p_iter = cv.points_begin(); 
         p_iter != cv.points_end(); p_iter++)
      cout<<p_iter->point()<<" ";
    
    cout<<std::endl;
  }
#endif // SWEEP_DEBUG_MODE

CGAL_END_NAMESPACE

#endif
