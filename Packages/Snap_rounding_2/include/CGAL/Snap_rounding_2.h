// ======================================================================
//
// Copyright (c) The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-24 $
// release_date  : $CGAL_Date: 2000/12/29 $
//
// file          : include/CGAL/Snap_rounding_2.h
// package       : arr (1.73)
// maintainer    : Eli Packer <elip@post.tau.ac.il>
// author(s)     : Eli Packer
// coordinator   : Dan Halperin Tel-Aviv University
//                 (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_SR_2_H
#define CGAL_SR_2_H

#include <CGAL/leda_rational.h> 

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Random.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/intersection_2.h>
//#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <list>
#include <set>
#include <CGAL/leda_real.h>
#include "Snap_rounding_kd_2.h"
#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE

enum Direction {SEG_UP_RIGHT,SEG_UP_LEFT,SEG_DOWN_RIGHT,SEG_DOWN_LEFT,
                SEG_UP,SEG_DOWN,SEG_LEFT,SEG_RIGHT,SEG_POINT_SEG};

template<class Rep_>
class Segment_data {

typedef Rep_                                Rep;
typedef typename Rep::FT                    NT;
typedef typename Rep::Segment_2             Segment_2;
typedef typename Rep::Point_2               Point_2;

private:
 Point_2 p;
 Point_2 q;

 Rep_   _gt;

public:
  Segment_data();
  Segment_data(const Point_2& p_inp,const Point_2& q_inp);
  Segment_2 segment() const {return(Segment_2(p,q));}
  Point_2 source() const {return(p);}
  Point_2 target() const {return(q);}
  inline void set_data(const Point_2& inp_p,const Point_2& inp_q);
  void determine_direction(Direction &seg_dir);
  bool equal(const Segment_2& s);
  Segment_data(const Segment_data &other);
};

template<class Rep_>
class Hot_Pixel {

typedef Rep_                       Rep;
typedef typename Rep::FT           NT;
typedef typename Rep::Segment_2    Segment_2;
typedef typename Rep::Point_2      Point_2;

private:
  // p is the center of the hot pixel
  Point_2 p;
  Point_2 p_left;
  Point_2 p_right;
  Point_2 p_down;
  Point_2 p_up;
  Rep_ _gt;
  NT pixel_size;
  Segment_2 *right_seg;
  Segment_2 *left_seg;
  Segment_2 *top_seg;
  Segment_2 *bot_seg;
  static Direction seg_dir;

public:
  template<class Out>
  void draw(Out &o) const;
  Hot_Pixel(const Point_2& inp_point,NT inp_pixel_size);
  ~Hot_Pixel();
  inline Point_2 get_center() const;
  inline Point_2 get_center(bool int_output) const;
  bool intersect_left(const Segment_2& seg,Direction seg_dir) const;
  bool intersect_right(const Segment_2& seg,Direction seg_dir) const;
  bool intersect_bot(const Segment_2& seg,Direction seg_dir) const;
  bool intersect_top(const Segment_2& seg,Direction seg_dir) const;
  bool intersect(Segment_data<Rep> &seg,Direction seg_dir) const;
  void set_direction(Direction inp_seg_dir) {seg_dir = inp_seg_dir;}
  Direction get_direction() const {return(seg_dir);}
};

template<class Rep_> Direction Hot_Pixel<Rep_>::seg_dir;

// a function for compare two hot pixels for the set of hot pixels
template<class Rep_>
struct hot_pixel_auclidian_cmp
{
  Rep_   _gt;

  bool operator()(const Hot_Pixel<Rep_> *h1,const Hot_Pixel<Rep_> *h2) const;
};

// a function for compare two hot pixels for the set of hot pixels a
// certain segment intersect
template<class Rep_>
struct hot_pixel_dir_cmp
{
  Rep_   _gt;

  bool operator ()(const Hot_Pixel<Rep_> *h1,const Hot_Pixel<Rep_> *h2);
};

template<class Rep_>
class Snap_rounding_2 {

typedef CGAL::Arr_segment_traits_2<Rep_ >            Traits;
typedef Rep_                                         Rep;
typedef typename Rep::FT                             NT;
typedef typename Traits::X_curve                     X_curve;
typedef typename Traits::Curve                       Curve;
typedef std::list<X_curve>                           CurveContainer;
typedef typename CurveContainer::iterator            CurveContainerIter;

public:
  friend class Segment_data<Rep>;
  friend class Hot_Pixel<Rep>;
  friend class hot_pixel_dir_cmp<Rep>;

  typedef typename Rep::Segment_2 Segment_2;
  typedef typename Rep::Point_2   Point_2;
  typedef std::list<Point_2> PointList;
  typedef typename std::list<std::list<Point_2> > Polylines_container;
  typedef typename Polylines_container::const_iterator Polyline_const_iterator;
  typedef typename Polylines_container::iterator Polyline_iterator;
  typedef typename std::list<typename Rep::Point_2> Points_container;
  typedef typename Points_container::const_iterator Point_const_iterator;

  void find_hot_pixels_and_create_kd_trees(
       NT pixel_size,
       unsigned int number_of_kd_trees,
       std::list<Segment_data<Rep> >& seg_list,
       Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> **mul_kd_tree);

  void iterate(
        Polylines_container& output_container,
        NT pixel_size,
        bool int_output,
        bool do_isr,
        std::list<Segment_data<Rep> >& seg_list,
        Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree);

private:
  Rep _gt;

  static const int default_number_of_kd_trees = 1;

  void find_intersected_hot_pixels(Segment_data<Rep> &seg,
                         std::set<Hot_Pixel<Rep> *,
                         hot_pixel_dir_cmp<Rep> > &hot_pixels_intersected_set,
                         int &number_of_intersections,
                         NT pixel_size,
                         Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree);
  void reroute_sr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                  &inp_hot_pixels_intersected_set,std::list<Point_2>
                  &seg_output,bool int_output);
  void reroute_isr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                   &inp_hot_pixels_intersected_set,std::list<Point_2>
                   &seg_output,int number_of_intersections,bool first_time,
                   NT pixel_size,bool int_output,
                   Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree);

  //  void list_copy(std::list<Point_2>& target,std::list<Point_2>& source);
};

// ctor
template<class Rep_>
Segment_data<Rep_>::Segment_data() {}
template<class Rep_>
Segment_data<Rep_>::Segment_data(const Point_2& p_inp,const Point_2& q_inp) :
                    p(p_inp), q(q_inp) {}

// cctor
template<class Rep_>
Segment_data<Rep_>::Segment_data(const Segment_data& other)
{
  p = other.p;
  q = other.q;
}

template<class Rep_>
inline void Segment_data<Rep_>::set_data(const Point_2& inp_p,
                                         const Point_2& inp_q)
{
  p = inp_p;
  q = inp_q;
}

template<class Rep_>
bool Segment_data<Rep_>::equal(const Segment_2& s)
{
  return(s.source() == p && s.target() == q);
}

template<class Rep_>
void Segment_data<Rep_>::determine_direction(Direction &seg_dir)
{
  Comparison_result cx = _gt.compare_x_2_object()(p,q);
  Comparison_result cy = _gt.compare_y_2_object()(p,q);

  if(cx == SMALLER) {
    if(cy == SMALLER)
      seg_dir = SEG_UP_RIGHT;
    else if(cy == EQUAL)
      seg_dir = SEG_RIGHT;
    else
      seg_dir = SEG_DOWN_RIGHT;
  } else if(cx == EQUAL) {
    if(cy == SMALLER)
      seg_dir = SEG_UP;
    else if(cy == EQUAL)
      seg_dir = SEG_POINT_SEG;
    else
      seg_dir = SEG_DOWN;
  } else {
    if(cy == SMALLER)
      seg_dir = SEG_UP_LEFT;
    else if(cy == EQUAL)
      seg_dir = SEG_LEFT;
    else
      seg_dir = SEG_DOWN_LEFT;
  }
}

template<class Rep_>
template<class Out> void Hot_Pixel<Rep_>::draw(Out &o) const
  {
    o << *right_seg;
    o << *left_seg;
    o << *top_seg;
    o << *bot_seg;
  }

// intersection pixel
template<class Rep_>
Hot_Pixel<Rep_>::Hot_Pixel(const Point_2& inp_point,NT inp_pixel_size) :
                           pixel_size(inp_pixel_size)
  {
    NT x,y;
    _gt.snap_2_object()(inp_point,pixel_size,x,y);

    p = Point_2(x,y);
    p_left = Point_2(x - pixel_size / 2.0,y);
    p_right = Point_2(x + pixel_size / 2.0,y);
    p_down = Point_2(x,y - pixel_size / 2.0);
    p_up = Point_2(x,y + pixel_size / 2.0);

    right_seg = new Segment_2(Point_2(x + pixel_size / 2.0,y -
                              pixel_size / 2.0),
                              Point_2(x + pixel_size / 2.0,y +
                              pixel_size / 2.0));
    left_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y -
                             pixel_size / 2.0),
                             Point_2(x - pixel_size / 2.0,y +
                             pixel_size / 2.0));
    top_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y +
                            pixel_size / 2.0),
                            Point_2(x + pixel_size / 2.0,y +
                            pixel_size / 2.0));
    bot_seg = new Segment_2(Point_2(x - pixel_size / 2.0,y -
                            pixel_size / 2.0),
                            Point_2(x + pixel_size / 2.0,y -
                            pixel_size / 2.0));
  }

template<class Rep_>
Hot_Pixel<Rep_>::~Hot_Pixel()
  {
    delete(right_seg);
    delete(left_seg);
    delete(top_seg);
    delete(bot_seg);
  }

template<class Rep_>
inline typename Hot_Pixel<Rep_>::Point_2 Hot_Pixel<Rep_>::get_center() const
    {return(p);}

template<class Rep_>
inline typename Hot_Pixel<Rep_>::Point_2 Hot_Pixel<Rep_>::get_center(
            bool int_output) const
  {
    if(int_output) {
      Point_2 out_p = _gt.integer_grid_point_2_object()(p,pixel_size);
      return(out_p);
    } else
      return(p);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_left(const Segment_2& seg,
        Direction seg_dir) const
  {
    Object result;    
    Point_2 p;
    Segment_2 s;

    result = intersection(seg,*left_seg);

    if(assign(p,result)) {
      Comparison_result c_p = _gt.compare_y_2_object()(p,p_up);
      Comparison_result c_target = _gt.compare_y_2_object()(seg.target(),p_up);
      Comparison_result c_source = _gt.compare_y_2_object()(seg.source(),p_up);

      return(c_p != EQUAL || seg_dir == SEG_UP_LEFT && c_source != EQUAL ||
             seg_dir == SEG_DOWN_RIGHT && c_target != EQUAL);
    } else if(assign(s,result))
      return(true);
    else
      return(false);
  }


template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_right(const Segment_2& seg,
         Direction seg_dir) const
  {
    Object result;    
    Point_2 p;
    Segment_2 s;
    
    result = intersection(seg,*right_seg);

    if(assign(p,result)) {
      // bottom right point was checked in intersect_bot
      Comparison_result c1 = _gt.compare_y_2_object()(p,p_up);
      Comparison_result c2 = _gt.compare_y_2_object()(p,p_down);
      Comparison_result c3 = _gt.compare_y_2_object()(seg.source(),p_up);
      Comparison_result c4 = _gt.compare_y_2_object()(seg.target(),p_up);

      if(c1 == EQUAL)
        return(seg_dir == SEG_UP_RIGHT && c3 != EQUAL ||
               seg_dir == SEG_DOWN_LEFT && c4 != EQUAL);
      else if(c2 == EQUAL)
        return(false);// was checked
      else {
        Comparison_result c_target =
             _gt.compare_x_2_object()(p_right,seg.target());
        Comparison_result c_source =
             _gt.compare_x_2_object()(p_right,seg.source());

        return((seg_dir == SEG_LEFT || seg_dir == SEG_DOWN_LEFT ||
                seg_dir == SEG_UP_LEFT) && c_target != EQUAL ||
                (seg_dir == SEG_RIGHT || seg_dir == SEG_DOWN_RIGHT ||
                seg_dir == SEG_UP_RIGHT) && c_source != EQUAL);
      }
    } else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_bot(const Segment_2& seg,
         Direction seg_dir) const
  {
    Object result;
    Point_2 p;
    Segment_2 s;

    result = intersection(seg,*bot_seg);

    if(assign(p,result)) {
      Comparison_result c_p = _gt.compare_x_2_object()(p,p_right);
      Comparison_result c_target = _gt.compare_x_2_object()
            (seg.target(),p_right);
      Comparison_result c_source = _gt.compare_x_2_object()
            (seg.source(),p_right);

      return(c_p != EQUAL || seg_dir == SEG_UP_LEFT &&
             c_target != EQUAL || seg_dir == SEG_DOWN_RIGHT &&
             c_source != EQUAL);
    } else if(assign(s,result))
      return(true);
    else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_top(const Segment_2& seg,
         Direction seg_dir) const
  {
    Object result;
    Point_2 p;
    Segment_2 s;
    
    result = intersection(seg,*top_seg);

    if(assign(p,result)) {
      // corner points was checked in intersect_bot
      Comparison_result c1 = _gt.compare_x_2_object()(p,p_left);
      Comparison_result c2 = _gt.compare_x_2_object()(p,p_right);
      Comparison_result c3 = _gt.compare_y_2_object()(seg.target(),p_up);
      Comparison_result c4 = _gt.compare_y_2_object()(seg.source(),p_up);

      if(c1 == EQUAL || c2 == EQUAL)
        return(false);// were checked
      else
        return((seg_dir == SEG_DOWN || seg_dir == SEG_DOWN_LEFT ||
               seg_dir == SEG_DOWN_RIGHT) && c3 != EQUAL ||
               (seg_dir == SEG_UP || seg_dir == SEG_UP_LEFT ||
               seg_dir == SEG_UP_RIGHT) && c4 != EQUAL);
    } else
    return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect(Segment_data<Rep_> &seg,
                                Direction seg_dir) const
  {
    Segment_2 s = seg.segment();

    return(intersect_bot(s,seg_dir) || intersect_left(s,seg_dir) ||
           intersect_right(s,seg_dir) || intersect_top(s,seg_dir));
  }

// a function for compare two hot pixels for the set of hot pixels
template<class Rep_>
bool hot_pixel_auclidian_cmp<Rep_>::operator()(const Hot_Pixel<Rep_> *h1,
     const Hot_Pixel<Rep_> *h2) const
  {
    Comparison_result cx = _gt.compare_x_2_object()(
        h1.get_center(),h2.get_center());
    Comparison_result cy = _gt.compare_y_2_object()(
        h1.get_center(),h2.get_center());

    return(cx == SMALLER ||
         cx == EQUAL && cy == SMALLER);
  }

// a function for compare two hot pixels for the set of hot pixels a certain
// segment intersect
template<class Rep_>
bool hot_pixel_dir_cmp<Rep_>::operator ()(const Hot_Pixel<Rep_> *h1,
     const Hot_Pixel<Rep_> *h2) 
{
  Comparison_result cx = _gt.compare_x_2_object()(
      h1->get_center(),h2->get_center());
  Comparison_result cy = _gt.compare_y_2_object()(
      h1->get_center(),h2->get_center());
  Direction seg_dir = h1->get_direction();

  return(
     // Point segment intersects only one pixel, thus ignored
    seg_dir == SEG_UP_RIGHT && (cx == SMALLER || cx == EQUAL &&
    cy == SMALLER) || seg_dir == SEG_UP_LEFT && (cx == LARGER || 
    cx == EQUAL && cy == SMALLER) || seg_dir == SEG_DOWN_RIGHT &&
    (cx == SMALLER || cx == EQUAL && cy == LARGER) ||
    seg_dir == SEG_DOWN_LEFT && (cx == LARGER || cx == EQUAL &&
    cy == LARGER) || seg_dir == SEG_UP && cy == SMALLER ||
    seg_dir == SEG_DOWN && cy == LARGER || seg_dir == SEG_LEFT &&
    cx == LARGER || seg_dir == SEG_RIGHT && cx == SMALLER);
}

template<class Rep_>
void Snap_rounding_2<Rep_>::find_hot_pixels_and_create_kd_trees(
       NT pixel_size,
       unsigned int number_of_kd_trees,
       std::list<Segment_data<Rep> >& seg_list,
       Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> **mul_kd_tree)
  {
    Hot_Pixel<Rep_> *hp;
    typename std::list<Segment_data<Rep_> >::iterator iter1;
    Object result;
    Point_2 p;
    std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> > hot_pixels_list;
    list<X_curve> segments;

    if(seg_list.empty())
      return;

    for(iter1 = seg_list.begin();iter1 != seg_list.end();++iter1)
      segments.push_back(X_curve(iter1->source(),iter1->target()));

    //    PM pm(new Pm_naive_point_location<PM>);
    // sweep_to_construct_planar_map(segments.begin(), segments.end(), pm);
    std::list<X_curve>  subcurves;

    /*    sweep_to_produce_planar_map_subcurves(segments.begin(), 
					  segments.end(),  
					  traits, 
					  subcurves);*/

    /*// get subcurves with overlapping ********** 
//    CGAL::Sweep_line_tight_2<CurveContainerIter, Traits, Event, SubCurve> sl;
    Sweep_line_2<CurveContainerIter, Traits> sl;
    sl.get_subcurves(segments.begin(), segments.end(),
    std::back_inserter(subcurves));*/

    // get intersection points (with endpoints)
    PointList mypointlist;
    //    CGAL::Sweep_line_tight_2<CurveContainerIter, Traits,
    //                             Event, SubCurve> sl;
    Sweep_line_2<CurveContainerIter, Traits> sl;
    sl.get_intersection_points(segments.begin(), segments.end(),
                             std::back_inserter(mypointlist));

    for(typename std::list<Point_2>::const_iterator
            v_iter = mypointlist.begin();
	v_iter != mypointlist.end();++v_iter) {
      hp = new Hot_Pixel<Rep_>(*v_iter,pixel_size);
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
				hp->get_center(),hp));
    }

    // create kd multiple tree
    // create simple_list from seg_list
    std::list<Segment_2> simple_seg_list;
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter)
      simple_seg_list.push_back(Segment_2(iter->source(),iter->target()));

    *mul_kd_tree = new Multiple_kd_tree<Rep,Hot_Pixel<Rep> *>(hot_pixels_list,
                  number_of_kd_trees,simple_seg_list);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::find_intersected_hot_pixels(Segment_data<Rep_>
                    &seg,
                    std::set<Hot_Pixel<Rep_> *,
                    hot_pixel_dir_cmp<Rep_> > &hot_pixels_intersected_set,
                    int &number_of_intersections,
                    NT pixel_size,
                    Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree)
  {
    typename std::list<Hot_Pixel<Rep_> *>::iterator iter;
    Direction seg_dir;

    hot_pixels_intersected_set.clear();
    seg.determine_direction(seg_dir);
    number_of_intersections = 0;

    std::list<Hot_Pixel<Rep_> *> hot_pixels_list;
    mul_kd_tree->get_intersecting_points(hot_pixels_list,
	   Segment_2(seg.segment()),pixel_size);

    for(iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
      if((*iter)->intersect(seg,seg_dir)) {
        (*iter)->set_direction(seg_dir);
        hot_pixels_intersected_set.insert(*iter);
      }
    }

    number_of_intersections = hot_pixels_intersected_set.size();
  }


template<class Rep_>
void Snap_rounding_2<Rep_>::reroute_sr(std::set<Hot_Pixel<Rep_> *,
     hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
     std::list<Point_2> &seg_output,bool int_output)
  {
    typename std::set<Hot_Pixel<Rep_> *,
    hot_pixel_dir_cmp<Rep_> >::iterator hot_pixel_iter =
    inp_hot_pixels_intersected_set.begin();
    ++hot_pixel_iter;

    while(hot_pixel_iter != inp_hot_pixels_intersected_set.end()) {
      seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
      ++hot_pixel_iter;
    }

  }

template<class Rep_>
void Snap_rounding_2<Rep_>::reroute_isr(
   std::set<Hot_Pixel<Rep_> *,
   hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
   std::list<Point_2> &seg_output,
   int number_of_intersections,
   bool first_time,
   NT pixel_size,
   bool int_output,
   Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree)
  {
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter,next_hot_pixel_iter,before_last_hot_pixel_iter;
    Segment_data<Rep_> seg;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;
    Direction seg_dir;

    if(number_of_intersections > 2 || first_time) {
      before_last_hot_pixel_iter = inp_hot_pixels_intersected_set.end();
      --before_last_hot_pixel_iter;

      for(hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
          hot_pixel_iter != before_last_hot_pixel_iter;++hot_pixel_iter) {
        next_hot_pixel_iter = hot_pixel_iter;
        ++next_hot_pixel_iter;
        seg.set_data((*hot_pixel_iter)->get_center(),
                     (*next_hot_pixel_iter)->get_center());
        seg.determine_direction(seg_dir);
        find_intersected_hot_pixels(seg,hot_pixels_intersected_set,
            number_of_intersections,pixel_size,mul_kd_tree);
        reroute_isr(hot_pixels_intersected_set,seg_output,
            number_of_intersections,false,pixel_size,
            int_output,mul_kd_tree);
      }
    } else {
      // insert second hot pixel
      hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
      ++hot_pixel_iter;
      seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
    }
  }

/*template<class Rep_>
void Snap_rounding_2<Rep_>::list_copy(
        std::list<Point_2>& target,std::list<Point_2>& source)
  {
    target = std::list<Point_2>();

    for(typename std::list<Point_2>::iterator i = source.begin();
        i != source.end();++i)
      target.push_back(*i);
      }*/

template<class Rep_>
void Snap_rounding_2<Rep_>::iterate(
        Polylines_container& output_container,
        NT pixel_size,
        bool int_output,
        bool do_isr,
        std::list<Segment_data<Rep> >& seg_list,
        Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree)
  {
    std::list<Point_2> seg_output;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter;
    int number_of_intersections;
    Hot_Pixel<Rep_> *hp;
    Direction seg_dir;

    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      seg_output.clear();
      iter->determine_direction(seg_dir);
      find_intersected_hot_pixels(*iter,hot_pixels_intersected_set,
        number_of_intersections,pixel_size,mul_kd_tree);
      // hot_pixels_intersected_set must have at least two hot pixels when the
      // segment is not in entirely inside a hot pixel enter first hot pixel
      hot_pixel_iter = hot_pixels_intersected_set.begin();
      if(hot_pixel_iter == hot_pixels_intersected_set.end()) {
        // segment entirely inside a pixel
        hp = new Hot_Pixel<Rep_>(iter->source(),pixel_size);
        seg_output.push_back(hp->get_center(int_output));
        delete(hp);
      } else {
        seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
        if(number_of_intersections > 1) {
          // segments that have at most one intersecting hot pixel are
          // done(it was inserted)
          if(do_isr)
            reroute_isr(hot_pixels_intersected_set,seg_output,
                        number_of_intersections,true,pixel_size,
                        int_output,mul_kd_tree);
          else
            reroute_sr(hot_pixels_intersected_set,seg_output,int_output);
	}
      }

      output_container.push_back(seg_output);
    }
  }

template<class Rep_>
void snap_rounding_2(
  typename std::list<typename Rep_::Segment_2>::const_iterator begin,
  typename std::list<typename Rep_::Segment_2>::const_iterator end,
  std::list<std::list<typename Rep_::Point_2> >& output_container,
  typename Rep_::FT pixel_size,
  bool do_isr,
  bool int_output,
  unsigned int number_of_kd_trees)
  {
    std::list<Segment_data<Rep_> > seg_list;
    Multiple_kd_tree<Rep_,Hot_Pixel<Rep_> *> *mul_kd_tree;

    // copy segments list
    while(begin != end) {
      seg_list.push_back(Segment_data<Rep_>(begin->source(),
                                            begin->target()));
      ++begin;
    }

    Snap_rounding_2<Rep_> s;

    s.find_hot_pixels_and_create_kd_trees(pixel_size,number_of_kd_trees,
         seg_list,&mul_kd_tree);
    s.iterate(output_container,pixel_size,int_output,do_isr,seg_list,
         mul_kd_tree);
  }

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_H
