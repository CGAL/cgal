// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// $Revision$
// $Name$
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>

#ifndef CGAL_SR_2_H
#define CGAL_SR_2_H

//#define DEBUG

#include <CGAL/leda_rational.h> 

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
// @@@@ special includes for pm
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/squared_distance_2.h>
// @@@@ end of special includes
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Iso_rectangle_2.h>

//#endif

//#include <CGAL/Sweep_line_tight_2.h>
#include <CGAL/Sweep_line_2.h>

#include <CGAL/Arr_segment_traits_2.h>
//#include <CGAL/Arr_polyline_traits.h>
#include <list>
#include <set>
#include <CGAL/leda_real.h>
#include "Snap_rounding_kd_2.h"

#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE

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
  Segment_data(Point_2 p_inp,Point_2 q_inp);
  Segment_2 segment() const {return(Segment_2(p,q));}
  Point_2 source() const {return(p);}
  Point_2 target() const {return(q);}
  inline void set_data(Point_2 inp_p,Point_2 inp_q);
  void determine_direction();
  bool equal(Segment_2 s);
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
  bool is_of_isrs;

public:
  void draw_hp(CGAL::Window_stream &w);
  template<class Out>
  void draw(Out &o) const;
  Hot_Pixel(Point_2 inp_point,NT inp_pixel_size,bool inp_isrs);
  ~Hot_Pixel();
  inline Point_2 get_center() const;
  inline Point_2 get_center(bool int_output) const;
  bool intersect_left(Segment_2 &seg) const;
  bool intersect_right(Segment_2 &seg) const;
  bool intersect_bot(Segment_2 &seg) const;
  bool intersect_top(Segment_2 &seg) const;
  bool intersect(Segment_data<Rep> &seg) const;
  bool belong_to_isrs() {return(is_of_isrs);}

  // @@@@ next function for isrs
  Point_2 get_ll() {return(Point_2(p_left.x(),p_down.y()));}
  Point_2 get_lr() {return(Point_2(p_right.x(),p_down.y()));}
  Point_2 get_ul() {return(Point_2(p_left.x(),p_up.y()));}
  Point_2 get_ur() {return(Point_2(p_right.x(),p_up.y()));}

  // @@@@ end of function for isrs

};

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

typedef CGAL::Arr_segment_traits_2<Rep_ >     Traits;
typedef Rep_                                         Rep;
typedef typename Rep::FT                             NT;
// @@@@ special typedefs for pm
//typedef CGAL::Quotient<CGAL::MP_Float>             NT2;
//typedef CGAL::Cartesian<NT>                        Kernel2;
//typedef CGAL::Arr_segment_traits_2<Kernel2> Traits2;
typedef CGAL::Pm_default_dcel<Traits>                      Dcel2;
typedef CGAL::Planar_map_2<Dcel2,Traits>                   Planar_map_2;
typedef CGAL::Planar_map_with_intersections_2<Planar_map_2> Pmwx;
typedef typename Traits::X_curve                                    X_curve2;
//typedef CGAL::squared_distance            squared_distance;
// @@@@ end special typedefs for pm
typedef typename Traits::X_curve                     X_curve;
typedef typename Traits::Curve                       Curve;
typedef std::list<X_curve>                           CurveContainer;
typedef typename CurveContainer::iterator            CurveContainerIter;
typedef CGAL::Iso_rectangle_2<Rep>                     Iso_rectangle_2;

public:
  friend class Segment_data<Rep>;
  friend class Hot_Pixel<Rep>;
  friend class hot_pixel_dir_cmp<Rep>;

  typedef typename Rep::Segment_2 Segment_2;
  typedef typename Rep::Point_2   Point_2;
  typedef std::list<Point_2> PointList;
  typedef typename std::list<std::list<Point_2> >
             Polylines_container;
  typedef typename Polylines_container::const_iterator Polyline_const_iterator;
  typedef typename Polylines_container::iterator Polyline_iterator;
  typedef typename std::list<typename Rep::Point_2> Points_container;
  typedef typename Points_container::const_iterator Point_const_iterator;
  typedef typename std::list<Segment_2> Segments_container;
  typedef typename Segments_container::const_iterator Segment_const_iterator;
  typedef typename Segments_container::iterator Segment_iterator;

  static bool erase_hp;
  static inline bool get_erase_hp() {return(erase_hp);}

  //! A constructor
  Snap_rounding_2(Segment_const_iterator begin,
                  Segment_const_iterator end,
                  NT inp_pixel_size,
                  bool inp_do_isr = true,
                  bool inp_int_output = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees,
                  // @@@@ next one is for isrs
                  NT inp_delta = 2);

  //! A constructor
  Snap_rounding_2(NT inp_pixel_size,
                  bool inp_do_isr = true,
                  bool inp_int_output = true,
                  int inp_number_of_kd_trees = default_number_of_kd_trees,
                  // @@@@ next one is for isrs
                  NT inp_delta = 2);

  //! A copy constructor
  Snap_rounding_2(const Snap_rounding_2& other);

  //! An operator =
  Snap_rounding_2& operator =(const Snap_rounding_2& other);

   //! Returns a constant iterator to the first input segment.
  inline Segment_const_iterator segments_begin() const {
             return(seg_2_list.begin());}

  //! Returns a constant iterator to the after-the-last of the input segments.
  inline Segment_const_iterator segments_end() const {
    return(seg_2_list.end());}

  //! Returns an iterator to the first input segment.
  inline Segment_iterator segments_begin() {return(seg_2_list.begin());}

  //! Returns an iterator to the after-the-last of the input segments.
  inline Segment_iterator segments_end() {return(seg_2_list.end());}

  //! Returns a constant iterator to the output of the first input segment.
  const Polyline_const_iterator polylines_begin();

  /*! Returns a constant iterator to the after-the-last of the output of
   *  the input segments.
   */
  const Polyline_const_iterator polylines_end();

  //! insert a segment
  bool insert(Segment_2 seg);

  //! The STL standard member function for insertion.
  bool push_back(Segment_2 seg);

  //! Insertion of an iterator range.
  template < class InputIterator >
    int insert(InputIterator first, InputIterator last);
  
  //! Remove a segment
  bool remove(Segment_2 seg);

  //! Remove all segments  .           .
  void clear();

  //! Change the number of kd-trees used for segment-hot pixel queries.
  bool change_number_of_kd_trees(int inp_number_of_kd_trees);

  //! Change the pixel size
  bool change_pixel_size(NT inp_pixel_size);

  // @@@@ isrs
  bool change_delta(NT inp_delta);

  /*! Determine whether to apply Iterated Snap Rounding (ISR)
      or Snap Rounding (SR).
   */
  void do_isr(bool inp_do_isr);

  // @@@@ next function
  void do_isrs(bool inp_do_isrs);

  /*! Determine the output mode: integer grid or closest grid point
      or Snap Rounding (SR).
   */
  void use_integer_mode(bool inp_int_output);

  template<class Out>
  void output(Out &o);
  void draw_hp(CGAL::Window_stream &w);

private:
  Rep _gt;
  enum Direction {UP_RIGHT,UP_LEFT,DOWN_RIGHT,DOWN_LEFT,UP,DOWN,LEFT,
                  RIGHT,POINT_SEG};

  std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> > hot_pixels_list;
  static Direction seg_dir;
  // the next variable is for lazy evaluation:
  // it determines whether an isr/sr work has
  // to be done (at the beginning, after insertion, etc.)
  bool need_sr;
  static const int default_number_of_kd_trees = 1;

  std::set<Hot_Pixel<Rep> *,hot_pixel_auclidian_cmp<Rep> > hp_set;
  NT pixel_size,min_x,max_x,min_y,max_y;
  Segments_container seg_2_list;
  std::list<Segment_data<Rep> > seg_list;
  std::list<std::list<Point_2> > segments_output_list;
  int number_of_segments,number_of_kd_trees;
  Multiple_kd_tree<Rep,Hot_Pixel<Rep> *> *mul_kd_tree;
  bool wheteher_to_do_isr;

  // @@@@ next one
  bool wheteher_to_do_isrs;

  bool int_output;
  // @@@@ next two for ISRS
  Pmwx pm;
  NT delta;
  // @@@@ end for isrs

  int NUMBER_OF_HP;

  void find_hot_pixels_and_create_kd_trees();

  //@@@@ next functions for isrs
  void produce_extra_hot_pixels(std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list,std::list<Segment_2>& segment_list);
  Point_2 next_center_to_right(const Point_2& p,const Segment_2& s);
  Point_2 next_center_to_left(const Point_2& p,const Segment_2& s);
  bool to_the_right(const Point_2& query_point,const Segment_2& s);
  bool to_the_left(const Point_2& query_point,const Segment_2& s);
  Segment_2 find_segment_to_right(const Point_2& p_center,const Point_2& query_point,std::list<Segment_2>& segment_list,NT& dis,bool& found);
  Segment_2 find_segment_to_left(const Point_2& p_center,const Point_2& query_point,std::list<Segment_2>& segment_list,NT& dis,bool& found);
  bool inside_bounding_box(const Point_2& query_point,const Segment_2& s);
  Point_2 find_point_on_segment_right(const Segment_2& s,const Point_2& p,bool& found);
  Point_2 find_point_on_segment_left(const Segment_2& s,const Point_2& p,bool& found);
  Point_2 find_point_on_segment_up(const Segment_2& s,const Point_2& p,bool& found);
  Point_2 find_point_on_segment_down(const Segment_2& s,const Point_2& p,bool& found);
  bool triangle_empty_up_right(const Segment_2& s,const Point_2& p,
        std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  bool triangle_empty_down_right(const Segment_2& s,const Point_2& p,
        std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  bool triangle_empty_up_left(const Segment_2& s,const Point_2& p,
        std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  bool triangle_empty_down_left(const Segment_2& s,const Point_2& p,
        std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  bool negative_slope(const Segment_2& s);
  bool positive_slope(const Segment_2& s);
  bool not_hot_pixel(Hot_Pixel<Rep_> *hp,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  void heat_pixel_up_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  void heat_pixel_down_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  void heat_pixel_up_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  void heat_pixel_down_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list);
  //@@@@ end functions for isrs

  void find_intersected_hot_pixels(Segment_data<Rep> &seg,
                         std::set<Hot_Pixel<Rep> *,
                         hot_pixel_dir_cmp<Rep> > &hot_pixels_intersected_set,
                         int &number_of_intersections);

  void reroute_sr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                  &inp_hot_pixels_intersected_set,std::list<Point_2>
                  &seg_output);
  void reroute_isr(std::set<Hot_Pixel<Rep> *,hot_pixel_dir_cmp<Rep> >
                   &inp_hot_pixels_intersected_set,std::list<Point_2>
                   &seg_output,int number_of_intersections,bool first_time);
  void iterate();
  void copy(const Snap_rounding_2<Rep_>& other);
  static inline Direction get_direction() {return(seg_dir);}
  static inline void set_direction(Direction dir) {seg_dir = dir;}
};

// ctor
template<class Rep_>
Segment_data<Rep_>::Segment_data() {}
template<class Rep_>
Segment_data<Rep_>::Segment_data(Point_2 p_inp,Point_2 q_inp) :
                    p(p_inp), q(q_inp) {}

// cctor
template<class Rep_>
Segment_data<Rep_>::Segment_data(const Segment_data& other)
{
  p = other.p;
  q = other.q;
}

template<class Rep_>
inline void Segment_data<Rep_>::set_data(Point_2 inp_p,Point_2 inp_q)
{
  p = inp_p;
  q = inp_q;
}

template<class Rep_>
bool Segment_data<Rep_>::equal(Segment_2 s)
{
  return(s.source() == p && s.target() == q);
}

template<class Rep_>
void Segment_data<Rep_>::determine_direction()
{
  Comparison_result cx = _gt.compare_x_2_object()(p,q);
  Comparison_result cy = _gt.compare_y_2_object()(p,q);

  if(cx == SMALLER) {
   if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_RIGHT);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::RIGHT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_RIGHT);
  } else if(cx == EQUAL) {
    if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::POINT_SEG);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN);
  } else {
    if(cy == SMALLER)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::UP_LEFT);
    else if(cy == EQUAL)
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::LEFT);
    else
      Snap_rounding_2<Rep_>::set_direction(Snap_rounding_2<Rep_>::DOWN_LEFT);
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
Hot_Pixel<Rep_>::Hot_Pixel(Point_2 inp_point,NT inp_pixel_size,bool is_isrs) :
                           pixel_size(inp_pixel_size),is_of_isrs(is_isrs)
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
    if(Snap_rounding_2<Rep_>::get_erase_hp()) {
      delete(right_seg);
      delete(left_seg);
      delete(top_seg);
      delete(bot_seg);
    }
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
bool Hot_Pixel<Rep_>::intersect_left(Segment_2 &seg) const
  {
    Object result;    
    Point_2 p;
    Segment_2 s;

    result = intersection(seg,*left_seg);

    if(assign(p,result)) {
      Comparison_result c_p = _gt.compare_y_2_object()(p,p_up);
      Comparison_result c_target = _gt.compare_y_2_object()(seg.target(),p_up);
      Comparison_result c_source = _gt.compare_y_2_object()(seg.source(),p_up);

      return(c_p != EQUAL || Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::UP_LEFT && c_source != EQUAL ||
             Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::DOWN_RIGHT && c_target != EQUAL);
    } else if(assign(s,result))
      return(true);
    else
      return(false);
  }


template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_right(Segment_2 &seg) const
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
        return(Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_RIGHT && c3 != EQUAL ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_LEFT && c4 != EQUAL);
      else if(c2 == EQUAL)
        return(false);// was checked
      else {
        Comparison_result c_target =
             _gt.compare_x_2_object()(p_right,seg.target());
        Comparison_result c_source =
             _gt.compare_x_2_object()(p_right,seg.source());

        return((Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::LEFT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::DOWN_LEFT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::UP_LEFT) &&
                c_target != EQUAL ||
                (Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::RIGHT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::DOWN_RIGHT ||
                Snap_rounding_2<Rep_>::get_direction() ==
                Snap_rounding_2<Rep_>::UP_RIGHT) &&
                c_source != EQUAL);
      }
    } else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_bot(Segment_2 &seg) const
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

      return(c_p != EQUAL || Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::UP_LEFT && c_target != EQUAL ||
             Snap_rounding_2<Rep_>::get_direction() ==
             Snap_rounding_2<Rep_>::DOWN_RIGHT && c_source != EQUAL);
    } else if(assign(s,result))
      return(true);
    else
      return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect_top(Segment_2 &seg) const
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
        return((Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_LEFT ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::DOWN_RIGHT) && c3 != EQUAL ||
               (Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_LEFT ||
               Snap_rounding_2<Rep_>::get_direction() ==
               Snap_rounding_2<Rep_>::UP_RIGHT) && c4 != EQUAL);
    } else
    return(false);
  }

template<class Rep_>
bool Hot_Pixel<Rep_>::intersect(Segment_data<Rep_> &seg) const
  {
    Segment_2 s = seg.segment();

    return(intersect_bot(s) || intersect_left(s) || intersect_right(s) ||
           intersect_top(s));
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

  return(
     // Point segment intersects only one pixel, thus ignored
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::UP_RIGHT &&
    (cx == SMALLER || 
    cx == EQUAL && cy == SMALLER) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::UP_LEFT &&
    (cx == LARGER || 
    cx == EQUAL && cy == SMALLER) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::DOWN_RIGHT &&
    (cx == SMALLER || 
    cx == EQUAL && cy == LARGER) ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::DOWN_LEFT &&
    (cx == LARGER || 
    cx == EQUAL && cy == LARGER) ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::UP &&
    cy == SMALLER ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::DOWN &&
    cy == LARGER ||
    Snap_rounding_2<Rep_>::get_direction() == Snap_rounding_2<Rep_>::LEFT &&
    cx == LARGER ||
    Snap_rounding_2<Rep_>::get_direction() ==
    Snap_rounding_2<Rep_>::RIGHT &&
    cx == SMALLER);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::to_the_right(const Point_2& query_point,const Segment_2& s)
{
  Point_2 so = s.source(),ta = s.target(),tmp;

  if(so.y() > ta.y()) {
    tmp = so;
    so = ta;
    ta = tmp;
  }

  #ifdef ISR_DEBUG
    std::cout << "so = " << so << std::endl;
    std::cout << "query_point = " << query_point << std::endl;
    std::cout << "ta = " << ta << std::endl;
  #endif

  return(right_turn(so,query_point,ta));
}

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::draw_hp(CGAL::Window_stream &w)
{
  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::iterator i = hot_pixels_list.begin();
      i != hot_pixels_list.end();++i) {
    Point_2 p(i->first);
    if(i->second->belong_to_isrs())
      w << CGAL::BLUE;
    else
      w << CGAL::GREEN;
 
    w << Iso_rectangle_2(Point_2(p.x() - pixel_size / 2.0,
		                 p.y() - pixel_size / 2.0),
                         Point_2(p.x() + pixel_size / 2.0,
				 p.y() + pixel_size / 2.0));
  }
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::to_the_left(const Point_2& query_point,const Segment_2& s)
{
  Point_2 so = s.source(),ta = s.target(),tmp;

  if(so.y() > ta.y()) {
    tmp = so;
    so = ta;
    ta = tmp;
  }

  #ifdef ISR_DEBUG
    std::cout << "so = " << so << std::endl;
    std::cout << "query_point = " << query_point << std::endl;
    std::cout << "ta = " << ta << std::endl;
  #endif

  return(left_turn(so,query_point,ta));
}


// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Segment_2 Snap_rounding_2<Rep_>::find_segment_to_right(
     const Point_2& p_center,const Point_2& query_point,
     std::list<Segment_2>& segment_list,NT& dis,bool& found)
{
  typename std::list<Segment_2>::const_iterator iter;
  NT curr_dis;
  bool first_time = true;
  Segment_2 result_segment;

  #ifdef ISR_DEBUG
    std::cout << "query_point = " << query_point << std::endl;
  #endif
  found = false;
  for(iter = segment_list.begin();iter != segment_list.end();++iter) {
    Segment_2 s = *iter;
    #ifdef ISR_DEBUG
      std::cout << "s = " << s << std::endl;
    #endif
    if(to_the_right(query_point,s)) { // both to the right and a x-ray shooting intersects it
      curr_dis = CGAL::squared_distance(p_center,s);
      if(first_time || curr_dis < dis) {
        first_time = false;
        result_segment = s;
        dis = curr_dis;
        found = true;
      }
    }
  }

  return(result_segment);
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Segment_2 Snap_rounding_2<Rep_>::find_segment_to_left(
     const Point_2& p_center,const Point_2& query_point,
     std::list<Segment_2>& segment_list,NT& dis,bool& found)
{
  typename std::list<Segment_2>::const_iterator iter;
  NT curr_dis;
  bool first_time = true;
  Segment_2 result_segment;

  #ifdef ISR_DEBUG
    std::cout << "query_point = " << query_point << std::endl;
  #endif
  found = false;
  for(iter = segment_list.begin();iter != segment_list.end();++iter) {
    Segment_2 s = *iter;
    #ifdef ISR_DEBUG
      std::cout << "s = " << s << std::endl;
    #endif
    if(to_the_left(query_point,s)) { // both to the LEFT and a x-ray shooting intersects it
      curr_dis = CGAL::squared_distance(p_center,s);
      if(first_time || curr_dis < dis) {
        first_time = false;
        result_segment = s;
        dis = curr_dis;
        found = true;
      }
    }
  }

  return(result_segment);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::inside_bounding_box(const Point_2& query_point,
                                                const Segment_2& s)
{
  NT x = query_point.x(),
     y = query_point.y(),
     s1x = s.source().x(),
     s1y = s.source().y(),
     s2x = s.target().x(),
     s2y = s.target().y();

  return(!(x < s1x - pixel_size && x < s2x - pixel_size || x > s1x + pixel_size && x > s2x + pixel_size ||
         y < s1y && y < s2y || y > s1y && y > s2y));
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::find_point_on_segment_right(const Segment_2& s,const Point_2& p,bool& found)
{
  Point_2 inter_p;
  Object result;

  Segment_2 crossing_s(p,Point_2(s.max().x(),p.y()));

  #ifdef ISR_DEBUG
    std::cout << "crossing_s = " << crossing_s << std::endl;
    std::cout << "s = " << s << std::endl;
  #endif

  result = intersection(s,Segment_2(p,Point_2(s.max().x(),p.y())));
  if(assign(inter_p,result))
    found = true;
  else
    found = false;

  return(inter_p);
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::find_point_on_segment_left(const Segment_2& s,const Point_2& p,bool& found)
{
  Point_2 inter_p;
  Object result;

  result = intersection(s,Segment_2(p,Point_2(s.min().x(),p.y())));
  if(assign(inter_p,result))
    found = true;
  else
    found = false;

  return(inter_p);
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::find_point_on_segment_up(const Segment_2& s,const Point_2& p,bool& found)
{
  Point_2 inter_p;
  Object result;

  result = intersection(s,Segment_2(p,Point_2(p.x(),s.source().y() > s.target().y() ?
                                    s.source().y() : s.target().y())));
  if(assign(inter_p,result))
    found = true;
  else
    found = false;

  return(inter_p);
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::find_point_on_segment_down(const Segment_2& s,const Point_2& p,bool& found)
{
  Point_2 inter_p;
  Object result;

  result = intersection(s,Segment_2(p,Point_2(p.x(),s.source().y() > s.target().y() ?
                                    s.target().y() : s.source().y())));
  if(assign(inter_p,result))
    found = true;
  else
    found = false;

  return(inter_p);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::negative_slope(const Segment_2& s)
{
  return(
	 s.source().x() < s.target().x() ?
         s.source().y() > s.target().y() :
         s.source().y() < s.target().y());
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::positive_slope(const Segment_2& s)
{
  return(
	 s.source().x() < s.target().x() ?
         s.source().y() < s.target().y() :
         s.source().y() > s.target().y());
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::triangle_empty_up_right(const Segment_2& s,const Point_2& p_inp,
      std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Point_2 p(p_inp.x() + pixel_size / 2,p_inp.y() + pixel_size / 2);
  bool found;
  Point_2 s1 = find_point_on_segment_right(s,p,found);
  if(!found)
    return(true);
  Point_2 s2 = find_point_on_segment_up(s,p,found);
  if(!found)
    return(true);
  Point_2 h1,h2,h3,h4;

  #ifdef ISR_DEBUG
    std::cout << "s2 = " << s << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;
  #endif

  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator i = hot_pixels_list.begin();i != hot_pixels_list.end();++i) {
    h1 = i->second->get_ll();
    h2 = i->second->get_lr();
    h3 = i->second->get_ul();
    h4 = i->second->get_ur();
    #ifdef ISR_DEBUG
      std::cout << "h1 = " << h1 << std::endl;
      std::cout << "h2 = " << h2 << std::endl;
      std::cout << "h3 = " << h3 << std::endl;
      std::cout << "h4 = " << h4 << std::endl;
    #endif
    if(left_turn(s1,h1,p) && left_turn(s2,h1,s1) && left_turn(p,h1,s2) ||
       left_turn(s1,h2,p) && left_turn(s2,h2,s1) && left_turn(p,h2,s2) ||
       left_turn(s1,h3,p) && left_turn(s2,h3,s1) && left_turn(p,h3,s2) ||
       left_turn(s1,h4,p) && left_turn(s2,h4,s1) && left_turn(p,h4,s2)) {
      #ifdef ISR_DEBUG
        std::cout << "inside TRIANGLE" << std::endl;
      #endif
      return(false);
    }
  }

  return(true);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::triangle_empty_down_right(const Segment_2& s,const Point_2& p_inp,
      std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Point_2 p(p_inp.x() + pixel_size / 2,p_inp.y() - pixel_size / 2);
  bool found;
  Point_2 s1 = find_point_on_segment_right(s,p,found);
  if(!found)
    return(true);
  Point_2 s2 = find_point_on_segment_down(s,p,found);
  if(!found)
    return(true);
  Point_2 h1,h2,h3,h4;

  #ifdef ISR_DEBUG
    std::cout << "s2 = " << s << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;
  #endif

  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator i = hot_pixels_list.begin();i != hot_pixels_list.end();++i) {
    h1 = i->second->get_ll();
    h2 = i->second->get_lr();
    h3 = i->second->get_ul();
    h4 = i->second->get_ur();
    #ifdef ISR_DEBUG
      std::cout << "h1 = " << h1 << std::endl;
      std::cout << "h2 = " << h2 << std::endl;
      std::cout << "h3 = " << h3 << std::endl;
      std::cout << "h4 = " << h4 << std::endl;
    #endif
    if(right_turn(s1,h1,p) && right_turn(s2,h1,s1) && right_turn(p,h1,s2) ||
       right_turn(s1,h2,p) && right_turn(s2,h2,s1) && right_turn(p,h2,s2) ||
       right_turn(s1,h3,p) && right_turn(s2,h3,s1) && right_turn(p,h3,s2) ||
       right_turn(s1,h4,p) && right_turn(s2,h4,s1) && right_turn(p,h4,s2)) {
      #ifdef ISR_DEBUG
        std::cout << "inside TRIANGLE" << std::endl;
      #endif
      return(false);
    }
  }

  return(true);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::triangle_empty_up_left(const Segment_2& s,const Point_2& p_inp,
      std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Point_2 p(p_inp.x() - pixel_size / 2,p_inp.y() + pixel_size / 2);
  bool found;
  Point_2 s1 = find_point_on_segment_left(s,p,found);
  if(!found)
    return(true);
  Point_2 s2 = find_point_on_segment_up(s,p,found);
  if(!found)
    return(false);
  Point_2 h1,h2,h3,h4;

  #ifdef ISR_DEBUG
    std::cout << "s2 = " << s << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;
  #endif

  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator i = hot_pixels_list.begin();i != hot_pixels_list.end();++i) {
    h1 = i->second->get_ll();
    h2 = i->second->get_lr();
    h3 = i->second->get_ul();
    h4 = i->second->get_ur();
    #ifdef ISR_DEBUG
      std::cout << "h1 = " << h1 << std::endl;
      std::cout << "h2 = " << h2 << std::endl;
      std::cout << "h3 = " << h3 << std::endl;
      std::cout << "h4 = " << h4 << std::endl;
    #endif
    if(right_turn(s1,h1,p) && right_turn(s2,h1,s1) && right_turn(p,h1,s2) ||
       right_turn(s1,h2,p) && right_turn(s2,h2,s1) && right_turn(p,h2,s2) ||
       right_turn(s1,h3,p) && right_turn(s2,h3,s1) && right_turn(p,h3,s2) ||
       right_turn(s1,h4,p) && right_turn(s2,h4,s1) && right_turn(p,h4,s2)) {
      #ifdef ISR_DEBUG
        std::cout << "inside TRIANGLE" << std::endl;
      #endif
      return(false);
    }
  }

  return(true);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::triangle_empty_down_left(const Segment_2& s,const Point_2& p_inp,
      std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Point_2 p(p_inp.x() - pixel_size / 2,p_inp.y() - pixel_size / 2);
  bool found;
  Point_2 s1 = find_point_on_segment_left(s,p,found);
  if(!found)
    return(true);
  Point_2 s2 = find_point_on_segment_down(s,p,found);
  if(!found)
    return(false);
  Point_2 h1,h2,h3,h4;

  #ifdef ISR_DEBUG
    std::cout << "s2 = " << s << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;
  #endif

  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator i = hot_pixels_list.begin();i != hot_pixels_list.end();++i) {
    h1 = i->second->get_ll();
    h2 = i->second->get_lr();
    h3 = i->second->get_ul();
    h4 = i->second->get_ur();
    #ifdef ISR_DEBUG
      std::cout << "h1 = " << h1 << std::endl;
      std::cout << "h2 = " << h2 << std::endl;
      std::cout << "h3 = " << h3 << std::endl;
      std::cout << "h4 = " << h4 << std::endl;
    #endif
    if(left_turn(s1,h1,p) && left_turn(s2,h1,s1) && left_turn(p,h1,s2) ||
       left_turn(s1,h2,p) && left_turn(s2,h2,s1) && left_turn(p,h2,s2) ||
       left_turn(s1,h3,p) && left_turn(s2,h3,s1) && left_turn(p,h3,s2) ||
       left_turn(s1,h4,p) && left_turn(s2,h4,s1) && left_turn(p,h4,s2)) {
      #ifdef ISR_DEBUG
        std::cout << "inside TRIANGLE" << std::endl;
      #endif
      return(false);
    }
  }

  return(true);
}

// @@@@ a function for ISRS
template<class Rep_>
bool Snap_rounding_2<Rep_>::not_hot_pixel(Hot_Pixel<Rep_> *hp,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  for(typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator iter = hot_pixels_list.begin();
      iter != hot_pixels_list.end();++iter) {
    if(hp->get_center() == iter->first)
      return(false);
  }

  return(true);
}

/*// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_up_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  #ifdef ISR_DEBUG
    std::cout << "** p = " << p << std::endl;
    std::cout << "** s = " << s << std::endl;
  #endif

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.max().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
    #endif
    // reach a point which is in the left pixel columns to the intersection
    Segment_2 r(Point_2(inter_p.x() - pixel_size,inter_p.y()),
                Point_2(inter_p.x() - pixel_size,s.target().y() > s.source().y() ?
                                                 s.target().y() : s.source().y()));
    result = intersection(r,s);
    #ifdef ISR_DEBUG
      std::cout << "r = " << r << std::endl;
      std::cout << "s = " << s << std::endl;
    #endif

    if(assign(inter_p,result)) {
      #ifdef ISR_DEBUG
        std::cout << "2) inter_p = " << inter_p << std::endl;
      #endif
      // heat the resective pixel and insert it to hot_pixels_list
      Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
      if(not_hot_pixel(hp,hot_pixels_list)) {
        hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
 	  			hp->get_center(),hp));
        ++NUMBER_OF_HP;
      } else
        delete hp;
    } else {
      std::cout << "error 312\n";
      exit(-1);
    }
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}*/

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_up_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  #ifdef ISR_DEBUG
    std::cout << "** p = " << p << std::endl;
    std::cout << "** s = " << s << std::endl;
  #endif

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.max().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
    #endif

    // heat the resective pixel and insert it to hot_pixels_list
    Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
    if(not_hot_pixel(hp,hot_pixels_list)) {
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
  			hp->get_center(),hp));
        ++NUMBER_OF_HP;
    } else
      delete hp;
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}

/*// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_down_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.max().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
    #endif
    // reach a point which is in the left pixel columns to the intersection
    Segment_2 r(Point_2(inter_p.x() - pixel_size,inter_p.y()),
                Point_2(inter_p.x() - pixel_size,s.target().y() < s.source().y() ?
                                                 s.target().y() : s.source().y()));
    result = intersection(r,s);
    #ifdef ISR_DEBUG
      std::cout << "r = " << r << std::endl;
      std::cout << "s = " << s << std::endl;
    #endif

    if(assign(inter_p,result)) {
      #ifdef ISR_DEBUG
        std::cout << "2) inter_p = " << inter_p << std::endl;
      #endif
      // heat the resective pixel and insert it to hot_pixels_list
      Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
      if(not_hot_pixel(hp,hot_pixels_list)) {
        hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
 	  			hp->get_center(),hp));
        ++NUMBER_OF_HP;
      } else
        delete hp;
    } else {
      std::cout << "error 311\n";
      exit(-1);
    }
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}*/

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_down_right(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.max().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
    #endif

    // heat the resective pixel and insert it to hot_pixels_list
    Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
    if(not_hot_pixel(hp,hot_pixels_list)) {
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
  			hp->get_center(),hp));
      ++NUMBER_OF_HP;
    } else
      delete hp;
  } else {
    std::cout << "error 311\n";
    exit(-1);
  }
}

/* @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_up_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.min().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
      std::cout << "t = " << t << std::endl;
    #endif
    // reach a point which is in the right pixel columns to the intersection
    Segment_2 r(Point_2(inter_p.x() + pixel_size,inter_p.y()),
                Point_2(inter_p.x() + pixel_size,s.target().y() > s.source().y() ?
                                                 s.target().y() : s.source().y()));
    result = intersection(r,s);
    #ifdef ISR_DEBUG
      std::cout << "r = " << r << std::endl;
      std::cout << "s = " << s << std::endl;
    #endif
    if(assign(inter_p,result)) {
      #ifdef ISR_DEBUG
        std::cout << "2) inter_p = " << inter_p << std::endl;
      #endif
      // heat the resective pixel and insert it to hot_pixels_list
      Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
      if(not_hot_pixel(hp,hot_pixels_list)) {
        hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
 				hp->get_center(),hp));
        ++NUMBER_OF_HP;
      } else
        delete hp;

    } else {
      std::cout << "error 321\n";
      exit(-1);
    }
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}*/

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_up_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the right on s
  Segment_2 t(p,Point_2(s.min().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
      std::cout << "t = " << t << std::endl;
    #endif

    // heat the resective pixel and insert it to hot_pixels_list
    Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
    if(not_hot_pixel(hp,hot_pixels_list)) {
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
			hp->get_center(),hp));
      ++NUMBER_OF_HP;
    } else
      delete hp;
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}

/* @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_down_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the left on s
  Segment_2 t(p,Point_2(s.min().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
      std::cout << "t = " << t << std::endl;
    #endif
    // reach a point which is in the right pixel columns to the intersection
    Segment_2 r(Point_2(inter_p.x() + pixel_size,inter_p.y()),
                Point_2(inter_p.x() + pixel_size,s.target().y() < s.source().y() ?
                                                 s.target().y() : s.source().y()));
    result = intersection(r,s);
    #ifdef ISR_DEBUG
      std::cout << "r = " << r << std::endl;
      std::cout << "s = " << s << std::endl;
    #endif
    if(assign(inter_p,result)) {
      #ifdef ISR_DEBUG
        std::cout << "2) inter_p = " << inter_p << std::endl;
      #endif
      // heat the resective pixel and insert it to hot_pixels_list
      Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
      if(not_hot_pixel(hp,hot_pixels_list)) {
        hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
 				hp->get_center(),hp));
        ++NUMBER_OF_HP;
      } else
        delete hp;

    } else {
      std::cout << "error 322\n";
      exit(-1);
    }
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}*/

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::heat_pixel_down_left(const Point_2& p,const Segment_2& s,std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list)
{
  Object result;
  Point_2 inter_p;  

  // find the intersction to the left on s
  Segment_2 t(p,Point_2(s.min().x(),p.y()));
  result = intersection(t,s);
  if(assign(inter_p,result)) {
    #ifdef ISR_DEBUG
      std::cout << "1) inter_p = " << inter_p << std::endl;
      std::cout << "t = " << t << std::endl;
    #endif

    // heat the resective pixel and insert it to hot_pixels_list
    Hot_Pixel<Rep_> *hp = new Hot_Pixel<Rep_>(inter_p,pixel_size,true);
    if(not_hot_pixel(hp,hot_pixels_list)) {
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
			hp->get_center(),hp));
      ++NUMBER_OF_HP;
    } else
      delete hp;
  } else {
    std::cout << "error 4\n";
    exit(-1);
  }
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::next_center_to_right(const Point_2& p,const Segment_2& s)
{
  Point_2 ass_p;
  Object result;
  Segment_2 inter_s(p,Point_2(s.max().x(),p.y()));
  #ifdef ISR_DEBUG
    std::cout << "s = " << s << std::endl;
    std::cout << "inter_s = " << inter_s << std::endl;
  #endif
  result = intersection(s,inter_s);
  if(assign(ass_p,result))
    return(ass_p);
  else {
    std::cout << "error 51\n";
    exit(-1);
  }
}

// @@@@ a function for ISRS
template<class Rep_>
typename Snap_rounding_2<Rep_>::Point_2 Snap_rounding_2<Rep_>::next_center_to_left(const Point_2& p,const Segment_2& s)
{
  Point_2 ass_p;
  Object result;
  Segment_2 inter_s(p,Point_2(s.min().x(),p.y()));
  #ifdef ISR_DEBUG
    std::cout << "s = " << s << std::endl;
    std::cout << "inter_s = " << inter_s << std::endl;
  #endif
  result = intersection(s,inter_s);
  if(assign(ass_p,result))
    return(ass_p);
  else {
    std::cout << "error 52\n";
    exit(-1);
  }
}

// @@@@ a function for ISRS
template<class Rep_>
void Snap_rounding_2<Rep_>::produce_extra_hot_pixels(std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >& hot_pixels_list,std::list<Segment_2>& segment_list)
{
  typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator iter;
  NT R = 0.70710678 * pixel_size;
  NT dis1 = delta * delta;
  NT dis2 = delta * delta + R * R + 2 * delta * R;

  #ifdef ISR_DEBUG
    std::cout << "delta = " << delta << std::endl;
  #endif

  for(iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
    // UP RIGHT
    Point_2 p_center = iter->first;
    Point_2 query_point = Point_2(p_center.x(),p_center.y() + pixel_size);
    bool done = false;

    // UP RIGHT
    while(!done) {
      bool found;
      NT sq_dis;
      Segment_2 first_s = find_segment_to_right(p_center,query_point,segment_list,sq_dis,found);
      if(found && //negative_slope(first_s) &&
         inside_bounding_box(query_point,first_s)) {
        if(sq_dis < dis1)
          query_point = next_center_to_right(query_point,first_s);
        else if(sq_dis < dis2 &&
                triangle_empty_up_right(first_s,query_point,hot_pixels_list) && negative_slope(first_s)) {
          heat_pixel_up_right(query_point,first_s,hot_pixels_list);// insert to the end of the list
          done = true;
	} else
          done = true;
      } else
	 done = true;
    }

    // UP LEFT
    p_center = iter->first;
    query_point = Point_2(p_center.x(),p_center.y() + pixel_size);
    done = false;
    while(!done) {
      bool found;
      NT sq_dis;
      Segment_2 first_s = find_segment_to_left(p_center,query_point,segment_list,sq_dis,found);
      if(found && //positive_slope(first_s) &&
         inside_bounding_box(query_point,first_s)) {
        if(sq_dis < dis1) {
          query_point = next_center_to_left(query_point,first_s);
        } else if(sq_dis < dis2 &&
                triangle_empty_up_left(first_s,query_point,hot_pixels_list) && positive_slope(first_s)) {
          heat_pixel_up_left(query_point,first_s,hot_pixels_list);// insert to the end of the list
          done = true;
	} else
          done = true;
      } else
	 done = true;
    }

    // DOWN RIGHT
    p_center = iter->first;
    #ifdef ISR_DEBUG
      std::cout << "p_center = " << p_center << std::endl;
    #endif
    query_point = Point_2(p_center.x(),p_center.y() - pixel_size);
    done = false;
    while(!done) {
      #ifdef ISR_DEBUG
        std::cout << "query_point = " << query_point << std::endl;
      #endif
      bool found;
      NT sq_dis;
      Segment_2 first_s = find_segment_to_right(p_center,query_point,segment_list,sq_dis,found);
      if(found && //negative_slope(first_s) &&
         inside_bounding_box(query_point,first_s)) {
        #ifdef ISR_DEBUG
          std::cout << "first_s = " << first_s << std::endl;
	  std::cout << "sq_dis = " << sq_dis << std::endl;
        #endif
        if(sq_dis < dis1)
          query_point = next_center_to_right(query_point,first_s);
        else if(sq_dis < dis2 &&
                triangle_empty_down_right(first_s,query_point,hot_pixels_list) && positive_slope(first_s)) {
          heat_pixel_down_right(query_point,first_s,hot_pixels_list);// insert to the end of the list
          done = true;
	} else
          done = true;
      } else
	 done = true;
    }

    // DOWN LEFT
    p_center = iter->first;
    query_point = Point_2(p_center.x(),p_center.y() - pixel_size);
    done = false;
    while(!done) {
      bool found;
      NT sq_dis;
      Segment_2 first_s = find_segment_to_left(p_center,query_point,segment_list,sq_dis,found);
      #ifdef ISR_DEBUG
          std::cout << "query_point = " << query_point << std::endl;
          std::cout << "first_s = " << first_s << std::endl;
	  std::cout << "sq_dis = " << sq_dis << std::endl;
      #endif
      if(found && //positive_slope(first_s) &&
         inside_bounding_box(query_point,first_s)) {
        if(sq_dis < dis1) {
          query_point = next_center_to_left(query_point,first_s);
        } else if(sq_dis < dis2 &&
                triangle_empty_down_left(first_s,query_point,hot_pixels_list) && negative_slope(first_s)) {
          heat_pixel_down_left(query_point,first_s,hot_pixels_list);// insert to the end of the list
          done = true;
	} else
          done = true;
      } else
	 done = true;
    }
  }

  std::cout << "number of new hp is " << NUMBER_OF_HP << std::endl;
  std::cout << "total number of hp is " << hot_pixels_list.size() << std::endl;

  // below is the ray shooting code
  //typename std::list<std::pair<Point_2,Hot_Pixel<Rep_> *> >::const_iterator iter;

  //for(iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
  // Point_2 p_center = iter->first;
  //  Planar_map_2::Locate_type lt;
  //  Planar_map_2::Halfedge_handle e = pm.vertical_ray_shoot(p_center, lt, true);
  //}

  // Shoot a vertical ray upward from p:
  //Point_2 p(95, 30);
  //Planar_map_2::Locate_type lt;

  //std::cout << "Upward vertical ray shooting from " << p << std::endl; 
  //Planar_map_2::Halfedge_handle e = pm.vertical_ray_shoot(p, lt, true);
  //std::cout << "returned the curve " << e->curve() <<  ", oriented toward " 
  //<< e->target()->point() << std::endl;
}


template<class Rep_>
void Snap_rounding_2<Rep_>::find_hot_pixels_and_create_kd_trees()
  {
    Hot_Pixel<Rep_> *hp;
    typename std::list<Segment_data<Rep_> >::iterator iter1;
    Object result;
    Point_2 p;

    NUMBER_OF_HP = 0;
    std::list<X_curve> segments;
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
      hp = new Hot_Pixel<Rep_>(*v_iter,pixel_size,false);
      hot_pixels_list.push_back(std::pair<Point_2,Hot_Pixel<Rep_> *>(
				hp->get_center(),hp));
    }

    std::cout << "number of hot pixels : " << hot_pixels_list.size() << std::endl;

    // create kd multiple tree
    // create simple_list from seg_list
    std::list<Segment_2> simple_seg_list;
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter)
      simple_seg_list.push_back(Segment_2(iter->source(),iter->target()));


    // @@@@ for ISRS : create new hot pixels
    if(wheteher_to_do_isrs)
      produce_extra_hot_pixels(hot_pixels_list,simple_seg_list);


    mul_kd_tree = new Multiple_kd_tree<Rep,Hot_Pixel<Rep> *>(hot_pixels_list,
                  number_of_kd_trees,simple_seg_list);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::find_intersected_hot_pixels(Segment_data<Rep_>
                    &seg,
                    std::set<Hot_Pixel<Rep_> *,
                    hot_pixel_dir_cmp<Rep_> > &hot_pixels_intersected_set,
                    int &number_of_intersections)
  {
    typename std::list<Hot_Pixel<Rep_> *>::iterator iter;

    hot_pixels_intersected_set.clear();
    seg.determine_direction();
    number_of_intersections = 0;

    std::list<Hot_Pixel<Rep_> *> hot_pixels_list;
    mul_kd_tree->get_intersecting_points(hot_pixels_list,
	   Segment_2(seg.segment()),pixel_size);

    for(iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
      if((*iter)->intersect(seg))
        hot_pixels_intersected_set.insert(*iter);
    }

    number_of_intersections = hot_pixels_intersected_set.size();
  }


template<class Rep_>
void Snap_rounding_2<Rep_>::reroute_sr(std::set<Hot_Pixel<Rep_> *,
     hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
     std::list<Point_2> &seg_output)
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
void Snap_rounding_2<Rep_>::reroute_isr(std::set<Hot_Pixel<Rep_> *,
   hot_pixel_dir_cmp<Rep_> > &inp_hot_pixels_intersected_set,
   std::list<Point_2> &seg_output,int number_of_intersections,bool first_time)
  {
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter,next_hot_pixel_iter,before_last_hot_pixel_iter;
    Segment_data<Rep_> seg;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;

    if(number_of_intersections > 2 || first_time) {
      before_last_hot_pixel_iter = inp_hot_pixels_intersected_set.end();
      --before_last_hot_pixel_iter;

      for(hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
          hot_pixel_iter != before_last_hot_pixel_iter;++hot_pixel_iter) {
        next_hot_pixel_iter = hot_pixel_iter;
        ++next_hot_pixel_iter;
        seg.set_data((*hot_pixel_iter)->get_center(),
                     (*next_hot_pixel_iter)->get_center());
        seg.determine_direction();
        find_intersected_hot_pixels(seg,hot_pixels_intersected_set,
            number_of_intersections);
        reroute_isr(hot_pixels_intersected_set,seg_output,
            number_of_intersections,false);
      }
    } else {
      // insert second hot pixel
      hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
      ++hot_pixel_iter;
      seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
    }
  }


template<class Rep_>
void Snap_rounding_2<Rep_>::iterate()
  {
    int number_of_links = 0;

    std::list<Point_2> seg_output;
    std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >
      hot_pixels_intersected_set;
    typename std::set<Hot_Pixel<Rep_> *,hot_pixel_dir_cmp<Rep_> >::
      iterator hot_pixel_iter;
    int number_of_intersections;
    Hot_Pixel<Rep_> *hp;
    segments_output_list.clear();
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      seg_output.clear();
      iter->determine_direction();
      find_intersected_hot_pixels(*iter,hot_pixels_intersected_set,
        number_of_intersections);
      // hot_pixels_intersected_set must have at least two hot pixels when the
      // segment is not in entirely inside a hot pixel enter first hot pixel
      hot_pixel_iter = hot_pixels_intersected_set.begin();
      if(hot_pixel_iter == hot_pixels_intersected_set.end()) {
        // segment entirely inside a pixel
        hp = new Hot_Pixel<Rep_>(iter->source(),pixel_size,false);
        seg_output.push_back(hp->get_center(int_output));
        erase_hp = true;
        delete(hp);
        erase_hp = false;
      } else {
        seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
        if(number_of_intersections > 1) {
          // segments that have at most one intersecting hot pixel are
          // done(it was inserted)
          if(wheteher_to_do_isr || wheteher_to_do_isrs) // @@@@ the second part of the condition
            reroute_isr(hot_pixels_intersected_set,seg_output,
                        number_of_intersections,true);
          else
            reroute_sr(hot_pixels_intersected_set,seg_output);
	}
      }

      number_of_links += seg_output.size() - 1;
      segments_output_list.push_back(seg_output);
    }

    std::cout << "total number of links = " << number_of_links << std::endl;
  }

template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(
  Segment_const_iterator begin,
  Segment_const_iterator end,
  NT inp_pixel_size,
  bool inp_do_isr,
  bool inp_int_output,
  int inp_number_of_kd_trees,
  // @@@@ next one is for isrs
  NT inp_delta)
  {
    // @@@@ for isrs
    if(inp_delta < 2 * inp_pixel_size) {
      std::cout << "small max-deviation " << inp_delta << ". exit.\n";
      exit(-1);
    } else
      delta = inp_delta;
    // initialize approximation angles map    
    erase_hp = false;
    wheteher_to_do_isr = inp_do_isr;
    // @@@@ next
    wheteher_to_do_isrs = true;
    int_output = inp_int_output;
    pixel_size = inp_pixel_size;
    number_of_segments = 0;
    number_of_kd_trees = inp_number_of_kd_trees;
    need_sr = true;
    // copy segments list
    while(begin != end) {
      seg_list.push_back(Segment_data<Rep_>(begin->source(),
                                            begin->target()));
      seg_2_list.push_back(*begin);
      ++number_of_segments;

      // @@@@ insert to planar map
      //std:: cout << "1st : " << begin->source() <<
      //   " snd : " << begin->target() << std::endl;
      //X_curve2 c(begin->source(),begin->target());
      //std::cout << "bbb\n";
      //pm.insert(c);
      // @@@@ end of pmwx for here

      ++begin;
    }
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::copy(const Snap_rounding_2<Rep_>& other)
{
  // @@@@ for isrs
  NUMBER_OF_HP = other.NUMBER_OF_HP;
  delta = other.delta;

  erase_hp = false;
  wheteher_to_do_isr = other.wheteher_to_do_isr;
  // @@@@ next
  wheteher_to_do_isrs = other.wheteher_to_do_isrs;
  int_output = other.int_output;
  pixel_size = other.pixel_size;
  number_of_segments = other.number_of_segments;
  number_of_kd_trees = other.number_of_kd_trees;
  need_sr = true;
  seg_list = other.seg_list;
  seg_2_list = other.seg_2_list;
}

// cctor
template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(const Snap_rounding_2<Rep_>& other)
  {
    copy(other);
  }

// operator =
template<class Rep_>
Snap_rounding_2<Rep_>&
Snap_rounding_2<Rep_>::operator =(const Snap_rounding_2<Rep_>& other)
  {
    copy(other);

    return(*this);
  }

template<class Rep_>
Snap_rounding_2<Rep_>::Snap_rounding_2(
  NT inp_pixel_size,
  bool inp_do_isr,
  bool inp_int_output,
  int inp_number_of_kd_trees,
  // @@@@ for isrs
  NT inp_delta)
  {
    // @@@@ for isrs
    if(inp_delta < 2 * inp_pixel_size) {
      std::cout << "small max-deviation " << inp_delta << ". exit.\n";
      exit(-1);
    } else
      delta = inp_delta;

    // initialize approximation angles map
    need_sr = true;
    erase_hp = false;
    wheteher_to_do_isr = inp_do_isr;
    // @@@@@ next
    wheteher_to_do_isrs = true;

    int_output = inp_int_output;
    pixel_size = inp_pixel_size;
    number_of_segments = 0;
    number_of_kd_trees = inp_number_of_kd_trees;
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::insert(Segment_2 seg)
  {
    need_sr = true;
    seg_list.push_back(Segment_data<Rep_>(
					  seg.source(),
                         seg.target()));

    seg_2_list.push_back(seg);
    ++number_of_segments;

    return(true);
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::push_back(Segment_2 seg)
  {
    return(insert(seg));
  }

template < class Rep_ >
template < class InputIterator >
int
Snap_rounding_2<Rep_>::insert(InputIterator first, InputIterator last)
  {
    need_sr = true;
    int n = 0;
    while(first != last){
      if(insert(*first)){
	n++;
      }
      ++first;
    }
    return n;
  }


template<class Rep_>
const typename Snap_rounding_2<Rep_>::Polyline_const_iterator
      Snap_rounding_2<Rep_>::polylines_begin()
{
  if(need_sr) {
    hot_pixels_list.clear();
    need_sr = false;
    find_hot_pixels_and_create_kd_trees();
    iterate();
  }    

  return(segments_output_list.begin());
}

template<class Rep_>
const typename Snap_rounding_2<Rep_>::Polyline_const_iterator
      Snap_rounding_2<Rep_>::polylines_end()
{
  if(need_sr) {
    hot_pixels_list.clear();
    need_sr = false;
    find_hot_pixels_and_create_kd_trees();
    iterate();
  }    

  return(segments_output_list.end());
}

template<class Rep_>
bool Snap_rounding_2<Rep_>::remove(Segment_2 seg)
  {
    need_sr = true;
    bool found = false;
    Segment_data<Rep> s;

    for(typename std::list<Segment_data<Rep> >::iterator i1 = seg_list.begin();
      i1 != seg_list.end();++i1) {
       s = *i1;  
      if(s.equal(seg)) {
        found = true;
        seg_list.erase(i1);
        --number_of_segments;
        break;
      }
    }

    if(found) {
      for(Segment_iterator i2 = seg_2_list.begin();
        i2 != seg_2_list.end();++i2) {
        if(seg == *i2) {
          seg_2_list.erase(i2);
          break;
        }
      }
    }

    return(found);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::clear()
  { 
    need_sr = true;
    seg_list.clear();
    seg_2_list.clear();
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::change_number_of_kd_trees(
          int inp_number_of_kd_trees)
  {
    if(inp_number_of_kd_trees > 0) {
      number_of_kd_trees = inp_number_of_kd_trees;
      need_sr = true;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
bool Snap_rounding_2<Rep_>::change_pixel_size(NT inp_pixel_size)
  {
    if(inp_pixel_size > 0) {
      pixel_size = inp_pixel_size;
      need_sr = true;
      return(true);
    } else
      return(false);
  }

// @@@@ isrs
template<class Rep_>
bool Snap_rounding_2<Rep_>::change_delta(NT inp_delta)
  {
    if(inp_delta > 0) {
      delta = inp_delta;
      need_sr = true;
      return(true);
    } else
      return(false);
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::do_isr(bool inp_do_isr)
  { 
    wheteher_to_do_isr = inp_do_isr;
    need_sr = true;
  }

// @@@@ next fiunction
template<class Rep_>
void Snap_rounding_2<Rep_>::do_isrs(bool inp_do_isrs)
  { 
    wheteher_to_do_isrs = inp_do_isrs;
    need_sr = true;
  }

template<class Rep_>
void Snap_rounding_2<Rep_>::use_integer_mode(bool inp_int_output)
  { 
    int_output = inp_int_output;
    need_sr = true;
  }

template<class Rep>
typename Snap_rounding_2<Rep>::Direction Snap_rounding_2<Rep>::seg_dir;

template<class Rep>
bool Snap_rounding_2<Rep>::erase_hp = false;

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_H
