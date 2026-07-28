// Copyright (c) 2001  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>

#ifndef CGAL_HOT_PIXEL_SNAP_ROUNDING_2_H
#define CGAL_HOT_PIXEL_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>


#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Named_function_parameters.h>
#include <list>
#include <set>
#include <CGAL/internal/Snap_rounding_kd_2.h>
#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/tss.h>

#include<CGAL/Hot_pixel_snap_rounding_traits_2.h>
#include <CGAL/internal/Snap_rounding_helpers.h>

namespace CGAL {

namespace Snap_Rounding_2::internal {
enum SEG_Direction {SEG_UP_RIGHT,SEG_UP_LEFT,SEG_DOWN_RIGHT,SEG_DOWN_LEFT,
                SEG_UP,SEG_DOWN,SEG_LEFT,SEG_RIGHT,SEG_POINT_SEG};

template<class Traits_>
class Segment_data {
private:
  typedef Traits_                               Traits;
  typedef typename Traits::FT                   NT;
  typedef typename Traits::Segment_2            Segment_2;
  typedef typename Traits::Point_2              Point_2;

private:
  Point_2 p;
  Point_2 q;

  Traits m_gt;

public:
  Segment_data();
  Segment_data(const Point_2 & p_inp,const Point_2 & q_inp);
  Segment_2 segment() const { return(Segment_2(p, q)); }
  const Point_2& source() const { return(p); }
  const Point_2& target() const { return(q); }
  inline void set_data(const Point_2 & inp_p,const Point_2 & inp_q);
  void determine_direction(SEG_Direction & seg_dir);
  bool equal(const Segment_2 & s);
  Segment_data(const Segment_data & other);
};

/*! */
template<class Traits_>
class Hot_pixel {
private:
  typedef Traits_                       Traits;
  typedef typename Traits::FT           NT;
  typedef typename Traits::Segment_2    Segment_2;
  typedef typename Traits::Point_2      Point_2;
  typedef CGAL::Snap_Rounding_2::internal::Segment_data<Traits>          Segment_data;

private:
  // p is the center of the hot pixel
  Point_2 p;
  Point_2 p_left;
  Point_2 p_right;
  Point_2 p_down;
  Point_2 p_up;
  Traits_ m_gt;
  NT pixel_size;
  Segment_2 *right_seg;
  Segment_2 *left_seg;
  Segment_2 *top_seg;
  Segment_2 *bot_seg;

  static SEG_Direction& direction()
  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(SEG_Direction, seg_dir,SEG_UP);
    return seg_dir;
  }



public:
  Hot_pixel(const Point_2 & inp_point, NT inp_pixel_size);
  ~Hot_pixel();
  inline Point_2 get_center() const;
  inline Point_2 get_center(bool int_output) const;
  bool intersect_left(const Segment_2 & seg, SEG_Direction seg_dir) const;
  bool intersect_right(const Segment_2 & seg, SEG_Direction seg_dir) const;
  bool intersect_bot(const Segment_2 & seg, SEG_Direction seg_dir) const;
  bool intersect_top(const Segment_2 & seg, SEG_Direction seg_dir) const;
  bool intersect(Segment_data & seg, SEG_Direction seg_dir) const;
  void set_direction(SEG_Direction inp_seg_dir) { direction() = inp_seg_dir; }
  SEG_Direction get_direction() const { return direction(); }
};

/*! */


// a function to compare two hot pixels for the set of hot pixels
template<class Traits_>
struct Hot_pixel_auclidian_cmp {
  typedef CGAL::Snap_Rounding_2::internal::Hot_pixel<Traits_>    Hot_pixel;

  Traits_ m_gt;

  bool operator()(const Hot_pixel * h1, const Hot_pixel * h2) const;
};

// a function to compare two hot pixels for the set of hot pixels a
// certain segment intersect
template<class Traits_>
struct Hot_pixel_dir_cmp {
  typedef CGAL::Snap_Rounding_2::internal::Hot_pixel<Traits_>    Hot_pixel;

  Traits_ m_gt;

  bool operator()(const Hot_pixel * h1, const Hot_pixel * h2) const;
};

#ifdef CGAL_SR_DEBUG
int number_of_false_hp;
#endif

template<class Traits_,class OutputContainer>
class Hot_pixel_snap_rounding_2 {
private:
  typedef Traits_                                               Traits;
  typedef typename Traits::FT                                   NT;
  typedef typename Traits::X_monotone_curve_2                   X_monotone_curve_2;
  typedef typename OutputContainer::value_type                  Polyline_type;
  typedef CGAL::Snap_Rounding_2::internal::Hot_pixel<Traits_>                    Hot_pixel;
  typedef CGAL::Snap_Rounding_2::internal::Segment_data<Traits>                  Segment_data;
  typedef CGAL::Snap_Rounding_2::internal::Multiple_kd_tree<Traits,Hot_pixel *>  Multiple_kd_tree;
  typedef std::list<Segment_data>                                                Segment_data_list;
  typedef CGAL::Snap_Rounding_2::internal::Hot_pixel_dir_cmp<Traits>             Hot_pixel_dir_cmp;
  typedef std::set<Hot_pixel *, Hot_pixel_dir_cmp>                               Hot_pixel_set;

public:
  // friend class Segment_data<Traits>;
  // friend class Hot_pixel<Traits>;
  // friend struct Hot_pixel_dir_cmp<Traits>;

  typedef typename Traits::Segment_2                    Segment_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef std::list<Point_2>                            Point_list;

  void find_hot_pixels_and_create_kd_trees(NT pixel_size,
                                           unsigned int number_of_kd_trees,
                                           Segment_data_list & seg_list,
                                           Multiple_kd_tree ** mul_kd_tree);

  void iterate(OutputContainer & output_container,
               NT pixel_size, bool int_output, bool do_isr,
               Segment_data_list & seg_list,
               Multiple_kd_tree * mul_kd_tree);

private:
  Traits m_gt;

  static const int default_number_of_kd_trees = 1;

  void find_intersected_hot_pixels(Segment_data & seg,
                                   Hot_pixel_set &hot_pixels_intersected_set,
                                   int &number_of_intersections,
                                   NT pixel_size,
                                   Multiple_kd_tree * mul_kd_tree);

  void reroute_sr(Hot_pixel_set & inp_hot_pixels_intersected_set,
                  Polyline_type & seg_output,
                  bool int_output);

  void reroute_isr(Hot_pixel_set & inp_hot_pixels_intersected_set,
                   Polyline_type & seg_output,
                   int number_of_intersections,
                   bool first_time,
                   NT pixel_size,
                   bool int_output,
                   Multiple_kd_tree * mul_kd_tree);
};

/*! Constructor */
template<class Traits_>
Segment_data<Traits_>::Segment_data() {}

/*! Constructor */
template<class Traits_>
Segment_data<Traits_>::Segment_data(const Point_2 & p_inp,
                                    const Point_2 & q_inp) :
  p(p_inp), q(q_inp) {}

/*! Constructor */
template<class Traits_>
Segment_data<Traits_>::Segment_data(const Segment_data & other)
{
  p = other.p;
  q = other.q;
}

/*! */
template<class Traits_>
inline void Segment_data<Traits_>::set_data(const Point_2 & inp_p,
                                            const Point_2 & inp_q)
{
  p = inp_p;
  q = inp_q;
}

/*! */
template<class Traits_>
bool Segment_data<Traits_>::equal(const Segment_2 & s)
{
  typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;
  Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();
  return(construct_vertex(s, 0) == p && construct_vertex(s, 1) == q);
}

/*! */
template<class Traits_>
void Segment_data<Traits_>::determine_direction(SEG_Direction & seg_dir)
{
  typedef typename Traits_::Compare_y_2         Compare_y_2;
  typedef typename Traits_::Compare_x_2         Compare_x_2;

  Compare_x_2 compare_x = m_gt.compare_x_2_object();
  Compare_y_2 compare_y = m_gt.compare_y_2_object();

  Comparison_result cx = compare_x(p, q);
  Comparison_result cy = compare_y(p, q);

  if (cx == SMALLER) {
    if (cy == SMALLER) seg_dir = SEG_UP_RIGHT;
    else if (cy == EQUAL) seg_dir = SEG_RIGHT;
    else seg_dir = SEG_DOWN_RIGHT;
  } else if (cx == EQUAL) {
    if (cy == SMALLER) seg_dir = SEG_UP;
    else if (cy == EQUAL) seg_dir = SEG_POINT_SEG;
    else seg_dir = SEG_DOWN;
  } else {
    if (cy == SMALLER) seg_dir = SEG_UP_LEFT;
    else if (cy == EQUAL) seg_dir = SEG_LEFT;
    else seg_dir = SEG_DOWN_LEFT;
  }
}

// intersection pixel
template<class Traits_>
Hot_pixel<Traits_>::Hot_pixel(const Point_2 & inp_point, NT inp_pixel_size) :
  pixel_size(inp_pixel_size)
{
  NT x,y;
  m_gt.snap_2_object()(inp_point, pixel_size, x, y);

  NT left_coord = x - pixel_size / NT(2.0);
  NT right_coord = x + pixel_size / NT(2.0);
  NT bottom_coord = y - pixel_size / NT(2.0);
  NT top_coord = y + pixel_size / NT(2.0);
  p = Point_2(x, y);
  p_left = Point_2(left_coord, y);
  p_right = Point_2(right_coord, y);
  p_down = Point_2(x, bottom_coord);
  p_up = Point_2(x, top_coord);

  typedef typename Traits::Construct_segment_2  Construct_segment_2;
  Construct_segment_2 construct_seg = m_gt.construct_segment_2_object();

  Point_2 lb(left_coord, bottom_coord);
  Point_2 rb(right_coord, bottom_coord);
  Point_2 lt(left_coord, top_coord);
  Point_2 rt(right_coord, top_coord);
  right_seg = new Segment_2(construct_seg(rb, rt));
  left_seg = new Segment_2(construct_seg(lb, lt));
  top_seg = new Segment_2(construct_seg(lt, rt));
  bot_seg = new Segment_2(construct_seg(lb, rb));
}

/*! */
template<class Traits_>
Hot_pixel<Traits_>::~Hot_pixel()
{
  delete right_seg;
  delete left_seg;
  delete top_seg;
  delete bot_seg;
}

/*! */
template<class Traits_>
inline typename Hot_pixel<Traits_>::Point_2
Hot_pixel<Traits_>::get_center() const
{
  return(p);
}

/*! */
template<class Traits_>
inline typename Hot_pixel<Traits_>::Point_2
Hot_pixel<Traits_>::get_center(bool int_output) const
{
  if (int_output) {
    Point_2 out_p = m_gt.integer_grid_point_2_object()(p,pixel_size);
    return(out_p);
  } else
    return(p);
}

/*! */
template<class Traits_>
bool Hot_pixel<Traits_>::intersect_left(const Segment_2 & seg,
                                        SEG_Direction seg_dir) const
{
  typedef typename Traits_::Compare_y_2         Compare_y_2;
  typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;

  Object result;
  Point_2 p;
  Segment_2 s;

  result = intersection(seg, *left_seg);

  if (assign(p, result)) {
    Compare_y_2 compare_y = m_gt.compare_y_2_object();
    Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();

    Comparison_result c_p = compare_y(p, p_up);
    Comparison_result c_target = compare_y(construct_vertex(seg, 1), p_up);
    Comparison_result c_source = compare_y(construct_vertex(seg, 0), p_up);

    return(c_p != EQUAL || (seg_dir == SEG_UP_LEFT && c_source != EQUAL) ||
           (seg_dir == SEG_DOWN_RIGHT && c_target != EQUAL));
  } else if (assign(s,result))
    return(true);
  else
    return(false);
}

/*! */
template<class Traits_>
bool Hot_pixel<Traits_>::intersect_right(const Segment_2 & seg,
                                         SEG_Direction seg_dir) const
{
  typedef typename Traits_::Compare_y_2         Compare_y_2;
  typedef typename Traits_::Compare_x_2         Compare_x_2;
  typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;

  Object result;
  Point_2 p;
  Segment_2 s;

  result = intersection(seg, *right_seg);

  if (assign(p,result)) {
    // bottom right point was checked in intersect_bot
    Compare_y_2 compare_y = m_gt.compare_y_2_object();
    Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();

    Comparison_result c1 = compare_y(p, p_up);
    Comparison_result c2 = compare_y(p, p_down);

    Point_2 src = construct_vertex(seg, 0);
    Point_2 trg = construct_vertex(seg, 1);
    Comparison_result c3 = compare_y(src, p_up);
    Comparison_result c4 = compare_y(trg, p_up);

    if (c1 == EQUAL)
      return ((seg_dir == SEG_UP_RIGHT && c3 != EQUAL) ||
              (seg_dir == SEG_DOWN_LEFT && c4 != EQUAL));
    else if (c2 == EQUAL)
      return false;// was checked

    Compare_x_2 compare_x = m_gt.compare_x_2_object();
    Comparison_result c_target = compare_x(p_right, trg);
    Comparison_result c_source = compare_x(p_right, src);

    return (((seg_dir == SEG_LEFT || seg_dir == SEG_DOWN_LEFT ||
              seg_dir == SEG_UP_LEFT) && c_target != EQUAL) ||
            ((seg_dir == SEG_RIGHT || seg_dir == SEG_DOWN_RIGHT ||
              seg_dir == SEG_UP_RIGHT) && c_source != EQUAL));
  }
  return false;
}

/*! */
template<class Traits_>
bool Hot_pixel<Traits_>::intersect_bot(const Segment_2 & seg,
                                       SEG_Direction seg_dir) const
{
  typedef typename Traits_::Compare_x_2         Compare_x_2;
  typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;

  Object result;
  Point_2 p;
  Segment_2 s;

  result = intersection(seg,*bot_seg);

  if (assign(p,result)) {
    Compare_x_2 compare_x = m_gt.compare_x_2_object();
    Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();

    Comparison_result c_p = compare_x(p, p_right);
    Comparison_result c_target = compare_x(construct_vertex(seg, 1), p_right);
    Comparison_result c_source = compare_x(construct_vertex(seg, 0), p_right);

    return(c_p != EQUAL || (seg_dir == SEG_UP_LEFT && c_target != EQUAL) ||
           (seg_dir == SEG_DOWN_RIGHT && c_source != EQUAL));
  } else if (assign(s,result))
    return(true);
  else
    return(false);
}

/*! */
template<class Traits_>
bool Hot_pixel<Traits_>::intersect_top(const Segment_2 & seg,
                                       SEG_Direction seg_dir) const
{
  typedef typename Traits_::Compare_x_2         Compare_x_2;
  typedef typename Traits_::Compare_y_2         Compare_y_2;
  typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;

  Object result;
  Point_2 p;
  Segment_2 s;

  result = intersection(seg, *top_seg);

  if (assign(p,result))
  {
    Compare_x_2 compare_x = m_gt.compare_x_2_object();
    Compare_y_2 compare_y = m_gt.compare_y_2_object();
    Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();

    // corner points was checked in intersect_bot
    Comparison_result c1 = compare_x(p, p_left);
    Comparison_result c2 = compare_x(p, p_right);
    Comparison_result c3 = compare_y(construct_vertex(seg, 1), p_up);
    Comparison_result c4 = compare_y(construct_vertex(seg, 0), p_up);

    if (c1 == EQUAL || c2 == EQUAL)
      return(false);// were checked
    else
      return(((seg_dir == SEG_DOWN || seg_dir == SEG_DOWN_LEFT ||
               seg_dir == SEG_DOWN_RIGHT) && c3 != EQUAL) ||
             ((seg_dir == SEG_UP || seg_dir == SEG_UP_LEFT ||
               seg_dir == SEG_UP_RIGHT) && c4 != EQUAL));
  }
  return(false);
}

/*! */
template<class Traits_>
bool
Hot_pixel<Traits_>::intersect(Segment_data & seg, SEG_Direction my_seg_dir) const
{
  Segment_2 s = seg.segment();

  return(intersect_bot(s,my_seg_dir) || intersect_left(s, my_seg_dir) ||
         intersect_right(s, my_seg_dir) || intersect_top(s, my_seg_dir));
}

// a function for compare two hot pixels for the set of hot pixels
template<class Traits_>
bool Hot_pixel_auclidian_cmp<Traits_>::
operator()(const Hot_pixel * h1, const Hot_pixel * h2) const
{
  typedef typename Traits_::Compare_x_2         Compare_x_2;
  typedef typename Traits_::Compare_y_2         Compare_y_2;

  Compare_x_2 compare_x = m_gt.compare_x_2_object();
  Compare_y_2 compare_y = m_gt.compare_y_2_object();

  Comparison_result cx = compare_x(h1->get_center(), h2->get_center());
  Comparison_result cy = compare_y(h1->get_center(), h2->get_center());

  return(cx == SMALLER || ( cx == EQUAL && cy == SMALLER));
}

// a function for compare two hot pixels for the set of hot pixels a certain
// segment intersect
template<class Traits_>
bool Hot_pixel_dir_cmp<Traits_>::operator ()(const Hot_pixel * h1,
                                             const Hot_pixel * h2) const
{
  typedef typename Traits_::Compare_x_2         Compare_x_2;
  typedef typename Traits_::Compare_y_2         Compare_y_2;

  Compare_x_2 compare_x = m_gt.compare_x_2_object();
  Compare_y_2 compare_y = m_gt.compare_y_2_object();

  Comparison_result cx = compare_x(h1->get_center(), h2->get_center());
  Comparison_result cy = compare_y(h1->get_center(), h2->get_center());
  SEG_Direction seg_dir = h1->get_direction();

  // Point segment intersects only one pixel, thus ignored
  return((seg_dir == SEG_UP_RIGHT &&  (cx == SMALLER || (cx == EQUAL &&
                                                         cy == SMALLER))) ||
         (seg_dir == SEG_UP_LEFT && (cx == LARGER || (cx == EQUAL &&
                                                      cy == SMALLER))) ||
         (seg_dir == SEG_DOWN_RIGHT && (cx == SMALLER || (cx == EQUAL &&
                                                          cy == LARGER))) ||
         (seg_dir == SEG_DOWN_LEFT && (cx == LARGER || (cx == EQUAL &&
                                                        cy == LARGER))) ||
         (seg_dir == SEG_UP && cy == SMALLER) ||
         (seg_dir == SEG_DOWN && cy == LARGER) ||
         (seg_dir == SEG_LEFT && cx == LARGER) ||
         (seg_dir == SEG_RIGHT && cx == SMALLER));
}

/*! */
template<class Traits,class OutputContainer>
void Hot_pixel_snap_rounding_2<Traits,OutputContainer>::
find_hot_pixels_and_create_kd_trees(NT pixel_size,
                                    unsigned int number_of_kd_trees,
                                    std::list<Segment_data> & seg_list,
                                    Multiple_kd_tree ** mul_kd_tree)
{
  typedef std::pair<Point_2, Hot_pixel *>               Point_hot_pixel_pair;
  typedef typename std::list<Segment_data>::iterator    Segment_data_iter;
  typedef std::list<Segment_2>                          Segment_list;
  typedef typename std::list<Point_2>::const_iterator   Point_const_iter;

  typedef typename Traits::Construct_segment_2  Construct_segment_2;
  Construct_segment_2 construct_seg = m_gt.construct_segment_2_object();

  Hot_pixel * hp;
  Segment_data_iter iter1;
  Object result;
  Point_2 p;
  std::list<Point_hot_pixel_pair> hot_pixels_list;
  Segment_list segments;

  if (seg_list.empty()) return;

  for (iter1 = seg_list.begin(); iter1 != seg_list.end(); ++iter1) {
    segments.push_back(construct_seg(iter1->source(), iter1->target()));
  }

  // get intersection points (with endpoints)
  Point_list mypointlist;

  CGAL::compute_intersection_points(segments.begin(), segments.end(),
                             std::back_inserter(mypointlist), true, m_gt);

  for (Point_const_iter v_iter = mypointlist.begin();
       v_iter != mypointlist.end(); ++v_iter)
  {
    hp = new Hot_pixel(*v_iter, pixel_size);
    hot_pixels_list.push_back(Point_hot_pixel_pair(hp->get_center(), hp));
  }

  // create kd multiple tree
  // create simple_list from seg_list
  Segment_list simple_seg_list;
  for (Segment_data_iter iter = seg_list.begin(); iter != seg_list.end();
       ++iter)
  {
    simple_seg_list.push_back(construct_seg(iter->source(), iter->target()));
  }
  *mul_kd_tree = new Multiple_kd_tree(hot_pixels_list, number_of_kd_trees,
                                      simple_seg_list);
}

/*! */
template<class Traits_,class OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_,OutputContainer>::
find_intersected_hot_pixels(Segment_data & seg,
                            std::set<Hot_pixel *,
                            Hot_pixel_dir_cmp> & hot_pixels_intersected_set,
                            int &number_of_intersections,
                            NT pixel_size,
                            Multiple_kd_tree * mul_kd_tree)
{
  typedef typename std::list<Hot_pixel *>::iterator     Hot_pixel_iter;

  Hot_pixel_iter iter;
  SEG_Direction seg_dir;

  hot_pixels_intersected_set.clear();
  seg.determine_direction(seg_dir);
  number_of_intersections = 0;
  std::list<Hot_pixel *> hot_pixels_list;

  mul_kd_tree->get_intersecting_points(hot_pixels_list,
                                       Segment_2(seg.segment()), pixel_size);

  for (iter = hot_pixels_list.begin();iter != hot_pixels_list.end();++iter) {
    if ((*iter)->intersect(seg,seg_dir)) {
      (*iter)->set_direction(seg_dir);
      hot_pixels_intersected_set.insert(*iter);
    }

#ifdef CGAL_SR_DEBUG
    else
      ++number_of_false_hp;
#endif

  }

  number_of_intersections = static_cast<int>(hot_pixels_intersected_set.size());
}

/*! */
template<class Traits_,class OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_,OutputContainer>::
reroute_sr(std::set<Hot_pixel *, Hot_pixel_dir_cmp>
           & inp_hot_pixels_intersected_set,
           Polyline_type & seg_output,
           bool int_output)
{
  typedef typename std::set<Hot_pixel *, Hot_pixel_dir_cmp>::iterator
    Hot_pixel_iter;

  Hot_pixel_iter hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
  ++hot_pixel_iter;

  while (hot_pixel_iter != inp_hot_pixels_intersected_set.end()) {
    seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
    ++hot_pixel_iter;
  }
}

/*! */
template<class Traits_,class OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_,OutputContainer>::
reroute_isr(std::set<Hot_pixel *, Hot_pixel_dir_cmp>
            & inp_hot_pixels_intersected_set,
            Polyline_type & seg_output,
            int number_of_intersections,
            bool first_time,
            NT pixel_size,
            bool int_output,
            Multiple_kd_tree * mul_kd_tree)
{
  typedef std::set<Hot_pixel *, Hot_pixel_dir_cmp>      Hot_pixel_set;
  typedef typename std::set<Hot_pixel *, Hot_pixel_dir_cmp>::iterator
    Hot_pixel_iter;

  Hot_pixel_iter hot_pixel_iter, next_hot_pixel_iter,
    before_last_hot_pixel_iter;
  Segment_data seg;
  Hot_pixel_set hot_pixels_intersected_set;
  SEG_Direction seg_dir;

  if (number_of_intersections > 2 || first_time) {
    before_last_hot_pixel_iter = inp_hot_pixels_intersected_set.end();
    --before_last_hot_pixel_iter;

    for (hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
         hot_pixel_iter != before_last_hot_pixel_iter; ++hot_pixel_iter)
    {
      next_hot_pixel_iter = hot_pixel_iter;
      ++next_hot_pixel_iter;
      seg.set_data((*hot_pixel_iter)->get_center(),
                   (*next_hot_pixel_iter)->get_center());
      seg.determine_direction(seg_dir);
      find_intersected_hot_pixels(seg, hot_pixels_intersected_set,
                                  number_of_intersections, pixel_size,
                                  mul_kd_tree);
      reroute_isr(hot_pixels_intersected_set, seg_output,
                  number_of_intersections,false, pixel_size,
                  int_output, mul_kd_tree);
    }
  } else {
    // insert second hot pixel
    hot_pixel_iter = inp_hot_pixels_intersected_set.begin();
    ++hot_pixel_iter;
    seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
  }
}

/*! */
template<class Traits_, class OutputContainer>
void Hot_pixel_snap_rounding_2<Traits_,OutputContainer>::
iterate(OutputContainer & output_container,
        NT pixel_size,
        bool int_output, bool do_isr,
        std::list<Segment_data> & seg_list,
        Multiple_kd_tree * mul_kd_tree)
{
  typedef std::set<Hot_pixel *, Hot_pixel_dir_cmp>      Hot_pixel_set;
  typedef typename std::set<Hot_pixel *, Hot_pixel_dir_cmp>::iterator
    Hot_pixel_iter;
  typedef typename std::list<Segment_data>::iterator    Segment_data_iter;

  Polyline_type seg_output;
  Hot_pixel_set hot_pixels_intersected_set;
  Hot_pixel_iter hot_pixel_iter;
  int number_of_intersections;
  Hot_pixel * hp;
  SEG_Direction seg_dir;

  for (Segment_data_iter iter = seg_list.begin(); iter != seg_list.end();
       ++iter)
  {
    seg_output.clear();
    iter->determine_direction(seg_dir);
    find_intersected_hot_pixels(*iter, hot_pixels_intersected_set,
                                number_of_intersections, pixel_size,
                                mul_kd_tree);
    // hot_pixels_intersected_set must have at least two hot pixels when the
    // segment is not in entirely inside a hot pixel enter first hot pixel
    hot_pixel_iter = hot_pixels_intersected_set.begin();
    if (hot_pixel_iter == hot_pixels_intersected_set.end()) {
      // segment entirely inside a pixel
      hp = new Hot_pixel(iter->source(), pixel_size);
      seg_output.push_back(hp->get_center(int_output));
      delete hp;
    } else {
      seg_output.push_back((*hot_pixel_iter)->get_center(int_output));
      if (number_of_intersections > 1) {
        // segments that have at most one intersecting hot pixel are
        // done(it was inserted)
        if (do_isr)
          reroute_isr(hot_pixels_intersected_set, seg_output,
                      number_of_intersections, true,pixel_size,
                      int_output, mul_kd_tree);
        else
          reroute_sr(hot_pixels_intersected_set, seg_output, int_output);
      }
    }

    output_container.push_back(seg_output);
  }
}

/*! */
template<class Traits, class InputIterator, class OutputContainer>
void hot_pixel_snap_rounding_2(InputIterator begin,
                               InputIterator end,
                               OutputContainer & output_container,
                               typename Traits::NT pixel_size,
                               bool do_isr = true,
                               bool int_output = true,
                               unsigned int number_of_kd_trees = 1)
{
#ifdef CGAL_SR_DEBUG
  number_of_false_hp = 0;
#endif

  typedef Hot_pixel<Traits>                     Hot_pixel;
  typedef Segment_data<Traits>                  Segment_data;
  typedef Multiple_kd_tree<Traits,Hot_pixel *>  Multiple_kd_tree;
  typedef std::list<Segment_data>               Segment_data_list;

  Segment_data_list seg_list;
  Multiple_kd_tree * mul_kd_tree = nullptr;

  output_container.clear();
  // copy segments list
  while (begin != end) {
    seg_list.push_back(Segment_data(begin->source(), begin->target()));
    ++begin;
  }

  Hot_pixel_snap_rounding_2<Traits,OutputContainer> s;

  s.find_hot_pixels_and_create_kd_trees(pixel_size, number_of_kd_trees,
                                        seg_list, &mul_kd_tree);
  s.iterate(output_container, pixel_size, int_output, do_isr,seg_list,
            mul_kd_tree);

  // hope that find_hot_pixels_and_create_kd_trees does not suddenly
  // new up an array of multiple kd_tree
  delete mul_kd_tree;

#ifdef CGAL_SR_DEBUG
  std::cout << "Overall number of false hot pixels in all the queries : "
            << number_of_false_hp << std::endl;
#endif

}

template <class PolygonRange, class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2_polygon(const PolygonRange  &polygons,
                                       OutputContainer     &out,
                                       const NamedParameters &np = parameters::default_values())
{
  using Polygon_2 = typename std::iterator_traits<typename PolygonRange::iterator>::value_type;
  using Kernel = typename Kernel_traits<typename Polygon_2::Point_2>::Kernel;
  using DefaultTraits = Hot_pixel_snap_rounding_traits_2<Kernel>;
  using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;
  using Point_2 = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector< Segment_2 > input_segments;
  std::vector< std::vector< Point_2 > > polylines;

  // Store the indices of segment of a new polygon, segments between [ polygon_index[i] and polygon_index[i+1] [ belong to polygon i
  std::vector< std::size_t > polygon_indices;

  polygon_indices.reserve(std::distance(polygons.begin(), polygons.end()));
  for(const Polygon_2 &P: polygons){
    polygon_indices.push_back(input_segments.size());
    for(std::size_t i=0; i<P.size()-1; ++i)
      input_segments.emplace_back(P[i], P[i+1]);
    input_segments.emplace_back(P[P.size()-1], P[0]);
  }
  polygon_indices.push_back(input_segments.size());

  // Main algorithm
  hot_pixel_snap_rounding(input_segments, polylines, np);

  // Reassemble the polygons
  for(std::size_t polygon_idx = 0; polygon_idx != polygon_indices.size()-1; ++polygon_idx){
    Polygon_2 P;
    std::size_t idx_start = polygon_indices[polygon_idx];
    std::size_t idx_end = polygon_indices[polygon_idx+1];
    Point_2 last_insert;
    for(std::size_t pl_idx = idx_start; pl_idx != idx_end; ++pl_idx){
      const auto &pl = polylines[pl_idx];
      // The first element is not add, it is identical to the last

      bool go_forward;
      if(pl_idx == idx_start)
        go_forward = (pl.back() == polylines[pl_idx+1].front()) || (pl.back() == polylines[pl_idx+1].back());
      else
        go_forward = (last_insert == pl.front());

      // Add the element in forward direction
      if(go_forward){
        for(std::size_t i = 1; i != pl.size(); ++i)
          P.push_back(pl[i]);
        last_insert = pl.back();

      // Add the element in backward direction
      } else {
        for(std::size_t i = pl.size()-1; i != 0; --i)
          P.push_back(pl[i-1]);
        last_insert = pl.front();
      }
    }
    out.push_back(P);
  }
}

} // namespace internal

#if DOXYGEN_RUNNING

/**
* \ingroup Snap_rounding_hot_pixel_grp
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
*
* By default, each polyline of the output corresponds to an input segment. Consequently, duplicate segments may appear in the output, for instance when multiple input segments collapse.
* When the parameter `output_unique_segments` is set to `true`, the polylines are decomposed into individual segments (represented as polylines with two points), and duplicates are removed.
*
* @tparam SegmentRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type is `geom_traits::Segment_2`, where the type of geom_traits is detailed by `np::geom_traits`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is itself a model of the concepts `DefaultConstructible` and `BackInsertionSequence` whose value type is `geom_traits::Point_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segment range
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{pixel_size}
*     \cgalParamDescription{The size of the pixel. The plane will be tiled with square pixels of that width such that the origin is the center of a pixel.}
*     \cgalParamType{`geom_traits::FT`}
*     \cgalParamDefault{1}
*   \cgalParamNEnd
*   \cgalParamNBegin{do_iterative_snap_rounding}
*     \cgalParamDescription{determines whether to apply Iterative Snap Rounding, see the user manual for more details.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
*   \cgalParamNBegin{use_grid_coordinates}
*     \cgalParamDescription{If set to true, the output coordinates are expressed in the integer grid (pixel indices).
*                           Otherwise, they are given in the input coordinate system.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
*   \cgalParamNBegin{output_unique_segments}
*     \cgalParamDescription{If set to true, the output polylines are unique pairs of distinct points represented a segment. As a result, the total number of output polylines may differ from the number of input segments.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{false}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must be a model of `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the segment type, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class SegmentRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const SegmentRange &segments,
                               OutputContainer    &out,
                               const NamedParameters &np = parameters::default_values());

/**
* \ingroup Snap_rounding_hot_pixel_grp
*
* \brief subdivides and rounds a range of polygons so that their boundary segments are pairwise disjoint in their interiors.
*
* If the input polygons are disjoint, the output polygons remain non-overlapping, although they may share vertices or edges.
* Each output polygon is free of self-intersections but may present pinched sections.
*
* @tparam PolygonRange model of a `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type is model of `CGAL::Polygon_2`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is model of `CGAL::Polygon_2`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param polygons the range of input polygons
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{pixel_size}
*     \cgalParamDescription{The size of the pixel. The plane will be tiled with square pixels of that width such that the origin is the center of a pixel.}
*     \cgalParamType{`geom_traits::FT`}
*     \cgalParamDefault{1}
*   \cgalParamNEnd
*   \cgalParamNBegin{do_iterative_snap_rounding}
*     \cgalParamDescription{determines whether to apply Iterative Snap Rounding, see the user manual for more details.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
*   \cgalParamNBegin{use_grid_coordinates}
*     \cgalParamDescription{If set to true, the output coordinates are expressed in the integer grid (pixel indices).
*                           Otherwise, they are given in the input coordinate system.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Hot_pixel_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the point type of the polygons, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* @warning a convex input polygon might no longer be convex after rounding.
*/
template <class PolygonRange, class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const PolygonRange  &polygons,
                               OutputContainer     &out,
                               const NamedParameters &np = parameters::default_values());

# else

template <class InputRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void hot_pixel_snap_rounding_2(const InputRange &inputs,
                               OutputContainer  &out,
                               const NamedParameters &np = parameters::default_values())
{
  using Input = std::remove_cv_t<typename std::iterator_traits<typename InputRange::iterator>::value_type>;

  if constexpr(internal::is_instance_of_Polygon_2< Input >){
    Snap_Rounding_2::internal::hot_pixel_snap_rounding_2_polygon(inputs, out, np);
    return;
  } else {
    using Polyline = std::remove_cv_t<typename OutputContainer::value_type>;

    using Kernel = typename Kernel_traits<std::remove_cv_t<typename std::iterator_traits<typename InputRange::iterator>::value_type>>::Kernel;
    using DefaultTraits = Hot_pixel_snap_rounding_traits_2<Kernel>;
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                              NamedParameters,
                                                              DefaultTraits>::type;

    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    auto pixel_size = choose_parameter(get_parameter(np, internal_np::pixel_size), 1.);
    bool do_isr = choose_parameter(get_parameter(np, internal_np::do_iterative_snap_rounding), true);
    bool int_output = choose_parameter(get_parameter(np, internal_np::use_grid_coordinates), true);
    bool unique_segments = choose_parameter(get_parameter(np, internal_np::output_unique_segments), false);
    unsigned int number_of_kd_trees = 1;

    if(unique_segments){
      // Output Segments while removing duplicate ones
      std::vector< std::vector< Point_2 > > output_container;
      Snap_Rounding_2::internal::hot_pixel_snap_rounding_2<Traits>(inputs.begin(), inputs.end(), output_container, pixel_size, do_isr, int_output, number_of_kd_trees);

      // Build a set with lexicographic order of the segments
      auto comp = Traits().compare_xy_2_object();
      auto less = Traits().less_xy_2_object();
      auto pred = [&](const Segment_2 &a, const Segment_2 &b){
        auto res = comp(a.source(), b.source());
        if(res != EQUAL)
          return res == SMALLER;
        return less(a.target(), b.target());
      };
      std::vector< Segment_2 > out_segs;
      for(auto &poly: output_container)
        for(std::size_t i=1; i<poly.size(); ++i)
          if(less(poly[i-1], poly[i]))
            out_segs.emplace_back(poly[i-1], poly[i]);
          else
            out_segs.emplace_back(poly[i], poly[i-1]);
      std::sort(out_segs.begin(), out_segs.end(), pred);
      auto last = std::unique(out_segs.begin(), out_segs.end());
      out_segs.erase(last, out_segs.end());
      for(Segment_2 &s: out_segs){
        Polyline pl;
        pl.push_back(std::move(s.source()));
        pl.push_back(std::move(s.target()));
        out.push_back(std::move(pl));
      }
    } else {
      // Output polylines
      Snap_Rounding_2::internal::hot_pixel_snap_rounding_2<Traits>(inputs.begin(), inputs.end(), out, pixel_size, do_isr, int_output, number_of_kd_trees);
    }
 }
}

#endif


} //namespace CGAL

#endif
