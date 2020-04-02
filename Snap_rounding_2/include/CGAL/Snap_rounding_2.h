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

#ifndef CGAL_SNAP_ROUNDING_2_H
#define CGAL_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>


#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <list>
#include <set>
#include <CGAL/Snap_rounding_kd_2.h>
#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/tss.h>

namespace CGAL {

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
  typedef CGAL::Segment_data<Traits>          Segment_data;

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


// a function for compare two hot pixels for the set of hot pixels
template<class Traits_>
struct Hot_pixel_auclidian_cmp {
  typedef CGAL::Hot_pixel<Traits_>    Hot_pixel;

  Traits_ m_gt;

  bool operator()(const Hot_pixel * h1, const Hot_pixel * h2) const;
};

// a function for compare two hot pixels for the set of hot pixels a
// certain segment intersect
template<class Traits_>
struct Hot_pixel_dir_cmp {
  typedef CGAL::Hot_pixel<Traits_>    Hot_pixel;

  Traits_ m_gt;

  bool operator()(const Hot_pixel * h1, const Hot_pixel * h2) const;
};

#ifdef CGAL_SR_DEBUG
int number_of_false_hp;
#endif

template<class Traits_,class OutputContainer>
class Snap_rounding_2 {
private:
  typedef Traits_                                       Traits;
  typedef typename Traits::FT                           NT;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename OutputContainer::value_type          Polyline_type;
  typedef CGAL::Hot_pixel<Traits_>                      Hot_pixel;
  typedef CGAL::Segment_data<Traits>                    Segment_data;
  typedef CGAL::Multiple_kd_tree<Traits,Hot_pixel *>    Multiple_kd_tree;
  typedef std::list<Segment_data>                       Segment_data_list;
  typedef CGAL::Hot_pixel_dir_cmp<Traits>               Hot_pixel_dir_cmp;
  typedef std::set<Hot_pixel *, Hot_pixel_dir_cmp>      Hot_pixel_set;

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
void Snap_rounding_2<Traits,OutputContainer>::
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
void Snap_rounding_2<Traits_,OutputContainer>::
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
void Snap_rounding_2<Traits_,OutputContainer>::
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
void Snap_rounding_2<Traits_,OutputContainer>::
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
void Snap_rounding_2<Traits_,OutputContainer>::
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
void snap_rounding_2(InputIterator begin,
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

  typedef CGAL::Hot_pixel<Traits>                     Hot_pixel;
  typedef CGAL::Segment_data<Traits>                  Segment_data;
  typedef CGAL::Multiple_kd_tree<Traits,Hot_pixel *>  Multiple_kd_tree;
  typedef std::list<Segment_data>                     Segment_data_list;

  Segment_data_list seg_list;
  Multiple_kd_tree * mul_kd_tree = nullptr;

  output_container.clear();
  // copy segments list
  while (begin != end) {
    seg_list.push_back(Segment_data(begin->source(), begin->target()));
    ++begin;
  }

  Snap_rounding_2<Traits,OutputContainer> s;

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

} //namespace CGAL

#endif
