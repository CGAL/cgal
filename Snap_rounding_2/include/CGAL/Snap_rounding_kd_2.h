// Copyright (c) 2001, 2009, 2014 Tel-Aviv University (Israel), Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// author(s)     : Eli Packer <elip@post.tau.ac.il>,
//                 Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_SNAP_ROUNDING_KD_2_H
#define CGAL_SNAP_ROUNDING_KD_2_H

#include <list>
#include <CGAL/basic.h>
#include <CGAL/predicates_on_points_2.h>
#include <iostream>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/utility.h>
#include <CGAL/assertions.h>
#include <CGAL/Dimension.h>

#include <boost/type_traits/is_pointer.hpp>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Fuzzy_iso_box.h>

namespace CGAL {

namespace internal {

//////////////////////
//////////////////////
//Point_with_hot_pixel_history
//////////////////////

template<class Traits, class SAVED_OBJECT>
class Point_with_hot_pixel_history : public Traits::Point_2 {

private:

  typedef typename Traits::Point_2                                Base;
  typedef typename Traits::Point_2                                Point_2;
  typedef typename Traits::FT                                     NT;

public:

  Point_2 orig;
  SAVED_OBJECT object;

  Point_with_hot_pixel_history(const Base& p, const Point_2& inp_orig, SAVED_OBJECT obj) : Base(p), orig(inp_orig), object(obj) {}

  Point_with_hot_pixel_history(const Base& p) : Base(p), orig(Point_2(0, 0)) {}

  Point_with_hot_pixel_history() : Base(), orig() {}

  Point_with_hot_pixel_history(NT x, NT y) : Base(x, y), orig(Point_2(0, 0)) {}

}; // Point_with_hot_pixel_history


//////////////////////
//////////////////////
//Search_traits_kd_tree_2
//
//(Search traits modified to be used by the Spacial Searching kd_trees for Snap rounding)
//////////////////////

template < class Traits_, class Point_ = typename Traits_::Point_2 >
class Search_traits_kd_tree_2 {

public:
  typedef Traits_                                                Traits;
  typedef Point_                                                 Point_d;

  typedef Dimension_tag<2>                                       Dimension;
  typedef typename Traits::Iso_rectangle_2                       Iso_box_d;
  typedef typename Traits::Cartesian_const_iterator_2            Cartesian_const_iterator_d;
  typedef typename Traits::Construct_cartesian_const_iterator_2  Construct_cartesian_const_iterator_d;

  typedef typename Traits::Construct_min_vertex_2                Construct_min_vertex_d;
  typedef typename Traits::Construct_max_vertex_2                Construct_max_vertex_d;

  typedef typename Traits::Construct_iso_rectangle_2             Construct_iso_box_d;
  typedef typename Traits::FT                                    FT;

  Construct_cartesian_const_iterator_d construct_cartesian_const_iterator_d_object() const
  {
    return Construct_cartesian_const_iterator_d();
  }

}; // Search_traits_kd_tree_2

} // namespace internal

/////////////////////
/////////////////////
//Multiple_kd_tree
/////////////////////

template<class Traits_, class SAVED_OBJECT>
class Multiple_kd_tree {
  CGAL_static_assertion_msg((boost::is_pointer<SAVED_OBJECT>::value), "SAVED_OBJECT is not a pointer.");
private:
  typedef Traits_                                       Traits;
  typedef typename Traits::FT                           NT;
  typedef typename Traits::Segment_2                    Segment_2;
  typedef typename Traits::Point_2                      Point_2;
  typedef typename Traits::Vector_2                     Vector_2;
  typedef typename Traits::Iso_rectangle_2              Iso_rectangle_2;
  typedef typename Traits::Direction_2                  Direction_2;
  typedef typename Traits::Line_2                       Line_2;
  typedef typename Traits::Aff_transformation_2         Transformation_2;

  typedef CGAL::internal::Point_with_hot_pixel_history<Traits, SAVED_OBJECT>     Point_with_hot_pixel_history_saved;
  typedef CGAL::internal::Search_traits_kd_tree_2<Traits, Point_with_hot_pixel_history_saved>
                                                                                 Search_traits;
  typedef CGAL::Kd_tree<Search_traits>                                           Kd_tree;
  typedef CGAL::Fuzzy_iso_box<Search_traits>                                     Box;

  typedef std::list<Point_with_hot_pixel_history_saved>   Points_List;
  typedef std::pair<Direction_2, NT>                      Direction_nt_pair;
  typedef std::pair<Kd_tree *,Direction_nt_pair>          Kd_triple;
  typedef std::pair<Kd_tree *, Direction_nt_pair>         Kd_direction_nt_pair;
  typedef std::list<Kd_direction_nt_pair>                 Kd_triple_list;

  typedef std::pair<Point_2, SAVED_OBJECT>              Point_saved_pair;
  typedef std::list<Point_saved_pair>                   Point_saved_pair_list;
  typedef typename Point_saved_pair_list::iterator      Point_saved_pair_iter;

  typedef typename std::list<Point_with_hot_pixel_history_saved>            Point_with_hot_pixel_history_saved_list;
  typedef typename Point_with_hot_pixel_history_saved_list::iterator        Point_with_hot_pixel_history_saved_iter;

  typedef std::list<Point_2>                            Point_list;
  typedef typename Point_list::iterator                 Point_iter;

  typedef std::list<Segment_2>                          Segment_list;
  typedef typename Segment_list::const_iterator         Segment_const_iter;

  typedef std::list<Direction_2>                        Direction_list;
  typedef typename Direction_list::const_iterator       Direction_const_iter;


private:
  Traits m_gt;
  const double pi, half_pi;
  int number_of_trees;

  Kd_triple_list    kd_trees_list;

  Point_saved_pair_list input_points_list;
  std::map<int, NT> angle_to_sines_appr; // was const int

  /*! */
  void rotate(Point_2& p, NT angle)
  {
    static const double rad_to_deg = 57.297;

    typename Traits::To_double to_dbl;
    int tranc_angle = int(to_dbl(angle) * rad_to_deg);

    NT cosine_val = angle_to_sines_appr[90 - tranc_angle],
       sine_val = angle_to_sines_appr[tranc_angle];

    Transformation_2 rotate(ROTATION, sine_val, cosine_val);
    p = rotate(p);
  }


  /*! */
  Kd_triple create_kd_tree(NT angle)
  {

    Kd_tree *tree = new Kd_tree();

    tree->reserve(input_points_list.size());

    for (Point_saved_pair_iter iter = input_points_list.begin();  iter != input_points_list.end();  ++iter)
    {
      Point_2 p(iter->first);
      rotate(p,angle);
      Point_with_hot_pixel_history_saved rotated_point(p,iter->first,iter->second);

      tree->insert(rotated_point);
    }

    tree->build();


    typename Traits::To_double to_dbl;
    double buffer_angle(to_dbl(angle) - half_pi / (2 * number_of_trees));

    if (buffer_angle < 0)
      buffer_angle = 0;

    Line_2 li(std::tan(buffer_angle), -1, 0);
    Direction_2 d(li);

    // rotate_by 180 degrees
    Transformation_2 t(ROTATION, 0, -1);
    d = d.transform(t);
    Direction_nt_pair kp(d, angle);

    Kd_triple kt(tree, kp);

    return(kt);
  }

  inline NT square(NT x) {return(x * x);}

  inline NT min BOOST_PREVENT_MACRO_SUBSTITUTION (NT x, NT y) {return((x < y) ? x : y);}
  inline NT max BOOST_PREVENT_MACRO_SUBSTITUTION (NT x, NT y) {return((x < y) ? y : x);}

  inline NT min BOOST_PREVENT_MACRO_SUBSTITUTION  (NT x1, NT x2, NT x3, NT x4, NT x5,
                NT x6)
  {return(min BOOST_PREVENT_MACRO_SUBSTITUTION (min BOOST_PREVENT_MACRO_SUBSTITUTION (min BOOST_PREVENT_MACRO_SUBSTITUTION (x1, x2),
                                                min BOOST_PREVENT_MACRO_SUBSTITUTION (x3, x4)),min BOOST_PREVENT_MACRO_SUBSTITUTION (x5, x6)));}

  inline NT max BOOST_PREVENT_MACRO_SUBSTITUTION (NT x1, NT x2, NT x3, NT x4, NT x5, NT x6)
  {return(max BOOST_PREVENT_MACRO_SUBSTITUTION (max BOOST_PREVENT_MACRO_SUBSTITUTION (max BOOST_PREVENT_MACRO_SUBSTITUTION (x1, x2),
                                                max BOOST_PREVENT_MACRO_SUBSTITUTION (x3, x4)),max BOOST_PREVENT_MACRO_SUBSTITUTION (x5, x6)));}

  /*! */
  Direction_2 get_direction(Segment_2 seg)
  {
    typedef typename Traits_::Construct_vertex_2  Construct_vertex_2;
    Construct_vertex_2 construct_vertex = m_gt.construct_vertex_2_object();

    // force the segment slope to [0-180)
    Point_2 s = construct_vertex(seg, 0);
    Point_2 t = construct_vertex(seg, 1);
    Comparison_result cx = m_gt.compare_x_2_object()(s, t);
    Comparison_result cy = m_gt.compare_y_2_object()(s, t);
    if (cy == LARGER || (cy == EQUAL && cx == LARGER)) {
      typedef typename Traits::Construct_segment_2  Construct_segment_2;
      Construct_segment_2 construct_seg = m_gt.construct_segment_2_object();
      seg = construct_seg(t, s);
    }

    // force the vector to [0-90)
    Vector_2 v(construct_vertex(seg, 0), construct_vertex(seg, 1));
    if (cx == EQUAL || (cx == LARGER && cy == SMALLER) ||
        (cx == SMALLER && cy == LARGER))
      v = v.perpendicular(RIGHT_TURN);

    Direction_2 d(v.direction());

    return(d);
  }

  /*! */
  int get_kd_num(Segment_2 seg, int n, Direction_list & directions)
  {
    Direction_2 d = get_direction(seg);
    int i = 0;
    bool found = false;
    Direction_const_iter iter = directions.begin();

    while (i < n && !found) {
      if (*iter > d) found = true;
      ++i;
      ++iter;
    }

    if (found) --i;
    else i = 0;

    return(i);
  }

  /*! */
  void check_kd(int *kd_counter,int number_of_trees,
                const Segment_list & seg_list,
                Direction_list & directions)
  {
    for (int i = 0;i < number_of_trees;++i)
      kd_counter[i] = 0;

    int kd_num;
    for (Segment_const_iter iter = seg_list.begin(); iter != seg_list.end();
         ++iter)
    {
      kd_num = get_kd_num(*iter, number_of_trees, directions);
      kd_counter[kd_num]++;
    }
  }

  /*! */
  void init_angle_to_sines_table()
  {
    angle_to_sines_appr[0] = NT(0);
    angle_to_sines_appr[1] = NT(115) / NT(6613);
    angle_to_sines_appr[2] = NT(57) / NT(1625);
    angle_to_sines_appr[3] = NT(39) / NT(761);
    angle_to_sines_appr[4] = NT(29) / NT(421);
    angle_to_sines_appr[5] = NT(23) / NT(265);
    angle_to_sines_appr[6] = NT(19) / NT(181);
    angle_to_sines_appr[7] = NT(32) / NT(257);
    angle_to_sines_appr[8] = NT(129) / NT(929);
    angle_to_sines_appr[9] = NT(100) / NT(629);
    angle_to_sines_appr[10] = NT(92) / NT(533);
    angle_to_sines_appr[11] = NT(93) / NT(485);
    angle_to_sines_appr[12] = NT(76) / NT(365);
    angle_to_sines_appr[13] = NT(156) / NT(685);
    angle_to_sines_appr[14] = NT(205) / NT(853);
    angle_to_sines_appr[15] = NT(69) / NT(269);
    angle_to_sines_appr[16] = NT(7) / NT(25);
    angle_to_sines_appr[17] = NT(120) / NT(409);
    angle_to_sines_appr[18] = NT(57) / NT(185);
    angle_to_sines_appr[19] = NT(12) / NT(37);
    angle_to_sines_appr[20] = NT(51) / NT(149);
    angle_to_sines_appr[21] = NT(135) / NT(377);
    angle_to_sines_appr[22] = NT(372) / NT(997);
    angle_to_sines_appr[23] = NT(348) / NT(877);
    angle_to_sines_appr[24] = NT(231) / NT(569);
    angle_to_sines_appr[25] = NT(36) / NT(85);
    angle_to_sines_appr[26] = NT(39) / NT(89);
    angle_to_sines_appr[27] = NT(300) / NT(661);
    angle_to_sines_appr[28] = NT(8) / NT(17);
    angle_to_sines_appr[29] = NT(189) / NT(389);
    angle_to_sines_appr[30] = NT(451) / NT(901);
    angle_to_sines_appr[31] = NT(180) / NT(349);
    angle_to_sines_appr[32] = NT(28) / NT(53);
    angle_to_sines_appr[33] = NT(432) / NT(793);
    angle_to_sines_appr[34] = NT(161) / NT(289);
    angle_to_sines_appr[35] = NT(228) / NT(397);
    angle_to_sines_appr[36] = NT(504) / NT(865);
    angle_to_sines_appr[37] = NT(3) / NT(5);
    angle_to_sines_appr[38] = NT(580) / NT(941);
    angle_to_sines_appr[39] = NT(341) / NT(541);
    angle_to_sines_appr[40] = NT(88) / NT(137);
    angle_to_sines_appr[41] = NT(48) / NT(73);
    angle_to_sines_appr[42] = NT(65) / NT(97);
    angle_to_sines_appr[43] = NT(429) / NT(629);
    angle_to_sines_appr[44] = NT(555) / NT(797);
    angle_to_sines_appr[45] = NT(697) / NT(985);
    angle_to_sines_appr[46] = NT(572) / NT(797);
    angle_to_sines_appr[47] = NT(460) / NT(629);
    angle_to_sines_appr[48] = NT(72) / NT(97);
    angle_to_sines_appr[49] = NT(55) / NT(73);
    angle_to_sines_appr[50] = NT(105) / NT(137);
    angle_to_sines_appr[51] = NT(420) / NT(541);
    angle_to_sines_appr[52] = NT(741) / NT(941);
    angle_to_sines_appr[53] = NT(4) / NT(5);
    angle_to_sines_appr[54] = NT(703) / NT(865);
    angle_to_sines_appr[55] = NT(325) / NT(397);
    angle_to_sines_appr[56] = NT(240) / NT(289);
    angle_to_sines_appr[57] = NT(665) / NT(793);
    angle_to_sines_appr[58] = NT(45) / NT(53);
    angle_to_sines_appr[59] = NT(299) / NT(349);
    angle_to_sines_appr[60] = NT(780) / NT(901);
    angle_to_sines_appr[61] = NT(340) / NT(389);
    angle_to_sines_appr[62] = NT(15) / NT(17);
    angle_to_sines_appr[63] = NT(589) / NT(661);
    angle_to_sines_appr[64] = NT(80) / NT(89);
    angle_to_sines_appr[65] = NT(77) / NT(85);
    angle_to_sines_appr[66] = NT(520) / NT(569);
    angle_to_sines_appr[67] = NT(805) / NT(877);
    angle_to_sines_appr[68] = NT(925) / NT(997);
    angle_to_sines_appr[69] = NT(352) / NT(377);
    angle_to_sines_appr[70] = NT(140) / NT(149);
    angle_to_sines_appr[71] = NT(35) / NT(37);
    angle_to_sines_appr[72] = NT(176) / NT(185);
    angle_to_sines_appr[73] = NT(391) / NT(409);
    angle_to_sines_appr[74] = NT(24) / NT(25);
    angle_to_sines_appr[75] = NT(260) / NT(269);
    angle_to_sines_appr[76] = NT(828) / NT(853);
    angle_to_sines_appr[77] = NT(667) / NT(685);
    angle_to_sines_appr[78] = NT(357) / NT(365);
    angle_to_sines_appr[79] = NT(476) / NT(485);
    angle_to_sines_appr[80] = NT(525) / NT(533);
    angle_to_sines_appr[81] = NT(621) / NT(629);
    angle_to_sines_appr[82] = NT(920) / NT(929);
    angle_to_sines_appr[83] = NT(255) / NT(257);
    angle_to_sines_appr[84] = NT(180) / NT(181);
    angle_to_sines_appr[85] = NT(264) / NT(265);
    angle_to_sines_appr[86] = NT(420) / NT(421);
    angle_to_sines_appr[87] = NT(760) / NT(761);
    angle_to_sines_appr[88] = NT(1624) / NT(1625);
    angle_to_sines_appr[89] = NT(6612) / NT(6613);
    angle_to_sines_appr[90] = NT(1);
  }

public:

  /*! */
  Multiple_kd_tree(const Point_saved_pair_list & inp_points_list,
                   int inp_number_of_trees,
                   const Segment_list & seg_list) :
    pi(3.1415), half_pi(1.57075),
    number_of_trees(inp_number_of_trees), input_points_list(inp_points_list)
  {

    Kd_triple kd;

    // check that there are at least two trees
    CGAL_precondition_msg(number_of_trees >= 1, "There must be at least one kd-tree" );

    init_angle_to_sines_table();

    // find the kd trees that have enough candidates  (segments with a close
    // slope)
    int * kd_counter = new int[number_of_trees];
    int number_of_segments = seg_list.size();

    // auxilary directions
    Direction_list directions;
    double buffer_angle;
    Line_2 li;
    Direction_2 d;

    int i = 0;
    for (double angle = 0; i < number_of_trees;
         angle += half_pi / number_of_trees,++i)
    {
      buffer_angle = angle - half_pi / (2 * number_of_trees);

      if (buffer_angle < 0)
        buffer_angle = 0;

      li = Line_2(std::tan(buffer_angle), -1, 0);
      d = Direction_2(li);

      // rotate_by 180 degrees
      Transformation_2 t(ROTATION, 0, -1);
      d = d.transform(t);

      directions.push_back(d);
    }

    check_kd(kd_counter, number_of_trees, seg_list,directions);
    int ind = 0;

#ifdef CGAL_SR_DEBUG
    int number_of_actual_kd_trees = 0;
#endif

    i = 0;

    for (NT angle = 0; i < number_of_trees;
         angle += NT(half_pi / number_of_trees),++i)
    {
      if (kd_counter[ind] >=
          (double)number_of_segments / (double)number_of_trees / 2.0)
      {
        kd = create_kd_tree(angle);

        kd_trees_list.push_back(kd);

#ifdef CGAL_SR_DEBUG
        ++number_of_actual_kd_trees;
#endif

      }

      ++ind;
    }

    delete[] kd_counter;

#ifdef CGAL_SR_DEBUG
    std::cout << "Actual number of kd-trees created : " <<
      number_of_actual_kd_trees << std::endl;
#endif

  }

  ~Multiple_kd_tree()
  {
    //delete all the kd_trees.
    for(typename Kd_triple_list::iterator it = kd_trees_list.begin();   it != kd_trees_list.end(); ++it)
      delete (it->first);

    //delete all the points.
    for(typename Point_saved_pair_list::iterator it = input_points_list.begin();
        it != input_points_list.end(); ++it) {
      delete (it->second);
    }

  }

  /*! */
  Point_2 small_x_point(const Point_2 & p1, const Point_2 & p2)
  {
    Comparison_result c = m_gt.compare_x_2_object()(p1, p2);
    return (c == SMALLER) ? p1 : p2;
  }

  /*! */
  Point_2 big_x_point(const Point_2 & p1, const Point_2 & p2)
  {
    Comparison_result c = m_gt.compare_x_2_object()(p1, p2);
    return (c == SMALLER) ? p2 : p1;
  }

  /*! */
  Point_2 small_y_point(const Point_2 & p1, const Point_2 & p2)
  {
    Comparison_result c = m_gt.compare_y_2_object()(p1, p2);
    return (c == SMALLER) ? p1 : p2;
  }

  /*! */
  Point_2 big_y_point(const Point_2 & p1, const Point_2 & p2)
  {
    Comparison_result c = m_gt.compare_y_2_object()(p1, p2);
    return (c == SMALLER) ? p2 : p1;
  }

  /*! */
  void get_intersecting_points(std::list<SAVED_OBJECT> & result_list,
                               Segment_2 s, NT unit_square)
  {
    // determine right kd-tree to work on, depending on the segment's slope
    Direction_2 d = get_direction(s);

    int i = 0;
    int n = kd_trees_list.size();
    bool found = false;
    typename Kd_triple_list::const_iterator iter = kd_trees_list.begin();

    while(i < n && !found)
    {
      if (iter->second.first > d)
        found = true;

      ++i;
      ++iter;
    }

    if (!found)
      iter = kd_trees_list.begin();

    else
      --iter;

    Point_list points_list;
    m_gt.minkowski_sum_with_pixel_2_object()(points_list, s, unit_square);

    Point_iter points_iter;

    for (points_iter = points_list.begin();   points_iter != points_list.end();   ++points_iter)
      rotate(*points_iter, iter->second.second);

    // query
    points_iter = points_list.begin();
    Point_2 point_left, point_right, point_bot, point_top;
    point_left = point_right = point_bot = point_top = *points_iter;

    for (++points_iter; points_iter != points_list.end(); ++points_iter)
    {
      point_left = small_x_point(point_left,*points_iter);
      point_right = big_x_point(point_right,*points_iter);
      point_bot = small_y_point(point_bot,*points_iter);
      point_top = big_y_point(point_top,*points_iter);
    }

    typedef typename Traits::Construct_iso_rectangle_2    Construct_iso_rectangle_2;

    Construct_iso_rectangle_2 construct_rec = m_gt.construct_iso_rectangle_2_object();

    Iso_rectangle_2 rec =  construct_rec(point_left, point_right, point_bot, point_top);

    Point_2 p1 = rec.vertex(0);
    Point_2 p2 = rec.vertex(2);

    Point_with_hot_pixel_history_saved point1(p1);
    Point_with_hot_pixel_history_saved point2(p2);

    Box b(point1, point2);

    // the kd-tree query
    Point_with_hot_pixel_history_saved_list result;

    iter->first->search(std::back_inserter(result), b);

    // create result
    result_list.empty();

    for( Point_with_hot_pixel_history_saved_iter my_point_iter = result.begin();    my_point_iter != result.end();   ++my_point_iter )
      result_list.push_back(my_point_iter->object);
  }
};

} //namespace CGAL

#endif
