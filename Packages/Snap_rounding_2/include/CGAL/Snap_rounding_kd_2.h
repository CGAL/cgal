// Copyright (c) 1997  Tel-Aviv University (Israel).
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
#ifndef CGAL_SR_KD_2_H
#define CGAL_SR_KD_2_H

#include <list>
#include <CGAL/config.h>
#include <CGAL/kdtree_d.h>
#include <CGAL/predicates_on_points_2.h>
#include <iostream>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/utility.h>

CGAL_BEGIN_NAMESPACE

template<class Rep,class SAVED_OBJECT>
class my_point : public Rep::Point_2 {

typedef typename Rep::Point_2               Point_2;
typedef typename Rep::FT                    NT;

public:
  Point_2 orig;
  SAVED_OBJECT object;
  my_point(Point_2 p,Point_2 inp_orig,SAVED_OBJECT obj) :
           Point_2(p),orig(inp_orig),object(obj) {}
  my_point(Point_2 p) : Point_2(p),orig(Point_2(0,0)) {}
  my_point() : Point_2(),orig() {}
  my_point(NT x,NT y) : Point_2(x,y),orig(Point_2(0,0)) {}
};

template<class Rep_,class SAVED_OBJECT>
class Multiple_kd_tree {

typedef Rep_                                          Rep;
typedef typename Rep::FT                              NT;
typedef typename Rep::Segment_2                       Segment_2;
typedef typename Rep::Point_2                         Point_2;
typedef typename Rep::Vector_2                        Vector_2;
typedef typename Rep::Iso_rectangle_2                 Iso_rectangle_2;
typedef typename Rep::Direction_2                     Direction_2;
typedef typename Rep::Line_2                          Line_2;
typedef typename Rep::Aff_transformation_2            Transformation_2;
typedef CGAL::Kdtree_interface_2d<my_point<Rep,SAVED_OBJECT> >
                                                      kd_interface;
typedef CGAL::Kdtree_d<kd_interface>                  kd_tree;
typedef typename kd_tree::Box                         Box;
typedef std::list<my_point<Rep,SAVED_OBJECT> >        Points_List; 
typedef std::pair<kd_tree *,std::pair<Direction_2,NT> >
                                                      kd_triple;
typedef std::list<std::pair<kd_tree *,std::pair<Direction_2,NT> > >
                                                      kd_triple_list;
private:
  Rep_   _gt;
  const double pi,half_pi;
  int number_of_trees;
  kd_triple_list kd_trees_list;
  std::list<std::pair<Point_2,SAVED_OBJECT > > input_points_list;
  std::map<int,typename Rep::FT> angle_to_sines_appr; // was const int

  void rotate(Point_2& p,NT angle)
  {
    static const double rad_to_deg = 57.297;
    int tranc_angle = int(to_double(angle) * rad_to_deg);

    NT cosine_val = angle_to_sines_appr[90 - tranc_angle],
       sine_val = angle_to_sines_appr[tranc_angle];

    Transformation_2 rotate(ROTATION, sine_val, cosine_val);
    p = rotate(p);
  }


  kd_triple create_kd_tree(NT angle)
  {
    Points_List l;
    kd_tree *tree = new kd_tree(2);

    for(typename std::list<std::pair<Point_2,SAVED_OBJECT> >::iterator
        iter = input_points_list.begin();
        iter != input_points_list.end();
        ++iter) {
      Point_2 p(iter->first);
      rotate(p,angle);
      my_point<Rep,SAVED_OBJECT> rotated_point(p,iter->first,iter->second);
      l.push_back(rotated_point);
    }

    tree->build(l);

    //checking validity
    if(!tree->is_valid() )
      tree->dump();
    assert(tree->is_valid());

    double buffer_angle(to_double(angle) - half_pi / (2 * number_of_trees));
    if(buffer_angle < 0)
      buffer_angle = 0;
    Line_2 li(tan(buffer_angle),-1,0);
    Direction_2 d(li);
    // rotate_by 180 degrees
    Transformation_2 t(ROTATION,0,-1);
    d = d.transform(t);
    std::pair<Direction_2,NT> kp(d,angle);
    kd_triple kt(tree,kp);

    return(kt);
  }

  inline NT square(NT x) {return(x * x);}
  inline NT min(NT x,NT y) {return((x < y) ? x : y);}
  inline NT max(NT x,NT y) {return((x < y) ? y : x);}
  inline NT min(NT x1,NT x2,NT x3,NT x4,NT x5,NT x6) 
       {return(min(min(min(x1,x2),
                min(x3,x4)),min(x5,x6)));}
  inline NT max(NT x1,NT x2,NT x3,NT x4,NT x5,NT x6) 
       {return(max(max(max(x1,x2),
                max(x3,x4)),max(x5,x6)));}

  Direction_2 get_direction(Segment_2 seg)
  {
    // force the segment slope to [0-180)
    Point_2 s = seg.source(),t = seg.target();
    Comparison_result cx = _gt.compare_x_2_object()(s,t);
    Comparison_result cy = _gt.compare_y_2_object()(s,t);
    if(cy == LARGER || cy == EQUAL && cx == LARGER)
      seg = Segment_2(t,s);

    // force the vector to [0-90)
    Vector_2 v(seg.source(),seg.target());
    if(cx == EQUAL || cx == LARGER && cy == SMALLER ||
       cx == SMALLER && cy == LARGER)
      v = v.perpendicular(RIGHT_TURN);

    Direction_2 d(v.direction());
 
    return(d);
  }

  int get_kd_num(Segment_2 seg,int n,std::list<Direction_2>& directions)
  {
    Direction_2 d = get_direction(seg);
    int i = 0;
    bool found = false;
    typename std::list<Direction_2>::const_iterator
      iter = directions.begin();

    while(i < n && !found) {
      if(*iter > d)
        found = true;
      ++i;
      ++iter;
    }

    if(found)
      --i;
    else
      i = 0;

    return(i);
  }

  void check_kd(int *kd_counter,int number_of_trees,
       const std::list<Segment_2> &seg_list,std::list<Direction_2>& directions)
  {
    for(int i = 0;i < number_of_trees;++i)
      kd_counter[i] = 0;

    int kd_num;
    for(typename std::list<Segment_2>::const_iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      kd_num = get_kd_num(*iter,number_of_trees,directions);
      kd_counter[kd_num]++;
    }
  }

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

  Multiple_kd_tree(const std::list<std::pair<Point_2,SAVED_OBJECT> > 
                   &inp_points_list,int inp_number_of_trees,
                   const std::list<Segment_2> &seg_list) : 
    pi(3.1415),half_pi(1.57075),
    number_of_trees(inp_number_of_trees),input_points_list(inp_points_list)
  {
    kd_triple kd;

    // check that there are at least two trees
    if(number_of_trees < 1) {
      std::cerr << "There must be at least one kd-tree\n";
      exit(1);
    }

    init_angle_to_sines_table();

    // find the kd trees that have enough candidates  (segments with a close
    // slope)
    int *kd_counter = new int[number_of_trees];
    int number_of_segments = seg_list.size();

    // auxilary directions
    std::list<Direction_2> directions;
    double buffer_angle;
    Line_2 li;
    Direction_2 d;

    int i = 0;
    for(double angle = 0;
        i < number_of_trees;
        angle += half_pi / number_of_trees,++i) {
      buffer_angle = angle - half_pi / (2 * number_of_trees);
      if(buffer_angle < 0)
        buffer_angle = 0;
      li = Line_2(tan(buffer_angle),-1,0);
      d = Direction_2(li);
      // rotate_by 180 degrees
      Transformation_2 t(ROTATION,0,-1);
      d = d.transform(t);

      directions.push_back(d);
    }

    check_kd(kd_counter,number_of_trees,seg_list,directions);
    int ind = 0;

#ifdef KD_DEBUG
    int number_of_actual_kd_trees = 0;
#endif
    i = 0;
    for(NT angle = 0;
        i < number_of_trees;
        angle += NT(half_pi / number_of_trees),++i) {
      if(kd_counter[ind] >= (double)number_of_segments /
	                    (double)number_of_trees / 2.0) {
        kd = create_kd_tree(angle);
        kd_trees_list.push_back(kd);

#ifdef KD_DEBUG
        ++number_of_actual_kd_trees;
#endif

      }

      ++ind;
    }

#ifdef KD_DEBUG
    std::cout << "Actual number of kd-trees created : " << 
      number_of_actual_kd_trees << std::endl;
#endif

  }

  Point_2 small_x_point(const Point_2& p1,const Point_2& p2)
  {
    Comparison_result c = _gt.compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_x_point(const Point_2& p1,const Point_2& p2)
  {
    Comparison_result c = _gt.compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  Point_2 small_y_point(const Point_2& p1,const Point_2& p2)
  {
    Comparison_result c = _gt.compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_y_point(const Point_2& p1,const Point_2& p2)
  {
    Comparison_result c = _gt.compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  void get_intersecting_points(std::list<SAVED_OBJECT> &result_list,
                               Segment_2 s,
                               NT unit_square)
  {
    // determine right kd-tree to work on, depending on the segment's slope
    Direction_2 d = get_direction(s);
    int i = 0;
    int n = kd_trees_list.size();
    bool found = false;
    typename kd_triple_list::const_iterator
      iter = kd_trees_list.begin();

    while(i < n && !found) {
      if(iter->second.first > d)
        found = true;
      ++i;
      ++iter;
    }

    if(!found)
      iter = kd_trees_list.begin();
    else
      --iter;

    std::list<Point_2> points_list;
    _gt.minkowski_sum_with_pixel_2_object()
      (points_list,s,unit_square);

    typename std::list<Point_2>::iterator points_iter;

    for(points_iter = points_list.begin();
        points_iter != points_list.end();
        ++points_iter)
      rotate(*points_iter,iter->second.second);

    // query
    points_iter = points_list.begin();
    Point_2 point_left,point_right,point_bot,point_top;
    point_left = point_right = point_bot = point_top = *points_iter;
    for(++points_iter;
        points_iter != points_list.end();
        ++points_iter) {
      point_left = small_x_point(point_left,*points_iter);
      point_right = big_x_point(point_right,*points_iter);
      point_bot = small_y_point(point_bot,*points_iter);
      point_top = big_y_point(point_top,*points_iter);
    }

    Iso_rectangle_2 rec(point_left,point_right,point_bot,point_top);

    Point_2 p1 = rec.vertex(0);
    Point_2 p2 = rec.vertex(2);

    my_point<Rep,SAVED_OBJECT> point1(p1); 
    my_point<Rep,SAVED_OBJECT> point2(p2);

    Box b(point1,point2,2);
 
    // the kd-tree query
    std::list<my_point<Rep,SAVED_OBJECT> > res;
    iter->first->search(std::back_inserter(res),b);

    // create result
    result_list.empty();
    for(typename std::list<my_point<Rep,SAVED_OBJECT> >::iterator 
        my_point_iter = res.begin();
        my_point_iter != res.end();
        ++my_point_iter)
      result_list.push_back(my_point_iter->object);
  }
};

CGAL_END_NAMESPACE

#endif
