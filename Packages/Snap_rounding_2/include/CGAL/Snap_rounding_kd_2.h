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
// file          : include/CGAL/Isr_kd_2.h
// package       : arr (1.73)
// maintainer    : Eli Packer <elip@post.tau.ac.il>
// author(s)     : Eli Packer
// coordinator   : ana Halperin Tel-Aviv University
//                 (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_SR_KD_2_H
#define CGAL_SR_KD_2_H

#include <list>
#include <CGAL/config.h>
#include <CGAL/kdtree_d.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/leda_rational.h> 
#include <iostream>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/utility.h>
#include "../../include/CGAL/Snap_rounding_2_utility.h"

CGAL_BEGIN_NAMESPACE

template<class NT,class SAVED_OBJECT>
class my_point : public CGAL::Point_2<CGAL::Cartesian<NT> > {

typedef CGAL::Point_2<CGAL::Cartesian<NT> >  Point;

public:
  Point orig;
  SAVED_OBJECT object;
  my_point(Point p,Point inp_orig,SAVED_OBJECT obj) : Point(p),orig(inp_orig),
           object(obj) {}
  my_point(Point p) : Point(p),orig(Point(0,0)) {}
  my_point() : Point(),orig() {}
  my_point(NT x,NT y) : Point(x,y),orig(Point(0,0)) {}
};

template<class Rep_,class SAVED_OBJECT>
class Multiple_kd_tree {

typedef Rep_                                Rep;
typedef typename Rep::FT                    NT;
typedef typename Rep::Segment_2             Segment_2;
typedef typename Rep::Point_2               Point_2;
typedef typename Rep::Vector_2              Vector_2;
typedef typename Rep::Iso_rectangle_2       Iso_rectangle_2;
typedef typename Rep::Direction_2           Direction_2;
typedef typename Rep::Line_2                Line_2;
typedef CGAL::Kdtree_interface_2d<my_point<NT,SAVED_OBJECT> >  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef typename kd_tree::Box Box;
typedef std::list<my_point<NT,SAVED_OBJECT> > Points_List; 
typedef std::pair<kd_tree *,std::pair<Direction_2,NT> > kd_triple;
typedef std::list<std::pair<kd_tree *,std::pair<Direction_2,NT> > >
  kd_triple_list;


private:
  Rep_   _gt;
  const double pi,half_pi,epsilon;
  int number_of_trees;
  kd_triple_list kd_trees_list;
  std::list<std::pair<Point_2,SAVED_OBJECT > > input_points_list;

  kd_triple create_kd_tree(NT angle)
  {
    Points_List l;
    kd_tree *tree = new kd_tree(2);

    for(typename std::list<std::pair<Point_2,SAVED_OBJECT> >::iterator
        iter = input_points_list.begin();
        iter != input_points_list.end();
        ++iter) {
      Point_2 p(iter->first);

      static Snap_rounding_rotation<Rep> r;
      r(p,angle);

      my_point<NT,SAVED_OBJECT> rotated_point(p,iter->first,iter->second);

      l.push_back(rotated_point);
    }

    tree->build(l);

    //checking validity
    if(!tree->is_valid() )
      tree->dump();
    assert(tree->is_valid());

    NT buffer_angle(angle - half_pi / (2 * number_of_trees));
    Line_2 li(tan(buffer_angle.to_double()),-1,0);
    Direction_2 d(li);
    std::pair<Direction_2,NT> kp(d,angle);
    kd_triple kt(tree,kp);

    return(kt);
  }

  inline NT squere(NT x) {return(x * x);}
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
      v = v.perpendicular(RIGHTTURN);

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
       std::list<Segment_2> &seg_list,std::list<Direction_2>& directions)
  {
    for(int i = 0;i < number_of_trees;++i)
      kd_counter[i] = 0;

    int kd_num;
    for(typename std::list<Segment_2>::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      kd_num = get_kd_num(*iter,number_of_trees,directions);
      kd_counter[kd_num]++;
    }
  }

public:

  Multiple_kd_tree(std::list<std::pair<Point_2,SAVED_OBJECT> > 
                   &inp_points_list,int inp_number_of_trees,
                   std::list<Segment_2> &seg_list) : 
    pi(3.1415),half_pi(1.57075),epsilon(0.001),
    number_of_trees(inp_number_of_trees),input_points_list(inp_points_list)
  {
    kd_triple kd;

    // check that there are at least two trees
    if(number_of_trees < 1) {
      std::cerr << "There must be at least one kd-tree\n";
      exit(1);
    }

    // find the kd trees that have enough candidates  (segments with a close
    // slope)
    int *kd_counter = new int[number_of_trees];
    int number_of_segments = seg_list.size();

    // auxilary directions
    std::list<Direction_2> directions;
    NT buffer_angle;
    Line_2 li;
    Direction_2 d;
    for(NT angle = 0;angle < half_pi;angle += half_pi / number_of_trees) {
      buffer_angle = angle - half_pi / (2 * number_of_trees);
      li = Line_2(tan(buffer_angle.to_double()),-1,0);
      d = Direction_2(li);
      directions.push_back(d);
    }

    check_kd(kd_counter,number_of_trees,seg_list,directions);
    int ind = 0;
    for(NT angle = 0;angle < half_pi;angle += half_pi / number_of_trees) {
      if(kd_counter[ind] >= (double)number_of_segments /
	                    (double)number_of_trees / 2.0) {
        kd = create_kd_tree(angle);
        kd_trees_list.push_back(kd);
      }

      ++ind;
    }
  }

  Point_2 small_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt.compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt.compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  Point_2 small_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt.compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt.compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  void get_intersecting_points(list<SAVED_OBJECT> &result_list,
                               Segment_2 s,
                               NT unit_squere)
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
    _gt.bounding_box_of_minkowski_sum_2_object()
      (points_list,s,unit_squere);

    static Snap_rounding_rotation<Rep> r;

    typename std::list<Point_2>::iterator points_iter;

    for(points_iter = points_list.begin();
        points_iter != points_list.end();
        ++points_iter)
      r(*points_iter,iter->second.second);

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

    my_point<NT,SAVED_OBJECT> point1(p1); 
    my_point<NT,SAVED_OBJECT> point2(p2);

    Box b(point1,point2,2);
 
    // the kd-tree query
    list<my_point<NT,SAVED_OBJECT> > res;
    iter->first->search(std::back_inserter(res),b);

    // create result
    result_list.empty();
    for(typename list<my_point<NT,SAVED_OBJECT> >::iterator 
        my_point_iter = res.begin();
        my_point_iter != res.end();
        ++my_point_iter)
      result_list.push_back(my_point_iter->object);
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_SR_KD_2_H
