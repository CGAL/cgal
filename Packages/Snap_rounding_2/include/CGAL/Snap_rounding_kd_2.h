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
// coordinator   : Dna Halperin Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
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
typedef CGAL::Segment_2<Rep>                Segment;// !!!! look at LER
typedef CGAL::Point_2<Rep>                  Point;
typedef CGAL::Iso_rectangle_2<Rep>          Iso_rectangle_2;
typedef CGAL::Kdtree_interface_2d<my_point<NT,SAVED_OBJECT> >  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef typename kd_tree::Box Box;
typedef std::list<my_point<NT,SAVED_OBJECT> > Points_List; 

private:
  Rep_   _gt;
  const double pi,half_pi,epsilon;
  int number_of_trees;
  std::list<std::pair<kd_tree *,NT> > kd_trees_list;
  std::list<std::pair<Point,SAVED_OBJECT > > input_points_list;

  std::pair<kd_tree *,NT> create_kd_tree(NT angle)
  {
    Points_List l;
    kd_tree *tree = new kd_tree(2);

    for(typename std::list<std::pair<Point,SAVED_OBJECT> >::iterator
        iter = input_points_list.begin();
        iter != input_points_list.end();
        ++iter) {

      Point p(iter->first);

      _gt.rotate_point(p,angle);

      my_point<NT,SAVED_OBJECT> rotated_point(p,iter->first,iter->second);

      l.push_back(rotated_point);
    }

    tree->build(l);

    //checking validity
    if(!tree->is_valid() )
      tree->dump();
    assert(tree->is_valid());

    return(std::pair<kd_tree *,NT>(tree,angle));
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

  int get_kd_num(Segment seg,int n)
  {
    double alpha = _gt.segment_direction(seg);

    int i;

    if(alpha < 0)
      alpha += pi / 2.0;

    if(alpha >= pi * (2 * n - 1) / (4 * n))
      i = 0;
    else {
      alpha += pi / (4 * n);
      i = int(2 * n * alpha / pi);
    }

    return(i);
  }

  void check_kd(int *kd_counter,int number_of_trees,
       std::list<Segment> &seg_list)
  {
    for(int i = 0;i < number_of_trees;++i)
      kd_counter[i] = 0;

    int kd_num;
    for(typename std::list<Segment>::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      kd_num = get_kd_num(*iter,number_of_trees);
      kd_counter[kd_num]++;
    }
  }

public:
  Multiple_kd_tree(std::list<std::pair<Point,SAVED_OBJECT> > 
                   &inp_points_list,int inp_number_of_trees,
                   std::list<Segment> &seg_list) : 
    pi(3.1415),half_pi(1.57075),epsilon(0.001),
    number_of_trees(inp_number_of_trees),input_points_list(inp_points_list)
  {
    std::pair<kd_tree *,NT> kd;

    // check that there are at least two trees
    if(number_of_trees < 1) {
      std::cerr << "There must be at least one kd-tree\n";
      exit(1);
    }

    // find the kd trees that have enough candidates  (segments with a close
    // slope
    int *kd_counter = new int[number_of_trees];
    int number_of_segments = seg_list.size();
    check_kd(kd_counter,number_of_trees,seg_list);
    int ind = 0;
    for(NT angle = 0;angle < half_pi;angle += half_pi / number_of_trees) {
      if(kd_counter[ind] >= (double)number_of_segments /
	                    (double) number_of_trees / 2.0) {
        kd = create_kd_tree(angle);
        kd_trees_list.push_back(kd);
      }

      ++ind;
    }
  }

  Point small_x_point(Point p1,Point p2)
  {
     Comparison_result c = _gt.compare_x_2_object()(p1,p2);
     if(c == SMALLER)
       return(p1);
     else
       return(p2);
  }

  Point small_y_point(Point p1,Point p2)
  {
     Comparison_result c = _gt.compare_y_2_object()(p1,p2);
     if(c == SMALLER)
       return(p1);
     else
       return(p2);
  }

  Point big_x_point(Point p1,Point p2)
  {
     Comparison_result c = _gt.compare_x_2_object()(p1,p2);
     if(c == SMALLER)
       return(p2);
     else
       return(p1);
  }

  Point big_y_point(Point p1,Point p2)
  {
     Comparison_result c = _gt.compare_y_2_object()(p1,p2);
     if(c == SMALLER)
       return(p2);
     else
       return(p1);
  }


  void get_intersecting_points(list<SAVED_OBJECT> &result_list,
                               Segment inp_s,
                               NT unit_squere)
  {
    Comparison_result cy = _gt.compare_y_2_object()(inp_s.source(),inp_s.target());
    Segment s(cy == SMALLER ?
              inp_s.source() : inp_s.target(),
              cy == SMALLER ?
              inp_s.target() : inp_s.source());

    // determine right kd-tree to work on, depending on the segment's slope
    double alpha_double = _gt.segment_direction(s);

    if(alpha_double < 0)
      alpha_double += pi / 2.0;

    NT alpha = alpha_double;

    bool found = false;
    NT last_dif;

    typename list<std::pair<kd_tree *,NT> >::iterator iter,right_iter;

    for(iter = kd_trees_list.begin();
        iter != kd_trees_list.end() && !found;
        ++iter) {
      if(iter->second > alpha) {
        right_iter = iter;
        if(iter != kd_trees_list.begin()) {
	  if(iter->second - alpha > last_dif)
            --right_iter;
	}
 
        found = true;
      } else
        last_dif = iter->second - alpha;
    }
    if(!found) {
      right_iter = kd_trees_list.end();
      --right_iter;
    }

    // query
    Point ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points
    list<my_point<NT,SAVED_OBJECT> > res;

    Comparison_result cx = _gt.compare_x_2_object()(s.source(),s.target());
    NT x1 = s.source().x(),y1 = s.source().y(),x2 =
       s.target().x(),y2 = s.target().y();
    if(cx == SMALLER) {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point(x1 - 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms3 = Point(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms4 = Point(x2 + 0.6 * unit_squere,y2 - 0.6 * unit_squere);
      ms5 = Point(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
    } else {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms3 = Point(x1 + 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms4 = Point(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms5 = Point(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point(x2 - 0.6 * unit_squere,y2 - 0.6 * unit_squere);
    }

    _gt.rotate_point(ms1,right_iter->second);
    _gt.rotate_point(ms2,right_iter->second);
    _gt.rotate_point(ms3,right_iter->second);
    _gt.rotate_point(ms4,right_iter->second);
    _gt.rotate_point(ms5,right_iter->second);
    _gt.rotate_point(ms6,right_iter->second);

    // query
    Point point_left,point_right,point_bot,point_top;

    point_left = small_x_point(ms1,ms2);
    point_left = small_x_point(point_left,ms3);
    point_left = small_x_point(point_left,ms4);
    point_left = small_x_point(point_left,ms5);
    point_left = small_x_point(point_left,ms6);

    point_right = big_x_point(ms1,ms2);
    point_right = big_x_point(point_right,ms3);
    point_right = big_x_point(point_right,ms4);
    point_right = big_x_point(point_right,ms5);
    point_right = big_x_point(point_right,ms6);

    point_bot = small_y_point(ms1,ms2);
    point_bot = small_y_point(point_bot,ms3);
    point_bot = small_y_point(point_bot,ms4);
    point_bot = small_y_point(point_bot,ms5);
    point_bot = small_y_point(point_bot,ms6);

    point_top = big_y_point(ms1,ms2);
    point_top = big_y_point(point_top,ms3);
    point_top = big_y_point(point_top,ms4);
    point_top = big_y_point(point_top,ms5);
    point_top = big_y_point(point_top,ms6);

    Iso_rectangle_2 rec(point_left,point_right,point_bot,point_top);

    Point p1 = rec.vertex(0);
    Point p2 = rec.vertex(2);// end of new code

    my_point<NT,SAVED_OBJECT> point1(p1); 
    my_point<NT,SAVED_OBJECT> point2(p2);

    Box b(point1,point2,2);
 
    // the kd-tree query
    right_iter->first->search(std::back_inserter(res),b);

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
