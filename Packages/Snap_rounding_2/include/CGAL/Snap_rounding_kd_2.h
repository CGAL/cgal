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
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
#ifndef CGAL_SR_KD_2_H
#define CGAL_SR_KD_2_H

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif

#include <list>

#ifndef  CGAL_KDTREE_D_H
#include <CGAL/kdtree_d.h>
#endif

#ifndef CGAL_AFF_TRANSFORMATION_2_H
#include <CGAL/Aff_transformation_2.h>
#endif

#ifndef CGAL_VECTOR_2_H
#include <CGAL/Vector_2.h>
#endif

#define PI 3.1415
#define HALF_PI NT(1.57075)
// Since kd does not output intersecting points with the BBOX we add EPSILON
#define EPSILON 0.5

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

template<class NT,class SAVED_OBJECT>
class Multiple_kd_tree {

typedef CGAL::Point_2<CGAL::Cartesian<NT> >  Point;
typedef CGAL::Segment_2<CGAL::Cartesian<NT> >      Segment;
typedef CGAL::Kdtree_interface_2d<my_point<NT,SAVED_OBJECT> >  kd_interface;
typedef CGAL::Kdtree_d<kd_interface>  kd_tree;
typedef typename kd_tree::Box Box;
typedef std::list<my_point<NT,SAVED_OBJECT> > Points_List; 
typedef CGAL::Aff_transformation_2<CGAL::Cartesian<NT> > Transformation;

private:
  int number_of_trees;
  std::list<std::pair<kd_tree *,NT> > kd_trees_list;
  std::list<std::pair<std::pair<NT,NT>,SAVED_OBJECT > > input_points_list;

  kd_tree *create_kd_tree(NT angle)
  {
    Points_List l;
    kd_tree *tree = new kd_tree(2);

    Transformation rotate(CGAL::ROTATION, sin(angle.to_double()),
                          cos(angle.to_double()));
    for(typename std::list<std::pair<std::pair<NT,NT>,SAVED_OBJECT> >::iterator
        iter = input_points_list.begin();
        iter != input_points_list.end();
        ++iter) {

      Point p(iter->first.first,iter->first.second);
      p = rotate(p);

      my_point<NT,SAVED_OBJECT> rotated_point(p,Point(iter->first.first,
                                iter->first.second),iter->second);

      l.push_back(rotated_point);
    }

    tree->build(l);

    //checking validity
    if(!tree->is_valid() )
      tree->dump();
    assert(tree->is_valid());

    return(tree);
  }

  inline NT squere(NT x) {return(x * x);}
  inline NT min(NT x,NT y) {return((x < y) ? x : y);}
  inline NT max(NT x,NT y) {return((x < y) ? y : x);}
  inline NT min(NT x1,NT x2,NT x3,NT x4,NT x5,NT x6) {return(min(min(min(x1,x2),
                min(x3,x4)),min(x5,x6)));}
  inline NT max(NT x1,NT x2,NT x3,NT x4,NT x5,NT x6) {return(max(max(max(x1,x2),
                max(x3,x4)),max(x5,x6)));}

  int get_kd_num(std::pair<std::pair<NT,NT>,std::pair<NT,NT> > &seg,int n)
  {
    int i;
    double x1 = seg.first.first.to_double();
    double y1 = seg.first.second.to_double();
    double x2 = seg.second.first.to_double();
    double y2 = seg.second.second.to_double();
    double alpha = atan((y2 -y1)/(x2 - x1));

    if(alpha < 0)
      alpha += PI / 2.0;

    if(alpha >= PI * (2 * n - 1) / (4 * n))
      i = 0;
    else {
      alpha += PI / (4 * n);
      i = int(2 * n * alpha / PI);
    }

    return(i);
  }

  void check_kd(int *kd_counter,int number_of_trees,
       std::list<std::pair<std::pair<NT,NT>,std::pair<NT,NT> > > &seg_list)
  {
    for(int i = 0;i < number_of_trees;++i)
      kd_counter[i] = 0;

    int kd_num;
    for(typename std::list<std::pair<std::pair<NT,NT>,std::pair<NT,NT> > >::
        iterator iter = seg_list.begin();iter != seg_list.end();++iter) {
      kd_num = get_kd_num(*iter,number_of_trees);
      kd_counter[kd_num]++;
    }
  }

public:
  Multiple_kd_tree(std::list<std::pair<std::pair<NT,NT>,SAVED_OBJECT> > 
                   &inp_points_list,int inp_number_of_trees,
                   std::list<std::pair<std::pair<NT,NT>,std::pair<NT,NT> > >
                   &seg_list) : 
    number_of_trees(inp_number_of_trees),input_points_list(inp_points_list)
  {
#ifdef TIMER
    CGAL::Timer t;
    t.start();
#endif

    kd_tree *kd;

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

#ifdef DEBUG
    int x = 0;
#endif
    int ind = 0;
    for(NT angle = 0;angle < HALF_PI;angle += HALF_PI / number_of_trees) {

#ifdef DEBUG
      std::cerr << "creating kd tree:  nu. " << ++x << endl;
#endif
      if(kd_counter[ind] >= (double)number_of_segments /
                            (double) number_of_trees / 2.0)
        kd = create_kd_tree(angle);
      else
        kd = NULL;

      kd_trees_list.push_back(std::pair<kd_tree *,NT>(kd,angle.to_double()));

      ++ind;
    }

#ifdef TIMER
    t.stop();

    std::cerr << endl << "Kd trees creation took " <<
                 t.time() << " seconds\n\n";
#endif
  }

  void get_intersecting_points(std::list<SAVED_OBJECT> &result_list,
                               Segment inp_s,NT unit_squere)
  {
    Segment s((inp_s.source().y() < inp_s.target().y()) ? inp_s.source() :
            inp_s.target(),(inp_s.source().y() < inp_s.target().y()) ?
            inp_s.target() : inp_s.source());
    NT x1 = s.source().x(),y1 = s.source().y(),
            x2 = s.target().x(),y2 = s.target().y();
    Point ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points
    Point op_ms1,op_ms2,op_ms3,op_ms4,op_ms5,op_ms6;// optimal min sums points

    // determine minkowski sum points
    NT half_unit_squere = unit_squere / 2.0;
    if(x1 < x2) {
      ms1 = Point(x1 - half_unit_squere - EPSILON,y1 -
            half_unit_squere - EPSILON);
      ms2 = Point(x1 - half_unit_squere - EPSILON,y1 +
            half_unit_squere + EPSILON);
      ms3 = Point(x1 + half_unit_squere + EPSILON,y1 -
            half_unit_squere - EPSILON);
      ms4 = Point(x2 + half_unit_squere + EPSILON,y2 -
            half_unit_squere - EPSILON);
      ms5 = Point(x2 + half_unit_squere + EPSILON,y2 +
            half_unit_squere + EPSILON);
      ms6 = Point(x2 - half_unit_squere - EPSILON,y2 +
            half_unit_squere + EPSILON);
    } else {
      ms1 = Point(x1 + half_unit_squere + EPSILON,y1 -
            half_unit_squere - EPSILON);
      ms2 = Point(x1 - half_unit_squere - EPSILON,y1 -
            half_unit_squere - EPSILON);
      ms3 = Point(x1 + half_unit_squere + EPSILON,y1 +
            half_unit_squere + EPSILON);
      ms4 = Point(x2 + half_unit_squere + EPSILON,y2 +
            half_unit_squere + EPSILON);
      ms5 = Point(x2 - half_unit_squere - EPSILON,y2 +
            half_unit_squere + EPSILON);
      ms6 = Point(x2 - half_unit_squere - EPSILON,y2 -
            half_unit_squere - EPSILON);
    }

    Transformation rotate(CGAL::ROTATION, sin((HALF_PI / number_of_trees).
        to_double()), cos((HALF_PI / number_of_trees).to_double()));
    Point p1,p2;
    std::list<my_point<NT,SAVED_OBJECT> > res;

    // find appropriate kd-tree
    NT size = -1,tmp_size,min_x,min_y,max_x,max_y;

    typename std::list<std::pair<kd_tree *,NT> >::iterator iter,right_iter;

    for(iter = kd_trees_list.begin();iter != kd_trees_list.end();++iter) {
      if(iter->first) {// skip nonexisting kd trees
        min_x = min(ms1.x(),ms2.x(),ms3.x(),ms4.x(),ms5.x(),ms6.x());
        min_y = min(ms1.y(),ms2.y(),ms3.y(),ms4.y(),ms5.y(),ms6.y());
        max_x = max(ms1.x(),ms2.x(),ms3.x(),ms4.x(),ms5.x(),ms6.x());
        max_y = max(ms1.y(),ms2.y(),ms3.y(),ms4.y(),ms5.y(),ms6.y());
        tmp_size = abs((max_x - min_x)*(max_y - min_y));

        if(size == -1 || tmp_size < size) {
          size = tmp_size;
          right_iter = iter;
          op_ms1 = ms1;
          op_ms2 = ms2;
          op_ms3 = ms3;
          op_ms4 = ms4;
          op_ms5 = ms5;
          op_ms6 = ms6;
        }
      }

      ms1 = rotate(ms1);
      ms2 = rotate(ms2);
      ms3 = rotate(ms3);
      ms4 = rotate(ms4);
      ms5 = rotate(ms5);
      ms6 = rotate(ms6);
    }

    // query
    p1 = Point(min(op_ms1.x(),op_ms2.x(),op_ms3.x(),op_ms4.x(),op_ms5.x(),
               op_ms6.x()),min(op_ms1.y(),op_ms2.y(),op_ms3.y(),op_ms4.y(),
               op_ms5.y(),op_ms6.y()));
    p2 = Point(max(op_ms1.x(),op_ms2.x(),op_ms3.x(),op_ms4.x(),op_ms5.x(),
               op_ms6.x()),max(op_ms1.y(),op_ms2.y(),op_ms3.y(),op_ms4.y(),
               op_ms5.y(),op_ms6.y()));
    my_point<NT,SAVED_OBJECT> point1(p1); 
    my_point<NT,SAVED_OBJECT> point2(p2);

    Box b(point1,point2,2);
    
    right_iter->first->search(std::back_inserter(res),b);
    // create result
    result_list.empty();
    for(typename std::list<my_point<NT,SAVED_OBJECT> >::iterator my_point_iter
        = res.begin();my_point_iter != res.end();++my_point_iter)
      result_list.push_back(my_point_iter->object);
  }
};

#endif // CGAL_SR_KD_2_H
