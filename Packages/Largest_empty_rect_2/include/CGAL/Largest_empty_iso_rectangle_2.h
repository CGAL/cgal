// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Largest_empty_iso_rectangle_2.h
// package       : Largest_empty_iso_rectangle_2
// maintainer    : Eli Packer
// author(s)     : Eli Packer (algorithm), Andreas Fabri (cgal conformance)
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
#define CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H

#include <set.h>
#include <list.h>

#include <CGAL/quadruple.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE



template<class T>
class Largest_empty_iso_rectangle_2 {
  enum Point_type{REG, BOT_RIGHT, BOT_LEFT, TOP_LEFT, TOP_RIGHT};
  typedef typename T::FT NT;
  typedef T::Point_2               Point_2;
  //  typedef T::Vector_2              Vector_2; 
  typedef T::Iso_rectangle_2       Iso_rectangle_2;

  typedef T                        Traits;
private:

  bool cache_valid;
  struct y_ptr_larger;
  struct x_ptr_larger;

  class Point_data {
  public:

    Point_2 p;

    set<Point_data *,y_ptr_larger> *right_tent;
    set<Point_data *,y_ptr_larger> *left_tent;
    Point_type type;

    Point_data(const Point_2& p);

    Point_data(const Point_2& p,
	       set<Point_data *,y_ptr_larger> *r_tent,
	       set<Point_data *,y_ptr_larger> *l_tent);

    Point_data(const Point_2& p,
	       set<Point_data *,y_ptr_larger> *r_tent,
	       set<Point_data *,y_ptr_larger> *l_tent,
	       Point_type i_type);

    Point_data (const Point_data &other) : 
      p(other.p), right_tent(other.right_tent), left_tent(other.left_tent) 
    {}

    ~Point_data() {
      delete right_tent;
      delete left_tent;
    }

    bool x_smaller(Point_data *second) {
      T::Less_x_2 lx;
      return lx(p, second->p);
    }

    bool y_smaller(Point_data *second) {
      T::Less_yx_2 lyx;
      return lyx(p, second->p);
    }

    bool x_larger(Point_data *second) {
      T::Compare_x_2 cx;
      return (cx(p, second->p) == LARGER);
    }
    
    bool y_larger(Point_data *second) {
      T::Compare_y_2 cy;
      Comparison_result c = cy(p, second->p);
      if(c == LARGER) {
	return true;
      } else if (c == EQUAL) {
	T::Less_x_2 lx;
	return lx(second->p, p);
      } 
      return false;
    }
  }; 

  struct y_ptr_larger
  {
    bool operator()(const Point_data *a, const Point_data *b) const
    {
      T::Less_yx_2 lyx;
      return lyx(a->p, b->p);
    }
  };

  struct x_ptr_larger
  {
    bool operator()(const Point_data *a, const Point_data *b) const
    {
      T::Less_xy_2 lxy;
      return lxy(a->p, b->p);
    }
  };

  typedef set<Point_data *,x_ptr_larger> Point_data_set_of_x;
  typedef set<Point_data *,y_ptr_larger> Point_data_set_of_y;

  Point_data_set_of_x x_sorted;
  Point_data_set_of_y y_sorted;

  Point_2 bl_p, tr_p;

  // the points that define the largest empty iso-rectangle
  Point_2 left_p, bottom_p, right_p ,top_p; 

  NT largest_rect_size;
  //  Polygon *polygon;

  void phase_1();
  void phase_1_on_x();
  void phase_1_on_y();
  void phase_2_on_bot();
  void phase_2_on_top();
  void phase_2_on_left();
  void phase_2_on_right();
  void phase_2();
  void phase_3();
  void check_for_larger(const Point_2& px0, 
			const Point_2& py0, 
			const Point_2& px1,
			const Point_2& py1);
  void tent(Point_data *first, Point_data *second);
  void tent(Point_data *first, Point_data *second, Point_data *third);
  void get_next_for_top(list<Point_data *>::iterator &iter,list<Point_data *>::iterator &beyond);
  void get_prev_for_top(list<Point_data *>::iterator &iter);
  void get_next_for_bot(list<Point_data *>::iterator &iter,list<Point_data *>::iterator &beyond);
  void get_prev_for_bot(list<Point_data *>::iterator &iter);
  void get_next_for_bot(Point_data_set_of_y::iterator &iter);
  void get_next_for_bot(Point_data_set_of_y::iterator &iter,Point_data_set_of_y::iterator &last);
  void get_prev_for_bot(Point_data_set_of_y::iterator &iter);
  void get_next_for_left(list<Point_data *>::iterator &iter,list<Point_data *>::iterator &beyond);
  void get_prev_for_left(list<Point_data *>::iterator &iter);
  void get_next_for_right(list<Point_data *>::iterator &iter,list<Point_data *>::iterator &beyond);
  void get_prev_for_right(list<Point_data *>::iterator &iter);
  void determine_first_two_iters(Point_data_set_of_y::iterator &iter1,Point_data_set_of_y::iterator &iter2,Point_data_set_of_y::iterator &iter3,bool &first_iter_is_right,bool &second_iter_is_right,bool &third_iter_is_right);
  void determine_next_iter(Point_data_set_of_y::iterator &iter,Point_data_set_of_y::iterator &right_iter,Point_data_set_of_y::iterator &left_iter,Point_data_set_of_y::const_iterator right_iter_end,Point_data_set_of_y::const_iterator left_iter_end,bool &iter_is_right,bool &exist);
  void calls_for_tents(Point_data_set_of_y::iterator iter1,Point_data_set_of_y::iterator iter2);
  void calls_for_tents(Point_data_set_of_y::iterator iter1,Point_data_set_of_y::iterator iter2,Point_data_set_of_y::iterator iter3);
  void phase_2_update_y_sorted_list();
  void phase_3_check_for_larger(Point_data_set_of_y::iterator iter,Point_data_set_of_y::iterator iter1,Point_data_set_of_y::iterator iter2,Point_data_set_of_y::iterator iter3,bool first_iter_is_right,bool second_iter_is_right,bool third_iter_is_right);
  void empty_tents();
  void update();
  void init(const Point_2& bl, const Point_2& tr);
  // Auxiliary iterators for convenience

  template < class Node>
  struct Proj_point {
    typedef Node                  argument_type;
    typedef Point_2                 result_type;
    Point_2&       operator()( Node& x)       const { return x->p; }
    const Point_2& operator()( const Node& x) const { return x->p; }
};

public:

  // The following is an iterator adapter that allows us to enumerate 
  // the points in a set where they are stored
  typedef Iterator_project<Point_data_set_of_x::const_iterator, 
                           Proj_point<Point_data*>, 
                           const Point_2&,
                           const Point_2*> const_iterator;


  const Traits & geom_traits(){return Traits();};

  // ctor
  Largest_empty_iso_rectangle_2(const Point_2& bl, const Point_2& tr);

  // ctor
  Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b);

  // ctor
  //  Largest_empty_iso_rectangle_2(Polygon &inp_polygon);

  // add a point to data
  void 
  insert(const Point_2& p, Point_type i_type = REG);

  // and the STL standard member function for insertion:
  void
  push_back(const Point_2& _p)
  {
    insert(_p);
  }
  
  // and the insertion of an iterator range:
  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
    int n = 0;
    while(first != last){
      insert(*first);
      n++;
      ++first;
    }
    return n;
  }

  // remove a point from data
  bool 
  remove(const Point_2& p);

  Iso_rectangle_2
  get_bounding_box();

  // retrieve largest rectangle
  Iso_rectangle_2 
  get_largest_empty_iso_rectangle();

  // retrieve four points from the input that define the largest rectangle
  quadruple<Point_2, Point_2, Point_2, Point_2>
  get_left_bottom_right_top();

  // clear data(remove points)
  void 
  clear();

  // get a begin iterator to points
  const_iterator 
  begin();

  // get a after-the-end iterator to points
  const_iterator 
  end();


  // dtor
  ~Largest_empty_iso_rectangle_2();
};


template<class T>
struct Delete {
  void operator()(T it){
    delete(*it);
  }
};

template<class T>
Largest_empty_iso_rectangle_2<T>::~Largest_empty_iso_rectangle_2()
{
  
  for(Point_data_set_of_x::iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter)
    delete(*iter);
  
  /*
  // Why does this not compile ????
  for_each(x_sorted.begin(), 
	   x_sorted.end(), 
	   Delete<Point_data_set_of_x::iterator>());
  */
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Point_data::Point_data(const Point_2& _p) : p(_p),type(REG)
{
  right_tent = 0;
  left_tent = 0;
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Point_data::Point_data(const Point_2& _p,
							 Point_data_set_of_y *r_tent,
							 Point_data_set_of_y *l_tent)
  : p(_p),right_tent(r_tent),left_tent(l_tent),type(REG) 
{}

template<class T>
Largest_empty_iso_rectangle_2<T>::Point_data::Point_data(const Point_2& _p,
							 Point_data_set_of_y *r_tent,
							 Point_data_set_of_y *l_tent,
							 Point_type i_type) 
  : p(_p),right_tent(r_tent),left_tent(l_tent),type(i_type) 
{}


template<class T>
void
Largest_empty_iso_rectangle_2<T>::insert(const Point_2& _p,
					 Point_type i_type)
{
  cache_valid = false;
  Point_data_set_of_y *right_tent = new Point_data_set_of_y;
  Point_data_set_of_y *left_tent = new Point_data_set_of_y;
  Point_data *po = new Point_data(_p,right_tent,left_tent,i_type);

  x_sorted.insert(po);
  y_sorted.insert(po);

}



template<class T>
bool
Largest_empty_iso_rectangle_2<T>::remove(const Point_2& _p)
{
  cache_valid = false;
  Point_data *po = new Point_data(_p);
  Point_data_set_of_x::iterator iter1 = x_sorted.find(po),
                                iter2 = y_sorted.find(po);

  // point does not exist or a corner point
  if(iter1 == x_sorted.end() || iter1->type != REG)
    return(false);

  delete(iter->right_tent);
  delete(iter->left_tent);

  x_sorted.erase(iter1);
  y_sorted.erase(iter2);

  return(true);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::check_for_larger(const Point_2& px0,
						   const Point_2& py0,
						   const Point_2& px1,
						   const Point_2& py1)
{
  bool do_check = true;
  /*
  if(polygon) {
    NT bw = x1 - x0;
    NT bh = y1 - y0;
    NT bx = x0 + bw/2;
    NT by = y0 + bh/2;

    Point_2 center(bx,by);
    do_check = polygon->has_on_bounded_side(center);
  }
  */
  // check if the rectangle represented by the parameters is larger 
  //than the current one
  NT rect_size = CGAL_NTS abs(px1.x() - px0.x()) * CGAL_NTS abs(py1.y() - py0.y());
  if(do_check && rect_size > largest_rect_size) {
    largest_rect_size = rect_size;
    left_p = px0;
    bottom_p = py0;
    right_p = px1;
    top_p = py1;
  }
}

template<class T>
void
Largest_empty_iso_rectangle_2<T>::phase_1_on_x()
{
  Point_data_set_of_x::const_iterator iter = x_sorted.begin(),
                                      last_iter = x_sorted.end(),
                                      prev_iter = iter;
  ++iter;

  // filter false points
  while((*iter)->type == TOP_RIGHT || (*iter)->type == TOP_LEFT) {
    ++iter;
    ++prev_iter;
  }

  // traverse over all possibilities for finding a larger empty rectangle
  // rectangles here touch the top and the buttom of the bounding box  
  while(iter != last_iter) {
    // filter false points
    if((*iter)->type != TOP_RIGHT && (*iter)->type != TOP_LEFT) {
      check_for_larger((*prev_iter)->p, bl_p, (*iter)->p, tr_p);
      prev_iter = iter;
    }
    ++iter;
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_1_on_y()
{
  Point_data_set_of_y::const_iterator iter = y_sorted.begin(),
                                             last_iter = y_sorted.end(),
                                             prev_iter = iter;
  ++iter;

  // filter false points
  while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT) {
    ++iter;
    ++prev_iter;
  }

  // traverse over all possibilities for finding a larger empty rectangle
  // rectangles here touch the left and the right of the bounding box  
  while(iter != last_iter) {
    // filter false points
    if((*iter)->type != BOT_RIGHT && (*iter)->type != TOP_RIGHT) {
      check_for_larger(bl_p, (*prev_iter)->p, tr_p, (*iter)->p);
      prev_iter = iter;
    }
    ++iter;
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_1()
{
  phase_1_on_x();
  phase_1_on_y();
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::tent(Point_data *first, Point_data *second)
{
  if(first->y_smaller(second))
    first->right_tent->insert(second);
  else
    second->left_tent->insert(first);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::tent(Point_data *first,
				       Point_data *second,
				       Point_data *third)
{
  first->right_tent->insert(second);
  third->left_tent->insert(second);
  if(first->y_smaller(third))
    first->right_tent->insert(third);
  else
    third->left_tent->insert(first);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_top(
  list<Point_data *>::iterator &iter,
  list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT 
			   || (*iter)->type == BOT_LEFT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_top(
  list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == BOT_LEFT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  list<Point_data *>::iterator &iter,
  list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == TOP_LEFT 
			   || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_bot(
  list<Point_data *>::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  Point_data_set_of_y::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  Point_data_set_of_y::iterator &iter,
  Point_data_set_of_y::iterator &last)
{
  while(iter != last && ((*iter)->type == TOP_LEFT 
			 || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_bot(Point_data_set_of_y::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_left(
  list<Point_data *>::iterator &iter, 
  list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT 
			   || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_left(
  list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_right(
  list<Point_data *>::iterator &iter,
  list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_LEFT 
			   || (*iter)->type == TOP_LEFT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_right(list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_LEFT || (*iter)->type == TOP_LEFT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2_on_bot()
{
  list<Point_data *> Point_data_list;
  copy(x_sorted.begin(), x_sorted.end(), back_inserter(Point_data_list));
  list<Point_data *>::iterator iter1 = Point_data_list.begin(),
                               iter2,iter3,first_iter,
                               beyond = Point_data_list.end();
  int points_removed = 0, 
      size = Point_data_list.size();

  get_next_for_bot(iter1,beyond);
  first_iter = iter1;
  iter2 = iter1;
  ++iter2;
  get_next_for_bot(iter2,beyond);
  iter3 = iter2;
  ++iter3;
  get_next_for_bot(iter3,beyond);

  while(size - 4 > points_removed && iter3 != Point_data_list.end()) {
    if((*iter1)->y_smaller(*iter2) && (*iter2)->y_larger(*iter3)) {
      // Rectangles in phase 2 should be ignored for polygon
      //if(!polygon)
        check_for_larger((*iter1)->p, bl_p, (*iter3)->p, (*iter2)->p);
      tent(*iter1,*iter2,*iter3);
      ++points_removed;
      Point_data_list.erase(iter2);
      if(iter1 != first_iter) { // move back
        iter2 = iter1;
        --iter1;
        get_prev_for_bot(iter1);
      } else {
        iter2 = iter3;
        ++iter3;
        get_next_for_bot(iter3,beyond);
      }
    } else {// iter3 can't be last
      iter1 = iter2;
      iter2 = iter3;
      ++iter3;
      get_next_for_bot(iter3,beyond);
    }
  }
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2_on_top()
{
  list<Point_data *> Point_data_list;
   copy(x_sorted.begin(), x_sorted.end(), back_inserter(Point_data_list));

  list<Point_data *>::iterator iter1 = Point_data_list.begin(),
    iter2, iter3, first_iter,
    beyond = Point_data_list.end();
  int points_removed = 0,
    size = Point_data_list.size();

  get_next_for_top(iter1,beyond);
  iter2 = iter1;
  first_iter = iter1;
  ++iter2;
  get_next_for_top(iter2,beyond);
  iter3 = iter2;
  ++iter3;
  get_next_for_top(iter3,beyond);

  while(size - 4 > points_removed && iter3 != Point_data_list.end()) {
    if((*iter1)->y_larger(*iter2) && (*iter2)->y_smaller(*iter3)) {
      check_for_larger((*iter1)->p,tr_p, (*iter3)->p,(*iter2)->p);
      ++points_removed;
      Point_data_list.erase(iter2);
      if(iter1 != first_iter) { // move back
        iter2 = iter1;
        --iter1;
        get_prev_for_top(iter1);
      } else {
        iter2 = iter3;
        ++iter3;
        get_next_for_top(iter3,beyond);
      }
    } else {// iter3 can't be last
      iter1 = iter2;
      iter2 = iter3;
      ++iter3;
      get_next_for_top(iter3,beyond);
    }
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2_on_left()
{
  list<Point_data *> Point_data_list;
  copy(y_sorted.begin(), y_sorted.end(), back_inserter(Point_data_list));
  list<Point_data *>::iterator iter1 = Point_data_list.begin(),
    iter2, iter3, first_iter,
    beyond = Point_data_list.end();
  int points_removed = 0,size = Point_data_list.size();

  get_next_for_left(iter1,beyond);
  first_iter = iter1;
  iter2 = iter1;
  ++iter2;
  get_next_for_left(iter2,beyond);
  iter3 = iter2;
  ++iter3;
  get_next_for_left(iter3,beyond);

  while(size - 4 > points_removed && iter3 != Point_data_list.end()) {
    if((*iter1)->x_smaller(*iter2) && (*iter2)->x_larger(*iter3)) {
      check_for_larger(bl_p, (*iter1)->p, (*iter2)->p, (*iter3)->p);
      ++points_removed;
      Point_data_list.erase(iter2);
      if(iter1 != first_iter) { // move back
        iter2 = iter1;
        --iter1;
        get_prev_for_left(iter1);
      } else {
       iter2 = iter3;
        ++iter3;
        get_next_for_left(iter3,beyond);
      }
    } else {// iter3 can't be last
      iter1 = iter2;
      iter2 = iter3;
      ++iter3;
      get_next_for_left(iter3,beyond);
    }
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2_on_right()
{
  list<Point_data *> Point_data_list;
  copy(y_sorted.begin(), y_sorted.end(), back_inserter(Point_data_list));
  list<Point_data *>::iterator iter1 = Point_data_list.begin(),
    iter2, iter3, first_iter, 
    beyond = Point_data_list.end();
  int points_removed = 0,size = Point_data_list.size();

  get_next_for_right(iter1,beyond);
  first_iter = iter1;
  iter2 = iter1;
  ++iter2;
  get_next_for_right(iter2,beyond);
  iter3 = iter2;
  ++iter3;
  get_next_for_right(iter3,beyond);

  while(size - 4 > points_removed && iter3 != Point_data_list.end()) {
    if((*iter1)->x_larger(*iter2) && (*iter2)->x_smaller(*iter3)) {
      check_for_larger((*iter2)->p, (*iter1)->p, tr_p, (*iter3)->p);
      ++points_removed;
      Point_data_list.erase(iter2);
      if(iter1 != first_iter) { // move back
        iter2 = iter1;
        --iter1;
        get_prev_for_right(iter1);
      } else {
        iter2 = iter3;
        ++iter3;
        get_next_for_right(iter3,beyond);
      }
    } else {// iter3 can't be last
      iter1 = iter2;
      iter2 = iter3;
      ++iter3;
      get_next_for_right(iter3,beyond);
    }
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2()
{
  // Rectangles in phase 2 should be ignored for polygon
  //if(!polygon) {
    phase_2_on_top();
    phase_2_on_left();
    phase_2_on_right();
    //}

  // Done only for building tents for phase 3
  phase_2_on_bot();
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_next_iter(
  Point_data_set_of_y::iterator &iter,
  Point_data_set_of_y::iterator &right_iter,
  Point_data_set_of_y::iterator &left_iter,
  Point_data_set_of_y::const_iterator right_iter_end,
  Point_data_set_of_y::const_iterator left_iter_end,
  bool &iter_is_right,
  bool &exist)
{
  if((Point_data_set_of_y::const_iterator)right_iter != right_iter_end) {
    if((Point_data_set_of_y::const_iterator)left_iter != left_iter_end) {
      if((*right_iter)->y_smaller(*left_iter)) {
        iter = right_iter;
        iter_is_right = true;
        ++right_iter;
      } else {
        iter = left_iter;
        iter_is_right = false;
        ++left_iter;
      }
    } else {
      iter = right_iter;
      iter_is_right = true;
      ++right_iter;
    }
  } else { 
    if((Point_data_set_of_y::const_iterator)left_iter != left_iter_end) {
      iter = left_iter;
      iter_is_right = false;
      ++left_iter;
     } else
      exist = false;
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_3_check_for_larger(
  Point_data_set_of_y::iterator iter,
  Point_data_set_of_y::iterator iter1,
  Point_data_set_of_y::iterator iter2,
  Point_data_set_of_y::iterator iter3,
  bool first_iter_is_right,
  bool second_iter_is_right,
  bool third_iter_is_right)
{
  if(first_iter_is_right) {
    if(!second_iter_is_right)
      check_for_larger((*iter2)->p, (*iter)->p, (*iter1)->p, (*iter3)->p);
  } else
    if(second_iter_is_right)
      check_for_larger((*iter1)->p,(*iter)->p,(*iter2)->p,(*iter3)->p);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::calls_for_tents(
  Point_data_set_of_y::iterator iter1,
  Point_data_set_of_y::iterator iter2,
  Point_data_set_of_y::iterator iter3)
{
  bool first_is_right_to_second = (*iter1)->x_smaller(*iter2);
  bool second_is_right_to_third = (*iter2)->x_smaller(*iter3);

  if(first_is_right_to_second) {
    if(second_is_right_to_third) {
      tent(*iter1,*iter2);
      tent(*iter2,*iter3);
    } else {
      tent(*iter1,*iter3,*iter2);
    }
  } else {
    if(second_is_right_to_third) {
      tent(*iter2,*iter3,*iter1);
    }
    else {
      tent(*iter2,*iter1);
      tent(*iter3,*iter2);
    }
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::calls_for_tents(
  Point_data_set_of_y::iterator iter1,
  Point_data_set_of_y::iterator iter2)
{
  if((*iter1)->x_smaller(*iter2))
    tent(*iter1,*iter2);
  else
    tent(*iter2,*iter1);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_3()
{
  bool first_iter_is_right, second_iter_is_right, third_iter_is_right,
    first_exist,second_exist,third_exist;
  Point_data_set_of_y::iterator iter, last_iter = y_sorted.end();
  Point_data_set_of_y::iterator iter1, iter2, iter3, 
                                right_iter, left_iter, last = last_iter;

  --last_iter;
  --last_iter;
  for(iter = y_sorted.begin();iter != last_iter;++iter) {
    get_next_for_bot(iter,last);
    if(iter == last)
      return;
    first_exist = true;
    second_exist = true;
    third_exist = true;

    right_iter = (*iter)->right_tent->begin();
    left_iter = (*iter)->left_tent->begin();
    determine_next_iter(iter1,
			right_iter,
			left_iter,
			(*iter)->right_tent->end(),
			(*iter)->left_tent->end(),
			first_iter_is_right,
			first_exist);
    determine_next_iter(iter2,right_iter,left_iter,
			(*iter)->right_tent->end(),
			(*iter)->left_tent->end(),
			second_iter_is_right,
			second_exist);
    determine_next_iter(iter3,
			right_iter,
			left_iter,
			(*iter)->right_tent->end(),
			(*iter)->left_tent->end(),
			third_iter_is_right,
			third_exist);
    bool had_three = false;

    while(third_exist) {
      had_three = true;
      phase_3_check_for_larger(iter,
			       iter1,
			       iter2,
			       iter3,
			       first_iter_is_right,
			       second_iter_is_right,
			       third_iter_is_right);
      calls_for_tents(iter1, iter2, iter3);
      determine_first_two_iters(iter1,
				iter2,
				iter3,
				first_iter_is_right,
				second_iter_is_right,
				third_iter_is_right);
      determine_next_iter(iter3,
			  right_iter,left_iter,
			  (*iter)->right_tent->end(),
			  (*iter)->left_tent->end(),
			  third_iter_is_right,
			  third_exist);
    }

    if(!had_three && second_exist)
      calls_for_tents(iter1, iter2);
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::empty_tents()
{
  for(Point_data_set_of_x::const_iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter) {
    (*iter)->right_tent->clear();
    (*iter)->left_tent->clear();
  }
}
template<class T>
Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_bounding_box()
{
  return(Iso_rectangle_2(bl_p, tr_p));
}

/* Performs the computation if the cache is invalid
 *
 */
template<class T>
void
Largest_empty_iso_rectangle_2<T>::update()
{
  if(! cache_valid){
    largest_rect_size = 0;

    // Rectangles in phase 1 should be ignored for polygon
    //if(!polygon)
    phase_1();

    phase_2();
    phase_3();
    empty_tents();
    cache_valid = true;
  }
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_largest_empty_iso_rectangle()
{
  if(x_sorted.size() == 4) {
    return(get_bounding_box());
  }
  update();
  return(Iso_rectangle_2(Point_2(left_p.x(), bottom_p.y()),
			 Point_2(right_p.x(),top_p.y())));
}

/* Some applications might be more interested in the four points
 * that are from the input and that define the empty rectangle
 */
template<class T>
quadruple<Largest_empty_iso_rectangle_2<T>::Point_2,
          Largest_empty_iso_rectangle_2<T>::Point_2,
          Largest_empty_iso_rectangle_2<T>::Point_2,
          Largest_empty_iso_rectangle_2<T>::Point_2>
Largest_empty_iso_rectangle_2<T>::get_left_bottom_right_top()
{
  if(x_sorted.size() == 4) {
    return(make_quadruple(bl_p, bl_p, tr_p, tr_p));
  }
  update();
  return(make_quadruple(left_p, bottom_p, right_p, top_p));
}

/*
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(Polygon &inp_polygon)
{
  polygon = new Polygon(inp_polygon);

  // determine extreme values of bounding box
  min_x2 = min_x = polygon->left_vertex()->x();
  min_y2 = min_y = polygon->bottom_vertex()->y();
  max_x2 = max_x = polygon->right_vertex()->x();
  max_y2 = max_y = polygon->top_vertex()->y();

  // add extreme points
  insert(Point_2(min_x - 0.000001,min_y - 0.000001),BOT_LEFT);
  insert(Point_2(max_x + 0.000001,min_y2 - 0.000001),BOT_RIGHT);
  insert(Point_2(min_x2 - 0.000001,max_y + 0.000001),TOP_LEFT);
  insert(Point_2(max_x2 + 0.000001,max_y2 + 0.000001),TOP_RIGHT);

  // insert the polygon 
  Polygon::Vertex_iterator it = polygon->vertices_begin();
  Polygon::Vertex_iterator next = it;

  Point_2 p,q;
  Point_2 p0 = *it;

  while(it != polygon->vertices_end()) {
    insert(*it);
    ++next;
    if(next == polygon->vertices_end())
      q = p0;
    else
      q = *next;
    p = *it;

    // add some points on the segment 
    Vector_2 v = (q - p)/6;
    for(int j = 1; j < 6; j++) {
      insert(p + j * v);
    }    

    ++it;
  }
}
*/  


template<class T>
void
Largest_empty_iso_rectangle_2<T>::init(const Point_2& bl, const Point_2& tr)
{
  //  polygon = NULL;

  // determine extreme values of bounding box
  bl_p = bl;
  tr_p = tr;

  // add extreme points
  insert(Point_2(bl.x() - 0.000001, bl.y() - 0.000001), BOT_LEFT);
  insert(Point_2(tr.x() + 0.000001, bl.y() - 0.000001), BOT_RIGHT);
  insert(Point_2(bl.x() - 0.000001, tr.y() + 0.000001), TOP_LEFT);
  insert(Point_2(tr.x() + 0.000001, tr.y() + 0.000001), TOP_RIGHT);
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(const Point_2& bl, const Point_2& tr)
  : cache_valid(false)
{
  // precondition: bl and tr
  init(bl, tr);
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b)
  : cache_valid(false)
{
  init(b.min(), b.max());
}


template<class T>
Largest_empty_iso_rectangle_2<T>::const_iterator 
Largest_empty_iso_rectangle_2<T>::begin()
{
  return const_iterator(x_sorted.begin());
}

template<class T>
Largest_empty_iso_rectangle_2<T>::const_iterator 
Largest_empty_iso_rectangle_2<T>::end()
{
    return const_iterator(x_sorted.end());
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::clear()
{
  cache_valid = false;
  for(Point_data_set_of_x::iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter)
    delete(*iter);

  x_sorted.clear();
  y_sorted.clear();

  insert(Point_2(bl_p.x() - 0.000001, bl_p.y() - 0.000001),BOT_LEFT);
  insert(Point_2(tr_p.x() + 0.000001, bl_p.y() - 0.000001),BOT_RIGHT);
  insert(Point_2(bl_p.x() - 0.000001, tr_p.y() + 0.000001),TOP_LEFT);
  insert(Point_2(tr_p.x() + 0.000001, tr_p.y() + 0.000001),TOP_RIGHT);
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_first_two_iters(
  Point_data_set_of_y::iterator &iter1,
  Point_data_set_of_y::iterator &iter2,
  Point_data_set_of_y::iterator &iter3,
  bool &first_iter_is_right,
  bool &second_iter_is_right,
  bool &third_iter_is_right)
{
  if(first_iter_is_right) {
    if(second_iter_is_right) {
      iter1 = iter2;
      iter2 = iter3;
      first_iter_is_right = second_iter_is_right;
      second_iter_is_right = third_iter_is_right;
    } else {
      if(third_iter_is_right) {
        iter1 = iter2;
        iter2 = iter3;
        first_iter_is_right = second_iter_is_right;
        second_iter_is_right = third_iter_is_right;  
      } else {
        iter2 = iter3;
        second_iter_is_right = third_iter_is_right;  
      }
    }
  } else {
    if(second_iter_is_right) {
      if(third_iter_is_right) {
        iter2 = iter3;
        second_iter_is_right = third_iter_is_right;  
      } else {
        iter1 = iter2;
        iter2 = iter3;
        first_iter_is_right = second_iter_is_right;
        second_iter_is_right = third_iter_is_right;  
      }
    } else {
      iter1 = iter2;
      iter2 = iter3;
      first_iter_is_right = second_iter_is_right;
      second_iter_is_right = third_iter_is_right;  
    }
  }
}


CGAL_END_NAMESPACE
    

#endif // CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
