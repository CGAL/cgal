
#ifndef CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
#define CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H

#include <set.h>
#include <list.h>

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

private:

  struct y_ptr_bigger;
  struct x_ptr_bigger;

  class Point_data {
  public:

    Point_2 p;

    set<Point_data *,y_ptr_bigger> *right_tent;
    set<Point_data *,y_ptr_bigger> *left_tent;
    Point_type type;

    Point_data(const Point_2& p);

    Point_data(const Point_2& p,
	       set<Point_data *,y_ptr_bigger> *r_tent,
	       set<Point_data *,y_ptr_bigger> *l_tent);

    Point_data(const Point_2& p,
	       set<Point_data *,y_ptr_bigger> *r_tent,
	       set<Point_data *,y_ptr_bigger> *l_tent,
	       Point_type i_type);

    Point_data(const Point_data &other) : 
      p(other.p),right_tent(other.right_tent),left_tent(other.left_tent) {}

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

    bool x_bigger(Point_data *second) {
      T::Compare_x_2 cx;
      return (cx(p, second->p) == LARGER);
    }
    
    bool y_bigger(Point_data *second) {
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

  struct y_ptr_bigger
  {
    bool operator()(const Point_data *a,const Point_data *b) const
    {
      T::Less_yx_2 lyx;
      return lyx(a->p, b->p);
    }
  };

  struct x_ptr_bigger
  {
    bool operator()(const Point_data *a,const Point_data *b) const
    {
      T::Less_xy_2 lxy;
      return lxy(a->p, b->p);
    }
  };

  struct point_bigger
  {
    bool operator()(const Point_2 &a,const Point_2 &b) const
    {
      T::Less_xy_2 lxy;
      return lxy(a, b);
    }
  };

  typedef set<Point_data *,x_ptr_bigger> Point_data_set_of_x;
  typedef set<Point_data *,y_ptr_bigger> Point_data_set_of_y;


  Point_data_set_of_x x_sorted;
  Point_data_set_of_y y_sorted;
  NT min_x, min_y, max_x, max_y, min_x2, min_y2, max_x2, max_y2;
  NT biggest_rect_x0, biggest_rect_y0, biggest_rect_x1 ,biggest_rect_y1, biggest_rect_size;
  //  Polygon *polygon;

  Iso_rectangle_2 get_biggest_rect();
  void phase_1();
  void phase_1_on_x();
  void phase_1_on_y();
  void phase_2_on_bot();
  void phase_2_on_top();
  void phase_2_on_left();
  void phase_2_on_right();
  void phase_2();
  void phase_3();
  void check_for_bigger(NT x0, NT y0, NT x1, NT y1);
  void copy_y_Point_data_to_list(list<Point_data *> &Point_data_list);
  void copy_x_Point_data_to_list(list<Point_data *> &Point_data_list);
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
  void phase_3_check_for_bigger(Point_data_set_of_y::iterator iter,Point_data_set_of_y::iterator iter1,Point_data_set_of_y::iterator iter2,Point_data_set_of_y::iterator iter3,bool first_iter_is_right,bool second_iter_is_right,bool third_iter_is_right);
  void empty_tents();


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


  //typedef set<Point_2,point_bigger>::const_iterator const_iterator;
  // ctor
  Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b);

  // ctor
  //  Largest_empty_iso_rectangle_2(Polygon &inp_polygon);

  // add a point to data
  void 
  insert(const Point_2& p, Point_type i_type = REG);

  // remove a point from data
  bool 
  remove(const Point_2& p);

  // retrieve biggest rectangle
  Iso_rectangle_2 
  get_largest_empty_rectangle();

  // clear data(remove points)
  void 
  clear();

  // get a begin iterator to points
  const_iterator 
  begin();

  // get a after-the-end iterator to points
  const_iterator 
  end();


  // get bounding box
  Iso_rectangle_2 
  get_bounding_box();

  // dtor
  ~Largest_empty_iso_rectangle_2();
};

template<class T>
Largest_empty_iso_rectangle_2<T>::~Largest_empty_iso_rectangle_2()
{
  for(Point_data_set_of_x::iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter)
    delete(*iter);
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
  Point_data *po = new Point_data(_p);
  Point_data_set_of_x::iterator iter1 = x_sorted.find(po), iter2 = y_sorted.find(po);

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
Largest_empty_iso_rectangle_2<T>::check_for_bigger(NT x0, NT y0, NT x1, NT y1)
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
  // check if the rectangle represented by the parameters is bigger than the current one
  if(do_check && abs(x1 - x0) * abs(y1 - y0) > biggest_rect_size) {
    biggest_rect_size = abs(x1 - x0) * abs(y1 - y0);
    biggest_rect_x0 = x0;
    biggest_rect_y0 = y0;
    biggest_rect_x1 = x1;
    biggest_rect_y1 = y1;
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

  // traverse over all possibilities for finding a bigger empty rectangle
  // rectangles here touch the top and the buttom of the bounding box  
  while(iter != last_iter) {
    if((*iter)->type != TOP_RIGHT && (*iter)->type != TOP_LEFT) {// filter false points
      check_for_bigger((*prev_iter)->p.x(), min_y, (*iter)->p.x(), max_y);
      prev_iter = iter;
    }
    ++iter;
  }
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_1_on_y()
{
  Point_data_set_of_y::const_iterator iter = y_sorted.begin(),last_iter = y_sorted.end(),prev_iter = iter;
  ++iter;

  // filter false points
  while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT) {
    ++iter;
    ++prev_iter;
  }

  // traverse over all possibilities for finding a bigger empty rectangle
  // rectangles here touch the left and the right of the bounding box  
  while(iter != last_iter) {
    if((*iter)->type != BOT_RIGHT && (*iter)->type != TOP_RIGHT) {// filter false points
      check_for_bigger(min_x, (*prev_iter)->p.y(), max_x, (*iter)->p.y());
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
Largest_empty_iso_rectangle_2<T>::copy_x_Point_data_to_list(list<Point_data *> &point_data_list)
{
  for(Point_data_set_of_x::const_iterator iter = x_sorted.begin();iter != x_sorted.end();++iter)
    point_data_list.push_back(*iter);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::copy_y_Point_data_to_list(list<Point_data *> &point_data_list)
{
  for(Point_data_set_of_y::const_iterator iter = y_sorted.begin();
      iter != (Point_data_set_of_y::const_iterator)y_sorted.end();
      ++iter)
    point_data_list.push_back(*iter);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::tent(Point_data *first,Point_data *second)
{
  if(first->y_smaller(second))
    /* if(first->y < second->y)*/
    first->right_tent->insert(second);
  else
    second->left_tent->insert(first);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::tent(Point_data *first,Point_data *second,Point_data *third)
{
  first->right_tent->insert(second);
  third->left_tent->insert(second);
  if(first->y_smaller(third))
    /* if(first->y < third->y) */
    first->right_tent->insert(third);
  else
    third->left_tent->insert(first);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_top(list<Point_data *>::iterator &iter,
						   list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT || (*iter)->type == BOT_LEFT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_top(list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == BOT_LEFT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(list<Point_data *>::iterator &iter,
						   list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_bot(list<Point_data *>::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(Point_data_set_of_y::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(Point_data_set_of_y::iterator &iter,
						   Point_data_set_of_y::iterator &last)
{
  while(iter != last && ((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT))
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
Largest_empty_iso_rectangle_2<T>::get_next_for_left(list<Point_data *>::iterator &iter, 
						    list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_left(list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_right(list<Point_data *>::iterator &iter,
						     list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_LEFT || (*iter)->type == TOP_LEFT))
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
  copy_x_Point_data_to_list(Point_data_list);
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
    if((*iter1)->y_smaller(*iter2) && (*iter2)->y_bigger(*iter3)) {
      // Rectangles in phase 2 should be ignored for polygon
      //if(!polygon)
        check_for_bigger((*iter1)->p.x(), min_y, (*iter3)->p.x(), (*iter2)->p.y());
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
  copy_x_Point_data_to_list(Point_data_list);
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
    if((*iter1)->y_bigger(*iter2) && (*iter2)->y_smaller(*iter3)) {
      check_for_bigger((*iter1)->p.x(),max_y, (*iter3)->p.x(),(*iter2)->p.y());
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
  copy_y_Point_data_to_list(Point_data_list);
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
    if((*iter1)->x_smaller(*iter2) && (*iter2)->x_bigger(*iter3)) {
      check_for_bigger(min_x, (*iter1)->p.y(), (*iter2)->p.x(), (*iter3)->p.y());
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
  copy_y_Point_data_to_list(Point_data_list);
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
    if((*iter1)->x_bigger(*iter2) && (*iter2)->x_smaller(*iter3)) {
      check_for_bigger((*iter2)->p.x(), (*iter1)->p.y(), max_x, (*iter3)->p.y());
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
Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_biggest_rect()
{
  return(Iso_rectangle_2(Point_2(biggest_rect_x0, biggest_rect_y0),
			 Point_2(biggest_rect_x1,biggest_rect_y1)));
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_next_iter(Point_data_set_of_y::iterator &iter,
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
	/* if((*right_iter)->y < (*left_iter)->y) { */
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
Largest_empty_iso_rectangle_2<T>::phase_3_check_for_bigger(Point_data_set_of_y::iterator iter,
							   Point_data_set_of_y::iterator iter1,
							   Point_data_set_of_y::iterator iter2,
							   Point_data_set_of_y::iterator iter3,
							   bool first_iter_is_right,
							   bool second_iter_is_right,
							   bool third_iter_is_right)
{
  if(first_iter_is_right) {
    if(!second_iter_is_right)
      check_for_bigger((*iter2)->p.x(), (*iter)->p.y(), (*iter1)->p.x(), (*iter3)->p.y());
  } else
    if(second_iter_is_right)
      check_for_bigger((*iter1)->p.x(),(*iter)->p.y(),(*iter2)->p.x(),(*iter3)->p.y());
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::calls_for_tents(Point_data_set_of_y::iterator iter1,
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
Largest_empty_iso_rectangle_2<T>::calls_for_tents(Point_data_set_of_y::iterator iter1,
						  Point_data_set_of_y::iterator iter2)
{
  if((*iter1)->x_smaller(*iter2))
    /* if((*iter1)->p.x() < (*iter2)->p.x()) */
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
  Point_data_set_of_y::iterator iter1, iter2, iter3, right_iter, left_iter, last = last_iter;

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
      phase_3_check_for_bigger(iter,
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
  for(Point_data_set_of_x::const_iterator iter = x_sorted.begin();iter != x_sorted.end();++iter) {
    (*iter)->right_tent->clear();
    (*iter)->left_tent->clear();
  }
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_largest_empty_rectangle()
{
  if(x_sorted.size() == 4)
    return(get_bounding_box());

  biggest_rect_size = 0;

  // Rectangles in phase 1 should be ignored for polygon
  //if(!polygon)
    phase_1();

  phase_2();
  phase_3();
  Iso_rectangle_2 b = get_biggest_rect();
  empty_tents();

  return(b);
}

template<class T>
Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_bounding_box()
{
  return(Iso_rectangle_2(Point_2(min_x,min_y),
			 Point_2(max_x,max_y)));
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
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b)
{
  //  polygon = NULL;


  NT x0 = b.min().x(),y0 = b.min().y();
  NT x1 = b.max().x(), y1 = b.max().y();
  // determine extreme values of bounding box
  min_x2 = min_x = x0;
  min_y2 = min_y = y0;
  max_x2 = max_x = x1;
  max_y2 = max_y = y1;

  // add extreme points
  insert(Point_2(min_x - 0.000001,min_y - 0.000001),BOT_LEFT);
  insert(Point_2(max_x + 0.000001,min_y2 - 0.000001),BOT_RIGHT);
  insert(Point_2(min_x2 - 0.000001,max_y + 0.000001),TOP_LEFT);
  insert(Point_2(max_x2 + 0.000001,max_y2 + 0.000001),TOP_RIGHT);
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
  for(Point_data_set_of_x::iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter)
    delete(*iter);

  x_sorted.clear();
  y_sorted.clear();


  insert(Point_2(min_x - 0.000001,min_y - 0.000001),BOT_LEFT);
  insert(Point_2(max_x + 0.000001,min_y2 - 0.000001),BOT_RIGHT);
  insert(Point_2(min_x2 - 0.000001,max_y + 0.000001),TOP_LEFT);
  insert(Point_2(max_x2 + 0.000001,max_y2 + 0.000001),TOP_RIGHT);
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_first_two_iters(Point_data_set_of_y::iterator &iter1,
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
