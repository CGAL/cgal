// =======================================================================
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
// coordinator   : Tel-Aviv University(Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
#define CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H

#include <set>
#include <list>

#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

CGAL_BEGIN_NAMESPACE


template<class T>
class Largest_empty_iso_rectangle_2 {
public: 
  enum Point_type{REG, BOT_RIGHT, BOT_LEFT, TOP_LEFT, TOP_RIGHT};
  typedef typename T::FT NT;
  typedef typename T::Point_2               Point;
  typedef typename T::Iso_rectangle_2       Iso_rectangle_2;
  typedef T                        Traits;

private:

  bool cache_valid;
  Traits _gt;
  class Less_yx;
  class Less_xy;

  class Point_data {
  public:

    Point p;

    std::set<Point_data *,Less_yx> *right_tent;
    std::set<Point_data *,Less_yx> *left_tent;
    Point_type type;

    Point_data(const Point& _p) 
      :  p(_p),type(REG)
    {
      right_tent = 0;
      left_tent = 0;
    }

    Point_data(const Point& _p,
	       std::set<Point_data *,Less_yx> *r_tent,
	       std::set<Point_data *,Less_yx> *l_tent)
      : p(_p),right_tent(r_tent),left_tent(l_tent),type(REG) 
    {}

    Point_data(const Point& _p,
	       std::set<Point_data *,Less_yx> *r_tent,
	       std::set<Point_data *,Less_yx> *l_tent,
	       Point_type i_type)
      : p(_p),right_tent(r_tent),left_tent(l_tent),type(i_type) 
    {}


    Point_data (const Point_data &other) : 
      p(other.p), right_tent(other.right_tent), left_tent(other.left_tent) 
    {}

    ~Point_data() {
      delete right_tent;
      delete left_tent;
    }

  }; 

  // was x_smaller
  bool less_xy(const Point_data *a, const Point_data *b) const
  {
    return traits().less_xy_2_object()(a->p, b->p);
  }

  // was y_smaller
  bool less_yx(const Point_data *a, const Point_data *b) {
    return traits().less_yx_2_object()(a->p, b->p);
    }

  // was x_larger
  bool larger_xy(const Point_data *a, const Point_data *b) {
    return traits().compare_xy_2_object()(a->p, b->p) == LARGER;
  }

  // was y_larger
  bool larger_yx(const Point_data *a, const Point_data *b) {

    Comparison_result c = traits().compare_y_2_object()(a->p, b->p);
    if(c == LARGER) {
      return true;
    } else if (c == EQUAL) {
      return traits().less_x_2_object()(b->p, a->p); 
    } 
    return false;
  }

  class Less_yx
  {
  private:
    Traits gt;

    Less_yx(){}

  public:
    Less_yx(const Traits& t)
      : gt(t)
    {}

    bool operator()(const Point_data *a, const Point_data *b) const
    {
      return gt.less_yx_2_object()(a->p, b->p);
    }
  };


  class Less_xy
  {
  private:
    Traits gt;

  public:
    Less_xy(const Traits& t)
      : gt(t)
    {}

    bool operator()(const Point_data *a, const Point_data *b) const
    {
      return gt.less_xy_2_object()(a->p, b->p);
    }
  };

  typedef std::set<Point_data *,Less_xy> Point_data_set_of_x;
  typedef std::set<Point_data *,Less_yx> Point_data_set_of_y;

  Point_data_set_of_x x_sorted;
  Point_data_set_of_y y_sorted;

  Point bl_p, tr_p;

  // the points that define the largest empty iso-rectangle
  Point left_p, bottom_p, right_p ,top_p; 

  NT largest_rect_size;


  bool insert(const Point& _p,Point_type i_type);
  void phase_1();
  void phase_1_on_x();
  void phase_1_on_y();
  void phase_2_on_bot();
  void phase_2_on_top();
  void phase_2_on_left();
  void phase_2_on_right();
  void phase_2();
  void phase_3();
  void check_for_larger(const Point& px0, 
			const Point& py0, 
			const Point& px1,
			const Point& py1);
  void tent(Point_data *first, Point_data *second);
  void tent(Point_data *first, Point_data *second, Point_data *third);
  void get_next_for_top(typename std::list<Point_data *>::iterator &iter,
			typename std::list<Point_data *>::iterator &beyond);
  void get_prev_for_top(typename std::list<Point_data *>::iterator &iter);
  void get_next_for_bot(typename std::list<Point_data *>::iterator &iter,
			typename std::list<Point_data *>::iterator &beyond);
  void get_prev_for_bot(typename std::list<Point_data *>::iterator &iter);
  void get_next_for_bot(typename Point_data_set_of_y::iterator &iter);
  void get_next_for_bot(typename Point_data_set_of_y::iterator &iter,
			typename Point_data_set_of_y::iterator &last);
  void get_prev_for_bot(typename Point_data_set_of_y::iterator &iter);
  void get_next_for_left(typename std::list<Point_data *>::iterator &iter,
			 typename std::list<Point_data *>::iterator &beyond);
  void get_prev_for_left(typename std::list<Point_data *>::iterator &iter);
  void get_next_for_right(typename std::list<Point_data *>::iterator &iter,
			  typename std::list<Point_data *>::iterator &beyond);
  void get_prev_for_right(typename std::list<Point_data *>::iterator &iter);
  void determine_first_two_iters(typename Point_data_set_of_y::iterator &iter1,
				 typename Point_data_set_of_y::iterator &iter2,
				 typename Point_data_set_of_y::iterator &iter3,
				 bool &first_iter_is_right,
				 bool &second_iter_is_right,
				 bool &third_iter_is_right);
  void determine_next_iter(
		    typename Point_data_set_of_y::iterator &iter,
		    typename Point_data_set_of_y::iterator &right_iter,
		    typename Point_data_set_of_y::iterator &left_iter,
		    typename Point_data_set_of_y::const_iterator right_iter_end,
		    typename Point_data_set_of_y::const_iterator left_iter_end,
		    bool &iter_is_right,
		    bool &exist);

  void calls_for_tents(typename Point_data_set_of_y::iterator iter1,
		       typename Point_data_set_of_y::iterator iter2);
  void calls_for_tents(typename Point_data_set_of_y::iterator iter1,
		       typename Point_data_set_of_y::iterator iter2,
		       typename Point_data_set_of_y::iterator iter3);
  void phase_2_update_y_sorted_list();
  void phase_3_check_for_larger(typename Point_data_set_of_y::iterator iter,
				typename Point_data_set_of_y::iterator iter1,
				typename Point_data_set_of_y::iterator iter2,
				typename Point_data_set_of_y::iterator iter3,
				bool first_iter_is_right,
				bool second_iter_is_right);
  void empty_tents();
  void update();
  void init(const Point& bl, const Point& tr);
  void copy_memory(const Largest_empty_iso_rectangle_2<T>& ler);
  void free_memory();

  // Auxiliary iterators for convenience

  template < class Node>
  struct Proj_point {
    typedef Node                  argument_type;
    typedef Point                 result_type;
    Point&       operator()( Node& x)       const { return x->p; }
    const Point& operator()( const Node& x) const { return x->p; }
};

public:

  // The following is an iterator adapter that allows us to enumerate 
  // the points in a set where they are stored
  typedef Iterator_project<typename Point_data_set_of_x::const_iterator, 
                           Proj_point<Point_data*>, 
                           const Point&,
                           const Point*> const_iterator;


  const Traits & traits() const {return _gt;};

  // ctor
  Largest_empty_iso_rectangle_2(const Point& bl, const Point& tr);

  // ctor
  Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b);

  // ctor
  Largest_empty_iso_rectangle_2();


  // add a point to data
  bool
  insert(const Point& p);

  // and the STL standard member function for insertion:
  void
  push_back(const Point& _p)
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
      if(insert(*first)){
	n++;
      }
      ++first;
    }
    return n;
  }

  // remove a point from data
  bool 
  remove(const Point& p);

  Iso_rectangle_2
  get_bounding_box();

  // retrieve largest rectangle
  Iso_rectangle_2 
  get_largest_empty_iso_rectangle();

  // retrieve four points from the input that define the largest rectangle
  Quadruple<Point, Point, Point, Point>
  get_left_bottom_right_top();

  // clear data(remove points)
  void 
  clear();

  // get a begin iterator to points
  const_iterator 
  begin()  const;

  // get a after-the-end iterator to points
  const_iterator
  end()  const;


  // dtor
  ~Largest_empty_iso_rectangle_2();

  // operator=
  Largest_empty_iso_rectangle_2<T>&
    operator =(const Largest_empty_iso_rectangle_2<T>& ler);

  // cctor
  Largest_empty_iso_rectangle_2<T>(
	       const Largest_empty_iso_rectangle_2<T>& ler);
  
};


template<class Ptr>
struct Delete {
  void operator()(Ptr ptr)const {
    delete(ptr);
  }
};

template<class T>
void Largest_empty_iso_rectangle_2<T>::free_memory()
{
  std::for_each(x_sorted.begin(), 
	   x_sorted.end(), 
	   Delete<Point_data*>());

  x_sorted.clear();
  y_sorted.clear();
}

template<class T>
Largest_empty_iso_rectangle_2<T>::~Largest_empty_iso_rectangle_2()
{
  free_memory();
}

template<class T>
void Largest_empty_iso_rectangle_2<T>::
  copy_memory(const Largest_empty_iso_rectangle_2<T>& ler)
{
  // copy bounding box
  bl_p = ler.bl_p;
  tr_p = ler.tr_p;

    // copy points
  for(typename Point_data_set_of_x::const_iterator iter = ler.x_sorted.begin();
        iter != ler.x_sorted.end();
	++iter) {
      if((*iter)->type == REG)
	insert((*iter)->p);
      else
	insert((*iter)->p,(*iter)->type);
  }
}

template<class T>
Largest_empty_iso_rectangle_2<T>&
Largest_empty_iso_rectangle_2<T>::operator =(
               const Largest_empty_iso_rectangle_2<T>& ler)
{
  if(this != &ler) {
    free_memory();
    copy_memory(ler);
  }

  return *this;
}

template<class T>
Largest_empty_iso_rectangle_2<T>::
Largest_empty_iso_rectangle_2(
               const Largest_empty_iso_rectangle_2<T>& ler)
: cache_valid(false), _gt(),
  x_sorted(Less_xy(traits())),
  y_sorted(Less_yx(traits()))
{
  copy_memory(ler);
}





template<class T>
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point& _p)
{
  // check that the point is inside the bounding box 
  if(_p.x() <= bl_p.x() || _p.x() >= tr_p.x() ||
     _p.y() <= bl_p.y() || _p.y() >= tr_p.y())
    return(false);

  // check that the point is not already inserted
  Point_data *po = new Point_data(_p);
  typename Point_data_set_of_x::iterator iter = x_sorted.find(po);
  delete(po);

  if(iter != x_sorted.end())
    return(false);

  cache_valid = false;
  Point_data_set_of_y *right_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data_set_of_y *left_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  po = new Point_data(_p,right_tent,left_tent,REG);

  x_sorted.insert(po);
  y_sorted.insert(po);
  return(true);
}



template<class T>
bool
Largest_empty_iso_rectangle_2<T>::remove(const Point& _p)
{
  cache_valid = false;
  Point_data *po = new Point_data(_p);
  typename Point_data_set_of_x::iterator iter1 = x_sorted.find(po);
  typename Point_data_set_of_y::iterator iter2 = y_sorted.find(po);

  // point does not exist or a corner point
  if(iter1 == x_sorted.end() || (*iter1)->type != REG)
    return(false);

  delete((*iter1)->right_tent);
  delete((*iter2)->left_tent);

  x_sorted.erase(iter1);
  y_sorted.erase(iter2);

  return(true);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::check_for_larger(const Point& px0,
						   const Point& py0,
						   const Point& px1,
						   const Point& py1)
{
  bool do_check = true;

  // check if the rectangle represented by the parameters is larger 
  //than the current one
  NT rect_size =
    CGAL_NTS abs(px1.x() - px0.x()) * CGAL_NTS abs(py1.y() - py0.y());
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
  typename Point_data_set_of_x::const_iterator iter = x_sorted.begin(),
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
  typename Point_data_set_of_y::const_iterator iter = y_sorted.begin(),
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
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point& _p,
					 Point_type i_type)
{
  // check that the point is inside the bounding box 
  if((i_type == REG) && (_p.x() <= bl_p.x() || _p.x() >= tr_p.x() ||
     _p.y() <= bl_p.y() || _p.y() >= tr_p.y()))
    return(false);

  // check that the point is not already inserted
  Point_data *po = new Point_data(_p);
  typename Point_data_set_of_x::iterator iter = x_sorted.find(po);
  delete(po);

  if(iter != x_sorted.end())
    return(false);

  cache_valid = false;
  Point_data_set_of_y *right_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data_set_of_y *left_tent = 
    new Point_data_set_of_y(Less_yx(traits()));
  po = new Point_data(_p,right_tent,left_tent,i_type);

  x_sorted.insert(po);
  y_sorted.insert(po);
  return(true);
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
  if(less_yx(first, second))
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
  if(less_yx(first, third))
    first->right_tent->insert(third);
  else
    third->left_tent->insert(first);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_top(
  typename std::list<Point_data *>::iterator &iter,
  typename std::list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT 
			   || (*iter)->type == BOT_LEFT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_top(
  typename std::list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == BOT_LEFT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  typename std::list<Point_data *>::iterator &iter,
  typename std::list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == TOP_LEFT 
			   || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_bot(
  typename std::list<Point_data *>::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  typename Point_data_set_of_y::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_bot(
  typename Point_data_set_of_y::iterator &iter,
  typename Point_data_set_of_y::iterator &last)
{
  while(iter != last && ((*iter)->type == TOP_LEFT 
			 || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_bot(
    typename Point_data_set_of_y::iterator &iter)
{
  while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_left(
  typename std::list<Point_data *>::iterator &iter, 
  typename std::list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_RIGHT 
			   || (*iter)->type == TOP_RIGHT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_left(
  typename std::list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_next_for_right(
  typename std::list<Point_data *>::iterator &iter,
  typename std::list<Point_data *>::iterator &beyond)
{
  while(iter != beyond && ((*iter)->type == BOT_LEFT 
			   || (*iter)->type == TOP_LEFT))
    ++iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::get_prev_for_right(
    typename std::list<Point_data *>::iterator &iter)
{
  while((*iter)->type == BOT_LEFT || (*iter)->type == TOP_LEFT)
    --iter;
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::phase_2_on_bot()
{
  std::list<Point_data *> Point_data_list;
  std::copy(x_sorted.begin(), 
	    x_sorted.end(), 
	    std::back_inserter(Point_data_list));
  typename std::list<Point_data *>::iterator iter1 = Point_data_list.begin(),
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
    if(less_yx(*iter1, *iter2) && larger_yx(*iter2, *iter3)) {

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
  std::list<Point_data *> Point_data_list;
   std::copy(x_sorted.begin(),
	     x_sorted.end(),
	     std::back_inserter(Point_data_list));

  typename std::list<Point_data *>::iterator iter1 = Point_data_list.begin(),
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
    if(larger_yx(*iter1, *iter2) && less_yx(*iter2, *iter3)) {
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
  std::list<Point_data *> Point_data_list;
  std::copy(y_sorted.begin(), 
	    y_sorted.end(), 
	    std::back_inserter(Point_data_list));
  typename std::list<Point_data *>::iterator iter1 = Point_data_list.begin(),
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
    if(less_xy(*iter1, *iter2) && larger_xy(*iter2, *iter3)) {
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
  std::list<Point_data *> Point_data_list;
  std::copy(y_sorted.begin(), 
	    y_sorted.end(), 
	    std::back_inserter(Point_data_list));
  typename std::list<Point_data *>::iterator iter1 = Point_data_list.begin(),
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
    if(larger_xy(*iter1, *iter2) && less_xy(*iter2, *iter3)) {
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
    phase_2_on_top();
    phase_2_on_left();
    phase_2_on_right();

  // Done only for building tents for phase 3
  phase_2_on_bot();
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_next_iter(
  typename Point_data_set_of_y::iterator &iter,
  typename Point_data_set_of_y::iterator &right_iter,
  typename Point_data_set_of_y::iterator &left_iter,
  typename Point_data_set_of_y::const_iterator right_iter_end,
  typename Point_data_set_of_y::const_iterator left_iter_end,
  bool &iter_is_right,
  bool &exist)
{
  if((typename Point_data_set_of_y::const_iterator)right_iter 
     != right_iter_end) {
    if((typename Point_data_set_of_y::const_iterator)left_iter 
       != left_iter_end) {
      if(less_yx(*right_iter, *left_iter)) {
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
    if((typename Point_data_set_of_y::const_iterator)left_iter 
       != left_iter_end) {
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
  typename Point_data_set_of_y::iterator iter,
  typename Point_data_set_of_y::iterator iter1,
  typename Point_data_set_of_y::iterator iter2,
  typename Point_data_set_of_y::iterator iter3,
  bool first_iter_is_right,
  bool second_iter_is_right)
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
  typename Point_data_set_of_y::iterator iter1,
  typename Point_data_set_of_y::iterator iter2,
  typename Point_data_set_of_y::iterator iter3)
{
  bool first_is_right_to_second = less_xy(*iter1, *iter2);
  bool second_is_right_to_third = less_xy(*iter2, *iter3);

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
  typename Point_data_set_of_y::iterator iter1,
  typename Point_data_set_of_y::iterator iter2)
{
  if(less_xy(*iter1, *iter2))
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
  typename Point_data_set_of_y::iterator iter, last_iter = y_sorted.end();
  typename Point_data_set_of_y::iterator iter1, iter2, iter3, 
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
			       second_iter_is_right);
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
  for(typename Point_data_set_of_x::const_iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter) {
    (*iter)->right_tent->clear();
    (*iter)->left_tent->clear();
  }
}
template<class T>
typename Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
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

    phase_1();

    phase_2();
    phase_3();
    empty_tents();
    cache_valid = true;
  }
}

template<class T>
typename Largest_empty_iso_rectangle_2<T>::Iso_rectangle_2
Largest_empty_iso_rectangle_2<T>::get_largest_empty_iso_rectangle()
{
  if(x_sorted.size() == 4) {
    return(get_bounding_box());
  }
  update();
  return(Iso_rectangle_2(Point(left_p.x(), bottom_p.y()),
			 Point(right_p.x(),top_p.y())));
}

/* Some applications might be more interested in the four points
 * that are from the input and that define the empty rectangle
 */
template<class T>
Quadruple<Largest_empty_iso_rectangle_2<T>::Point,
          Largest_empty_iso_rectangle_2<T>::Point,
          Largest_empty_iso_rectangle_2<T>::Point,
          Largest_empty_iso_rectangle_2<T>::Point>
Largest_empty_iso_rectangle_2<T>::get_left_bottom_right_top()
{
  if(x_sorted.size() == 4) {
    return(make_quadruple(bl_p, bl_p, tr_p, tr_p));
  }
  update();
  return(make_quadruple(left_p, bottom_p, right_p, top_p));
}


template<class T>
void
Largest_empty_iso_rectangle_2<T>::init(const Point& bl, const Point& tr)
{
  // determine extreme values of bounding box
  bl_p = bl;
  tr_p = tr;

  // add extreme points
  insert(Point(bl.x(),bl.y()), BOT_LEFT);
  insert(Point(tr.x(),bl.y()), BOT_RIGHT);
  insert(Point(bl.x(),tr.y()), TOP_LEFT);
  insert(Point(tr.x(),tr.y()), TOP_RIGHT);
}

// ctor
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(
  const Point& bl,
  const Point& tr)
  : cache_valid(false), _gt(),
    x_sorted(Less_xy(traits())),
    y_sorted(Less_yx(traits()))
{
  // precondition: bl and tr
  init(bl, tr);
}

// ctor
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(
  const Iso_rectangle_2 &b)
  : cache_valid(false), _gt(),
    x_sorted(Less_xy(traits())),
    y_sorted(Less_yx(traits()))
{
  init(b.min(), b.max());
}

// ctor
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2()
  : cache_valid(false), _gt(),
    x_sorted(Less_xy(traits())),
    y_sorted(Less_yx(traits()))
{
  Point bl(0,0);
  Point tr(1,1);

  init(bl,tr);
}


template<class T>
typename Largest_empty_iso_rectangle_2<T>::const_iterator 
Largest_empty_iso_rectangle_2<T>::begin() const
{
  typename Point_data_set_of_x::const_iterator i = x_sorted.begin();
  while(i != x_sorted.end() && (*i)->type != REG) 
    ++i;
 
  return const_iterator(i);
}

template<class T>
typename Largest_empty_iso_rectangle_2<T>::const_iterator 
Largest_empty_iso_rectangle_2<T>::end() const
{
   typename Point_data_set_of_x::const_iterator i = x_sorted.end();
   while(--i != x_sorted.begin() && (*i)->type != REG);
   if((*i)->type != REG)
     // The points list is actually empty. Point to end() to make
     // begin() == end()
     i =  x_sorted.end();
   else
     // increment i to make it point to a corner point
     ++i;

   return const_iterator(i);
}

template<class T>
void 
Largest_empty_iso_rectangle_2<T>::clear()
{
  cache_valid = false;
  for(typename Point_data_set_of_x::iterator iter = x_sorted.begin();
      iter != x_sorted.end();
      ++iter)
    delete(*iter);

  x_sorted.clear();
  y_sorted.clear();

  insert(Point(bl_p.x(),bl_p.y()),BOT_LEFT);
  insert(Point(tr_p.x(),bl_p.y()),BOT_RIGHT);
  insert(Point(bl_p.x(),tr_p.y()),TOP_LEFT);
  insert(Point(tr_p.x(),tr_p.y()),TOP_RIGHT);
}


template<class T>
void 
Largest_empty_iso_rectangle_2<T>::determine_first_two_iters(
  typename Point_data_set_of_y::iterator &iter1,
  typename Point_data_set_of_y::iterator &iter2,
  typename Point_data_set_of_y::iterator &iter3,
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

