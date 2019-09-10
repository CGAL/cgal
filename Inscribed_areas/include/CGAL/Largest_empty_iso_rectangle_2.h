// Copyright (c) 1999  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Eli Packer (algorithm), Andreas Fabri (cgal conformance)

#ifndef CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
#define CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H

#include <CGAL/license/Inscribed_areas.h>


/*! \file
 * The implementation of the Largest_empty_iso_rectangle_2<Traits> class.
 */

#include <iostream>
#include <set>
#include <list>

#include <CGAL/utility.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

namespace CGAL {

/*!
  Largest_empty_iso_rectangle_2 is a class that implements the largest
  empty rectangle algorithm based on the algorithm described in

  M. Orlowski. A new algorithm for the largest empty rectangle problem.
  Algorithmica, 5:65-73, 1990.

  The problem is the following. Given a set of points and a bounding box that
  include the points, find the largest axis parallel rectangle that contains
  no points inside it.

  The algorithm is extended to support the degenerate case in which two points
  have the same x or y coordinates.

  The algorithm checks all the empty rectangles that are bounded by either
  points or edges of the bounding box (other empty rectangles can be enlarged
  and remain empty). There are O(n^2) such rectangles. It is done in three
  phases. In the first one empty rectangles that are bounded by two opposite
  edges of the bounding box are checked. In the second one, other empty
  rectangles that are bounded by one or two edges of the bounding box are
  checked. In the third phase, all other empty rectangles, namely the ones
  that are bounded by four points, are checked.
*/

template<class T>
class Largest_empty_iso_rectangle_2 {
public:

  struct Internal_point;

  typedef typename T::FT                NT;
  typedef typename T::Point_2           Point_2;
  typedef Internal_point                Point;
  typedef typename T::Iso_rectangle_2   Iso_rectangle_2;
  typedef T                             Traits;

  class Point_data;

  class Less_yx
  {
  private:
    class Less_xy_internal_point;
    class Less_yx_internal_point;

    Traits _gt;

    const Traits & traits() const {return _gt;};

    Less_yx(){}

  public:
    Less_yx(const Traits& t)
      : _gt(t)
    {}

    bool operator()(const Point_data *a, const Point_data *b) const
    {
      Comparison_result c = traits().compare_y_2_object()
             (b->p.y_part, a->p.y_part);
      if(c == LARGER) {
        return true;
      } else if (c == EQUAL) {
        return traits().less_x_2_object()(a->p.x_part, b->p.x_part); 
      } 
      return false;
    }
  };

  class Less_xy
  {
  private:
    Traits _gt;

    const Traits & traits() const {return _gt;};

  public:
    Less_xy(const Traits& t)
      : _gt(t)
    {}

    bool operator()(const Point_data *a, const Point_data *b) const
    {
      Comparison_result c = traits().compare_x_2_object()
             (b->p.x_part, a->p.x_part);
      if(c == LARGER) {
        return true;
      } else if (c == EQUAL) {
        return traits().less_y_2_object()(a->p.y_part, b->p.y_part); 
      } 
      return false;
    }
  };

  template < class Node>
  struct Proj_point {
    typedef Node                  argument_type;
    typedef Point_2                 result_type;
    Point_2&       operator()( Node& x)       const { return x->p.original; }
    const Point_2& operator()( const Node& x) const { return x->p.original; }
  };

  typedef std::set<Point_data *,Less_xy> Point_data_set_of_x;
  typedef std::set<Point_data *,Less_yx> Point_data_set_of_y;

  // The following is an iterator adapter that allows us to enumerate 
  // the points in a set where they are stored
  typedef Iterator_project<typename Point_data_set_of_x::const_iterator, 
                           Proj_point<Point_data*>, 
                           const Point_2&,
                           const Point_2*> const_iterator;


  enum Point_type{REG, BOT_RIGHT, BOT_LEFT, TOP_LEFT, TOP_RIGHT};

  const Traits & traits() const {return _gt;};

  //! A constructor given two points parameters. The parameters are two
  //! opposite corners of the bounding box.
  Largest_empty_iso_rectangle_2(const Point_2& bl, const Point_2& tr);

  //! Constructor given an Iso Rectangle parameter. The parameter is
  //! the bounding box.
  Largest_empty_iso_rectangle_2(const Iso_rectangle_2 &b);

  //! A parameter-less constructor
  Largest_empty_iso_rectangle_2();

  //! Add a point to the data.
  bool
  insert(const Point_2& p);

  //! The STL standard member function for insertion.
  void
  push_back(const Point& _p)
  {
    insert(_p);
  }
  
  //! Insertion of an iterator range.
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

  //! Remove a point from data.
  bool
  remove(const Point& p);

  //! Get the bounding box of the instance.
  Iso_rectangle_2
  get_bounding_box();

  //! Retrieve largest empty iso rectangle.

  /*! get_largest_empty_iso_rectangle() retrieves the largest empty iso
   * rectangle which lies inside the bounding box of the instance. An
   * empty rectangle is defined as a rectangle that contains no points
   * inside its interior.
   * \return the largest empty iso rectangle.
   */
  Iso_rectangle_2 
  get_largest_empty_iso_rectangle();

  //! Retrieve four points from the input that define the largest
  //! empty iso rectangle.
  Quadruple<Point_2, Point_2, Point_2, Point_2>
  get_left_bottom_right_top()
  {
    if(x_sorted.size() == 4) {
      return(make_quadruple(bl_p.original, bl_p.original,
                            tr_p.original, tr_p.original));
    }
    update();
    return(make_quadruple(left_p.original, bottom_p.original, 
                          right_p.original, top_p.original));
  }

  //! Clear the data(remove the points).
  void 
  clear();

  //! Get a begin iterator to points
  const_iterator 
  begin()  const;

  //! Get a after-the-end iterator to points
  const_iterator
  end()  const;


  //! A destructor
  ~Largest_empty_iso_rectangle_2();

  //! An operator=
  Largest_empty_iso_rectangle_2<T>&
    operator =(const Largest_empty_iso_rectangle_2<T>& ler);

  //! A copy constructor
  Largest_empty_iso_rectangle_2<T>(
	       const Largest_empty_iso_rectangle_2<T>& ler);

  struct Internal_point {
    Point_2 x_part;// the x coordinate of the point
    Point_2 y_part;// the y coordinate of the point
    Point_2 original;

    Internal_point &
     operator=(const Internal_point &other)
    {
      x_part = other.x_part;
      y_part = other.y_part;
      original = other.original;

      return(*this);
    }

    Internal_point() // no real value - just to allow construction of LER
      : x_part(Point_2(0,0)), y_part(Point_2(0,0)), original(Point_2(0,0)) {}

    Internal_point(int x,int y)
      : x_part(Point_2(x,y)), y_part(Point_2(x,y)), original(Point_2(x,y)) {}

    Internal_point(const Point_2 &p)
       : x_part(p), y_part(p), original(p) {}
  };

  class Point_data {
  public:

    Point p;

    /*! the next two members save set of points that are needed
        for the third phase.
    */
    std::set<Point_data *,Less_yx> *right_tent;
    std::set<Point_data *,Less_yx> *left_tent;
    
    /* detemine whether the point is a bounding box corner
       (thus not implicitely inserted as a point, or not.
    */

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
      if(right_tent != NULL){
        delete right_tent;
      }
      if (left_tent != NULL){
        delete left_tent;
      }
    }
  }; 


private:

  /* this struct is the point representation. It is composed of two points
   * such that one holds the x coordinate and the other holds the y coordinate
   */

  /*! false if no points were inserted or removed, thus the previous
     results hold. Otherwise there is a need to find the new largest
     empty rectangle..
  */

  class Less_xy_internal_point
  {
  private:
    Traits _gt;

    const Traits & traits() const {return _gt;};

  public:
    Less_xy_internal_point(const Traits& t)
      : _gt(t)
    {}

    bool operator()(const Point &a, const Point &b) const
    {
      Comparison_result c = traits().compare_x_2_object()
             (b.x_part, a.x_part);
      if(c == LARGER) {
        return true;
      } else if (c == EQUAL) {
        return traits().less_y_2_object()(a.y_part, b.y_part); 
      } 
      return false;
    }
  };

  class Less_yx_internal_point
  {
  private:
    Traits _gt;

    const Traits & traits() const {return _gt;};

    Less_yx_internal_point(){}

  public:
    Less_yx_internal_point(const Traits& t)
      : _gt(t)
    {}

    bool operator()(const Point &a, const Point &b) const
    {
      Comparison_result c = traits().compare_y_2_object()
             (b.y_part, a.y_part);
      if(c == LARGER) {
        return true;
      } else if (c == EQUAL) {
        return traits().less_x_2_object()(a.x_part, b.x_part); 
      } 
      return false;
    }
  };

  bool cache_valid;

  Traits _gt;

  /*! this class holds points' data as needed in the LER process.
   */

  bool less_xy(const Point_data *a, const Point_data *b) const
  {
    return(!larger_xy(a,b));
  }

  bool less_yx(const Point_data *a, const Point_data *b) const
  {
    return(!larger_yx(a,b));
  }

  bool larger_xy(const Point_data *a, const Point_data *b) const
  {
    Comparison_result c = traits().compare_x_2_object()
           (a->p.x_part, b->p.x_part);
    if(c == LARGER) {
      return true;
    } else if (c == EQUAL) {
      return traits().less_y_2_object()(b->p.y_part, a->p.y_part); 
    } 
    return false;
  }

  bool larger_yx(const Point_data *a, const Point_data *b) const
  {
    Comparison_result c = traits().compare_y_2_object()
             (a->p.y_part, b->p.y_part);
    if(c == LARGER) {
      return true;
    } else if (c == EQUAL) {
      return traits().less_x_2_object()(b->p.x_part, a->p.x_part); 
    } 
    return false;
  }

  // the next sets store the points sorted
  Point_data_set_of_x x_sorted;
  Point_data_set_of_y y_sorted;

  Less_xy_internal_point less_xy_point;
  Less_yx_internal_point less_yx_point;  

  // bottom left and top right points of the bounding box
  Point bl_p, tr_p;

  // the bounding box of the points
  Iso_rectangle_2 bbox_p;

  // the points that define the largest empty iso-rectangle
  Point left_p, bottom_p, right_p ,top_p; 

  // save the largest empty rectangle size found by now
  NT largest_rect_size;

  // insert points
  bool insert(const Point& _p,Point_type i_type);
  bool insert(const Point_2& _p,Point_type i_type);
  
  /* the phases of the algorithm as described in the paper
  */

  void phase_1();
  // the first phase is divided to work on the x-axis and work on the y-axis
  void phase_1_on_x();
  void phase_1_on_y();
  // the second phase is divided to four parts,
  // one for each edge of the bounding box
  void phase_2_on_bot();
  void phase_2_on_top();
  void phase_2_on_left();
  void phase_2_on_right();
  void phase_2();
  void phase_3();

  // the next functions are used by the functions of the three phases 
  void check_for_larger(const Point& px0, 
			const Point& py0, 
			const Point& px1,
			const Point& py1);
  void tent(Point_data *first, Point_data *second);
  void tent(Point_data *first, Point_data *second, Point_data *third);

  void get_next_for_top(typename std::list<Point_data *>::iterator &iter,
			typename std::list<Point_data *>::iterator &beyond)
  {
    while(iter != beyond && ((*iter)->type == BOT_RIGHT 
	  		   || (*iter)->type == BOT_LEFT))
      ++iter;
  }
 

  void get_prev_for_top(typename std::list<Point_data *>::iterator &iter)
  {
    while((*iter)->type == BOT_RIGHT || (*iter)->type == BOT_LEFT)
      --iter;
  }

  void get_next_for_bot(typename std::list<Point_data *>::iterator &iter,
			typename std::list<Point_data *>::iterator &beyond)
  {
    while(iter != beyond && ((*iter)->type == TOP_LEFT 
			     || (*iter)->type == TOP_RIGHT))
      ++iter;
  }

  void get_prev_for_bot(typename std::list<Point_data *>::iterator &iter)
  {
    while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
      --iter;
  }

  void get_next_for_bot(typename Point_data_set_of_y::iterator &iter)
  {
    while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
      ++iter;
  }
 
  void get_next_for_bot(typename Point_data_set_of_y::iterator &iter,
			typename Point_data_set_of_y::iterator &last)
  {
   while(iter != last && ((*iter)->type == TOP_LEFT 
			  || (*iter)->type == TOP_RIGHT))
     ++iter;
  }

  void get_prev_for_bot(typename Point_data_set_of_y::iterator &iter)
  {
    while((*iter)->type == TOP_LEFT || (*iter)->type == TOP_RIGHT)
      --iter;
  }


  void get_next_for_left(typename std::list<Point_data *>::iterator &iter,
			 typename std::list<Point_data *>::iterator &beyond)
  {
    while(iter != beyond && ((*iter)->type == BOT_RIGHT 
			     || (*iter)->type == TOP_RIGHT))
      ++iter;
  }

  void get_prev_for_left(typename std::list<Point_data *>::iterator &iter)
  {
    while((*iter)->type == BOT_RIGHT || (*iter)->type == TOP_RIGHT)
      --iter;
  }

  void get_next_for_right(typename std::list<Point_data *>::iterator &iter,
			  typename std::list<Point_data *>::iterator &beyond)
  {
    while(iter != beyond && ((*iter)->type == BOT_LEFT 
			   || (*iter)->type == TOP_LEFT))
      ++iter;
  }


  void get_prev_for_right(typename std::list<Point_data *>::iterator &iter)
  {
    while((*iter)->type == BOT_LEFT || (*iter)->type == TOP_LEFT)
      --iter;
  }

  void determine_first_two_iters(typename Point_data_set_of_y::iterator& iter1,
				 typename Point_data_set_of_y::iterator& iter2,
				 typename Point_data_set_of_y::iterator& iter3,
				 bool& first_iter_is_right,
				 bool& second_iter_is_right,
				 bool& third_iter_is_right)
  {
    if (first_iter_is_right) {
      if (second_iter_is_right) {
        iter1 = iter2;
        iter2 = iter3;
        first_iter_is_right = second_iter_is_right;
        second_iter_is_right = third_iter_is_right;
      } else {
        if (third_iter_is_right) {
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
      if (second_iter_is_right) {
        if (third_iter_is_right) {
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

  void determine_next_iter(
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

  void calls_for_tents(typename Point_data_set_of_y::iterator iter1,
		       typename Point_data_set_of_y::iterator iter2)
  {
    if(less_xy(*iter1, *iter2))
      tent(*iter1,*iter2);
    else
      tent(*iter2,*iter1);
  }


  void calls_for_tents(typename Point_data_set_of_y::iterator iter1,
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

  void phase_2_update_y_sorted_list();
  void phase_3_check_for_larger(typename Point_data_set_of_y::iterator iter,
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
  
  void empty_tents();

  // call the computation of the largest empty rectangle .
  void update();
  // init class.
  void init(const Point_2& bl, const Point_2& tr);
  void copy_memory(const Largest_empty_iso_rectangle_2<T>& ler);
  void free_memory();

  // add a point to data
  bool
  insert(const Point& p);

  // Auxiliary iterators for convenience
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
  bbox_p = ler.bbox_p;
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
  y_sorted(Less_yx(traits())),
  less_xy_point(traits()),
  less_yx_point(traits())
{
  copy_memory(ler);
}

template<class T>
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point_2& _p)
{
  // check that the point is inside the bounding box 
  if(bbox_p.has_on_unbounded_side(_p)) {
    return(false);
  }

  return(insert(Point(_p)));
}

template<class T>
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point& _p)
{
  // check that the point is not already inserted
  Point_data po(_p);
  typename Point_data_set_of_x::iterator iter = x_sorted.find(&po);

  if(iter != x_sorted.end())
    return(false);

  cache_valid = false;
  Point_data_set_of_y *right_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data_set_of_y *left_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data * ppo = new Point_data(_p,right_tent,left_tent,REG);

  x_sorted.insert(ppo);
  y_sorted.insert(ppo);
  return(true);
}

template<class T>
bool
Largest_empty_iso_rectangle_2<T>::remove(const Point& _p)
{
  cache_valid = false;
  Point_data po(_p);
  typename Point_data_set_of_x::iterator iter1 = x_sorted.find(&po);
  typename Point_data_set_of_y::iterator iter2 = y_sorted.find(&po);

  // point does not exist or a corner point
  if(iter1 == x_sorted.end() || (*iter1)->type != REG)
    return(false);

  Point_data* ptr  = *iter1;
  x_sorted.erase(iter1);
  y_sorted.erase(iter2);
  delete ptr;

  return(true);
}

template<class T>
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point_2& _p,
					 Point_type i_type)
{
  // check that the point is inside the bounding box 
  if((i_type == REG) && bbox_p.has_on_unbounded_side(_p)) {
    return false;
  }

  return(insert(Point(_p),i_type));
}


template<class T>
bool
Largest_empty_iso_rectangle_2<T>::insert(const Point& _p,
					 Point_type i_type)
{
  // check that the point is not already inserted
  Point_data po(_p);
  typename Point_data_set_of_x::iterator iter = x_sorted.find(&po);

  if(iter != x_sorted.end())
    return(false);

  cache_valid = false;
  Point_data_set_of_y *right_tent =
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data_set_of_y *left_tent = 
    new Point_data_set_of_y(Less_yx(traits()));
  Point_data *ppo = new Point_data(_p,right_tent,left_tent,i_type);

  x_sorted.insert(ppo);
  y_sorted.insert(ppo);
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
  Iso_rectangle_2 rec(less_xy_point(px0,px1) ? px0.x_part : px1.x_part,
                      less_xy_point(px0,px1) ? px1.x_part : px0.x_part,
                      less_yx_point(py0,py1) ? py0.y_part : py1.y_part,
                      less_yx_point(py0,py1) ? py1.y_part : py0.y_part);
  NT rect_size = rec.area();

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
    size = static_cast<int>(Point_data_list.size());

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
      if(iter1 != first_iter) {
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
    size = static_cast<int>(Point_data_list.size());

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
      if(iter1 != first_iter) {
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
  int points_removed = 0,
    size = static_cast<int>(Point_data_list.size());

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
      if(iter1 != first_iter) {
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
  int points_removed = 0,size = static_cast<int>(Point_data_list.size());

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
Largest_empty_iso_rectangle_2<T>::phase_3()
{
  // init for compiler warning
  bool first_iter_is_right(true);
  bool second_iter_is_right(true);
  bool third_iter_is_right(true);
  bool first_exist(true);
  bool second_exist(true);
  bool third_exist(true);
  
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
  return bbox_p;
}

/* Performs the computation if the cache is invalid.
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

  return(Iso_rectangle_2(
           less_xy_point(left_p.x_part,right_p.x_part) ?
                 left_p.x_part : right_p.x_part,
           less_xy_point(left_p.x_part,right_p.x_part) ?
                 right_p.x_part : left_p.x_part,
           less_yx_point(bottom_p.x_part,top_p.x_part) ?
                 bottom_p.y_part : top_p.y_part,
           less_yx_point(bottom_p.x_part,top_p.x_part) ?
                 top_p.y_part : bottom_p.y_part));
}

template<class T>
void
Largest_empty_iso_rectangle_2<T>::init(const Point_2& bl, const Point_2& tr)
{
  // determine extreme values of bounding box
  bbox_p = Iso_rectangle_2(bl,tr);
  // add extreme points
  
  insert(bbox_p.vertex(0), BOT_LEFT);
  insert(bbox_p.vertex(1), BOT_RIGHT);
  insert(bbox_p.vertex(3), TOP_LEFT);
  insert(bbox_p.vertex(2), TOP_RIGHT);
}

// ctor
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2(
  const Point_2& bl,
  const Point_2& tr)
  : cache_valid(false), _gt(),
     x_sorted(Less_xy(traits())),
     y_sorted(Less_yx(traits())),
     less_xy_point(traits()),
     less_yx_point(traits()),
     bl_p(bl),
     tr_p(tr)
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
     y_sorted(Less_yx(traits())),
     less_xy_point(traits()),
     less_yx_point(traits()),
     bl_p((b.min)()),
     tr_p((b.max)())
{
  init((b.min)(), (b.max)());
}

// ctor
template<class T>
Largest_empty_iso_rectangle_2<T>::Largest_empty_iso_rectangle_2()
  : cache_valid(false), _gt(),
    x_sorted(Less_xy(traits())),
    y_sorted(Less_yx(traits())),
    less_xy_point(traits()),
    less_yx_point(traits())
{
  Point bl(0,0);
  Point tr(1,1);

  init(bl.x_part,tr.x_part);
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
   while(--i != x_sorted.begin() && (*i)->type != REG) {}
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

  // add extreme points
  insert(bbox_p.vertex(0), BOT_LEFT);
  insert(bbox_p.vertex(1), BOT_RIGHT);
  insert(bbox_p.vertex(3), TOP_LEFT);
  insert(bbox_p.vertex(2), TOP_RIGHT);
}



} //namespace CGAL

#endif // CGAL_LARGEST_EMPTY_ISO_RECTANGLE_2_H
