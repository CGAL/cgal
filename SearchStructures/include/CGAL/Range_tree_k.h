// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Gabriele Neyer

#ifndef CGAL_RANGE_TREE_K_H
#define CGAL_RANGE_TREE_K_H

#include <CGAL/license/SearchStructures.h>


// Predefined k-dimensional Range Trees (k=1..4) 
// The trees can either be templated with d arbitrary types
// (e.g., Range_tree_3) 
// or with an unary type for each dimension
// (e.g., Range_tree_uni_4).
// The container class and sequence container class as well as the 
// interfaces are defined in these classes.

#include <iostream>
#include <iterator>
#include <CGAL/Tree_base.h>
#include <CGAL/Tree_traits.h>
#include <CGAL/Range_tree_d.h>

//-------------------------------------------------------------------
// A one dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

namespace CGAL {


template <class C_Traits_1>
class Range_tree_1
{ 

public:
  typedef C_Traits_1 Traits;
  typedef typename C_Traits_1::Key Key;
  typedef typename C_Traits_1::Interval Interval;
  typedef typename C_Traits_1::Key_1 Key_1;
  typedef typename C_Traits_1::key_1 key_1;
  typedef typename C_Traits_1::low_1 low_1;
  typedef typename C_Traits_1::high_1 high_1;
  typedef typename C_Traits_1::compare_1 compare_1;


  typedef tree_point_traits<Key, Interval, Key_1, 
                            key_1, low_1, high_1, compare_1> I1;

  typedef Tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_1_type;
  Range_tree_1_type * range_tree_1;


  Range_tree_1()
  {
    anchor = new Tree_anchor_type;
    range_tree_1 = new Range_tree_1_type(*anchor);
  }

  template <class T>
  Range_tree_1(const T& first, 
	       const T& last)  
    : anchor(new Tree_anchor_type), range_tree_1(new Range_tree_1_type(*anchor))
  {
   range_tree_1->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete range_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    range_tree_1 = new Range_tree_1_type(*anchor);
    return range_tree_1->make_tree(first,last);
  }

  template <class T>
  T  window_query(Interval const &win,  
		  const T& result)
  {
    return range_tree_1->window_query(win, result);
  }

  ~Range_tree_1()
  {
    if (range_tree_1!=0)
      delete range_tree_1;
    if (anchor!=0)
      delete anchor;
  }

};


//-------------------------------------------------------------------
// A two dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

template <class C_Traits_2>
class Range_tree_2
{ 

public:
  typedef C_Traits_2 Traits;
  typedef typename C_Traits_2::Key Key;
  typedef typename C_Traits_2::Interval Interval;
  typedef typename C_Traits_2::Key_2 Key_2;
  typedef typename C_Traits_2::Key_1 Key_1;
  typedef typename C_Traits_2::key_1 key_1;
  typedef typename C_Traits_2::key_2 key_2;
  typedef typename C_Traits_2::low_1 low_1;
  typedef typename C_Traits_2::high_1 high_1;
  typedef typename C_Traits_2::low_2 low_2;
  typedef typename C_Traits_2::high_2 high_2;
  typedef typename C_Traits_2::compare_1 compare_1;
  typedef typename C_Traits_2::compare_2 compare_2;


  typedef tree_point_traits<Key, Interval, 
                            Key_1, key_1, low_1, high_1, compare_1> I1;

  typedef tree_point_traits<Key, Interval, 
                            Key_2, key_2, low_2, high_2, compare_2> I2;


  typedef Tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_1_type;
  Range_tree_1_type * range_tree_1;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_2_type;
  Range_tree_2_type *range_tree_2;


  Range_tree_2()
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1))
  {}

  template <class T>
  Range_tree_2(const T& first, 
	       const T& last)
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1))
  {
    range_tree_2->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete range_tree_2;
    delete range_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    range_tree_1 = new Range_tree_1_type(*anchor);
    range_tree_2 = new Range_tree_2_type(*range_tree_1);
    return range_tree_2->make_tree(first,last);
  }
  
  template <class T>
  T window_query(Interval const &win,  
		 const T& result)
  {
    return range_tree_2->window_query(win, result);
  }

  ~Range_tree_2()
  {
    if (range_tree_2!=0)
      delete range_tree_2;
    if (range_tree_1!=0)
      delete range_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};

//-------------------------------------------------------------------
// A three dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.
template <class C_Traits_3>
class Range_tree_3
{ 
public:
  typedef C_Traits_3 Traits;
  typedef typename C_Traits_3::Key Key;
  typedef typename C_Traits_3::Interval Interval;
  typedef typename C_Traits_3::Key_1 Key_1;
  typedef typename C_Traits_3::Key_2 Key_2;
  typedef typename C_Traits_3::Key_3 Key_3;
  typedef typename C_Traits_3::key_1 key_1;
  typedef typename C_Traits_3::key_2 key_2;
  typedef typename C_Traits_3::key_3 key_3;
  typedef typename C_Traits_3::low_1 low_1;
  typedef typename C_Traits_3::high_1 high_1;
  typedef typename C_Traits_3::low_2 low_2;
  typedef typename C_Traits_3::high_2 high_2;
  typedef typename C_Traits_3::low_3 low_3;
  typedef typename C_Traits_3::high_3 high_3;
  typedef typename C_Traits_3::compare_1 compare_1;
  typedef typename C_Traits_3::compare_2 compare_2;
  typedef typename C_Traits_3::compare_3 compare_3;

  typedef tree_point_traits<Key, Interval, Key_1,
                            key_1, low_1, high_1, compare_1> I1;

  typedef tree_point_traits<Key, Interval, Key_2,
                            key_2, low_2, high_2, compare_2> I2;

  typedef tree_point_traits<Key, Interval, Key_3, 
                            key_3, low_3, high_3, compare_3> I3;

  typedef Tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Range_tree_d<Key, Interval, I3> Range_tree_1_type;
  Range_tree_1_type * range_tree_1;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_2_type;
  Range_tree_2_type * range_tree_2;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_3_type;
  Range_tree_3_type *range_tree_3;

  Range_tree_3()
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1)),
      range_tree_3(new Range_tree_3_type(*range_tree_2))
  {}
  
  template <class T>
  Range_tree_3(const T& first, 
	       const T& last)
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1)),
      range_tree_3(new Range_tree_3_type(*range_tree_2))
  {
    range_tree_3->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete range_tree_3;
    delete range_tree_2;
    delete range_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    range_tree_1 = new Range_tree_1_type(*anchor);
    range_tree_2 = new Range_tree_2_type(*range_tree_1);
    range_tree_3 = new Range_tree_3_type(*range_tree_2);
    return range_tree_3->make_tree(first,last);
  }
  
  template <class T>
  T  window_query(Interval const &win,  
		  const T& result)
  {
    return range_tree_3->window_query(win, result);
  }

  ~Range_tree_3()
  {
    if (range_tree_3!=0)
      delete range_tree_3;
    if (range_tree_2!=0)
      delete range_tree_2;
    if (range_tree_1!=0)
      delete range_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};


//-------------------------------------------------------------------
// A three dimensional Range Tree is defined in this class.
// Ti is the type of each dimension of the tree.

template <class C_Traits_4>
class Range_tree_4
{ 
public:
  typedef C_Traits_4 Traits;
  typedef typename C_Traits_4::Key Key;
  typedef typename C_Traits_4::Interval Interval;
  typedef typename C_Traits_4::Key_1 Key_1;
  typedef typename C_Traits_4::Key_2 Key_2;
  typedef typename C_Traits_4::Key_3 Key_3;
  typedef typename C_Traits_4::Key_4 Key_4;
  typedef typename C_Traits_4::key_1 key_1;
  typedef typename C_Traits_4::key_2 key_2;
  typedef typename C_Traits_4::key_4 key_4;
  typedef typename C_Traits_4::key_3 key_3;
  typedef typename C_Traits_4::low_1 low_1;
  typedef typename C_Traits_4::high_1 high_1;
  typedef typename C_Traits_4::low_2 low_2;
  typedef typename C_Traits_4::high_2 high_2;
  typedef typename C_Traits_4::low_3 low_3;
  typedef typename C_Traits_4::high_3 high_3;
  typedef typename C_Traits_4::low_4 low_4;
  typedef typename C_Traits_4::high_4 high_4;
  typedef typename C_Traits_4::compare_1 compare_1;
  typedef typename C_Traits_4::compare_2 compare_2;
  typedef typename C_Traits_4::compare_3 compare_3;
  typedef typename C_Traits_4::compare_4 compare_4;

  typedef tree_point_traits<Key, Interval, Key_1,
                            key_1, low_1, high_1, compare_1> I1;  

  typedef tree_point_traits<Key, Interval, Key_2,
                            key_2, low_2, high_2, compare_2> I2;  

  typedef tree_point_traits<Key, Interval, Key_3,
                            key_3, low_3, high_3, compare_3> I3;  

  typedef tree_point_traits<Key, Interval, Key_4, 
                            key_4, low_4, high_4, compare_4> I4;  

  typedef Tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Range_tree_d<Key, Interval, I4> Range_tree_1_type;
  Range_tree_1_type * range_tree_1;

  typedef Range_tree_d<Key, Interval, I3> Range_tree_2_type;
  Range_tree_2_type * range_tree_2;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_3_type;
  Range_tree_3_type *range_tree_3;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_4_type;
  Range_tree_4_type *range_tree_4;

  Range_tree_4()
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1)),
      range_tree_3(new Range_tree_3_type(*range_tree_2)),
      range_tree_4(new Range_tree_4_type(*range_tree_3))
  {}

  template <class T>
  Range_tree_4(const T& first, 
	       const T& last)
    : anchor(new Tree_anchor_type),
      range_tree_1(new Range_tree_1_type(*anchor)),
      range_tree_2(new Range_tree_2_type(*range_tree_1)),
      range_tree_3(new Range_tree_3_type(*range_tree_2)),
      range_tree_4(new Range_tree_4_type(*range_tree_3))
  {
    range_tree_4->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete range_tree_4;
    delete range_tree_3;
    delete range_tree_2;
    delete range_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    range_tree_1 = new Range_tree_1_type(*anchor);
    range_tree_2 = new Range_tree_2_type(*range_tree_1);
    range_tree_3 = new Range_tree_3_type(*range_tree_2);
    range_tree_4 = new Range_tree_4_type(*range_tree_3);
    return range_tree_4->make_tree(first,last);
  }

  template <class T>
  T  window_query(Interval const &win,  
		  const T& result)
  {
    return range_tree_4->window_query(win, result);
  }

  ~Range_tree_4()
  {
    if (range_tree_4!=0)
      delete range_tree_4;
    if (range_tree_3!=0)
      delete range_tree_3;
    if (range_tree_2!=0)
      delete range_tree_2;
    if (range_tree_1!=0)
      delete range_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};

} //namespace CGAL

#endif
