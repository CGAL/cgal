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

#ifndef CGAL_SEGMENT_TREE_K_H
#define CGAL_SEGMENT_TREE_K_H

#include <CGAL/license/SearchStructures.h>


// Predefined k-dimensional Segment Trees (k=1..4) 
// The trees can either be templated with d arbitrary types
// (e.g., Segment_tree_3) 
// or with an unary type for each dimension
// (e.g., Segment_tree_uni_4).
// The container class and sequence container class as well as the 
// data accessors are defined in these classes.

#include <iostream>
#include <iterator>
#include <list>
#include <CGAL/Tree_base.h>
#include <CGAL/Tree_traits.h>
#include <CGAL/Segment_tree_d.h>

namespace CGAL {

template <class C_Traits_1>
class Segment_tree_1
{ 
public:
  typedef  C_Traits_1 Traits;
  typedef typename C_Traits_1::Key Key;
  typedef typename C_Traits_1::Interval Interval;
  typedef typename C_Traits_1::Key_1 Key_1;
  typedef typename C_Traits_1::key_1 key_1;
  typedef typename C_Traits_1::low_1 low_1;
  typedef typename C_Traits_1::high_1 high_1;
  typedef typename C_Traits_1::compare_1 compare_1;

  typedef tree_interval_traits<Interval, Interval, Key_1,  
                               low_1, high_1, low_1, 
                               high_1, compare_1> I1; 


  typedef Tree_anchor<Interval, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Segment_tree_d<Interval, Interval, I1> Segment_tree_1_type;
  Segment_tree_1_type * segment_tree_1;


  Segment_tree_1()
    : anchor(new Tree_anchor_type), segment_tree_1(new Segment_tree_1_type(*anchor))
  {}
  
  template <class T>
  Segment_tree_1(const T& first, 
		 const T& last) 
    : anchor(new Tree_anchor_type), segment_tree_1(new Segment_tree_1_type(*anchor))
 {
   segment_tree_1->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete segment_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    segment_tree_1 = new Segment_tree_1_type(*anchor);
    return segment_tree_1->make_tree(first,last);
  }

  template <class T>  
  T  window_query(Interval const &win, const T& result)
  {
    return segment_tree_1->window_query(win, result);
  }

  template <class T>  
  T  enclosing_query(Interval const &win, const T& result)
  {
    return segment_tree_1->enclosing_query(win, result);
  }

  ~Segment_tree_1()
  {
    if (segment_tree_1!=0)
      delete segment_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};


//-------------------------------------------------------------------
// A two dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.

template <class C_Traits_2>
class Segment_tree_2
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

  typedef typename std::list<Interval>::iterator l_iterator;
  typedef typename std::vector<Interval>::iterator v_iterator;

  typedef tree_interval_traits<Interval, Interval, Key_1,  
                               low_1,  high_1, 
                               low_1,  high_1, compare_1> I1; 

  typedef tree_interval_traits<Interval, Interval, Key_2,  
                               low_2,  high_2, 
                               low_2,  high_2, compare_2> I2; 

  typedef Tree_anchor<Interval, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Segment_tree_d<Interval, Interval, I2> Segment_tree_1_type;
  Segment_tree_1_type * segment_tree_1;

  typedef Segment_tree_d<Interval, Interval, I1> Segment_tree_2_type;
  Segment_tree_2_type *segment_tree_2;

  Segment_tree_2()
    : anchor( new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1))
  {}

  template <class T>
  Segment_tree_2(const T& first, 
		 const T& last)
    : anchor( new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1))  
  {
    segment_tree_2->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete segment_tree_2;
    delete segment_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    segment_tree_1 = new Segment_tree_1_type(*anchor);
    segment_tree_2 = new Segment_tree_2_type(*segment_tree_1);
    return segment_tree_2->make_tree(first,last);
  }
  
  template <class T>
  T  window_query(Interval const &win, const T& result)
  {
    return segment_tree_2->window_query(win, result);
  }

  template <class T>
  T enclosing_query(Interval const &win, const T& result)
  {
    return segment_tree_2->enclosing_query(win, result);
  }

  ~Segment_tree_2()
  {
    if (segment_tree_2!=0)
      delete segment_tree_2;
    if (segment_tree_1!=0)
      delete segment_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};


//-------------------------------------------------------------------
// A three dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.
template <class C_Traits_3>
class Segment_tree_3
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

  typedef typename std::list<Interval>::iterator l_iterator;
  typedef typename std::vector<Interval>::iterator v_iterator;

  typedef tree_interval_traits<Interval, Interval, Key_1,  
                               low_1, high_1, low_1, high_1,  compare_1> I1; 

  typedef tree_interval_traits<Interval, Interval, Key_2,  
                               low_2, high_2, low_2, high_2,  compare_2> I2;


  typedef tree_interval_traits<Interval, Interval, Key_3, 
                               low_3, high_3, low_3, high_3,  compare_3> I3;


  typedef Tree_anchor<Interval, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Segment_tree_d<Interval, Interval, I3> Segment_tree_1_type;
  Segment_tree_1_type * segment_tree_1;

  typedef Segment_tree_d<Interval, Interval, I2> Segment_tree_2_type;
  Segment_tree_2_type *segment_tree_2;

  typedef Segment_tree_d<Interval, Interval, I1> Segment_tree_3_type;
  Segment_tree_3_type *segment_tree_3;

  Segment_tree_3()
    : anchor(new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1)),
      segment_tree_3(new Segment_tree_3_type(*segment_tree_2))
  {}

  template <class T>
  Segment_tree_3(const T& first, 
		 const T& last)
    : anchor(new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1)),
      segment_tree_3(new Segment_tree_3_type(*segment_tree_2))
  {
    segment_tree_3->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete segment_tree_3;
    delete segment_tree_2;
    delete segment_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    segment_tree_1 = new Segment_tree_1_type(*anchor);
    segment_tree_2 = new Segment_tree_2_type(*segment_tree_1);
    segment_tree_3 = new Segment_tree_3_type(*segment_tree_2);
    return segment_tree_3->make_tree(first,last);
  }

  template <class T>  
  T window_query(Interval const &win, const T& result)
  {
    return (*segment_tree_3).window_query(win, result);
  }

  template <class T>  
  T  enclosing_query(Interval const &win, const T& result)
  {
    return (*segment_tree_3).enclosing_query(win, result);
  }

  ~Segment_tree_3()
  {
    if (segment_tree_3!=0)
      delete segment_tree_3;
    if (segment_tree_2!=0)
      delete segment_tree_2;
    if (segment_tree_1!=0)
      delete segment_tree_1;
    if (anchor!=0)
      delete anchor;
  }
};


//-------------------------------------------------------------------
// A three dimensional Segment Tree is defined in this class.
// Ti is the type of each dimension of the tree.
template <class C_Traits_4>
class Segment_tree_4
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

  typedef typename std::list<Interval>::iterator l_iterator;
  typedef typename std::vector<Interval>::iterator v_iterator;

  typedef tree_interval_traits<Interval, Interval, Key_1, 
                               low_1, high_1, low_1, high_1, compare_1> I1;

  typedef tree_interval_traits<Interval, Interval, Key_2,  
                               low_2, high_2, low_2, high_2, compare_2> I2; 

  typedef tree_interval_traits<Interval, Interval, Key_3, 
                               low_3, high_3, low_3, high_3, compare_3> I3; 

  typedef tree_interval_traits<Interval, Interval, Key_4,  
                               low_4, high_4, low_4, high_4, compare_4> I4; 


  typedef Tree_anchor<Interval, Interval> Tree_anchor_type;
  Tree_anchor_type *anchor;

  typedef Segment_tree_d<Interval, Interval, I4> Segment_tree_1_type;
  Segment_tree_1_type * segment_tree_1;

  typedef Segment_tree_d<Interval, Interval, I3> Segment_tree_2_type;
  Segment_tree_2_type *segment_tree_2;

  typedef Segment_tree_d<Interval, Interval, I2> Segment_tree_3_type;
  Segment_tree_3_type *segment_tree_3;

  typedef Segment_tree_d<Interval, Interval, I1> Segment_tree_4_type;
  Segment_tree_4_type *segment_tree_4;

  Segment_tree_4()
    : anchor(new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1)),
      segment_tree_3(new Segment_tree_3_type(*segment_tree_2)),
      segment_tree_4(new Segment_tree_4_type(*segment_tree_3))
  {}

  template <class T>
  Segment_tree_4(const T& first, 
		 const T& last)
    : anchor(new Tree_anchor_type),
      segment_tree_1(new Segment_tree_1_type(*anchor)),
      segment_tree_2(new Segment_tree_2_type(*segment_tree_1)),
      segment_tree_3(new Segment_tree_3_type(*segment_tree_2)),
      segment_tree_4(new Segment_tree_4_type(*segment_tree_3))
  {
   segment_tree_4->make_tree(first,last);
  }

  template <class T>
  bool make_tree(const T& first, 
		 const T& last)
  {
    delete segment_tree_4;
    delete segment_tree_3;
    delete segment_tree_2;
    delete segment_tree_1;
    delete anchor;
    anchor = new Tree_anchor_type;
    segment_tree_1 = new Segment_tree_1_type(*anchor);
    segment_tree_2 = new Segment_tree_2_type(*segment_tree_1);
    segment_tree_3 = new Segment_tree_3_type(*segment_tree_2);
    segment_tree_4 = new Segment_tree_4_type(*segment_tree_3);
    return segment_tree_4->make_tree(first,last);
  }

  template <class T>
  T window_query(Interval const &win, const T& result)
  {
    return (*segment_tree_4).window_query(win, result);
  }

  template <class T>
  T  enclosing_query(Interval const &win, const T& result)
  {
    return (*segment_tree_4).enclosing_query(win, result);
  }

  ~Segment_tree_4()
  {
    if (segment_tree_4!=0)
      delete segment_tree_4;

    if (segment_tree_3!=0)
      delete segment_tree_3;

    if (segment_tree_2!=0)
      delete segment_tree_2;

    if (segment_tree_1!=0)
      delete segment_tree_1;

    if (anchor!=0)
      delete anchor;
  }
};

} //namespace CGAL

#endif
