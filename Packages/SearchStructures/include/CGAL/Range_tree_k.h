// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 2000, August 09
//
// file          : include/CGAL/Range_tree_k.h
// package       : SearchStructures (2.54)
// maintainer    : Philipp Kramer <kramer@inf.ethz.ch>
// source        : include/CGAL/Range_tree_k.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Gabriele Neyer
//
// coordinator   : Peter Widmayer, ETH Zurich
//
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#ifndef __CGAL_Range_tree_k__
#define __CGAL_Range_tree_k__

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

CGAL_BEGIN_NAMESPACE


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


  typedef tree_point_traits<Key, Interval, Key_1, key_1, low_1, 
  high_1, compare_1> I1;

  typedef tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *Tree_anchor;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_1_type;
  Range_tree_1_type * CRange_tree_1;


  Range_tree_1()
  {
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
  }

  template <class T>
  Range_tree_1(T& first, 
	       T& last)  {
   Tree_anchor = new Tree_anchor_type;
   CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
   (*CRange_tree_1).make_tree(first,last);
  }

  template <class T>
  bool make_tree(T& first, 
		 T& last)
  {
    delete CRange_tree_1;
    delete Tree_anchor;
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    return (*CRange_tree_1).make_tree(first,last);
  }

  template <class T>
  T  window_query(Interval const &win,  
		  T result)
  {
    return (*CRange_tree_1).window_query(win, result);
  }

  ~Range_tree_1()
  {
    if (CRange_tree_1!=0)  //(Range_tree_1_type *)
      delete CRange_tree_1;
    CRange_tree_1=0;  //(Range_tree_1_type *)
    if (Tree_anchor!=0) //(Tree_anchor_type *)
      delete Tree_anchor;
    Tree_anchor=0; //(Tree_anchor_type *)
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
  Key_1, key_1, low_1, high_1, 
  compare_1> I1;

  typedef tree_point_traits<Key, Interval, 
  Key_2,   key_2, low_2, high_2, 
  compare_2> I2;


  typedef tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *Tree_anchor;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_1_type;
  Range_tree_1_type * CRange_tree_1;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_2_type;
  Range_tree_2_type *CRange_tree_2;


  Range_tree_2()
  {
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
  }

  template <class T>
  Range_tree_2(T& first, 
	       T& last)  {
   Tree_anchor = new Tree_anchor_type;
   CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
   CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
   (*CRange_tree_2).make_tree(first,last);
  }

  template <class T>
  bool make_tree(T& first, 
		 T& last)
  {
    delete CRange_tree_2;
    delete CRange_tree_1;
    delete Tree_anchor;
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
    return (*CRange_tree_2).make_tree(first,last);
  }
  
  template <class T>
  T window_query(Interval const &win,  
		 T result)
  {
    return (*CRange_tree_2).window_query(win, result);
  }

  ~Range_tree_2()
  {
    if (CRange_tree_2!=0) //(Range_tree_2_type *)
      delete CRange_tree_2;
    CRange_tree_2=0; //(Range_tree_2_type *)
    if (CRange_tree_1!=0) //(Range_tree_2_type *)
      delete CRange_tree_1;
    CRange_tree_1=0;
    if (Tree_anchor!=0)
      delete Tree_anchor;
    Tree_anchor=0;
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

  typedef tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *Tree_anchor;

  typedef Range_tree_d<Key, Interval, I3> Range_tree_1_type;
  Range_tree_1_type * CRange_tree_1;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_2_type;
  Range_tree_2_type * CRange_tree_2;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_3_type;
  Range_tree_3_type *CRange_tree_3;

  Range_tree_3()
  {
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
    CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
  }
  template <class T>
  Range_tree_3(T& first, 
	       T& last)  {
   Tree_anchor = new Tree_anchor_type;
   CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
   CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
   CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
   (*CRange_tree_3).make_tree(first,last);
  }

  template <class T>
  bool make_tree(T& first, 
		 T& last)
  {
    delete CRange_tree_3;
    delete CRange_tree_2;
    delete CRange_tree_1;
    delete Tree_anchor;
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
    CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
    return (*CRange_tree_3).make_tree(first,last);
  }
  
  template <class T>
  T  window_query(Interval const &win,  
		  T result)
  {
    return (*CRange_tree_3).window_query(win, result);
  }

  ~Range_tree_3()
  {
    if (CRange_tree_3!=0)
      delete CRange_tree_3;
    CRange_tree_3=0;
    if (CRange_tree_2!=0)
      delete CRange_tree_2;
    CRange_tree_2=0;
    if (CRange_tree_1!=0)
      delete CRange_tree_1;
    CRange_tree_1=0;
    if (Tree_anchor!=0)
      delete Tree_anchor;
    Tree_anchor=0;
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

  typedef tree_anchor<Key, Interval> Tree_anchor_type;
  Tree_anchor_type *Tree_anchor;

  typedef Range_tree_d<Key, Interval, I4> Range_tree_1_type;
  Range_tree_1_type * CRange_tree_1;

  typedef Range_tree_d<Key, Interval, I3> Range_tree_2_type;
  Range_tree_2_type * CRange_tree_2;

  typedef Range_tree_d<Key, Interval, I2> Range_tree_3_type;
  Range_tree_3_type *CRange_tree_3;

  typedef Range_tree_d<Key, Interval, I1> Range_tree_4_type;
  Range_tree_4_type *CRange_tree_4;

  Range_tree_4()
  {
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
    CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
    CRange_tree_4 = new Range_tree_4_type(*CRange_tree_3);
  }

  template <class T>
  Range_tree_4(T& first, 
	       T& last)  {
   Tree_anchor = new Tree_anchor_type;
   CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
   CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
   CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
   CRange_tree_4 = new Range_tree_4_type(*CRange_tree_3);
   (*CRange_tree_4).make_tree(first,last);
  }

  template <class T>
  bool make_tree(T& first, 
		 T& last)
  {
    delete CRange_tree_4;
    delete CRange_tree_3;
    delete CRange_tree_2;
    delete CRange_tree_1;
    delete Tree_anchor;
    Tree_anchor = new Tree_anchor_type;
    CRange_tree_1 = new Range_tree_1_type(*Tree_anchor);
    CRange_tree_2 = new Range_tree_2_type(*CRange_tree_1);
    CRange_tree_3 = new Range_tree_3_type(*CRange_tree_2);
    CRange_tree_4 = new Range_tree_4_type(*CRange_tree_3);
    return (*CRange_tree_4).make_tree(first,last);
  }

  template <class T>
  T  window_query(Interval const &win,  
		  T result)
  {
    return (*CRange_tree_4).window_query(win, result);
  }

  ~Range_tree_4()
  {
    if (CRange_tree_4!=0)
      delete CRange_tree_4;
    CRange_tree_4=0;
    if (CRange_tree_3!=0)
      delete CRange_tree_3;
    CRange_tree_3=0;
    if (CRange_tree_2!=0)
      delete CRange_tree_2;
    CRange_tree_2=0;
    if (CRange_tree_1!=0)
      delete CRange_tree_1;
    CRange_tree_1=0;
    if (Tree_anchor!=0)
      delete Tree_anchor;
    Tree_anchor=0;
  }
};

CGAL_END_NAMESPACE
#endif



