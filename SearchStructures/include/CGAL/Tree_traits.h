// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Gabriele Neyer

#ifndef CGAL_TREE_TRAITS_H
#define CGAL_TREE_TRAITS_H

#include <CGAL/license/SearchStructures.h>


// Implementation of a minimal tree traits for CGAL trees, derived
// from Tree_base.h
// (e.g. Range_tree_d.h and Segment_tree_d.h).
// Any other interface must at least provide the interface of this class.

namespace CGAL {

// Interface as it is expected for one layer of a Range Tree
template<class C_Data, class C_Window, class C_Key,
         class C_Data_func, class C_Window_left_func,
         class C_Window_right_func, class C_Compare>
class tree_point_traits{
 public:
  typedef  C_Key Key;
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef  C_Compare Compare;

   C_Key get_key(const  C_Data&  dt)
    const {return  C_Data_func()(dt);}
   C_Key get_left(const  C_Window& dt)
    const {return  C_Window_left_func()(dt);}
   C_Key get_right(const  C_Window& dt)
    const {return  C_Window_right_func()(dt);}

  static bool comp( C_Key const& key1,  C_Key const& key2)
  {return  C_Compare()(key1, key2);}
  static bool key_comp( C_Data const& data1,  C_Data const& data2)
  {
    const Key d1= C_Data_func()(data1);
    const Key d2= C_Data_func()(data2);
    return  C_Compare()(d1,d2);
  }
};


// Interface as it is expected for one layer of a Segment Tree
template<class C_Data, class C_Window, class C_Key,
         class C_Data_left_func, class C_Data_right_func,
         class C_Window_left_func, class C_Window_right_func, class C_Compare>
class tree_interval_traits{
 public:
  typedef  C_Data Data;
  typedef  C_Window Window;
  typedef  C_Key Key;
  typedef  C_Compare Compare;

   C_Key get_left(const  C_Data&  dt)
    const {return  C_Data_left_func()(dt);}
   C_Key get_right(const  C_Data& dt)
    const {return  C_Data_right_func()(dt);}
   C_Key get_left_win(const  C_Window& dt)
    const {return  C_Window_left_func()(dt);}
   C_Key get_right_win(const  C_Window& dt)
    const {return  C_Window_right_func()(dt);}

  static bool comp( C_Key const& key1,  C_Key const& key2)
    {return  C_Compare()(key1, key2);}
};

} //namespace CGAL

#endif
