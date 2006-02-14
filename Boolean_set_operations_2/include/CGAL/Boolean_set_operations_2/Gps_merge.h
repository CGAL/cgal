// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_MERGE_H
#define GPS_MERGE_h

#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_xor_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_intersection_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Arrangement>
class Join_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_join_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }

};


template <class Arrangement>
class Intersection_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_intersection_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }
};

template <class Arrangement>
class Xor_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_xor_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }
};

CGAL_END_NAMESPACE

#endif
