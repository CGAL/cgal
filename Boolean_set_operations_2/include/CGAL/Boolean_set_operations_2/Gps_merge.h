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

#ifndef CGAL_GPS_MERGE_H
#define CGAL_GPS_MERGE_H

#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_xor_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_intersection_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Join_merge
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef std::pair<Arrangement_2 *,
                    std::vector<Vertex_handle> *>     Arr_entry;

public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arr_entry>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement_2::Traits_2*  tr = arr_vec[i].first->traits();
    Arrangement_2              *res = new Arrangement_2(tr);
    std::vector<Vertex_handle> *verts = new std::vector<Vertex_handle>;

    Gps_agg_op<Arrangement_2, Gps_bfs_join_visitor<Arrangement_2> >
      agg_op(*res, *verts, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete (arr_vec[count].first);
      delete (arr_vec[count].second);
    }
    
    arr_vec[i].first = res;
    arr_vec[i].second = verts;
  }

};


template <class Arrangement_>
class Intersection_merge
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef std::pair<Arrangement_2 *,
                    std::vector<Vertex_handle> *>     Arr_entry;

public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arr_entry>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement_2::Traits_2*  tr = arr_vec[i].first->traits();
    Arrangement_2              *res = new Arrangement_2 (tr);
    std::vector<Vertex_handle> *verts = new std::vector<Vertex_handle>;

    Gps_agg_op<Arrangement_2, Gps_bfs_intersection_visitor<Arrangement_2> >
      agg_op(*res, *verts, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete (arr_vec[count].first);
      delete (arr_vec[count].second);
    }
    
    arr_vec[i].first = res;
    arr_vec[i].second = verts;
  }
};

template <class Arrangement_>
class Xor_merge
{
  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef std::pair<Arrangement_2 *,
                    std::vector<Vertex_handle> *>     Arr_entry;

public:

  // Temporarily defined to see if this avoids a warning on SunPro CC
  Xor_merge()
  {}

   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arr_entry>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement_2::Traits_2*  tr = arr_vec[i].first->traits();
    Arrangement_2              *res = new Arrangement_2(tr);
    std::vector<Vertex_handle> *verts = new std::vector<Vertex_handle>;

    Gps_agg_op<Arrangement_2, Gps_bfs_xor_visitor<Arrangement_2> >
      agg_op(*res, *verts, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete (arr_vec[count].first);
      delete (arr_vec[count].second);
    }
    
    arr_vec[i].first = res;
    arr_vec[i].second = verts;
  }
};

CGAL_END_NAMESPACE

#endif
