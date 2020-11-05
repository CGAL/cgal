// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_GPS_MERGE_H
#define CGAL_GPS_MERGE_H

#include <CGAL/license/Boolean_set_operations_2.h>


#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_xor_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_intersection_visitor.h>
#include <vector>

namespace CGAL {

/*!
  \file   Gps_merge.h
  \brief  This file contains classes that are responsible for merging
          two sets of polygons in the divide-and-conquer algorithm.
          The file contains 3 mergers: Join_merge, Intersection_merge and
          Xor_merge. Join_merge is used when we want to merge the two sets,
          Intersection_merge is used for intersection, and Xor_merge is used
          for symmetric difference.
*/

//! Base_merge
/*! Base_merge is the base class for all merger classes.
    All merges used BFS algorithm with a different visitor when discovering
    a new face.
 */
template <class Arrangement_, class Visitor_>
class Base_merge
{
  typedef Arrangement_                                Arrangement_2;
  typedef Visitor_                                    Visitor;
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

    const typename Arrangement_2::Geometry_traits_2 * tr =
      arr_vec[i].first->geometry_traits();
    Arrangement_2              *res = new Arrangement_2(tr);
    std::vector<Vertex_handle> *verts = new std::vector<Vertex_handle>;

    Gps_agg_op<Arrangement_2, Visitor>
      agg_op(*res, *verts, *(res->traits_adaptor()));
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

//! Join_merge
/*! Join_merge is used to join two sets of polygons together in the D&C
    algorithm. It is a base merge with a visitor that joins faces.
 */
template <class Arrangement_>
class Join_merge : public Base_merge<Arrangement_,
  Gps_bfs_join_visitor<Arrangement_> >
{};


//! Intersection_merge
/*! Intersection_merge is used to merge two sets of polygons creating their
    intersection.
 */
template <class Arrangement_>
class Intersection_merge : public Base_merge<Arrangement_,
  Gps_bfs_intersection_visitor<Arrangement_> >
{};

//! Xor_merge
/*! Xor_merge is used to merge two sets of polygons creating their
    symmetric difference.
 */
template <class Arrangement_>
class Xor_merge : public Base_merge<Arrangement_,
  Gps_bfs_xor_visitor<Arrangement_> >
{
};

} //namespace CGAL

#endif
