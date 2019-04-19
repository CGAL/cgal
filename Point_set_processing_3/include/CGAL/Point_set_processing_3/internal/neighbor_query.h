// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H
#define CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/point_set_processing_assertions.h>

#include <boost/function_output_iterator.hpp>

namespace CGAL {
namespace Point_set_processing_3 {
namespace internal {

struct Maximum_points_reached_exception : public std::exception { };

template <typename Point,
          typename TreeTraits,
          typename TreeSplitter,
          typename TreeUseExtendedNode,
          typename FT,
          typename PointContainer>
void neighbor_query (const Point& query,
                     const CGAL::Kd_tree<TreeTraits, TreeSplitter, TreeUseExtendedNode>& tree,
                     unsigned int k,
                     FT neighbor_radius,
                     PointContainer& points)
{
  typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;
  typedef CGAL::Fuzzy_sphere<TreeTraits> Sphere;
  
  if (neighbor_radius != FT(0))
  {
    Sphere fs (query, neighbor_radius, 0, tree.traits());

    // if k=0, no limit on the number of neighbors returned
    if (k == 0)
      k = std::numeric_limits<unsigned int>::max();
    
    try
    {
      std::function<void(const Point&)> back_insert_with_limit
        = [&](const Point& point) -> void
        {
          points.push_back (point);
          if (points.size() == k)
            throw Maximum_points_reached_exception();
        };
      
      auto function_output_iterator
        = boost::make_function_output_iterator (back_insert_with_limit);

      tree.search (function_output_iterator, fs);
    }
    catch (const Maximum_points_reached_exception&)
    { }

    // Fallback, if less than 3 points are return, search for the 3
    // first points
    if (points.size() < 3)
      k = 3;
    // Else, no need to search for K nearest neighbors
    else
      k = 0;
  }
  
  if (k != 0)
  {
    // Gather set of (k+1) neighboring points.
    // Perform k+1 queries (as in point set, the query point is
    // output first). Search may be aborted if k is greater
    // than number of input points.
    points.reserve(k+1);
    Neighbor_search search(tree,query,k+1);
    Search_iterator search_iterator = search.begin();
    unsigned int i;
    for(i=0;i<(k+1);i++)
    {
      if(search_iterator == search.end())
        break; // premature ending
      points.push_back(search_iterator->first);
      search_iterator++;
    }
    CGAL_point_set_processing_precondition(points.size() >= 1);
  }
}

} } } // namespace CGAL::Point_set_processing_3::internal

#endif // CGAL_PSP_INTERNAL_NEIGHBOR_QUERY_H
