// Copyright (c) 2009-2013 GeometryFactory (France).
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
// 
//
// Author(s)     : Laurent Rineau and Ilker O. Yaz


#ifndef CGAL_ORIENT_POLYGON_SOUP
#define CGAL_ORIENT_POLYGON_SOUP

#include <set>
#include <stack>
#include <algorithm>
#include <iostream>

#include <CGAL/array.h>

namespace CGAL {

namespace internal {

template<class Point_3, class Polygon_3>
class Polygon_soup_orienter
{
  typedef typename std::iterator_traits<typename Polygon_3::iterator>::value_type Index;
  typedef std::vector<Point_3> Points;
  typedef std::map<std::pair<Index, Index>, std::set<std::size_t> > Edges_map;
  typedef cpp11::array<Index, 2> Edge;
  typedef std::vector<Polygon_3> Polygons;
  typedef std::set<Edge> Edges;
  typedef typename Polygons::size_type size_type;

  const Points& points;
  Polygons&     polygons;

  Edges_map edges;
  Edges non_manifold_edges;

public:
  Polygon_soup_orienter(const Points& points, Polygons& polygons)
    : points(points), polygons(polygons)
  {
    fill_edges();
  }

private:
  Edge canonical_edge(Index i, Index j)
  {
    return i<j ? CGAL::make_array(i,j):CGAL::make_array(j,i);
  }

  void fill_edges() {
    // Fill edges
    edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const Index& i0 = polygons[i][j];
        const Index& i1 = polygons[i][ j+1 < size ? j+1: 0];
        edges[std::make_pair(i0, i1)].insert(i);
      }
    }

    // Fill non-manifold edges
    non_manifold_edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const Index& i0 = polygons[i][j];
        const Index& i1 = polygons[i][ j+1 < size ? j+1: 0];

        if( edges[std::make_pair(i0, i1)].size() +
            edges[std::make_pair(i1, i0)].size() > 2 )
        {
          non_manifold_edges.insert(canonical_edge(i0,i1));
        }
      }
    }
  }

  void inverse_orientation(const std::size_t index) {
    std::reverse(polygons[index].begin(), polygons[index].end());
  }

public:
  bool orient()
  {
    std::vector<bool> oriented;
    std::stack<std::size_t> stack;
    using std::make_pair;

    // no polygon is oriented
    oriented.resize(polygons.size());

    size_type polygon_index = 0;
    bool success = true;

    while (polygon_index != polygons.size()) 
    {
      while ( polygon_index != polygons.size() && oriented[polygon_index] ) {
        ++polygon_index;
      }
      if(polygon_index == polygons.size()) break;

      oriented[polygon_index] = true;
      stack.push(polygon_index);
      while(! stack.empty() )
      {
        const size_type to_be_oriented_index = stack.top();
        stack.pop();
        const size_type size = polygons[to_be_oriented_index].size();
        for(size_type ih = 0 ; ih < size ; ++ih) {
          size_type ihp1 = ih+1;
          if(ihp1>=size) ihp1 = 0;
          const Index& i1 = polygons[to_be_oriented_index][ih];
          const Index& i2 = polygons[to_be_oriented_index][ihp1];

          if(non_manifold_edges.count(canonical_edge(i1,i2)) > 0) {
            continue;
          }

          // edge (i1,i2)
          typename Edges_map::iterator it_same_orient = edges.find(make_pair(i1, i2));
          // edges (i2,i1)
          typename Edges_map::iterator it_other_orient = edges.find(make_pair(i2, i1));

          CGAL_assertion(it_same_orient != edges.end());
          if(it_same_orient->second.size() > 1) {
            if((it_other_orient != edges.end() && it_other_orient->second.size() > 0) ||
              it_same_orient->second.size() > 2) {
                // three polygons at the edge
                success = false; // non-orientable
            }
            {
              // one neighbor polyhedron, opposite orientation
              size_type index = *(it_same_orient->second.begin());
              if(index == to_be_oriented_index)
                index = *(++it_same_orient->second.begin());
              if(oriented[index]) {
                //  "neighbor polygon #%1 is already oriented, but in opposite orientation").arg(index);
                success = false; // non-orientable
                continue; // next edge
              }

              // reverse the orientation
              const size_type size = polygons[index].size();
              for(size_type j = 0; j < size; ++j) {
                const Index& i0 = polygons[index][j];
                const Index& i1 = polygons[index][ j+1 < size ? j+1: 0];
                CGAL_assertion_code(const bool r = )
                  edges[std::make_pair(i0, i1)].erase(index)
                CGAL_assertion_code(!= 0);
                CGAL_assertion(r);
              }
              inverse_orientation(index);
              for(size_type j = 0; j < size; ++j) {
                const Index& i0 = polygons[index][j];
                const Index& i1 = polygons[index][ j+1 < size ? j+1: 0];
                edges[std::make_pair(i0, i1)].insert(index);
              }
              // "inverse the orientation of polygon #%1\n").arg(index);
              oriented[index] = true;
              stack.push(index);
            }
          }
          else if(it_other_orient != edges.end() && it_other_orient->second.size() == 1) {
            // one polygon, same orientation
            const size_type index = *(it_other_orient->second.begin());
            if(oriented[index])
              continue;
            oriented[index] = true;
            // "keep the orientation of polygon #%1\n").arg(index);
            stack.push(index);
          }
          else {
            if(it_same_orient->second.size() != 1 || 
              (it_other_orient != edges.end() && it_other_orient->second.size() > 0)) 
            {
              success = false; // non-orientable
            }
          }
        } // end for on all edges of one 
      } // end while loop on the polygons of the connected component
    } // end while loop on all non-oriented polygons remaining 

    return success;
  }
};
} // namespace internal

/** 
 * Tries to consistently orient a soup of polygons in 3D space.  
 * If a consistent orientation has been found, `true` is returned.
 * In any case `polygons` is updated.
 * @tparam Point_3 the point type
 *
 * @param points points of the soup of polygons.
 * @param[in, out] polygons each element in the vector describes a polygon using the index of the points in the vector.
 *
 * @return true if a consistent orientation has been found
 *
 * \TODO code: there is no check for duplicate points, yet it can be implemented as separate filter function
 */ 
template <class Point_3, class Polygon_3>
bool orient_polygon_soup(const std::vector<Point_3>& points,
                         std::vector< Polygon_3 >& polygons)
{
  internal::Polygon_soup_orienter<Point_3, Polygon_3> orienter(points, polygons);
  return orienter.orient();
}

}// namespace CGAL

#endif // CGAL_ORIENT_POLYGON_SOUP
