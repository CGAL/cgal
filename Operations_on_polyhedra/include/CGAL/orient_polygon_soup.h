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


#ifndef CGAL_ORIENT_POLYHEDRON_3
#define CGAL_ORIENT_POLYHEDRON_3
#include <CGAL/IO/generic_print_polyhedron.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

#include <set>
#include <stack>
#include <algorithm>
#include <iostream>

#include <boost/array.hpp>

namespace CGAL {

namespace internal {

template<class Point_3>
class Polygon_soup_orienter
{
  typedef std::vector<std::size_t> Polygon_3;
  typedef std::vector<Point_3> Points;
  typedef std::map<std::pair<std::size_t, std::size_t>, std::set<std::size_t> > Edges_map;
  typedef boost::array<std::size_t, 2> Edge;
  typedef std::vector<Polygon_3> Polygons;
  typedef std::set<Edge> Edges;
  typedef Polygons::size_type size_type;

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
  void fill_edges() {
    // Fill edges
    edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        edges[std::make_pair(i0, i1)].insert(i);
      }
    }

    // Fill non-manifold edges
    non_manifold_edges.clear();
    for(size_type i = 0; i < polygons.size(); ++i)
    {
      const size_type size = polygons[i].size();
      for(size_type j = 0; j < size; ++j) {
        const std::size_t& i0 = polygons[i][j];
        const std::size_t& i1 = polygons[i][ j+1 < size ? j+1: 0];
        if( (i0 < i1) && 
            (edges[std::make_pair(i0, i1)].size() +
             edges[std::make_pair(i1, i0)].size() > 2) )
        {
          Edge edge;
          edge[0] = i0;
          edge[1] = i1;
          if(i0 > i1) std::swap(edge[0], edge[1]);
          non_manifold_edges.insert(edge);
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
          const std::size_t& i1 = polygons[to_be_oriented_index][ih];
          const std::size_t& i2 = polygons[to_be_oriented_index][ihp1];

          Edge edge;
          edge[0] = i1;
          edge[1] = i2;
          if(i1 > i2) std::swap(edge[0], edge[1]);

          if(non_manifold_edges.count(edge) > 0) {
            continue;
          }

          // edge (i1,i2)
          Edges_map::iterator it_same_orient = edges.find(make_pair(i1, i2));
          // edges (i2,i1)
          Edges_map::iterator it_other_orient = edges.find(make_pair(i2, i1));

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
                const std::size_t& i0 = polygons[index][j];
                const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
                CGAL_assertion_code(const bool r = )
                  edges[std::make_pair(i0, i1)].erase(index);
                CGAL_assertion(r);
              }
              inverse_orientation(index);
              for(size_type j = 0; j < size; ++j) {
                const std::size_t& i0 = polygons[index][j];
                const std::size_t& i1 = polygons[index][ j+1 < size ? j+1: 0];
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
 * \TODO code: support fixed size arrays for polygons, or creating a concept which provides .size and .operator[]
 */ 
template <class Point_3>
bool orient_polygon_soup(const std::vector<Point_3>& points,
                         std::vector< std::vector<std::size_t> >& polygons)
{
  internal::Polygon_soup_orienter<Point_3> orienter(points, polygons);
  return orienter.orient();
}

/**
  * Modifier to build a polyhedron from a soup of polygons.
  */
template <class HDS, class Point_3>
class Polygon_soup_to_polyhedron_3: public CGAL::Modifier_base<HDS>
{
  typedef std::vector<std::size_t> Polygon_3;

  const std::vector<Point_3>& points;
  const std::vector<std::vector<std::size_t> >& polygons;
public:
  /** 
   * The constructor for modifier object.
   * @param points points of the soup of polygons.
   * @param polygons each element in the vector describes a polygon using the index of the points in the vector.
   */
  Polygon_soup_to_polyhedron_3(const std::vector<Point_3>& points, 
    const std::vector<std::vector<std::size_t> >& polygons) 
    : points(points), polygons(polygons)
  { }

  void operator()(HDS& out_hds)
  {
    Polyhedron_incremental_builder_3<HDS> builder(out_hds);

    builder.begin_surface(points.size(), polygons.size());

    for(std::size_t i = 0, end = points.size(); i < end; ++i)
    { builder.add_vertex(points[i]); }

    for(std::size_t i = 0, end = polygons.size(); i < end; ++i)
    {
      const Polygon_3& polygon = polygons[i]; 
      const std::size_t size = polygon.size();

      builder.begin_facet();
      for(std::size_t j = 0; j < size; ++j) {
        builder.add_vertex_to_facet(polygon[j]);
      }
      builder.end_facet();
    }
    builder.end_surface();
  }
};

}// namespace CGAL

#endif // CGAL_ORIENT_POLYHEDRON_3
