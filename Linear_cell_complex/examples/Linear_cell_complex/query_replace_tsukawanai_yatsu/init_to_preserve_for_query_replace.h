// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef INIT_TO_PRESERVE_FOR_QUERY_REPLACE_H
#define INIT_TO_PRESERVE_FOR_QUERY_REPLACE_H

#include"lcc_geometrical_tests.h"
#include <CGAL/Bbox_3.h>

///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool is_square_extremal_point(const Point& p, const CGAL::Bbox_3& b)
{
  return (b.xmin()==b.xmax() || (p.x()==b.xmin() || p.x()==b.xmax())) &&
      (b.ymin()==b.ymax() || (p.y()==b.ymin() || p.y()==b.ymax())) &&
      (b.zmin()==b.zmax() || (p.z()==b.zmin() || p.z()==b.zmax()));
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool is_edge_on_square_border(const Point& p1, const Point& p2,
                              const CGAL::Bbox_3& b)
{
  return (b.xmin()!=b.xmax() && p1.x()==p2.x()) ||
      (b.ymin()!=b.ymax() && p1.y()==p2.y()) ||
      (b.zmin()!=b.zmax() && p1.z()==p2.z());
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool is_hexa_extremal_point(const Point& p, const CGAL::Bbox_3& b)
{
  return (p.x()==b.xmin() || p.x()==b.xmax()) &&
      (p.y()==b.ymin() || p.y()==b.ymax()) &&
      (p.z()==b.zmin() || p.z()==b.zmax());
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool is_edge_on_hexa_border(const Point& p1, const Point& p2)
{
  return ((p1.x()==p2.x() && p1.y()==p2.y()) ||
          (p1.x()==p2.x() && p1.z()==p2.z()) ||
          (p1.y()==p2.y() && p1.z()==p2.z()));
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_fpattern_corners(LCC& pattern, typename LCC::size_type mark_to_preserve)
{
  double epsilon=.1;
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend; ++it)
  {
    if(pattern.template is_free<2>(it))
    {
      typename LCC::Dart_handle prev=pattern.template beta<0>(it);
      while(!pattern.template is_free<2>(prev))
      { prev=pattern.template beta<2,0>(prev); }
      if(!lcc_tests::aligned(pattern.point(prev),
                             pattern.point(it),
                             pattern.point(pattern.other_extremity(it)), epsilon))
      { pattern.mark(it, mark_to_preserve); }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_face_corners(LCC& lcc, typename LCC::Dart_handle dh,
                       typename LCC::size_type mark_to_preserve)
{
  double epsilon=.1;
  typename LCC::Dart_handle cur=dh;
  do
  {
    if(!lcc_tests::aligned(lcc.point(lcc.beta(cur, 0)),
                           lcc.point(cur),
                           lcc.point(lcc.beta(cur, 1)), epsilon))
    { lcc.mark(cur, mark_to_preserve); }
    cur=lcc.template beta<1>(cur);
  }
  while(cur!=dh);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_spattern_edges(LCC& pattern, typename LCC::size_type mark_border)
{
  // Mark edges between non coplanar faces
  double epsilon=.1;
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend; ++it)
  {
    assert(!pattern.template is_free<2>(it) && pattern.template is_free<3>(it));
    if(!pattern.is_marked(it, mark_border) && it<pattern.template beta<2>(it))
    {
      typename LCC::Vector n1=CGAL::compute_normal_of_cell_2(pattern, it);
      typename LCC::Vector n2=CGAL::compute_normal_of_cell_2
                              (pattern, pattern.template beta<2>(it));
      if(!lcc_tests::coplanar(n1, n2, epsilon))
      { pattern.template mark_cell<1>(it, mark_border); }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_vpattern_corners(LCC& pattern, typename LCC::size_type mark_to_preserve)
{
  // 1) mark edges between non coplanar faces
  // 2) mark origin of these edges having its previous marked edge non align
  double epsilon=.1;
  typename LCC::Dart_handle dh2;
  typename LCC::size_type cube_border=pattern.get_new_mark();
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend; ++it)
  {
    if(pattern.template is_free<3>(it) && !pattern.is_marked(it, cube_border))
    {
      dh2=pattern.template beta<2>(it);
      while(!pattern.template is_free<3>(dh2))
      { dh2=pattern.template beta<3,2>(dh2); }
      typename LCC::Vector n1=CGAL::compute_normal_of_cell_2(pattern, it);
      typename LCC::Vector n2=CGAL::compute_normal_of_cell_2(pattern, dh2);
      if(!lcc_tests::coplanar(n1, n2, epsilon))
      { pattern.template mark_cell<1>(it, cube_border); }
    }
  }
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend; ++it)
  {
    if(pattern.is_marked(it, cube_border) && pattern.template is_free<3>(it))
    {
      dh2=pattern.template beta<0>(it);
      while(!pattern.is_marked(dh2, cube_border))
      {
        dh2=pattern.template beta<2>(dh2);
        while(!pattern.template is_free<3>(dh2))
        { dh2=pattern.template beta<3,2>(dh2); }
        dh2=pattern.template beta<0>(dh2);
      }
      if(!lcc_tests::aligned(pattern.point(dh2), pattern.point(it),
                             pattern.point(pattern.other_extremity(it)), epsilon))
      { pattern.mark(it, mark_to_preserve); }
    }
  }
  for(auto it=pattern.darts().begin(), itend=pattern.darts().end(); it!=itend;
      ++it)
  { pattern.unmark(it, cube_border); }
  pattern.free_mark(cube_border);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_volume_corners(LCC& lcc, typename LCC::Dart_handle dh,
                         typename LCC::size_type mark_to_preserve)
{
  // 1) mark edges between non coplanar faces
  // 2) mark origin of these edges having its previous marked edge non align
  double epsilon=.1;
  typename LCC::Dart_handle dh2;
  typename LCC::size_type amark=lcc.get_new_mark();
  typename LCC::size_type cube_border=lcc.get_new_mark();
  for(auto it=lcc.template darts_of_cell_basic<3>(dh, amark).begin(),
        itend=lcc.template darts_of_cell_basic<3>(dh, amark).end(); it!=itend; ++it)
  {
    lcc.mark(it, amark);
    if(!lcc.is_marked(it, cube_border))
    {
      dh2=lcc.template beta<2>(it);
      typename LCC::Vector n1=CGAL::compute_normal_of_cell_2(lcc, it);
      typename LCC::Vector n2=CGAL::compute_normal_of_cell_2(lcc, dh2);
      if(!lcc_tests::coplanar(n1, n2, epsilon))
      { lcc.template mark_cell<1>(it, cube_border); }
    }
  }
  lcc.negate_mark(amark);
  for(auto it=lcc.template darts_of_cell_basic<3>(dh, amark).begin(),
        itend=lcc.template darts_of_cell_basic<3>(dh, amark).end(); it!=itend; ++it)
  {
    lcc.mark(it, amark);
    if(lcc.is_marked(it, cube_border))
    {
      dh2=lcc.template beta<0>(it);
      while(!lcc.is_marked(dh2, cube_border))
      { dh2=lcc.template beta<2,0>(dh2); }
      if(!lcc_tests::aligned(lcc.point(dh2), lcc.point(it),
                             lcc.point(lcc.other_extremity(it)), epsilon))
      { lcc.mark(it, mark_to_preserve); }
    }
  }
  lcc.negate_mark(amark);

  for(auto it=lcc.template darts_of_cell<3>(dh).begin(),
        itend=lcc.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(lcc.is_marked(it, cube_border))
    { lcc.template unmark_cell<1>(it, cube_border); }
  }

  assert(lcc.is_whole_map_unmarked(amark));
  assert(lcc.is_whole_map_unmarked(cube_border));
  lcc.free_mark(amark);
  lcc.free_mark(cube_border);
}
////////////////////////////////////////////////////////////////////////////////
#endif // INIT_TO_PRESERVE_FOR_QUERY_REPLACE_H
