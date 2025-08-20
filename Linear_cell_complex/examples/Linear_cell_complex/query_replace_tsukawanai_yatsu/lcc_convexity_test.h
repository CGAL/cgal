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
#ifndef LCC_CONVEXITY_TEST_H
#define LCC_CONVEXITY_TEST_H
///////////////////////////////////////////////////////////////////////////////
#include <CGAL/squared_distance_3.h>
#include <CGAL/Linear_cell_complex_operations.h>
///////////////////////////////////////////////////////////////////////////////
/// @return true iff face(df) is convex
template<class LCC>
bool is_face_convex(LCC& lcc,
                    typename LCC::Dart_handle dh)
{
  return true;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void mark_vertices_of_face(LCC& lcc, typename LCC::Dart_handle dh,
                           typename LCC::size_type mark)
{
  typename LCC::Dart_handle dh2=dh;
  do
  {
    if(!lcc.is_marked(dh2, mark)) // Already marked => dangling or inner edge
    { lcc.template mark_cell<0>(dh2, mark); }
    dh2=lcc.next(dh2);
  }
  while(dh2!=dh);
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
void unmark_vertices_of_face(LCC& lcc, typename LCC::Dart_handle dh,
                             typename LCC::size_type mark)
{
  lcc.negate_mark(mark);
  mark_vertices_of_face(lcc, dh, mark);
  lcc.negate_mark(mark);
}
///////////////////////////////////////////////////////////////////////////////
/// @return true iff volume(df) is convex
/// Use EPSILON to test if a point belongs to a plane.
template<class LCC>
bool is_volume_convex(LCC& lcc,
                      typename LCC::Dart_handle dh,
                      const typename LCC::FT EPSILON=0.0001)
{ // TODO we can improve the test: in the current code, we test face f1
  // and its adjacent faces f2; and reciprocally for f2 we test again f1.
  bool first=true;
  bool positive=true;
  CGAL::Oriented_side side;
  auto mark=lcc.get_new_mark();
  bool on_plane=false;

  for(auto it=lcc.template one_dart_per_incident_cell<2,3>(dh).begin(),
      itend=lcc.template one_dart_per_incident_cell<2,3>(dh).end(); it!=itend; ++it)
  {
    typename LCC::Traits::Plane_3 plane(lcc.point(it),
                                        CGAL::compute_normal_of_cell_2(lcc, it));
    mark_vertices_of_face(lcc, it, mark);
    typename LCC::Dart_handle dh2=it;
    do
    {
      auto dh3=lcc.opposite2(dh2);
      do
      {
        if (!lcc.is_marked(dh3, mark))
        {
          on_plane=(CGAL::squared_distance(lcc.point(dh3), plane)<EPSILON);
          if (!on_plane)
          {
            side=plane.oriented_side(lcc.point(dh3));
            assert(side!=CGAL::ON_ORIENTED_BOUNDARY);
            if (first)
            {
              if (side==CGAL::ON_POSITIVE_SIDE) { first=false; positive=true; }
              else if (side==CGAL::ON_NEGATIVE_SIDE) { first=false; positive=false; }
            }
            else
            {
              if ((side==CGAL::ON_POSITIVE_SIDE && !positive) ||
                  (side==CGAL::ON_NEGATIVE_SIDE && positive))
              {
                unmark_vertices_of_face(lcc, it, mark);
                lcc.free_mark(mark);
                return false;
              }
            }
          }
        }
        dh3=lcc.next(dh3);
      }
      while(dh3!=lcc.opposite2(dh2));

      dh2=lcc.next(dh2);
    }
    while(dh2!=it);
    unmark_vertices_of_face(lcc, it, mark);
  }
  lcc.free_mark(mark);
  return true;
}
///////////////////////////////////////////////////////////////////////////////
/// @return true iff face(df) union face(beta2(dh)) is convex
template<class LCC>
bool adjacent_faces_are_convex(LCC& lcc,
                               typename LCC::Dart_handle dh)
{
  if(lcc.template is_free<2>(dh)) { return false; }
  return true;
}
///////////////////////////////////////////////////////////////////////////////
/// @return true iff volume(df) union volume(beta3(dh)) is convex
/// Use EPSILON to test if a point belongs to a plane.
template<class LCC>
bool adjacent_volume_are_convex(LCC& lcc,
                                typename LCC::Dart_handle dh,
                                const typename LCC::FT EPSILON=0.0001)
{ // TODO we can improve the test: in the current code, we test face f1
  // and its adjacent faces f2; and reciprocally for f2 we test again f1.
  if(lcc.template is_free<3>(dh)) { return false; }
  bool first=true;
  bool positive=true;
  CGAL::Oriented_side side;
  auto mark=lcc.get_new_mark();
  auto markface=lcc.get_new_mark();
  typename LCC::Dart_handle sd=dh;
  bool on_plane=false;

  lcc.template mark_cell<2>(dh, markface);
  for(int i=0; i<2; ++i)
  {
    for(auto it=lcc.template one_dart_per_incident_cell<2,3>(sd).begin(),
        itend=lcc.template one_dart_per_incident_cell<2,3>(sd).end(); it!=itend; ++it)
    {
      if (!lcc.is_marked(it, markface))
      {
        typename LCC::Traits::Plane_3 plane(lcc.point(it),
                                            CGAL::compute_normal_of_cell_2(lcc, it));
        mark_vertices_of_face(lcc, it, mark);
        typename LCC::Dart_handle dh2=it;
        do
        {
          typename LCC::Dart_handle sdh3=lcc.opposite2(dh2);
          while(lcc.is_marked(sdh3, markface))
          { sdh3=lcc.opposite2(lcc.template opposite<3>(sdh3)); }
          if(sdh3!=lcc.null_handle)
          {
            if(lcc.template beta<3>(sdh3)==dh2)
            { // Case of future dangling face => return false
              unmark_vertices_of_face(lcc, it, mark);
              lcc.template unmark_cell<2>(dh, markface);
              lcc.free_mark(mark);
              lcc.free_mark(markface);
              return false;
            }
            typename LCC::Dart_handle dh3=sdh3;
            do
            {
              if (!lcc.is_marked(dh3, mark))
              {
                on_plane=(CGAL::squared_distance(lcc.point(dh3), plane)<EPSILON);
                if (!on_plane)
                {
                  side=plane.oriented_side(lcc.point(dh3));
                  assert(side!=CGAL::ON_ORIENTED_BOUNDARY);
                  if (first)
                  {
                    if (side==CGAL::ON_POSITIVE_SIDE) { first=false; positive=true; }
                    else if (side==CGAL::ON_NEGATIVE_SIDE) { first=false; positive=false; }
                  }
                  else
                  {
                    if ((side==CGAL::ON_POSITIVE_SIDE && !positive) ||
                        (side==CGAL::ON_NEGATIVE_SIDE && positive))
                    {
                      unmark_vertices_of_face(lcc, it, mark);
                      lcc.template unmark_cell<2>(dh, markface);
                      lcc.free_mark(mark);
                      lcc.free_mark(markface);
                      return false;
                    }
                  }
                }
              }
              dh3=lcc.next(dh3);
            }
            while(dh3!=sdh3);
          }
          dh2=lcc.next(dh2);
        }
        while(dh2!=it);
        unmark_vertices_of_face(lcc, it, mark);
      }
    }
    sd=lcc.template opposite<3>(dh);
  }
  lcc.template unmark_cell<2>(dh, markface);
  lcc.free_mark(mark);
  lcc.free_mark(markface);

  return true;
}
///////////////////////////////////////////////////////////////////////////////
#endif // LCC_CONVEXITY_TEST_H
