// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef LCC_SEW3_SIMILAR_FACETS_H
#define LCC_SEW3_SIMILAR_FACETS_H

#include <vector>
#include <map>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL::internal
{
///////////////////////////////////////////////////////////////////////////////
template<class FT>
bool are_length_similar(FT d1, FT d2, double epsilon=1e-9)
{ return CGAL::abs(d1-d2)<=epsilon*epsilon; }
///////////////////////////////////////////////////////////////////////////////
template<class Point>
bool are_point_similar(const Point& p1, const Point& p2, double epsilon=1e-9)
{ return are_length_similar(CGAL::squared_distance(p1, p2), 0., epsilon); }
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
void print_face(LCC& lcc,
                typename LCC::Dart_descriptor d)
{
  std::cout<<"[  ";
  typename LCC::Dart_descriptor cur=d;
  do
  {
    std::cout<<"("<<lcc.point(d)<<")  ";
    d=lcc.next(d);
  }
  while(cur!=d);
  std::cout<<"]"<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
bool are_facets_opposite_and_similar_geometry(LCC& lcc,
                                              typename LCC::Dart_descriptor d1,
                                              typename LCC::Dart_descriptor d2,
                                              double epsilon=1e-9)
{
  if(lcc.other_extremity(d2)!=lcc.null_descriptor &&
     !are_point_similar(lcc.point(d1), lcc.point(lcc.other_extremity(d2)), epsilon))
  { return false; }

  using Dart_descriptor=typename LCC::Dart_descriptor;
  Dart_descriptor s1=d1;
  Dart_descriptor s2=d2;
  while (lcc.is_previous_exist(d1) && lcc.previous(s1)!=d1)
  {
    s1=lcc.previous(s1);
    if (!lcc.is_next_exist(d2)) { return false; }
    s2=lcc.next(s2);
  }

  d1=s1;
  d2=s2;
  do
  {
    if (lcc.is_next_exist(d1)!=lcc.is_previous_exist(d2))
    { return false; }

    if (lcc.other_extremity(d2)!=lcc.null_descriptor &&
        !are_point_similar(lcc.point(d1), lcc.point(lcc.other_extremity(d2)), epsilon))
    { return false; }

    // The only case where d2 could have no other_extremity
    // is the end of an open path. In this case d1 is the
    // beginning of an open path and we do not compare points
    // but this is the correct thing to do.

    d1=lcc.next(d1);
    d2=lcc.previous(d2);
  }
  while(lcc.is_next_exist(d1) && d1!=s1);

  if (lcc.is_next_exist(d1)!=lcc.is_previous_exist(d2))
  { return false; }

  if (d1==s1 && d2!=s2) { return false; }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
/// @return sew3 all pairs of facets that are "similar", given epsilon
/// (a facet is considered marked if ALL its darts are marked).
/// Only marked faces are proceed, but they can be 3-sewn with non
/// marked faces.
template<class LCC>
std::size_t sew3_similar_facets(LCC& lcc,
                                typename LCC::size_type AMark,
                                double epsilon=1e-9)
{
  std::size_t res=0;
  using Dart_descriptor=typename LCC::Dart_descriptor;
  using FT=typename LCC::FT;
  std::map<std::size_t, std::vector<std::tuple<Dart_descriptor, FT, FT>>>
      one_dart_per_facet;
  auto mymark=lcc.get_new_mark();
  typename LCC::Point p0=CGAL::ORIGIN;

  // First we fill the std::map by one dart per facet, and by using
  // the number of edges as index.
  for (auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
  {
    if (!lcc.is_marked(it, mymark))
    {
      if(!lcc.template is_opposite_exist<3>(it))
      {
        FT min_lg=CGAL::squared_distance(lcc.point(it), p0);
        FT max_lg=min_lg;
        std::size_t nb=1;
        auto it2=lcc.template darts_of_cell_basic<2>(it, mymark).begin();
        lcc.mark(it2, mymark);
        for(++it2; it2.cont(); ++it2, ++nb )
        {
          lcc.mark(it2, mymark);
          FT cur_lg=CGAL::squared_distance(lcc.point(it2), p0);
          if (cur_lg<min_lg) { min_lg=cur_lg; }
          if (cur_lg>max_lg) { max_lg=cur_lg; }
        }
        one_dart_per_facet[nb].push_back(std::make_tuple(it, min_lg, max_lg));
      }
      else
      { lcc.mark(it, mymark); }
    }
  }

  // Second we run through the map: candidates for sew3 have necessary the
  // same minimal and maximal points.
  for (auto itmap=one_dart_per_facet.begin(),
       itmapend=one_dart_per_facet.end(); itmap!=itmapend; ++itmap)
  {
    for (auto it1=(itmap->second).begin(),
         it1end=(itmap->second).end(); it1!=it1end; ++it1)
    {
      // We only proceed 3-free marked faces for it1
      if (!lcc.template is_opposite_exist<3>(std::get<0>(*it1)) &&
          lcc.is_marked(std::get<0>(*it1), AMark))
      {
        //std::cout<<"**********************************"<<std::endl;
        //std::cout<<"F1->  "; print_face(lcc, std::get<0>(*it1));
        auto it2=(itmap->second).begin();// DEBUG it1; // it2 iterate through the std::vector, to test all
                      // possible candidates to sew with it1
        for (/* ++it2 */; it2!=it1end; )
        {
          //std::cout<<"F2->  "; print_face(lcc, std::get<0>(*it2));
          if (it1!=it2 &&
              !lcc.template is_opposite_exist<3>(std::get<0>(*it2)) &&
              are_length_similar(std::get<1>(*it1), std::get<1>(*it2), epsilon) &&
              are_length_similar(std::get<2>(*it1), std::get<2>(*it2), epsilon))
          {
            Dart_descriptor curdh=std::get<0>(*it2);
            do
            {
              CGAL_assertion(lcc.is_previous_exist(curdh));
              if(are_facets_opposite_and_similar_geometry
                  (lcc, std::get<0>(*it1), lcc.previous(curdh), epsilon))
              {
                ++res;
                lcc.template sew<3>(std::get<0>(*it1),
                                    lcc.other_orientation(lcc.previous(curdh)));
              }
              else
              { curdh=lcc.previous(curdh); }
            }
            while(!lcc.template is_opposite_exist<3>(std::get<0>(*it1)) &&
                  curdh!=std::get<0>(*it2));
          }
          if(lcc.template is_opposite_exist<3>(std::get<0>(*it1)))
          { it2=it1end; } // to leave the "for loop" since it1 is no more 3-free
          else { ++it2; }
        }
      }
    }
  }
  CGAL_assertion(lcc.is_whole_map_marked(mymark));
  lcc.free_mark(mymark);
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Sew3 the facets having similar geometry
/// (all the facets of the map are considered)
template<typename LCC>
std::size_t sew3_similar_facets(LCC& lcc,
                                double epsilon=1e-9)
{
  auto mark=lcc.get_new_mark();
  lcc.negate_mark(mark);
  std::size_t res=sew3_similar_facets(lcc, mark, epsilon);
  lcc.free_mark(mark);
  return res;
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_SEW3_SIMILAR_FACETS_H
