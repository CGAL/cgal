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
#ifndef CMAP_ISOMORPHISMS_H
#define CMAP_ISOMORPHISMS_H

#include <CGAL/Handle_hash_function.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Combinatorial_map/internal/Combinatorial_map_internal_functors.h>

#include <queue>
#include <functional>

/// @pre both faces are a cycle of darts (no 1-boundary).
template<typename LCC1, typename LCC2>
bool are_faces_isomorphic_from_darts(LCC1& lcc1, typename LCC1::Dart_handle sd1,
                                     LCC2& lcc2, typename LCC2::Dart_handle sd2,
                                     typename LCC1::size_type marktopreserve1,
                                     typename LCC2::size_type marktopreserve2,
                                     std::function<typename LCC1::Dart_handle
                                     (typename LCC1::Dart_handle)> next1,
                                     std::function<typename LCC2::Dart_handle
                                     (typename LCC2::Dart_handle)> next2,
                                     bool testDartInfo=true,
                                     bool testAttributes=true,
                                     bool testPoint=true)
{
  if(marktopreserve1!=LCC1::INVALID_MARK &&
     marktopreserve2!=LCC2::INVALID_MARK &&
     lcc1.is_marked(sd1, marktopreserve1)!=lcc2.is_marked(sd2, marktopreserve2))
  { return false; }

  bool match=true;
  typename LCC1::Dart_handle dh1=sd1;
  typename LCC2::Dart_handle dh2=sd2;
  do
  {
    if(marktopreserve1!=LCC1::INVALID_MARK &&
       marktopreserve2!=LCC2::INVALID_MARK &&
       lcc1.is_marked(dh1, marktopreserve1)!=lcc2.is_marked(dh2, marktopreserve2))
    { match=false; }

    // We first test info of darts
    if(match && testDartInfo)
    { match=CGAL::internal::template Test_is_same_dart_info_functor<LCC1, LCC2>::
          run(lcc1, lcc2, dh1, dh2); }

    // We need to test in both direction because
    // Foreach_enabled_attributes only test non void attributes
    // of Self. Functor Test_is_same_attribute_functor will modify
    // the value of match to false if attributes do not match
    if(testAttributes)
    {
      if (match)
      { LCC1::Helper::template Foreach_enabled_attributes
            <CGAL::internal::template Test_is_same_attribute_functor
            <LCC1, LCC2>>::run(lcc1, lcc2, dh1, dh2, match); }
      if (match)
      { LCC2::Helper::template Foreach_enabled_attributes
            <CGAL::internal::template Test_is_same_attribute_functor
            <LCC2, LCC1>>::run(lcc2, lcc1, dh2, dh1, match); }
    }

    if(match && testPoint)
    {
      // Only point of 0-attributes are tested. TODO test point of all
      // attributes ?
      match=CGAL::internal::template Test_is_same_attribute_point_functor
            <LCC1, LCC2, 0>::run(lcc1, lcc2, dh1, dh2);
    }

    if(match)
    {
      dh1=next1(dh1);
      dh2=next2(dh2);
      assert(dh1!=lcc1.null_dart_handle);
      assert(dh2!=lcc2.null_dart_handle);
    }
  }
  while(match && dh1!=sd1 && dh2!=sd2);
  if(dh1!=sd1 || dh2!=sd2) { match=false; }

  return match;
}

/// Test if face(sd) is isomorphic with the border of the fpattern
/// @pre both faces are a cycle of darts (no 1-boundary).
/// @return the dart of vpattern satisfying the bijection, nullprt if there is
/// no isomorphism.
template<typename LCC1, typename LCC2>
typename LCC2::Dart_handle is_face_isomorphic_to_fpattern
(LCC1& lcc, typename LCC1::Dart_handle sd,
 LCC2& fpattern,
 typename LCC1::size_type marktopreserve1=LCC1::INVALID_MARK,
 typename LCC2::size_type marktopreserve2=LCC2::INVALID_MARK,
 bool testDartInfo=true,
 bool testAttributes=true,
 bool testPoint=true)
{
  typename LCC2::Dart_handle res=nullptr;
  for(auto it=fpattern.darts().begin(), itend=fpattern.darts().end();
      res==nullptr && it!=itend; ++it)
  {
    if(fpattern.template is_free<2>(it))
    {
      if(are_faces_isomorphic_from_darts(lcc, sd, fpattern, it,
                                         marktopreserve1, marktopreserve2,

                                         [&lcc](typename LCC1::Dart_handle dh)
                                         -> typename LCC1::Dart_handle
                                         { return lcc.template beta<1>(dh); },

                                         [&fpattern](typename LCC2::Dart_handle dh)
                                         -> typename LCC2::Dart_handle
                                         { typename LCC2::Dart_handle other=
                                         fpattern.template beta<1>(dh);
                                         while(!fpattern.template is_free<2>(other))
                                         { other=fpattern.template beta<2,1>(other); }
                                         return other;
                                         },
                                         testDartInfo,
                                         testAttributes,
                                         testPoint))
      { res=it; }
    }
  }
  return res;
}

template<typename LCC1, typename LCC2>
bool are_volumes_isomorphic_from_darts(LCC1& lcc1, typename LCC1::Dart_handle sd1,
                                       LCC2& lcc2, typename LCC2::Dart_handle sd2,
                                       typename LCC1::size_type marktopreserve1,
                                       typename LCC2::size_type marktopreserve2,
                                       std::function<typename LCC1::Dart_handle
                                       (typename LCC1::Dart_handle)> next1,
                                       std::function<typename LCC1::Dart_handle
                                       (typename LCC1::Dart_handle)> opposite1,
                                       std::function<typename LCC2::Dart_handle
                                       (typename LCC2::Dart_handle)> next2,
                                       std::function<typename LCC2::Dart_handle
                                       (typename LCC2::Dart_handle)> opposite2,
                                       bool testDartInfo=true,
                                       bool testAttributes=true,
                                       bool testPoint=true)
{
  if(marktopreserve1!=LCC1::INVALID_MARK &&
     marktopreserve2!=LCC2::INVALID_MARK &&
     lcc1.is_marked(sd1, marktopreserve1)!=lcc2.is_marked(sd2, marktopreserve2))
  { return false; }

  bool match = true;

  // Two stacks used to run through the two maps.
  std::deque<typename LCC1::Dart_handle> toTreat1;
  std::deque<typename LCC2::Dart_handle> toTreat2;

   // A dart of this map is marked with m1 if its bijection was set
  // (and similarly for mark m2 and darts of map2)
  typename LCC1::size_type m1=lcc1.get_new_mark();
  typename LCC2::size_type m2=lcc2.get_new_mark();

  // A dart of this map is marked with markpush if it was already pushed
  // in the queue toTreat1.
  typename LCC1::size_type markpush=lcc1.get_new_mark();

  toTreat1.push_back(sd1);
  toTreat2.push_back(sd2);
  lcc1.mark(sd1, markpush);

  typename LCC1::Dart_handle dh1, other1;
  typename LCC2::Dart_handle dh2, other2;

  CGAL::Unique_hash_map<typename LCC1::Dart_handle,
                        typename LCC2::Dart_handle,
                        typename CGAL::Handle_hash_function> bijection;

  while (match && !toTreat1.empty())
  {
    // Next dart
    dh1=toTreat1.front();
    toTreat1.pop_front();
    dh2=toTreat2.front();
    toTreat2.pop_front();

    if(lcc1.is_marked(dh1, m1)!=lcc2.is_marked(dh2, m2))
    { match=false; }
    else if(!lcc1.is_marked(dh1, m1))
    {
      bijection[dh1]=dh2;
      lcc1.mark(dh1, m1);
      lcc2.mark(dh2, m2);

      // We first test info of darts
      if(match && testDartInfo)
      { match=CGAL::internal::template Test_is_same_dart_info_functor<LCC1, LCC2>::
            run(lcc1, lcc2, dh1, dh2); }

      // We need to test in both direction because
      // Foreach_enabled_attributes only test non void attributes
      // of Self. Functor Test_is_same_attribute_functor will modify
      // the value of match to false if attributes do not match
      if(testAttributes)
      {
        if (match)
        { LCC1::Helper::template Foreach_enabled_attributes
              <CGAL::internal::template Test_is_same_attribute_functor
              <LCC1, LCC2>>::run(lcc1, lcc2, dh1, dh2, match); }
        if (match)
        { LCC2::Helper::template Foreach_enabled_attributes
              <CGAL::internal::template Test_is_same_attribute_functor
              <LCC2, LCC1>>::run(lcc2, lcc1, dh2, dh1, match); }
      }

      if(match && testPoint)
      {
        // Only point of 0-attributes are tested. TODO test point of all
        // attributes ?
        match=CGAL::internal::template Test_is_same_attribute_point_functor
              <LCC1, LCC2, 0>::run(lcc1, lcc2, dh1, dh2);
      }

      if(match &&
         marktopreserve1!=LCC1::INVALID_MARK &&
         marktopreserve2!=LCC2::INVALID_MARK &&
         lcc1.is_marked(dh1, marktopreserve1)!=
         lcc2.is_marked(dh2, marktopreserve2))
      { match=false; }

      // We test if the injection is valid with its neighboors.
      // We go out as soon as it is not satisfied.
      // Process next then opposite
      other1=next1(dh1); other2=next2(dh2);
      for(int i:{0,1})
      {
        if(other1==lcc1.null_dart_handle && other2!=lcc2.null_dart_handle)
        { match=false; }
        else if(other1!=lcc1.null_dart_handle && other2==lcc2.null_dart_handle)
        { match=false; }
        else if(other1!=lcc1.null_dart_handle && other2!=lcc2.null_dart_handle)
        {
          if(lcc1.is_marked(other1, m1)!=lcc2.is_marked(other2, m2))
          { match=false; }
          else
          {
            if(!lcc1.is_marked(other1, m1))
            {
              if (!lcc1.is_marked(other1, markpush))
              {
                toTreat1.push_back(other1);
                toTreat2.push_back(other2);
                lcc1.mark(other1, markpush);
              }
            }
            else
            {
              if (bijection[other1]!=other2)
              { match=false; }
            }
          }
        }
        if(i==0) { other1=opposite1(dh1); other2=opposite2(dh2); }
      }
    }
  }

  // Here we test if both queue are empty
  if(!toTreat1.empty() || !toTreat2.empty())
  { match=false; }

  // Here we unmark all the marked darts.
  toTreat1.clear();
  toTreat2.clear();

  toTreat1.push_back(sd1);
  toTreat2.push_back(sd2);
  lcc1.unmark(sd1, markpush);

  while (!toTreat1.empty())
  {
    dh1=toTreat1.front();
    toTreat1.pop_front();
    dh2=toTreat2.front();
    toTreat2.pop_front();

    lcc1.unmark(dh1, m1);
    lcc2.unmark(dh2, m2);

    other1=next1(dh1); other2=next2(dh2);
    for(int i:{0,1})
    {
      if(other1!=lcc1.null_dart_handle && other2!=lcc2.null_dart_handle)
      {
        if (lcc1.is_marked(other1, markpush))
        {
          toTreat1.push_back(other1);
          toTreat2.push_back(other2);
          lcc1.unmark(other1, markpush);
        }
      }
      if(i==0) { other1=opposite1(dh1); other2=opposite2(dh2); }
    }
  }

  assert(lcc1.is_whole_map_unmarked(m1));
  assert(lcc1.is_whole_map_unmarked(markpush));
  assert(lcc2.is_whole_map_unmarked(m2));
  lcc1.free_mark(m1);
  lcc1.free_mark(markpush);
  lcc2.free_mark(m2);

  return match;
}

/// Test if surface(sd) is isomorphic with the border of the spattern
/// @return the dart of spattern satisfying the bijection, nullprt if there is
/// no isomorphism.
template<typename LCC1, typename LCC2>
typename LCC2::Dart_handle is_surface_isomorphic_to_spattern
(LCC1& lcc, typename LCC1::Dart_handle sd,
 LCC2& spattern,
 typename LCC2::size_type faceborder,
 typename LCC1::size_type marktopreserve1=LCC1::INVALID_MARK,
 typename LCC2::size_type marktopreserve2=LCC2::INVALID_MARK,
 bool testDartInfo=true,
 bool testAttributes=true,
 bool testPoint=true)
{
  typename LCC2::Dart_handle res=nullptr;
  for(auto it=spattern.darts().begin(), itend=spattern.darts().end();
      res==nullptr && it!=itend; ++it)
  {
    if(spattern.is_marked(it, faceborder))
    {
      if(are_volumes_isomorphic_from_darts(lcc, sd, spattern, it,
                                           marktopreserve1, marktopreserve2,

                                           [&lcc](typename LCC1::Dart_handle dh)
                                           -> typename LCC1::Dart_handle
                                           { return lcc.template beta<1>(dh); },
                                           [&lcc](typename LCC1::Dart_handle dh)
                                           -> typename LCC1::Dart_handle
                                           { return lcc.template beta<2>(dh); },

                                           [&spattern, faceborder](typename LCC2::Dart_handle dh)
                                           -> typename LCC2::Dart_handle
                                           { typename LCC2::Dart_handle other=
                                           spattern.template beta<1>(dh);
                                           while(!spattern.is_marked(other, faceborder))
                                           { other=spattern.template beta<2,1>(other); }
                                           return other;
                                           },

                                           [&spattern](typename LCC2::Dart_handle dh)
                                           -> typename LCC2::Dart_handle
                                            { return spattern.template beta<2>(dh); },

                                           testDartInfo,
                                           testAttributes,
                                           testPoint))
      { res=it; }
    }
  }
  return res;
}

/// Test if volume(sd) is isomorphic with the border of the vpattern
/// @return the dart of vpattern satisfying the bijection, nullprt if there is
/// no isomorphism.
template<typename LCC1, typename LCC2>
typename LCC2::Dart_handle is_volume_isomorphic_to_vpattern
(LCC1& lcc, typename LCC1::Dart_handle sd,
 LCC2& vpattern,
 typename LCC1::size_type marktopreserve1=LCC1::INVALID_MARK,
 typename LCC2::size_type marktopreserve2=LCC2::INVALID_MARK,
 bool testDartInfo=true,
 bool testAttributes=true,
 bool testPoint=true)
{
  typename LCC2::Dart_handle res=nullptr;
  for(auto it=vpattern.darts().begin(), itend=vpattern.darts().end();
      res==nullptr && it!=itend; ++it)
  {
    if(vpattern.template is_free<3>(it))
    {
      if(are_volumes_isomorphic_from_darts(lcc, sd, vpattern, it,
                                           marktopreserve1, marktopreserve2,

                                           [&lcc](typename LCC1::Dart_handle dh)
                                           -> typename LCC1::Dart_handle
                                           { return lcc.template beta<1>(dh); },
                                           [&lcc](typename LCC1::Dart_handle dh)
                                           -> typename LCC1::Dart_handle
                                           { return lcc.template beta<2>(dh); },

                                           [&vpattern](typename LCC2::Dart_handle dh)
                                           -> typename LCC2::Dart_handle
                                           { return vpattern.template beta<1>(dh); },
                                           [&vpattern](typename LCC2::Dart_handle dh)
                                           -> typename LCC2::Dart_handle
                                           { typename LCC2::Dart_handle other=
                                           vpattern.template beta<2>(dh);
                                           while(!vpattern.template is_free<3>(other))
                                           { other=vpattern.template beta<3,2>(other); }
                                           return other;
                                           },

                                         testDartInfo,
                                         testAttributes,
                                         testPoint))
      { res=it; }
    }
  }
  return res;
}

#endif // CMAP_ISOMORPHISMS_H
