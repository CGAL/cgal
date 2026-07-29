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
#ifndef CMAP_QUERY_REPLACE_H
#define CMAP_QUERY_REPLACE_H

#include <queue>
#include <utility>
#include <vector>

// #include "cmap_copy.h" // FOR DEBUGGING
// #include <CGAL/draw_linear_cell_complex.h> // FOR DEBUGGING
#include <CGAL/Linear_cell_complex/cmap_3close_cc.h>
#include <CGAL/Linear_cell_complex/query_replace/cmap_query_replace_geometry.h>
#include <CGAL/Linear_cell_complex/cmap_isomorphisms.h>

namespace CGAL::internal
{
////////////////////////////////////////////////////////////////////////////////
/// Compute the bijection between the external edges of the pattern and
/// the face isomorphic to this external boundary.
/// Mark the dart of the external faces of the pattern.
/// @input dh1 is a copy of a dart of the pattern into lcc
/// @input dh2 is a dart of the face into lcc
/// @pre the two elements must be isomorphic
template<typename LCC>
void compute_face_bijection_from_pattern_to_dart(LCC& lcc,
                                                 typename LCC::Dart_descriptor dh1,
                                                 typename LCC::Dart_descriptor dh2,
                                                 typename LCC::size_type markexternal,
                                                 Dart_mapping<LCC>& pattern_to_face)
{
  assert(lcc.template is_free<2>(dh1));
  pattern_to_face.clear();
  typename LCC::Dart_descriptor cur1=dh1;
  typename LCC::Dart_descriptor cur2=dh2;
  do
  {
    pattern_to_face[cur1]=cur2;
    lcc.mark(cur1, markexternal);
    cur1=lcc.template beta<1>(cur1);
    while(!lcc.template is_free<2>(cur1))
    { cur1=lcc.template beta<2,1>(cur1); }
    cur2=lcc.template beta<1>(cur2);
  }
  while(cur1!=dh1);
}
////////////////////////////////////////////////////////////////////////////////
/// Compute the bijection between the edges of the pattern and
/// the faces isomorphic to these external boundaries.
/// 2-unsew all darts marked face border.
/// @input dh1 is a dart of the pattern (WARNING: and NOT a copy in LCC contrary to similar methods for face and volume)
/// @input dh2 is a dart of the face
/// @pre the two elements must be isomorphic
template<typename LCC>
void compute_surface_bijection_from_pattern_to_dart(LCC& lcc,
                                                    Pattern<LCC, 2>& pattern,
                                                    typename LCC::Dart_descriptor dh1,
                                                    typename LCC::Dart_descriptor dh2,
                                                    Dart_mapping<LCC>&
                                                        pattern_to_global,
                                                    Dart_mapping<LCC>&
                                                        pattern_to_surface)
{
  assert(!pattern.lcc().template is_free<2>(dh1));
  assert(pattern.lcc().is_marked(dh1, pattern.mark_faceborder()));
  std::queue<std::pair<typename LCC::Dart_descriptor,
                       typename LCC::Dart_descriptor>> to_treat;
  typename LCC::size_type treated=pattern.lcc().get_new_mark();
  typename LCC::Dart_descriptor other1, other2;
  pattern_to_surface.clear();
  to_treat.push(std::make_pair(dh1, dh2));
  pattern.lcc().mark(dh1, treated);
  while(!to_treat.empty())
  {
    auto cur=to_treat.front();
    to_treat.pop();
    pattern_to_surface[pattern_to_global[cur.first]]=cur.second;

    // Process beta1
    other1=pattern.lcc().template beta<1>(cur.first);
    while(!pattern.lcc().is_marked(other1, pattern.mark_faceborder()))
    { other1=pattern.lcc().template beta<2,1>(other1); }

    other2=lcc.template beta<1>(cur.second);
    assert(other1!=lcc.null_descriptor && other2!=lcc.null_descriptor);
    if(!pattern.lcc().is_marked(other1, treated))
    {
      to_treat.push(std::make_pair(other1, other2));
      pattern.lcc().mark(other1, treated);
    }

    // Process beta2
    other1=pattern.lcc().template beta<2>(cur.first);
    other2=lcc.template beta<2>(cur.second);
    assert(other1!=lcc.null_descriptor && other2!=lcc.null_descriptor);
    if(!pattern.lcc().is_marked(other1, treated))
    {
      to_treat.push(std::make_pair(other1, other2));
      pattern.lcc().mark(other1, treated);
    }
  }

  for(auto it=pattern.lcc().darts().begin(),
       itend=pattern.lcc().darts().end(); it!=itend; ++it)
  {
    if(pattern.lcc().is_marked(it, pattern.mark_faceborder()))
    {
      pattern.lcc().unmark(it, treated);
      if(!lcc.template is_free<2>(pattern_to_global[it]))
      { lcc.template unsew<2>(pattern_to_global[it]); }
    }
    assert(!pattern.lcc().is_marked(it, treated));
  }
  pattern.lcc().free_mark(treated);
}
////////////////////////////////////////////////////////////////////////////////
/// Compute the bijection between the external faces of the pattern and
/// the faces of the volume isomorphic to these external boundaries.
/// Mark the dart of the external boundary of the pattern.
/// @input dh1 is a copy of a dart of the pattern into lcc
/// @input dh2 is a dart of the volume into lcc
/// @pre the two elements must be isomorphic
///
template<typename LCC>
void compute_volume_bijection_from_pattern_to_dart(LCC& lcc,
                                                   typename LCC::Dart_descriptor dh1,
                                                   typename LCC::Dart_descriptor dh2,
                                                   typename LCC::size_type markexternal,
                                                   Dart_mapping<LCC>&
                                                       pattern_to_volume)
{
  assert(lcc.template is_free<3>(dh1));
  std::queue<std::pair<typename LCC::Dart_descriptor, typename LCC::Dart_descriptor>> to_treat;
  typename LCC::Dart_descriptor other1, other2;
  pattern_to_volume.clear();
  to_treat.push(std::make_pair(dh1, dh2));
  lcc.mark(dh1, markexternal);
  while(!to_treat.empty())
  {
    auto cur=to_treat.front();
    to_treat.pop();
    assert(lcc.template is_free<3>(cur.first));
    pattern_to_volume[cur.first]=cur.second;

    // Process beta1
    other1=lcc.template beta<1>(cur.first);
    other2=lcc.template beta<1>(cur.second);
    assert(other1!=lcc.null_descriptor && other2!=lcc.null_descriptor);
    if(!lcc.is_marked(other1, markexternal))
    {
      to_treat.push(std::make_pair(other1, other2));
      lcc.mark(other1, markexternal);
    }

    // Process beta2
    other1=lcc.template beta<2>(cur.first);
    while(!lcc.template is_free<3>(other1))
    { other1=lcc.template beta<3,2>(other1); }
    other2=lcc.template beta<2>(cur.second);
    assert(other1!=lcc.null_descriptor && other2!=lcc.null_descriptor);
    if(!lcc.is_marked(other1, markexternal))
    {
      to_treat.push(std::make_pair(other1, other2));
      lcc.mark(other1, markexternal);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/// Replace volume(dh1) by the vpattern, knowing that the surface of
/// vpattern is isomorphic with volume(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the surface of vpattern is isomorphic with volume(dh1)
template<typename LCC>
void replace_one_volume_from_dart(LCC& lcc,
                                  typename LCC::Dart_descriptor dh1,
                                  Pattern<LCC, 3>& vpattern,
                                  typename LCC::Dart_descriptor dh2)
{
  assert(is_volume_isomorphic_to_vpattern_from_dart(lcc, dh1,
                                                    vpattern.lcc(), dh2,
                                                    LCC::INVALID_MARK,
                                                    LCC::INVALID_MARK,
                                                    false, false, false));

  Dart_mapping<LCC> pattern_to_global;
  Dart_mapping<LCC> links_from_pattern_to_volume;
  auto amark=lcc.get_new_mark();
  // 1) Copy pattern into lcc. New darts are not marked
  lcc.copy(vpattern.lcc(), &pattern_to_global);
  // 2) Compute old_3sew to store 3-links of darts in volume(dh2)
  assert(vpattern.lcc().darts().owns(dh2));
  compute_volume_bijection_from_pattern_to_dart
      (lcc, pattern_to_global[dh2], dh1, amark,
       links_from_pattern_to_volume);
  transform_geometry_of_vpattern(lcc, amark, links_from_pattern_to_volume,
                                 pattern_to_global, vpattern);

  // 3) Remove all the external faces of the copy of the pattern, and 2-sew
  //    the internal faces of the copy of the pattern with the boundary of
  //    the volume.
  std::vector<std::pair<typename LCC::Dart_descriptor,
                        typename LCC::Dart_descriptor>> tosew;
  typename LCC::Dart_descriptor otherdh;
  tosew.reserve(links_from_pattern_to_volume.size());
  for(auto curdh: links_from_pattern_to_volume)
  {
    otherdh=lcc.template beta<2>(curdh.first);
    if(lcc.is_dart_used(otherdh) && !lcc.is_marked(otherdh, amark))
    {
      lcc.template topo_unsew<2>(curdh.first);
      if(!lcc.template is_free<2>(curdh.second))
      { lcc.template unsew<2>(curdh.second); }
      tosew.push_back(std::make_pair(curdh.second, otherdh));
    }
  }
  for(auto curdh: tosew)
  { lcc.template sew<2>(curdh.first, curdh.second); }

  for(auto curdh: links_from_pattern_to_volume)
  { lcc.erase_dart(curdh.first); }
  // assert(lcc.is_valid());
  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
}
////////////////////////////////////////////////////////////////////////////////
/// Replace face(dh1) by the fpattern, knowing that the border of
/// fpattern is isomorphic with face(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the border of fpattern is isomorphic with face(dh1)
template<typename LCC>
void replace_one_face_from_dart(LCC& lcc,
                                typename LCC::Dart_descriptor dh1,
                                Pattern<LCC, 1>& fpattern,
                                typename LCC::Dart_descriptor dh2)
{
  Dart_mapping<LCC> pattern_to_global;
  Dart_mapping<LCC> links_from_pattern_to_face;
  bool with_beta3=false;
  // 1) Copy pattern into lcc.
  lcc.copy(fpattern.lcc(), &pattern_to_global);
  if(!lcc.template is_free<3>(dh1))
  {
    close_cc_for_beta3(lcc, pattern_to_global[dh2]);
    with_beta3=true;
  }
  // 2) Compute mapping from the boundary of the pattern and
  //    the face isomorphic to this external boundary
  auto amark=lcc.get_new_mark();
  compute_face_bijection_from_pattern_to_dart
      (lcc, pattern_to_global[dh2], dh1, amark, links_from_pattern_to_face);
  transform_geometry_of_fpattern(lcc, links_from_pattern_to_face,
                                 pattern_to_global, fpattern);

  // 3) Remove all the external edges of the copy of the pattern, and 1-sew
  //    the internal edges of the copy of the pattern with the boundary of
  //    the face.
  std::vector<std::pair<typename LCC::Dart_descriptor,
                        typename LCC::Dart_descriptor>> tosew0, tosew1;
  tosew0.reserve(links_from_pattern_to_face.size());
  tosew1.reserve(links_from_pattern_to_face.size());
  typename LCC::Dart_descriptor otherdh;
  for(auto curdh: links_from_pattern_to_face)
  {
    otherdh=lcc.template beta<0>(curdh.first);
    if(lcc.is_dart_used(otherdh) && !lcc.is_marked(otherdh, amark))
    {
      lcc.template topo_unsew<0>(curdh.first);
      //lcc.template unsew<0>(curdh.first);
      if(!lcc.template is_free<0>(curdh.second))
      { lcc.template unsew<0>(curdh.second); }
      //lcc.template sew<0>(curdh.second, otherdh);
      tosew0.push_back(std::make_pair(curdh.second, otherdh));
    }
    otherdh=lcc.template beta<1>(curdh.first);
    if(lcc.is_dart_used(otherdh) && !lcc.is_marked(otherdh, amark))
    {
      lcc.template topo_unsew<1>(curdh.first);
      //lcc.template unsew<1>(curdh.first);
      if(!lcc.template is_free<1>(curdh.second))
      { lcc.template unsew<1>(curdh.second); }
      //lcc.template sew<1>(curdh.second, otherdh);
      tosew1.push_back(std::make_pair(curdh.second, otherdh));
    }
  }
  for(auto curdh: tosew0)
  { lcc.template sew<0>(curdh.first, curdh.second); }
  for(auto curdh: tosew1)
  { lcc.template sew<1>(curdh.first, curdh.second); }
  for(auto curdh: links_from_pattern_to_face)
  {
    if(with_beta3)
    { lcc.erase_dart(lcc.template beta<3>(curdh.first)); }
    lcc.erase_dart(curdh.first);
  }
  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
  // assert(lcc.is_valid());
}
////////////////////////////////////////////////////////////////////////////////
/// Replace surface(dh1) by the spattern, knowing that the border of
/// spattern is isomorphic with surface(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the border of spattern is isomorphic with surface(dh1)
template<typename LCC>
void replace_one_surface_from_dart(LCC& lcc,
                                   typename LCC::Dart_descriptor dh1,
                                   Pattern<LCC, 2>& spattern,
                                   typename LCC::Dart_descriptor dh2)
{
  Dart_mapping<LCC> pattern_to_global;
  Dart_mapping<LCC> links_from_pattern_to_surface;
  // 1) Copy pattern into lcc.
  lcc.copy(spattern.lcc(), &pattern_to_global);

  // 2) Compute mapping from the boundary of the pattern and
  //    the surface isomorphic to this external boundary, and 2-unsew
  //    each face border of the pattern
  auto amark=lcc.get_new_mark();
  compute_surface_bijection_from_pattern_to_dart(lcc, spattern,
                                                 dh2, dh1,
                                                 pattern_to_global,
                                                 links_from_pattern_to_surface);

  // Transform the geometry of all faces (same method than for faces)
  transform_geometry_of_spattern(lcc, links_from_pattern_to_surface,
                                 pattern_to_global, spattern);

  // 3) Remove all the external edges of the copy of the pattern, and 1-sew
  //    the internal edges of the copy of the pattern with the boundary of
  //    the face.
  for(auto curdh: links_from_pattern_to_surface)
  {
    if(!lcc.template is_free<3>(curdh.second) &&
        lcc.template is_free<3>(curdh.first))
    { close_cc_for_beta3(lcc, curdh.first); }

    dh2=lcc.template beta<0>(curdh.first);
    if(lcc.is_dart_used(dh2) && !lcc.is_marked(dh2, amark))
    {
      lcc.template unsew<0>(curdh.first);
      if(!lcc.template is_free<0>(curdh.second))
      { lcc.template unsew<0>(curdh.second); }
      lcc.template sew<0>(curdh.second, dh2);
    }
    dh2=lcc.template beta<1>(curdh.first);
    if(lcc.is_dart_used(dh2) && !lcc.is_marked(dh2, amark))
    {
      lcc.template unsew<1>(curdh.first);
      if(!lcc.template is_free<1>(curdh.second))
      { lcc.template unsew<1>(curdh.second); }
      lcc.template sew<1>(curdh.second, dh2);
    }
  }
  for(auto curdh: links_from_pattern_to_surface)
  {
    if(!lcc.template is_free<3>(curdh.first))
    { lcc.erase_dart(lcc.template beta<3>(curdh.first)); }
    lcc.erase_dart(curdh.first);
  }
  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
}
////////////////////////////////////////////////////////////////////////////////
}
#endif // CMAP_QUERY_REPLACE_H
