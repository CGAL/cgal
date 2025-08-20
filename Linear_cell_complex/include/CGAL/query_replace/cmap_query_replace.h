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
#ifndef CMAP_QUERY_REPLACE_H
#define CMAP_QUERY_REPLACE_H

#include <filesystem>
#include <limits>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cmap_3close_cc.h"
#include "cmap_query_replace_geometry.h"
#include "cmap_isomorphisms.h"
#include "cmap_signature.h"
#include "lcc_read_depending_extension.h"

///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern_substituer
{
public:
  using Dart_handle=typename LCC::Dart_handle;
  using size_type=typename LCC::size_type;

  using Signature_mapping=std::unordered_map<Signature,
  std::pair<Dart_handle, std::size_t>>;

  template<unsigned int type> // type==1 for face, 2 for surface and 3 for volume
  using Pattern_set=std::vector<Pattern<LCC, type>>;

  std::size_t number_of_fpatterns() const
  { return m_fpatterns.size(); }
  std::size_t number_of_spatterns() const
  { return m_spatterns.size(); }
  std::size_t number_of_vpatterns() const
  { return m_vpatterns.size(); }

  LCC& fpattern(std::size_t i)
  { return m_fpatterns[i].lcc(); }
  LCC& spattern(std::size_t i)
  { return m_spatterns[i].lcc(); }
  LCC& vpattern(std::size_t i)
  { return m_vpatterns[i].lcc(); }

   Signature_mapping& fsignatures()
  { return m_fsignatures; }
   Signature_mapping& ssignatures()
  { return m_ssignatures; }
   Signature_mapping& vsignatures()
  { return m_vsignatures; }

  typename Signature_mapping::const_iterator find_fpattern(const Signature& s) const
  { return m_fsignatures.find(s); }
  typename Signature_mapping::const_iterator find_spattern(const Signature& s) const
  { return m_ssignatures.find(s); }
  typename Signature_mapping::const_iterator find_vpattern(const Signature& s) const
  { return m_vsignatures.find(s); }

  typename Signature_mapping::const_iterator fpattern_end() const
  { return m_fsignatures.end(); }
  typename Signature_mapping::const_iterator spattern_end() const
  { return m_ssignatures.end(); }
  typename Signature_mapping::const_iterator vpattern_end() const
  { return m_vsignatures.end(); }

  void load_fpatterns(const std::string& directory_name,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    load_all_patterns<1>(directory_name, m_fpatterns);
    Signature signature;
    Dart_handle dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    std::size_t nb=0;
    m_fsignatures.clear();
    for(auto& pattern: m_fpatterns)
    {
      if(init_topreserve!=nullptr) // true iff the std::function is not empty
      {
        mark_to_preserve=pattern.reserve_mark_to_preserve();
        init_topreserve(pattern.lcc(), mark_to_preserve);
      }
      dh=fsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
      auto res=m_fsignatures.find(signature);
      if(res==m_fsignatures.end())
      {
        pattern.compute_barycentric_coord();
        m_fsignatures[signature]=std::make_pair(dh, nb);
      }
      else
      {
        std::cout<<"[ERROR] load_fpatterns: two patterns have same signature "
                 <<nb<<" and "<<res->second.second<<std::endl;
      }
      ++nb;
      // std::cout<<"[Pattern] Signature "<<nb<<": "; print_signature(signature);
    }
  }

  void load_additional_fpattern(const std::string& file_name,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    auto [success, id] = load_one_additional_pattern<1>(file_name, m_fpatterns);
    if (!success) {
      std::cerr << "load_additional_fpattern: file not found or format not readable" << std::endl;
      return;
    };

    Signature signature;
    Dart_handle dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;

    auto& pattern = m_fpatterns[id];

    if(init_topreserve!=nullptr) // true iff the std::function is not empty
    {
      mark_to_preserve=pattern.reserve_mark_to_preserve();
      init_topreserve(pattern.lcc(), mark_to_preserve);
    }
    dh=fsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
    auto res=m_fsignatures.find(signature);
    if(res==m_fsignatures.end())
    {
      pattern.compute_barycentric_coord();
      m_fsignatures[signature]=std::make_pair(dh, id);
    }
    else
    {
      std::cout<<"[ERROR] load_fpatterns: two patterns have same signature "
                <<id<<" and "<<res->second.second<<std::endl;
    }
  }

  void load_spatterns(const std::string& directory_name,
                      std::function<void(LCC&, size_type)> init_faceborder,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    load_all_patterns<2>(directory_name, m_spatterns);
    Signature signature;
    Dart_handle dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    std::size_t nb=0;
    m_ssignatures.clear();
    for(auto& pattern: m_spatterns)
    {
      init_faceborder(pattern.lcc(), pattern.m_mark_faceborder);
      if(init_topreserve!=nullptr) // true iff the std::function is not empty
      {
        mark_to_preserve=pattern.reserve_mark_to_preserve();
        init_topreserve(pattern.lcc(), mark_to_preserve);
      }
      dh=ssignature_of_pattern(pattern.lcc(), pattern.m_mark_faceborder,
                               mark_to_preserve, signature, false);
      auto res=m_ssignatures.find(signature);
      if(res==m_ssignatures.end())
      {
        pattern.compute_barycentric_coord();
        assert(pattern.lcc().is_marked(dh, pattern.m_mark_faceborder));
        m_ssignatures[signature]=std::make_pair(dh, nb);
      }
      else
      {
        std::cout<<"[ERROR] load_spatterns: two patterns have same signature "
                 <<nb<<" and "<<res->second.second<<std::endl;
      }
      ++nb;
      // std::cout<<"[Pattern] Signature "<<nb<<": "; print_signature(signature);
    }
  }
  void load_vpatterns(const std::string& directory_name,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    load_all_patterns<3>(directory_name, m_vpatterns);
    Signature signature;
    Dart_handle dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    std::size_t nb=0;
    m_vsignatures.clear();
    for(auto& pattern: m_vpatterns)
    {
      if(init_topreserve!=nullptr) // true iff the std::function is not empty
      {
        mark_to_preserve=pattern.reserve_mark_to_preserve();
        init_topreserve(pattern.lcc(), mark_to_preserve);
      }
      dh=vsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
      auto res=m_vsignatures.find(signature);
      if(res==m_vsignatures.end())
      {
        pattern.compute_barycentric_coord();
        m_vsignatures[signature]=std::make_pair(dh, nb);
      }
      else
      {
        std::cout<<"[ERROR] load_vpatterns: two patterns have same signature "
                 <<nb<<" and "<<res->second.second<<std::endl;
      }
      ++nb;
      // std::cout<<"[Pattern] Signature "<<nb<<": "; print_signature(signature);
    }
  }

  void load_additional_vpattern(const std::string& directory_name,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    auto [success, id] = load_one_additional_pattern<3>(directory_name, m_vpatterns);
    if (!success) {
      std::cerr << "load_additional_vpattern: file not found or format not readable" << std::endl;
      return;
    };

    Signature signature;
    Dart_handle dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    auto& pattern = m_vpatterns[id];


    if(init_topreserve!=nullptr) // true iff the std::function is not empty
    {
      mark_to_preserve=pattern.reserve_mark_to_preserve();
      init_topreserve(pattern.lcc(), mark_to_preserve);
    }
    dh=vsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
    auto res=m_vsignatures.find(signature);
    if(res==m_vsignatures.end())
    {
      pattern.compute_barycentric_coord();
      m_vsignatures[signature]=std::make_pair(dh, id);
    }
    else
    {
      std::cout<<"[ERROR] load_vpatterns: two patterns have same signature "
                <<id<<" and "<<res->second.second<<std::endl;
    }
  }

protected:
  template<unsigned int type>
  void load_all_patterns(const std::string& directory_name,
                         Pattern_set<type>& patterns)
  {
    patterns.clear();
    std::size_t nb=0;
    const std::filesystem::path dir(directory_name);
    if(!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir))
    { return; }

    for(auto const& dir_entry: std::filesystem::directory_iterator{dir})
    {
      if(dir_entry.is_regular_file() &&
         is_lcc_readable_file(dir_entry.path().string()))
      { ++nb; }
    }

    patterns.resize(nb);
    nb=0;
    // std::cout<<"##############################"<<std::endl;
    for(auto const& dir_entry: std::filesystem::directory_iterator{dir})
    {
      if(dir_entry.is_regular_file() &&
         is_lcc_readable_file(dir_entry.path().string()))
      {
        // std::cout<<"pattern "<<nb<<": "<<dir_entry.path().string()<<std::endl;
        read_depending_extension(dir_entry.path().string(),
                                 patterns[nb].lcc());
        ++nb;
      }
    }
  }

  // Returns the id of the loaded pattern, -1 if it couldn't be loaded
  template<unsigned int type>
  std::pair<bool, std::size_t> load_one_additional_pattern(const std::string& file_name,
                         Pattern_set<type>& patterns)
  {
    const std::filesystem::path file(file_name);
    if (!std::filesystem::exists(file)
      || !std::filesystem::is_regular_file(file)
      || !is_lcc_readable_file(file.string()))
    { return {false, 0}; }

    std::size_t id = patterns.size();
    patterns.push_back(Pattern<LCC,type>());

    read_depending_extension(file.string(),
                                 patterns[id].lcc());

    return {true, id};
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// Compute the bijection between the external edges of the pattern and
  /// the face isomorphic to this external boundary.
  /// Mark the dart of the external faces of the pattern.
  /// @input dh1 is a copy of a dart of the pattern into lcc
  /// @input dh2 is a dart of the face into lcc
  /// @pre the two elements must be isomorphic
  void compute_face_bijection_from_pattern_to_dart(LCC& lcc,
                                                   Dart_handle dh1,
                                                   Dart_handle dh2,
                                                   size_type markexternal,
                                                   Dart_mapping<LCC>&
                                                   pattern_to_face)
  {
    assert(lcc.template is_free<2>(dh1));
    pattern_to_face.clear();
    Dart_handle cur1=dh1;
    Dart_handle cur2=dh2;
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
  void compute_surface_bijection_from_pattern_to_dart(LCC& lcc,
                                                      Pattern<LCC, 2>& pattern,
                                                      Dart_handle dh1,
                                                      Dart_handle dh2,
                                                      Dart_mapping<LCC>&
                                                      pattern_to_global,
                                                      Dart_mapping<LCC>&
                                                      pattern_to_surface)
  {
    assert(!pattern.lcc().template is_free<2>(dh1));
    assert(pattern.lcc().is_marked(dh1, pattern.m_mark_faceborder));
    std::queue<std::pair<Dart_handle, Dart_handle>> to_treat;
    size_type treated=pattern.lcc().get_new_mark();
    Dart_handle other1, other2;
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
      while(!pattern.lcc().is_marked(other1, pattern.m_mark_faceborder))
      { other1=pattern.lcc().template beta<2,1>(other1); }

      other2=lcc.template beta<1>(cur.second);
      assert(other1!=lcc.null_handle && other2!=lcc.null_handle);
      if(!pattern.lcc().is_marked(other1, treated))
      {
        to_treat.push(std::make_pair(other1, other2));
        pattern.lcc().mark(other1, treated);
      }

      // Process beta2
      other1=pattern.lcc().template beta<2>(cur.first);
      other2=lcc.template beta<2>(cur.second);
      assert(other1!=lcc.null_handle && other2!=lcc.null_handle);
      if(!pattern.lcc().is_marked(other1, treated))
      {
        to_treat.push(std::make_pair(other1, other2));
        pattern.lcc().mark(other1, treated);
      }
    }

    for(auto it=pattern.lcc().darts().begin(),
        itend=pattern.lcc().darts().end(); it!=itend; ++it)
    {
      if(pattern.lcc().is_marked(it, pattern.m_mark_faceborder))
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
  void compute_volume_bijection_from_pattern_to_dart(LCC& lcc,
                                                     Dart_handle dh1,
                                                     Dart_handle dh2,
                                                     size_type markexternal,
                                                     Dart_mapping<LCC>&
                                                     pattern_to_volume)
  {
    assert(lcc.template is_free<3>(dh1));
    std::queue<std::pair<Dart_handle, Dart_handle>> to_treat;
    Dart_handle other1, other2;
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
      assert(other1!=lcc.null_handle && other2!=lcc.null_handle);
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
      assert(other1!=lcc.null_handle && other2!=lcc.null_handle);
      if(!lcc.is_marked(other1, markexternal))
      {
        to_treat.push(std::make_pair(other1, other2));
        lcc.mark(other1, markexternal);
      }
    }
  }

public:

////////////////////////////////////////////////////////////////////////////////
/// Replace volume(dh1) by the vpattern, knowing that the surface of
/// vpattern is isomorphic with volume(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the surface of vpattern is isomorphic with volume(dh1)
void replace_one_volume_from_dart(LCC& lcc,
                                  Dart_handle dh1,
                                  Pattern<LCC, 3>& vpattern,
                                  Dart_handle dh2)
{
  Dart_mapping<LCC> pattern_to_global;
  Dart_mapping<LCC> links_from_pattern_to_volume;
  auto amark=lcc.get_new_mark();
  // 1) Copy pattern into lcc. New darts are not be marked
  lcc.copy(vpattern.lcc(), &pattern_to_global);
  // 2) Compute old_3sew to store 3-links of darts in volume(dh2)
  assert(vpattern.lcc().darts().owns(dh2));
  compute_volume_bijection_from_pattern_to_dart
      (lcc, pattern_to_global[dh2], dh1, amark,
      links_from_pattern_to_volume);
  transform_geometry_of_vpattern(lcc, links_from_pattern_to_volume,
                                 pattern_to_global, vpattern);

  // 3) Remove all the external faces of the copy of the pattern, and 2-sew
  //    the internal faces of the copy of the pattern with the boundary of
  //    the volume.
  std::vector<std::pair<Dart_handle, Dart_handle>> tosew;
  Dart_handle otherdh;
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
std::size_t query_replace_one_volume_from_dart(LCC& lcc,
                                               Dart_handle dh,
                                               size_type marktopreserve)
{
  std::size_t replaced=std::numeric_limits<std::size_t>::max();
  Signature word_signature;
  Dart_handle
      dh2=vsignature_of_volume_for_dart(lcc, dh, marktopreserve, word_signature); //, true);

  if (dh2==nullptr) { return replaced; }

  Signature signature;
  vsignature_of_volume(lcc, dh, marktopreserve, signature);
  if(signature!=word_signature) { return replaced; }

  //std::cout<<"Source: "; print_signature(signature);
  auto res=m_vsignatures.find(signature);
  if(res!=m_vsignatures.end())
  {
    replace_one_volume_from_dart(lcc, dh2, m_vpatterns[res->second.second],
                                 res->second.first);
    replaced=res->second.second;
  }
  // else { std::cout<<"NOT found"<<std::endl; }
  return replaced;
}

bool replace_one_volume_from_signature(LCC& lcc,
                                Dart_handle dh1,
                                Signature& signature,
                                Dart_handle dh2)
{
  //std::cout<<"Source: "; print_signature(signature);
  std::size_t replaced=std::numeric_limits<std::size_t>::max();
  auto res=m_vsignatures.find(signature);
  if(res!=m_vsignatures.end())
  {
    // std::cout<<"FOUND Pattern "<<res->second.second+1<<std::endl;
    replace_one_volume_from_dart(lcc, dh2, m_vpatterns[res->second.second],
        res->second.first);
    replaced=res->second.second;
  }
  // else { std::cout<<"face NOT found"<<std::endl; }
  return replaced;
}

////////////////////////////////////////////////////////////////////////////////
/// Query volume(dh) in the set of patterns, and if one pattern matches,
/// replace volume(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_volume(LCC& lcc,
                                     Dart_handle dh,
                                     size_type marktopreserve=LCC::INVALID_MARK)
{
  Signature signature;
  Dart_handle
      dh2=vsignature_of_volume(lcc, dh, marktopreserve, signature); //, true);
  return replace_one_volume_from_signature(lcc, dh, signature, dh2);
}


////////////////////////////////////////////////////////////////////////////////
/// Query volume(dh) but without using signatures. If one pattern matches,
/// replace volume(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_volume_without_signature(LCC& lcc,
                                                       Dart_handle dh,
                                                       size_type marktopreserve=LCC::INVALID_MARK)
{
  Dart_handle res=nullptr, sd=dh;

  if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
  {
    auto it=lcc.template darts_of_cell<3>(dh).begin(),
          itend=lcc.template darts_of_cell<3>(dh).end();
    while(it!=itend && !lcc.is_marked(it, marktopreserve))
    { ++it; }
    if(it!=itend) { sd=it; }
  }

  std::size_t i=0;
  while(res==nullptr && i<number_of_vpatterns())
  {
    res=is_volume_isomorphic_to_vpattern(lcc, sd, vpattern(i), marktopreserve,
                                         m_vpatterns[i].mark_to_preserve(),
                                         false, false, false);
    if(res==nullptr) { ++i; }
  }

  if(res!=nullptr)
  { replace_one_volume_from_dart(lcc, sd, m_vpatterns[i], res); }
  else
  { i=std::numeric_limits<std::size_t>::max(); }
  // std::cout<<"NOT found"<<std::endl;
  return i;
}
////////////////////////////////////////////////////////////////////////////////
/// Replace face(dh1) by the fpattern, knowing that the border of
/// fpattern is isomorphic with face(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the border of fpattern is isomorphic with face(dh1)
void replace_one_face_from_dart(LCC& lcc,
                                Dart_handle dh1,
                                Pattern<LCC, 1>& fpattern,
                                Dart_handle dh2)
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
  std::vector<std::pair<Dart_handle, Dart_handle>> tosew0, tosew1;
  tosew0.reserve(links_from_pattern_to_face.size());
  tosew1.reserve(links_from_pattern_to_face.size());
  Dart_handle otherdh;
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
/// Query face(dh) in the set of patterns, and if one pattern matches,
/// replace face(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_face(LCC& lcc,
                                   Dart_handle dh,
                                   size_type marktopreserve=LCC::INVALID_MARK)
{
  Signature signature;
  Dart_handle dh2=fsignature_of_face(lcc, dh, marktopreserve, signature); //, true);
  return replace_one_face_from_signature(lcc, dh, signature, dh2);
}

bool replace_one_face_from_signature(LCC& lcc,
                                Dart_handle dh1,
                                Signature& signature,
                                Dart_handle dh2)
{
  //std::cout<<"Source: "; print_signature(signature);
  std::size_t replaced=std::numeric_limits<std::size_t>::max();
  auto res=m_fsignatures.find(signature);
  if(res!=m_fsignatures.end())
  {
    // std::cout<<"FOUND Pattern "<<res->second.second+1<<std::endl;
    replace_one_face_from_dart(lcc, dh2, m_fpatterns[res->second.second],
        res->second.first);
    replaced=res->second.second;
  }
  // else { std::cout<<"face NOT found"<<std::endl; }
  return replaced;
}

////////////////////////////////////////////////////////////////////////////////
/// Query face(dh) but without using signatures. If one pattern matches,
/// replace face(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_face_without_signature(LCC& lcc,
                                                     Dart_handle dh,
                                                     size_type marktopreserve=LCC::INVALID_MARK)
{
  Dart_handle res=nullptr, sd=dh;

  if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
  {
    auto it=lcc.template darts_of_cell<2,2>(dh).begin(),
          itend=lcc.template darts_of_cell<2,2>(dh).end();
    while(it!=itend && !lcc.is_marked(it, marktopreserve))
    { ++it; }
    if(it!=itend) { sd=it; }
  }

  std::size_t i=0;
  while(res==nullptr && i<number_of_fpatterns())
  {
    res=is_face_isomorphic_to_fpattern(lcc, sd, fpattern(i), marktopreserve,
                                       m_fpatterns[i].mark_to_preserve(),
                                       false, false, false);
    if(res==nullptr) { ++i; }
  }

  if(res!=nullptr)
  { replace_one_face_from_dart(lcc, sd, m_fpatterns[i], res); }
  else
  { i=std::numeric_limits<std::size_t>::max(); }
  // std::cout<<"NOT found"<<std::endl;
  return i;
}
////////////////////////////////////////////////////////////////////////////////
/// Replace surface(dh1) by the spattern, knowing that the border of
/// spattern is isomorphic with surface(dh1) starting from the pair of darts
/// (dh1, dh2).
/// @pre the border of spattern is isomorphic with surface(dh1)
void replace_one_surface_from_dart(LCC& lcc,
                                   Dart_handle dh1,
                                   Pattern<LCC, 2>& spattern,
                                   Dart_handle dh2)
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
/// Query surface(dh) in the set of patterns, and if one pattern matches,
/// replace surface(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_surface(LCC& lcc,
                                      Dart_handle dh,
                                      size_type marktopreserve=LCC::INVALID_MARK)
{
  Signature signature;
  Dart_handle
      dh2=ssignature_of_surface(lcc, dh, marktopreserve, signature); //, true);
  typename LCC::Vector v1, v2;
  // std::cout<<"Source: "; print_signature(signature);
  std::size_t replaced=std::numeric_limits<std::size_t>::max();
  auto res=m_ssignatures.find(signature);
  if(res!=m_ssignatures.end())
  {
    // std::cout<<"FOUND Pattern "<<res->second.second+1<<std::endl;
    replace_one_surface_from_dart(lcc, dh2, m_spatterns[res->second.second],
        res->second.first);
    replaced=res->second.second;
    // assert(lcc.is_valid());
  }
  // else { std::cout<<"NOT found"<<std::endl; }
  return replaced;
}
////////////////////////////////////////////////////////////////////////////////
/// Query surface(dh) but without using signatures. If one pattern matches,
/// replace surface(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_surface_without_signature(LCC& lcc,
                                                        Dart_handle dh,
                                                        size_type marktopreserve=LCC::INVALID_MARK)
{
  Dart_handle res=nullptr, sd=dh;

  if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
  {
    auto it=lcc.template darts_of_cell<3>(dh).begin(),
          itend=lcc.template darts_of_cell<3>(dh).end();
    while(it!=itend && !lcc.is_marked(it, marktopreserve))
    { ++it; }
    if(it!=itend) { sd=it; }
  }

  std::size_t i=0;
  while(res==nullptr && i<number_of_spatterns())
  {
    res=is_surface_isomorphic_to_spattern(lcc, sd, spattern(i),
                                          m_spatterns[i].m_mark_faceborder,
                                          marktopreserve,
                                          m_spatterns[i].mark_to_preserve(),
                                          false, false, false);
    if(res==nullptr) { ++i; }
  }

  if(res!=nullptr)
  { replace_one_surface_from_dart(lcc, sd, m_spatterns[i], res); }
  else
  { i=std::numeric_limits<std::size_t>::max(); }
  // std::cout<<"NOT found"<<std::endl;
  return i;
}
////////////////////////////////////////////////////////////////////////////////
std::size_t replace_vpatterns(LCC& lcc,
                              size_type marktopreserve,
                              bool nosignature=false,
                              bool all=true,
                              bool trace=false)
{
  auto amark=lcc.get_new_mark();
  lcc.negate_mark(amark); // All darts are marked
  std::size_t res=0;
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      lcc.template unmark_cell<3>(it, amark);
      // New darts will not be marked
      std::size_t replaced=
          (nosignature?query_replace_one_volume_without_signature
                       (lcc, it, marktopreserve):
           query_replace_one_volume(lcc, it, marktopreserve));
      if(replaced!=std::numeric_limits<std::size_t>::max())
      {
        ++res;
        if(!all)
        {
          lcc.free_mark(amark);
          return true;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  lcc.free_mark(amark);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
std::size_t replace_spatterns(LCC& lcc,
                              size_type marktopreserve,
                              bool nosignature=false,
                              bool all=true,
                              bool trace=false)
{
  auto amark=lcc.get_new_mark();
  lcc.negate_mark(amark); // All darts are marked
  std::size_t res=0;
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      lcc.template unmark_cell<3>(it, amark);
      // New darts will not be marked
      std::size_t replaced=
          (nosignature?query_replace_one_surface_without_signature
                       (lcc, it, marktopreserve):
           query_replace_one_surface(lcc, it, marktopreserve));
      if(replaced!=std::numeric_limits<std::size_t>::max())
      {
        ++res;
        if(!all)
        {
          lcc.free_mark(amark);
          return true;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  lcc.free_mark(amark);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
std::size_t replace_fpatterns(LCC& lcc,
                              size_type marktopreserve,
                              bool nosignature=false,
                              bool all=true,
                              bool trace=false)
{
  auto amark=lcc.get_new_mark();
  lcc.negate_mark(amark); // All darts are marked
  std::size_t res=0;
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, amark))
    {
      lcc.template unmark_cell<2>(it, amark);
      // New darts will not be marked
      std::size_t replaced=
          (nosignature?query_replace_one_face_without_signature
                       (lcc, it, marktopreserve):
           query_replace_one_face(lcc, it, marktopreserve));
      if(replaced!=std::numeric_limits<std::size_t>::max())
      {
        ++res;
        if(!all)
        {
          lcc.free_mark(amark);
          return true;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  lcc.free_mark(amark);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
void generate_all_face_replacement(LCC& lcc,
                                   size_type marktopreserve,
                                   std::list<LCC>& reslccs)
{
  std::list<LCC> totreat; // list to avoid copy of LCC
  totreat.push_back(LCC());
  totreat.back()=lcc; // copy
  while(!totreat.empty())
  {
    LCC current=std::move(totreat.front());
    totreat.pop_front();

    bool replaced=false;
    for(auto it=current.darts().begin(), itend=current.darts().end();
        !replaced && it!=itend; ++it)
    {
      Signature signature;
      fsignature_of_face_for_dart(current, it, marktopreserve, signature); //, true);
      auto res=m_fsignatures.find(signature);
      if(res!=m_fsignatures.end())
      {
        Dart_handle dh2=it;
        do
        {
          totreat.push_back(LCC());
          std::unordered_map<Dart_handle, Dart_handle> origin_to_copy;
          totreat.back().copy(current, &origin_to_copy, nullptr);
          replace_one_face_from_dart(totreat.back(), origin_to_copy[dh2],
                                     m_fpatterns[res->second.second],
                                     res->second.first);
          replaced=true;
          do
          {
            dh2=current.template beta<1>(dh2);
            fsignature_of_face_for_dart(current, dh2, marktopreserve, signature); //, true);
            res=m_fsignatures.find(signature);
          }
          while(res==m_fsignatures.end() && dh2!=it);
        }
        while(dh2!=it);
      }
    }
    if(!replaced)
    {
      reslccs.push_back(LCC());
      current.swap(reslccs.back());
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void generate_all_surface_replacement(LCC& lcc,
                                      size_type marktopreserve,
                                      std::list<LCC>& reslccs)
{
  std::list<LCC> totreat; // list to avoid copy of LCC
  totreat.push_back(LCC());
  totreat.back()=lcc;
  while(!totreat.empty())
  {
    LCC current;
    current->swap(totreat.front());
    totreat.pop_front();

    std::size_t replaced=std::numeric_limits<std::size_t>::max();
    for(auto it=current.darts().begin(), itend=current.darts().end();
        it!=itend; ++it)
    {
      Signature signature;
      Dart_handle
          dh2=ssignature_of_surface_for_dart(current, it, marktopreserve, signature); //, true);
      auto res=m_ssignatures.find(signature);
      if(res!=m_ssignatures.end())
      {
        totreat.push_back(LCC());
        std::unordered_map<Dart_handle, Dart_handle> origin_to_copy;
        totreat.back().copy(current, &origin_to_copy, nullptr);
        replace_one_surface_from_dart(totreat.back(), origin_to_copy[dh2],
                                      m_spatterns[res->second.second],
                                      res->second.first);
        replaced=res->second.second;
      }
    }
    if(replaced==std::numeric_limits<std::size_t>::max())
    {
      reslccs.push_back(LCC());
      current->swap(reslccs.back());
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void generate_all_volume_replacement(LCC& lcc,
                                     size_type marktopreserve,
                                     std::list<LCC>& reslccs)
{
  std::list<LCC> totreat; // list to avoid copy of LCC
  totreat.push_back(LCC());
  totreat.back()=lcc;
  while(!totreat.empty())
  {
    LCC current;
    current->swap(totreat.front());
    totreat.pop_front();

    std::size_t replaced=std::numeric_limits<std::size_t>::max();
    for(auto it=current.darts().begin(), itend=current.darts().end();
        it!=itend; ++it)
    {
      Signature signature;
      Dart_handle
          dh2=vsignature_of_volume_for_dart(current, it, marktopreserve, signature); //, true);
      auto res=m_vsignatures.find(signature);
      if(res!=m_vsignatures.end())
      {
        totreat.push_back(LCC());
        std::unordered_map<Dart_handle, Dart_handle> origin_to_copy;
        totreat.back().copy(current, &origin_to_copy, nullptr);
        replace_one_volume_from_dart(totreat.back(), origin_to_copy[dh2],
                                     m_vpatterns[res->second.second],
                                     res->second.first);
        replaced=res->second.second;
      }
    }
    if(replaced==std::numeric_limits<std::size_t>::max())
    {
      reslccs.push_back(LCC());
      current->swap(reslccs.back());
    }
  }
}

public:
  Pattern_set<1> m_fpatterns;
  Signature_mapping m_fsignatures;

  Pattern_set<2> m_spatterns;
  Signature_mapping m_ssignatures;

  Pattern_set<3> m_vpatterns;
  Signature_mapping m_vsignatures;
};
////////////////////////////////////////////////////////////////////////////////
#endif // CMAP_QUERY_REPLACE_H
