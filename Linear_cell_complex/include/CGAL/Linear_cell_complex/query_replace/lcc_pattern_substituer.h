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
#ifndef LCC_PATTERN_SUBSTITUER_H
#define LCC_PATTERN_SUBSTITUER_H

#include <filesystem>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <CGAL/Combinatorial_map/internal/cmap_isomorphisms.h>
#include <CGAL/Linear_cell_complex/query_replace/cmap_signature.h>
#include <CGAL/Linear_cell_complex/query_replace/lcc_pattern.h>
#include <CGAL/Linear_cell_complex/IO/lcc_read_write_depending_extension.h>

namespace CGAL::internal
{
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern_substituer
{
public:
  using Dart_descriptor=typename LCC::Dart_descriptor;
  using size_type=typename LCC::size_type;

  using Signature_mapping=
      std::unordered_map<Signature, std::pair<Dart_descriptor, std::size_t>>;

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
    Dart_descriptor dh;
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

  void load_fpattern(const std::string& filename,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    load_pattern<1>(filename, m_fpatterns);
    Signature signature;
    Dart_descriptor dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    m_fsignatures.clear();
    auto& pattern = m_fpatterns[0];
    if(init_topreserve!=nullptr) // true iff the std::function is not empty
    {
      mark_to_preserve=pattern.reserve_mark_to_preserve();
      init_topreserve(pattern.lcc(), mark_to_preserve);
    }
    dh=fsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
    pattern.compute_barycentric_coord();
    m_fsignatures[signature]=std::make_pair(dh, 0);
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
    Dart_descriptor dh;
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
    Dart_descriptor dh;
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
        CGAL_assertion(pattern.lcc().is_marked(dh, pattern.m_mark_faceborder));
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
    Dart_descriptor dh;
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

  void load_vpattern(const std::string& filename,
                      std::function<void(LCC&, size_type)> init_topreserve=nullptr)
  {
    load_pattern<3>(filename, m_vpatterns);
    Signature signature;
    Dart_descriptor dh;
    size_type mark_to_preserve=LCC::INVALID_MARK;
    m_vsignatures.clear();
    auto& pattern = m_vpatterns[0];
    if(init_topreserve!=nullptr) // true iff the std::function is not empty
    {
      mark_to_preserve=pattern.reserve_mark_to_preserve();
      init_topreserve(pattern.lcc(), mark_to_preserve);
    }
    dh=vsignature_of_pattern(pattern.lcc(), mark_to_preserve, signature, false);
    pattern.compute_barycentric_coord();
    m_vsignatures[signature]=std::make_pair(dh, 0);
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
    Dart_descriptor dh;
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
         IO::is_an_lcc_known_extension(dir_entry.path().string()))
      { ++nb; }
    }

    patterns.resize(nb);
    nb=0;
    // std::cout<<"##############################"<<std::endl;
    for(auto const& dir_entry: std::filesystem::directory_iterator{dir})
    {
      if(dir_entry.is_regular_file() &&
         IO::is_an_lcc_known_extension(dir_entry.path().string()))
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
      || !IO::is_an_lcc_known_extension(file.string()))
    { return {false, 0}; }

    std::size_t id = patterns.size();
    patterns.push_back(Pattern<LCC,type>());

    IO::read_depending_extension(file.string(), patterns[id].lcc());

    return {true, id};
  }

  template<unsigned int type>
  void load_pattern(const std::string& filename,
                    Pattern_set<type>& patterns)
  {
    patterns.clear();

    patterns.resize(1);
    IO::read_depending_extension(filename, patterns[0].lcc());
  }

public:

  ////////////////////////////////////////////////////////////////////////////////
  std::size_t query_replace_one_volume_from_dart(LCC& lcc,
                                                 Dart_descriptor dh,
                                                 size_type marktopreserve)
  {
    std::size_t replaced=std::numeric_limits<std::size_t>::max();
    Signature word_signature;
    Dart_descriptor
      dh2=vsignature_of_volume_for_dart(lcc, dh, marktopreserve, word_signature); //, true);

    if (dh2==lcc.null_descriptor) { return replaced; }

    Signature signature;
    vsignature_of_volume(lcc, dh, marktopreserve, signature);
    if(signature!=word_signature) { return replaced; }

    return replace_one_volume_from_signature(lcc, dh, signature, dh2);
  }
  ////////////////////////////////////////////////////////////////////////////////
  bool replace_one_volume_from_signature(LCC& lcc,
                                         Dart_descriptor dh1,
                                         Signature& signature,
                                         Dart_descriptor dh2)
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
                                       Dart_descriptor dh,
                                       size_type marktopreserve=LCC::INVALID_MARK)
  {
    Signature signature;
    Dart_descriptor
      dh2=vsignature_of_volume(lcc, dh, marktopreserve, signature); //, true);
    return replace_one_volume_from_signature(lcc, dh, signature, dh2);
  }
  ////////////////////////////////////////////////////////////////////////////////
  /// Query volume(dh) but without using signatures. If one pattern matches,
  /// replace volume(dh).
  /// @return the index of the replaced pattern, max(std::size_t) if no match.
  std::size_t query_replace_one_volume_without_signature(LCC& lcc,
                                                         Dart_descriptor dh,
                                                         size_type marktopreserve=LCC::INVALID_MARK)
  {
    Dart_descriptor res=lcc.null_descriptor, sd=dh;

    if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
    {
      auto it=lcc.template darts_of_cell<3>(dh).begin(),
        itend=lcc.template darts_of_cell<3>(dh).end();
      while(it!=itend && !lcc.is_marked(it, marktopreserve))
      { ++it; }
      if(it!=itend) { sd=it; }
    }

    std::size_t i=0;
    while(res==lcc.null_descriptor && i<number_of_vpatterns())
    {
      res=is_volume_isomorphic_to_vpattern(lcc, sd, vpattern(i), marktopreserve,
                                           m_vpatterns[i].mark_to_preserve(),
                                           false, false, false);
      if(res==lcc.null_descriptor) { ++i; }
    }

    if(res!=lcc.null_descriptor)
    { replace_one_volume_from_dart(lcc, sd, m_vpatterns[i], res); }
    else
    { i=std::numeric_limits<std::size_t>::max(); }
    // std::cout<<"NOT found"<<std::endl;
    return i;
  }
  ////////////////////////////////////////////////////////////////////////////////
  bool replace_one_face_from_signature(LCC& lcc,
      Dart_descriptor dh1,
      Signature& signature,
      Dart_descriptor dh2)
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
/// Query face(dh) in the set of patterns, and if one pattern matches,
/// replace face(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_face(LCC& lcc,
                                   Dart_descriptor dh,
                                   size_type marktopreserve=LCC::INVALID_MARK)
{
  Signature signature;
  Dart_descriptor
      dh2=fsignature_of_face(lcc, dh, marktopreserve, signature); //, true);
   return replace_one_face_from_signature(lcc, dh, signature, dh2);
}
////////////////////////////////////////////////////////////////////////////////
/// Query face(dh) but without using signatures. If one pattern matches,
/// replace face(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_face_without_signature(LCC& lcc,
                                                     Dart_descriptor dh,
                                                     size_type marktopreserve=LCC::INVALID_MARK)
{
  Dart_descriptor res=lcc.null_descriptor, sd=dh;

  if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
  {
    auto it=lcc.template darts_of_cell<2,2>(dh).begin(),
          itend=lcc.template darts_of_cell<2,2>(dh).end();
    while(it!=itend && !lcc.is_marked(it, marktopreserve))
    { ++it; }
    if(it!=itend) { sd=it; }
  }

  std::size_t i=0;
  while(res==lcc.null_descriptor && i<number_of_fpatterns())
  {
    res=is_face_isomorphic_to_fpattern(lcc, sd, fpattern(i), marktopreserve,
                                       m_fpatterns[i].mark_to_preserve(),
                                       false, false, false);
    if(res==lcc.null_descriptor) { ++i; }
  }

  if(res!=lcc.null_descriptor)
  { replace_one_face_from_dart(lcc, sd, m_fpatterns[i], res); }
  else
  { i=std::numeric_limits<std::size_t>::max(); }
  // std::cout<<"NOT found"<<std::endl;
  return i;
}
////////////////////////////////////////////////////////////////////////////////
/// Query surface(dh) in the set of patterns, and if one pattern matches,
/// replace surface(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_surface(LCC& lcc,
                                      Dart_descriptor dh,
                                      size_type marktopreserve=LCC::INVALID_MARK)
{
  Signature signature;
  Dart_descriptor
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
    // CGAL_assertion(lcc.is_valid());
  }
  // else { std::cout<<"NOT found"<<std::endl; }
  return replaced;
}
////////////////////////////////////////////////////////////////////////////////
/// Query surface(dh) but without using signatures. If one pattern matches,
/// replace surface(dh).
/// @return the index of the replaced pattern, max(std::size_t) if no match.
std::size_t query_replace_one_surface_without_signature(LCC& lcc,
                                                        Dart_descriptor dh,
                                                        size_type marktopreserve=LCC::INVALID_MARK)
{
  Dart_descriptor res=lcc.null_descriptor, sd=dh;

  if(marktopreserve!=LCC::INVALID_MARK && !lcc.is_marked(dh, marktopreserve))
  {
    auto it=lcc.template darts_of_cell<3>(dh).begin(),
          itend=lcc.template darts_of_cell<3>(dh).end();
    while(it!=itend && !lcc.is_marked(it, marktopreserve))
    { ++it; }
    if(it!=itend) { sd=it; }
  }

  std::size_t i=0;
  while(res==lcc.null_descriptor && i<number_of_spatterns())
  {
    res=is_surface_isomorphic_to_spattern(lcc, sd, spattern(i),
                                          m_spatterns[i].m_mark_faceborder,
                                          marktopreserve,
                                          m_spatterns[i].mark_to_preserve(),
                                          false, false, false);
    if(res==lcc.null_descriptor) { ++i; }
  }

  if(res!=lcc.null_descriptor)
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
      if(all)
      { lcc.template unmark_cell<3>(it, amark); }
      else
      { lcc.negate_mark(amark); } // All darts are not marked
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
          CGAL_assertion(lcc.is_whole_map_unmarked(amark));
          lcc.free_mark(amark);
          return 1;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  CGAL_assertion(lcc.is_whole_map_unmarked(amark));
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
      if(all)
      { lcc.template unmark_cell<3>(it, amark); }
      else
      { lcc.negate_mark(amark); } // All darts are not marked
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
          CGAL_assertion(lcc.is_whole_map_unmarked(amark));
          lcc.free_mark(amark);
          return 1;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  CGAL_assertion(lcc.is_whole_map_unmarked(amark));
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
      if(all)
      { lcc.template unmark_cell<2>(it, amark); }
      else
      { lcc.negate_mark(amark); } // All darts are not marked
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
          CGAL_assertion(lcc.is_whole_map_unmarked(amark));
          lcc.free_mark(amark);
          return 1;
        }
        if(trace) { std::cout<<replaced+1<<" "; }
      }
    }
  }

  CGAL_assertion(lcc.is_whole_map_unmarked(amark));
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
        Dart_descriptor dh2=it;
        do
        {
          totreat.push_back(LCC());
          std::unordered_map<Dart_descriptor, Dart_descriptor> origin_to_copy;
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
      Dart_descriptor
          dh2=ssignature_of_surface_for_dart(current, it, marktopreserve, signature); //, true);
      auto res=m_ssignatures.find(signature);
      if(res!=m_ssignatures.end())
      {
        totreat.push_back(LCC());
        std::unordered_map<Dart_descriptor, Dart_descriptor> origin_to_copy;
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
      Dart_descriptor
          dh2=vsignature_of_volume_for_dart(current, it, marktopreserve, signature); //, true);
      auto res=m_vsignatures.find(signature);
      if(res!=m_vsignatures.end())
      {
        totreat.push_back(LCC());
        std::unordered_map<Dart_descriptor, Dart_descriptor> origin_to_copy;
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
}
#endif // LCC_PATTERN_SUBSTITUER_H
