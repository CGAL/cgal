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
#ifndef CMAP_SIGNATURE_H
#define CMAP_SIGNATURE_H

#include <cassert>
#include <unordered_map>
#include <queue>
#include <functional>
#include <iostream>
#include <boost/container_hash/hash.hpp>

// We have 3 type of words:
//      fword for faces; sword for surfaces; vword for volumes
// 3 types of patterns:
//      fpattern: a connected set of faces
//      spattern: a closed surface (a set of connected faces without boundary)
//      vpattern: a connected sef of volumes
// and 3 types of signatures:
//      fsignature; ssignature; vsignature
///////////////////////////////////////////////////////////////////////////////
using MyInt=std::uint16_t;
using Signature=std::vector<MyInt>;
////////////////////////////////////////////////////////////////////////////////
namespace  std {
template<>
class hash<Signature>
{
public:
  size_t operator() (const Signature& s) const
  {
    std::size_t seed=0;
    for(auto n: s)
    { boost::hash_combine(seed, n); }
    return seed;
  }
};
}
////////////////////////////////////////////////////////////////////////////////
void print_signature(const Signature& s)
{
  bool first=true;
  std::cout<<"[";
  for(auto n: s)
  { if(!first) { std::cout<<" "; } else { first=false; } std::cout<<(int)n; }
  std::cout<<"]  "<<std::hash<Signature>()(s)<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
//// Signature for faces //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// Compute the face word of the given map starting from a dart.
/// Takes one function next as parameter, allowing to change
/// the object considered (face or border of the fpattern)
/// If signature is non empty, compare the current word with signature, and
/// stop as soon as the word becomes bigger than signature.
/// @return true iff the computed word is the new minimal one
template<class CMAP>
bool compute_fword_from_dart(CMAP& cmap,
                             typename CMAP::Dart_handle dh,
                             typename CMAP::size_type marktopreserve,
                             Signature& word,
                             const Signature& signature,
                             std::function<typename CMAP::Dart_handle
                             (typename CMAP::Dart_handle)> next,
                             bool trace=false)
{
  word.clear();
  if(marktopreserve!=CMAP::INVALID_MARK && !cmap.is_marked(dh, marktopreserve))
  { return false; }

  if(!signature.empty()) { word.reserve(signature.size()); }
  MyInt nb=0;
  typename CMAP::Dart_handle cur=dh;
  bool same_prefix=true, bigger=false;
  do
  {
    nb=0;
    do
    {
      ++nb;
      cur=next(cur);
    }
    while(cur!=dh &&
          (marktopreserve==CMAP::INVALID_MARK ||
           !cmap.is_marked(cur, marktopreserve)));
    word.push_back(nb);
    if(same_prefix && !signature.empty())
    {
      if(word.back()!=signature[word.size()-1])
      {
        same_prefix=false;
        if(word.back()>signature[word.size()-1])
        {
          bigger=true;
          assert(word>signature);
        }
      }
    }
  }
  while(!bigger && cur!=dh);

  if(trace)
  {
    bool first=true;
    std::cout<<"[";
    for(auto n: word)
    { if(!first) { std::cout<<" "; } else { first=false; } std::cout<<(int)n; }
    std::cout<<"]  "<<std::endl;
  }

  if(signature.empty() || (!bigger && !same_prefix))
  {
    assert(signature.empty() || word<signature);
    return true; // word<signature
  }
  assert(!signature.empty() && word>=signature);
  return false; // word>=signature

}
///////////////////////////////////////////////////////////////////////////////
/// Compute the face signature of the given pattern.
/// @pre cmap is a fpattern, i.e. a connected set of faces.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle fsignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.template is_free<2>(it))
    {
      if(compute_fword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                 [&cmap](typename CMAP::Dart_handle dh)
                                 -> typename CMAP::Dart_handle
                                 { typename CMAP::Dart_handle other=
                                 cmap.template beta<1>(dh);
                                 while(!cmap.template is_free<2>(other))
                                 { other=cmap.template beta<2,1>(other); }
                                 return other;
                                 },

                                 trace))
      {
        res=it;
        std::swap(current_word, signature);
        if(signature.size()==1)
        { return res; } // No need to test all starting darts if we have only one value
      }
    }
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the face signature of the given face.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle fsignature_of_face_for_dart
(CMAP& cmap, typename CMAP::Dart_handle dh,
 typename CMAP::size_type marktopreserve, Signature& signature, bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  if(compute_fword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_handle dh)
                             -> typename CMAP::Dart_handle
                             { return cmap.template beta<1>(dh); },

                             trace))
  {
    res=dh;
    std::swap(current_word, signature);
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the face signature of the given face.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle fsignature_of_face(CMAP& cmap,
                                              typename CMAP::Dart_handle dh,
                                              typename CMAP::size_type marktopreserve,
                                              Signature& signature,
                                              bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  typename CMAP::Dart_handle cur=dh;
  do
  {
    if(compute_fword_from_dart(cmap, cur, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                               { return cmap.template beta<1>(dh); },

                               trace))
    {
      res=cur;
      std::swap(current_word, signature);
      if(signature.size()==1)
      { return res; } // No need to test all starting darts if we have only one value
    }
    cur=cmap.template beta<1>(cur);
  }
  while(cur!=dh);
  return res;
}
///////////////////////////////////////////////////////////////////////////////
//// Signature for volumes/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// Compute the volume word of the given map starting from a dart.
/// Takes two functions next and opposite as parameter, allowing to change
/// the object considered (volume, surface...)
/// If signature is non empty, compare the current word with signature, and
/// stop as soon as the word becomes bigger than signature.
/// @return true iff the computed word is the new minimal one
template<class CMAP>
bool compute_vword_from_dart(CMAP& cmap,
                             typename CMAP::Dart_handle dh,
                             typename CMAP::size_type marktopreserve,
                             Signature& word,
                             const Signature& signature,
                             std::function<typename CMAP::Dart_handle
                             (typename CMAP::Dart_handle)> next,
                             std::function<typename CMAP::Dart_handle
                             (typename CMAP::Dart_handle)> opposite,
                             bool trace=false)
{
  word.clear();
  if(marktopreserve!=CMAP::INVALID_MARK && !cmap.is_marked(dh, marktopreserve))
  { return false; }

  if(!signature.empty()) { word.reserve(signature.size()); }
  typename CMAP::size_type amark=cmap.get_new_mark();
  std::unordered_map<typename CMAP::Dart_handle, MyInt> indices;
  std::queue<typename CMAP::Dart_handle> to_treat;
  std::vector<typename CMAP::Dart_handle> to_unmark;
  typename CMAP::Dart_handle cur, other;
  bool same_prefix=true, bigger=false;
  MyInt nb=1;

  to_treat.push(dh);
  cmap.mark(dh, amark);
  to_unmark.push_back(dh);
  indices[dh]=nb++;
  while(!bigger && !to_treat.empty())
  {
    cur=to_treat.front();
    to_treat.pop();

    if(marktopreserve==CMAP::INVALID_MARK || cmap.is_marked(cur, marktopreserve))
    {
      word.push_back(0);
      if(same_prefix && !signature.empty())
      {
        if(word.back()!=signature[word.size()-1])
        {
          same_prefix=false;
          if(word.back()>signature[word.size()-1])
          {
            bigger=true;
            assert(word>signature);
          }
        }
      }
    }

    // Process next then opposite
    for(auto f: {next, opposite})
    {
      other=f(cur);
      assert(other!=cmap.null_handle);
      if(!cmap.is_marked(other, amark))
      {
        to_treat.push(other);
        cmap.mark(other, amark);
        to_unmark.push_back(other);
        assert(nb!=std::numeric_limits<MyInt>::max());
        indices[other]=nb++;
      }
      assert(indices.count(other)==1);
      word.push_back(indices[other]);
      if(same_prefix && !signature.empty())
      {
        if(word.back()!=signature[word.size()-1])
        {
          same_prefix=false;
          if(word.back()>signature[word.size()-1])
          {
            bigger=true;
            assert(word>signature);
          }
        }
      }
    }
  }

  if(trace)
  {
    bool first=true;
    std::cout<<"[";
    for(auto n: word)
    { if(!first) { std::cout<<" "; } else { first=false; } std::cout<<(int)n; }
    std::cout<<"]  "<<std::endl;
  }

  for(auto dhtou: to_unmark)
  { cmap.unmark(dhtou, amark); }
  cmap.free_mark(amark);
  if(signature.empty() || (!bigger && !same_prefix))
  {
    assert(signature.empty() || word<signature);
    return true; // word<signature
  }
  assert(!signature.empty() && word>=signature);
  return false; // word>=signature
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the volume signature of the given pattern.
/// @pre cmap is a connected set of volumes.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle vsignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.template is_free<3>(it))
    {
      if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                 [&cmap](typename CMAP::Dart_handle dh)
                                 -> typename CMAP::Dart_handle
                                 { return cmap.template beta<1>(dh); },

                                 [&cmap](typename CMAP::Dart_handle dh)
                                 -> typename CMAP::Dart_handle
                                 { typename CMAP::Dart_handle other=
                                 cmap.template beta<2>(dh);
                                 while(!cmap.template is_free<3>(other))
                                 { other=cmap.template beta<3,2>(other); }
                                 return other;
                                 },

                                 trace))
      {
        res=it;
        std::swap(current_word, signature);
      }
    }
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the vsignature of one volume but only for the given dart.
/// Function used only in order to validate all the patterns. Use
/// vsignature_of_volume instead.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle vsignature_of_volume_for_dart
(CMAP& cmap, typename CMAP::Dart_handle dh,
 typename CMAP::size_type marktopreserve, Signature& signature, bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  if(compute_vword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_handle dh)
                             -> typename CMAP::Dart_handle
                             { return cmap.template beta<1>(dh); },

                             [&cmap](typename CMAP::Dart_handle dh)
                             -> typename CMAP::Dart_handle
                             { return cmap.template beta<2>(dh); },

                             trace))
  {
    res=dh;
    std::swap(current_word, signature);
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the vsignature of one volume.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle vsignature_of_volume(CMAP& cmap,
                                                typename CMAP::Dart_handle dh,
                                                typename CMAP::size_type marktopreserve,
                                                Signature& signature,
                                                bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
        itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                                { return cmap.template beta<1>(dh); },

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                                { return cmap.template beta<2>(dh); },

                               trace))
    {
      res=it;
      std::swap(current_word, signature);
    }
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the number of automorphisms of the vpattern of the given CMap,
/// knowing its vsignature.
/// @pre cmap is a vpattern, i.e. a connected set of volumes.
template<class CMAP>
std::size_t number_of_automorphisms_of_vpattern(CMAP& cmap,
                                                typename CMAP::size_type marktopreserve,
                                                const Signature& signature,
                                                bool trace=false)
{
  Signature current_word;
  std::size_t nb=0;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.template is_free<3>(it))
    {
      if(!compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                  [&cmap](typename CMAP::Dart_handle dh)
                                  -> typename CMAP::Dart_handle
                                  { return cmap.template beta<1>(dh); },

                                  [&cmap](typename CMAP::Dart_handle dh)
                                  -> typename CMAP::Dart_handle
                                  { typename CMAP::Dart_handle other=
                                  cmap.template beta<2>(dh);
                                  while(!cmap.template is_free<3>(other))
                                  { other=cmap.template beta<3,2>(other); }
                                  return other;
                                  },

                                  trace))
      {
        // std::cout<<&*it<<"  "; print_signature(current_word);
        if(current_word==signature) { ++nb; }
      }
    }
  }
  return nb;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the number of automorphisms of the volume, knowing its signature.
template<class CMAP>
std::size_t number_of_automorphisms_of_volume
(CMAP& cmap, typename CMAP::Dart_handle dh,
 const Signature& signature, typename CMAP::size_type marktopreserve,
 bool trace=false)
{
  Signature current_word;
  std::size_t nb=0;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
      itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(!compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                [&cmap](typename CMAP::Dart_handle dh)
                                -> typename CMAP::Dart_handle
                                  { return cmap.template beta<1>(dh); },

                                [&cmap](typename CMAP::Dart_handle dh)
                                -> typename CMAP::Dart_handle
                                 { return cmap.template beta<2>(dh); },

                                trace))
    {
      if(current_word==signature) { ++nb; }
    }
  }
  return nb;
}
///////////////////////////////////////////////////////////////////////////////
//// Signature for surfaces////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// Note: compute_vword_from_dart is reused to compute sword, only changing
/// the next and opposite methods.
///////////////////////////////////////////////////////////////////////////////
/// Compute the ssignature of the given pattern.
/// @pre cmap is a spattern, i.e. a closed connected set of faces.
/// @return the initial dart of the signature
/// @note Be careful: contrary to f and vsignature, a ssignature does not exist
///       without some darts marked as face borders. Indeed, without marked darts,
///       it is not possible to know which edges are border of faces and which
///       ones are not. Note that this mark is only used for the pattern, not
///       for the target.
template<class CMAP>
typename CMAP::Dart_handle ssignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type faceborder,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.is_marked(it, faceborder) &&
       compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap, faceborder](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                               { typename CMAP::Dart_handle other=
                               cmap.template beta<1>(dh);
                               while(!cmap.is_marked(other, faceborder))
                               { other=cmap.template beta<2,1>(other); }
                               return other;
                               },

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                                { return cmap.template beta<2>(dh); },

                               trace))
    {
      res=it;
      std::swap(current_word, signature);
      if(signature.size()==1)
      { return res; } // No need to test all starting darts if we have only one value
    }
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the ssignature of one surface.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle ssignature_of_surface_for_dart
(CMAP& cmap, typename CMAP::Dart_handle dh,
 typename CMAP::size_type marktopreserve, Signature& signature, bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  if(compute_vword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_handle dh)
                             -> typename CMAP::Dart_handle
                             { return cmap.template beta<1>(dh); },

                             [&cmap](typename CMAP::Dart_handle dh)
                             -> typename CMAP::Dart_handle
                             { return cmap.template beta<2>(dh); },

                             trace))
  {
    res=dh;
    std::swap(current_word, signature);
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the ssignature of one surface.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_handle ssignature_of_surface(CMAP& cmap,
                                                 typename CMAP::Dart_handle dh,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_handle res=nullptr;
  Signature current_word;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
        itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                                { return cmap.template beta<1>(dh); },

                               [&cmap](typename CMAP::Dart_handle dh)
                               -> typename CMAP::Dart_handle
                                { return cmap.template beta<2>(dh); },

                               trace))
    {
      res=it;
      std::swap(current_word, signature);
    }
  }
  return res;
}
///////////////////////////////////////////////////////////////////////////////
#endif // CMAP_SIGNATURE_H
