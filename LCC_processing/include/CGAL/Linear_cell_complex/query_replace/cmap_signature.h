// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//

#ifndef CGAL_CMAP_SIGNATURE_H
#define CGAL_CMAP_SIGNATURE_H

#include <CGAL/license/LCC_processing.h>

#include <boost/container_hash/hash.hpp>
#include <functional>
#include <iostream>
#include <queue>
#include <ranges>
#include <unordered_map>
#include <vector>

#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Random.h>
#include <CGAL/Combinatorial_map/internal/cmap_copy.h>

// We have 3 type of words:
//      fword for faces; sword for surfaces; vword for volumes
// 3 types of patterns:
//      fpattern: a connected set of faces
//      spattern: a closed surface (a set of connected faces without boundary)
//      vpattern: a connected sef of volumes
// and 3 types of signatures:
//      fsignature; ssignature; vsignature
///////////////////////////////////////////////////////////////////////////////
namespace CGAL::internal
{
using MyInt=std::int32_t;
using Signature=std::vector<MyInt>;
const MyInt EPSILON = MyInt(0); // must be i-free
const MyInt ASTERISK = MyInt(-1); // can be anything
const MyInt NOTEPSILON = MyInt(-2); // must be i-sewn with any other dart
const MyInt NOMATCH = MyInt(-3); // can not be i-sewn with a dart of the pattern
}
////////////////////////////////////////////////////////////////////////////////
namespace  std
{
template<>
class hash<CGAL::internal::Signature>
{
public:
  size_t operator() (const CGAL::internal::Signature& s) const
  {
    std::size_t seed=0;
    for(auto n: s)
    { boost::hash_combine(seed, n); }
    return seed;
  }
};
}
////////////////////////////////////////////////////////////////////////////////
namespace CGAL::internal
{
////////////////////////////////////////////////////////////////////////////////
inline
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
                             typename CMAP::Dart_descriptor dh,
                             typename CMAP::size_type marktopreserve,
                             Signature& word,
                             const Signature& signature,
                             std::function<typename CMAP::Dart_descriptor
                             (typename CMAP::Dart_descriptor)> next,
                             bool trace=false)
{
  word.clear();
  if(marktopreserve!=CMAP::INVALID_MARK && !cmap.is_marked(dh, marktopreserve))
  { return false; }

  if(!signature.empty()) { word.reserve(signature.size()); }
  MyInt nb=0;
  typename CMAP::Dart_descriptor cur=dh;
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
          CGAL_assertion(word>signature);
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
    CGAL_assertion(signature.empty() || word<signature);
    return true; // word<signature
  }
  CGAL_assertion(!signature.empty() && word>=signature);
  return false; // word>=signature

}
///////////////////////////////////////////////////////////////////////////////
/// Compute the face signature of the given pattern.
/// @pre cmap is a fpattern, i.e. a connected set of faces.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_descriptor fsignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.template is_free<2>(it))
    {
      if(compute_fword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                 [&cmap](typename CMAP::Dart_descriptor dh)
                                 -> typename CMAP::Dart_descriptor
                                 { typename CMAP::Dart_descriptor other=
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
typename CMAP::Dart_descriptor fsignature_of_face_for_dart
(CMAP& cmap, typename CMAP::Dart_descriptor dh,
 typename CMAP::size_type marktopreserve, Signature& signature, bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  if(compute_fword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_descriptor dh)
                             -> typename CMAP::Dart_descriptor
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
typename CMAP::Dart_descriptor fsignature_of_face(CMAP& cmap,
                                              typename CMAP::Dart_descriptor dh,
                                              typename CMAP::size_type marktopreserve,
                                              Signature& signature,
                                              bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  typename CMAP::Dart_descriptor cur=dh;
  do
  {
    if(compute_fword_from_dart(cmap, cur, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
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
                             typename CMAP::Dart_descriptor dh,
                             typename CMAP::size_type marktopreserve,
                             Signature& word,
                             const Signature& signature,
                             std::function<typename CMAP::Dart_descriptor
                             (typename CMAP::Dart_descriptor)> next,
                             std::function<typename CMAP::Dart_descriptor
                             (typename CMAP::Dart_descriptor)> opposite,
                             bool trace=false)
{
  word.clear();
  if(marktopreserve!=CMAP::INVALID_MARK && !cmap.is_marked(dh, marktopreserve))
  { return false; }

  if(!signature.empty()) { word.reserve(signature.size()); }
  typename CMAP::size_type amark=cmap.get_new_mark();
  std::unordered_map<typename CMAP::Dart_descriptor, MyInt> indices;
  std::queue<typename CMAP::Dart_descriptor> to_treat;
  std::vector<typename CMAP::Dart_descriptor> to_unmark;
  typename CMAP::Dart_descriptor cur, other;
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

    if(marktopreserve==CMAP::INVALID_MARK ||
        cmap.is_marked(cur, marktopreserve))
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
            CGAL_assertion(word>signature);
          }
        }
      }
    }

    // Process next then opposite
    for(auto f: {next, opposite})
    {
      other=f(cur);
      CGAL_assertion(other!=cmap.null_descriptor);
      if(!cmap.is_marked(other, amark))
      {
        to_treat.push(other);
        cmap.mark(other, amark);
        to_unmark.push_back(other);
        CGAL_assertion(nb!=std::numeric_limits<MyInt>::max());
        indices[other]=nb++;
      }
      CGAL_assertion(indices.count(other)==1);
      word.push_back(indices[other]);
      if(same_prefix && !signature.empty())
      {
        if(word.back()!=signature[word.size()-1])
        {
          same_prefix=false;
          if(word.back()>signature[word.size()-1])
          {
            bigger=true;
            CGAL_assertion(word>signature);
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
  CGAL_assertion(cmap.is_whole_map_unmarked(amark));
  cmap.free_mark(amark);
  if(signature.empty() || (!bigger && !same_prefix))
  {
    CGAL_assertion(signature.empty() || word<signature);
    return true; // word<signature
  }
  CGAL_assertion(!signature.empty() && word>=signature);
  return false; // word>=signature
}
///////////////////////////////////////////////////////////////////////////////
/// Compute the volume signature of the given pattern.
/// @pre cmap is a connected set of volumes.
/// @return the initial dart of the signature
template<class CMAP>
typename CMAP::Dart_descriptor vsignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.template is_free<3>(it))
    {
      if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                 [&cmap](typename CMAP::Dart_descriptor dh)
                                 -> typename CMAP::Dart_descriptor
                                 { return cmap.template beta<1>(dh); },

                                 [&cmap](typename CMAP::Dart_descriptor dh)
                                 -> typename CMAP::Dart_descriptor
                                 { typename CMAP::Dart_descriptor other=
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
typename CMAP::Dart_descriptor vsignature_of_volume_for_dart
    (CMAP& cmap, typename CMAP::Dart_descriptor dh,
     typename CMAP::size_type marktopreserve, Signature& signature,
     bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  if(compute_vword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_descriptor dh)
                             -> typename CMAP::Dart_descriptor
                             { return cmap.template beta<1>(dh); },

                             [&cmap](typename CMAP::Dart_descriptor dh)
                             -> typename CMAP::Dart_descriptor
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
typename CMAP::Dart_descriptor vsignature_of_volume(CMAP& cmap,
                                                typename CMAP::Dart_descriptor dh,
                                                typename CMAP::size_type marktopreserve,
                                                Signature& signature,
                                                bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
        itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
                                { return cmap.template beta<1>(dh); },

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
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

                                  [&cmap](typename CMAP::Dart_descriptor dh)
                                  -> typename CMAP::Dart_descriptor
                                  { return cmap.template beta<1>(dh); },

                                  [&cmap](typename CMAP::Dart_descriptor dh)
                                  -> typename CMAP::Dart_descriptor
                                  { typename CMAP::Dart_descriptor other=
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
(CMAP& cmap, typename CMAP::Dart_descriptor dh,
 const Signature& signature, typename CMAP::size_type marktopreserve,
 bool trace=false)
{
  Signature current_word;
  std::size_t nb=0;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
      itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(!compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                                [&cmap](typename CMAP::Dart_descriptor dh)
                                -> typename CMAP::Dart_descriptor
                                  { return cmap.template beta<1>(dh); },

                                [&cmap](typename CMAP::Dart_descriptor dh)
                                -> typename CMAP::Dart_descriptor
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
typename CMAP::Dart_descriptor ssignature_of_pattern(CMAP& cmap,
                                                 typename CMAP::size_type faceborder,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  for(auto it=cmap.darts().begin(), itend=cmap.darts().end(); it!=itend; ++it)
  {
    if(cmap.is_marked(it, faceborder) &&
       compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap, faceborder](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
                               { typename CMAP::Dart_descriptor other=
                               cmap.template beta<1>(dh);
                               while(!cmap.is_marked(other, faceborder))
                               { other=cmap.template beta<2,1>(other); }
                               return other;
                               },

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
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
typename CMAP::Dart_descriptor ssignature_of_surface_for_dart
(CMAP& cmap, typename CMAP::Dart_descriptor dh,
 typename CMAP::size_type marktopreserve, Signature& signature, bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  if(compute_vword_from_dart(cmap, dh, marktopreserve, current_word, signature,

                             [&cmap](typename CMAP::Dart_descriptor dh)
                             -> typename CMAP::Dart_descriptor
                             { return cmap.template beta<1>(dh); },

                             [&cmap](typename CMAP::Dart_descriptor dh)
                             -> typename CMAP::Dart_descriptor
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
typename CMAP::Dart_descriptor ssignature_of_surface(CMAP& cmap,
                                                 typename CMAP::Dart_descriptor dh,
                                                 typename CMAP::size_type marktopreserve,
                                                 Signature& signature,
                                                 bool trace=false)
{
  signature.clear();
  typename CMAP::Dart_descriptor res=cmap.null_descriptor;
  Signature current_word;

  for(auto it=cmap.template darts_of_cell<3>(dh).begin(),
        itend=cmap.template darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(compute_vword_from_dart(cmap, it, marktopreserve, current_word, signature,

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
                                { return cmap.template beta<1>(dh); },

                               [&cmap](typename CMAP::Dart_descriptor dh)
                               -> typename CMAP::Dart_descriptor
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
inline
bool are_signatures_equal(const Signature& signature1, const Signature& signature2)
{
  if(signature1.size() != signature2.size()) return false;
  for(unsigned int i = 0; i < signature1.size(); ++i)
  {
    if(signature1[i] != signature2[i]) return false;
  }
  return true;
}
///////////////////////////////////////////////////////////////////////////////
/// Gets the indice of beta_i of a dart from a word of a full combinatorial map
/// @param word is a word of a combinatorial map
/// @param indice is the indice of a dart. Must be between 1 and
///        (word.size()/dimension) to be valid
/// @param i is the beta we want. Must be between 1 and dimension.
/// @param dimension is the dimension of the combinatorial map
/// @return the indice of the beta_i of the dart
template<unsigned int i, unsigned int dimension>
MyInt beta_from_word(const Signature& word, const MyInt& indice)
{ return word.at((indice-1) * dimension + i - 1); }
///////////////////////////////////////////////////////////////////////////////
/// Gets the indice of beta_i of a dart from a word of a full combinatorial map
/// @param word is a word of a combinatorial map
/// @param indice is the indice of a dart. Must be between 1 and
///        (word.size()/dimension) to be valid
/// @param i is the beta we want. Must be between 1 and dimension.
/// @param dimension is the dimension of the combinatorial map
/// @return the indice of the beta_i of the dart
template<unsigned int dimension>
MyInt beta_from_word(const Signature& word, const MyInt& indice, int i)
{ return word.at((indice-1) * dimension + i - 1); }
#ifdef REFACTOR_VINCENT
///////////////////////////////////////////////////////////////////////////////
/// Checks whether a dart is compatible with an expected label

template<typename CMAP>
struct MarkConstraint {

  enum Policy {
    EQUAL,
    SUBSET,
    CONTAIN
  } ;

  typename CMAP::size_type mark ;
  Policy policy ;
  bool expected ;

  MarkConstraint(
      typename CMAP::size_type mark,
      bool expected = false,
      Policy policy = EQUAL
      ) : mark(mark), policy(policy), expected(expected) {}

  bool check(bool mark) const {
    switch(policy) {
      case EQUAL:
        return mark == expected ;
      case SUBSET:
        return mark <= expected ;
      case CONTAIN:
        return mark >= expected ;
    }
    return false ;
  }
} ;

template<
  class CMAP,
  class DartMapping,
  class IndexMapping
  >
bool check_dart_compatibility_and_mark(
    const CMAP& cmap,
    DartMapping& dart_mapping,
    IndexMapping& index_mapping,
    typename CMAP::Dart_descriptor dh,
    MyInt expected,
    std::ranges::input_range auto&& mark_constraints
    )
{
  using index_it = typename DartMapping::iterator ;

  //any dart is compatible with asterisk, forbidden or not
  if(expected == ASTERISK) return true ;

  //at this point, expected can be either :
  // a/ epsilon
  // b/ notepsilon
  // c/ nomatch
  // d/ an index
  //
  //the dart can be :
  // 1/ epsilon (therefore with index 0)
  // 2/ assigned with an index != 0
  // 3/ not assigned with an index
  // 4/ 2 or 3 and marked by a constraint
  //
  //the desired behaviour is then :
  // * a1 epsilon/epsilon => OK
  // * a2 epsilon/indexed => KO
  // * a3 epsilon/not indexed => KO
  // * a4 epsilon/marked => KO
  // * b1 notepsilon/epsilon => KO
  // * b2 notepsilon/indexed => OK
  // * b3 notepsilon/not indexed => OK
  // * b4 notepsilon/marked => OK
  // * c1 nomatch/epsilon => OK
  // * c2 nomatch/indexed => KO unless the index is notmatch
  // * c3 nomatch/not indexed => OK, mark the dart to prevent further index attempt
  // * c4 nomatch/marked => mark irrelevant, do as c3
  // * d1 index/epsilon => KO
  // * d2 index/indexed => KO unless the indices match
  // * d3 index/not_indexed => OK unless the index is already assigned to another dart
  // * d4 index/marked => fail on mark incompatibility, then do as d2 or d3

  //cases a1 a2 a3 a4
  if(expected==EPSILON) return dh == cmap.null_dart_descriptor ;
  //case b1 b2 b3 b4
  if(expected==NOTEPSILON) return dh != cmap.null_dart_descriptor ;
  //case c1
  if(expected==NOMATCH && dh==cmap.null_dart_descriptor) return true;

  //case d4
  if(expected != NOMATCH) {
    for(const MarkConstraint<CMAP>& constraint : mark_constraints) {
      //check whether the dart is marked
      bool marked = constraint.mark != cmap.INVALID_MARK && cmap.is_marked(dh, constraint.mark) ;

      //check the compatibility with the expectation
      if(!constraint.check(marked)) return false ;
    }
  }

  //check whether this dart already has an index
  index_it dh_it = dart_mapping.find(dh) ;
  if(dh_it != dart_mapping.end()) {
    //cases c2 d1 d2
    //when expected is nomatch or dh is epsilon this works as well
    return dh_it->second == expected ;
  } else {
    //cases c3 d3
    //ensure that the expected label was not previously assigned
    if(expected == NOMATCH || index_mapping[expected] == cmap.null_descriptor) {
      //the dart receives an index
      dart_mapping[dh] = expected ;
      if(expected != NOMATCH) {
        //if the index is an actual index assign the dart to the index
        index_mapping[expected] = dh ;
      }
      return true ;
    } else {
      return false ;
    }
  }
}

template<
  class CMAP,
  class DartMapping,
  class IndexMapping
  >
bool check_dart_compatibility_and_mark(
    const CMAP& cmap,
    DartMapping& dart_mapping,
    IndexMapping& index_mapping,
    typename CMAP::Dart_descriptor dh,
    MyInt expected
    ) {
  return check_dart_compatibility_and_mark(
      cmap, dart_mapping, index_mapping,
      dh, expected,
      std::views::empty<MarkConstraint<CMAP>>()
      ) ;
}

template<
  class CMAP,
  class DartMapping,
  class IndexMapping
  >
bool check_dart_compatibility_and_mark(
    const CMAP& cmap,
    DartMapping& dart_mapping,
    IndexMapping& index_mapping,
    typename CMAP::Dart_descriptor dh,
    MyInt expected,
    typename CMAP::size_type to_ignore
    ) {
  return check_dart_compatibility_and_mark(
      cmap, dart_mapping, index_mapping,
      dh, expected,
      std::views::single(MarkConstraint<CMAP>(to_ignore))
      ) ;
}
#endif
///////////////////////////////////////////////////////////////////////////////
/// Checks whether there is a match between a pattern and a part of an input
/// map by trying to recreate a partial labelling of the map, starting with a
/// dart and following a word of the pattern.
template<unsigned int dimension, class CMAP>
bool check_word_from_dart(CMAP& cmap,
                          const Signature& word,
                          typename CMAP::Dart_descriptor dh,
                          std::unordered_map<typename CMAP::Dart_descriptor, MyInt>& indices,
                          typename CMAP::size_type to_ignore = CMAP::INVALID_MARK,
                          const MyInt& starting_indice = 1
                          )
{
#ifndef REFACTOR_VINCENT
  indices.clear();
  std::queue<typename CMAP::Dart_descriptor> to_treat;
  to_treat.push(dh);
  indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>(dh, starting_indice));
  indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>(cmap.null_dart_descriptor, EPSILON));
  std::unordered_set<MyInt> used_indices;
  used_indices.insert(starting_indice);
  while(!to_treat.empty())
  {
    typename CMAP::Dart_descriptor dk = to_treat.front();
    to_treat.pop();
    MyInt k = indices.at(dk);
    for(unsigned int i = 1; i <= dimension; i++)
    {
      MyInt expected_label=beta_from_word<dimension>(word, k, i);
      if(expected_label!=ASTERISK)
      {
        typename CMAP::Dart_descriptor dn=cmap.beta(dk, i);

        if(expected_label==EPSILON && dn!=cmap.null_dart_descriptor)
        { return false; }
        if(expected_label==NOTEPSILON && dn==cmap.null_dart_descriptor)
        { return false; }
        if(expected_label==NOMATCH && dn!=cmap.null_dart_descriptor)
        {
          if(indices.find(dn)!=indices.end())
          { return false; }
        }
          if(to_ignore != cmap.INVALID_MARK && dn != cmap.null_dart_descriptor &&
              cmap.is_marked(dn, to_ignore))
          { return false; }

        if(expected_label!=EPSILON && expected_label!=NOTEPSILON &&
           expected_label!=NOMATCH)
        {
          if(indices.find(dn)==indices.end())
          {
              // check that there is no other dart with the same label already
              if(used_indices.find(expected_label)!=used_indices.end())
              { return false; }

              used_indices.insert(expected_label);
            to_treat.push(dn);

            indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>(dn, expected_label));
          }
          else
          {
            if(indices.at(dn) != expected_label)
            { return false; }
          }
        }
      }
    }
  }
  return true;
#else
  using Dart = typename CMAP::Dart_descriptor ;
  using index_it = typename std::unordered_map<Dart, MyInt>::iterator;

  //initialize indices
  indices.clear();
  //the starting dart has index epsilon
  indices[dh]=starting_indice;
  //epsilon has index epsilon
  indices[cmap.null_dart_descriptor]=EPSILON;

  //reverse index mapping to check used indices
  //can be replaced by a simple index if the word was built using a BFS
  MyInt nb_darts = word.size()/dimension ;
  std::vector<Dart> reverse_indices(nb_darts+1, cmap.null_descriptor) ;
  reverse_indices[starting_indice]=dh ;

  //traversal is only necessary here because of starting_index, when
  //starting_index is 1 a simple loop on the number of darts is enough
  //using a BFS here
  std::queue<Dart> to_treat ;
  to_treat.push(dh) ;
  while(!to_treat.empty())
  {
    //get the next dart to check which should have been indexed previously
    Dart dk=to_treat.front() ;
    to_treat.pop() ;

    //get the index of the dart
    //safety check for debug : only one dart got indexed and dn has an index
    CGAL_assertion(indices.find(dk) != indices.end() && indices[dk]> 0) ;
    MyInt k = indices[dk] ;

    //iterate over its neighboring darts in increasing dimension order
    for(unsigned int i=1; i<=dimension; ++i)
    {
      //get the neighboring dart and the label it should match
      Dart dn=cmap.beta(dk, i);
      MyInt expected=beta_from_word<dimension>(word, k, i);

      //record indices size to check whether dn receives a mark
      std::size_t nb_indexed = indices.size() ;

      //check the compatibility and mark
      if(!check_dart_compatibility_and_mark(
            cmap, indices, reverse_indices,
            dn, expected,
            to_ignore
            )
          ) return false ;

      //push dn to queue if it received an index
      if(indices.size()>nb_indexed) {
        //safety check for debug : only one dart got indexed and dn has an index
        CGAL_assertion(indices.size()==nb_indexed+1 && indices.find(dn)!=indices.end());
        if(indices[dn]!=NOMATCH) {
          //dn received a true index, push it to the queue to treat it
          to_treat.push(dn);
        }
      }
    }
  }

  //filter out the nomatch indices
  for(index_it it = indices.begin(); it != indices.end();) {
    if(it->second == NOMATCH) {
      //erase increments the iterator
      it = indices.erase(it) ;
    } else {
      ++it ;
    }
  }
  return true ;
#endif
}
///////////////////////////////////////////////////////////////////////////////
/// Checks whether there is a match between a pattern and a part of an input
/// map, including a particular dart. Checks all possible positions of the dart
/// in the pattern, by testing all the possible labels it could have.
template<unsigned int dimension, class CMAP>
void check_word_from_dart_all_positions(CMAP& cmap, const Signature& word,
                               typename CMAP::Dart_descriptor dh,
                               std::vector<typename CMAP::Dart_descriptor>& possible_matches,
                               typename CMAP::size_type to_ignore = CMAP::INVALID_MARK)
{
  for(unsigned int i=0; i<word.size()/dimension; ++i)
  {
    std::unordered_map<typename CMAP::Dart_descriptor, MyInt> indices;
    if(check_word_from_dart<dimension>(cmap, word, dh, indices, to_ignore, i+1))
    {
      for(auto pair : indices)
      {
        if(pair.second == 1)
        { possible_matches.push_back(pair.first); }
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Computes a word of a full combinatorial map. By default, i-border darts are
/// marked with NOMATCH for betai (i.e. they can not be linked with a dart of
/// the pattern), but they can be marked by epsilon, not epsilon or asterisk.
template<unsigned int dimension, class CMAP>
void word_of_n_map(CMAP& cmap, typename CMAP::Dart_descriptor starting_dart,
                   Signature& word,
                   typename CMAP::size_type epsilon=CMAP::INVALID_MARK,
                   typename CMAP::size_type notepsilon=CMAP::INVALID_MARK,
                   typename CMAP::size_type asterisk=CMAP::INVALID_MARK)
{
  word.clear();

  std::unordered_map<typename CMAP::Dart_descriptor, MyInt> indices;
  std::queue<typename CMAP::Dart_descriptor> to_treat;
  to_treat.push(starting_dart);
  indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>(starting_dart, 1));
  MyInt next_l = 2;

  indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>
                 (cmap.null_dart_descriptor, 0));

  while(!to_treat.empty())
  {
    typename CMAP::Dart_descriptor dk = to_treat.front();
    to_treat.pop();
    for(unsigned int i = 1; i <= dimension; i++)
    {
      typename CMAP::Dart_descriptor dn=cmap.beta(dk, i);

      if(indices.find(dn) == indices.end())
      {
        indices.insert(std::pair<typename CMAP::Dart_descriptor, MyInt>
                       (dn, next_l));
        next_l++;
        to_treat.push(dn);
      }

      if(dn!=cmap.null_dart_descriptor) { word.push_back(indices.at(dn)); }
      else
      {
        if(epsilon!=cmap.INVALID_MARK && cmap.is_marked(dk, epsilon))
        { word.push_back(EPSILON); }
        else if(notepsilon!=cmap.INVALID_MARK && cmap.is_marked(dk, notepsilon))
        { word.push_back(NOTEPSILON); }
        else if(asterisk!=cmap.INVALID_MARK && cmap.is_marked(dk, asterisk))
        { word.push_back(ASTERISK); }
        else { word.push_back(NOMATCH); }
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Computes the signature (smallest word in lexicographical order) of a
/// full combinatorial map.
template<unsigned int dimension, class CMAP>
typename CMAP::Dart_descriptor signature_of_n_map(CMAP& cmap, Signature& signature,
                                                  typename CMAP::size_type epsilon,
                        typename CMAP::size_type notepsilon,
                                                  typename CMAP::size_type asterisk)
{
  signature.clear();
  Signature word;
  typename CMAP::Dart_descriptor dh;
  for(auto d=cmap.darts().begin(), itend=cmap.darts().end(); d!=itend; ++d)
  {
    // TODO improve word_of_n_map complexity, taking the current minimal, and testing
    // during the computation if the word is bigger to stop the algorithm
    word_of_n_map<dimension, CMAP>(cmap, d, word, epsilon, notepsilon, asterisk);
    if(signature.empty() || word<signature)
    {
      signature=word;
      dh=d;
    }
  }
  return dh;
}
///////////////////////////////////////////////////////////////////////////////
/// builds the indices map based on signature and starting dart
/// @pre the pattern described by the word exists in lcc starting from dart dh
///      (ie check_word_from_dart returns true)
template<unsigned int dimension, class LCC>
void find_indices_from_word(LCC& lcc,
                            typename LCC::Dart_descriptor dh,
                            const Signature& word,
                            std::unordered_map<typename LCC::Dart_descriptor, MyInt>&
                                indices,
                            std::unordered_map<MyInt, typename LCC::Dart_descriptor>*
                                indices_to_darts=nullptr)
{
  indices.clear();
  indices.reserve(word.size()/dimension + 1);
  std::queue<typename LCC::Dart_descriptor> to_treat;
  to_treat.push(dh);
  indices[dh] = 1;
  indices[lcc.null_dart_descriptor] = EPSILON;
  MyInt max_label=0;
  while(!to_treat.empty())
  {
    typename LCC::Dart_descriptor dk = to_treat.front();
    to_treat.pop();
    MyInt k = indices.at(dk);
    for(unsigned int i = 1; i <= dimension; i++)
    {
      MyInt expected_label = beta_from_word<dimension>(word, k, i);
      if(    expected_label != ASTERISK
          && expected_label != NOTEPSILON
          && expected_label != NOMATCH
          && expected_label != EPSILON)
      {
        typename LCC::Dart_descriptor dn=lcc.beta(dk, i);
        if(indices.find(dn) == indices.end())
        {
          indices[dn] = expected_label;
          to_treat.push(dn);
          if(expected_label>max_label) { max_label=expected_label; }
        }
      }
    }
  }

  if(indices_to_darts!=nullptr)
  {
    indices_to_darts->clear();
    for(auto& i: indices)
    { (*indices_to_darts)[i.second]=i.first; }
  }
}
///////////////////////////////////////////////////////////////////////////////
/// @brief Copy all the query matched from lcc to lcc_copy.
template<unsigned int dimension, typename LCC>
void copy_match(LCC& lcc, const Signature& word,
                const std::unordered_map<typename LCC::Dart_descriptor, MyInt>& indices,
                LCC& lcc_copy,
                std::unordered_map
                <typename LCC::Dart_descriptor, typename LCC::Dart_descriptor>*
                    origin_to_copy=nullptr,
                std::unordered_map
                <typename LCC::Dart_descriptor, typename LCC::Dart_descriptor>*
                    copy_to_origin=nullptr)
{
  // Construct another map with reversed key and value to find descriptor from indice
  std::unordered_map<MyInt, typename LCC::Dart_descriptor> descriptor_from_indice;
  for(auto pair : indices)
  { descriptor_from_indice[pair.second] = pair.first; }

  std::vector<typename LCC::Dart_descriptor> to_copy; //contains one dart of each cell to copy

  // mark all darts of the cells to copy so we don't add more than one dart
  // of each cell in the vector
  typename LCC::size_type to_be_copied= lcc.get_new_mark();

  for(unsigned int i=1; i<=word.size()/dimension; ++i)
  {
    typename LCC::Dart_descriptor dh = descriptor_from_indice.at(i);
    if(!lcc.is_marked(dh, to_be_copied))
    {
      to_copy.push_back(dh);
      lcc.template mark_cell<dimension>(dh, to_be_copied);
    }
  }

  CGAL::CMap_copy::copy_cells<dimension>(lcc, to_copy, lcc_copy,
                                         origin_to_copy, copy_to_origin);

  /// Unmark darts
  for(unsigned int i=1; i<=word.size()/dimension; ++i)
  {
    typename LCC::Dart_descriptor dh = descriptor_from_indice.at(i);
    if(lcc.is_marked(dh, to_be_copied))
    { lcc.template unmark_cell<dimension>(dh, to_be_copied); }
  }
  CGAL_assertion(lcc.is_whole_map_unmarked(to_be_copied));
  lcc.free_mark(to_be_copied);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
struct Deterministic_volume_colors:
  public CGAL::Graphics_scene_options<LCC,
    typename LCC::Dart_const_descriptor,
    typename LCC::Dart_const_descriptor,
    typename LCC::Dart_const_descriptor,
    typename LCC::Dart_const_descriptor
  >
{
private:
  std::vector<CGAL::IO::Color> m_volume_colors; // Store the list of all volume colors
  std::unordered_map<typename LCC::Dart_const_descriptor, std::size_t>
      m_volume_id; // Give the volume id associated with each dart.

public:
  Deterministic_volume_colors(LCC& lcc)
  {
    auto treated=lcc.get_new_mark();
    for(auto it=lcc.darts().begin(), itend=lcc.darts().end(); it!=itend; ++it)
    {
      if(!lcc.is_marked(it, treated))
      {
        // Compute the signature of the volume
        Signature sig;
        vsignature_of_volume(lcc, it, LCC::INVALID_MARK, sig) ;

        // Initialize a random engine from the signature
        unsigned int seed=std::hash<Signature>()(sig) ;
        CGAL::Random alea(seed) ;

        //randomly sample a bright color
        CGAL::IO::Color res(alea.get_double(0,255),
                            alea.get_double(0,255),
                            alea.get_double(0,255));
        /* C'est vraiment très flash !!
           res.set_hsv(alea.get_double(0,360),
                    alea.get_double(50,100),
                    alea.get_double(90,100));*/
        m_volume_colors.push_back(res);

        for(auto itv=lcc.template darts_of_cell_basic<3>(it, treated).begin(),
             itvend=lcc.template darts_of_cell_basic<3>(it, treated).end();
             itv!=itvend; ++itv)
        {
          m_volume_id[itv]=m_volume_colors.size()-1;
          lcc.mark(itv, treated);
        }
      }
    }

    CGAL_assertion(lcc.is_whole_map_marked(treated));
    lcc.free_mark(treated);
  }

  // All volumes are colored.
  bool colored_volume(const LCC&, typename LCC::Dart_const_descriptor) const
  { return true; }

  // Determine the color of a volume from its signature.
  CGAL::IO::Color volume_color(const LCC& /*lcc*/,
                               typename LCC::Dart_const_descriptor d) const
  { return m_volume_colors[m_volume_id.at(d)]; }
};
///////////////////////////////////////////////////////////////////////////////
}
#endif // CGAL_CMAP_SIGNATURE_H
