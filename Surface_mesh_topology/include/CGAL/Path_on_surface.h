// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_PATH_ON_SURFACE_H
#define CGAL_PATH_ON_SURFACE_H 1

#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Random.h>
#include <boost/algorithm/searching/knuth_morris_pratt.hpp>
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace CGAL {

template<typename Map_>
class Path_on_surface_with_rle;

template<typename Map_>
class Path_on_surface
{
  friend class Path_on_surface_with_rle<Map_>;

public:
  typedef Path_on_surface<Map_>           Self;
  typedef Map_                            Map;
  typedef typename Map::Dart_const_handle Dart_const_handle;

  Path_on_surface(const Map& amap) : m_map(amap), m_is_closed(false)
  {}

  Path_on_surface(const Path_on_surface_with_rle<Map>& apath) :
    m_map(apath.get_map()),
    m_is_closed(apath.is_closed())
  {
    for (auto it=apath.m_path.begin(), itend=apath.m_path.end(); it!=itend; ++it)
    {
      push_back(it->begin, false);
      if (it->length>0)
      { extend_straight_positive(it->length, false); }
      else if (it->length<0)
      { extend_straight_negative(-(it->length), false); }
    }
    CGAL_assertion(is_valid());
  }

  void swap(Self& p2)
  {
    if (this==&p2) { return; }
    
    CGAL_assertion(&m_map==&(p2.m_map));
    m_path.swap(p2.m_path);
    std::swap(m_is_closed, p2.m_is_closed);
  }

  Self& operator=(const Self& other)
  {
    CGAL_assertion(&m_map==&(other.m_map));
    if (this!=&other)
    {
      m_path=other.m_path;
      m_is_closed=other.m_is_closed;
    }
    return *this;
  }

  /// @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  /// @return the length of the path, i.e. its number of darts.
  std::size_t length() const
  { return m_path.size(); }

  /// @return true iff the path is closed.
  ///  (m_is_closed is updated after each path modification).
  bool is_closed() const
  { return m_is_closed; }

  /// @return the combinatorial map supporting this path.
  const Map& get_map() const
  { return m_map; }

  /// clear the path.
  void clear()
  {
    m_path.clear();
    m_is_closed=false;
  }
  
  /// @return true iff the next index exists
  bool next_index_exists(std::size_t i) const
  { return is_closed() || i<(m_path.size()-1); }

  /// @return the index after index i.
  std::size_t next_index(std::size_t i) const
  { return (is_closed() && i==m_path.size()-1?0:i+1); }

  /// @return the index before index i.
  std::size_t prev_index(std::size_t i) const
  { return (is_closed() && i==0?m_path.size()-1:i-1); }

  /// @return the ith dart of the path.
  Dart_const_handle get_ith_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return m_path[i];
  }
  
  /// @return the ith dart of the path.
  Dart_const_handle operator[] (std::size_t i) const
  { return get_ith_dart(i); }

  /// @return the dart before the ith dart of the path,
  ///          NULL if such a dart does not exist.
  Dart_const_handle get_prev_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==0 && !is_closed()) return NULL;
    return m_path[prev_index(i)];
  }

  /// @return the dart after the ith dart of the path,
  ///          NULL if such a dart does not exist.
  Dart_const_handle get_next_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==m_path.size()-1 && !is_closed()) return NULL;
    return m_path[next_index(i)];
  }

  /// @return the first dart of the path.
  /// @pre !is_empty()
  Dart_const_handle front() const
  {
    CGAL_assertion(!is_empty());
    return m_path.front();
  }
  
  /// @return the last dart of the path.
  /// @pre !is_empty()
  Dart_const_handle back() const
  {
    CGAL_assertion(!is_empty());
    return m_path.back();
  }

  /// @return true iff df can be added at the end of the path.
  bool can_be_pushed(Dart_const_handle dh) const
  {
    // This assert is too long CGAL_assertion(m_map.darts().owns(dh));

    if (is_empty()) return true;

    return CGAL::template belong_to_same_cell<Map, 0>
      (m_map, m_map.other_extremity(back()), dh);
  }
  
  /// Add the given dart at the end of this path.
  /// @pre can_be_pushed(dh)
  void push_back(Dart_const_handle dh, bool update_isclosed=true)
  {
    CGAL_assertion(dh!=Map::null_handle);
    /* This assert is too long, it is tested in the is_valid method.
       CGAL_assertion(can_be_pushed(dh)); */

    m_path.push_back(dh);
    if (update_isclosed) { update_is_closed(); }
  }

  Self& operator+=(const Self& other)
  {
    m_path.reserve(m_path.size()+other.m_path.size());
    // Be careful to the special case when *this==other
    // this is the reason of the iend.
    for (std::size_t i=0, iend=other.length(); i<iend; ++i)
    { push_back(other[i], false); }
    update_is_closed();
    return *this;
  }

  Self operator+(const Self& other) const
  {
    Self res=*this;
    res+=other;
    return res;
  }

  void cut(std::size_t n, bool update_isclosed=true)
  {
    if (n>=length()) return;
    m_path.resize(n);
    if (update_isclosed) { update_isclosed; }
  }

  /// copy all darts starting from begin and going to the dart before end
  /// from this path to new_path.
  void copy_rest_of_path(std::size_t begin, std::size_t end,
                         Self& new_path)
  {
    CGAL_assertion(begin<=end);
    CGAL_assertion(end<=length());
    new_path.m_path.reserve(new_path.m_path.size()+end-begin+1);
    while(begin!=end)
    {
      new_path.push_back(get_ith_dart(begin));
      ++begin;
    }
  }

  void extend_straight_positive(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty() || nb==0)
    { return; }

    m_path.reserve(m_path.size()+nb);
    Dart_const_handle d2;
    for (std::size_t i=0; i<nb; ++i)
    {
      d2=get_map().template beta<1,2,1>(back());
      if (d2!=get_map().null_dart_handle)
      { push_back(d2, false); }
    }
    if (update_isclosed) { update_is_closed(); }
  }

  void extend_straight_negative(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty() || nb==0)
    { return; }

    m_path.reserve(m_path.size()+nb);
    Dart_const_handle d2;
    for (std::size_t i=0; i<nb; ++i)
    {
      d2=get_map().template beta<2,0,2,0,2>(back());
      if (d2!=get_map().null_dart_handle)
      { push_back(d2, false); }
    }
    if (update_isclosed) { update_is_closed(); }
  }

  void extend_straight_positive_until(Dart_const_handle dend,
                                      bool update_isclosed=true)
  {
    if (is_empty() || back()==dend)
    { return; }

    Dart_const_handle d2=get_map().template beta<1,2,1>(back());
    while(d2!=dend)
    {
      push_back(d2, false);
      d2=get_map().template beta<1,2,1>(d2);
    }
    if (update_isclosed) { update_is_closed(); }
  }

  void extend_straight_negative_until(Dart_const_handle dend,
                                      bool update_isclosed=true)
  {
    if (is_empty() || back()==dend)
    { return; }

    Dart_const_handle d2=get_map().template beta<2,0,2,0,2>(back());
    while(d2!=dend)
    {
      push_back(d2, false);
      d2=get_map().template beta<2,0,2,0,2>(d2);
    }
    if (update_isclosed) { update_is_closed(); }
  }

  void extend_positive_turn(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty()) { return; }

    if (nb==0)
    {
      if (!get_map().template is_free<2>(back()))
      { push_back(get_map().template beta<2>(back())); }
      return;
    }

    Dart_const_handle d2=get_map().template beta<1>(back());
    for (std::size_t i=1; i<nb; ++i)
    { d2=get_map().template beta<2, 1>(d2); }

    if (d2!=get_map().null_dart_handle)
    { push_back(d2, update_isclosed); }
  }

  void extend_negative_turn(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty()) { return; }

    if (nb==0)
    {
      if (!get_map().template is_free<2>(back()))
      { push_back(get_map().template beta<2>(back())); }
      return;
    }

    Dart_const_handle d2=get_map().template beta<2>(back());
    for (std::size_t i=0; i<nb; ++i)
    { d2=get_map().template beta<0, 2>(d2); }

    if (d2!=get_map().null_dart_handle)
    { push_back(d2, update_isclosed); }
  }

  /// Replace edge [i] by the path of darts along the face.
  /// Problem of complexity when used many times (like in update_path_randomly).
  void push_around_face(std::size_t i, bool update_isclosed=true)
  {
    CGAL_assertion(i<length());

    Self p2(get_map());
    std::size_t begin=i;
    Dart_const_handle dh=get_map().template beta<0>(get_ith_dart(begin));
    do
    {
      p2.push_back(get_map().template beta<2>(dh));
      dh=get_map().template beta<0>(dh);
    }
    while(dh!=get_ith_dart(begin));

    p2.m_path.reserve(p2.m_path.size()+length()-begin);
    for (std::size_t j=begin+1; j<length(); ++j)
    { p2.push_back(get_ith_dart(j), false); }

    cut(begin, false);
    m_path.reserve(m_path.size()+p2.length());
    for (std::size_t j=0; j<p2.length(); ++j)
    { push_back(p2[j], false); }

    if (update_isclosed) { update_isclosed; }
  }

  /// Push back a random dart, if the path is empty.
  bool initialize_random_starting_dart(CGAL::Random& random,
                                       bool update_isclosed=true)
  {
    if (!is_empty() || get_map().is_empty()) { return false; }

    unsigned int index=random.get_int(0, get_map().darts().capacity());
    while (!get_map().darts().is_used(index))
    {
      ++index;
      if (index==get_map().darts().capacity()) index=0;
    }
    push_back(get_map().darts().iterator_to(get_map().darts()[index]),
              update_isclosed);
    return true;
  }

  bool initialize_random_starting_dart(bool update_isclosed=true)
  {
    CGAL::Random random;
    return initialize_random_starting_dart(random, update_isclosed);
  }

  bool extend_path_randomly(CGAL::Random& random,
                            bool allow_half_turn=true,
                            bool update_isclosed=true)
  {
    if (is_empty())
    { return initialize_random_starting_dart(random, update_isclosed); }

    Dart_const_handle pend=get_map().template beta<2>(back());
    if (pend==Map::null_handle)
    {
      if (!get_map().template is_free<1>(back()))
      { // Here there is no other possibility to extend the path !
        push_back(get_map().template beta<1>(back()), update_isclosed);
        return true;
      }
      else { return false; }
    }

    unsigned int index=random.get_int  //get_int(a,b) returns an int in {a,...,b-1}
        ((allow_half_turn?0:1),
         get_map().template darts_of_cell<0>(pend).size());

    Dart_const_handle res=pend;
    for(unsigned int i=0; i<index; ++i)
    { res=get_map().template beta<2,1>(res); }

    CGAL_assertion(allow_half_turn || res!=pend);

    push_back(res, update_isclosed);
    return true;
  }

  bool extend_path_randomly(bool allow_half_turn=false,
                            bool update_isclosed=true)
  {
    CGAL::Random random;
    extend_path_randomly(random, allow_half_turn, update_isclosed);
  }

  void generate_random_path(std::size_t length, CGAL::Random& random,
                            bool update_isclosed=true)
  {
    m_path.reserve(m_path.size()+length);
    for (std::size_t i=0; i<length; ++i)
    { extend_path_randomly(random, true, false); }
    if (update_isclosed) { update_is_closed(); }
  }

  template<typename Path>
  void generate_random_path(CGAL::Random& random,
                            bool update_isclosed=true)
  { generate_random_path(random.get_int(1, 10000), random, update_isclosed); }

  template<typename Path>
  void generate_random_path(std::size_t length,
                            bool update_isclosed=true)
  {
    CGAL::Random random;
    generate_random_path(length, random, update_isclosed);
  }

  template<typename Path>
  void generate_random_path(bool update_isclosed=true)
  {
    CGAL::Random random;
    generate_random_path(random, update_isclosed);
  }

  void generate_random_closed_path(std::size_t length, CGAL::Random& random)
  {
    m_path.reserve(m_path.size()+length);
    std::size_t i=0;
    while(i<length || !is_closed())
    {
      extend_path_randomly(random, true, true);
      ++i;
    }
  }
  void generate_random_closed_path(std::size_t length)
  {
    CGAL::Random random;
    generate_random_closed_path(random, length);
  }

  void generate_random_closed_path(CGAL::Random& random)
  { generate_random_closed_path(random.get_int(1, 10000), random); }

  void generate_random_closed_path()
  {
    CGAL::Random random;
    generate_random_closed_path(random.get_int(1, 10000), random);
  }

  /// Transform the current path by pushing some dart around faces.
  /// At the end, the new path is homotopic to the original one.
  void update_path_randomly(std::size_t nb, CGAL::Random& random,
                            bool update_isclosed=true)
  {
    if (is_empty()) return;

    for (unsigned int i=0; i<nb; ++i)
    {
      push_around_face(random.get_int(0, length()), false);
    }
    if (update_isclosed) { update_isclosed; }
  }

  void update_path_randomly(CGAL::Random& random,
                            bool update_isclosed=true)
  { update_path_randomly(random.get_int(0, 10000), update_isclosed); }

  void update_path_randomly(std::size_t nb, bool update_isclosed=true)
  {
    CGAL::Random random;
    update_path_randomly(nb, random, update_isclosed);
  }

  void update_path_randomly(bool update_isclosed=true)
  {
    CGAL::Random random;
    update_path_randomly(random, update_isclosed);
  }

  /// @Return true if this path is equal to other path, identifying dart 0 of
  ///          this path with dart start in other path.
  bool are_same_paths_from(const Self& other, std::size_t start) const
  {
    CGAL_assertion(start==0 || start<length());
    CGAL_assertion(is_closed() || start==0);
    CGAL_assertion(length()==other.length() && is_closed()==other.is_closed());

    for(std::size_t i=0; i<length(); ++i)
    {
      if (get_ith_dart(i)!=other.get_ith_dart(start))
      { return false; }
      start=next_index(start);
    }
    return true;
  }

  /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts. Old quadratic version, new version
  ///         (operator==) use linear version based on Knuth, Morris, Pratt
  bool are_paths_equals(const Self& other) const
  {
    if (length()!=other.length() || is_closed()!=other.is_closed())
    { return false; }

    if (!is_closed())
    { return are_same_paths_from(other, 0); }

    for(std::size_t start=0; start<length(); ++start)
    {
      if (are_same_paths_from(other, start))
      { return true; }
    }
    return false;
  }

  /// @return true if this path is equal to other path. For closed paths,
  ///         equality is achieved whatever the first dart.
  bool operator==(const Self& other) const
  {
    if (length()!=other.length() || is_closed()!=other.is_closed())
    { return false; }

    if (!is_closed())
    { return are_same_paths_from(other, 0); }

    Self p2(*this); p2+=p2;
    // Now we search if other is a sub-motif of p2 <=> *this==other

    return boost::algorithm::knuth_morris_pratt_search(p2.m_path.begin(),
                                                       p2.m_path.end(),
                                                       other.m_path.begin(),
                                                       other.m_path.end()).first
        !=p2.m_path.end();
  }
  bool operator!=(const Self&  other) const
  { return !(operator==(other)); }

  /// @Return true if this path is equal to other path, identifying dart 0 of
  ///          this path with dart start in other path. other path is given
  ///          by index of its darts, in text format.
  bool are_same_paths_from(const char* other, std::size_t start) const
  {
    CGAL_assertion(start==0 || start<length());
    CGAL_assertion(is_closed() || start==0);

    std::string sother(other);
    std::istringstream iss(sother);
    uint64_t nb;

    for(std::size_t i=0; i<length(); ++i)
    {
      if (!iss.good())
      { return false; }
      iss>>nb;
      if (nb!=m_map.darts().index(get_ith_dart(start)))
      { return false; }
      start=next_index(start);
    }
    iss>>nb;
    if (iss.good())
    { return false; } // There are more elements in other than in this path

    return true;
  }
  /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts. other path is given by index of its
  ///         darts, in text format.
  bool operator==(const char* other) const
  {
    if (!is_closed())
    { return are_same_paths_from(other, 0); }

    for(std::size_t start=0; start<length(); ++start)
    {
      if (are_same_paths_from(other, start))
      { return true; }
    }
    return false;
  }
  bool operator!=(const char*  other) const
  { return !(operator==(other)); }


  /// @return true iff the path is valid; i.e. a sequence of edges two by
  ///              two adjacent.
  bool is_valid() const
  {
    if (is_empty()) { return !is_closed(); } // an empty past is not closed

    for (unsigned int i=1; i<m_path.size(); ++i)
    {
      /* This assert is long if (!m_map.darts().owns(m_path[i]))
      { return false; } */

      if (m_path[i]==NULL || m_path[i]==m_map.null_dart_handle)
      { return false; }

      Dart_const_handle pend=m_map.other_extremity(m_path[i-1]);
      if (pend==Map::null_handle) { return false; }

      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, m_path[i], pend))
      { return false; }
    }
    if (is_closed())
    {
      Dart_const_handle pend=m_map.other_extremity(back());
      if (pend==Map::null_handle) { return false; }
      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, pend, front()))
      { return false; }
    }
    else
    {
      Dart_const_handle pend=m_map.other_extremity(back());
      if (pend==Map::null_handle) { return true; }
      if (CGAL::template belong_to_same_cell<Map,0>(m_map, pend, front()))
      { return false; }
    }

    return true;
  }

  /// Update m_is_closed to true iff the path is closed (i.e. the second
  ///   extremity of the last dart of the path is the same vertex than the one
  ///   of the first dart of the path).
  void update_is_closed()
  {
    // CGAL_assertion(is_valid());
    if (is_empty()) { m_is_closed=false; }
    else
    {
      Dart_const_handle pend=m_map.other_extremity(back());
      if (pend==Map::null_handle) { m_is_closed=false; }
      else
      { m_is_closed=CGAL::belong_to_same_cell<Map,0>(m_map, m_path[0], pend); }
    }
  }

  /// @return true iff the path does not pass twice through a same edge
  ///              or a same vertex.
  bool is_simple() const
  {
    typename Map::size_type markvertex=m_map.get_new_mark();
    typename Map::size_type markedge=m_map.get_new_mark();

    bool res=true;
    unsigned int i=0;
    for (i=0; res && i<m_path.size(); ++i)
    {
      if (m_map.is_marked(m_path[i], markvertex)) { res=false; }
      else { CGAL::mark_cell<Map, 0>(m_map, m_path[i], markvertex); }

      if (m_map.is_marked(m_path[i], markedge)) { res=false; }
      else  { CGAL::mark_cell<Map, 1>(m_map, m_path[i], markedge); }
    }

    i=0;
    while(m_map.number_of_marked_darts(markedge)>0 ||
          m_map.number_of_marked_darts(markvertex)>0)
    {
      CGAL_assertion(i<m_path.size());
      if (m_map.is_marked(m_path[i], markvertex))
      { CGAL::unmark_cell<Map, 0>(m_map, m_path[i], markvertex); }
      if (m_map.is_marked(m_path[i], markedge))
      { CGAL::unmark_cell<Map, 1>(m_map, m_path[i], markedge); }
      ++i;
    }

    m_map.free_mark(markvertex);
    m_map.free_mark(markedge);

    return res;
  }

  /// Reverse the path (i.e. negate its orientation).
  void reverse()
  {
    std::vector<Dart_const_handle> new_path(m_path.size());
    for (std::size_t i=0; i<m_path.size()/2; ++i)
    {
      m_path[m_path.size()-1-i]=
          m_map.template beta<2>(m_path[m_path.size()-1-i]);
      m_path[i]=m_map.template beta<2>(m_path[i]);
      std::swap(m_path[i], m_path[m_path.size()-1-i]);
    }
    if (m_path.size()%2==1)
    {
      m_path[m_path.size()/2+1]=
          m_map.template beta<2>(m_path[m_path.size()/2+1]);
    }
  }

  /// If the given path is opened, close it by doing the same path that the
  /// first one in reverse direction.
  void close()
  { // TODO follow shortest path ?
    if (!is_closed())
    {
      for (int i=m_path.size()-1; i>=0; --i)
      { m_path.push_back(m_map.template beta<2>(get_ith_dart(i)), false); }
      m_is_closed=true;
    }
  }

  /// @return the turn between dart number i and dart number i+1.
  ///         (turn is position of the second edge in the cyclic ordering of
  ///          edges starting from the first edge around the second extremity
  ///          of the first dart)
  std::size_t next_positive_turn(std::size_t i) const
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(i<m_path.size());
    CGAL_assertion (is_closed() || i<length()-1);

    return m_map.positive_turn(get_ith_dart(i), get_next_dart(i));
  }

  /// Same than next_positive_turn but turning in reverse orientation around vertex.
  std::size_t next_negative_turn(std::size_t i) const
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(i<m_path.size());
    CGAL_assertion (is_closed() || i<length()-1);

    return m_map.negative_turn(get_ith_dart(i), get_next_dart(i));
  }

  std::vector<std::size_t> compute_positive_turns() const
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;

    std::size_t i;
    for (i=0; i<m_path.size()-1; ++i)
    { res.push_back(next_positive_turn(i)); }
    if (is_closed())
    { res.push_back(next_positive_turn(i)); }
    return res;
  }

  std::vector<std::size_t> compute_negative_turns() const
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;

    std::size_t i;
    for (i=0; i<m_path.size()-1; ++i)
    { res.push_back(next_negative_turn(i)); }
    if (is_closed())
    { res.push_back(next_negative_turn(i)); }
    return res;
  }

  std::vector<std::size_t> compute_turns(bool positive) const
  { return (positive?compute_positive_turns():compute_negative_turns()); }

  bool same_turns_from(const char* turns,
                       const std::vector<std::size_t>& resplus,
                       const std::vector<std::size_t>& resmoins,
                       std::size_t start) const
  {
    CGAL_assertion(start==0 || start<resplus.size());
    CGAL_assertion(resplus.size()==resmoins.size());

    std::string sturns(turns);
    std::istringstream iss(sturns);
    int64_t nb;

    for(std::size_t i=0; i<resplus.size(); ++i)
    {
      if (!iss.good())
      { return false; }
      iss>>nb;
      if ((nb>=0 && resplus[start]!=nb) ||
          (nb<0 && resmoins[start]!=-nb))
      { return false; }

      start=next_index(start);
    }
    iss>>nb;
    if (iss.good())
    { return false; } // There are more elements in turns than in res

    return true;
  }

  bool same_turns(const char* turns) const
  {
    std::vector<std::size_t> resplus=compute_positive_turns();
    std::vector<std::size_t> resmoins=compute_negative_turns();

    if (!is_closed())
    { return same_turns_from(turns, resplus, resmoins, 0); }

    for (std::size_t start=0; start<length(); ++start)
    {
      if (same_turns_from(turns, resplus, resmoins, start))
      { return true; }
    }

    return false;
  }

  void display_positive_turns() const
  {
    std::cout<<"+(";
    std::vector<std::size_t> res=compute_positive_turns();
    for (std::size_t i=0; i<res.size(); ++i)
    { std::cout<<res[i]<<(i<res.size()-1?" ":""); }
    std::cout<<")";
  }

  void display_negative_turns() const
  {
    std::cout<<"-(";
    std::vector<std::size_t> res=compute_negative_turns();
    for (std::size_t i=0; i<res.size(); ++i)
    { std::cout<<res[i]<<(i<res.size()-1?" ":""); }
    std::cout<<")";
  }

  void display_pos_and_neg_turns() const
  {
    display_positive_turns();
    std::cout<<"  ";
    display_negative_turns();
  }

  void display() const
  {
    for (std::size_t i=0; i<length(); ++i)
    {
      std::cout<<m_map.darts().index(get_ith_dart(i));
      if (i<length()-1) { std::cout<<" "; }
    }
     if (is_closed())
     { std::cout<<" c "; } //<<m_map.darts().index(get_ith_dart(0)); }
  }

  friend std::ostream& operator<<(std::ostream& os, const Self& p)
  {
    p.display();
    return os;
  }

protected:
  const Map& m_map;                      /// The underlying map
  std::vector<Dart_const_handle> m_path; /// The sequence of darts
  bool m_is_closed;                      /// True iff the path is a cycle
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_H //
// EOF //
