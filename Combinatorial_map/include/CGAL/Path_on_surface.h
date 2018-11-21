// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
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

#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <CGAL/Path_generators.h>
#include <CGAL/Combinatorial_map_operations.h>

namespace CGAL {

template<typename Map_>
class Path_on_surface_with_rle;

template<typename Map_>
class Path_on_surface
{
  friend class Path_on_surface_with_rle<Map_>;

public:
  typedef Map_ Map;
  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::Dart_const_handle Dart_const_handle;

  typedef Path_on_surface<Map> Self;

  Path_on_surface(const Map& amap) : m_map(amap), m_is_closed(false)
  {}

  Path_on_surface(const Path_on_surface_with_rle<Map>& apath) :
    m_map(apath.get_map()),
    m_is_closed(apath.is_closed())
  {
    for (auto it=apath.m_path.begin(), itend=apath.m_path.end(); it!=itend; ++it)
    {
      push_back(it->first, false);
      if (it->second>0)
      { CGAL::extend_straight_positive(*this, it->second-1, false); }
      else if (it->second<0)
      { CGAL::extend_straight_negative(*this, -(it->second)-1, false); }
    }
  }

  void swap(Self& p2)
  {
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
  ///         all possible starting darts.
  bool operator==(const Self& other) const
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
  bool operator==(const char*  other) const
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

  // @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  std::size_t length() const
  { return m_path.size(); }

  // @return true iff the path is closed (update after each path modification).
  bool is_closed() const
  { return m_is_closed; }

  const Map& get_map() const
  { return m_map; }

  void clear()
  {
    m_path.clear();
    m_is_closed=false;
  }
  
  void cut(std::size_t n, bool update_isclosed=true)
  {
    if (n>=length()) return;
    m_path.resize(n);
    if (update_isclosed) { update_is_closed(); }
  }

  std::size_t next_index(std::size_t i) const
  { return (is_closed() && i==m_path.size()-1?0:i+1); }

  std::size_t prev_index(std::size_t i) const
  { return (is_closed() && i==0?m_path.size()-1:i-1); }

  Dart_const_handle get_ith_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return m_path[i];
  }
  
  Dart_const_handle operator[] (std::size_t i) const
  { return get_ith_dart(i); }

  Dart_const_handle get_prev_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==0 && !is_closed()) return NULL;
    return m_path[prev_index(i)];
  }

  Dart_const_handle get_next_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==m_path.size()-1 && !is_closed()) return NULL;
    return m_path[next_index(i)];
  }

  Dart_const_handle front() const
  {
    CGAL_assertion(!is_empty());
    return m_path.front();
  }
  
  Dart_const_handle back() const
  {
    CGAL_assertion(!is_empty());
    return m_path.back();
  }
  
  void push_back(Dart_const_handle dh, bool update_isclosed=true)
  {
    CGAL_assertion(dh!=NULL && dh!=m_map.null_dart_handle);
    /* This assert is too long...
     CGAL_assertion((is_empty() ||
           CGAL::template belong_to_same_cell<Map, 0>
           (m_map, m_map.other_extremity(back()), dh))); */

    m_path.push_back(dh);
    if (update_isclosed) { update_is_closed(); }
  }

  // @return true iff the path is valid; i.e. a sequence of edges two by
  //              two adjacent.
  bool is_valid() const
  {
    for (unsigned int i=1; i<m_path.size(); ++i)
    {
      if (!m_map.darts().owns(m_path[i]))
      { return false; }

      if (m_path[i]==NULL || m_path[i]==m_map.null_dart_handle)
      { return false; }

      Dart_const_handle pend=m_map.other_extremity(m_path[i-1]);
      if (pend==Map::null_handle) { return false; }

      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, m_path[i], pend))
      { return false; }
    }
    if (is_closed())
    {
      Dart_const_handle pend=m_map.other_extremity(m_path[m_path.size()-1]);
      if (pend==Map::null_handle) { return false; }
      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, pend, m_path[0]))
      { return false; }
    }

    return true;
  }

  // Update m_is_closed to true iff the path is closed (i.e. the second
  //   extremity of the last dart of the path is the same vertex than the one
  //   of the first dart of the path).
  void update_is_closed()
  {
    CGAL_assertion(is_valid());
    if (is_empty()) { m_is_closed=false; }
    else
    {
      Dart_const_handle pend=m_map.other_extremity(back());
      if (pend==Map::null_handle) { m_is_closed=false; }
      else
      { m_is_closed=CGAL::belong_to_same_cell<Map,0>(m_map, m_path[0], pend); }
    }
  }

  // @return true iff the path does not pass twice through a same edge
  //              or a same vertex.
  bool is_simple() const
  {
    typename Map::size_type markvertex=m_map.get_new_mark();
    typename Map::size_type markedge=m_map.get_new_mark();

    bool res=true;
    unsigned int i=0;
    for (i=0; res && i<m_path.size(); ++i)
    {
      if (m_map.is_marked(m_path[i], markvertex)) res=false;
      if (m_map.is_marked(m_path[i], markedge)) res=false;

      CGAL::mark_cell<Map, 0>(m_path[i], markvertex);
      CGAL::mark_cell<Map, 1>(m_path[i], markedge);
    }

    i=0;
    while(m_map.number_of_marked_darts(markedge)>0)
    {
      CGAL_assertion(i<m_path.size());
      CGAL::unmark_cell<Map, 0>(m_path[i], markvertex);
      CGAL::unmark_cell<Map, 1>(m_path[i], markedge);
      ++i;
    }

    m_map.free_mark(markvertex);
    m_map.free_mark(markedge);

    return res;
  }

  void reverse()
  {
    std::vector<Dart_const_handle> new_path(m_path.size());
    for (std::size_t i=0; i<m_path.size(); ++i)
    {
      new_path[m_path.size()-1-i]=m_map.template beta<2>(m_path[i]);
    }
    new_path.swap(m_path);
  }

  /// If the given path is opened, close it by doing the same path that the
  /// first one in reverse direction.
  void close()
  {
    if (!is_closed())
    {
      for (int i=m_path.size()-1; i>=0; --i)
      { m_path.push_back(m_map.template beta<2>(get_ith_dart(i))); }
      m_is_closed=true;
    }
  }

  // copy all darts starting from begin and going to the dart before end
  // from this path to new_path.
  void copy_rest_of_path(std::size_t begin, std::size_t end,
                         Self& new_path)
  {
    CGAL_assertion(end<=length());
    CGAL_assertion(begin<=end);
    while(begin!=end)
    {
      new_path.push_back(get_ith_dart(begin));
      ++begin;
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
  const Map& m_map; // The underlying map
  std::vector<Dart_const_handle> m_path; // The sequence of darts
  bool m_is_closed; // True iff the path is a cycle
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_H //
// EOF //
