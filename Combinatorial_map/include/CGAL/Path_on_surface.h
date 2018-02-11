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

#include <stack>
#include <CGAL/Union_find.h>
#include <boost/unordered_map.hpp>
#include<CGAL/Random.h>

namespace CGAL {

template<typename Map>
class Path_on_surface
{
public:
  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::Dart_const_handle Dart_const_handle;

  Path_on_surface(const Map& amap) : m_map(amap)
  {}

  unsigned int length() const
  { return m_path.size(); }

  Dart_const_handle get_ith_dart(std::size_t i) const
  {
    assert(i<=m_path.size());
    return m_path[(i==m_path.size()?0:i)];
  }
  
  Dart_const_handle operator[] (std::size_t i) const
  {
    assert(i<=m_path.size());
    return m_path[(i==m_path.size()?0:i)];
  }

  void push_back(Dart_const_handle dh)
  { m_path.push_back(dh); }

  // @return true iff the path is valid; i.e. a sequence of edges two by
  //              two adjacent.
  bool is_valid() const
  {
    for (unsigned int i=1; i<m_path.size(); ++i)
    {
      Dart_const_handle pend=m_map.other_extremity(m_path[i-1]);
      if (pend==Map::null_handle) { return false; }

      if (!CGAL::template belong_to_same_cell<Map,0>(m_map, m_path[i], pend))
      { return false; }
    }
    return true;
  }

  // @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  // @return true iff the path is closed (i.e. the second extremity of the
  //              last dart of the path is the same vertex than the one of the
  //              first dart of the path.
  bool is_closed() const
  {
    if (is_empty()) { return false; } // or true by vacuity ?
    if (!is_valid()) { return false; } // Interest ??

    Dart_const_handle pend=m_map.other_extremity(m_path.back());
    if (pend==Map::null_handle) { return false; }

    return CGAL::belong_to_same_cell<Map,0>(m_map, m_path[0], pend);
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
      assert(i<m_path.size());
      CGAL::unmark_cell<Map, 0>(m_path[i], markvertex);
      CGAL::unmark_cell<Map, 1>(m_path[i], markedge);
      ++i;
    }

    m_map.free_mark(markvertex);
    m_map.free_mark(markedge);

    return res;
  }

  // Generate a random path with about percent % of edge
  void generate_random_path(unsigned int percent, CGAL::Random& random)
  {
    m_path.clear();
    unsigned int length=((m_map.number_of_darts()/2)*percent)/100;
    for (unsigned int i=0; i<length; ++i)
    { extend_randomly(random); }
  }
  void generate_random_path(unsigned int percent)
  {
    CGAL::Random random;
    generate_random_path(percent, random);
  }
  
  bool extend_randomly(CGAL::Random& random, bool allow_half_turn=false)
  {
    if (m_path.empty())
    {
      unsigned int index=random.get_int(0, m_map.darts().capacity());
      while (!m_map.darts().is_used(index))
      {
         ++index;
        if (index==m_map.darts().capacity()) index=0;
      }
      m_path.push_back(m_map.darts().iterator_to(m_map.darts()[index]));
      return true;
    }

    Dart_const_handle pend=m_map.opposite(m_path.back());
    if (pend==Map::null_handle)
    {
      if (!m_map.template is_free<1>(m_path.back()))
      {
        m_path.push_back(m_map.template beta<1>(m_path.back()));
        return true;
      }
      else { return false; }
    }

    typename Map::template Dart_of_cell_range<0>::const_iterator
        it=m_map.template darts_of_cell<0>(pend).begin();

     unsigned int index=random.get_int
         ((allow_half_turn?0:1), m_map.template darts_of_cell<0>(pend).size());
     for(unsigned int i=0; i<index; ++i)
     { ++it; }

     assert(allow_half_turn || it!=pend);
     
     m_path.push_back(it);
     return true;
  }

  /// @return the turn between dart number i and dart number i+1.
  ///         (turn is position of the second edge in the cyclic ordering of
  ///          edges starting from the first edge around the second extremity
  ///          of the first dart)
  std::size_t next_turn(std::size_t i) const
  {
    assert(is_valid());
    assert(i<m_path.size());

    Dart_const_handle d1=m_path[i];
    Dart_const_handle d2=m_path[i+1]; // Work also for the last dart.
    assert(d1!=d2);

    std::size_t res=1;
    while (m_map.template beta<1>(d1)!=d2)
    {
      ++res;
      d1=m_map.template beta<1, 2>(d1);
    }
    // std::cout<<"next_turn="<<res<<std::endl;
    return res;
  }

  /// Same than nex_turn but turning in reverse orientation around vertex.
  std::size_t next_negative_turn(std::size_t i) const
  {
    assert(is_valid());
    assert(i<m_path.size());

    Dart_const_handle d1=m_map.template beta<2>(m_path[i]);
    Dart_const_handle d2=m_map.template beta<2>(m_path[i+1]); // Work also for the last dart.
    assert(d1!=d2);

    std::size_t res=1;
    while (m_map.template beta<0>(d1)!=d2)
    {
      ++res;
      d1=m_map.template beta<0, 2>(d1);
    }
    // std::cout<<"next_negative_turn="<<res<<std::endl;
    return res;
  }


  std::size_t find_end_of_braket(std::size_t begin) const
  {
    assert(next_turn(begin)==1 || next_negative_turn(begin)==1);
    std::size_t end=begin+1;
    bool positive=(next_turn(begin)==1);

    while ((positive && next_turn(end)==2) ||
           (!positive && next_negative_turn(end)==2))
    {
      ++end;
      if (end==m_path.size()) { end=0; }
    }
    ++end;
    if (end==m_path.size()) { end=0; }
    return end;
  }

  void transform_brancket(std::size_t begin, std::size_t end,
                          std::vector<Dart_const_handle>& new_path)
  {
    Dart_const_handle d1=m_map.template beta<0>(get_ith_dart(begin));
    Dart_const_handle d2=m_map.template beta<2,0>(get_ith_dart(end));
    do
    {
      new_path.push_back(m_map.template beta<2>(d1));
      d1=m_map.template beta<0,2,0>(d1);
    }
    while(d1!=d2);
  }

  void bracket_flattening()
  {
    std::vector<Dart_const_handle> new_path;

    for (std::size_t i=0; i<m_path.size()-1; )
    {
      if (next_turn(i)!=1 /* && next_negative_turn(i)!=1 */)
      {
        new_path.push_back(m_path[i]); // We copy this dart
        ++i;
      }
      else
      {
        // i is maybe the beginning of a bracket
        std::size_t begin=i;
        std::size_t end=find_end_of_braket(begin);
        std::cout<<"Bracket: ["<<begin<<"; "<<end<<"]"<<std::endl;
        if (end<begin)
        { i=m_path.size(); }
        else
        { i=end; }
        transform_brancket(begin, end, new_path);
      }
    }
    new_path.swap(m_path);
  }

protected:
  const Map& m_map;
  std::vector<Dart_const_handle> m_path;
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_H //
// EOF //
