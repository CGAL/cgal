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
#ifndef CGAL_PATH_ON_SURFACE_WITH_RLE_H
#define CGAL_PATH_ON_SURFACE_WITH_RLE_H 1

#include <list>
#include <utility>
#include <iostream>
#include <iterator>

namespace CGAL {

template<typename Map_>
class Path_on_surface;

template<typename Map_>
class Path_on_surface_with_rle
{
  friend class Path_on_surface<Map_>;

public:
  typedef Map_ Map;
  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::Dart_const_handle Dart_const_handle;

  typedef Path_on_surface_with_rle<Map> Self;

  Path_on_surface_with_rle(const Map& amap) : m_map(amap),
                                              m_is_closed(false),
                                              m_length(0)
  {}

  Path_on_surface_with_rle(const Path_on_surface<Map>& apath) : m_map(apath.get_map()),
                                                                m_is_closed(apath.is_closed()),
                                                                m_length(apath.length())
  {
    std::size_t i=0, j=0, starti=0, length=0;
    bool positive_flat=false;
    bool negative_flat=false;    
    
    if (apath.is_closed())
    {
      while (apath.next_positive_turn(i)==2 || apath.next_negative_turn(i)==2)
      {
        i=apath.next_index(i);
        if (i==0) // Case of a closed path, made of only one flat part.
        {
          m_path.push_back(std::make_pair(apath.front(),
                                          apath.next_positive_turn(0)==2?
                                              apath.length():
                                             -apath.length()));
          return;
        }
      }
    }

    starti=i;
    do
    {
      // Here dart i is the beginning of a flat part (maybe of length 0)
      if (apath.next_positive_turn(i)==2)
      { positive_flat=true; negative_flat=false; }
      else if (apath.next_negative_turn(i)==2)
      { positive_flat=false; negative_flat=true; }
      else
      { positive_flat=false; negative_flat=false; }

      if (!positive_flat && !negative_flat)
      {
        m_path.push_back(std::make_pair(apath[i], 0));
        i=apath.next_index(i);
      }
      else
      {
        j=i;
        length=0;
        while ((positive_flat && apath.next_positive_turn(j)==2) ||
               (negative_flat && apath.next_negative_turn(j)==2))
        {
          j=apath.next_index(j);
          ++length;
        }
        assert(length>0);
        m_path.push_back(std::make_pair(apath[i], positive_flat?length:-length)); // begining of the flat part
        i=j;
      }
    }
    while(i<apath.length() && i!=starti);
  }
  
  void swap(Self& p2)
  {
    assert(&m_map==&(p2.m_map));
    m_path.swap(p2.m_path);
    std::swap(m_is_closed, p2.m_is_closed);
    std::swap(m_length, p2.m_length);
  }

  /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts.
  bool operator==(const Self& other) const
  {
    if (is_closed()!=other.is_closed() ||
        length()!=other.length())
      return false;
    
    // TODO TEST THE TWO PATHS
    return true;
  }
  bool operator!=(const Self&  other) const
  { return !(operator==(other)); }

  // @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  std::size_t length() const
  { return m_length; }

  std::size_t size_of_list() const
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
    m_length=0;
  }
  
  // @return true iff the path is valid; i.e. a sequence of edges two by
  //              two adjacent.
  bool is_valid() const
  {
    // TODO
    return true;
  }

  void display() const
  {
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      std::cout<<m_map.darts().index(it->first)<<"("<<it->second<<")";
      if (std::next(it)!=itend) { std::cout<<" "; }
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
  std::list<std::pair<Dart_const_handle, int> > m_path; // The sequence of turning darts, plus the length of the flat part after the dart (a flat part is a sequence of dart with positive turn == 2). If negative value k, -k is the length of the flat part, for negative turns (-2).
  bool m_is_closed; // True iff the path is a cycle
  std::size_t m_length;
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_WITH_RLE_H //
// EOF //
