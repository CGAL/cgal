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
#include <vector>

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
  typedef std::list<std::pair<Dart_const_handle, int> > List_of_dart_length;
  typedef typename List_of_dart_length::iterator List_iterator;
  typedef typename List_of_dart_length::const_iterator List_const_iterator;

  typedef Path_on_surface_with_rle<Map> Self;

  Path_on_surface_with_rle(const Map& amap) : m_map(amap),
                                              m_is_closed(false),
                                              m_length(0)
  {}

  Path_on_surface_with_rle(const Path_on_surface<Map>& apath) : m_map(apath.get_map()),
                                                                m_is_closed(apath.is_closed()),
                                                                m_length(apath.length())
  {
    if (apath.is_empty()) return;

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
        CGAL_assertion(length>0);
        m_path.push_back(std::make_pair(apath[i], positive_flat?length:-length)); // begining of the flat part
        i=j;
      }
    }
    while(i<apath.length() && i!=starti);
  }
  
  void swap(Self& p2)
  {
    CGAL_assertion(&m_map==&(p2.m_map));
    m_path.swap(p2.m_path);
    std::swap(m_is_closed, p2.m_is_closed);
    std::swap(m_length, p2.m_length);
  }

  Self& operator=(const Self& other)
  {
    CGAL_assertion(&m_map==&(other.m_map));
    m_path=other.m_path;
    m_is_closed=other.m_is_closed;
    m_length=other.m_length;
    return *this;
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

  void advance_iterator(List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    it=std::next(it);
    if (is_closed() && it==m_path.end())
    { it=m_path.begin(); } // Here the path is closed, and it is the last element of the list
  }

  void advance_iterator(List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    it=std::next(it);
    if (is_closed() && it==m_path.end())
    { it=m_path.begin(); } // Here the path is closed, and it is the last element of the list
  }

  void retreat_iterator(List_iterator& it)
  {
    CGAL_assertion(it!=m_path.begin() || is_closed());
    if (is_closed() && it==m_path.begin())
    { it=m_path.end(); }
    it=std::prev(it);
  }

  void retreat_iterator(List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(it!=m_path.begin() || is_closed());
    if (is_closed() && it==m_path.begin())
    { it=m_path.end(); }
    it=std::prev(it);
  }

  List_iterator next_iterator(const List_iterator& it)
  {
    List_iterator res=it;
    advance_iterator(res);
    return res;
  }

  List_const_iterator next_iterator(const List_const_iterator& it) const
  {
    List_const_iterator res=it;
    advance_iterator(res);
    return res;
  }

  List_iterator prev_iterator(const List_iterator& it)
  {
    List_iterator res=it;
    retreat_iterator(res);
    return res;
 }

  List_const_iterator prev_iterator(const List_const_iterator& it) const
  {
    List_const_iterator res=it;
    retreat_iterator(res);
    return res;
 }

  // @return true iff there is a dart after it
  bool next_dart_exist(const List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return is_closed() || std::next(it)!=m_path.end();
  }

  // @return true iff 'it' is the beginning of a flat part (possibly of null length)
  // In fact, return false only if 'it' is the second dart of a flat part of non
  // null length.
  bool is_beginning_of_flat(const List_const_iterator& it) const
  {
    if (it->second!=0)
    { return true; } // Only the beginning of a flat part has a non null length

    if (!is_closed() && it==m_path.begin())
    { return true; }

    return prev_iterator(it)->second==0;
  }

  // @return true iff 'it' is the beginning of a non null flat part
  bool is_beginning_of_non_null_flat(const List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(it->second==0 || next_dart_exist(it));

    return (it->second!=0);
  }

  // @return true iff 'it' is the end of a flat part (possibly of null length)
  // In fact, return false only if 'it' is the first dart of a flat part of non
  // null length.
  bool is_end_of_flat(const List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return (it->second==0);
  }

  // @return the second dart of the given flat.
  Dart_const_handle second_dart_of_a_flat(const List_const_iterator& it) const
  {
    CGAL_assertion(is_beginning_of_non_null_flat(it));

    if (it->second>0)
    { return m_map.template beta<1, 2, 1>(it->first); }

    // here it->second<0
    return m_map.template beta<2, 0, 2, 0, 2>(it->first);
  }

  // @return the dart before the last dart of the given flat.
  Dart_const_handle before_last_dart_of_a_flat(const List_const_iterator& it) const
  {
    CGAL_assertion(is_beginning_of_non_null_flat(it));

    if (it->second>0)
    { return m_map.template beta<0, 2, 0>(next_iterator(it)->first); }

    // here it->second<0
    return m_map.template beta<2, 1, 2, 1, 2>(next_iterator(it)->first);
  }

  // Reduce the length of the flat part starting at 'it' from its beginning
  // 'it' moves to the end of the previous flat part if the current flat part
  // If the flat was empty, remove 'it'.
  // If the flat part becomes empty, remove its second element.
  // The path could be not valid after this operation (consistency with next
  // element should be ensure, by possibly updating the next flat part).
  void reduce_flat_from_beginning(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    if (it->second==0)
    {
      it=m_path.erase(it);
      if (!m_path.empty()) { retreat_iterator(it); } // TODO if !is_closed && first dart
    }
    else
    {
      it->first=second_dart_of_a_flat(it);

      if (it->second>0) { --(it->second); }
      else              { ++(it->second); }

      if (it->second==0)
      {
        m_path.erase(next_iterator(it)); // erase the second element of the flat
      }
      else { it=next_iterator(it); }
    }
  }

  // Reduce the length of the flat part starting at 'it' from its end
  // 'it' moves to the end of the previous flat part if the current flat part
  // becomes empty (and thus was removed) or to the end of this part otherwise
  // If the flat was empty, remove 'it'.
  // If the flat part becomes empty, remove its second element.
  // The path could be not valid after this operation (consistency with next
  // element should be ensure, by possibly updating the next flat part).
  void reduce_flat_from_end(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    if (it->second==0)
    {
      it=m_path.erase(it);
      if (!m_path.empty()) { retreat_iterator(it); } // TODO if !is_closed && first dart
    }
    else
    {
      next_iterator(it)->first=before_last_dart_of_a_flat(it);

      if (it->second>0) { --(it->second); }
      else              { ++(it->second); }

      if (it->second==0)
      {
        m_path.erase(next_iterator(it)); // erase the second element of the flat
      }
      else { it=next_iterator(it); }
    }
  }

  void merge_adjacent_flats_if_possible(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));
    bool positive1=false, negative1=false;
    bool positive2=false, negative2=false;
    bool positive3=false, negative3=false;

    List_iterator it2=next_iterator(it);
    if (!is_beginning_of_flat(it2)) { advance_iterator(it2); }

    if (it==it2) { return; } // only one flat in the path

    if (it->second>0) positive1=true;
    else if (it->second<0) negative1=true;

    List_iterator it1_end=it;
    if (!is_end_of_flat(it)) { advance_iterator(it1_end); }

    if (next_positive_turn(it1_end)==2) { positive2=true; }
    if (next_negative_turn(it1_end)==2) { negative2=true; }

    if (!positive2 && !negative2) { return; } // the middle turn is not a flat

    if (it2->second>0) positive3=true;
    else if (it2->second<0) negative3=true;

    if (!positive1 && !negative1)
    { // First flat is empty (length 0)
      if (!positive3 && !negative3)
      { // and second flat too
        it->second=(positive2?+1:-1); // we only initialize the length of the beginning of the flat
      }
      else
      { // here second flat is not empty => we have 3 darts (1 for first flat, 2 for second flat)
        if ((positive3 && !positive2) || (negative3 && !negative2))
        { return; } // the two flats cannot be merged

        it->second=it2->second+(positive3?1:-1);
        m_path.erase(it2);
      }
    }
    else
    { // here first flat is not empty
      if (!positive3 && !negative3)
      { // Second flat is empty =>  we have 3 darts (2 for first flat, 1 for second flat)
        if ((positive1 && !positive2) || (negative1 && !negative2))
        { return; } // the two flats cannot be merged

        it->second+=(positive1?1:-1);
        m_path.erase(it1_end);
      }
      else
      {
        // Second flat is not empty =>  we have 4 darts (2 for first flat, 2 for second flat)
        if ((positive1!=positive3) || (negative1!=negative3) ||
            (positive1 && !positive2) || (negative1 && !negative2))
        { return; } // the two flats cannot be merged

        it->second+=(it2->second+(positive1?1:-1));
        m_path.erase(it2);
        m_path.erase(it1_end);
      }
    }
  }

  // Return true if it is the beginning of a spur.
  bool is_spur(const List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return it->second==0 &&
        next_dart_exist(it) &&
        m_map.template beta<2>(it->first)==next_iterator(it)->first;
  }

  // Remove the spur given by 'it'.
  // Either 'it' stay on the same element (if the flat still exist),
  // or 'it' move to the previous element if the flat containing the spur is removed.
  // ('it' will be equal to m_path.end() if the path becomes empty).
  // 'it' is necessary either the beginning of a null length flat, or the end
  // of a non null flat.
  void remove_spur(List_iterator& it)
  {
    CGAL_assertion(is_spur(it));

    List_iterator it2=next_iterator(it); // it2 is the beginning of the next flat

    if (!is_beginning_of_flat(it)) { retreat_iterator(it); }
    reduce_flat_from_end(it);

    CGAL_assertion(is_beginning_of_flat(it2));
    reduce_flat_from_beginning(it2);

    // Here it is on the flat before the

    // Now move it to the element before the removed spur
    // except if the path has become empty, or if it is not closed
    // and we are in its first element.
    if (m_path.empty())
    {
      it=m_path.end();
      m_is_closed=false;
    }
    else
    {
      CGAL_assertion(it!=m_path.end());
      // here 'it' is the end of the flat
      if (!is_beginning_of_flat(it)) { retreat_iterator(it); }
      merge_adjacent_flats_if_possible(it);
    }
  }

  // Move it to the next spur after it. Go to m_path.end() if there is no
  // spur in the path.
  void move_to_next_spur(List_iterator& it)
  {
     CGAL_assertion(it!=m_path.end());
    List_iterator itend=(is_closed()?it:m_path.end());
    do
    {
      advance_iterator(it);
      if (is_spur(it)) { return; }
    }
    while(it!=itend);
    it=m_path.end(); // Here there is no spur in the whole path
  }

  // Simplify the path by removing all spurs
  // @return true iff at least one spur was removed
  bool remove_spurs()
  {
    bool res=false;
    List_iterator it=m_path.begin();
    while(it!=m_path.end())
    {
      if (is_spur(it)) { remove_spur(it); res=true; }
      else { move_to_next_spur(it); }
    }
    return res;
  }

  /// @return the turn between dart it and next dart.
  ///         (turn is position of the second edge in the cyclic ordering of
  ///          edges starting from the first edge around the second extremity
  ///          of the first dart)
  std::size_t next_positive_turn(const List_const_iterator& it) const
  {
    CGAL_assertion(is_valid());
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(is_closed() || std::next(it)!=m_path.end());

    Dart_const_handle d1=it->first;
    if (it->second!=0)
    {
      if (it->second<0) { return -(it->second); }
      else { return it->second; }
    }

    Dart_const_handle d2=next_iterator(it)->first;
    CGAL_assertion(d1!=d2);

    if (d2==m_map.template beta<2>(d1))
    { return 0; }

    std::size_t res=1;
    while (m_map.template beta<1>(d1)!=d2)
    {
      ++res;
      d1=m_map.template beta<1, 2>(d1);
    }
    // std::cout<<"next_positive_turn="<<res<<std::endl;
    return res;
  }

  /// Same than next_positive_turn but turning in reverse orientation around vertex.
  std::size_t next_negative_turn(const List_const_iterator& it) const
  {
    CGAL_assertion(is_valid());
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(is_closed() || std::next(it)!=m_path.end());

    Dart_const_handle d1=m_map.template beta<2>(it->first);
    if (it->second!=0)
    {
      if (it->second<0) { return -(it->second); }
      else { return it->second; }
    }

    Dart_const_handle d2=m_map.template beta<2>(next_iterator(it)->first);
    CGAL_assertion(d1!=d2);

    if (d2==m_map.template beta<2>(d1))
    { return 0; }

    std::size_t res=1;
    while (m_map.template beta<0>(d1)!=d2)
    {
      ++res;
      d1=m_map.template beta<0, 2>(d1);
    }
    // std::cout<<"next_negative_turn="<<res<<std::endl;
    return res;
  }

  std::vector<std::size_t> compute_positive_turns() const
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      if (is_closed() || std::next(it)!=m_path.end())
      { res.push_back(next_positive_turn(it)); }
    }
    return res;
  }

  std::vector<std::size_t> compute_negative_turns() const
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      if (is_closed() || std::next(it)!=m_path.end())
      { res.push_back(next_negative_turn(it)); }
    }
    return res;
  }

  std::vector<std::size_t> compute_turns(bool positive) const
  { return (positive?compute_positive_turns():compute_negative_turns()); }

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
  List_of_dart_length m_path; // The sequence of turning darts, plus the length of the flat part after the dart (a flat part is a sequence of dart with positive turn == 2). If negative value k, -k is the length of the flat part, for negative turns (-2).
  bool m_is_closed; // True iff the path is a cycle
  std::size_t m_length;
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_WITH_RLE_H //
// EOF //
