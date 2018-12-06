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
#ifndef CGAL_PATH_ON_SURFACE_WITH_RLE_H
#define CGAL_PATH_ON_SURFACE_WITH_RLE_H 1

#include <list>
#include <utility>
#include <iostream>
#include <iterator>
#include <vector>

#define CGAL_PWRLE_TURN_V1  // Compute turns by turning (method of CMap)
// #define CGAL_PWRLE_TURN_V2  // Compute turns by using an id of darts, given by an hash-table (built and given by Surface_mesh_curve_topology)
// #define CGAL_PWRLE_TURN_V3  // Compute turns by using an id of darts, associated in Info of Darts (build by Surface_mesh_curve_topology)

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
  // typedef typename List_of_dart_length::const_iterator List_const_iterator;

  typedef Path_on_surface_with_rle<Map> Self;

#ifdef CGAL_PWRLE_TURN_V2
    typedef boost::unordered_map<Dart_const_handle, std::size_t> TDartIds;
#endif //CGAL_PWRLE_TURN_V2

  /// Constructor
  Path_on_surface_with_rle(const Map& amap
                         #ifdef CGAL_PWRLE_TURN_V2
                           , const TDartIds & darts_ids
                         #endif //CGAL_PWRLE_TURN_V2
                           )
    : m_map(amap),
      m_is_closed(false),
      m_length(0)
  #ifdef CGAL_PWRLE_TURN_V2
    , m_darts_ids(darts_ids)
  #endif //CGAL_PWRLE_TURN_V2
  {}

  /// Copy constructor
  Path_on_surface_with_rle(const Path_on_surface<Map>& apath
                         #ifdef CGAL_PWRLE_TURN_V2
                           , const TDartIds & darts_ids
                         #endif //CGAL_PWRLE_TURN_V2
                           ) :
    m_map(apath.get_map()),
    m_is_closed(apath.is_closed()),
    m_length(apath.length())
#ifdef CGAL_PWRLE_TURN_V2
  , m_darts_ids(darts_ids)
#endif //CGAL_PWRLE_TURN_V2
  {
    if (apath.is_empty()) return;

    std::size_t i=0, j=0, starti=0, length=0;
    bool positive_flat=false;
    bool negative_flat=false;    
    
    if (apath.is_closed())
    {
      if (apath.next_positive_turn(i)==2)
      { positive_flat=true; negative_flat=false; }
      else if (apath.next_negative_turn(i)==2)
      { positive_flat=false; negative_flat=true; }

      while ((positive_flat && apath.next_positive_turn(i)==2) ||
             (negative_flat && apath.next_negative_turn(i)==2))
      {
        i=apath.next_index(i);
        if (i==0) // Case of a closed path, made of only one flat part.
        {
          m_path.push_back(std::make_pair(apath.front(),
                                          positive_flat?(apath.length()-1):
                                                        -(apath.length())-1));
          if (apath.length()>1)
          { m_path.push_back(std::make_pair(apath.back(),0)); }
          return;
        }
      }
      // Here i is the last dart of a flat
      i=apath.next_index(i); // Now we are sure that i is the beginning of a flat
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
        m_path.push_back(std::make_pair(apath[j], 0));
        i=apath.next_index(j);
      }
    }
    while(i<apath.length() && i!=starti);
    CGAL_assertion(is_valid());
  }

  void swap(Self& p2)
  {
    if (this!=&p2)
    {
      CGAL_assertion(&m_map==&(p2.m_map));
#ifdef CGAL_PWRLE_TURN_V2
      CGAL_assertion(&m_darts_ids==&(p2.m_darts_ids));
#endif // CGAL_PWRLE_TURN_V2
      m_path.swap(p2.m_path);
      std::swap(m_is_closed, p2.m_is_closed);
      std::swap(m_length, p2.m_length);
    }
  }

  Self& operator=(const Self& other)
  {
    CGAL_assertion(&m_map==&(other.m_map));
#ifdef CGAL_PWRLE_TURN_V2
    CGAL_assertion(&m_darts_ids==&(other.m_darts_ids));
#endif // CGAL_PWRLE_TURN_V2
    if (this!=&other)
    {
      m_path=other.m_path;
      m_is_closed=other.m_is_closed;
      m_length=other.m_length;
    }
    return *this;
  }

  Self& operator+=(Self& other)
  {
    // Be careful to the special case when *this==other
    // this is the reason of the iend.
    for (List_iterator lit=other.m_path.begin(),
         litend=std::prev(other.m_path.end()); lit!=litend; ++lit)
    { m_path.push_back(*lit); }
    m_path.push_back(other.m_path.back()); // Last element
    update_is_closed();
    return *this;
  }

  Self operator+(const Self& other) const
  {
    Self res=*this;
    res+=other;
    return res;
  }  
  /// @Return true if this path is equal to other path, identifying dart begin
  ///         of this path with dart 'itother' in other path.
  bool are_same_paths_from(Self& other, List_iterator itother)
  {
    CGAL_assertion(itother!=other.m_path.end());
    CGAL_assertion(is_closed() || itother==other.m_path.begin());
    CGAL_assertion(is_closed()==other.is_closed());

    for(auto it=m_path.begin(); it!=m_path.end(); ++it)
    {
      if (it->first!=itother->first || it->second!=itother->second)
      { return false; }
      other.advance_iterator(itother);
    }
    return true;
  }

  /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts. Old quadratic version, new version
  ///         (operator==) use linear version based on Knuth, Morris, Pratt
  bool are_paths_equals(const Self& other)
  {
    if (is_closed()!=other.is_closed() ||
        length()!=other.length() ||
        size_of_list()!=other.size_of_list())
      return false;
    
    if (!is_closed())
    { return are_same_paths_from(other, other.m_path.begin()); }

    for(auto itother=other.m_path.begin(); itother!=other.m_path.end(); ++itother)
    {
      if (are_same_paths_from(other, itother))
      { return true; }
    }

    return false;
  }

    /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts.
  bool operator==(Self& other)
  {
    if (is_closed()!=other.is_closed() ||
        length()!=other.length() ||
        size_of_list()!=other.size_of_list())
      return false;

    if (is_empty() && other.is_empty())
    { return true; }

    if (!is_closed())
    { return are_same_paths_from(other, other.m_path.begin()); }

    // TODO: is it possible to avoid the transformations into Path_on_surface
    //  and use directly knuth_morris_pratt_search with Path_on_surface_with_rle ?
    Path_on_surface<Map> p1(*this);
    Path_on_surface<Map> p2(other);
    return p1==p2;
  }
  
  bool operator!=(Self&  other)
  { return !(operator==(other)); }

  /// @return true iff the path is empty
  bool is_empty() const
  { return m_path.empty(); }

  /// @return the length of the path, i.e. its number of darts
  std::size_t length() const
  { return m_length; }

  /// @return the number of darts in the double linked list.
  /// note that size_of_list()<=length().
  std::size_t size_of_list() const
  { return m_path.size(); }

  /// @return true iff the path is closed (update after each path modification).
  bool is_closed() const
  { return m_is_closed; }

  /// @return the underlying map.
  const Map& get_map() const
  { return m_map; }

  /// clear the path.
  void clear()
  {
    m_path.clear();
    m_is_closed=false;
    m_length=0;
  }

  /// @return true iff the path is valid; i.e. a sequence of edges two by
  ///              two adjacent.
  bool is_valid()
  {
    if (is_empty()) { return !is_closed(); } // an empty past is not closed

    unsigned int nbdarts=0,i=0;
    Dart_const_handle prev=NULL;
    bool loopend=false;
    auto it=m_path.begin();
    if (!is_beginning_of_flat(it)) { advance_iterator(it); }

    do
    {
      if (prev!=NULL &&
          !CGAL::belong_to_same_cell<Map, 0>
          (m_map, m_map.template beta<1>(prev), it->first))
      {
        std::cerr<<"ERROR: The path is not valid: dart in position "<<i
                 <<" does not follow the previous dart."<<std::endl;
        return false;
      }
      if (it->second!=0 && next_iterator(it)==m_path.end())
      {
        std::cerr<<"ERROR: The path is not valid: a non null flat does not "
                <<"have a second dart in position "<<i<<"."<<std::endl;
        return false;
      }
      if (it->second!=0 && next_iterator(it)->second!=0)
      {
        std::cerr<<"ERROR: The path is not valid: two non null flat are "
                <<"consecutive in position "<<i<<"."<<std::endl;
        return false;
      }
      if (it->second!=0)
      {
        Dart_const_handle dhend=it->first;
        int nb=0;
        while(nb!=it->second)
        {
          if (it->second>0)
          { dhend=m_map.template beta<1, 2, 1>(dhend); ++nb; }
          else
          { dhend=m_map.template beta<2, 0, 2, 0, 2>(dhend); --nb; }
          ++nbdarts;
        }
        advance_iterator(it); ++i; ++nbdarts;
        if (it==m_path.end() || it==m_path.begin()) loopend=true;

        if (dhend!=it->first)
        {
          std::cout<<"ERROR: The path is not valid: flat at position "<<i-1
                  <<" with length "<<nb<<" is not correct: its end does not"
                 <<" correspond to the flat."<<std::endl;
          return false;
        }
      }
      else
      { ++nbdarts; }
      prev=it->first; // end of the previous flat
      ++i; advance_iterator(it);
      if (it==m_path.end() || it==m_path.begin()) loopend=true;
    }
    while(!loopend);

    if (is_closed())
    {
      if (prev==NULL)
      {
        std::cout<<"ERROR: The path is not valid: it is empty and closed."
                <<std::endl;
        return false;
      }
      if (!CGAL::belong_to_same_cell<Map, 0>
          (m_map, m_map.template beta<1>(prev), it->first))
      {
        std::cerr<<"ERROR: The path is not valid: dart in position "<<i
                 <<" does not follow the previous dart."<<std::endl;
        return false;
      }
    }

    if (length()!=nbdarts)
    {
      std::cerr<<"ERROR: The path is not valid: store length=="<<length()
               <<" different of real length=="<<nbdarts<<std::endl;
      return false;
    }

    return true;
  }

  /// For debugging purpose, test if 'it' is a valid iterator.
  bool is_valid_iterator(const List_iterator& ittotest)
  {
    if (ittotest==m_path.end()) { return true; }
    for (auto it=m_path.begin(); it!=m_path.end(); ++it)
    {
      if (it==ittotest) { return true; }
    }
    return false;
  }

  Dart_const_handle front()
  {
    CGAL_assertion(!is_empty());
    return m_path.front().first;
  }

  Dart_const_handle back()
  {
    CGAL_assertion(!is_empty());
    return m_path.back().first;
  }

  /* TODO REMOVE
   void push_back(Dart_const_handle dh, bool update_isclosed=true)
  {
    CGAL_assertion(dh!=NULL && dh!=m_map.null_dart_handle);
    // This assert is too long...
    // CGAL_assertion((is_empty() ||
    //      CGAL::template belong_to_same_cell<Map, 0>
    //       (m_map, m_map.other_extremity(back()), dh)));

    m_path.push_back(dh);
    if (update_isclosed) { update_is_closed(); }
  }
  */
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
      { m_is_closed=CGAL::belong_to_same_cell<Map,0>(m_map, front(), pend); }
    }
  }

  /// @return true iff there is a dart after it
  bool next_dart_exist(const List_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return is_closed() || std::next(it)!=m_path.end();
  }

  /// @return true iff there is a dart before it
  bool prev_dart_exist(const List_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return is_closed() || it!=m_path.begin();
  }

  /// @return true iff there is a flat after the flat given by 'it'
  bool next_flat_exist(const List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    return next_dart_exist(it) &&
      (is_beginning_of_flat(next_iterator(it)) ||
       next_dart_exist(next_iterator(it)));
  }
  
  /// @return true iff there is a flat before the flat given by 'it'
  bool previous_flat_exist(const List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    return prev_dart_exist(it) &&
      (is_beginning_of_flat(prev_iterator(it)) ||
       prev_dart_exist(prev_iterator(it)));
  }
  
  /// @return true iff 'it' is the beginning of a flat part (possibly of null
  ///   length). In fact, returns false only if 'it' is the second dart of a
  ///   non length null flat.
  bool is_beginning_of_flat(const List_iterator& it)
  {
    if (it->second!=0)
    { return true; } // Only the beginning of a flat part has a non null length

    if (!is_closed() && it==m_path.begin())
    { return true; }

    return prev_iterator(it)->second==0;
  }

  /// @return true iff 'it' is the beginning of a non length null flat
  bool is_beginning_of_non_null_flat(const List_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(it->second==0 || next_dart_exist(it));

    return (it->second!=0);
  }

  /// @return true iff 'it' is the end of a flat part (possibly of null length).
  /// In fact, return false only if 'it' is the first dart of a non
  /// length null flat.
  bool is_end_of_flat(const List_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return (it->second==0);
  }

/*  void advance_iterator(List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    it=std::next(it);
    if (is_closed() && it==m_path.end())
    { it=m_path.begin(); } // Here the path is closed, and it is the last element of the list
  } */

  void advance_iterator(List_iterator& it)
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

  /*void retreat_iterator(List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(it!=m_path.begin() || is_closed());
    if (is_closed() && it==m_path.begin())
    { it=m_path.end(); }
    it=std::prev(it);
  } */

  void advance_to_next_flat(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));
    advance_iterator(it);
    if (it!=m_path.end() && !is_beginning_of_flat(it))
    { advance_iterator(it); }
  }

  /* void advance_to_next_flat(List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    advance_iterator(it);
    if (it!=m_path.end() && !is_beginning_of_flat(it))
    { advance_iterator(it); }
  } */

  void retreat_to_prev_flat(List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    retreat_iterator(it);
    if (!is_beginning_of_flat(it))
    { retreat_iterator(it); }
  }

  /* void retreat_to_prev_flat(List_const_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    retreat_iterator(it);
    if (!is_beginning_of_flat(it))
    { retreat_iterator(it); }
  } */

  List_iterator next_iterator(const List_iterator& it)
  {
    List_iterator res=it;
    advance_iterator(res);
    return res;
  }

  /* List_const_iterator next_iterator(const List_const_iterator& it) const
  {
    List_const_iterator res=it;
    advance_iterator(res);
    return res;
  } */

  List_iterator prev_iterator(const List_iterator& it)
  {
    List_iterator res=it;
    retreat_iterator(res);
    return res;
 }

  /* List_const_iterator prev_iterator(const List_const_iterator& it) const
  {
    List_const_iterator res=it;
    retreat_iterator(res);
    return res;
 } */

  List_iterator next_flat(const List_iterator& it)
  {
    List_iterator res=it;
    advance_to_next_flat(res);
    return res;
  }

  /* List_const_iterator next_flat(const List_const_iterator& it) const
  {
    List_iterator res=it;
    advance_to_next_flat(res);
    return res;
  } */

  List_iterator prev_flat(const List_iterator& it)
  {
    List_iterator res=it;
    retreat_to_prev_flat(res);
    return res;
  }

  /* List_const_iterator prev_flat(const List_const_iterator& it) const
  {
    List_iterator res=it;
    retreat_to_prev_flat(res);
    return res;
  } */

  /// @return the length of the given flat.
  int flat_length(const List_iterator& it) const
  {
    CGAL_assertion(is_beginning_of_flat(it));
    return it->second;
  }

  /// @return the second dart of the given flat.
  Dart_const_handle second_dart_of_a_flat(const List_iterator& it) const
  {
    CGAL_assertion(is_beginning_of_non_null_flat(it));

    if (it->second>0)
    { return m_map.template beta<1, 2, 1>(it->first); }

    CGAL_assertion(it->second<0);
    return m_map.template beta<2, 0, 2, 0, 2>(it->first);
  }

  /// @return the dart before the last dart of the given flat.
  Dart_const_handle before_last_dart_of_a_flat(const List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_non_null_flat(it));

    if (it->second>0)
    { return m_map.template beta<0, 2, 0>(next_iterator(it)->first); }

    CGAL_assertion(it->second<0);
    return m_map.template beta<2, 1, 2, 1, 2>(next_iterator(it)->first);
  }

  /// @return the positive turn given two ids of darts (unsed for CGAL_PWRLE_TURN_V1)
  std::size_t compute_positive_turn_given_ids(std::size_t id1,
                                               std::size_t id2) const
  {
    std::size_t number_of_edges=m_map.number_of_darts()/2;
    if (id1>=number_of_edges)
    {
      id1-=number_of_edges; // id of the first dart in its own vertex
      CGAL_assertion(id2>=number_of_edges);
      id2-=number_of_edges; // id of the second dart in its own vertex
    }

    return (id1<id2?id2-id1:number_of_edges-id1+id2);
  }

  /// @return the negative turn given two ids of darts (unused for CGAL_PWRLE_TURN_V1)
  std::size_t compute_negative_turn_given_ids(std::size_t id1,
                                              std::size_t id2) const
  {
    std::size_t number_of_edges=m_map.number_of_darts()/2;
    if (id1>=number_of_edges)
    {
      id1-=number_of_edges; // id of the first dart in its own vertex
      CGAL_assertion(id2>=number_of_edges);
      id2-=number_of_edges; // id of the second dart in its own vertex
    }

    return (id1<=id2?number_of_edges-id2+id1:id1-id2);
  }

  std::size_t get_dart_id(Dart_const_handle dh) const
  {
#ifdef CGAL_PWRLE_TURN_V2
    return m_darts_ids.at(dh);
#else // CGAL_PWRLE_TURN_V2
#ifdef CGAL_PWRLE_TURN_V3
    return m_map.info(dh);
#else //  CGAL_PWRLE_TURN_V3
    std::cerr<<"Error: impossible to get dart id without method V2 or V3."<<std::endl;
    return -1;
#endif  CGAL_PWRLE_TURN_V3
#endif //  CGAL_PWRLE_TURN_V2
  }

  /// @return the turn between dart it and next dart.
  ///         (turn is position of the second edge in the cyclic ordering of
  ///          edges starting from the first edge around the second extremity
  ///          of the first dart)
  std::size_t next_positive_turn(const List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(is_closed() || std::next(it)!=m_path.end());

    if (it->second!=0)
    {
      if (it->second<0) { return -2; }
      else { return 2; }
    }

#ifdef CGAL_PWRLE_TURN_V1
    return m_map.positive_turn(it->first, next_iterator(it)->first);
#else // CGAL_PWRLE_TURN_V1
    return compute_positive_turn_given_ids(get_dart_id(m_map.template beta<2>(it->first)),
                                           get_dart_id(next_iterator(it)->first));
#endif // CGAL_PWRLE_TURN_V1
  }

  /// Same than next_positive_turn but turning in reverse orientation around vertex.
  std::size_t next_negative_turn(const List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(it!=m_path.end());
    CGAL_assertion(is_closed() || std::next(it)!=m_path.end());

    if (it->second!=0)
    {
      if (it->second<0) { return -2; }
      else { return 2; }
    }

#ifdef CGAL_PWRLE_TURN_V1
    return m_map.negative_turn(it->first, next_iterator(it)->first);
#else // CGAL_PWRLE_TURN_V1
    return compute_negative_turn_given_ids(get_dart_id(m_map.template beta<2>(it->first)),
                                           get_dart_id(next_iterator(it)->first));
#endif // CGAL_PWRLE_TURN_V1
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

  /// Reduce the length of the flat part starting at 'it' from its beginning
  /// 'it' moves to the end of the previous flat part if the current flat part
  /// If the flat was empty, remove 'it'.
  /// If the flat part becomes empty, remove its second element.
  /// The path could be not valid after this operation (consistency with next
  /// element should be ensure, by possibly updating the next flat part).
  void reduce_flat_from_beginning(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    if (it->second==0)
    {
      it=m_path.erase(it);
      if (!m_path.empty()) { retreat_iterator(it); } // TODO ? if !is_closed && first dart
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
    CGAL_assertion(m_length>0);
    --m_length;
  }

  /// Reduce the length of the flat part starting at 'it' from its end
  /// 'it' moves to the end of the previous flat part if the current flat part
  /// becomes empty (and thus was removed) or to the end of this part otherwise
  /// If the flat was empty, remove 'it'.
  /// If the flat part becomes empty, remove its second element.
  /// The path could be not valid after this operation (consistency with next
  /// element should be ensure, by possibly updating the next flat part).
  void reduce_flat_from_end(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    if (it->second==0)
    {
      it=m_path.erase(it);
      if (!m_path.empty()) { retreat_iterator(it); } // TODO ? if !is_closed && first dart
    }
    else
    {
      if (it->second>0) { --(it->second); }
      else              { ++(it->second); }

      List_iterator it2=next_iterator(it);
      if (it->second==0)
      {
        m_path.erase(it2); // erase the second element of the flat
      }
      else
      {
        it2->first=before_last_dart_of_a_flat(it);
        it=it2;
      }
    }
    CGAL_assertion(m_length>0);
    --m_length;
  }

  void merge_adjacent_flats_if_possible(List_iterator& it)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    bool positive1=false, negative1=false;
    bool positive2=false, negative2=false;
    bool positive3=false, negative3=false;

    List_iterator it1_end=it;
    if (!is_end_of_flat(it1_end)) { advance_iterator(it1_end); }
    List_iterator it2=next_iterator(it1_end);

    if (it==it2) { return; } // only one flat in the path

    if (next_positive_turn(it1_end)==2) { positive2=true; }
    if (next_negative_turn(it1_end)==2) { negative2=true; }

    if (!positive2 && !negative2) { return; } // the middle turn is not a flat

    if (it->second>0) positive1=true;
    else if (it->second<0) negative1=true;

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
        m_path.erase(it1_end);
        m_path.erase(it2);
      }
    }

    // CGAL_assertion(is_valid());
  }

  /// Return true if it is the beginning of a spur.
  bool is_spur(const List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());
    return it->second==0 &&
        next_dart_exist(it) &&
        m_map.template beta<2>(it->first)==next_iterator(it)->first;
  }

  /// Remove the spur given by 'it'.
  /// Either 'it' stay on the same element (if the flat still exist),
  /// or 'it' move to the previous element if the flat containing the spur is removed.
  /// ('it' will be equal to m_path.end() if the path becomes empty).
  /// 'it' is necessary either the beginning of a null length flat, or the end
  /// of a non null flat.
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

  /// Move it to the next spur after it. Go to m_path.end() if there is no
  /// spur in the path.
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

  /// Simplify the path by removing all spurs, if all==true; otherwise
  /// remove only one spur.
  /// @return true iff at least one spur was removed
  bool remove_spurs(bool all=true)
  {
    bool res=false;
    List_iterator it=m_path.begin();
    while(it!=m_path.end())
    {
      if (is_spur(it))
      {
        remove_spur(it); res=true;
        if (!all) { return true; }
      }
      else { move_to_next_spur(it); }
    }
    CGAL_assertion(is_valid());
    return res;
  }

  /// Return true if 'it' is the beginning of a positive bracket.
  /// If true, 'itend' is updated to be the end of the bracket
  bool is_positive_bracket(const List_iterator& it,
                           List_iterator& itend)
  {
    CGAL_assertion(it!=m_path.end());
    if (it->second!=0 || !next_dart_exist(it) || next_positive_turn(it)!=1)
    { return false; }
    // Here it is the end of a first flat, and first turn is 1

    itend=next_iterator(it);
    // Here itend is the beginning of the second flat
    if (itend->second<0)
    { return false; } // This is not a positive bracket, we have some -2

    if (!is_end_of_flat(itend)) { advance_iterator(itend); }
    // Here itend is the end of the second flat
    if (!next_dart_exist(itend) || next_positive_turn(itend)!=1)
    { return false; }

    // Now we move itend to the beginning of the third flat.
    advance_iterator(itend);
    return true;
  }
  
  /// Return true if it is the beginning of a negative bracket.
  /// If true, itend is updated to be the end of the bracket
  bool is_negative_bracket(const List_iterator& it,
                           List_iterator& itend)
  {
    CGAL_assertion(it!=m_path.end());
    if (it->second!=0 || !next_dart_exist(it) || next_negative_turn(it)!=1)
    { return false; }
    // Here it is the end of a first flat, and first negative turn is 1

    itend=next_iterator(it);
    // Here itend is the beginning of the second flat
    if (itend->second>0)
    { return false; } // This is not a negative bracket, we have some +2

    if (!is_end_of_flat(itend)) { advance_iterator(itend); }
    // Here itend is the end of the second flat
    if (!next_dart_exist(itend) || next_negative_turn(itend)!=1)
    { return false; }

    // Now we move itend to the beginning of the third flat.
    advance_iterator(itend);
    return true;
  }

  /// Move 'it' to the next bracket after 'it'. Go to m_path.end() if there is no
  /// bracket in the path.
  void move_to_next_bracket(List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    CGAL_assertion(it!=m_path.end());
    List_iterator itend=(is_closed()?it:m_path.end());
    List_iterator it2;
    do
    {
      advance_iterator(it);
      if (is_positive_bracket(it, it2) || is_negative_bracket(it, it2))
      { return; }
    }
    while(it!=itend);
    it=m_path.end(); // Here there is no spur in the whole path
  }

  /// Remove the given negative bracket. it1 is the end of a flat beginning
  /// of the bracket; it3 is the beginning of the flat end of the bracket.
  void remove_negative_bracket(List_iterator& it1, List_iterator& it3)
  {
    List_iterator it2;
    CGAL_assertion(is_negative_bracket(it1, it2)); // Here it2 is unused

    if (next_iterator(it1)==it3)
    { // Case of cyclic bracket
      CGAL_assertion(it3->second<0);
      CGAL_assertion(next_iterator(it3)==it1);
      it3->first=m_map.template  beta<2,0,2,1>(it3->first);
      it1->first=m_map.template  beta<2,1,2,0>(it1->first);
      it3->second=(-it3->second)-2;
      CGAL_assertion(m_length>1);
      m_length-=2;
      return;
    }
    
    it2=next_iterator(it1); // it2 is the beginning of the next flat after it1
    if (!is_beginning_of_flat(it1)) { retreat_iterator(it1); }    

    assert(is_beginning_of_flat(it1));
    assert(is_beginning_of_flat(it2));
    assert(is_beginning_of_flat(it3));
    
    reduce_flat_from_end(it1); // decrease also m_length
    reduce_flat_from_beginning(it3);

    it2->first=m_map.template  beta<2,1,1>(it2->first);
    if (it2->second<0)
    {
      it2->second=-it2->second;
      advance_iterator(it2);
      it2->first=m_map.template  beta<2,1,1>(it2->first);
    }

    if (is_beginning_of_non_null_flat(it1)) { advance_iterator(it1); }

    // CGAL_assertion(is_valid());
  }

  /// Remove the given positive bracket. it1 is the end of a flat beginning
  /// of the bracket; it3 is the beginning of the flat end of the bracket.
  void remove_positive_bracket(List_iterator& it1,
                               List_iterator& it3)
  {
    List_iterator it2;
    CGAL_assertion(is_positive_bracket(it1, it2)); // Here it2 is unused

    if (next_iterator(it1)==it3)
    { // Case of cyclic bracket
      CGAL_assertion(it3->second>0);
      CGAL_assertion(next_iterator(it3)==it1);
      it3->first=m_map.template  beta<1,2,0,2>(it3->first);
      it1->first=m_map.template  beta<0,2,1,2>(it1->first);
      it3->second=-(it3->second-2);
      m_length-=2;
      return;
    }

    it2=next_iterator(it1); // it2 is the beginning of the next flat after it1
    if (!is_beginning_of_flat(it1)) { retreat_iterator(it1); }    

    assert(is_beginning_of_flat(it1));
    assert(is_beginning_of_flat(it2));
    assert(is_beginning_of_flat(it3));
    
    reduce_flat_from_end(it1); // decrease also m_length
    reduce_flat_from_beginning(it3);

    it2->first=m_map.template  beta<1,1,2>(it2->first);
    if (it2->second>0)
    {
      it2->second=-it2->second;
      advance_iterator(it2);
      it2->first=m_map.template  beta<1,1,2>(it2->first);
    }

    if (is_beginning_of_non_null_flat(it1)) { advance_iterator(it1); }

    // CGAL_assertion(is_valid());
  }
  
  /// Simplify the path by removing all brackets, if all==true (default),
  /// or by removing only one bracket, if all==false.
  /// @return true iff at least one bracket was removed
  bool remove_brackets(bool all=true)
  {
    bool res=false;
    List_iterator it1=m_path.begin();
    List_iterator it2;
    while(it1!=m_path.end())
    {
      // CGAL_assertion(is_valid());
      if (is_positive_bracket(it1, it2))
      { remove_positive_bracket(it1, it2); res=true; }
      else if (is_negative_bracket(it1, it2))
      { remove_negative_bracket(it1, it2); res=true; }
      else { move_to_next_bracket(it1); }
      if (!all && res) { return true; }
      // CGAL_assertion(is_valid());
    }
    CGAL_assertion(is_valid());
    return res;
  }

  /// @return true iff 'it' is the beginning of a l-shape.
  bool is_l_shape(const List_iterator& it)
  {
    CGAL_assertion(it!=m_path.end());

    if (!is_beginning_of_flat(it) || it->second>0)
    { return false; }
    
    // Here it->second<=0
    List_iterator it2=(it->second==0?it:next_iterator(it));
    // Here it2 is the end of the first flat
    if (next_negative_turn(it2)!=1)
    { return false; }

    advance_iterator(it2); // Here it2 is the beginning of the second flat
    if (it2->second>0) { return false; }
    
    return true;
  }

  /// Move it to the next l_shape after it. Go to m_path.end() if there is no
  /// l_shape in the path.
  void move_to_next_l_shape(List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(it!=m_path.end());
    List_iterator itend=(is_closed()?it:m_path.end());
    do
    {
      advance_iterator(it);
      if (is_l_shape(it)) { return; }
    }
    while(it!=itend);
    it=m_path.end(); // Here there is no spur in the whole path
  }

  /// @return true iff the flat before flat 'it' can be extended by adding
  ///              dart 'dh' to its end.
  bool is_prev_flat_can_be_extended_at_end(const List_iterator& it,
                                           Dart_const_handle dh,
                                           bool& positive_flat,
                                           bool& negative_flat)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    positive_flat=false; negative_flat=false;
    if (!prev_dart_exist(it)) { return false; }
    List_iterator ittemp=prev_iterator(it);
    if (m_map.positive_turn(ittemp->first, dh)==2) { positive_flat=true; }
    if (m_map.negative_turn(ittemp->first, dh)==2) { negative_flat=true; }

    if (is_beginning_of_flat(ittemp)) // Case of flat lengh 0
    { return positive_flat || negative_flat; }

    retreat_iterator(ittemp);
    return (ittemp->second>0 && positive_flat) ||
        (ittemp->second<0 && negative_flat);
  }

  /// @return true iff the flat after flat 'it' can be extended by adding
  ///              dart 'dh' to its beginning.
  bool is_next_flat_can_be_extended_at_beginning(const List_iterator& it,
                                                 Dart_const_handle dh,
                                                 bool& positive_flat,
                                                 bool& negative_flat)
  {
     CGAL_assertion(is_end_of_flat(it));

     positive_flat=false; negative_flat=false;
     if (!next_dart_exist(it)) { return false; }
     List_iterator ittemp=next_iterator(it);
     if (m_map.positive_turn(dh, ittemp->first)==2) { positive_flat=true; }
     if (m_map.negative_turn(dh, ittemp->first)==2) { negative_flat=true; }

     if (ittemp->second==0)  // Case of flat lengh 0
     { return positive_flat || negative_flat; }

     return (ittemp->second>0 && positive_flat) ||
         (ittemp->second<0 && negative_flat);
  }
  
  /// Add the given dart 'dh' before the flat 'it'.
  void add_dart_before(const List_iterator& it, Dart_const_handle dh)
  {
    CGAL_assertion(is_beginning_of_flat(it));

    bool positive_flat, negative_flat;
    if (is_prev_flat_can_be_extended_at_end(it, dh,
                                            positive_flat, negative_flat))
    {
      List_iterator ittemp=prev_iterator(it);
      if (is_beginning_of_flat(ittemp))
      {
        ittemp->second=(negative_flat?-1:+1);
        m_path.insert(it, std::make_pair(dh, 0));
        // insert the new element before 'it', thus after 'ittemp'
      }
      else
      {
        ittemp->first=dh; // Move the last dart of the previous flat
        retreat_iterator(ittemp);
        CGAL_assertion(ittemp->second!=0);
        ittemp->second+=(ittemp->second>0?+1:-1); // Increment the length of the flat
      }
    }
    else
    {
      m_path.insert(it, std::make_pair(dh, 0));
      // insert the new element before 'it'
    }
    ++m_length;
  }
  
  /// Add the given dart 'dh' after the flat 'it' (given by its end).
  void add_dart_after(const List_iterator& it, Dart_const_handle dh)
  {
    CGAL_assertion(is_end_of_flat(it));
    bool positive_flat, negative_flat;
    List_iterator ittemp=next_iterator(it);
    if (is_next_flat_can_be_extended_at_beginning(it, dh,
                                                  positive_flat, negative_flat))
    {
      if (ittemp->second==0)
      {
        m_path.insert(ittemp, std::make_pair(dh, (negative_flat?-1:+1)));
        // insert the new dart before 'ittemp'
      }
      else
      {
        List_iterator ittemp=next_iterator(it);
        CGAL_assertion(ittemp->second!=0);
        ittemp->first=dh; // Move the first dart of the flat
        ittemp->second+=(ittemp->second>0?+1:-1); // Increment the length of the flat
      }
    }
    else
    {
      m_path.insert(ittemp, std::make_pair(dh, 0));
      // insert the new element before 'ittemp', thus after 'it'
    }
    ++m_length;
  }
  
  /// Right push the given l-shape.
  void right_push_l_shape(List_iterator& it)
  {
    CGAL_assertion(is_l_shape(it));

    if (m_path.size()==2 && it->second==0 &&
        next_iterator(it)->second==0)
    { // Special case of a path with only 2 darts
      it->first=m_map.template beta<2,1>(it->first);
      List_iterator ittemp=next_iterator(it);
      ittemp->first=m_map.template beta<2,0>(ittemp->first);
      return;
    }

    List_iterator it2=next_flat(it);  // Beginning of the second flat
    List_iterator it3=(is_end_of_flat(it2)?it2:next_iterator(it2)); // End of the second flat
     
    // 1) Add the first dart before flat 'it'
    add_dart_before(it, m_map.template beta<2,1>(it->first));
  
    // 2) Add the last dart after flat 'it3'
    add_dart_after(it3, m_map.template beta<2,0>(it3->first));

    // 3) Move the first flat
    if (it->second==0)
    { it=m_path.erase(it); }
    else
    {
      CGAL_assertion(it->second<0);
      it->first=m_map.template beta<2,1,1>(it->first);
      it->second=(-it->second)-1;
      if (it->second==0)
      { it=m_path.erase(next_iterator(it)); }
      else
      {
        List_iterator ittemp=next_iterator(it);
        ittemp->first=m_map.template beta<2,1,2,0>(ittemp->first);
      }
    }
   
    // 4) Move the second flat
    if (it2->second==0)
    {
      it=m_path.erase(it2);
    }
    else
    {
      CGAL_assertion(it2->second<0);
      it2->first=m_map.template beta<2,0,2,1>(it2->first);
      it2->second=(-it2->second)-1;
      if (it2->second==0)
      { it=m_path.erase(next_iterator(it2)); }
      else
      { it3->first=m_map.template beta<2,1,1>(it3->first); }
    }
    m_length-=2;
    if (it==m_path.end() && !is_empty())
    { it=m_path.begin(); }
  }

  /// Right push the path, if all all l-shape are pushed, otherwise only one.
  /// @return true iff the path was pushed
  bool right_push(bool all=true)
  {
    bool res=false;
    List_iterator it=m_path.begin();
    // unsigned int TOTO=0;
    while(it!=m_path.end())
    {
      if (is_l_shape(it))
      {
        // std::cout<<"right_push "<<TOTO++<<std::endl;
        right_push_l_shape(it); res=true;
        //CGAL_assertion(is_valid());
        //CGAL_assertion(is_valid_iterator(it));
        if (!all) { return true; }
      }
      else { move_to_next_l_shape(it); }
    }
    CGAL_assertion(is_valid());
    return res;
  }

  /// Canonize the path
  void canonize()
  {
     // ?? if (!is_closed()) { return; }

     bool modified=false;
     remove_spurs();

     do
     {
       do
       {
         modified=remove_brackets();
         modified=modified || remove_spurs();
       }
       while(modified);
       modified=right_push();
     }
     while(modified);

     CGAL_assertion(is_valid());
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
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      std::cout<<m_map.darts().index(it->first)<<"("<<it->second<<")";
      if (std::next(it)!=itend) { std::cout<<" "; }
    }
     if (is_closed())
     { std::cout<<" c "; } //<<m_map.darts().index(get_ith_dart(0)); }
     std::cout<<" length("<<length()<<")";
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

#ifdef CGAL_PWRLE_TURN_V2
  TDartIds m_darts_ids;
#endif //CGAL_PWRLE_TURN_V2
};

} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_WITH_RLE_H //
// EOF //
