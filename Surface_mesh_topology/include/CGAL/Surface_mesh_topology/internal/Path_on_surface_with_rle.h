// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_PATH_ON_SURFACE_WITH_RLE_H
#define CGAL_PATH_ON_SURFACE_WITH_RLE_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <list>
#include <utility>
#include <iostream>
#include <iterator>
#include <vector>
#include <set>
#include <CGAL/assertions.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

namespace CGAL {
namespace Surface_mesh_topology {

template<typename Map_>
class Path_on_surface;

namespace internal {

// A flat is a sequence of darts given by its two extremities: begin and end,
// with +2 turns (if length>0) or -2 turns (if length<0).
// length==0 => begin==end.
template<typename Map_>
class CFlat
{
  typedef Map_                            Map;
  typedef typename Map::Dart_const_handle Dart_const_handle;
  typedef CFlat<Map>                      Self;

public:
  CFlat(Dart_const_handle dh) : begin(dh), end(dh), length(0)
  {}

  CFlat(Dart_const_handle dh1, Dart_const_handle dh2, int l) :
    begin(dh1), end(dh2), length(l)
  {}

  bool operator==(const Self& other) const
  { return begin==other.begin && end==other.end && length==other.length; }

  bool operator!=(const Self& other) const
  { return !(operator==(other)); }

  friend std::ostream& operator<< (std::ostream& os, const Self& flat)
  {
    os<<"["<<&*(flat.begin)<<", "<<&*(flat.end)<<" (l="<<flat.length<<")]";
    return os;
  }

  Dart_const_handle begin, end;
  int length; // Length of the flat, positive flat if >0, negative flat if <0
};

template<typename Mesh>
class Path_on_surface_with_rle
{
  friend class Path_on_surface<Mesh>;

public:
  typedef Path_on_surface_with_rle<Mesh>     Self;
  typedef typename Get_map<Mesh, Mesh>::type Map;
  typedef typename Map::Dart_handle          Dart_handle;
  typedef typename Map::Dart_const_handle    Dart_const_handle;
  typedef CFlat<Map>                         Flat;
  typedef std::list<Flat>                    List_of_flats;
  typedef typename List_of_flats::iterator   List_iterator;
  // TODO typedef typename List_of_dart_length::const_iterator List_const_iterator;

  struct List_iterator_hash
  {
    std::size_t operator() (const List_iterator& lit) const
    {
      return boost::hash<const void*>()(&*lit);
    }
  };

  typedef boost::unordered_set<List_iterator, List_iterator_hash> Set_of_it;

#ifdef CGAL_PWRLE_TURN_V2
    typedef boost::unordered_map<Dart_const_handle, std::size_t> TDartIds;
#endif //CGAL_PWRLE_TURN_V2

  /// Constructor
  Path_on_surface_with_rle(const Mesh& amap
                         #ifdef CGAL_PWRLE_TURN_V2
                           , const TDartIds & darts_ids
                         #endif //CGAL_PWRLE_TURN_V2
                           )
    : m_map(amap),
      m_is_closed(false),
      m_length(0),
      m_use_only_positive(false),
      m_use_only_negative(false)
  #ifdef CGAL_PWRLE_TURN_V2
    , m_darts_ids(darts_ids)
  #endif //CGAL_PWRLE_TURN_V2
  {}

  /// Creates a Path_on_surface_with_rle from a Path_on_surface.
  /// If use_only_positive, consider only positive flats and not negative ones.
  /// If use_only_negative, consider only negative flats and not positive ones.
  /// If both are false, consider both positive and negative flats.
  /// Both cannot be true at the same time.
  /// Note that for a minimal surface of genus>=2, we cannot have both -2 and
  /// +2 as turn, and thus these parameters are useless.
  /// However, this case can occured for our unit tests on the cube, this is
  /// the reason of these parameters.
  Path_on_surface_with_rle(const Path_on_surface<Mesh>& apath,
                         #ifdef CGAL_PWRLE_TURN_V2
                           const TDartIds & darts_ids,
                         #endif //CGAL_PWRLE_TURN_V2
                           bool use_only_positive=false,
                           bool use_only_negative=false
                           ) :
    m_map(apath.get_map()),
    m_is_closed(apath.is_closed()),
    m_length(0),
#ifdef CGAL_PWRLE_TURN_V2
    m_darts_ids(darts_ids),
#endif //CGAL_PWRLE_TURN_V2
    m_use_only_positive(use_only_positive),
    m_use_only_negative(use_only_negative)
  {
    if (apath.is_empty()) return;

    std::size_t i=0, j=0, starti=0, length=0;
    bool positive_flat=false;
    bool negative_flat=false;
    
    if (apath.is_closed())
    {
      /// Tempo pour debug
      // i=get_default_random().get_int(0, apath.length());
      /// END Tempo pour debug
      if (!use_only_negative && apath.next_positive_turn(i)==2)
      { positive_flat=true; negative_flat=false; }
      else if (!use_only_positive && apath.next_negative_turn(i)==2)
      { positive_flat=false; negative_flat=true; }

      while ((positive_flat && apath.next_positive_turn(i)==2) ||
             (negative_flat && apath.next_negative_turn(i)==2))
      {
        i=apath.next_index(i);
        if (i==0) // Case of a closed path, made of only one flat part.
        {
          m_path.push_back(Flat(apath.front(), apath.back(),
                                (positive_flat?(apath.length()-1):
                                               -(apath.length()-1))));
          m_length=apath.length();
          CGAL_assertion(is_valid());
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
        push_back(apath[i], false);
        i=apath.next_index(i);
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
      std::swap(m_use_only_negative, p2.m_use_only_negative);
      std::swap(m_use_only_positive, p2.m_use_only_positive);
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
      m_use_only_negative=other.m_use_only_negative;
      m_use_only_positive=other.m_use_only_positive;
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

  /// @return true if this path is equal to other path. For closed paths, test
  ///         all possible starting darts.
  bool operator==(Self& other)
  {
    if (is_closed()!=other.is_closed() ||
        length()!=other.length())
      return false;

    if (is_empty() && other.is_empty())
    { return true; }

    // Note that we need to transform the Path_on_surface_with_rle into
    // the correspondings Path_on_surface, because a same path can have two
    // different rle representations.
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

  /// @return true if the given flat is valid
  bool is_flat_valid(const List_iterator& flat) const
  {
    CGAL_assertion(is_valid_iterator(flat));
    Dart_const_handle dhend=flat->begin;
    int nb=0;
    while(nb!=flat->length)
    {
      if (flat->length>0)
      { dhend=m_map.template beta<1, 2, 1>(dhend); ++nb; }
      else
      { dhend=m_map.template beta<2, 0, 2, 0, 2>(dhend); --nb; }
    }

    return dhend==flat->end;
  }

  /// Display the given flat
  void display_flat(std::ostream& os, const List_iterator& flat)
  {
    CGAL_assertion(is_valid_iterator(flat));
    os<<"["<<m_map.darts().index(flat->begin)<<", "
      <<m_map.darts().index(flat->end)<<" (l="<<flat->length<<")]";
  }

  /// @return true iff the path is valid; i.e. a sequence of flats two by
  ///              two adjacent.
  /// if test_minimal is true, test that there is no two consecutive flats
  /// that can be merged. If test_minimal is false, this test is not done.
  bool is_valid(bool test_minimal=true)
  {
    if (is_empty()) { return !is_closed(); } // an empty past is not closed

    unsigned int nbdarts=0,i=0;
    Dart_const_handle prev=(is_closed()?back():NULL); // Last dart of the path
    for (auto it=m_path.begin(); it!=m_path.end(); ++it)
    {
      if (prev!=NULL && !m_map.template belong_to_same_cell<0>
          (m_map.template beta<1>(prev), begin_of_flat(it)))
      {
        std::cerr<<"ERROR: The path is not valid: flat in position "<<i
                 <<" does not follow the previous dart."<<std::endl;
        return false;
      }

      if (!is_flat_valid(it))
      {
        display();
        std::cout<<"ERROR: The path is not valid: flat at position "<<i
                <<" with length "<<flat_length(it)
               <<" is not correct: its end does not"
              <<" correspond to the flat."<<std::endl;
        return false;
      }

      if (test_minimal && can_merge_with_next_flat(it))
      {
       /* std::cout<<"************* ERROR ************"<<std::endl;
        display();
        Path_on_surface<Map>(*this).display_pos_and_neg_turns(); std::cout<<std::endl;
        display_pos_and_neg_turns(); std::cout<<std::endl;
        std::cout<<"FLAT ERROR "<<*it<<std::endl; */

        std::cout<<"ERROR: The path is not valid: flat at position "<<i
                <<" with length "<<flat_length(it)<<" is not correct: it can be merged "
               <<"with the next flat with length "
              <<flat_length(next_iterator(it))<<" - turns between the two flats = +"
             <<next_positive_turn(it)<<" -"<<next_negative_turn(it)<<std::endl;
        return false;
      }
      prev=end_of_flat(it); // end of the previous flat
      nbdarts+=std::abs(flat_length(it))+1; // Number of darts of the flat
      ++i;
    }

    if (is_closed())
    {
      if (prev==NULL)
      {
        std::cout<<"ERROR: The path is not valid: it is empty and closed."
                <<std::endl;
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
  bool is_valid_iterator(const List_iterator& ittotest) const
  {
    if (ittotest==m_path.end()) { return false; }
    return true;
    // Assert too long; uncomment in case of bug.
    /* for (auto it=m_path.begin(); it!=m_path.end(); ++it)
    {
      if (it==ittotest) { return true; }
    }
    return false; */
  }

  /// @return the first dart of the path
  /// @pre !is_empty()
  Dart_const_handle front()
  {
    CGAL_assertion(!is_empty());
    return begin_of_flat(m_path.begin());
  }

  /// @return the last dart of the path
  /// @pre !is_empty()
  Dart_const_handle back()
  {
    CGAL_assertion(!is_empty());
    return end_of_flat(std::prev(m_path.end()));
  }

  /// @return true iff df can be added at the end of the path.
  bool can_be_pushed(Dart_const_handle dh)
  {
    // This assert is too long CGAL_assertion(m_map.darts().owns(dh));

    if (is_empty()) return true;

    return m_map.template belong_to_same_cell<0>
      (m_map.other_extremity(back()), dh);
  }

  /// Add the given dart at the end of this path.
  /// @pre can_be_pushed(dh)
  void push_back(Dart_const_handle dh, bool update_isclosed=true)
  {
    CGAL_assertion(dh!=Map::null_handle);
    // This assert is too long, it is tested in the is_valid method.
    CGAL_assertion(can_be_pushed(dh));

    List_iterator itlast=m_path.end();

    bool positive_flat, negative_flat;
    if (is_empty() ||
        !is_prev_flat_can_be_extended_at_end(itlast, dh,
                                             positive_flat, negative_flat))
    { m_path.push_back(Flat(dh)); } // Create a new flat
    else
    {
      --itlast;
      set_end_of_flat(itlast, dh); // Move the last dart of the last flat
      set_flat_length(itlast, flat_length(itlast)+(positive_flat?+1:-1)); // Increment the length of the flat
    }

    ++m_length;
    if (update_isclosed) { update_is_closed(); }
  }

  // Update m_is_closed to true iff the path is closed (i.e. the second
  //   extremity of the last dart of the path is the same vertex than the one
  //   of the first dart of the path).
  void update_is_closed()
  {
    if (is_empty()) { m_is_closed=false; }
    else
    {
      Dart_const_handle pend=m_map.other_extremity(back());
      if (pend==Map::null_handle) { m_is_closed=false; }
      else
      { m_is_closed=m_map.template belong_to_same_cell<0>(front(), pend); }
    }
  }

  /// @return true iff there is a flat after it
  bool next_flat_exist(const List_iterator& it) const
  {
    CGAL_assertion(it!=m_path.end());
    return is_closed() || std::next(it)!=m_path.end();
  }

  /// @return true iff there is a flat before it
  bool prev_flat_exist(const List_iterator& it) const
  {
    return is_closed() || it!=m_path.begin();
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
    if (it==m_path.end())
    {
      if (is_closed())
      { it=m_path.begin(); } // Here the path is closed, and it is the last element of the list
    }
  }

  void retreat_iterator(List_iterator& it)
  {
    if (it==m_path.begin())
    {
      if (is_closed())
      { it=m_path.end(); }
      else { return; }
    }
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

  /// @return the beginning of the flat
  Dart_const_handle begin_of_flat(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    return it->begin;
  }

  /// @return the end of the flat
  Dart_const_handle end_of_flat(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    return it->end;
  }

  /// @return the length of the flat
  int flat_length(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    return it->length;
  }

  void set_flat_length(const List_iterator& it, int i)
  {
    CGAL_assertion(is_valid_iterator(it));
    it->length=i;
  }

  void decrease_flat_length(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    if (flat_length(it)>0) { --(it->length); }
    else                   { ++(it->length); }
  }

  void increase_flat_length(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    if (flat_length(it)>0) { ++(it->length); }
    else                   { --(it->length); }
  }

  void set_begin_of_flat(const List_iterator& it, Dart_const_handle dh)
  {
    CGAL_assertion(is_valid_iterator(it));
    it->begin=dh;
  }

  void set_end_of_flat(const List_iterator& it, Dart_const_handle dh)
  {
    CGAL_assertion(is_valid_iterator(it));
    it->end=dh;
  }

  /// @return the second dart of the given flat.
  Dart_const_handle second_dart_of_a_flat(const List_iterator& it)
  {
    if (flat_length(it)>0)
    { return m_map.template beta<1, 2, 1>(begin_of_flat(it)); }

    CGAL_assertion(flat_length(it)<0);
    return m_map.template beta<2, 0, 2, 0, 2>(begin_of_flat(it));
  }

  /// @return the dart before the last dart of the given flat.
  Dart_const_handle before_last_dart_of_a_flat(const List_iterator& it)
  {
    if (flat_length(it)>0)
    { return m_map.template beta<0, 2, 0>(end_of_flat(it)); }

    CGAL_assertion(flat_length(it)<0);
    return m_map.template beta<2, 1, 2, 1, 2>(end_of_flat(it));
  }

  /// @return the positive turn given two ids of darts (unsed for CGAL_PWRLE_TURN_V1)
  std::size_t compute_positive_turn_given_ids(std::size_t id1,
                                              std::size_t id2) const
  {
    if (id1==id2) { return 0; }

    std::size_t number_of_edges=m_map.number_of_darts()/2;
    if (id1>=number_of_edges)
    {
      id1-=number_of_edges; // id of the first dart in its own vertex
      CGAL_assertion(id2>=number_of_edges);
      id2-=number_of_edges; // id of the second dart in its own vertex
    }

    return (id1<id2?(id2-id1):(number_of_edges-id1+id2));
  }

  /// @return the negative turn given two ids of darts (unused for CGAL_PWRLE_TURN_V1)
  std::size_t compute_negative_turn_given_ids(std::size_t id1,
                                              std::size_t id2) const
  {
    if (id1==id2) { return 0; }

    std::size_t number_of_edges=m_map.number_of_darts()/2;
    if (id1>=number_of_edges)
    {
      id1-=number_of_edges; // id of the first dart in its own vertex
      CGAL_assertion(id2>=number_of_edges);
      id2-=number_of_edges; // id of the second dart in its own vertex
    }

    return (id1<=id2?(number_of_edges-id2+id1):(id1-id2));
  }

  std::size_t get_dart_id(Dart_const_handle dh) const
  {
#ifdef CGAL_PWRLE_TURN_V3
    return m_map.info(dh);
#else // CGAL_PWRLE_TURN_V3
#ifdef CGAL_PWRLE_TURN_V2
    return m_darts_ids.at(dh);
#else //  CGAL_PWRLE_TURN_V2
    std::cerr<<"Error: impossible to get dart id without method V2 or V3."<<std::endl;
    return -1;
#endif // CGAL_PWRLE_TURN_V2
#endif // CGAL_PWRLE_TURN_V3
  }

  /// @return the turn between the two given darts.
  std::size_t positive_turn(Dart_const_handle dh1,
                            Dart_const_handle dh2)
  {
#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    unsigned int toto1=m_map.positive_turn(dh1, dh2);
    unsigned int toto2=compute_positive_turn_given_ids(get_dart_id(m_map.template beta<2>(dh1)),
                                                       get_dart_id(dh2));
   // std::cout<<"toto1="<<toto1<<"  toto2="<<toto2<<std::endl;
    CGAL_assertion(toto1==toto2);

    return compute_positive_turn_given_ids(get_dart_id(m_map.template beta<2>(dh1)),
                                           get_dart_id(dh2));
#else // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    return m_map.positive_turn(dh1, dh2);
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
  }

  /// @return the turn between the two given darts.
  std::size_t negative_turn(Dart_const_handle dh1,
                            Dart_const_handle dh2)
  {
#if defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    // TEMPO POUR DEBUG
    CGAL_assertion(m_map.negative_turn(dh1, dh2)==
                   compute_negative_turn_given_ids(get_dart_id(m_map.template beta<2>(dh1)),
                                                   get_dart_id(dh2)));

    return compute_negative_turn_given_ids(get_dart_id(m_map.template beta<2>(dh1)),
                                           get_dart_id(dh2));
#else // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
    return m_map.negative_turn(dh1, dh2);
#endif // defined(CGAL_PWRLE_TURN_V2) || defined(CGAL_PWRLE_TURN_V3)
  }

  /// @return the turn between dart it and next dart.
  ///         (turn is position of the second edge in the cyclic ordering of
  ///          edges starting from the first edge around the second extremity
  ///          of the first dart)
  std::size_t next_positive_turn(const List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    // CGAL_assertion(is_valid_iterator(it));
    CGAL_assertion(next_flat_exist(it));

    /* if (flat_length(it)!=0)
    {
      if (flat_length(it)<0) { return (m_map.number_of_darts()/2)-2; }
      else { return 2; }
    } */

    return positive_turn(end_of_flat(it), begin_of_flat(next_iterator(it)));
  }

  /// Same than next_positive_turn but turning in reverse orientation around vertex.
  std::size_t next_negative_turn(const List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    // CGAL_assertion(is_valid_iterator(it));
    CGAL_assertion(next_flat_exist(it));

    /*if (flat_length(it)!=0)
    {
      if (flat_length(it)<0) { return 2; }
      else { return (m_map.number_of_darts()/2)-2; }
    }*/

    return negative_turn(end_of_flat(it), begin_of_flat(next_iterator(it)));
  }

  std::vector<std::size_t> compute_positive_turns()
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      if (next_flat_exist(it))
      { res.push_back(next_positive_turn(it)); }
    }
    return res;
  }

  std::vector<std::size_t> compute_negative_turns()
  {
    std::vector<std::size_t> res;
    if (is_empty()) return res;
    for (auto it=m_path.begin(), itend=m_path.end(); it!=itend; ++it)
    {
      if (next_flat_exist(it))
      { res.push_back(next_negative_turn(it)); }
    }
    return res;
  }

  std::vector<std::size_t> compute_turns(bool positive)
  { return (positive?compute_positive_turns():compute_negative_turns()); }

  void flat_modified(const List_iterator& it,
                     Set_of_it& modified_flats)
  {
    modified_flats.insert(it);
    if (prev_flat_exist(it))
    { modified_flats.insert(prev_iterator(it)); }
    if (next_flat_exist(it))
    { modified_flats.insert(next_iterator(it)); }
  }

  void erase_flat(List_iterator& it, Set_of_it& modified_flats)
  {
    if (next_flat_exist(it))
    { modified_flats.insert(next_iterator(it)); }
    if (prev_flat_exist(it))
    { modified_flats.insert(prev_iterator(it)); }

    modified_flats.erase(it);

    it=m_path.erase(it); // after the erase, it is the element after the erased one
    if (prev_flat_exist(it))
    {
      retreat_iterator(it); // this is why we move it backward here
      // modified_flats.insert(it);
    }
  }

  List_iterator merge_modified_flats_when_possible(Set_of_it& modified_flats)
  {
    if (m_path.empty())
    {
      m_is_closed=false;
      return m_path.end();
    }

    List_iterator smallest_it=m_path.end();
    for (auto it=modified_flats.begin(), itend=modified_flats.end();
         it!=itend; ++it)
    {
      if (prev_flat_exist(*it) && modified_flats.count(prev_iterator(*it))==0)
      { smallest_it=*it; }
    }

    List_iterator itcur;
    if (smallest_it==m_path.end())
    { // Case of a closed path, with all paths modified
      for (auto it=modified_flats.begin(); it!=modified_flats.end(); ++it)
      {
        itcur=*it;
        while(merge_with_next_flat_if_possible(itcur, modified_flats))
        {}
      }
      return m_path.begin();
    }

    for (auto it=modified_flats.begin(); it!=modified_flats.end(); ++it)
    {
      itcur=*it;
      while(merge_with_next_flat_if_possible(itcur, modified_flats))
      {}
    }

    // Note that smallest is not necessarily the smallest, this is a flat
    // such that the previous flat is not modified.
    return smallest_it;
  }

  /// Reduce the length of the flat part starting at 'it' from its beginning
  /// 'it' moves to the previous flat if the current flat disapeared.
  /// The path could be not valid after this operation (consistency with next
  /// element should be ensure, by possibly updating the next flat part).
  /// @return true iff the flat disapeared after its reduction.
  bool reduce_flat_from_beginning(List_iterator& it,
                                  Set_of_it& modified_flats)
  {
    CGAL_assertion(is_valid_iterator(it));

    CGAL_assertion(m_length>0);
    --m_length;

    flat_modified(it, modified_flats);

    if (flat_length(it)==0)
    {
      erase_flat(it, modified_flats);
      return true;
    }

    //else
    set_begin_of_flat(it, second_dart_of_a_flat(it));
    decrease_flat_length(it);
    return false;
  }

  /// Reduce the length of the flat part starting at 'it' from its end.
  /// 'it' moves to the previous flat if the current flat disapeared.
  /// The path could be not valid after this operation (consistency with next
  /// element should be ensure, by possibly updating the next flat part).
  /// @return true iff the flat disapeared after its reduction.
  bool reduce_flat_from_end(List_iterator& it,
                            Set_of_it& modified_flats)
  {
    CGAL_assertion(is_valid_iterator(it));

    CGAL_assertion(m_length>0);
    --m_length;

    flat_modified(it, modified_flats);

    if (flat_length(it)==0)
    {
      erase_flat(it, modified_flats);
      return true;
    }

    // else
    set_end_of_flat(it, before_last_dart_of_a_flat(it));
    decrease_flat_length(it);
    return false;
  }

  /// @return true iff the flat 'it' can be merged with its next flat.
  bool can_merge_with_next_flat(const List_iterator& it,
                                bool& positive2, bool& negative2)
  {
    if (m_path.size()<=1) { return false; }

    /* if (!is_valid_iterator(it))
    {
      std::cout<<"ERROR"<<std::endl;
      CGAL_assertion(is_valid_iterator(it));
    } */

    if (!next_flat_exist(it)) { return false; }

    List_iterator it2=next_iterator(it);
    CGAL_assertion(m_path.size()>1 || it==it2);
    if (it==it2) { return false; } // only one flat in the path

    positive2=(!m_use_only_negative && (next_positive_turn(it)==2));
    negative2=(!m_use_only_positive && (next_negative_turn(it)==2));
    CGAL_assertion(!(positive2 && negative2));

    if (!positive2 && !negative2) { return false; } // the middle turn is not a flat

    bool positive1=false, negative1=false;
    bool positive3=false, negative3=false;

    if (flat_length(it)>0) positive1=true;
    else if (flat_length(it)<0) negative1=true;

    if (flat_length(it2)>0) positive3=true;
    else if (flat_length(it2)<0) negative3=true;

    if (!positive1 && !negative1)
    { // First flat is empty (length 0)
      if (!positive3 && !negative3)
      { // and second flat too
        return true;
      }
      else
      { // here second flat is not empty
        if ((positive3 && !positive2) || (negative3 && !negative2))
        { return false; } // the two flats cannot be merged
      }
    }
    else
    { // here first flat is not empty
      if (!positive3 && !negative3)
      { // and second flat is empty
        if ((positive1 && !positive2) || (negative1 && !negative2))
        { return false; } // the two flats cannot be merged
      }
      else
      {
        // Second flat is not empty
        if ((positive1!=positive3) || (negative1!=negative3) ||
            (positive1 && !positive2) || (negative1 && !negative2))
        { return false; } // the two flats cannot be merged
      }
    }
    return true;
  }

  /// @return true iff the flat 'it' can be merged with its next flat.
  bool can_merge_with_next_flat(const List_iterator& it)
  {
    bool dummy1, dummy2;
    return can_merge_with_next_flat(it, dummy1, dummy2);
  }

  /// Merge flat 'it' with its next flat if it is possible.
  /// @return true if a merging was done.
  bool merge_with_next_flat_if_possible(const List_iterator& it,
                                        Set_of_it& modified_flats)
  {
    bool positive2=false, negative2=false;
    if (!can_merge_with_next_flat(it, positive2, negative2))
    { return false; }

#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"Merge with next flat: ";
    display_flat(std::cout, it); std::cout<<" and ";
    display_flat(std::cout, next_iterator(it)); std::cout<<std::endl;
#endif

    //Path_on_surface<Map> pbefore(*this);

    List_iterator it2=next_iterator(it);
    set_flat_length(it, flat_length(it)+flat_length(it2)+
                    (positive2?1:-1));
    set_end_of_flat(it, end_of_flat(it2));
    m_path.erase(it2);
    if (modified_flats.count(it2)==1)
      modified_flats.erase(it2);

    //Path_on_surface<Map> pafter(*this);
    //CGAL_assertion(pbefore==pafter);

    // CGAL_assertion(is_flat_valid(it));
    return true;
  }

  /// Merge flat 'it' with its next flat if it is possible.
  /// @return true if a merging was done.
  bool merge_with_next_flat_if_possible(const List_iterator& it)
  {
    Set_of_it dummy;
    return merge_with_next_flat_if_possible(it, dummy);
  }

  /// Merge the last flat of this path with its next flat if it is possible.
  /// @return true if a merging was done.
  bool merge_last_flat_with_next_if_possible()
  {
    List_iterator lastit=m_path.end(); --lastit;
    return merge_with_next_flat_if_possible(lastit);
  }

  /// Return true if it is the beginning of a spur.
  bool is_spur(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    return next_flat_exist(it) &&
        m_map.template beta<2>(end_of_flat(it))==
        begin_of_flat(next_iterator(it));
  }

  /// Remove the spur given by 'it'.
  /// After the operation, either 'it' stay on the same element (if the flat
  /// still exist), or 'it' move to the previous element if the flat containing
  /// the spur is removed.
  /// ('it' will be equal to m_path.end() if the path becomes empty).
  void remove_spur(List_iterator& it)
  {
#ifdef CGAL_TRACE_PATH_TESTS
        std::cout<<"Remove spur between flats: ";
        display_flat(std::cout, it); std::cout<<" and ";
        display_flat(std::cout, next_iterator(it)); std::cout<<std::endl;
#endif

    CGAL_assertion(is_spur(it));
    Set_of_it modified_flats;

    //std::cout<<"*********************************** ";
   /* static int TOTO=0;
    std::cout<<"remove_spur "<<TOTO++<<std::endl; */

    /*if (TOTO==984)
    {
      int dh=TOTO;
     CGAL_assertion(is_valid());
      display();
     // Path_on_surface<Map>(*this).display_pos_and_neg_turns(); std::cout<<std::endl;
      display_pos_and_neg_turns(); std::cout<<std::endl;
    }*/

    List_iterator it2=next_iterator(it); // it2 is the next flat

    // We reduce the first flat
    reduce_flat_from_end(it, modified_flats); // If flat 'it' is erased, it is moved to the previous flat

    // And we reduce the second flat
    reduce_flat_from_beginning(it2, modified_flats);

    /* std::cout<<*prev_iterator(it)<<std::endl
            <<*it<<std::endl<<*it2<<std::endl<<*next_iterator(it2)<<std::endl;*/

    it=merge_modified_flats_when_possible(modified_flats);

    // CGAL_assertion(is_valid_iterator(it));
    // CGAL_assertion(is_valid());
  }

  /// Move 'it' to the next spur after it. Go to m_path.end() if there is no
  /// spur in the path.
  void move_to_next_spur(List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    List_iterator itend=(is_closed()?it:m_path.end());
    do
    {
      advance_iterator(it);
      if (it!=m_path.end() && is_spur(it)) { return; }
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
        // CGAL_assertion(is_valid_iterator(it));
        // CGAL_assertion(is_valid());
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
    CGAL_assertion(is_valid_iterator(it));
    if (m_use_only_negative || !next_flat_exist(it) ||
        next_positive_turn(it)!=1)
    { return false; }
    // Here first positive turn is 1

    itend=next_iterator(it);
    // Here itend is the beginning of the second flat
    if (flat_length(itend)<0)
    { return false; } // This is not a positive bracket, we have some -2

    if (!next_flat_exist(itend) || next_positive_turn(itend)!=1)
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
    CGAL_assertion(is_valid_iterator(it));
    if (m_use_only_positive || !next_flat_exist(it) ||
        next_negative_turn(it)!=1)
    { return false; }
    // Here first negative turn is 1

    itend=next_iterator(it);
    // Here itend is the beginning of the second flat
    if (flat_length(itend)>0)
    { return false; } // This is not a negative bracket, we have some +2

    if (!next_flat_exist(itend) || next_negative_turn(itend)!=1)
    { return false; }

    // Now we move itend to the beginning of the third flat.
    advance_iterator(itend);
    return true;
  }

  /// Move 'it' to the next bracket after 'it'. Go to m_path.end() if there
  /// is no bracket in the path.
  void move_to_next_bracket(List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    List_iterator itend=(is_closed()?it:m_path.end());
    List_iterator it2;
    do
    {
      advance_iterator(it);
      if (it!=m_path.end() &&
          (is_positive_bracket(it, it2) || is_negative_bracket(it, it2)))
      { return; }
    }
    while(it!=itend);
    it=m_path.end(); // Here there is no spur in the whole path
  }

  /// Remove the given negative bracket. it1 is the flat beginning
  /// of the bracket; it3 is the flat end of the bracket.
  void remove_negative_bracket(List_iterator& it1, List_iterator& it3)
  {
#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"Remove negative bracket: ";
    display_flat(std::cout, it1); std::cout<<", ";
    display_flat(std::cout, next_iterator(it1)); std::cout<<" and ";
    display_flat(std::cout, it3); std::cout<<std::endl;
#endif

    Set_of_it modified_flats;
    List_iterator it2;
    CGAL_assertion(is_negative_bracket(it1, it2)); // Here it2 is unused

    if (size_of_list()==1) // it1==it3)
    { // Case of cyclic bracket
      CGAL_assertion(size_of_list()==1);
      CGAL_assertion(flat_length(it1)<0);
      set_begin_of_flat(it1, m_map.template  beta<2,0,2,1>(begin_of_flat(it1)));
      set_end_of_flat(it1, m_map.template  beta<2,1,2,0>(end_of_flat(it1)));
      set_flat_length(it1, (-flat_length(it1))-2);
      flat_modified(it1, modified_flats);
      it1=merge_modified_flats_when_possible(modified_flats);
      CGAL_assertion(m_length>1);
      m_length-=2;
      return;
    }
    
    it2=next_iterator(it1); // it2 is the the next flat after it1
    
    reduce_flat_from_end(it1, modified_flats); // decrease also m_length
    reduce_flat_from_beginning(it3, modified_flats);

    it1=it2; // Beginning of the second flat, we are sure it still exists

    set_begin_of_flat(it2, m_map.template  beta<2,1,1>(begin_of_flat(it2)));
    set_end_of_flat(it2, m_map.template  beta<2,1,1>(end_of_flat(it2)));
    set_flat_length(it2, -flat_length(it2));

    flat_modified(it2, modified_flats);

    it1=merge_modified_flats_when_possible(modified_flats);

    // CGAL_assertion(is_valid());
  }

  /// Remove the given positive bracket. 'it1' is the flat beginning
  /// of the bracket; 'it3' is the flat end of the bracket.
  void remove_positive_bracket(List_iterator& it1,
                               List_iterator& it3)
  {
#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"Remove positive bracket: ";
    display_flat(std::cout, it1); std::cout<<", ";
    display_flat(std::cout, next_iterator(it1)); std::cout<<" and ";
    display_flat(std::cout, it3); std::cout<<std::endl;
#endif

    Set_of_it modified_flats;
    List_iterator it2;
    CGAL_assertion(is_positive_bracket(it1, it2)); // Here it2 is unused

    if (size_of_list()==1) // it1==it3)
    { // Case of cyclic bracket
      CGAL_assertion(size_of_list()==1);
      CGAL_assertion(flat_length(it1)>0);
      set_begin_of_flat(it1, m_map.template  beta<1,2,0,2>(begin_of_flat(it1)));
      set_end_of_flat(it1, m_map.template  beta<0,2,1,2>(end_of_flat(it1)));
      set_flat_length(it1, -(flat_length(it1)-2));
      flat_modified(it1, modified_flats);
      it1=merge_modified_flats_when_possible(modified_flats);
      CGAL_assertion(m_length>1);
      m_length-=2;
      return;
    }

    it2=next_iterator(it1); // it2 is the next flat after it1

    reduce_flat_from_end(it1, modified_flats); // decrease also m_length
    reduce_flat_from_beginning(it3, modified_flats);

    it1=it2; // Beginning of the second flat, we are sure it exists

    set_begin_of_flat(it2, m_map.template  beta<1,1,2>(begin_of_flat(it2)));
    set_end_of_flat(it2, m_map.template  beta<1,1,2>(end_of_flat(it2)));
    set_flat_length(it2, -flat_length(it2));

    flat_modified(it2, modified_flats);

    it1=merge_modified_flats_when_possible(modified_flats);

    // CGAL_assertion(is_valid());
  }
  
  /// Simplify the path by removing all brackets, if all==true (default),
  /// or by removing only one bracket, if all==false.
  /// @return true iff at least one bracket was removed
  bool remove_brackets(bool all=true)
  {
    // CGAL_assertion(is_valid());
    bool res=false;
    List_iterator it1=m_path.begin();
    List_iterator it2;
    while(it1!=m_path.end())
    {
      // CGAL_assertion(is_valid());
      if (is_positive_bracket(it1, it2))
      {
        remove_positive_bracket(it1, it2); res=true;
        // CGAL_assertion(is_valid_iterator(it1));
        // CGAL_assertion(is_valid());
      }
      else if (is_negative_bracket(it1, it2))
      {
        remove_negative_bracket(it1, it2); res=true;
        // CGAL_assertion(is_valid_iterator(it1));
        // CGAL_assertion(is_valid());
      }
      else { move_to_next_bracket(it1); }
      if (!all && res) { return true; }
    }
    // CGAL_assertion(is_valid());
    return res;
  }

  /// @return true iff 'it' is the beginning of a l-shape.
  bool is_l_shape(const List_iterator& it)
  {
    CGAL_assertion(is_valid_iterator(it));
    if (!next_flat_exist(it)) { return false; }
    std::size_t t=next_negative_turn(it);
    return (t==1 ||
            (flat_length(it)<0 && size_of_list()==1 && t==2));
  }

  /// Move it to the next l_shape after it. Go to m_path.end() if there is no
  /// l_shape in the path.
  void move_to_next_l_shape(List_iterator& it)
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(is_valid_iterator(it));
    List_iterator itend=(is_closed()?it:m_path.end());
    do
    {
      advance_iterator(it);
      if (it!=m_path.end() && is_l_shape(it)) { return; }
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
    positive_flat=false; negative_flat=false;
    if (!prev_flat_exist(it)) { return false; }
    List_iterator ittemp=prev_iterator(it);
    if (!m_use_only_negative && positive_turn(end_of_flat(ittemp), dh)==2)
    { positive_flat=true; }
    if (!m_use_only_positive && negative_turn(end_of_flat(ittemp), dh)==2)
    { negative_flat=true; }

    if (flat_length(ittemp)==0) // Case of flat lengh 0
    { return positive_flat || negative_flat; }

    return (flat_length(ittemp)>0 && positive_flat) ||
      (flat_length(ittemp)<0 && negative_flat);
  }

  /// @return true iff the flat before flat 'it' can be extended by adding
  ///              dart 'dh' to its end.
  bool is_prev_flat_can_be_extended_at_end(const List_iterator& it,
                                           Dart_const_handle dh)
  {
    bool dummy1, dummy2;
    return is_prev_flat_can_be_extended_at_end(it, dh, dummy1, dummy2);
  }

  /// @return true iff the flat after flat 'it' can be extended by adding
  ///              dart 'dh' to its beginning.
  bool is_next_flat_can_be_extended_at_beginning(const List_iterator& it,
                                                 Dart_const_handle dh,
                                                 bool& positive_flat,
                                                 bool& negative_flat)
  {
     CGAL_assertion(is_valid_iterator(it));

     positive_flat=false; negative_flat=false;
     if (!next_flat_exist(it)) { return false; }
     List_iterator ittemp=next_iterator(it);
     if (!m_use_only_negative && positive_turn(dh, begin_of_flat(ittemp))==2)
     { positive_flat=true; }
     if (!m_use_only_positive && negative_turn(dh, begin_of_flat(ittemp))==2)
     { negative_flat=true; }

     if (flat_length(ittemp)==0)  // Case of flat lengh 0
     { return positive_flat || negative_flat; }

     return (flat_length(ittemp)>0 && positive_flat) ||
       (flat_length(ittemp)<0 && negative_flat);
  }
  
  /// @return true iff the flat after flat 'it' can be extended by adding
  ///              dart 'dh' to its beginning.
  bool is_next_flat_can_be_extended_at_beginning(const List_iterator& it,
                                                 Dart_const_handle dh)
  {
    bool dummy1, dummy2;
    return is_next_flat_can_be_extended_at_beginning(it, dh, dummy1, dummy2);
  }

  /// Add the given dart 'dh' before the flat 'it'.
  void add_dart_before(const List_iterator& it, Dart_const_handle dh,
                       Set_of_it& modified_flats)
  {
    CGAL_assertion(is_valid_iterator(it));

    bool positive_flat, negative_flat;
    if (is_prev_flat_can_be_extended_at_end(it, dh,
                                            positive_flat, negative_flat))
    {
      List_iterator ittemp=prev_iterator(it);
      set_end_of_flat(ittemp, dh); // Move the last dart of the previous flat
      set_flat_length(ittemp, flat_length(ittemp)+(positive_flat?+1:-1)); // Increment the length of the flat
      flat_modified(ittemp, modified_flats);
    }
    else
    {
      // insert the new element before 'it'
      flat_modified(m_path.insert(it, Flat(dh)), modified_flats);
    }
    ++m_length;
  }
  
  /// Add the given dart 'dh' after the flat 'it'.
  void add_dart_after(const List_iterator& it, Dart_const_handle dh,
                      Set_of_it& modified_flats)
  {
    CGAL_assertion(is_valid_iterator(it));
    bool positive_flat, negative_flat;
    List_iterator ittemp=next_iterator(it);
    if (is_next_flat_can_be_extended_at_beginning(it, dh,
                                                  positive_flat, negative_flat))
    {
      set_begin_of_flat(ittemp, dh); // Move the first dart of the flat
      set_flat_length(ittemp, flat_length(ittemp)+(positive_flat?+1:-1)); // Increment the length of the flat
      flat_modified(ittemp, modified_flats);
    }
    else
    {
      // insert the new element before 'ittemp', thus after 'it'
      flat_modified(m_path.insert(ittemp, Flat(dh)), modified_flats);
    }
    ++m_length;
  }
  
  /// Right push the given l-shape.
  void right_push_l_shape(List_iterator& it1)
  {
#ifdef CGAL_TRACE_PATH_TESTS
        std::cout<<"Right push l-shape: ";
        display_flat(std::cout, it1); std::cout<<" and ";
        display_flat(std::cout, next_iterator(it1)); std::cout<<std::endl;
#endif

    // it is the beginning of a flat
    CGAL_assertion(is_l_shape(it1));
    // CGAL_assertion(is_valid());
    Set_of_it modified_flats;

    /* static int TOTO=0;
    std::cout<<"right_push_l_shape "<<TOTO++<<std::endl;
    if (TOTO==1)
    { std::cout<<"BUG"<<std::endl;} */

    //Path_on_surface<Map> pbefore(*this);

    if (m_path.size()==1)
    {
      if (next_negative_turn(it1)==2)
      { // Special case of a global shift
        set_flat_length(it1, -flat_length(it1));
        set_begin_of_flat(it1, m_map.template beta<2,1,1>(begin_of_flat(it1)));
        set_end_of_flat(it1, m_map.template beta<2,1,1>(end_of_flat(it1)));
        // CGAL_assertion(is_valid());
        return;
      }
      else // Here negative turn is 1
      {
        if (flat_length(it1)>0) // Case where the first flat is positive
        { // We split the flat in two parts
          CGAL_assertion(flat_length(it1)>=1);
          Dart_const_handle dh1=begin_of_flat(it1);
          Dart_const_handle dh2=end_of_flat(it1);
          if (flat_length(it1)==1) // only one flat with two darts
          {
            reduce_flat_from_end(it1, modified_flats);
            it1=m_path.insert(it1, Flat(dh2)); // insert dh1 before it1
            flat_modified(it1, modified_flats);
            ++m_length;
          }
          else
          {
            reduce_flat_from_beginning(it1, modified_flats);
            reduce_flat_from_end(it1, modified_flats);
            flat_modified(m_path.insert(it1, Flat(dh1)), modified_flats); // insert dh1 before it1
            it1=m_path.insert(next_iterator(it1), Flat(dh2)); // insert dh2 after it1
            flat_modified(it1, modified_flats);
            m_length+=2;
          }
          // Now we can continue with the normal case because we have 3 flats
          CGAL_assertion(flat_length(it1)==0);
          CGAL_assertion(flat_length(next_iterator(it1))==0);
        }
        else
        {
          // Here we have only one flat, with some -2
          // Is this case possible ??
          CGAL_assertion(false);
          return;
        }
      }
    }

    List_iterator it2=next_iterator(it1);  // Beginning of the second flat
    CGAL_assertion(it1!=it2);

    if (flat_length(it1)>0) // Case where the first flat is positive
    { // We split the flat in two parts
      Dart_const_handle dh=end_of_flat(it1); // last dart of the first flat
      reduce_flat_from_end(it1, modified_flats);
      it1=m_path.insert(it2, Flat(dh)); // insert dh before it2
      ++m_length;
      flat_modified(it1, modified_flats);
    }

    if (flat_length(it2)>0) // Case where the second flat is positive, we need to
    { // split it into two parts
      Dart_const_handle dh=begin_of_flat(it2); // first dart of the second flat
      reduce_flat_from_beginning(it2, modified_flats);
      it2=m_path.insert(it2, Flat(dh)); // insert dh before the second flat
      ++m_length;
      flat_modified(it2, modified_flats);
    }

    //Path_on_surface<Map> pafter(*this);
    //CGAL_assertion(pbefore==pafter);

    // CGAL_assertion(is_valid(false));
    // CGAL_assertion(is_l_shape(it1));

    Dart_const_handle dh1=m_map.template beta<2,1>(begin_of_flat(it1));
    Dart_const_handle dh2=m_map.template beta<1>(dh1);
    Dart_const_handle dh3=m_map.template beta<2,1,2,0>(end_of_flat(it1));
    Dart_const_handle dh4=m_map.template beta<2,0,2,1>(begin_of_flat(it2));
    Dart_const_handle dh5=m_map.template beta<2,0,0>(end_of_flat(it2));
    Dart_const_handle dh6=m_map.template beta<1>(dh5);

    bool first_flat_zero=(flat_length(it1)==0);
    bool second_flat_zero=(flat_length(it2)==0);

    if (first_flat_zero)
    {
      if (second_flat_zero)
      { // Special case of a l-shape with only 2 darts
        set_begin_of_flat(it1, dh1);
        set_end_of_flat(it1, dh1);
        set_begin_of_flat(it2, dh6);
        set_end_of_flat(it2, dh6);

        flat_modified(it1, modified_flats);
        flat_modified(it2, modified_flats);

        it1=merge_modified_flats_when_possible(modified_flats);

        // CGAL_assertion(is_valid());
        return;
      }
      else
      { // Here first flat length is 0, while second flat not
        set_begin_of_flat(it1, dh1);
        set_end_of_flat(it1, dh5);
        set_begin_of_flat(it2, dh6);
        set_end_of_flat(it2, dh6);
        set_flat_length(it1, -flat_length(it2));
        set_flat_length(it2, 0);

        flat_modified(it1, modified_flats);
        flat_modified(it2, modified_flats);

        it1=merge_modified_flats_when_possible(modified_flats);

        // CGAL_assertion(is_valid());
        return;
      }
    }
    else
    {
      if (second_flat_zero)
      { // Here first flat length is non zero, while second flat length is 0
        set_begin_of_flat(it1, dh1);
        set_end_of_flat(it1, dh1);
        set_begin_of_flat(it2, dh2);
        set_end_of_flat(it2, dh6); // TODO WHY THIS DOES NOT WORK ????
        //set_end_of_flat(it2, m_map.template beta<1, 2, 1>(dh3)); // TODO VERY STRANGE
        set_flat_length(it2, -flat_length(it1));
        set_flat_length(it1, 0);

        flat_modified(it1, modified_flats);
        flat_modified(it2, modified_flats);

        it1=merge_modified_flats_when_possible(modified_flats);

        // CGAL_assertion(is_valid());
        return;
      }
    }

    // General case, with two flats with non zero length

/*    List_iterator it_insert1=m_path.insert(it, std::make_pair
                                           (first_dart_toinsert, 0));
    List_iterator it_insert2=m_path.insert(next_iterator(it3),
                                           std::make_pair
                                           (second_dart_toinsert, 0)); */

    if (next_iterator(it2)!=it1 || next_negative_turn(it2)!=3)
    {
      // 1) Add the first dart before flat 'it'
      add_dart_before(it1, dh1, modified_flats);

      // 2) Add the last dart after flat 'it2'
      add_dart_after(it2, dh6, modified_flats);
    }
    else
    { // Case "-2^s -1 -2^t -3" is a special case
      dh2=m_map.template beta<2,0>(dh1);
      dh5=m_map.template beta<2,1>(dh6);
      increase_flat_length(it1);
      increase_flat_length(it2);
      m_length+=2;
      flat_modified(it1, modified_flats);
      flat_modified(it2, modified_flats);
    }

    // 3) Move the first flat
    CGAL_assertion(flat_length(it1)<0);
    set_begin_of_flat(it1, dh2);
    set_flat_length(it1, -(flat_length(it1))-1);
    if (flat_length(it1)==0) { set_end_of_flat(it1, dh2); }
    else { set_end_of_flat(it1, dh3); } // End of the moved flat
    flat_modified(it1, modified_flats);
   
    // 4) Move the second flat
    CGAL_assertion(flat_length(it2)<0);
    set_begin_of_flat(it2, dh4);
    set_flat_length(it2, -(flat_length(it2))-1);
    if (flat_length(it2)==0) { set_end_of_flat(it2, dh4); }
    else { set_end_of_flat(it2, dh5); } // End of the moved flat
    flat_modified(it2, modified_flats);

    CGAL_assertion(m_length>1);
    m_length-=2;

    it1=merge_modified_flats_when_possible(modified_flats);
  }

  /// Right push the path, if all all l-shape are pushed, otherwise only one.
  /// @return true iff the path was pushed
  bool right_push(bool all=true)
  {
    bool res=false;
    List_iterator it=m_path.begin();
    while(it!=m_path.end())
    {
      if (is_l_shape(it))
      {
        right_push_l_shape(it); res=true;
        // CGAL_assertion(is_valid_iterator(it));
        // CGAL_assertion(is_valid());
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
    CGAL_assertion(is_valid());

    if (is_empty()) { return; }

    CGAL_assertion(is_closed());

     bool modified=false;
     remove_spurs();

     do
     {
       modified=remove_brackets();
       modified=modified || remove_spurs();
     }
     while(modified);

     right_push();

     CGAL_assertion(remove_brackets()==false);
     CGAL_assertion(remove_spurs()==false);
     CGAL_assertion(is_valid());
  }

  void display_positive_turns()
  {
    std::cout<<"+(";
    for (List_iterator it=m_path.begin(), itend=m_path.end();
         it!=itend; ++it)
    { std::cout<<next_positive_turn(it)<<" "; }
    std::cout<<")";
  }

  void display_negative_turns()
  {
    std::cout<<"-(";
    for (List_iterator it=m_path.begin(), itend=m_path.end();
         it!=itend; ++it)
    { std::cout<<next_negative_turn(it)<<" "; }
   std::cout<<")";
  }

  void display_pos_and_neg_turns()
  {
    display_positive_turns();
    std::cout<<"  ";
    display_negative_turns();
  }

  void display()
  {
    if (!is_empty())
    {
      for (List_iterator it=m_path.begin(), itend=m_path.end();
           it!=itend; ++it)
      {
        std::cout<<"[ ";
        std::cout<<m_map.darts().index(begin_of_flat(it))
                 <<", "<<m_map.darts().index(end_of_flat(it))
                 <<"("<<flat_length(it)<<") ] ";
      }
      if (is_closed())
      { std::cout<<" c "; }
      std::cout<<" length("<<length()<<")";
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Self& p)
  {
    p.display();
    return os;
  }

protected:
  const typename Get_map<Mesh, Mesh>::storage_type m_map; // The underlying map
  List_of_flats m_path; // The sequence of flats (a flat part is a pair of darts with positive or negative turn == 2). If negative value k, -k is the length of the flat part, for negative turns (-2).
  bool m_is_closed; // True iff the path is a cycle
  std::size_t m_length;
  bool m_use_only_positive;
  bool m_use_only_negative;

#ifdef CGAL_PWRLE_TURN_V2
  const TDartIds& m_darts_ids;
#endif //CGAL_PWRLE_TURN_V2
};

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_WITH_RLE_H //
// EOF //
