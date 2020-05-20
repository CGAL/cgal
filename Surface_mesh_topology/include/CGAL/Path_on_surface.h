// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_PATH_ON_SURFACE_H
#define CGAL_PATH_ON_SURFACE_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Random.h>
#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h>
#include <boost/algorithm/searching/knuth_morris_pratt.hpp>
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <initializer_list>

// A Path_on_surface contains two vectors of equal length n
// The first one is a vector of darts called m_path and the second one a vector
// of booleans called m_flip.
// If n = 0, the path represented by those vectors is the empty path.
// Else, it is the path represented by the n-1 first elements of both vectors,
// at the one we add the m_path[n-1] dart if m_flip[n-1] is false and the
// opposite of this dart if m_flip[n-1] is true i.e. if m_flip[i] is true means
// that the i-th dart m_path[i] has to be flipped.
// We use flips because sometimes opposite darts doesn't exist on surfaces with
// boundaries. But if m_flip[i] is true doesn't necesary mean that
// m_path[i] is 2-free

namespace CGAL {
namespace Surface_mesh_topology {

template<typename Mesh_>
class Path_on_surface
{
public:
  typedef Path_on_surface<Mesh_>             Self;
  typedef Mesh_                              Mesh;
  typedef typename Get_map<Mesh, Mesh>::type Map; // Mesh seen as a 2-map
  typedef typename Map::Dart_const_handle    Dart_const_handle;

  typedef Dart_const_handle halfedge_descriptor; // To be compatible with BGL

  Path_on_surface(const Mesh& amesh) : m_map(amesh), m_is_closed(false)
  {}

  template<class COST>
  Path_on_surface(const internal::Path_on_surface_with_rle<COST>& apath) :
    m_map(apath.get_map()),
    m_is_closed(apath.is_closed())
  {
    for (auto it=apath.m_path.begin(), itend=apath.m_path.end(); it!=itend; ++it)
    {
      push_back(it->begin, false, false);
      if (it->length>0)
      { extend_straight_positive(it->length, false); }
      else if (it->length<0)
      { extend_straight_negative(-(it->length), false); }
    }
    update_is_closed();
    CGAL_assertion(is_valid(true));
  }

  Path_on_surface(const Self& apath) : m_map(apath.m_map),
                                       m_path(apath.m_path),
                                       m_is_closed(apath.m_is_closed),
                                       m_flip(apath.m_flip)
  {}

  void swap(Self& p2)
  {
    if (this==&p2) { return; }

    CGAL_assertion(&get_mesh()==&(p2.get_mesh()));
    m_path.swap(p2.m_path);
    std::swap(m_is_closed, p2.m_is_closed);
    m_flip.swap(p2.m_flip);
  }

  Self& operator=(const Self& other)
  {
    CGAL_assertion(&get_mesh()==&(other.get_mesh()));
    if (this!=&other)
    {
      m_path=other.m_path;
      m_is_closed=other.m_is_closed;
      m_flip=other.m_flip;
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

  /// @return the combinatorial map supporting this path.
  const Mesh& get_mesh() const
  { return Get_map<Mesh, Mesh>::get_mesh(m_map); }

  const std::vector<bool>& get_flip() const
  { return m_flip; }

  /// clear the path.
  void clear()
  {
    m_path.clear();
    m_flip.clear();
    m_is_closed=false;
  }

  /// @return true iff the prev index exists
  bool prev_index_exists(std::size_t i) const
  { return is_closed() || i>0; }

  /// @return true iff the next index exists
  bool next_index_exists(std::size_t i) const
  { return is_closed() || i<(m_path.size()-1); }

  /// @return the index after index i.
  std::size_t next_index(std::size_t i) const
  { return ((is_closed() && i==(m_path.size()-1))?0:(i+1)); }

  /// @return the index before index i.
  std::size_t prev_index(std::size_t i) const
  { return ((is_closed() && i==0)?(m_path.size()-1):(i-1)); }

  /// @return the ith dart of the path.
  Dart_const_handle get_ith_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return m_path[i];
  }

  /// @return true iff the ith dart is flipped
  bool get_ith_flip(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return m_flip[i];
  }

  /// @return the ith dart of the path.
  Dart_const_handle operator[] (std::size_t i) const
  { return get_ith_dart(i); }

  /// @return the dart before the ith dart of the path,
  ///          Map::null_handle if such a dart does not exist.
  Dart_const_handle get_prev_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==0 && !is_closed()) return Map::null_handle;
    return m_path[prev_index(i)];
  }

  /// @return the dart after the ith dart of the path,
  ///          Map::null_handle if such a dart does not exist.
  Dart_const_handle get_next_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==m_path.size()-1 && !is_closed()) return Map::null_handle;
    return m_path[next_index(i)];
  }

  /// @return the flip before the ith flip of the path,
  ///          false if such a flip does not exist.
  bool get_prev_flip(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==0 && !is_closed()) return false;
    return m_flip[prev_index(i)];
  }

  /// @return the flip after the ith flip of the path,
  ///          false if such a flip does not exist.
  bool get_next_flip(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    if (i==m_path.size()-1 && !is_closed()) return false;
    return m_flip[next_index(i)];
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

  /// @return the first flip of the path.
  /// @pre !is_empty()
  bool front_flip() const
  {
    CGAL_assertion(!is_empty());
    return m_flip.front();
  }

  /// @return the last flip of the path.
  /// @pre !is_empty()
  bool back_flip() const
  {
    CGAL_assertion(!is_empty());
    return m_flip.back();
  }

  /// @return the index of the first dart of the path.
  /// @pre !is_empty()
  std::size_t front_index() const
  { return get_map().darts().index(front()); }

  /// @return the index of the last dart of the path.
  /// @pre !is_empty()
  std::size_t back_index() const
  { return get_map().darts().index(back()); }

  /// @return the ith dart of the path taking into account flip.
  /// return null_handle if flip and there is no beta2
  Dart_const_handle get_ith_real_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return (get_ith_flip(i)?get_map().opposite2(get_ith_dart(i)):
                            get_ith_dart(i));
  }

  /// @return the opposite of the ith dart of the path taking into account flip.
  /// return null_handle if !flip and there is no beta2
  Dart_const_handle get_opposite_ith_real_dart(std::size_t i) const
  {
    CGAL_assertion(i<m_path.size());
    return (get_ith_flip(i)?get_ith_dart(i):
                            get_map().opposite2(get_ith_dart(i)));
  }

  /// @return the first dart of the path, taking into account flip.
  /// @pre !is_empty()
  Dart_const_handle real_front() const
  {
    CGAL_assertion(!is_empty());
    return get_ith_real_dart(0);
  }

  /// @return the last dart of the path, taking into account flip.
  /// @pre !is_empty()
  Dart_const_handle real_back() const
  {
    CGAL_assertion(!is_empty());
    return get_ith_real_dart(length()-1);
  }

  /// @return true iff df can be added at the end of the path.
  bool can_be_pushed(Dart_const_handle dh, bool flip=false) const
  {
    // This assert is too long CGAL_assertion(m_map.darts().owns(dh));

    if (is_empty()) return true;

    return m_map.template belong_to_same_cell<0>
        (m_flip.back() ? back() : m_map.other_extremity(back()),
         flip ? m_map.other_extremity(dh) : dh);
  }

  /// Add the given dart at the end of this path.
  /// @pre can_be_pushed(dh)
  void push_back(Dart_const_handle dh, bool flip=false,
                 bool update_isclosed=true)
  {
    CGAL_assertion(dh!=Map::null_handle);
    /* This assert is too long, it is tested in the is_valid method. */
    //  CGAL_assertion(can_be_pushed(dh, flip));

    m_path.push_back(dh);
    m_flip.push_back(flip);
    if (update_isclosed) { update_is_closed(); }
  }

  /// @return true iff the ith dart can be added at the end of the path.
  bool can_be_pushed_by_index(typename Map::size_type i, bool flip=false,
                              bool update_isclosed=true) const
  { return can_be_pushed(get_map().dart_handle(i), flip, update_isclosed); }

  /// Add the given ith dart at the end of this path.
  void push_back_by_index(typename Map::size_type i, bool flip=false,
                          bool update_isclosed=true)
  { push_back(get_map().dart_handle(i), flip, update_isclosed); }

  void push_back_by_index(std::initializer_list<typename Map::size_type> l,
                          bool update_isclosed=true)
  {
    for (std::size_t i : l)
    { push_back_by_index(i, false, update_isclosed); }
  }

  /// @return true iff the dart labeled e can be added at the end of the path.
  bool can_be_pushed_by_label(const std::string& e, bool flip=false) const
  {
    Dart_const_handle dh=get_map().get_dart_labeled(e);
    if (dh==Map::null_handle) { return false; }
    return can_be_pushed(dh, flip);
  }

  /// Add the dart having the given labels at the end of this path.
  /// Each label is a word, possibly starting by -, words are separated by spaces
  void push_back_by_label(const std::string& s, bool update_isclosed=true)
  {
    std::istringstream iss(s);
    for (std::string e; std::getline(iss, e, ' '); )
    {
      Dart_const_handle dh=get_map().get_dart_labeled(e);
      if (dh!=Map::null_handle) { push_back(dh, false, update_isclosed); }
    }
  }

  void push_back_by_label(std::initializer_list<const char*> l,
                          bool update_isclosed=true)
  {
    for (const char* e : l)
    { push_back_by_label(e, false, update_isclosed); }
  }

  Self& operator+=(const Self& other)
  {
    m_path.reserve(m_path.size()+other.m_path.size());
    // Be careful to the special case when *this==other
    // this is the reason of the iend.
    for (std::size_t i=0, iend=other.length(); i<iend; ++i)
    { push_back(other[i], other.m_flip[i], false); }
    update_is_closed();
    return *this;
  }

  Self operator+(const Self& other) const
  {
    Self res=*this;
    res+=other;
    return res;
  }

  /// change m_path and m_flip in order to get the lower number of flips possible
  void simplify_flips(bool show_flips_left=false)
  {
    if (show_flips_left)
    { std::cout<<"Flips left (maybe none) : "<<std::flush; }
    for(unsigned int i=0; i<length(); ++i)
    {
      if (m_flip[i] && !get_map().template is_free<2>(m_path[i]))
      {
        m_path[i]=get_map().opposite2(m_path[i]);
        m_flip[i]=!m_flip[i];
      }
      else if (show_flips_left)
      { std::cout<<i<<" "<<std::flush; }
    }
    if (show_flips_left)
    { std::cout<<std::endl; }
  }

  /// @return the number of flips of this path.
  unsigned int nb_flips()
  {
    unsigned int res=0;
    for (unsigned int i=0; i<length(); ++i)
    { if (m_flip[i]) ++res; }
    return res;
  }

  /// Cut this path to keep only the n first darts.
  void cut(std::size_t n, bool update_isclosed=true)
  {
    if (n>=length()) return;
    m_path.resize(n);
    m_flip.resize(n);
    if (update_isclosed) { update_is_closed(); }
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
      new_path.push_back(get_ith_dart(begin), get_ith_flip(begin), false);
      ++begin;
    }
    update_is_closed();
  }

  /// Debugging method.
  void display_failed_extention(const std::string& /*name_of_function*/)
  {
    // std::cout<<"Cant extend the path this way ("<<name_of_function<<")"
    //          <<std::endl;
  }

  /// Extend the path straight positive.
  /// @pre must be non empty.
  void extend_straight_positive(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty() || nb==0)
    { display_failed_extention("extend_straight_positive"); return; }

    Dart_const_handle dh=back();
    if(back_flip())
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_straight_positive"); return; }
      else
      { dh=get_map().opposite2(dh); }
    }

    for (unsigned int i=0; i<nb; ++i)
    {
      dh=get_map().next(dh);

      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_straight_positive"); return; }
      dh=get_map().next(get_map().opposite2(dh));

      push_back(dh, false, false);
    }

    if (update_isclosed) { update_is_closed(); }
  }

  /// Extend the path straight negative.
  /// @pre must be non empty.
  void extend_straight_negative(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty() || nb==0)
    { display_failed_extention("extend_straight_negative"); return; }

    Dart_const_handle dh=back();
    if(!back_flip())
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_straight_positive"); return; }
      else
      { dh=get_map().opposite2(dh); }
    }

    for (unsigned int i=0; i<nb; ++i)
    {
      dh=get_map().previous(dh);

      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_straight_negative"); return; }
      dh=get_map().previous(get_map().opposite2(dh));

      push_back(dh, true, false);
    }

    if (update_isclosed) { update_is_closed(); }
  }

  /// Extend the path given a positive turn.
  /// @pre must be non empty.
  void extend_positive_turn(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty())
    { display_failed_extention("extend_positive_turn"); return; }

    if (nb==0)
    {
      push_back(back(), !back_flip(), update_isclosed);
      return;
    }

    Dart_const_handle dh=back();
    if(back_flip())
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_positive_turn"); return; }
      else
      { dh=get_map().opposite2(dh); }
    }
    dh=get_map().next(dh);

    for (unsigned int i=1; i<nb; ++i)
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_positive_turn"); return; }
      dh=get_map().next(get_map().opposite2(dh));
    }

    push_back(dh, false, update_isclosed);
  }

  /// Extend the path given a negative turn.
  /// @pre must be non empty.
  void extend_negative_turn(std::size_t nb=1, bool update_isclosed=true)
  {
    if (is_empty()) { display_failed_extention("extend_negative_turn"); return; }

    if (nb==0)
    {
      push_back(back(), !back_flip(), update_isclosed);
      return;
    }

    Dart_const_handle dh=back();
    if(!back_flip())
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_negative_turn"); return; }
      else
      { dh=get_map().opposite2(dh); }
    }
    dh=get_map().previous(dh);

    for (unsigned int i=1; i<nb; ++i)
    {
      if (get_map().template is_free<2>(dh))
      { display_failed_extention("extend_negative_turn"); return; }
      dh=get_map().previous(get_map().opposite2(dh));
    }

    push_back(dh, true, update_isclosed);
  }

  /// Initializes this path to a random starting path.
  /// @pre must be empty.
  bool initialize_random_starting_dart(CGAL::Random& random,
                                       bool update_isclosed=true)
  {
    if (!is_empty() || get_map().is_empty()) { return false; }

    // first select a random edge by taking the lower index of
    // the two darts when it is not a boundary
    typename Map::size_type index=static_cast<typename Map::size_type>
      (random.get_int(0, static_cast<int>(get_map().darts().capacity())));
    while (!get_map().darts().is_used(index) ||
          (!get_map().template is_free<2>(get_map().dart_handle(index)) &&
           get_map().dart_handle(index)>get_map().
           opposite2(get_map().dart_handle(index))))
    {
      ++index;
      if (index==get_map().darts().capacity()) index=0;
    }

    // second we take randomly one of the two darts of this edge
    // (potentially with the help of a flip)
    bool heads_or_tails=random.get_bool();
    if (get_map().template is_free<2>(get_map().dart_handle(index)))
    {
      push_back(get_map().dart_handle(index), heads_or_tails, update_isclosed);
    }
    else
    {
      if (heads_or_tails)
      { push_back(get_map().dart_handle(index), false, update_isclosed); }
      else
      { push_back(get_map().opposite2(get_map().dart_handle(index)),
                  false, update_isclosed); }
    }
    return true;
  }

  /// Initializes this path to a random starting path.
  /// @pre must be empty.
  bool initialize_random_starting_dart(bool update_isclosed=true)
  {
    CGAL::Random& random=get_default_random();
    return initialize_random_starting_dart(random, update_isclosed);
  }

  /// Extends this path with a random dart.
  /// @pre must be non empty.
  bool extend_path_randomly(CGAL::Random& random,
                            bool allow_half_turn=true,
                            bool update_isclosed=true)
  {
    if (is_empty())
    { return initialize_random_starting_dart(random, update_isclosed); }

    if(get_map().template is_free<1>(back()))
    { return false; }

    Dart_const_handle next_vertex;
    if (back_flip())
    { next_vertex=back(); }
    else if (get_map().template is_free<2>(back()))
    { next_vertex=get_map().next(back()); }
    else
    { next_vertex=get_map().opposite2(back()); }

    std::vector<std::pair<Dart_const_handle, bool> > candidats;
    for (auto it=get_map().template darts_of_cell<0>(next_vertex).begin(),
           itend=get_map().template darts_of_cell<0>(next_vertex).end();
           it!=itend; ++it )
    {
      if (back_flip() || !get_map().template is_free<2>(back()))
      {
        candidats.push_back(std::make_pair(it, false));
        if (get_map().template is_free<2>(get_map().previous(it)))
        { candidats.push_back
              (std::make_pair(get_map().previous(it), true)); }
      }
      else
      {
        if (get_map().template is_free<2>(get_map().previous(it)))
        { candidats.push_back
              (std::make_pair(get_map().previous(it), true)); }
        candidats.push_back(std::make_pair(it, false));
      }
    }
    //candidats is now the list of all the darts that can be pushed back to
    // the path (maybe with a flip) the first of them in the list is the
    // opposite of back(), or back() itself if it is 2-free

    std::size_t i=static_cast<std::size_t>
      (random.get_int(allow_half_turn?0:1,static_cast<int>(candidats.size())));
    auto it=candidats.begin();
    for (std::size_t nb=0; nb<i; ++nb, ++it) {}
    push_back(it->first, it->second, update_isclosed);
    return true;
  }

  /// Extends this path with a random dart.
  /// @pre must be non empty.
  bool extend_path_randomly(bool allow_half_turn=false,
                            bool update_isclosed=true)
  {
    CGAL::Random& random=get_default_random();
    return extend_path_randomly(random, allow_half_turn, update_isclosed);
  }

  /// Generates a random path, with a number of darts >= length.
  void generate_random_path(std::size_t length,
                            CGAL::Random& random=get_default_random(),
                            bool allow_half_turns=true,
                            bool update_isclosed=true)
  {
    m_path.reserve(m_path.size()+length);
    for (std::size_t i=0; i<length; ++i)
    { extend_path_randomly(random, allow_half_turns, true); }
    if (update_isclosed) { update_is_closed(); }
  }

  /// Generates a random path.
  template<typename Path>
  void generate_random_path(CGAL::Random& random,
                            bool update_isclosed=true)
  { generate_random_path(random.get_int(1, 10000),
                         random, true, update_isclosed); }

  /// Generates a random path.
  template<typename Path>
  void generate_random_path(std::size_t length,
                            bool update_isclosed=true)
  {
    CGAL::Random& random=get_default_random();
    generate_random_path(length, random, true, update_isclosed);
  }

  /// Generates a random path.
  template<typename Path>
  void generate_random_path(bool update_isclosed=true)
  {
    CGAL::Random& random=get_default_random();
    generate_random_path(random, update_isclosed);
  }

  /// Generates a random closed path.
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

  /// Generates a random closed path.
  void generate_random_closed_path(std::size_t length)
  {
    CGAL::Random& random=get_default_random();
    generate_random_closed_path(length, random);
  }

  /// Generates a random closed path.
  void generate_random_closed_path(CGAL::Random& random)
  { generate_random_closed_path(random.get_int(1, 10000), random); }

  /// Generates a random closed path.
  void generate_random_closed_path()
  {
    CGAL::Random& random=get_default_random();
    generate_random_closed_path(random.get_int(1, 10000), random);
  }

  /// Replace edge [i] by the path of darts along the face.
  /// If this face does not exist (if it is a boundary) then replace the edge
  /// by the face on the other side. Problem of complexity when used many times
  /// (like in update_path_randomly).
  bool push_around_face(std::size_t i, bool update_isclosed=true)
  {
    CGAL_assertion(i<length());

    // It is not possible to push around a perforated face since it changes
    // the homotopy of the path.
    if (get_map().is_perforated(get_ith_dart(i))) { return false; }

    Self p2(get_mesh());

    // 1) We add in p2 the part of the path which is pushed.
    if (get_ith_flip(i))
    {
      Dart_const_handle dh=get_map().next(get_ith_dart(i));
      do
      {
        p2.push_back(dh, false, false);
        dh=get_map().next(dh);
      }
      while(dh!=get_ith_dart(i));
    }
    else
    {
      Dart_const_handle dh=get_map().previous(get_ith_dart(i));
      do
      {
        p2.push_back(dh, true, false);
        dh=get_map().previous(dh);
      }
      while(dh!=get_ith_dart(i));
    }

    // 2) We copy the end of the path.
    p2.m_path.reserve(p2.length()+length()-i);
    for (std::size_t j=i+1; j<length(); ++j)
    { p2.push_back(get_ith_dart(j), get_ith_flip(j), false); }

    // 3) We cut this path to keep the first i darts.
    cut(i, false);
    m_path.reserve(length()+p2.length());
    for (std::size_t j=0; j<p2.length(); ++j)
    { push_back(p2[j], p2.get_ith_flip(j), false); }

    if (update_isclosed) { update_is_closed(); }
    return true;

    //CGAL_assertion(is_valid());
  }

  /// Transform the current path by pushing some dart around faces.
  /// At the end, the new path is homotopic to the original one.
  void update_path_randomly(std::size_t nb, CGAL::Random& random,
                            bool update_isclosed=true)
  {
    if (is_empty()) return;

    for (unsigned int i=0; i<nb; ++i)
    {
      std::size_t dartn=static_cast<std::size_t>
        (random.get_int(0, static_cast<int>(length())));
      std::size_t j=dartn;
      while(!push_around_face(dartn, false) && dartn!=j)
      { ++dartn; }
    }
    if (update_isclosed) { update_is_closed(); }
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
    CGAL::Random& random=get_default_random();
    update_path_randomly(random, update_isclosed);
  }

  /// @return true iff the i-th dart of the path and the j-th dart of the other
  /// are the same (taking into account the flips !)
  bool are_same_step(std::size_t i, const Self& other, std::size_t j) const
  {
    if (get_ith_flip(i)==other.get_ith_flip(j))
    { return get_ith_dart(i)==other[j]; }

    if (get_map().template is_free<2>(get_ith_dart(i)) ||
        get_map().template is_free<2>(other[j]))
    { return false; }

    return get_ith_dart(i)==get_map().opposite2(other[j]);
  }

  /// @return true if this path is equal to other path, identifying dart 0 of
  ///          this path with dart start in other path.
  bool are_same_paths_from(const Self& other, std::size_t start) const
  {
    CGAL_assertion(start==0 || start<length());
    CGAL_assertion(is_closed() || start==0);
    CGAL_assertion(length()==other.length() && is_closed()==other.is_closed());

    for(std::size_t i=0; i<length(); ++i)
    {
      if (!are_same_step(i, other, start))
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

    Self pp1=*this;
    pp1.simplify_flips();
    Self pp2=other;
    pp2.simplify_flips();
    pp2+=pp2;
    // Now we search if pp1 is a sub-motif of pp2 <=> *this==other

    return boost::algorithm::knuth_morris_pratt_search(pp2.m_path.begin(),
                                                       pp2.m_path.end(),
                                                       pp1.m_path.begin(),
                                                       pp1.m_path.end())
#if BOOST_VERSION>=106200
      .first
#endif
        !=pp2.m_path.end();
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
  bool is_valid(bool display_error=false) const
  {
    if (is_empty()) { return !is_closed(); } // an empty past is not closed
    Dart_const_handle last_vertex;

    for (unsigned int i=1; i<m_path.size(); ++i)
    {
      /* This assert is long if (!m_map.darts().owns(m_path[i]))
      { return false; } */

      if (m_path[i]==Map::null_handle || m_path[i]==m_map.null_dart_handle)
      { return false; }

      last_vertex=m_flip[i-1]?m_path[i-1]:get_map().next(m_path[i-1]);
      if (last_vertex==Map::null_handle)
      {
        if (display_error)
        { std::cout<<"Invalid path: one of the vertices doesn't exist"
                   <<std::endl; }
        return false;
      }

      if (!m_map.template belong_to_same_cell<0>
          (m_flip[i]?get_map().next(m_path[i]):m_path[i], last_vertex))
      {
        if (display_error)
        { std::cout<<"Invalid path: dart "<<i-1<<" and dart "<<i
                   <<" are not adjacents"<<std::endl; }
        return false;
      }
    }
    last_vertex=back_flip()?back():get_map().next(back());
    if (is_closed())
    {
      if (last_vertex==Map::null_handle)
      {
        if (display_error)
        { std::cout<<"Invalid path: one of the vertices doesn't exist"
                   <<std::endl; }
        return false;
      }
      if (!m_map.template belong_to_same_cell<0>
          (front_flip()?get_map().next(front()):front(), last_vertex))
      {
        if (display_error)
        { std::cout<<"Invalid path: m_is_closed is true but the path is "
                   <<"not closed"<<std::endl; }
        return false;
      }
    }
    else
    {
      if (last_vertex==Map::null_handle)
      {
        if (display_error)
        { std::cout<<"Invalid path: one of the vertices doesn't exist"
                   <<std::endl; }
        return false;
      }
      if (m_map.template belong_to_same_cell<0>
          (front_flip()?get_map().next(front()):front(), last_vertex))
      {
        if (display_error)
        { std::cout<<"Invalid path: m_is_closed is false but the path "
                   <<"is closed"<<std::endl; }
        return false;
      }
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
      Dart_const_handle
          pend=m_flip.back()?back():m_map.other_extremity(back());
      if (pend==Map::null_handle) { m_is_closed=false; }
      else
      {
        Dart_const_handle
            pbegin=m_flip[0]?m_map.other_extremity(m_path[0]):m_path[0];
        m_is_closed=m_map.template belong_to_same_cell<0>(pbegin, pend);
      }
    }
  }

  /// @return true iff the path does not pass twice through a same edge
  ///              or a same vertex.
  bool is_simple() const
  {
    typename Map::size_type markvertex=m_map.get_new_mark();
    typename Map::size_type markedge=m_map.get_new_mark();

    bool res=true;
    Dart_const_handle dh_vertex;
    unsigned int i=0;
    for (i=0; res && i<m_path.size(); ++i)
    {
      dh_vertex=m_flip[i]?get_map().next(m_path[i]):m_path[i];
      if (m_map.is_marked(dh_vertex, markvertex)) { res=false; }
      else { CGAL::mark_cell<Map, 0>(m_map, dh_vertex, markvertex); }

      if (m_map.is_marked(m_path[i], markedge)) { res=false; }
      else  { CGAL::mark_cell<Map, 1>(m_map, m_path[i], markedge); }
    }

    i=0;
    while(m_map.number_of_marked_darts(markedge)>0 ||
          m_map.number_of_marked_darts(markvertex)>0)
    {
      CGAL_assertion(i<m_path.size());
      dh_vertex=m_flip[i]?get_map().next(m_path[i]):m_path[i];
      if (m_map.is_marked(dh_vertex, markvertex))
      { CGAL::unmark_cell<Map, 0>(m_map, dh_vertex, markvertex); }
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
    bool tmpbool;
    for (unsigned int i=0; i<length()/2; ++i)
    {
      std::swap(m_path[i], m_path[length()-1-i]);
      tmpbool=m_flip[i]; // Cannot swap in vector bool
      m_flip[i]=m_flip[length()-1-i];
      m_flip[length()-1-i]=tmpbool;
    }

    for (unsigned int i=0; i<length(); ++i)
    { m_flip[i]=!m_flip[i]; }
  }

  /// If the given path is opened, close it by doing the same path that the
  /// first one in reverse direction.
  void close()
  { // TODO follow shortest path ?
    if (!is_closed())
    {
      for (int i=m_path.size()-1; i>=0; --i)
      { m_path.push_back(m_path[i], !m_flip[i], false); }
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

    if ((get_ith_flip(i) && get_map().template is_free<2>(get_ith_dart(i))) ||
        (get_next_flip(i) && get_map().template is_free<2>(get_next_dart(i))))
    { return (std::numeric_limits<std::size_t>::max)(); }

    return m_map.positive_turn(get_ith_real_dart(i),
                               get_ith_real_dart(next_index(i)));
  }

  /// Same than next_positive_turn but turning in reverse orientation
  /// around vertex.
  std::size_t next_negative_turn(std::size_t i) const
  {
    // CGAL_assertion(is_valid());
    CGAL_assertion(i<m_path.size());
    CGAL_assertion (is_closed() || i<length()-1);

    if ((!get_ith_flip(i) && get_map().template is_free<2>(get_ith_dart(i))) ||
        (!get_next_flip(i) && get_map().template is_free<2>(get_next_dart(i))))
    { return (std::numeric_limits<std::size_t>::max)(); }

    return m_map.positive_turn(get_opposite_ith_real_dart(next_index(i)),
                               get_opposite_ith_real_dart(i));
  }

  /// Computes all positive turns of this path.
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

  /// Computes all negative turns of this path.
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

  /// Computes all positive or negative turns of this path, depending on p.
  std::vector<std::size_t> compute_turns(bool p) const
  { return (p?compute_positive_turns():compute_negative_turns()); }

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
      if ((nb>=0 && resplus[start]!=static_cast<std::size_t>(nb)) ||
          (nb<0 && resmoins[start]!=static_cast<std::size_t>(-nb)))
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
      if (m_flip[i])
      { std::cout<<"f"; }
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
  const typename Get_map<Mesh, Mesh>::storage_type m_map; // The underlying map
  std::vector<Dart_const_handle> m_path; /// The sequence of darts
  bool m_is_closed;                      /// True iff the path is a cycle
  std::vector<bool> m_flip;              /// The sequence of flips
};

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_PATH_ON_SURFACE_H //
// EOF //
