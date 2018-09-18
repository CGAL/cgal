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
#ifndef CGAL_PATH_GENERATORS_H
#define CGAL_PATH_GENERATORS_H 1

#include<CGAL/Random.h>
#include<unordered_set>
#include<unordered_map>

namespace CGAL {

template<typename Path>
void generate_random_path(Path& p, std::size_t length, CGAL::Random& random,
                          bool update_isclosed=true)
{
  for (unsigned int i=0; i<length; ++i)
  { extend_path_randomly(p, random, true, false); }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void generate_random_path(Path& p, CGAL::Random& random,
                          bool update_isclosed=true)
{ generate_random_path(p, random.get_int(0, 10000), random, update_isclosed); }

template<typename Path>
void generate_random_path(Path& p, std::size_t length,
                          bool update_isclosed=true)
{
  CGAL::Random random;
  generate_random_path(p, length, random, update_isclosed);
}

template<typename Path>
void generate_random_path(Path& p, bool update_isclosed=true)
{
  CGAL::Random random;
  generate_random_path(p, random, update_isclosed);
}

template<typename Path>
bool initialize_path_random_starting_dart(Path& p, CGAL::Random& random,
                                          bool update_isclosed=true)
{
  p.clear();
  
  unsigned int index=random.get_int(0, p.get_map().darts().capacity());
  while (!p.get_map().darts().is_used(index))
  {
    ++index;
    if (index==p.get_map().darts().capacity()) index=0;
  }
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[index]),
              update_isclosed);
  return true;
}
  
template<typename Path>
bool initialize_path_random_starting_dart(Path& p, bool update_isclosed=true)
{
  CGAL::Random random;
  return initialize_path_random_starting_dart(p, random, update_isclosed);
}
  
template<typename Path>
bool extend_path_randomly(Path& p, CGAL::Random& random,
                          bool allow_half_turn=true,
                          bool update_isclosed=true)
{
  if (p.is_empty())
  { return initialize_path_random_starting_dart(p, random, update_isclosed); }

  typename Path::Dart_const_handle pend=p.get_map().opposite(p.back());
  if (pend==Path::Map::null_handle)
  {
    if (!p.get_map().template is_free<1>(p.back()))
    {
      p.push_back(p.get_map().template beta<1>(p.back()), update_isclosed);
      return true;
    }
    else { return false; }
  }
  
  typename Path::Map::template Dart_of_cell_range<0>::const_iterator
    it=p.get_map().template darts_of_cell<0>(pend).begin();
  
  unsigned int index=random.get_int
      ((allow_half_turn?0:1),
       p.get_map().template darts_of_cell<0>(pend).size());
  for(unsigned int i=0; i<index; ++i)
  { ++it; }
  
  assert(allow_half_turn || it!=pend);
  
  p.push_back(it, update_isclosed);
  return true;
}
  
template<typename Path>
bool extend_path_randomly(Path& p, bool allow_half_turn=false,
                          bool update_isclosed=true)
{
  CGAL::Random random;
  extend_path_randomly(p, random, allow_half_turn, update_isclosed);
}

template<typename Path>  
void extend_straight_positive(Path& p, std::size_t nb=1,
                              bool update_isclosed=true)
{
  if (p.is_empty() || nb==0)
  { return; }

  typename Path::Dart_const_handle d2;
  for (std::size_t i=0; i<nb; ++i)
  {
    d2=p.get_map().template beta<1,2,1>(p.back());
    if (d2!=p.get_map().null_dart_handle)
    { p.push_back(d2, false); }
  }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>  
void extend_straight_negative(Path& p, std::size_t nb=1,
                              bool update_isclosed=true)
{
  if (p.is_empty() || nb==0)
  { return; }
  
  typename Path::Dart_const_handle d2;
  for (std::size_t i=0; i<nb; ++i)
  {
    d2=p.get_map().template beta<2,0,2,0,2>(p.back());
    if (d2!=p.get_map().null_dart_handle)
    { p.push_back(d2, false); }
  }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void extend_straight_positive_until(Path& p,
                                    typename Path::Dart_const_handle dend,
                                    bool update_isclosed=true)
{
  if (p.is_empty() || p.back()==dend)
  { return; }

  typename Path::Dart_const_handle
      d2=p.get_map().template beta<1,2,1>(p.back());
  while(d2!=dend)
  {
    p.push_back(d2, false);
    d2=p.get_map().template beta<1,2,1>(d2);
  }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void extend_straight_negative_until(Path& p,
                                    typename Path::Dart_const_handle dend,
                                    bool update_isclosed=true)
{
  if (p.is_empty() || p.back()==dend)
  { return; }

  typename Path::Dart_const_handle
      d2=p.get_map().template beta<2,0,2,0,2>(p.back());
  while(d2!=dend)
  {
    p.push_back(d2, false);
    d2=p.get_map().template beta<2,0,2,0,2>(d2);
  }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void extend_uturn_positive(Path& p, std::size_t nb=1,
                           bool update_isclosed=true)
{
  if (p.is_empty() || nb==0)
  { return; }

  typename Path::Dart_const_handle d2=p.get_map().template beta<1>(p.back());
  for (std::size_t i=1; i<nb; ++i)
  { d2=p.get_map().template beta<2, 1>(d2); }

  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2, update_isclosed); }
}

template<typename Path>  
void extend_uturn_negative(Path& p, std::size_t nb=1, bool update_isclosed=true)
{
  if (p.is_empty())
  { return; }
  
  typename Path::Dart_const_handle d2=p.get_map().template beta<2>(p.back());
  for (std::size_t i=0; i<nb; ++i)
  { d2=p.get_map().template beta<0, 2>(d2); }

  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2, update_isclosed); }
}
  
template<typename Path>
void extend_uturn_half_turn(Path& p, bool update_isclosed=true)
{
  if (p.is_empty())
  { return; }

  typename Path::Dart_const_handle d2=p.get_map().template beta<2>(p.back());
  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2, update_isclosed); }
}

template<typename Path>
void create_braket_positive(Path& p, std::size_t length, CGAL::Random& random,
                            bool update_isclosed=true)
{
  if (p.is_empty())
  { initialize_path_random_starting_dart(p, random, false); }
  
  extend_uturn_positive(p, 1, false);
  extend_straight_positive(p, length, false);
  extend_uturn_positive(p, 1, false);
  if (update_isclosed) { p.update_is_closed(); }
}
  
template<typename Path>  
void create_braket_positive(Path& p, std::size_t length,
                            bool update_isclosed=true)
{
  CGAL::Random random;
  create_braket_positive(p, length, random, update_isclosed);
}
  
template<typename Path>  
void create_braket_negative(Path& p, std::size_t length, CGAL::Random& random,
                            bool update_isclosed=true)
{
  if (p.is_empty())
  { initialize_path_random_starting_dart(p, random, false); }
  
  extend_uturn_negative(p, 1, false);
  extend_straight_negative(p, length, false);
  extend_uturn_negative(p, 1, false);
  if (update_isclosed) { p.update_is_closed(); }
}
  
template<typename Path>  
void create_braket_negative(Path& p, std::size_t length,
                            bool update_isclosed=true)
{
  CGAL::Random random;
  create_braket_negative(p, length, random, update_isclosed);
}

template<typename Path>
void push_around_face(Path& p, std::size_t i, bool update_isclosed=true)
{
  std::size_t begin=i, end=i;
  while (p.get_map().template beta<1>(p.get_prev_dart(begin))==
         p.get_ith_dart(begin))
  {
    begin=p.prev_index(begin);
    if (begin==i)
    { return; } // Case of a path that is equal to a face
  }

  while (p.get_map().template beta<1>(p.get_ith_dart(end))==
         p.get_next_dart(end))
  {
    end=p.next_index(end);
    assert(end!=i);
  }

  Path p2(p.get_map());
  typename Path::Dart_const_handle
      dh=p.get_map().template beta<0>(p.get_ith_dart(begin));
  do
  {
    p2.push_back(p.get_map().template beta<2>(dh));
    dh=p.get_map().template beta<0>(dh);
  }
  while(dh!=p.get_ith_dart(end));
  for (std::size_t i=end+1; i<p.length(); ++i)
  { p2.push_back(p.get_ith_dart(i), false); }

  p.cut(begin, false);
  for (std::size_t i=0; i<p2.length(); ++i)
  { p.push_back(p2[i], false); }

  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void update_path_randomly(Path& p, std::size_t nb, CGAL::Random& random,
                          bool update_isclosed=true)
{
  for (unsigned int i=0; i<nb; ++i)
  {
    push_around_face(p, random.get_int(0, p.length()), false);
  }
  if (update_isclosed) { p.update_is_closed(); }
}

template<typename Path>
void update_path_randomly(Path& p, CGAL::Random& random,
                          bool update_isclosed=true)
{ update_path_randomly(p, random.get_int(0, 10000), random, update_isclosed); }

template<typename Path>
void update_path_randomly(Path& p, std::size_t nb, bool update_isclosed=true)
{
  CGAL::Random random;
  update_path_randomly(p, nb, random, update_isclosed);
}

template<typename Path>
void update_path_randomly(Path& p, bool update_isclosed=true)
{
  CGAL::Random random;
  update_path_randomly(p, random, update_isclosed);
}

template<typename LCC>
typename LCC::Dart_const_handle
generate_random_connected_set_of_faces(const LCC& lcc, std::size_t nb,
                                       CGAL::Random& random,
                                       std::unordered_set<typename LCC::Dart_const_handle>& set,
                                       typename LCC::size_type amark)
{
  set.clear();
  if (lcc.is_empty()) { return NULL; }

  std::unordered_map<std::size_t, typename LCC::Dart_const_handle> border_faces;
  
  unsigned int index=random.get_int(0, lcc.darts().capacity());
  while (!lcc.darts().is_used(index))
  {
    ++index;
    if (index==lcc.darts().capacity()) { index=0; }
  }

  typename LCC::Dart_const_handle dh1=lcc.darts().iterator_to(lcc.darts()[index]);
  border_faces[0]=dh1;
  set.insert(dh1);
  CGAL::mark_cell<LCC, 2>(lcc, dh1, amark);
  
  for (std::size_t i=1; i<nb; ++i)
  {
    std::size_t facenumber=(std::size_t)(random.get_int(0, border_faces.size()));
    std::size_t nbborder=0;
    
    typename LCC::Dart_const_handle dh1_init=border_faces[facenumber];
    dh1=dh1_init;
    do
    {
      if (!lcc.template is_free<2>(dh1) &&
          !lcc.is_marked(lcc.template beta<2>(dh1), amark))
      { ++nbborder; }
      dh1=lcc.template beta<1>(dh1);
    }
    while (dh1!=dh1_init);

    while(lcc.template is_free<2>(dh1) ||
          lcc.is_marked(lcc.template beta<2>(dh1), amark))
    { dh1=lcc.template beta<1>(dh1); }

    std::size_t dartnumber=(std::size_t)(random.get_int(0, nbborder));
    for (std::size_t j=0; j<dartnumber;)
    {
      if (!lcc.template is_free<2>(dh1) &&
          !lcc.is_marked(lcc.template beta<2>(dh1), amark))
      { ++j; }
      dh1=lcc.template beta<1>(dh1);
      while(lcc.template is_free<2>(dh1) ||
            lcc.is_marked(lcc.template beta<2>(dh1), amark))
      { dh1=lcc.template beta<1>(dh1); }
    }

    // Here we have a new face
    set.insert(lcc.template beta<2>(dh1));
    CGAL::mark_cell<LCC, 2>(lcc, lcc.template beta<2>(dh1), amark);

    // We add it in the list of borders faces
    border_faces[border_faces.size()]=lcc.template beta<2>(dh1);

    // Then we update the list of border faces (because some of them could be
    // no more border due to the adding of the new face)
    std::unordered_map<std::size_t, typename LCC::Dart_const_handle> border_faces_new;
    for (typename std::unordered_map<std::size_t, typename LCC::Dart_const_handle>::iterator
         it=border_faces.begin(), itend=border_faces.end(); it!=itend; ++it)
    {
      bool isborder=false;
      dh1=it->second;
      do
      {
        if (!lcc.template is_free<2>(dh1) &&
            !lcc.is_marked(lcc.template beta<2>(dh1), amark))
        { isborder=true; }
        else
        { dh1=lcc.template beta<1>(dh1); }
      }
      while(!isborder && dh1!=it->second);
      if (isborder)
      { border_faces_new[border_faces_new.size()]=dh1; }
    }
    std::swap(border_faces, border_faces_new);

    if (border_faces.size()==0)
    { return NULL; }
  }

  assert (border_faces.size()!=0);
  typename LCC::Dart_const_handle dhres=border_faces[0];
  while(lcc.template is_free<2>(dhres) ||
        lcc.is_marked(lcc.template beta<2>(dhres), amark))
  { dhres=lcc.template beta<1>(dhres); }

  return dhres;
}

template<typename Path>
void generate_random_closed_path(Path& p, std::size_t nb,
                                 CGAL::Random& random)
{
  std::unordered_set<typename Path::Map::Dart_const_handle> faces;
  typename Path::Map::size_type amark=p.get_map().get_new_mark();

  typename Path::Map::Dart_const_handle dhi=
      generate_random_connected_set_of_faces(p.get_map(), nb, random,
                                             faces, amark);

  if (dhi==NULL)
  {
    p.get_map().free_mark(amark);
    return;  // We have selected all the faces.
  }

  typename Path::Map::Dart_const_handle dh=dhi;
  do
  {
    assert(p.get_map().template is_free<2>(dh) ||
           !p.get_map().is_marked(p.get_map().template beta<2>(dh), amark));
    p.push_back(dh, false);
    dh=p.get_map().template beta<1>(dh);
    while(!p.get_map().template is_free<2>(dh) &&
          p.get_map().is_marked(p.get_map().template beta<2>(dh), amark))
    { dh=p.get_map().template beta<2, 1>(dh); }
  }
  while(dh!=dhi);

  for (typename std::template unordered_set<typename Path::Map::Dart_const_handle>::iterator
       it=faces.begin(), itend=faces.end(); it!=itend; ++it)
  { CGAL::unmark_cell<typename Path::Map, 2>(p.get_map(), *it, amark); }

  p.get_map().free_mark(amark);

  p.update_is_closed(); // TODO we can avoid that because we know that we generated a closed path (to do so, we need a method that put p.is_closed to true)
  assert(p.is_closed());
}

} // namespace CGAL

#endif // CGAL_PATH_GENERATORS_H //
// EOF //
