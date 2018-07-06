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

} // namespace CGAL

#endif // CGAL_PATH_GENERATORS_H //
// EOF //
