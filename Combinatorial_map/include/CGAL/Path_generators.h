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
void generate_random_path(Path& p, std::size_t length, CGAL::Random& random)
{
  for (unsigned int i=0; i<length; ++i)
  { extend_path_randomly(p, random); }
}

template<typename Path>
void generate_random_path(Path& p, std::size_t length)
{
  CGAL::Random random;
  generate_random_path(p, length, random);
}
  
template<typename Path>
bool initialize_path_random_starting_dart(Path& p, CGAL::Random& random)
{
  p.clear();
  
  unsigned int index=random.get_int(0, p.get_map().darts().capacity());
  while (!p.get_map().darts().is_used(index))
  {
    ++index;
    if (index==p.get_map().darts().capacity()) index=0;
  }
  p.push_back(p.get_map().darts().iterator_to(p.get_map().darts()[index]));
  return true;
}
  
template<typename Path>
bool initialize_path_random_starting_dart(Path& p)
{
  CGAL::Random random;
  return initialize_path_random_starting_dart(p, random);
}
  
template<typename Path>
bool extend_path_randomly(Path& p, CGAL::Random& random, bool allow_half_turn=false)
{
  if (p.is_empty())
  { return initialize_path_random_starting_dart(p, random); }

  typename Path::Dart_const_handle pend=p.get_map().opposite(p.back());
  if (pend==Path::Map::null_handle)
  {
    if (!p.get_map().template is_free<1>(p.back()))
    {
      p.push_back(p.get_map().template beta<1>(p.back()));
      return true;
    }
    else { return false; }
  }
  
  typename Path::Map::template Dart_of_cell_range<0>::const_iterator
    it=p.get_map().template darts_of_cell<0>(pend).begin();
  
  unsigned int index=random.get_int
    ((allow_half_turn?0:1), p.get_map().template darts_of_cell<0>(pend).size());
  for(unsigned int i=0; i<index; ++i)
  { ++it; }
  
  assert(allow_half_turn || it!=pend);
  
  p.push_back(it);
  return true;
}
  
template<typename Path>
bool extend_path_randomly(Path& p, bool allow_half_turn=false)
{
  CGAL::Random random;
  extend_path_randomly(p, random, allow_half_turn);
}

template<typename Path>  
void extend_straight_positive(Path& p, std::size_t nb=1)
{
  if (p.is_empty() || nb==0)
  { return; }

  typename Path::Dart_const_handle d2;
  for (std::size_t i=0; i<nb; ++i)
  {
    d2=p.get_map().template beta<1,2,1>(p.back());
    if (d2!=p.get_map().null_dart_handle)
    { p.push_back(d2); }
  }
}

template<typename Path>  
void extend_straight_negative(Path& p, std::size_t nb=1)
{
  if (p.is_empty() || nb==0)
  { return; }
  
  typename Path::Dart_const_handle d2;
  for (std::size_t i=0; i<nb; ++i)
  {
    d2=p.get_map().template beta<2,0,2,0,2>(p.back());
    if (d2!=p.get_map().null_dart_handle)
    { p.push_back(d2); }
  }
}

template<typename Path>
void extend_straight_positive_until(Path& p,
                                    typename Path::Dart_const_handle dend)
{
  if (p.is_empty() || p.back()==dend)
  { return; }

  typename Path::Dart_const_handle
      d2=p.get_map().template beta<1,2,1>(p.back());
  while(d2!=dend)
  {
    p.push_back(d2);
    d2=p.get_map().template beta<1,2,1>(d2);
  }
}

template<typename Path>
void extend_straight_negative_until(Path& p,
                                    typename Path::Dart_const_handle dend)
{
  if (p.is_empty() || p.back()==dend)
  { return; }

  typename Path::Dart_const_handle
      d2=p.get_map().template beta<2,0,2,0,2>(p.back());
  while(d2!=dend)
  {
    p.push_back(d2);
    d2=p.get_map().template beta<2,0,2,0,2>(d2);
  }
}

template<typename Path>
void extend_uturn_positive(Path& p, std::size_t nb=1)
{
  if (p.is_empty() || nb==0)
  { return; }

  typename Path::Dart_const_handle d2=p.get_map().template beta<1>(p.back());
  for (std::size_t i=1; i<nb; ++i)
  { d2=p.get_map().template beta<2, 1>(d2); }

  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2); }
}

template<typename Path>  
void extend_uturn_negative(Path& p, std::size_t nb=1)
{
  if (p.is_empty())
  { return; }
  
  typename Path::Dart_const_handle d2=p.get_map().template beta<2>(p.back());
  for (std::size_t i=0; i<nb; ++i)
  { d2=p.get_map().template beta<0, 2>(d2); }

  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2); }
}
  
template<typename Path>
void extend_uturn_half_turn(Path& p)
{
  if (p.is_empty())
  { return; }

  typename Path::Dart_const_handle d2=p.get_map().template beta<2>(p.back());
  if (d2!=p.get_map().null_dart_handle)
  { p.push_back(d2); }
}

template<typename Path>
void create_braket_positive(Path& p, std::size_t length, CGAL::Random& random)
{
  if (p.is_empty())
  { initialize_path_random_starting_dart(p, random); }
  
  extend_uturn_positive(p);
  for (std::size_t i=0; i<length; ++i)
  { extend_straight_positive(p); }
  extend_uturn_positive(p);
}
  
template<typename Path>  
void create_braket_positive(Path& p, std::size_t length)
{
  CGAL::Random random;
  create_braket_positive(p, length, random);
}
  
template<typename Path>  
void create_braket_negative(Path& p, std::size_t length, CGAL::Random& random)
{
  if (p.is_empty())
  { initialize_path_random_starting_dart(p, random); }
  
  extend_uturn_negative(p);
  for (std::size_t i=0; i<length; ++i)
  { extend_straight_negative(p); }
  extend_uturn_negative(p);
}
  
template<typename Path>  
void create_braket_negative(Path& p, std::size_t length)
{
  CGAL::Random random;
  create_braket_negative(p, length, random);
}
  
} // namespace CGAL

#endif // CGAL_PATH_GENERATORS_H //
// EOF //
