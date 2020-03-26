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
#ifndef CGAL_PATH_GENERATORS_H
#define CGAL_PATH_GENERATORS_H 1

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Random.h>
#include <unordered_set>
#include <unordered_map>
#include <cstdlib>

namespace CGAL {
namespace Surface_mesh_topology {
namespace internal {

template<typename Path>
void create_braket_positive(Path& p, std::size_t length, CGAL::Random& random,
                            bool update_isclosed=true)
{
  if (p.is_empty())
  { p.initialize_random_starting_dart(random, false); }

  p.extend_positive_turn(1, false);
  p.extend_straight_positive(length, false);
  p.extend_positive_turn(1, false);
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
  { p.initialize_random_starting_dart(random, false); }

  p.extend_negative_turn(1, false);
  p.extend_straight_negative(length, false);
  p.extend_negative_turn(1, false);
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
void generate_random_positive_bracket(Path& path,
                                      std::size_t nb1,
                                      std::size_t nb2, std::size_t nb3,
                                      CGAL::Random& random)
{
  path.clear();
  path.initialize_random_starting_dart(random);
  path.extend_straight_positive(nb1-1);
  create_braket_positive(path, nb2);
  path.extend_straight_positive(nb3);
  path.generate_random_path(random.get_int(0, 15), random);
}

template<typename LCC>
typename LCC::Dart_const_handle
generate_random_connected_set_of_faces(const LCC& lcc, std::size_t nb,
                                       CGAL::Random& random,
                                       std::unordered_set
                                       <typename LCC::Dart_const_handle>& set,
                                       typename LCC::size_type amark)
{
  set.clear();
  if (lcc.is_empty()) { return NULL; }

  std::unordered_map<std::size_t, typename LCC::Dart_const_handle> border_faces;

  std::size_t index=static_cast<std::size_t>
    (random.get_int(0, static_cast<int>(lcc.darts().capacity())));
  while (!lcc.darts().is_used(index))
  {
    ++index;
    if (index==lcc.darts().capacity()) { index=0; }
  }

  typename LCC::Dart_const_handle dh1=lcc.darts().iterator_to(lcc.darts()[index]);
  border_faces[0]=dh1;
  set.insert(dh1);
  lcc.template mark_cell<2>(dh1, amark);

  for (std::size_t i=1; i<nb; ++i)
  {
    std::size_t facenumber=static_cast<std::size_t>
                           (random.get_int
                            (0, static_cast<int>(border_faces.size())));
    int nbborder=0;

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

    std::size_t dartnumber=static_cast<std::size_t>(random.get_int(0, nbborder));
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
    lcc.template mark_cell<2>(lcc.template beta<2>(dh1), amark);

    // We add it in the list of borders faces
    border_faces[border_faces.size()]=lcc.template beta<2>(dh1);

    // Then we update the list of border faces (because some of them could be
    // no more border due to the adding of the new face)
    std::unordered_map<std::size_t, typename LCC::Dart_const_handle>
        border_faces_new;
    for (typename std::unordered_map<std::size_t,
         typename LCC::Dart_const_handle>::iterator it=border_faces.begin(),
         itend=border_faces.end(); it!=itend; ++it)
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
    CGAL_assertion(p.get_map().template is_free<2>(dh) ||
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
  { p.get_map().template unmark_cell<2>(*it, amark); }

  p.get_map().free_mark(amark);

  p.update_is_closed();
  CGAL_assertion(p.is_closed());
}

} // namespace internal
} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_PATH_GENERATORS_H //
// EOF //
