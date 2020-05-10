// Copyright (c) 2020 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Thien Hoang <thienvhoang99@gmail.com>
//
#ifndef CGAL_EDGE_WEIGHT_FUNCTOR_H
#define CGAL_EDGE_WEIGHT_FUNCTOR_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL {
namespace Surface_mesh_topology {

// Unit weight, all edges have the same weight: 1.
struct Unit_weight_functor
{
  using Weight_t=unsigned int;
  template <class T>
  Weight_t operator() (T) const { return 1; }
};

// Euclidean distance weight functor: each edge has as weight its Euclidean length.
template<typename Mesh>
struct Euclidean_length_weight_functor
{
  using Weight_t=double;
  using Dart_const_handle=typename Get_map<Mesh, Mesh>::type::Dart_const_handle;

  Euclidean_length_weight_functor(const Mesh& m) : m_mesh(m), m_map(m)
  {}

  Weight_t operator() (Dart_const_handle dh) const
  {
    return CGAL::sqrt(CGAL::squared_distance
                      (Get_traits<Mesh>::get_point(m_mesh, dh),
                       Get_traits<Mesh>::get_point(m_mesh, m_map.other_extremity(dh))));
  }

protected:
  const Mesh& m_mesh;
  const typename Get_map<Mesh, Mesh>::storage_type m_map;
};

} // namespace Surface_mesh_topology
} // namespace CGAL

#endif // CGAL_EDGE_WEIGHT_FUNCTOR_H
