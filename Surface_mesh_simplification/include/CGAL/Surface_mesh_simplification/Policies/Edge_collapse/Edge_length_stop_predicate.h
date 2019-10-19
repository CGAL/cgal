// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGTH_STOP_PREDICATE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGTH_STOP_PREDICATE_H 1

#include <CGAL/license/Surface_mesh_simplification.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

template <class FT>
class Edge_length_stop_predicate
{
  FT m_edge_sq_length_threshold;
public:
  Edge_length_stop_predicate( double edge_length_threshold )
    : m_edge_sq_length_threshold(edge_length_threshold*edge_length_threshold)
  {}

  template <typename F, typename Profile>
  bool operator()( F const&
                 , Profile const&  profile
                 , std::size_t
                 , std::size_t
                 ) const
  {
    return  CGAL::squared_distance(profile.p0(), profile.p1()) >
            m_edge_sq_length_threshold;
  }
};

} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_EDGE_LENGTH_STOP_PREDICATE_H //

