// Copyright (c) 2014 GeometryFactory
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
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_EXPANSION_EVENT_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_EXPANSION_EVENT_H

#include <CGAL/license/Surface_mesh_shortest_path.h>


namespace CGAL {

namespace internal {

template<class Traits>
class Cone_tree_node;

template <class Traits>
struct Cone_expansion_event
{
public:
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::FT FT;

public:
  enum Expansion_type
  {
    LEFT_CHILD,
    RIGHT_CHILD,
    PSEUDO_SOURCE
  };

public:
  Cone_tree_node<Traits>* m_parent;
  FT m_distanceEstimate;
  Expansion_type m_type;
  Segment_2 m_windowSegment;
  bool m_cancelled;

public:
  Cone_expansion_event(Cone_tree_node<Traits>* parent, const FT& distanceEstimate, Expansion_type type)
    : m_parent(parent)
    , m_distanceEstimate(distanceEstimate)
    , m_type(type)
    , m_cancelled(false)
  {
  }

  Cone_expansion_event(Cone_tree_node<Traits>* parent, const FT& distanceEstimate, Expansion_type type, const Segment_2& windowSegment)
    : m_parent(parent)
    , m_distanceEstimate(distanceEstimate)
    , m_type(type)
    , m_windowSegment(windowSegment)
    , m_cancelled(false)
  {
  }
};

// Does the opposite of what you would expect in order to implement a min-priority queue
template <class Traits>
struct Cone_expansion_event_min_priority_queue_comparator
{
public:
  bool operator () (const Cone_expansion_event<Traits>* lhs, const Cone_expansion_event<Traits>* rhs) const
  {
    return rhs->m_distanceEstimate < lhs->m_distanceEstimate;
  }
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_EXPANSION_EVENT_H
