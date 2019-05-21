// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_KSR_3_EVENT_H
#define CGAL_KSR_3_EVENT_H

//#include <CGAL/license/Kinetic_shape_reconstruction.h>

#include <CGAL/KSR/utils.h>

namespace CGAL
{

namespace KSR_3
{

template <typename Data>
class Event_queue;

template <typename Data>
class Event
{
public:
  typedef typename Data::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Data::PVertex PVertex;
  typedef typename Data::PEdge PEdge;
  typedef typename Data::PFace PFace;
  typedef typename Data::IVertex IVertex;
  typedef typename Data::IEdge IEdge;
  
  typedef Event_queue<Data> Queue;
  friend Queue;

  enum Type
  {
    FREE_VERTEX_TO_INTERSECTION_EDGE,
    CONSTRAINED_VERTEX_TO_FREE_VERTEX,
    CONSTRAINED_VERTEX_TO_INTERSECTION_VERTEX,
    CONSTRAINED_VERTEX_TO_CONSTRAINED_VERTEX,
    EDGE_TO_INTERSECTION_EDGE
  };

private:

  PVertex m_pvertex;
  IEdge m_iedge;
  FT m_time;

public:

  Event () { }

  Event (PVertex pvertex, IEdge iedge, FT time)
    : m_pvertex (pvertex)
    , m_iedge (iedge), m_time (time)
  { }

  PVertex pvertex() const { return m_pvertex; }
  IEdge iedge() const { return m_iedge; }
  FT time() const { return m_time; }
  
  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    os << "Event at t=" << ev.m_time << " between vertex ("
       << ev.m_pvertex.first << ":" << ev.m_pvertex.second
       << ") and intersection edge " << ev.m_iedge;
    return os;
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_EVENT_H
