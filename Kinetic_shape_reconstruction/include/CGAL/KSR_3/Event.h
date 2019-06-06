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
  
private:

  PVertex m_pvertex;

  PVertex m_pother;
  IEdge m_iedge;
  IVertex m_ivertex;

  FT m_time;

public:

  Event ()
    : m_pvertex (Data::null_pvertex())
    , m_pother (Data::null_pvertex())
    , m_iedge (Data::null_iedge())
    , m_ivertex (Data::null_ivertex())
    , m_time (0)
  { }

  Event (PVertex pvertex, PVertex pother, FT time)
    : m_pvertex (pvertex)
    , m_pother (pother)
    , m_iedge (Data::null_iedge())
    , m_ivertex (Data::null_ivertex())
    , m_time (time)
  { }

  Event (PVertex pvertex, IEdge iedge, FT time)
    : m_pvertex (pvertex)
    , m_pother (Data::null_pvertex())
    , m_iedge (iedge)
    , m_ivertex (Data::null_ivertex())
    , m_time (time)
  { }

  Event (PVertex pvertex, IVertex ivertex, FT time)
    : m_pvertex (pvertex)
    , m_pother (Data::null_pvertex())
    , m_iedge (Data::null_iedge())
    , m_ivertex (ivertex)
    , m_time (time)
  { }

  Event (PVertex pvertex, PVertex pother, IVertex ivertex, FT time)
    : m_pvertex (pvertex)
    , m_pother (pother)
    , m_iedge (Data::null_iedge())
    , m_ivertex (ivertex)
    , m_time (time)
  { }

  PVertex pvertex() const { return m_pvertex; }
  PVertex pother() const { return m_pother; }
  IEdge iedge() const { return m_iedge; }
  IVertex ivertex() const { return m_ivertex; }
  FT time() const { return m_time; }

  bool is_pvertex_to_pvertex() const { return (m_pother != Data::null_pvertex()); }
  bool is_pvertex_to_iedge() const { return (m_iedge != Data::null_iedge()); }
  
  bool is_pvertex_to_ivertex() const { return (m_pother == Data::null_pvertex()
                                               && m_ivertex != Data::null_ivertex()); }
  bool is_pvertices_to_ivertex() const { return (m_pother != Data::null_pvertex()
                                                 && m_ivertex != Data::null_ivertex()); }
  
  friend std::ostream& operator<< (std::ostream& os, const Event& ev)
  {
    if (ev.is_pvertex_to_pvertex())
      os << "Event at t=" << ev.m_time << " between PVertex("
         << ev.m_pvertex.first << ":" << ev.m_pvertex.second
         << ") and PVertex(" << ev.m_pother.first << ":" << ev.m_pother.second << ")";
    else if (ev.is_pvertex_to_iedge())
      os << "Event at t=" << ev.m_time << " between PVertex("
         << ev.m_pvertex.first << ":" << ev.m_pvertex.second
         << ") and IEdge" << ev.m_iedge;
    else if (ev.is_pvertex_to_ivertex())
      os << "Event at t=" << ev.m_time << " between PVertex("
         << ev.m_pvertex.first << ":" << ev.m_pvertex.second
         << ") and IVertex(" << ev.m_ivertex << ")";
    else if (ev.is_pvertices_to_ivertex())
      os << "Event at t=" << ev.m_time << " between PVertex("
         << ev.m_pvertex.first << ":" << ev.m_pvertex.second
         << "), PVertex(" << ev.m_pother.first << ":" << ev.m_pother.second
         << " and IVertex(" << ev.m_ivertex << ")";
    else
      os << "Invalid event at t=" << ev.m_time;
    return os;
  }

};


}} // namespace CGAL::KSR_3


#endif // CGAL_KSR_3_EVENT_H
