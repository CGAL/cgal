// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel        <seel@mpi-sb.mpg.de>
//                 Miguel Granados     <granados@mpi-sb.mpg.de>
//                 Susan Hert          <hert@mpi-sb.mpg.de>
//                 Lutz Kettner        <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_ITEMS_H
#define CGAL_SNC_ITEMS_H

#include <CGAL/license/Nef_3.h>

#include <CGAL/Nef_3/Vertex.h>
#include <CGAL/Nef_3/Halfedge.h>
#include <CGAL/Nef_3/Halffacet.h>
#include <CGAL/Nef_3/Volume.h>
#include <CGAL/Nef_3/SHalfedge.h>
#include <CGAL/Nef_3/SHalfloop.h>
#include <CGAL/Nef_3/SFace.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename K, typename I, typename M> class SNC_sphere_map;
template <typename R> class SM_decorator;

class SNC_items {
 public:
  template <class Refs> class Vertex    :    public Vertex_base<Refs> {};
  template <class Refs> class SVertex   :    public Halfedge_base<Refs> {};
  template <class Refs> class Halffacet :    public Halffacet_base<Refs> {};
  template <class Refs> class Volume    :    public Volume_base<Refs> {};
  template <class Refs> class SHalfedge :    public SHalfedge_base<Refs> {};
  template <class Refs> class SHalfloop :    public SHalfloop_base<Refs> {};
  template <class Refs> class SFace     :    public SFace_base<Refs> {};
}; // SNC_items



} //namespace CGAL
#endif //CGAL_SNC_ITEMS_H
