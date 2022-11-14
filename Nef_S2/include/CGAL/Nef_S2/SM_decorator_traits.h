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
// Author(s)     : Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_SM_DECORATOR_TRAITS_H
#define CGAL_NEF_SM_DECORATOR_TRAITS_H

#include <CGAL/license/Nef_S2.h>


namespace CGAL {

template <class Refs_>
class SM_decorator_traits {
  typedef Refs_ Refs;
 public:
  typedef typename Refs::SVertex_handle SVertex_handle;
  typedef typename Refs::SHalfedge_handle SHalfedge_handle;
  typedef typename Refs::SHalfloop_handle SHalfloop_handle;
  typedef typename Refs::SFace_handle SFace_handle;

  typedef typename Refs::SVertex_iterator SVertex_iterator;
  typedef typename Refs::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Refs::SHalfloop_iterator SHalfloop_iterator;
  typedef typename Refs::SFace_iterator SFace_iterator;

  typedef typename Refs::SVertex_const_handle SVertex_const_handle;
  typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Refs::SFace_const_handle SFace_const_handle;

  typedef typename Refs::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Refs::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Refs::SHalfloop_const_iterator SHalfloop_const_iterator;
  typedef typename Refs::SFace_const_iterator SFace_const_iterator;

  typedef typename Refs::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;
  typedef typename Refs::SHalfedge_around_sface_circulator
    SHalfedge_around_sface_circulator;
  typedef typename Refs::SFace_cycle_iterator SFace_cycle_iterator;
};

template <class Refs_>
class SM_decorator_const_traits {
  typedef Refs_ Refs;
 public:
  typedef typename Refs::SVertex_const_handle SVertex_handle;
  typedef typename Refs::SHalfedge_const_handle SHalfedge_handle;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_handle;
  typedef typename Refs::SFace_const_handle SFace_handle;

  typedef typename Refs::SVertex_const_iterator SVertex_iterator;
  typedef typename Refs::SHalfedge_const_iterator SHalfedge_iterator;
  typedef typename Refs::SHalfloop_const_iterator SHalfloop_iterator;
  typedef typename Refs::SFace_const_iterator SFace_iterator;

  typedef typename Refs::SHalfedge_around_svertex_const_circulator
    SHalfedge_around_svertex_circulator;
  typedef typename Refs::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_circulator;
  typedef typename Refs::SFace_cycle_const_iterator SFace_cycle_iterator;

  typedef typename Refs::SVertex_const_handle SVertex_const_handle;
  typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Refs::SFace_const_handle SFace_const_handle;

  typedef typename Refs::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Refs::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Refs::SHalfloop_const_iterator SHalfloop_const_iterator;
  typedef typename Refs::SFace_const_iterator SFace_const_iterator;

  typedef typename Refs::SHalfedge_around_svertex_const_circulator
    SHalfedge_around_svertex_const_circulator;
  typedef typename Refs::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;
  typedef typename Refs::SFace_cycle_const_iterator SFace_cycle_const_iterator;
};

} //namespace CGAL
#endif // CGAL_NEF_SM_DECORATOR_TRAITS_H
