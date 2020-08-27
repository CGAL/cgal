// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_2_offset_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_triangulation_remove_traits_2.h>

namespace CGAL {

template <class Gt_,
          class Off_ = typename CGAL::Periodic_2_offset_2>
class Periodic_2_Delaunay_triangulation_remove_traits_2
  : public Periodic_2_triangulation_remove_traits_2<Gt_>
{
  typedef Periodic_2_Delaunay_triangulation_remove_traits_2<Gt_, Off_>  Self;
  typedef Periodic_2_triangulation_remove_traits_2<Gt_>                 Base;

public:
  typedef Gt_                                                           Geom_traits;
  typedef Off_                                                          Offset;

  typedef typename Geom_traits::RT                                      RT;
  typedef typename Geom_traits::FT                                      FT;
  typedef std::pair<typename Geom_traits::Point_2, Offset>              Point_2;

  // not allowing a default value for `gt` because we need to have
  // an initialized domain in `gt`
  Periodic_2_Delaunay_triangulation_remove_traits_2(const Geom_traits& gt) : Base(gt) { }

};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_REMOVE_TRAITS_2_H
