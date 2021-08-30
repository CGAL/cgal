// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_REGULAR_TRIANGULATION_ADAPTATION_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_ADAPTATION_TRAITS_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Regular_triangulation_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>

namespace CGAL {

template<class RT2>
struct Regular_triangulation_adaptation_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Adaptation_traits_base_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor
   <typename RT2::Geom_traits::Weighted_point_2,RT2,Tag_false>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_Voronoi_point_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;
};


} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_ADAPTATION_TRAITS_2_H
