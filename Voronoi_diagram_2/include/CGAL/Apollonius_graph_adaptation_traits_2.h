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

#ifndef CGAL_APOLLONIUS_GRAPH_ADAPTATION_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_ADAPTATION_TRAITS_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Apollonius_graph_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_base_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>


namespace CGAL {

template<class AG2>
struct Apollonius_graph_adaptation_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Adaptation_traits_base_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<typename AG2::Site_2,AG2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_Voronoi_point_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
{
  typedef typename AG2::Point_2                   Point_2;
  typedef typename AG2::Site_2                    Site_2;
};


} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_ADAPTATION_TRAITS_2_H
