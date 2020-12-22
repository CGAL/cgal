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

#ifndef CGAL_APOLLONIUS_GRAPH_ADAPTATION_POLICIES_2_H
#define CGAL_APOLLONIUS_GRAPH_ADAPTATION_POLICIES_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Apollonius_graph_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Policy_base.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>

#include <CGAL/Identity_policy_2.h>

namespace CGAL {

//=========================================================================
//=========================================================================


template<class AG2>
struct Apollonius_graph_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Policy_base
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_rejector<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename AG2::Site_2,AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_remover<AG2> >
{
  typedef typename AG2::Site_2                    Site_2;
};


//=========================================================================
//=========================================================================

template<class AG2>
struct Apollonius_graph_caching_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Caching_policy_base
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_rejector<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename AG2::Site_2,AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_remover<AG2> >
{
  typedef typename AG2::Site_2                    Site_2;
};

//=========================================================================
//=========================================================================

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_ADAPTATION_POLICIES_2_H
