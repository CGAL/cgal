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


#ifndef CGAL_IDENTITY_POLICY_2_H
#define CGAL_IDENTITY_POLICY_2_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Policy_base.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>

namespace CGAL {


template<class DG, class AT>
struct Identity_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Policy_base
  <DG,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_edge_rejector<DG>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_rejector<DG>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename AT::Site_2,DG>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_remover<DG> >
{
  typedef typename AT::Site_2     Site_2;
};


} //namespace CGAL


#endif // CGAL_IDENTITY_POLICY_2_H
