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

#ifndef CGAL_ADAPTATION_POLICY_CONCEPT_H
#define CGAL_ADAPTATION_POLICY_CONCEPT_H 1

#include <CGAL/Voronoi_diagram_2/Policy_base.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>

namespace CGAL {

//=========================================================================

template<class DG, class AT>
struct Adaptation_policy_concept
  : public CGAL_VORONOI_DIAGRAM_2_INS::Policy_base
  <DG,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_edge_rejector<DG>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_rejector<DG>,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor>
{
  typedef typename AT::Site_2   Site_2;
};


} //namespace CGAL


#endif // CGAL_ADAPTATION_POLICY_CONCEPT_H
