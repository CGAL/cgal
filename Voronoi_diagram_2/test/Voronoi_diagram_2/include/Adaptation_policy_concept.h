// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_ADAPTATION_POLICY_CONCEPT_H
#define CGAL_ADAPTATION_POLICY_CONCEPT_H 1

#include <CGAL/basic.h>
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
