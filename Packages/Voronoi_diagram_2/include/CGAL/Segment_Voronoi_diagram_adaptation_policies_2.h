// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_ADAPTATION_POLICIES_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_ADAPTATION_POLICIES_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Segment_Voronoi_diagram_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Policy_base.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_traits_functors.h>

#include <CGAL/Identity_policy_2.h>

CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class SVD2>
struct Segment_Voronoi_diagram_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Policy_base
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename SVD2::Site_2,
						     SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_remover<SVD2> >
{
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;
};


//=========================================================================
//=========================================================================

template<class SVD2>
struct Segment_Voronoi_diagram_caching_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Caching_policy_base
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor>
{
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;
};

//=========================================================================
//=========================================================================

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_ADAPTATION_POLICIES_2_H
