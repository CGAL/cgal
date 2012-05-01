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

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_ADAPTATION_POLICIES_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_ADAPTATION_POLICIES_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Segment_Delaunay_graph_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Policy_base.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>

#include <CGAL/Identity_policy_2.h>

namespace CGAL {

//=========================================================================
//=========================================================================

template<class SDG2>
struct Segment_Delaunay_graph_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Policy_base
  <SDG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Delaunay_graph_edge_tester_2<SDG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Delaunay_graph_face_tester_2<SDG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename SDG2::Site_2,
						     SDG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_remover<SDG2> >
{
  typedef typename SDG2::Point_2                   Point_2;
  typedef typename SDG2::Site_2                    Site_2;
};


//=========================================================================
//=========================================================================

template<class SDG2>
struct Segment_Delaunay_graph_caching_degeneracy_removal_policy_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Caching_policy_base
  <SDG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Delaunay_graph_edge_tester_2<SDG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Delaunay_graph_face_tester_2<SDG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor,
   CGAL_VORONOI_DIAGRAM_2_INS::Null_functor>
{
  typedef typename SDG2::Point_2                   Point_2;
  typedef typename SDG2::Site_2                    Site_2;
};

//=========================================================================
//=========================================================================

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_ADAPTATION_POLICIES_2_H
