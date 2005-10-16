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
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
#define CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Regular_triangulation_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Regular_triangulation_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_traits_functors.h>

CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class RT2>
struct Regular_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Voronoi_traits_base_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor
   <typename RT2::Geom_traits::Point_2,RT2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_dual_point_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter
   <typename RT2::Geom_traits::Point_2,RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;

  typedef Tag_true                                Has_remove;
};


//=========================================================================
//=========================================================================


template<class RT2>
struct Regular_triangulation_caching_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Caching_Voronoi_traits_base_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor
   <typename RT2::Geom_traits::Point_2,RT2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_dual_point_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter
   <typename RT2::Geom_traits::Point_2,RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;

  typedef Tag_true                                Has_remove;
};

//=========================================================================
//=========================================================================

template<class RT2>
struct Regular_triangulation_identity_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Voronoi_traits_base_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_edge_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Identity_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor
   <typename RT2::Geom_traits::Point_2,RT2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_dual_point_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter
   <typename RT2::Geom_traits::Point_2,RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;

  typedef Tag_true                                Has_remove;
};

//=========================================================================
//=========================================================================

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
