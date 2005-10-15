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


#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS
#include <CGAL/Voronoi_diagram_2/Identity_Voronoi_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class RT2>
struct Regular_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor
   <typename RT2::Geom_traits::Point_2,RT2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_dual_point_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter
   <typename RT2::Geom_traits::Point_2,RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
  typedef typename RT2::Geom_traits::Point_2           Point_2;
  typedef typename RT2::Geom_traits::Weighted_point_2  Site_2;

  typedef Tag_true                                Has_get_conflicts;
  typedef Tag_true                                Has_insert;
  typedef Tag_true                                Has_remove;
};


//=========================================================================
//=========================================================================


template<class RT2>
class Regular_triangulation_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

//=========================================================================
//=========================================================================

template<class RT2>
class Regular_triangulation_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <RT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_edge_tester_2<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<RT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Regular_triangulation_nearest_site_2<RT2> >
  Base;

  typedef Regular_triangulation_cached_Voronoi_traits_2<RT2>  Self;
};

//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt, class TDS> class Regular_triangulation_2;

template<class Gt, class TDS>
class Identity_Voronoi_traits_2< Regular_triangulation_2<Gt,TDS> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Regular_triangulation_2<Gt,TDS>,
    Regular_triangulation_Voronoi_traits_2
    < Regular_triangulation_2<Gt,TDS> >
  >
{};

template<class Gt, class TDS>
class Identity_Voronoi_traits_2
<  Triangulation_hierarchy_2< Regular_triangulation_2<Gt,TDS> >  >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Triangulation_hierarchy_2< Regular_triangulation_2<Gt,TDS> >,
    Regular_triangulation_Voronoi_traits_2
    <  Triangulation_hierarchy_2< Regular_triangulation_2<Gt,TDS> >  >
  >
{};

#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_VORONOI_TRAITS_2_H
