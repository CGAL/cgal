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

#ifndef CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Delaunay_triangulation_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Delaunay_triangulation_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>


CGAL_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================


template<class DT2>
class Delaunay_triangulation_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
{
 private:

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
  Base;

  typedef Delaunay_triangulation_Voronoi_traits_2<DT2>  Self;

 public:
  typedef typename DT2::Geom_traits::Point_2      Point_2;
  typedef Point_2                                 Site_2;
  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Face_handle              Face_handle;

  typedef Tag_true                                Has_get_conflicts;
  typedef Tag_true                                Has_insert;
  typedef Tag_true                                Has_remove;

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Point_accessor<Point_2,Vertex_handle,Tag_true>
  Get_site_2;

  Get_site_2 get_site_2_object() const { return Get_site_2(); }

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_dual_point_2<DT2>
  Get_point_2;

  Get_point_2 get_point_2_object() const { return Get_point_2(); }
};


//=========================================================================
//=========================================================================

template<class DT2>
class Delaunay_triangulation_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
{
private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
  Base;

  typedef Delaunay_triangulation_cached_Voronoi_traits_2<DT2>  Self;
};

//=========================================================================
//=========================================================================

template<class DT2>
class Delaunay_triangulation_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
{
 private:

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <DT2,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_edge_tester_2<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<DT2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Delaunay_triangulation_nearest_site_2<DT2> >
  Base;

  typedef Delaunay_triangulation_ref_counted_Voronoi_traits_2<DT2>  Self;
};

//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt,class TDS> class Delaunay_triangulation_2;

template<class Gt, class TDS>
class Identity_Voronoi_traits_2< Delaunay_triangulation_2<Gt,TDS> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Delaunay_triangulation_2<Gt,TDS>,
    Delaunay_triangulation_Voronoi_traits_2
    < Delaunay_triangulation_2<Gt,TDS> >
  >
{};

template<class Gt, class TDS>
class Identity_Voronoi_traits_2
<  Triangulation_hierarchy_2< Delaunay_triangulation_2<Gt,TDS> >  >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Triangulation_hierarchy_2< Delaunay_triangulation_2<Gt,TDS> >,
    Delaunay_triangulation_Voronoi_traits_2
    <  Triangulation_hierarchy_2< Delaunay_triangulation_2<Gt,TDS> >  >
  >
{};

#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_VORONOI_TRAITS_2_H
