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

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Segment_Voronoi_diagram_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Segment_Voronoi_diagram_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Default_Voronoi_traits_2.h>
#include <CGAL/Voronoi_diagram_2/Site_accessors.h>
#include <CGAL/Voronoi_diagram_2/Construct_dual_points.h>

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS
#include <CGAL/Voronoi_diagram_2/Identity_Voronoi_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE


//=========================================================================
//=========================================================================

template<class SVD2>
class Segment_Voronoi_diagram_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
  Base;

  typedef Segment_Voronoi_diagram_Voronoi_traits_2<SVD2>  Self;

 public:
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;
  typedef typename SVD2::Vertex_handle             Vertex_handle;
  typedef typename SVD2::Face_handle               Face_handle;

  typedef Tag_false                                Has_get_conflicts;
  typedef Tag_false                                Has_insert;
  typedef Tag_false                                Has_remove;

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<Site_2,Vertex_handle,Tag_false>
  Get_site_2;

  Get_site_2 get_site_2_object() const { return Get_site_2(); }

  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_dual_point_2<SVD2>
  Get_point_2;

  Get_point_2 get_point_2_object() const { return Get_point_2(); }
};



//=========================================================================
//=========================================================================

template<class SVD2>
class Segment_Voronoi_diagram_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;
};

//=========================================================================
//=========================================================================


template<class SVD2>
class Segment_Voronoi_diagram_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <SVD2,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_edge_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_face_tester_2<SVD2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Segment_Voronoi_diagram_nearest_site_2<SVD2> >
  Base;

  typedef Segment_Voronoi_diagram_cached_Voronoi_traits_2<SVD2>  Self;
};

//=========================================================================
//=========================================================================


//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt, class STag, class DS, class LTag >
class Segment_Voronoi_diagram_hierarchy_2;

template<class Gt, class DS, class LTag>
class Segment_Voronoi_diagram_2;


template<class Gt, class DS, class LTag>
class Identity_Voronoi_traits_2< Segment_Voronoi_diagram_2<Gt,DS,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Segment_Voronoi_diagram_2<Gt,DS,LTag>,
    Segment_Voronoi_diagram_Voronoi_traits_2
    < Segment_Voronoi_diagram_2<Gt,DS,LTag> >
  >
{};

template<class Gt, class STag, class DS, class LTag>
class Identity_Voronoi_traits_2<
  Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag>,
    Segment_Voronoi_diagram_Voronoi_traits_2
    < Segment_Voronoi_diagram_hierarchy_2<Gt,STag,DS,LTag> >
  >
{};

#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VORONOI_TRAITS_2_H
