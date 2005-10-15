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

#ifndef CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
#define CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Apollonius_graph_nearest_site_2.h>
#include <CGAL/Voronoi_diagram_2/Apollonius_graph_degeneracy_testers.h>
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

template<class AG2>
struct Apollonius_graph_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_Voronoi_traits_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Site_accessor<typename AG2::Site_2,AG2,Tag_true>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_dual_point_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_site_inserter<typename AG2::Site_2,
						     AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
{
  typedef typename AG2::Point_2                   Point_2;
  typedef typename AG2::Site_2                    Site_2;

  typedef Tag_true                                Has_get_conflicts;
  typedef Tag_true                                Has_insert;
  typedef Tag_true                                Has_remove;
};


//=========================================================================
//=========================================================================

template<class AG2>
class Apollonius_graph_cached_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_cached_Voronoi_traits_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
  Base;

  typedef Apollonius_graph_cached_Voronoi_traits_2<AG2>  Self;
};

//=========================================================================
//=========================================================================

template<class AG2>
class Apollonius_graph_ref_counted_Voronoi_traits_2
  : public CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>, 
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
{
 private:
  typedef
  CGAL_VORONOI_DIAGRAM_2_INS::Default_ref_counted_Voronoi_traits_2
  <AG2,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_edge_tester_2<AG2>, 
   CGAL_VORONOI_DIAGRAM_2_INS::Default_face_degeneracy_tester<AG2>,
   CGAL_VORONOI_DIAGRAM_2_INS::Apollonius_graph_nearest_site_2<AG2> >
  Base;

  typedef Apollonius_graph_ref_counted_Voronoi_traits_2<AG2>  Self;
};

//=========================================================================
//=========================================================================

#ifdef VDA_USE_IDENTITY_VORONOI_TRAITS

template<class Gt,class Agds,class LTag> class Apollonius_graph_2;
template<class Gt,class Agds,class LTag> class Apollonius_graph_hierarchy_2;

template<class Gt, class Agds, class LTag>
class Identity_Voronoi_traits_2< Apollonius_graph_2<Gt,Agds,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Apollonius_graph_2<Gt,Agds,LTag>,
    Apollonius_graph_Voronoi_traits_2< Apollonius_graph_2<Gt,Agds,LTag> >
  >
{};

template<class Gt, class Agds, class LTag>
class Identity_Voronoi_traits_2< Apollonius_graph_hierarchy_2<Gt,Agds,LTag> >
  : public CGAL_VORONOI_DIAGRAM_2_INS::Identity_Voronoi_traits_2_base
  < Apollonius_graph_hierarchy_2<Gt,Agds,LTag>,
    Apollonius_graph_Voronoi_traits_2
    < Apollonius_graph_hierarchy_2<Gt,Agds,LTag> >
  >
{};


#endif // VDA_USE_IDENTITY_VORONOI_TRAITS

CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_VORONOI_TRAITS_2_H
