// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H
#define CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H

#include <CGAL/license/Convex_decomposition_3.h>


#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_external_structure.h>
#include <CGAL/Nef_3/SNC_intersection.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 43
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Nef_>
class External_structure_builder : public Modifier_base<typename Nef_::SNC_and_PL> {

  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename SNC_structure::Items          Items;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>  SNC_intersection;
  typedef CGAL::SNC_external_structure<Items, SNC_structure>
    SNC_external_structure;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>         SM_decorator;
  typedef CGAL::SM_point_locator<SM_decorator>   SM_point_locator;

  typedef typename Base::Segment_3               Segment_3;
  typedef typename Base::Point_3                 Point_3;
  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Sphere_circle           Sphere_circle;
  typedef typename Base::Sphere_segment          Sphere_segment;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::SFace_iterator          SFace_iterator;
  typedef typename Base::SHalfedge_iterator      SHalfedge_iterator;
  typedef typename Base::SFace_cycle_iterator    SFace_cycle_iterator;
  typedef typename Base::SHalfedge_around_sface_circulator
    SHalfedge_around_sface_circulator;

  Halfedge_handle ein;
  Vector_3 dir;

 public:
  External_structure_builder() {}

  void operator()(SNC_and_PL& sncpl) {
    //    CGAL_NEF_TRACEN(43);

    SNC_structure* sncp(sncpl.sncp);
    SNC_point_locator* pl(sncpl.pl);
    SNC_point_locator* old_pl = pl;
    pl = pl->clone();
    sncpl.pl = pl;
    SNC_external_structure C(*sncp,pl);
    C.clear_external_structure();
    C.build_external_structure();
    delete old_pl;
  }
};

} //namespace CGAL
#endif //CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H
