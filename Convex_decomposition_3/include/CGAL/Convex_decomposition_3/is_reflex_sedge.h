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
#ifndef CGAL_CD3_IS_REFLEX_SEDGE_H
#define CGAL_CD3_IS_REFLEX_SEDGE_H

#include <CGAL/license/Convex_decomposition_3.h>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 239
#include <CGAL/Nef_2/debug.h>
#include <CGAL/use.h>

namespace CGAL {

/// @cond SKIP

/*
  int is_reflex_edge(Halfedge_handle e) {
    SHalfedge_around_svertex_circulator
      svc(e->out_sedge()), send(svc);
    int isrse = 0;
    CGAL_For_all(svc, send)
      isrse |= is_reflex_sedge(svc, dir);
    return isrse;
  }
*/

/*
  int is_reflex_sedge(SHalfedge_handle se) {
    return is_reflex_sedge(se, dir);
  }
*/

template<class SNC_structure>
bool is_reflex_sedge_in_any_direction
(typename SNC_structure::SHalfedge_const_handle se) {

  typename SNC_structure::SHalfedge_const_handle
    se2(se->sprev()->twin());
  if(se2 == se) {
    CGAL_NEF_TRACEN("isolated sedge "
                    << se->source()->source()->point()
                    << ": "
                    << se->source()->point()
                    << "->"
                    << se->twin()->source()->point());
    return true;
  }
  typename SNC_structure::Vector_3 vec1
    (se->source()->point() - CGAL::ORIGIN);
  typename SNC_structure::Vector_3 vec2
    (se->circle().orthogonal_vector());
  typename SNC_structure::Sphere_point sp1
    (CGAL::ORIGIN + cross_product(vec2,vec1));
  if(se2->circle().oriented_side(sp1) == ON_POSITIVE_SIDE) {
    CGAL_NEF_TRACEN("reflex sedge "
                    << se->source()->source()->point()
                    << ": "
                    << se->source()->point()
                    << "->"
                    << se->twin()->source()->point());
    return true;
  }
  return false;
}

template<class SNC_structure>
int is_reflex_sedge(typename SNC_structure::SHalfedge_handle se,
                    typename SNC_structure::Sphere_point dir,
                    bool only_small_to_large = true)
{
  typename SNC_structure::Halfedge_handle e = se->source();
  CGAL_NEF_TRACEN("is reflex edge?");
  CGAL_NEF_TRACEN("  e " << e->source()->point()
                  << "->" << e->twin()->source()->point()
                  << " (" << e->point() << ")");

  if(e->point() == dir || e->point() == CGAL::ORIGIN - dir)
    return 0;
  if(only_small_to_large &&
     e->source()->point() > e->twin()->source()->point())
    return 0;


  typename SNC_structure::Sphere_circle cp(e->point(), dir);
  typename SNC_structure::SHalfedge_handle se2 = se->sprev()->twin();
  CGAL_assertion(se2->source() == se->source());

  if(se2 == se) {
    typename SNC_structure::Sphere_segment
      seg(se->source()->point(), se->twin()->source()->point(), se->circle());
    CGAL_NEF_TRACEN("  only one sedge pair " << se->source()->point() <<
                    "->" << se->twin()->source()->point() <<
                    " circle " << se->circle() << " is_short " << seg.is_short());
    if(seg.sphere_circle() == cp)
      return 2;
    if(seg.sphere_circle() == cp.opposite())
      return 1;
    return 3;
  }

  CGAL_NEF_TRACEN(" se1 " << se->circle()
                  << ":" << se->source()->point()
                  << "->" << se->twin()->source()->point());
  CGAL_NEF_TRACEN(" se2 " << se2->circle()
                  << ":" << se2->source()->point()
                  << "->" << se2->twin()->source()->point());
  typename SNC_structure::Vector_3 vec1 = e->point() - CGAL::ORIGIN;
  typename SNC_structure::Vector_3 vec2 = se->circle().orthogonal_vector();
  typename SNC_structure::Sphere_point sp1 = CGAL::ORIGIN + cross_product(vec2,vec1);
  if(se2->circle().oriented_side(sp1) != ON_POSITIVE_SIDE) {
    CGAL_NEF_TRACEN("  too short");
    return 0;
  }

  int result = 0;
  CGAL_NEF_TRACEN(" cp " << cp);
  typename SNC_structure::Vector_3 vec3 = cp.orthogonal_vector();
  typename SNC_structure::Sphere_point sp3 = CGAL::ORIGIN + cross_product(vec3,vec1);
  CGAL_NEF_TRACEN(" sp3 " << sp3);
  CGAL::Oriented_side os1 = se->circle().oriented_side(sp3);
  CGAL_NEF_TRACEN(" os1 " << (int) os1);
  CGAL::Oriented_side os2 = se2->circle().oriented_side(sp3);
  CGAL_NEF_TRACEN(" os2 " << (int) os2);

  if(os1 == ON_POSITIVE_SIDE ||
     os2 == ON_NEGATIVE_SIDE)
    result |= 1;

  if(os1 == ON_NEGATIVE_SIDE ||
     os2 == ON_POSITIVE_SIDE)
    result |= 2;

  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  CGAL_USE_TYPE(Sphere_segment);
  if(os1 == ON_POSITIVE_SIDE &&
     se2->twin()->source()->point() == dir)
    CGAL_assertion(Sphere_segment(se2->source()->point(), se2->twin()->source()->point(), se2->circle()).is_long());

  if(os2 == ON_NEGATIVE_SIDE &&
     se->twin()->source()->point() == dir)
    CGAL_assertion(Sphere_segment(se->source()->point(), se->twin()->source()->point(), se->circle()).is_long());

  if(os1 == ON_NEGATIVE_SIDE &&
     se2->twin()->source()->point() == dir.antipode())
    CGAL_assertion(Sphere_segment(se2->source()->point(), se2->twin()->source()->point(), se2->circle()).is_long());

  if(os2 == ON_POSITIVE_SIDE &&
     se->twin()->source()->point() == dir.antipode())
    CGAL_assertion(Sphere_segment(se->source()->point(), se->twin()->source()->point(), se->circle()).is_long());

  CGAL_NEF_TRACEN("  result " << result);

  return result;
}

/// \endcond

} //namespace CGAL
#endif // CGAL_CD3_IS_REFLEX_SEDGE_H
