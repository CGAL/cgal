// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CONVEX_DECOMPOSITION_3_H
#define CGAL_CONVEX_DECOMPOSITION_3_H

#include <CGAL/license/Convex_decomposition_3.h>


#include <CGAL/Convex_decomposition_3/Single_wall_creator.h>
#include <CGAL/Convex_decomposition_3/Single_wall_creator2.h>
#include <CGAL/Convex_decomposition_3/Reflex_edge_searcher.h>
#include <CGAL/Convex_decomposition_3/Reflex_vertex_searcher.h>
#include <CGAL/Convex_decomposition_3/YVertical_wall_builder.h>
#include <CGAL/Convex_decomposition_3/Ray_hit_generator2.h>
#include <CGAL/Convex_decomposition_3/External_structure_builder.h>
#include <CGAL/Convex_decomposition_3/SFace_separator.h>
#include <CGAL/Convex_decomposition_3/Edge_sorter.h>
#include <CGAL/Convex_decomposition_3/is_reflex_sedge.h>

/*! 
  \file convex_decomposition_3.h 
*/


/// The CGAL namespace.
namespace CGAL {

/*!
\ingroup PkgConvexDecomposition3

The function `convex_decomposition_3()` inserts additional facets 
into the given `Nef_polyhedron_3` `N`, such that each bounded 
marked volume (the outer volume is unbounded) is subdivided into convex 
pieces. The modified polyhedron represents a decomposition into 
\f$ O(r^2)\f$ convex pieces, where \f$ r\f$ is the number of edges that have two 
adjacent facets that span an angle of more than 180 degrees with 
respect to the interior of the polyhedron. 

The worst-case running time of our implementation is 
\f$ O(n^2r^4\sqrt[3]{nr^2}\log{(nr)})\f$, where \f$ n\f$ is the complexity of 
the polyhedron (the complexity of a `Nef_polyhedron_3` is the sum 
of its `Vertices`, `Halfedges` and `SHalfedges`) and \f$ r\f$ 
is the number of reflex edges. 

\pre The polyhedron `N` is bounded. Otherwise, the outer volume is ignored. 

\post If the polyhedron `N` is non-convex, it is modified to represent the 
convex decomposition. If `N` is convex, it is not modified. 

\sa `CGAL::Nef_polyhedron_3<Traits>` 

*/
template<typename Nef_polyhedron>
void convex_decomposition_3(Nef_polyhedron& N) 
{
  typedef typename Nef_polyhedron::SNC_structure  SNC_structure;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename Nef_polyhedron::Point_3            Point_3;
  typedef typename Nef_polyhedron::Vector_3           Vector_3;
  typedef typename Nef_polyhedron::Sphere_point   Sphere_point;
  typedef typename Nef_polyhedron::FT FT;

  typedef typename CGAL::Single_wall_creator<Nef_polyhedron>  Single_wall;
  typedef typename CGAL::YVertical_wall_builder<Nef_polyhedron> YVertical_wall_builder;
  typedef typename CGAL::Reflex_vertex_searcher<Nef_polyhedron>  Reflex_vertex_searcher;
  typedef typename CGAL::Ray_hit_generator2<Nef_polyhedron> Ray_hit2;
  typedef typename CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;
  typedef typename CGAL::SFace_separator<Nef_polyhedron> SFace_separator;

  typedef Compare_halfedges_in_reflex_edge_sorter<Halfedge_handle, std::less<Point_3> >
	  Less_edges;
  typedef Compare_halfedges_in_reflex_edge_sorter<Halfedge_handle, std::greater<Point_3> >
	  Greater_edges;

  typedef typename std::multiset<Halfedge_handle, Less_edges> Negatively_sorted_set;
  typedef typename std::multiset<Halfedge_handle, Greater_edges> Positively_sorted_set;
  typedef typename Positively_sorted_set::const_iterator  Positive_reflex_edge_iterator;
  typedef typename Negatively_sorted_set::const_iterator  Negative_reflex_edge_iterator;

  typedef typename CGAL::Reflex_edge_searcher<Nef_polyhedron, Positively_sorted_set, Negatively_sorted_set>
	  Reflex_edge_searcher;

  typedef typename CGAL::Edge_sorter<Nef_polyhedron, std::less<FT>, Negatively_sorted_set> Edge_sorter;
  typedef typename CGAL::Edge_sorter<Nef_polyhedron, std::greater<FT>, Positively_sorted_set> Edge_sorter2;

  External_structure_builder esb;
  SFace_separator sf_sep;
  N.delegate(sf_sep,false, false);

  Reflex_edge_searcher res(Sphere_point(1,0,0));
  N.delegate(res,false,false);
  
  Edge_sorter es(res.get_negative_redges());
  N.delegate(es);
  
  Negative_reflex_edge_iterator nrei;
  for(nrei=res.negative_redges_begin(); nrei!=res.negative_redges_end(); ++nrei) {
    Halfedge_handle e = (*nrei);
    
    Single_wall W(e,Vector_3(-1,0,0));
    if(!W.need_to_create_wall()) continue;
    
    Reflex_vertex_searcher rvs(Sphere_point(1,0,0));
    if(rvs.need_to_shoot(e, true)) {
      Ray_hit2 rh2a(Vector_3(-1,0,0), e->source());
      N.delegate(rh2a);
    }
    if(rvs.need_to_shoot(e->twin(), true)) {
      Ray_hit2 rh2a(Vector_3(-1,0,0), e->twin()->source());
      N.delegate(rh2a);
    }  
  }
  
  // int i=0;
  for(nrei=res.negative_redges_begin(); nrei!=res.negative_redges_end(); ++nrei) {
    Halfedge_handle e = (*nrei);
    if(e->point().hx() > 0)
      e = e->twin();
    Single_wall W(e,Vector_3(-1,0,0));
    if(!W.need_to_create_wall()) continue;    
    N.delegate(W);
  }
  
  N.delegate(esb);
  N.delegate(res, false, false);
  
  CGAL_assertion(N.is_valid(0,0));
  
  Reflex_edge_searcher& res2 = res;
  
  Edge_sorter2 es2(res2.get_positive_redges());
  N.delegate(es2);
  
  Positive_reflex_edge_iterator prei;
  for(prei=res2.positive_redges_begin(); prei!=res2.positive_redges_end(); ++prei) {
    Halfedge_handle e = (*prei);
    
    CGAL_assertion(e->source()->point() >
		   e->twin()->source()->point());
    Single_wall W(e,Vector_3(1,0,0));
    if(!W.need_to_create_wall()) continue;
    
    Reflex_vertex_searcher rvs(Sphere_point(1,0,0));
    if(rvs.need_to_shoot(e, false)) {
      Ray_hit2 rh2a(Vector_3(1,0,0), e->source());
      N.delegate(rh2a);
    }
    if(rvs.need_to_shoot(e->twin(), false)) {
      Ray_hit2 rh2a(Vector_3(1,0,0), e->twin()->source());
      N.delegate(rh2a);
    }
  }
  
  // i=0;
  for(prei=res2.positive_redges_begin(); prei!=res2.positive_redges_end(); ++prei) {
    Halfedge_handle e = (*prei);
    Single_wall W(e,Vector_3(1,0,0));
    if(!W.need_to_create_wall()) continue;
    N.delegate(W);
  }
  
  N.delegate(esb);    
  CGAL_assertion(N.is_valid(0,0));
  
  YVertical_wall_builder Y;
  N.delegate(Y,false,false);
  
  N.delegate(esb);

  CGAL_assertion_code(typename Nef_polyhedron::SHalfedge_const_iterator cse);
  CGAL_assertion_code(CGAL_forall_shalfedges(cse, N)
    if(cse->incident_sface()->mark())
      CGAL_assertion(!CGAL::is_reflex_sedge_in_any_direction<Nef_polyhedron>(cse)));
}

} //namespace CGAL
#endif // CGAL_CONVEX_DECOMPOSITION_3_H
