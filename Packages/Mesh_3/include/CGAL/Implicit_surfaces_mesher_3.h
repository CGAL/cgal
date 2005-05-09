// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IMPLICIT_SURFACES_MESHER_3_H
#define CGAL_IMPLICIT_SURFACES_MESHER_3_H

#include <CGAL/Mesh_3/Refine_tets.h>
#include <CGAL/Chew_4_surfaces/Chew_4_surfaces.h>
#include <CGAL/Chew_4_surfaces/Chew_4_surfaces_visitor.h>

namespace CGAL {

  namespace Mesh_3 
  {
  }

template <
  typename Tr,
  typename Oracle,
  typename Facets_criteria,
  typename Tets_criteria
  >
class Implicit_surfaces_mesher_3
{
public:
  typedef typename Tr::Point Point;
  
  typedef typename
     Chew_4_surfaces::Chew_4_surfaces<Tr,
                                      Oracle,
                                      Facets_criteria> Facets_level;

  typedef typename Mesh_3::Refine_tets<Tr,
                                       Tets_criteria,
                                       Oracle,
                                       Mesh_3::Refine_tets_with_oracle_base<Tr,
                                        Tets_criteria, Oracle>, Facets_level>
                                                     Tets_level;

  typedef typename Mesh_3::tets::Refine_facets_visitor<Tr,
     Tets_level> Tets_facets_visitor;
  typedef typename Chew_4_surfaces::Visitor<Tr, Facets_level>
                 Chew_facets_visitor;
  
  typedef Combine_mesh_visitor<Chew_facets_visitor, Tets_facets_visitor>
                 Facets_visitor;

  typedef Null_mesh_visitor_level<Facets_visitor> Tets_visitor;

  typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;

private:
  Null_mesher_level null_mesher_level;
  Null_mesh_visitor null_visitor;

  C2t3 c2t3;
  Oracle& oracle;
  Facets_level facets;
  Tets_level tets;

  Chew_facets_visitor chew_facets_visitor;
  Tets_facets_visitor tets_facets_visitor;
  Tets_visitor tets_visitor;

  bool initialized;

public:
  Implicit_surfaces_mesher_3(Tr& t, Oracle& o,
                             Facets_criteria& c,
                             Tets_criteria tets_crit)
    : c2t3(t), oracle(o), 
      facets(t, c2t3, oracle, c), tets(t, tets_crit, oracle, facets),
      chew_facets_visitor(&facets),
      tets_facets_visitor(&tets),
      tets_visitor(Facets_visitor(&chew_facets_visitor,
				  &tets_facets_visitor)),
      initialized(false)
  {}

  void init()
  {
    Tr& tr = tets.triangulation_ref_impl();

    typedef typename Tr::Geom_traits Geom_traits;
    typename Geom_traits::Construct_circumcenter_3 circumcenter = 
      tr.geom_traits().construct_circumcenter_3_object();

    for(typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
	cit != tr.finite_cells_end();
	++cit)
      {
	const Point& p = cit->vertex(0)->point();
	const Point& q = cit->vertex(1)->point();
	const Point& r = cit->vertex(2)->point();
	const Point& s = cit->vertex(3)->point();

	cit->set_in_domain(oracle.surf_equation(circumcenter(p,q,r,s))<0.);
      }
    
    facets.scan_triangulation();

    for(typename Tr::Finite_vertices_iterator vit = 
      tr.finite_vertices_begin();
      vit != tr.finite_vertices_end();
      ++vit)
      vit->info()=true;
    std::cerr << "Restore infos.\n";

    tets.scan_triangulation();
    initialized = true;
  }

  void refine_mesh()
  {
    if(!initialized)
      init();
    tets.refine(tets_visitor);
  }

  void refine_surface()
  {
    if(!initialized)
      init();
    facets.refine(tets_visitor.previous_level());
  }

  void step_by_step()
  {
    if(!initialized)
      init();
    tets.try_to_insert_one_point(tets_visitor);
  }

  bool done()
  {
    return tets.no_longer_element_to_refine();
    
  }
}; // end Implicit_surfaces_mesher_3

} // end namespace CGAL

#endif // CGAL_IMPLICIT_SURFACES_MESHER_3_H
