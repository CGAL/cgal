// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008       GeometryFactory, Sophia Antipolis (France)
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
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_SURFACE_MESHER_MANIFOLD_H
#define CGAL_SURFACE_MESHER_MANIFOLD_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>

namespace CGAL {

  namespace Surface_mesher {

  template <
    class C2T3,
    class Surface,
    class SurfaceMeshTraits,
    class Criteria,
    typename SMMBB // ("SMM_base_base")
  >
  class Surface_mesher_manifold_base
    : public SMMBB
  {
    public:
      typedef C2T3 C2t3;
      typedef typename C2T3::Triangulation Tr;
      typedef typename Tr::Geom_traits GT;
      typedef typename Tr::Point Point;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;

  private:
    // because of the caching, these two members have to be mutable,
    // because they are actually updated in the const method
    // 'no_longer_element_to_refine_impl()'
    typedef std::set<Vertex_handle> Bad_vertices;
    mutable Bad_vertices bad_vertices;
    mutable bool bad_vertices_initialized;

  private:
      typedef typename Tr::Cell_handle Cell_handle;
      typedef std::list<Facet> Facets;
      typedef std::list<Vertex_handle> Vertices;
      typedef typename Facets::iterator Facets_iterator;
      typedef typename Vertices::iterator Vertices_iterator;
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

      Facet canonical_facet(const Facet& f) const {
	const Cell_handle& c = f.first;
	const Cell_handle& c2 = c->neighbor(f.second);
	return (c2 < c) ? std::make_pair(c2,c2->index(c)) : f;
      }

      // Action to perform on a facet on the boundary of the conflict zone
      void 
      before_insertion_handle_facet_on_boundary_of_conflict_zone(const Facet& f)
      {
	const Cell_handle& c = f.first;
	const int i = f.second;

       	// for each v of f
	for (int j = 0; j < 4; j++)
	  if (i != j)
            bad_vertices.erase(c->vertex(j));
      }

      Facet 
      biggest_incident_facet_in_complex(const Vertex_handle sommet) const
      {

	std::list<Facet> facets;
	SMMBB::c2t3.incident_facets(sommet, std::back_inserter(facets));

	typename std::list<Facet>::iterator it = facets.begin();
	Facet biggest_facet = *it;
	CGAL_assertion(it!=facets.end());

	for (++it;
	     it != facets.end();
	     ++it) {
	  Facet current_facet = *it;
	  // is the current facet bigger than the current biggest one
	  if ( SMMBB::compute_distance_to_facet_center(current_facet, sommet) >
	       SMMBB::compute_distance_to_facet_center(biggest_facet, sommet) )
	  {
	    biggest_facet = current_facet;
	  }
	}
	return biggest_facet;
      }

      void initialize_bad_vertices() const
      {
#ifdef CGAL_SURFACE_MESHER_VERBOSE
	std::cerr << "scanning vertices" << std::endl;
#endif
	int n = 0;
	for (Finite_vertices_iterator vit = SMMBB::tr.finite_vertices_begin();
	     vit != SMMBB::tr.finite_vertices_end();
	     ++vit) {
	  if ( (SMMBB::c2t3.face_status(vit)  // @TODO: appeler is_regular
		== C2t3::SINGULAR) ) {
	    bad_vertices.insert( vit );
	    ++n;
	  }
	}
	bad_vertices_initialized = true;
#ifdef CGAL_SURFACE_MESHER_VERBOSE
	std::cerr << "   -> found " << n << " bad vertices\n";
#endif
      }

    public:
      Surface_mesher_manifold_base (C2T3& c2t3,
                                    const Surface& surface,
                                    const SurfaceMeshTraits& mesh_traits,
                                    const Criteria& criteria)
	: SMMBB(c2t3, surface, mesh_traits, criteria),
          bad_vertices_initialized(false)
      {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
        std::cerr << "CONS: Surface_mesher_manifold_base\n";
#endif
      }

    public:

      // Tells whether there remain elements to refine
      bool no_longer_element_to_refine_impl() const {
	if(SMMBB::no_longer_element_to_refine_impl()){
	  if(! bad_vertices_initialized){
	    initialize_bad_vertices();
	  }
	  return bad_vertices.empty();
	}
	return false;
      }

      // Lazy initialization function
      void scan_triangulation_impl() {
	SMMBB::scan_triangulation_impl();
#ifdef CGAL_SURFACE_MESHER_VERBOSE
	std::cerr << "scanning vertices (lazy)" << std::endl;
#endif
      }

      // Returns the next element to refine
      Facet get_next_element_impl() {
	if ( !SMMBB::no_longer_element_to_refine_impl() ) {
	  return SMMBB::get_next_element_impl();
	}
	else {
	  CGAL_assertion(bad_vertices_initialized);
	  Vertex_handle first_bad_vertex = *(bad_vertices.begin());
	  return biggest_incident_facet_in_complex(first_bad_vertex);
	}
      }

      void before_insertion_impl(const Facet& f, const Point& s,
				 Zone& zone) {
        if(bad_vertices_initialized)
        {
          for (typename Zone::Facets_iterator fit =
                 zone.boundary_facets.begin(); fit !=
                 zone.boundary_facets.end(); ++fit)
	    before_insertion_handle_facet_on_boundary_of_conflict_zone (*fit); 
        }
	SMMBB::before_insertion_impl(f, s, zone);
      }

    void after_insertion_impl(const Vertex_handle v) {
      SMMBB::after_insertion_impl(v);

      if(bad_vertices_initialized){
	// foreach v' in star of v
	Vertices vertices;
	SMMBB::tr.incident_vertices(v, std::back_inserter(vertices));

	// is_regular_or_boundary_for_vertices
	// is used here also incident edges are not known to be
	// REGULAR which may cause some singular vertices to be forgotten
	// This causes no problem because 
	// those SINGULAR incident SINGULAR edges are going to be handled
	for (Vertices_iterator vit = vertices.begin();
	     vit != vertices.end();
	     ++vit)
	  if ( SMMBB::c2t3.is_in_complex(*vit)  &&
	       !SMMBB::c2t3.is_regular_or_boundary_for_vertices(*vit)) {
	    bad_vertices.insert(*vit);
	  }

	if ( SMMBB::c2t3.is_in_complex(v) &&
	     !SMMBB::c2t3.is_regular_or_boundary_for_vertices(v)) {
	  bad_vertices.insert(v);
	}
      }
    }
      
    std::string debug_info() const
    {
      std::stringstream s;
      s << SMMBB::debug_info() << ",";
      if(bad_vertices_initialized) 
	s << bad_vertices.size();
      else
	s << "non manifold vertices not initialized";
      return s.str();
    }

    std::string debug_info_header() const
    {
      std::stringstream s;
      s << SMMBB::debug_info_header() << "," << "#bad vertices";
      return s.str();
    }
  };  // end Surface_mesher_manifold_base

  }  // end namespace Surface_mesher

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESHER_MANIFOLD_H
