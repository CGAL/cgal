// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_SURFACE_MESHER_MANIFOLD_H
#define CGAL_SURFACE_MESHER_MANIFOLD_H

#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>

namespace CGAL {

  namespace Surface_mesher {

  template <
    class C2T3,
    class Surface,
    class SurfaceMeshTraits,
    class Criteria,
    typename SMREB = 
      Surface_mesher_regular_edges_base<C2T3, Surface, 
                                        SurfaceMeshTraits, Criteria>
  >
  class Surface_mesher_manifold_base
    : public SMREB 
  {
    public:
      typedef C2T3 C2t3;
      typedef typename C2T3::Triangulation Tr;
      typedef typename Tr::Geom_traits GT;
      typedef typename GT::FT FT;
      typedef typename Tr::Point Point;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
      typedef std::list<Facet> Facets;
      typedef std::list<Vertex_handle> Vertices;
      typedef typename Facets::iterator Facets_iterator;
      typedef typename Vertices::iterator Vertices_iterator;
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

    protected:
      mutable std::set<Vertex_handle> bad_vertices; // @TODO, BEURK: mais
      mutable bool bad_vertices_initialized;        // pourquoi mutable???!!!

    private:
      Facet canonical_facet(const Facet& f) const {
	Cell_handle c = f.first;
	Cell_handle c2 = c->neighbor(f.second);
	return (c2 < c) ? std::make_pair(c2,c2->index(c)) : f;
      }

      // Action to perform on a facet on the boundary of the conflict zone
      void handle_facet_on_boundary_of_conflict_zone (const Facet& f) {
	Facet f1 = canonical_facet(f);
	Cell_handle c = f1.first;
	int i = f1.second;

       	// for each v of f
	for (int j = 0; j < 4; j++) {
	  if (i != j) {
	    Vertex_handle v = c->vertex(j);

	    if(bad_vertices_initialized){
	      //if ( SMREB::c2t3.is_in_complex(v) ) { // no need to test here
		bad_vertices.erase(v);              // il faut tester les
                                                    // facets, avant
		//}
	    }
	  }
	}
      }
      /*
      Facet biggest_incident_facet_in_complex(const Vertex_handle sommet)
      const {
	Graph& g = sommet->get_umbrellas_dual();
	Nodes_map& nodes = g.get_nodes();
	Nodes_map_iterator nit = nodes.begin();
	Facet first_facet = (*nit).first;
	Facet biggest_facet = first_facet;

	for (++nit;
	     nit != nodes.end();
	     ++nit) {
	  Facet current_facet = (*nit).first;
	  // is the current facet bigger than the current biggest one
	  if ( SMREB::c2t3.compute_distance_to_facet_center(current_facet, sommet) >
	       SMREB::c2t3.compute_distance_to_facet_center(biggest_facet, sommet) ) {
	    biggest_facet = current_facet;
	  }
	}
	return biggest_facet;
      }
      */

      Facet biggest_incident_facet_in_complex(const Vertex_handle sommet) const {

	std::list<Facet> facets;
	SMREB::c2t3.incident_facets(sommet, std::back_inserter(facets));

	typename std::list<Facet>::iterator it = facets.begin();
	Facet first_facet = *it;
	Facet biggest_facet = first_facet;

	for (++it;
	     it != facets.end();
	     ++it) {
	  Facet current_facet = *it;
	  // is the current facet bigger than the current biggest one
	  if ( SMREB::compute_distance_to_facet_center(current_facet, sommet) >
	       SMREB::compute_distance_to_facet_center(biggest_facet, sommet) ) {
	    biggest_facet = current_facet;
	  }
	}
	return biggest_facet;
      }
    public:
      Surface_mesher_manifold_base (C2T3& c2t3,
                                    Surface& surface,
                                    SurfaceMeshTraits mesh_traits,
                                    Criteria& criteria)
	: SMREB(c2t3, surface, mesh_traits, criteria),
          bad_vertices_initialized(false)
      {}


    public:

      // Tells whether there remain elements to refine
      bool no_longer_element_to_refine_impl() const {
	if(SMREB::no_longer_element_to_refine_impl()){
	  if(! bad_vertices_initialized){
	    initialize_bad_vertices();
	  }
	  return bad_vertices.empty();
	}
	return false;
      }

      void initialize_bad_vertices() const
      {
	std::cout << "scanning vertices" << std::endl;
	int n = 0;
	for (Finite_vertices_iterator vit = SMREB::tr.finite_vertices_begin();
	     vit != SMREB::tr.finite_vertices_end();
	     ++vit) {
	  if ( (SMREB::c2t3.face_status(vit)  // @TODO: appeler is_regular
		== C2t3::SINGULAR) ) {
	    bad_vertices.insert( vit );
	    ++n;
	  }
	}
	bad_vertices_initialized = true;
	std::cout << "   -> found " << n << " bad vertices\n";
      }

      // Lazy initialization function
      void scan_triangulation_impl() {
	SMREB::scan_triangulation_impl();
	std::cout << "scanning vertices (lazy)" << std::endl;
      }

      // Returns the next element to refine
      Facet get_next_element_impl() {
	if ( !SMREB::no_longer_element_to_refine_impl() ) {
	  return SMREB::get_next_element_impl();
	}
	else {
	  assert(bad_vertices_initialized);
	  Vertex_handle first_bad_vertex = *(bad_vertices.begin());
	  return biggest_incident_facet_in_complex(first_bad_vertex);
	}
      }

      void before_insertion_impl(const Facet&, const Point& s,
				 Zone& zone) {
	for (typename Zone::Facets_iterator fit =
	       zone.boundary_facets.begin(); fit !=
	       zone.boundary_facets.end(); ++fit)
	  if (SMREB::c2t3.is_in_complex(*fit)) {
	    handle_facet_on_boundary_of_conflict_zone (*fit); 
	  }
	SMREB::before_insertion_impl(Facet(), s, zone);
      }

    void after_insertion_impl(const Vertex_handle v) {
      SMREB::after_insertion_impl(v);

      if(bad_vertices_initialized){
	// foreach v' in star of v
	Vertices vertices;
	SMREB::tr.incident_vertices(v, std::back_inserter(vertices));

	// is_regular_or_boundary_for_vertices
	// is used here also incident edges are not known to be
	// REGULAR which may cause some singular vertices to be forgotten
	// This causes no problem because 
	// those SINGULAR incident SINGULAR edges are going to be handled
	for (Vertices_iterator vit = vertices.begin();
	     vit != vertices.end();
	     ++vit)
	  if ( SMREB::c2t3.is_in_complex(*vit)  &&
	       !SMREB::c2t3.is_regular_or_boundary_for_vertices(*vit)) {
	    bad_vertices.insert(*vit);
	  }

	if ( SMREB::c2t3.is_in_complex(v) &&
	     !SMREB::c2t3.is_regular_or_boundary_for_vertices(v)) {
	  bad_vertices.insert(v);
	}
      }
    }
      

      std::string debug_info() const
      {
        std::stringstream s;
        s << SMREB::debug_info() << "," << bad_vertices.size();
        return s.str();
      }
    };  // end Surface_mesher_manifold_base

  template <
    typename C2T3,
    typename Surface,
    typename SurfaceMeshTraits,
    typename Criteria,
    typename SMREB = 
      Surface_mesher_regular_edges_base<C2T3, Surface, 
                                        SurfaceMeshTraits, Criteria>
    >
  class Surface_mesher_manifold
    : public Surface_mesher<C2T3, Surface, SurfaceMeshTraits, Criteria,
                            Surface_mesher_manifold_base<C2T3,
                                                         Surface,
                                                         SurfaceMeshTraits,
                                                         Criteria,
                                                         SMREB> >
  {
    typedef Surface_mesher<C2T3, Surface, SurfaceMeshTraits, Criteria,
                           Surface_mesher_manifold_base<C2T3,
                                                        Surface,
                                                        SurfaceMeshTraits,
                                                        Criteria,
                                                        SMREB> > SM;
  public:
    Surface_mesher_manifold(C2T3& c2t3,
                            Surface& surface,
                            SurfaceMeshTraits mesh_traits,
                            Criteria& criteria)
      : SM(c2t3, surface, mesh_traits, criteria)
    {}
  };  // end Surface_mesher_manifold

  }  // end namespace Surface_mesher

}  // end namespace CGAL


#endif // CGAL_SURFACE_MESHER_MANIFOLD_H

