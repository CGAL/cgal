// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_SURFACE_MESHER_MANIFOLD_H
#define CGAL_SURFACE_MESHER_MANIFOLD_H

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3_surface_mesh.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>
#include <CGAL/Node.h>
#include <CGAL/Graph.h>

namespace CGAL {

  namespace Surface_mesher {

  template < typename Tr, 
    typename Surface,
    typename Criteria,
    typename SMREB = Surface_mesher_regular_edges_base 
    <Tr, Surface, Criteria> >
    class Surface_mesher_manifold_base
    : public SMREB {
      
    public:
      typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;
      typedef typename Tr::Point Point;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
      typedef typename C2t3::Face_type Face_type;
      typedef std::list<Facet> Facets;
      typedef std::list<Vertex_handle> Vertices;
      typedef typename Facets::iterator Facets_iterator;
      typedef typename Vertices::iterator Vertices_iterator;
      typedef CGAL::Node<Facet> Node;
      typedef CGAL::Graph<Facet> Graph;
      typedef std::map<Facet, Node*> Nodes_map;
      typedef typename Nodes_map::iterator Nodes_map_iterator;
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

    protected:
      mutable std::set<Vertex_handle> bad_vertices;
      mutable bool bad_vertices_initialized;
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
	      if ( SMREB::c2t3.is_in_complex(v) ) {
		bad_vertices.erase(v);
	      }
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

	std::list<Facet>::iterator it = facets.begin();
	Facet first_facet = *it;
	Facet biggest_facet = first_facet;
	
	for (++it;
	     it != facets.end();
	     ++it) {
	  Facet current_facet = *it;
	  // is the current facet bigger than the current biggest one	  
	  if ( SMREB::c2t3.compute_distance_to_facet_center(current_facet, sommet) > 
	       SMREB::c2t3.compute_distance_to_facet_center(biggest_facet, sommet) ) {
	    biggest_facet = current_facet;
	  }
	}
	return biggest_facet;
      }
    public:
      Surface_mesher_manifold_base (Tr& t, 
						     C2t3& co, 
						     Surface& s, 
						     Criteria& c)
	: SMREB(t, co, s, c), bad_vertices_initialized(false)
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
	  if ( (SMREB::c2t3.face_type(vit)
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
	  handle_facet_on_boundary_of_conflict_zone (*fit);
	
	SMREB::before_insertion_impl(Facet(), s, zone);
      }

      void after_insertion_impl(const Vertex_handle v) {
	SMREB::after_insertion_impl(v);

	if(bad_vertices_initialized){
	  // foreach v' in star of v
	  Vertices vertices;
	  SMREB::c2t3.incident_vertices(v, std::back_inserter(vertices));
	  
	  for (Vertices_iterator vit = vertices.begin();
	       vit != vertices.end();
	       ++vit) {
	    if ( SMREB::c2t3.face_type(*vit) == 
		 SMREB::C2t3::SINGULAR ) {
	      bad_vertices.insert(*vit);
	    }
	  }
	  
	  if ( SMREB::c2t3.is_in_complex(v) ) {
	    if ( SMREB::c2t3.face_type(v) == 
		 SMREB::C2t3::SINGULAR ) {
	      bad_vertices.insert(v);
	    }
	  }
	}
      }
    };  // end Surface_mesher_manifold

  template <typename Tr,
    typename Surface,
    typename Criteria,
    typename Base = Surface_mesher_manifold_base<
      Tr, 
      Surface, 
      Criteria> >
  class Surface_mesher_manifold : 
      public Base, 
      public Mesher_level <
    Tr,
    Surface_mesher_manifold<Tr, Surface, Criteria, Base>,
    typename Tr::Facet,
    Null_mesher_level,
    Triangulation_mesher_level_traits_3<Tr>
  >
  {
  public:
    typedef Surface_mesher_manifold<
      Tr, 
      Surface, 
      Criteria, 
      Base
    > Self;
    typedef Mesher_level <
      Tr,
      Self,
      typename Tr::Facet,
      Null_mesher_level,
       Triangulation_mesher_level_traits_3<Tr>
    > Mesher;
    
    typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;

    using Mesher::scan_triangulation;
    using Mesher::refine;
    using Mesher::is_algorithm_done;
    using Mesher::one_step;
    using Base::check_restricted_delaunay;
    

  protected:
    Null_mesher_level null_mesher_level;
    Null_mesh_visitor null_visitor;
    bool initialized;
    
  public:
    Surface_mesher_manifold(Tr& t, 
					     C2t3& co, 
					     Surface& s, 
					     Criteria& c): 
      Base(t, co, s, c), 
      Mesher(null_mesher_level),
      initialized(false)
      {}
    
    
    // Initialization
    void init(bool debug = false)
      {
	scan_triangulation();
	initialized = true;
	if (debug)
	  check_restricted_delaunay();
      }
    
    void refine_mesh (bool verbose = false, bool debug = false) {
      if(!initialized)
	init (debug);
      
      
      if (!verbose)
	refine (null_visitor);
      else {  
	std::cerr << "Refining...\n";
	int nbsteps = 0;
	std::cerr << "(" << nbsteps << "," 
		  << this->facets_to_refine.size() << ")";
	while (!is_algorithm_done()) {
	  one_step (null_visitor);
	  std::cerr << "\r             \r"
		    << "(" << ++nbsteps << "," 
		    << this->facets_to_refine.size() 
		    << "," 
		    << this->bad_edges.size() 
		    << "," 
		    << this->bad_vertices.size() 
		    << ")";
	}
	std::cerr << "\ndone.\n";
      }
      
      if (debug)
	check_restricted_delaunay();
      
      initialized = false;
    }
    
  };  // end Surface_mesher_manifold
      
  }  // end namespace Surface_mesher
  
}  // end namespace CGAL

  
#endif // CGAL_SURFACE_MESHER_MANIFOLD_H
  
