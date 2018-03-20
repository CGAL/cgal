// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008-2009  GeometryFactory (France)
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
// Author(s) : Steve Oudot,
//             David Rey,
//             Mariette Yvinec,
//             Laurent Rineau,
//             Andreas Fabri



#ifndef CGAL_SURFACE_MESHER_SURFACE_MESHER_H
#define CGAL_SURFACE_MESHER_SURFACE_MESHER_H

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesher_level.h>
#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>
#include <CGAL/Double_map.h>
#include <CGAL/Timer.h>
#include <CGAL/Object.h>
#include <list>
#include <string>
#include <sstream>
#include <boost/format.hpp>

#include <CGAL/Surface_mesher/Verbose_flag.h>
#include <CGAL/Surface_mesher/Types_generators.h>
#include <CGAL/Surface_mesher/Profile_counter.h>
#include <CGAL/Surface_mesher/Profile_timer.h>

namespace CGAL {

  namespace Surface_mesher {

  // NB: by convention, the priority queue is sorted with respect to the
  // first element of the list of criteria

  template <
    class C2T3,
    class Surface_,
    class SurfaceMeshTraits,
    class Criteria_
    >
  class Surface_mesher_base
    : public Triangulation_mesher_level_traits_3<typename C2T3::Triangulation>
  {
  public:
    typedef C2T3 Complex_2_in_triangulation_3;
    typedef Surface_ Surface;
    typedef SurfaceMeshTraits Surface_mesh_traits;
    typedef Criteria_ Criteria;

    typedef typename C2T3::Triangulation Tr;

    typedef typename Tr::Point Point;
    typedef typename Surface_mesh_traits::Intersection_point Intersection_point;

    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Edge Edge;
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Cell_handle Cell_handle;

    typedef typename Tr::Geom_traits GT;

    typedef Triangulation_mesher_level_traits_3<Tr> Triangulation_traits;
    typedef typename Triangulation_traits::Zone Zone;

    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;

    typedef typename Criteria::Quality Quality;

    typedef Double_map<Facet, Quality> Bad_facets;

    // Constructor
    Surface_mesher_base (C2T3& co, 
                         const Surface& s, 
                         const Surface_mesh_traits& mesh_traits,
                         const Criteria& c) :
      Triangulation_mesher_level_traits_3<Tr>(co.triangulation()),
      c2t3(co),
      tr(co.triangulation()),
      surf(s),
      meshtraits(mesh_traits),
      criteria(c)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Surface_mesher_base\n";
#endif
    }

  protected:
    C2T3& c2t3;
    Tr& tr;     // Associated triangulation reference
    const Surface& surf;  // Surface
    const Surface_mesh_traits& meshtraits; // Surface mesh traits
    const Criteria& criteria;  // Meshing criteria
    Bad_facets facets_to_refine;  // Set of facets to refine

  public:

    // Helper functions
    Facet mirror_facet(const Facet& f) const
    {
      return tr.mirror_facet(f);
    }

    static void set_facet_visited(Facet f)
    {
      f.first->set_facet_visited(f.second);
    }

    static void set_facet_surface_center(Facet f, const Point& p)
    {
      f.first->set_facet_surface_center(f.second, p);
    }

    static const Point& get_facet_surface_center(const Facet& f)
    {
      return f.first->get_facet_surface_center(f.second);
    }

    static void reset_visited(Facet f)
    {
      f.first->reset_visited(f.second);
    }

    static bool is_facet_visited(const Facet& f)
    {
      return f.first->is_facet_visited(f.second);
    }

    // Remains unchanged
    Tr& triangulation_ref_impl()
      {
	return tr;
      }

    // Remains unchanged
    const Tr& triangulation_ref_impl() const
      {
	return tr;
      }

    // Initialization function
    void scan_triangulation_impl() {
      Point center;
      // We test only the finite Delaunay facets

#ifdef CGAL_SURFACE_MESHER_VERBOSE
      std::cout << "scanning facets..." << std::endl;
#endif
      for (Finite_facets_iterator fit = tr.finite_facets_begin(); fit !=
	     tr.finite_facets_end(); ++fit) {
	if (tr.dimension() == 3) {
          new_facet<true>(*fit);
          // see definition of
          // template <bool> new_facet()
	}
	else {
	  CGAL_assertion (tr.dimension() == 2);
	  if (is_facet_on_surface(*fit, center)) {
	    //Cell_handle c;
	    c2t3.set_in_complex(*fit);
	    set_facet_surface_center((*fit),center);
	    //c=(*fit).first;
	    //c->set_facet_on_surface((*fit).second,true);
	    //c->set_facet_surface_center((*fit).second,center);

            Quality a_r;
	    if (criteria.is_bad (*fit, a_r)) {
#ifdef CGAL_SURFACE_MESHER_TAG_BAD
              fit->first->set_bad(fit->second);
              const Facet mirror_facet = tr.mirror_facet(*fit);
              mirror_facet.first->set_bad(mirror_facet.second);
#endif // CGAL_SURFACE_MESHER_TAG_BAD
	      facets_to_refine.insert(*fit,a_r);
	    }
	  }
	  else
	    c2t3.remove_from_complex(*fit);
	  //	    (*fit).first->set_facet_on_surface((*fit).second,false);
	}
      }
    } // scan_triangulation_impl end

    // Tells whether there remain elements to refine
    bool no_longer_element_to_refine_impl() const {
      return facets_to_refine.empty();
    }

    // Returns the next element to refine
    Facet get_next_element_impl() {
      return facets_to_refine.front()->second;
    }

    // deletes the next element from the set of elements to refine
    // NB: it is useless here, since the update of the restricted
    // Delaunay triangulation automatically deletes the next element
    void pop_next_element_impl() {
//       facets_to_refine.pop_front();
    }

    // From the element to refine, gets the point to insert
    Point refinement_point_impl(const Facet& f) const
      {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
        std::cerr << "point from Surface_mesher: ";
#endif
	CGAL_assertion (c2t3.face_status(f) == C2T3::REGULAR);
	//CGAL_assertion (f.first->is_facet_on_surface (f.second));
	return get_facet_surface_center (f);
	//return f.first->get_facet_surface_center (f.second);
      }

    // Useless here
    void before_conflicts_impl(const Facet&, const Point&
#ifdef CGAL_SURFACE_MESHER_DEBUG_BEFORE_CONFLICTS
			                                s)
    {
      std::cerr << "Refine_facets: before conflicts of " << s << " ";
#else
                                                        )
    {
#endif
    }

  // Useless here
    Mesher_level_conflict_status private_test_point_conflict_impl(const Point& p,
								  Zone& zone)
    {
      if( zone.locate_type == Tr::VERTEX )
      {
        std::stringstream sstr;
        sstr << "(" << p << ") is already inserted on surface.\n";
        CGAL_error_msg(sstr.str().c_str());
	return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
      }
      else
	return NO_CONFLICT;
    }

    // Useless here
    void after_no_insertion_impl(const Facet&, const Point&, const Zone&)
    {
    }

    ///////////////////////////////////////////////////////////////////////////
    // Tests if the point p encroaches one facet of the restricted Delaunay
    // triangulation.

    Mesher_level_conflict_status
    test_point_conflict_from_superior_impl(const Point& p,
					   Zone& zone)
      {
	for (typename Zone::Facets_iterator fit =
	       zone.internal_facets.begin();
	     fit != zone.internal_facets.end();
	     ++fit)
	  if( test_if_facet_is_encroached(*fit, p) )
	    return CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;

	for (typename Zone::Facets_iterator fit =
	       zone.boundary_facets.begin();
	     fit != zone.boundary_facets.end();
	     ++fit)
	  if( test_if_facet_is_encroached(*fit, p) )
	    return CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;

	return NO_CONFLICT;
      }

    bool test_if_facet_is_encroached(const Facet f, const Point& p)
    {
      if( tr.is_infinite(f.first) )
	return false;

      if (c2t3.face_status(f) == C2T3::REGULAR)
	{
	  const Cell_handle& c = f.first;
	  const int& index = f.second;

	  typename GT::Compute_squared_distance_3 distance =
	    tr.geom_traits().compute_squared_distance_3_object();

	  // test Delaunay surfacique
	  Point center = get_facet_surface_center(f);

          for(bool exit = false; ; exit = true)
          {
            // this for loop is a trick to pass in the following "if" once
            // with center="surface center", and once with
            // center="circumcenter"

            if( distance(center, p) <
                distance(center, c->vertex((index+1)&3)->point()) )
	    {
              Quality q;
              criteria.is_bad(f, q); // to get q (passed as reference)
              Facet other_side = mirror_facet(f);
#ifdef CGAL_SURFACE_MESHER_TAG_BAD
              f.first->set_bad(f.second);
              other_side.first->set_bad(other_side.second);
#endif // CGAL_SURFACE_MESHER_TAG_BAD
              if(f.first < other_side.first)
                facets_to_refine.insert(f, q);
              else
                facets_to_refine.insert(other_side, q);
	      return true;
	    }

//             if(exit)
              return false;

            // test Gabriel
            center = tr.geom_traits().construct_circumcenter_3_object()
              (c->vertex((index+1)&3)->point(),
               c->vertex((index+2)&3)->point(),
               c->vertex((index+3)&3)->point());
          }
	}
      return false;
    }


  /* returns the conflicts zone */
    Zone conflicts_zone_impl(const Point& p, Facet f) const {
	Zone zone;

       	// TODO may be avoid the locate here
	zone.cell =
	  tr.locate(p, zone.locate_type, zone.i, zone.j, f.first);

        if(zone.locate_type != Tr::VERTEX)
          tr.find_conflicts(p, zone.cell,
                            std::back_inserter(zone.boundary_facets),
                            std::back_inserter(zone.cells),
                            std::back_inserter(zone.internal_facets));
	return zone;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Deletes old facets from the restricted Delaunay triangulation

    void before_insertion_impl(const Facet& source_facet,
                               const Point& p,
                               Zone& zone) {
      if (tr.dimension() == 3) {
        bool source_facet_is_in_conflict = false;
	// On s'occupe des facettes de la zone de conflit
	for (typename Zone::Facets_iterator fit =
	       zone.internal_facets.begin();
	     fit != zone.internal_facets.end();
	     ++fit) {
          if(before_insertion_handle_facet_inside_conflict_zone (*fit,
                                                                 source_facet))
            source_facet_is_in_conflict = true;
        }

	for (typename Zone::Facets_iterator fit =
	       zone.boundary_facets.begin(); fit !=
	       zone.boundary_facets.end(); ++fit) {
          if(before_insertion_handle_facet_on_boundary_of_conflict_zone (*fit,
                                                                         source_facet))
            source_facet_is_in_conflict = true;
        }

        // source_facet == Facet() when this->before_insertion_impl is
        // called from a Mesh_3 visitor.
        if(source_facet != Facet() && !source_facet_is_in_conflict)
        {
          const Facet source_other_side = mirror_facet(source_facet);
          std::stringstream error_msg;
          error_msg << 
            boost::format("Surface_mesher ERROR: "
                          "A facet is not in conflict with its refinement point!\n"
                          "Debugging informations:\n"
                          "  Facet: (%1%, %2%) = (%6%, %7%, %8%)\n"
                          "  Dual: (%3%, %4%)\n"
                          "  Refinement point: %5%\n")
            % (&*source_facet.first) % source_facet.second
            % source_facet.first->circumcenter()
            % source_other_side.first->circumcenter()
            % p
            % source_facet.first->vertex((source_facet.second + 1)&3)->point()
            % source_facet.first->vertex((source_facet.second + 2)&3)->point()
            % source_facet.first->vertex((source_facet.second + 3)&3)->point();
          CGAL_error_msg(error_msg.str().c_str());
        }
      } // end if dimension() == 3

      // If dim < 3, then the triangulation has only one facet and the
      // complex has no facet, generically
      else {
	CGAL_assertion (tr.dimension() == 2);
	facets_to_refine.clear();
      }
//       CGAL_assertion_code(const Facet& source_other_side = mirror_facet(source_facet));
//       CGAL_assertion(
//         facets_to_refine.erase( (source_facet.first < source_other_side.first ) ?
//                                 source_facet : source_other_side) == false);
    }


    ///////////////////////////////////////////////////////////////////////////
    // Adds new facets from the restricted Delaunay triangulation

    void after_insertion_impl(const Vertex_handle& v) {
      CGAL_SURFACE_MESHER_PROFILER("inserted point")
#ifdef CGAL_SURFACE_MESHER_DEBUG_AFTER_INSERTION
      std::cerr << "Inserted\n";
#endif
      restore_restricted_Delaunay(v);
    }

    void restore_restricted_Delaunay(const Vertex_handle& v)
    {
      // On met a jour les flags des nouvelles facettes (celles a
      // l'interieur du trou et celles sur le bord du trou)

      std::vector<Cell_handle> cellules;
      tr.incident_cells (v, std::back_inserter(cellules));

      for(typename std::vector<Cell_handle>::iterator cellule =
	    cellules.begin();
	  cellule != cellules.end();
	  ++cellule) 
      {
	// Look at all four facets of the cell, starting with the
	// facet opposite to the new vertex
	int indice = (*cellule)->index (v);
	after_insertion_handle_opposite_facet (Facet (*cellule, indice));

	for (int i = 1; i <= 3; ++i)
	  after_insertion_handle_incident_facet (Facet (*cellule, (indice+i)&3));
      }


    } // restore_restricted_Delaunay end




    ///////////////////////////////////////////////////////////////////////////

    // Useless here
    void after_no_insertion_impl(const Facet&, const Point&, Zone&) {}




    ///////////////////////////////////////////////////////////////////////////
    // Private functions
  private:


    ///////////////////////
    // For before_insertion

    // Actions to perform on a facet inside the conflict zone
    bool before_insertion_handle_facet_inside_conflict_zone (const Facet& f,
                                                             const Facet& source_facet) 
    {
      const Facet other_side = mirror_facet(f);

      if(tr.is_infinite(f.first) && tr.is_infinite(other_side.first))
        return (f==source_facet) || (other_side == source_facet);

      // On enleve la facette de la liste des mauvaises facettes
#ifdef CGAL_SURFACE_MESHER_TAG_BAD
      if(f.first->is_bad(f.second))
#endif // CGAL_SURFACE_MESHER_TAG_BAD
        if(f.first < other_side.first)
          facets_to_refine.erase(f);
        else
          facets_to_refine.erase(other_side);
      // Le compteur des visites est remis a zero
       reset_visited(f);
       reset_visited(other_side);

      // On retire la facette du complexe (car on doit etre
      // independant de l'implementation du complexe)
      c2t3.remove_from_complex (f);
      return (f==source_facet) || (other_side == source_facet);
    }

    // Action to perform on a facet on the boundary of the conflict zone
    bool before_insertion_handle_facet_on_boundary_of_conflict_zone (const Facet& f,
                                                                     const Facet& source_facet) {
      // perform the same operations as for an internal facet
      return before_insertion_handle_facet_inside_conflict_zone (f, source_facet);
    }



    ///////////////////////
    // For after_insertion

    // Action to perform on a facet incident to the new vertex
    void after_insertion_handle_incident_facet (const Facet& f) {

      // If the facet is infinite or has been already visited,
      // then there is nothing to do as for it or its edges
      if (tr.is_infinite(f) || is_facet_visited(f))
	return;

      new_facet<false>(f);
    }

    // Action to perform on any new facet
    template <bool remove_from_complex_if_not_in_restricted_Delaunay>
    void new_facet (const Facet& f) 
    {
      const Facet other_side = mirror_facet(f);

      if(tr.is_infinite(f.first) && tr.is_infinite(other_side.first)) return;

      // NB: set_facet_visited() is implementation dependent
      // and each side of the real facet has to be considered
      // separately
      set_facet_visited(f);
      set_facet_visited(other_side);

      // On regarde d'abord si la facette est dans le Delaunay
      // restreint
      Point center;
      if (is_facet_on_surface(f, center)) {

	// NB: set_in_complex() is implementation independent and
	// consider both sides of the facet at once
	c2t3.set_in_complex(f);

	// NB: set_facet_surface_center() is implementation dependent
	// and each side of the real facet has to be considered
	// separately
	set_facet_surface_center(f, center);
	set_facet_surface_center(other_side, center);

	// On regarde alors si la facette est bonne
        Quality a_r;
	if (criteria.is_bad (f, a_r)) {
#ifdef CGAL_SURFACE_MESHER_TAG_BAD
          f.first->set_bad(f.second);
          other_side.first->set_bad(other_side.second);
#endif // CGAL_SURFACE_MESHER_TAG_BAD

          if(f.first < other_side.first)
            facets_to_refine.insert (f, a_r);
          else
            facets_to_refine.insert (other_side, a_r);
	}
      }
      else
        if( remove_from_complex_if_not_in_restricted_Delaunay )
          c2t3.remove_from_complex(f);
    }

    // Action to perform on a facet opposite to the new vertex
    void after_insertion_handle_opposite_facet (const Facet& f) {
      // perform the same operations as for a facet incident to the
      // new vertex
      after_insertion_handle_incident_facet (f);
    }



    ///////////////////////
    // Predicate to test restricted Delaunay membership


    // Tests whether a given facet is restricted or not
    bool is_facet_on_surface(const Facet& f, Point& center) {
      typedef typename Surface_mesh_traits::Intersect_3 Intersect_3;
      Intersect_3 intersect = meshtraits.intersect_3_object();
      typename GT::Is_degenerate_3 is_degenerate;
        
      Object dual = tr.dual(f);

      Object intersection;
      // If the dual is a segment
      if (const typename GT::Segment_3* segment_ptr=object_cast<typename GT::Segment_3>(&dual)) {
	if(is_degenerate(*segment_ptr)) return false;
        intersection = intersect(surf, *segment_ptr);
      }
      // If the dual is a ray
      else if(const typename GT::Ray_3* ray_ptr=object_cast<typename GT::Ray_3>(&dual)) {
        // If a facet is on the convex hull, and if its finite incident
        // cell has a very bid Delaunay ball, then the dual of the facet is
        // a ray constructed with a point with very big coordinates, and a
        // vector with small coordinates. Its can happen than the
        // constructed ray is degenerate (the point(1) of the ray is
        // point(0) plus a vector whose coordinates are espilon).
        if(is_degenerate(*ray_ptr)) return false;

	intersection = intersect(surf, *ray_ptr);
      }
      // If the dual is a line
      else if(const typename GT::Line_3* line_ptr = object_cast<typename GT::Line_3>(&dual)) {
	intersection = intersect(surf, *line_ptr);
      }
      // Else there is a problem with the dual
      else {
	std::cerr << "In is_facet_on_surface(const Facet& f, Point& center)\n"
		  << "file " << __FILE__ << ", line " << __LINE__ << "\n";
	std::cerr << "Incorrect object type: " << dual.type().name() << "\n";
        CGAL_error();
      }

      if(const Intersection_point* point_ptr = object_cast<Intersection_point>(&intersection))
      {
        center = static_cast<Point>(*point_ptr);
        return true;
      }
      return false;
    }

  public:
    std::string debug_info() const
    {
      std::stringstream s;
      s << facets_to_refine.size();
      return s.str();
    }

    std::string debug_info_header() const
    {
      return "#facets";
    }
  };  // end Surface_mesher_base

  template <
    typename Base,
    typename Element = typename details::Facet_generator<Base>::type,
    typename PreviousLevel = Null_mesher_level,
    Verbose_flag verbose = NOT_VERBOSE
    >
  struct Surface_mesher
    : public Base,
      public details::Mesher_level_generator<
        Base,
        Surface_mesher<Base, Element, PreviousLevel, verbose>,
        Element,
        PreviousLevel
      >::Type
  {
  public:
    typedef typename Base::Complex_2_in_triangulation_3 C2T3;
    typedef typename C2T3::Triangulation Triangulation;
    typedef Triangulation Tr;
    typedef typename Base::Surface Surface;
    typedef typename Base::Criteria Criteria;
    typedef typename Base::Surface_mesh_traits Surface_mesh_traits;

    typedef Surface_mesher<Base> Self;

    typedef typename details::Mesher_level_generator<
      Base,
      Surface_mesher<Base, Element, PreviousLevel, verbose>,
      Element,
      PreviousLevel
      >::Type Mesher_lvl;

    using Mesher_lvl::scan_triangulation;
    using Mesher_lvl::refine;
    using Mesher_lvl::is_algorithm_done;
    using Mesher_lvl::one_step;

    typedef C2T3 Complex_2_in_triangulation_3;

  private:
    Null_mesher_level null_mesher_level;
    Null_mesh_visitor null_visitor;
    bool initialized;

  public:
    Surface_mesher(C2T3& c2t3,
                   const Surface& surface,
                   const Surface_mesh_traits& mesh_traits,
                   const Criteria& criteria)
      : Base(c2t3, surface, mesh_traits, criteria), 
        Mesher_lvl(null_mesher_level),
        initialized(false)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Surface_mesher\n";
#endif
    }

    Surface_mesher(C2T3& c2t3,
                   const Surface& surface,
                   const Surface_mesh_traits& mesh_traits,
                   const Criteria& criteria,
		   PreviousLevel& previous)
      : Base(c2t3, surface, mesh_traits, criteria), 
        Mesher_lvl(previous),
        initialized(false)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Surface_mesher\n";
#endif
    }

    std::string debug_info() const
    {
      std::stringstream s;
      s << Base::debug_info() 
	<< "," << Mesher_lvl::previous().debug_info();
      return s.str();
    }

    std::string debug_info_header() const
    {
      std::stringstream s;
      s << Base::debug_info_header()
	<< "," << Mesher_lvl::previous().debug_info_header();
      return s.str();
    }

    // Initialization
    void init()
      {
	scan_triangulation();
	initialized = true;
      }

    void refine_mesh () {
      refine_mesh(null_visitor);
    }

    template <typename Visitor>
    void refine_mesh(Visitor visitor) {
      Tr& tr = this->triangulation_ref_impl();

      if(!initialized)
	init();

      if (verbose == NOT_VERBOSE)
	refine (visitor);
      else {
	std::cerr << "Refining...\n";
	int nbsteps = 0;
	CGAL::Timer timer;
        std::cerr << "Legende of the following line: "
                  << "(#vertices,#steps," << this->debug_info_header()
                  << ")\n";
	std::cerr 
	  << boost::format("\r             \r"
			   "(%1%,%2%,%3%) (%|4$.1f| vertices/s)")
	  % tr.number_of_vertices()
	  % nbsteps % this->debug_info()
	  % (nbsteps / timer.time());
	++nbsteps;
	timer.start();
	while (!is_algorithm_done()) {
	  {
	    CGAL_SURFACE_MESHER_TIME_PROFILER("Surface_mesher::one_step()");
	    one_step (visitor);
	  }
	  std::cerr 
	    << boost::format("\r             \r"
			     "(%1%,%2%,%3%) (%|4$.1f| vertices/s)")
	    % tr.number_of_vertices()
	    % nbsteps % this->debug_info()
	    % (nbsteps / timer.time());
	  ++nbsteps;
	}
	std::cerr << "\ndone.\n";
      }

      initialized = false;
    }

  };  // end Surface_mesher


  }  // end namespace Surface_mesher

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESHER_SURFACE_MESHER_H
