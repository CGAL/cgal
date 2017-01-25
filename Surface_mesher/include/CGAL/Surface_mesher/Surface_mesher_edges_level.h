// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
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
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_SURFACE_MESHER_EDGES_LEVEL_H
#define CGAL_SURFACE_MESHER_EDGES_LEVEL_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Mesh_2/Output_stream.h>

#include <CGAL/Mesher_level.h>
#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>
#include <CGAL/Meshes/Simple_map_container.h>
#include <CGAL/iterator.h> // for CGAL::inserter
#include <CGAL/circulator.h> // for CGAL::Circulator_from_container<C>
#include <CGAL/use.h>
#include <CGAL/number_utils.h> // for CGAL::sqrt

#include <sstream>
#include <set>

#include <boost/format.hpp>


namespace CGAL {
namespace Surface_mesher {

//   template <class Tr>
//   typename Tr::Edge
//   canonical_edge(const Tr& tr, 
// 		 const typename Tr::Edge e) 
//   {
//     typedef typename Tr::Cell_handle Cell_handle;
//     typedef typename Tr::Vertex_handle Vertex_handle;

//     const Vertex_handle& va = e.first->vertex(e.second);
//     const Vertex_handle& vb = e.first->vertex(e.third);

//     typename Tr::Cell_circulator
//       ccirc = tr.incident_cells(e),
//       begin(ccirc);
//     Cell_handle cell = ccirc;
//     do {
//       if(Cell_handle(++ccirc) < cell)
// 	cell = ccirc;
//     } while(ccirc != begin);
//     return typename Tr::Edge(cell,
// 			     cell->index(va),
// 			     cell->index(vb));
//   }

//   template <class Tr,
// 	    typename ForwardIterator, typename OutputIterator>
//   void canonical_edges(const Tr& tr,
// 		       ForwardIterator begin,
// 		       ForwardIterator end,
// 		       OutputIterator it)
//   {
//     typedef typename Tr::Edge Edge;

//     for(ForwardIterator cit = begin; cit != end; ++cit)
//       for(int i = 1; i < 4; ++i)
// 	for(int j = 0 ; j < i; ++j)
// 	  *it++=canonical_edge(tr, Edge(*cit, i, j));
//   }

  template <typename Tr, 
	    typename Surface_mesh_traits>
  Object 
  compute_edge_intersection_curve(const Tr& tr,
				  const typename Surface_mesh_traits::
				    Surface_3& surf,
				  const Surface_mesh_traits& meshtraits,
				  const typename Tr::Edge& e)
  {
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
    std::cerr << 
      boost::format("compute_edge_intersection_curve(Edge(%1%, %2%, %3%)="
		    "(%4%, %5%))\n")
      % &*e.first % e.second % e.third
      % e.first->vertex(e.second)->point()
      % e.first->vertex(e.third)->point();
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
    typedef typename Tr::Geom_traits GT;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Cell_handle Cell_handle;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::Triangle_3 Triangle_3;
    typedef typename GT::Vector_3 Vector_3;
    typedef typename GT::FT FT;
    typedef typename Surface_mesh_traits::Intersection_point Intersection_point;

    typename GT::Construct_midpoint_3 midpoint = 
      tr.geom_traits().construct_midpoint_3_object();
    typename GT::Construct_ray_3 create_ray = 
      tr.geom_traits().construct_ray_3_object();
    typename GT::Construct_vector_3 create_vector = 
      tr.geom_traits().construct_vector_3_object();
    typename GT::Construct_cross_product_vector_3 cross_product = 
      tr.geom_traits().construct_cross_product_vector_3_object();
    typename GT::Construct_point_on_3 point_on = 
      tr.geom_traits().construct_point_on_3_object();
    typename GT::Construct_translated_point_3 translate = 
      tr.geom_traits().construct_translated_point_3_object();
    typename GT::Construct_divided_vector_3 divide = 
      tr.geom_traits().construct_divided_vector_3_object();
    typename GT::Construct_scaled_vector_3 scale = 
      tr.geom_traits().construct_scaled_vector_3_object();
    typename GT::Construct_sum_of_vectors_3 sum = 
      tr.geom_traits().construct_sum_of_vectors_3_object();
    typename GT::Compute_squared_length_3 sq_length =
      tr.geom_traits().compute_squared_length_3_object();

    typedef std::vector<Point_3> Edge_dual;
    Edge_dual edge_dual;
    edge_dual.reserve(12);

    // va and vb are the extremities of the edge. 
    const Vertex_handle& va = e.first->vertex(e.second);
    const Vertex_handle& vb = e.first->vertex(e.third);

    // m is the middle of [va->point(), vb->point()]
    const Point_3 edge_middle = midpoint(va->point(), vb->point());

    // square of (four times the radius of the bounding sphere)
    const FT sq_radius = surf.bounding_sphere_squared_radius() * FT(16);

    // If there is an infinite cell in the incident cells of e, it means
    // that e in on the convex hull of the triangulation.
    // In that case, there will be exactly two infinite cells, and they
    // will be consecutive. We do not want to start

    typename Tr::Cell_circulator circ = tr.incident_cells(e, e.first);
    while(tr.is_infinite(circ)) ++circ; // Start at the first finite
					// incident cell.
    typename Tr::Cell_circulator begin = circ;

    typedef typename Edge_dual::size_type size_type;
    CGAL_USE_TYPE(size_type);
    CGAL_assertion_code(size_type number_of_finite_incident_cells = 0);
    CGAL_assertion_code(size_type number_of_infinite_incident_cells = 0);

    do {
      const Cell_handle current_cell = circ;
      const Cell_handle next_cell = ++circ; // increment circ once

      // index_a and index_b are the indices of va and vb in the cell *circ.
      int index_a = -1;
      int index_b = -1;
      CGAL_assertion_code(bool va_is_vertex =) current_cell->has_vertex(va, index_a);
      CGAL_assertion_code(bool vb_is_vertex =) current_cell->has_vertex(vb, index_b);
      CGAL_assertion( va_is_vertex && vb_is_vertex );

      // index_c and index_d are the two others indices in *current_cell.
      const int index_c = tr.next_around_edge(index_a, index_b);
      const int index_d = tr.next_around_edge(index_b, index_a);
      CGAL_assertion(index_c == current_cell->index(next_cell));
      // vertices c and d are the corresponding two vertices, in *circ.
      const Vertex_handle& vc = current_cell->vertex(index_c);
      const Vertex_handle& vd = current_cell->vertex(index_d);

      // If 'current_cell' is infinite, then that is vd which is infinite.
      CGAL_assertion(!tr.is_infinite(vc));
      if(tr.is_infinite(vd))
      {
	const Cell_handle cell_after = ++circ; // increment cell once again.
	CGAL_assertion_code(++++number_of_infinite_incident_cells);

	CGAL_assertion(tr.is_infinite(current_cell));
	CGAL_assertion(tr.is_infinite(next_cell));
	CGAL_assertion(cell_after != current_cell);
	CGAL_assertion(!tr.is_infinite(cell_after));

	// last_circumcenter is the circumcenter of the cell before
	// 'current_cell'. That cell is finite.
	const Point_3& last_circumcenter = edge_dual[edge_dual.size()-1];
	// next_circumcenter is the circumcenter of 'cell_after'. That cell
	// is finite.
	Point_3 next_circumcenter;
	// 'cell_after' can be 'begin'.
	if(cell_after != begin)
	  next_circumcenter = tr.dual(cell_after);
	else
	  next_circumcenter = edge_dual[0];

	// remember vc is finite.
	const Vector_3 first_vector = 
	  cross_product(create_vector(va->point(),
				      vb->point()),
			create_vector(va->point(),
				      vc->point()));

	// 'index_d_in_next' is the index of the other finite vertex (with
	// va and vd), in 'next_cell'.
	const int index_d_in_next = next_cell->index(current_cell);
	const Vertex_handle& vd_in_next = next_cell->vertex(index_d_in_next);
	CGAL_assertion_code(tr.is_infinite(vd_in_next));

	const Vector_3 second_vector = 
	  cross_product(create_vector(vb->point(), // special
				      va->point()),// order
			create_vector(vb->point(),
				      vd_in_next->point()));

	const Vector_3 middle_vector = sum(first_vector, second_vector);

	const FT sq_norm_first_vector = sq_length(first_vector);
	const FT sq_norm_second_vector = sq_length(second_vector);
	const FT sq_norm_middle_vector = sq_length(middle_vector);

#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
	std::cerr << ::boost::format("lengths=%1%, %2%, %3%\n")
	  % sq_norm_first_vector
	  % sq_norm_middle_vector
	  % sq_norm_second_vector;
#endif

	edge_dual.push_back(translate(last_circumcenter, 
				      scale(first_vector,
					    CGAL::sqrt(sq_radius / sq_norm_first_vector))));
	edge_dual.push_back(translate(edge_middle,
				      scale(middle_vector,
					    CGAL::sqrt(sq_radius / sq_norm_middle_vector))));
	edge_dual.push_back(translate(next_circumcenter,
				      scale(second_vector,
					    CGAL::sqrt(sq_radius / sq_norm_second_vector))));
	if(cell_after != begin) 
	{
	  edge_dual.push_back(next_circumcenter);
	  ++circ; // again (third time), increment circ
	  CGAL_assertion_code(++number_of_finite_incident_cells);
	}
      }
      else // the cell 'current_cell' is finite
      { 
	CGAL_assertion_code(++number_of_finite_incident_cells);
	CGAL_assertion(!tr.is_infinite(current_cell));

        // Use tr.dual(), which is optimized, when the cell base class as
        // circumcenter().
	edge_dual.push_back(tr.dual(current_cell));

// 	const int index_of_next_cell = current_cell->index(next_cell);

// 	{ // center of the (finite) facet between current_cell and next_cell
// 	  const int i = index_of_next_cell(index_of_next_cell, 0);
// 	  const int j = index_of_next_cell(index_of_next_cell, 1);
// 	  const int k = index_of_next_cell(index_of_next_cell, 2);
// 	  edge_dual.push_back(circumcenter(current_cell->vertex(i),
// 					   current_cell->vertex(j),
// 					   current_cell->vertex(k)));
// 	}
      } // end of the case when 'current_cell' is finite
    } while(circ != begin);

#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
    CGAL_assertion_code(
    std::cerr << 
      boost::format("  number of finite/infinite incident cells: %1%/%2%\n"
		    "  edge_dual.size(): %3%\n")
      % number_of_finite_incident_cells
      % number_of_infinite_incident_cells
      % edge_dual.size();
    )
#endif
    CGAL_assertion(number_of_infinite_incident_cells == 0 ||
		   number_of_infinite_incident_cells == 2);
    CGAL_assertion(edge_dual.size() == number_of_finite_incident_cells +
		   number_of_infinite_incident_cells +
		   (number_of_infinite_incident_cells!=0 ? 1 : 0) );
    CGAL_assertion(edge_dual.size() >= 3);

    Circulator_from_container<Edge_dual> 
      vor_circ(&edge_dual), vor_begin(vor_circ); //Circulator in the list of
						//points of the edge dual
						//(a Voronoi facet).

    const Point_3& first = *vor_circ++;
    --vor_begin;
    while(vor_circ != vor_begin) {
      const Point_3& current = *vor_circ;
      const Point_3& next = *++vor_circ;

      Object o = 
	meshtraits.intersect_curves_with_triangle(surf,
						  Triangle_3(current,
							     next,
							     first));
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
      CGAL_assertion_code(
      std::cerr << boost::format("  triangle test (%1%, %2%, %3%) = %4%\n")
      % current % next % first % !o.is_empty();
      )
        std::cerr << boost::format("intersecion type: %1%\n") % o.type().name();
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
      if(const Intersection_point* point = object_cast<Intersection_point>(&o))
      {
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
        CGAL_assertion_code(
        std::cerr << boost::format("  result=(%1%)\n")
        % (*point);
			    )
#endif
        return make_object(static_cast<Point_3>(*point));
      }
    }

    return Object();
  } // end of compute_edge_intersection_curve(Tr, Surface_3, SMTraits, Edge)

  template <typename Tr, 
	    typename Surface_mesh_traits>
  typename Tr::Point 
  lineic_center(const Tr& tr,
		  const typename Surface_mesh_traits::Surface_3& surf,
		  const Surface_mesh_traits& meshtraits,
		  const typename Tr::Edge& e)
  {
    typedef typename Surface_mesh_traits::Intersection_point Intersection_point;
    // the following object cast can throw a 'Bad_object_cast' exception 
    Intersection_point point = 
      object_cast<Intersection_point>(compute_edge_intersection_curve(tr,
					     surf,
					     meshtraits,
					     e));
    return static_cast<typename Tr::Point>(point);
  }

  namespace details {

    /** This class defines several auxiliary types. */
    template <typename C2T3>
    struct Surface_mesher_edges_base_types
    {
      typedef C2T3 Complex_2_in_triangulation_3;
      typedef typename C2T3::Triangulation Tr;
      typedef typename Tr::Edge Edge;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Point Point_3;
      typedef std::pair<Vertex_handle, Vertex_handle> Pair_of_vertices;
      typedef std::pair<Edge, Point_3> Edge_info;
      typedef ::CGAL::Meshes::Simple_map_container<Pair_of_vertices,
						   Edge_info> Default_container;
    }; // end struct Surface_mesher_edges_base_types

  } // end namespace details

  template <
    class C2T3,
    class Surface_,
    class SurfaceMeshTraits,
    class EdgesCriteria,
    class Container = 
      typename details::Surface_mesher_edges_base_types<C2T3>::Default_container
    >
  class Surface_mesher_edges_level_base :
    public No_after_no_insertion,
    public No_before_conflicts,
    public Triangulation_mesher_level_traits_3<typename C2T3::Triangulation>,
    public Container
  {
  public:
    typedef C2T3 Complex_2_in_triangulation_3;
    typedef typename C2T3::Edge_info Edge_info;
    typedef Surface_ Surface;
    typedef SurfaceMeshTraits Surface_mesh_traits;
    typedef EdgesCriteria Edges_criteria;
    typedef Edges_criteria Criteria;

    typedef typename C2T3::Triangulation Tr;
    typedef typename Tr::Point Point_3;
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Edge Edge;
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Cell_handle Cell_handle;

    typedef std::pair<Vertex_handle, Vertex_handle> Pair_of_vertices;

    typedef typename Tr::Geom_traits GT;
    typedef typename GT::Triangle_3 Triangle_3;


    typedef Triangulation_mesher_level_traits_3<Tr> Triangulation_traits;
    typedef typename Triangulation_traits::Zone Zone;

    typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;

  protected:
    C2T3& c2t3;
    Tr& tr;     // Associated triangulation reference
    const Surface& surf;  // Surface
    const Surface_mesh_traits& meshtraits; // Surface mesh traits
    const Edges_criteria& criteria;  // Meshing criteria
    Point_3 refinement_point_cache;

  private:
    template <typename T>
    static
    void order_pair(T& a, T& b)
    {
      if( b < a )
	std::swap(a, b);
    }

    template <typename T>
    static
    std::pair<T, T>
    make_ordered_pair(const T& a, const T& b)
    {
      if( a < b )
	return std::make_pair(a, b);
      else
	return std::make_pair(b, a);
    }

    static
    Pair_of_vertices make_pair_of_vertices(const Vertex_handle& va,
					   const Vertex_handle& vb)
    {
      return make_ordered_pair(va, vb);
    }

    static
    Pair_of_vertices make_pair_of_vertices(const Edge e)
    {
      return make_pair_of_vertices(e.first->vertex(e.second),
				   e.first->vertex(e.third));
    }

    template <typename ForwardIterator, typename OutputIterator>
    static
    void all_pair_of_vertices(ForwardIterator begin,
			      ForwardIterator end,
			      OutputIterator it)
    {
      for(ForwardIterator cit = begin; cit != end; ++cit)
	for(int i = 1; i < 4; ++i)
	  for(int j = 0 ; j < i; ++j)
	    *it++=make_pair_of_vertices(Edge(*cit, i, j));
    }

    template <typename ForwardIterator, typename OutputIterator>
    static
    void all_edges(ForwardIterator begin,
		   ForwardIterator end,
		   OutputIterator it)
    {
      std::set<Pair_of_vertices> edges;
      for(ForwardIterator cit = begin; cit != end; ++cit)
	for(int i = 1; i < 4; ++i)
	  for(int j = 0 ; j < i; ++j)
	  {
	    const Pair_of_vertices pair = make_pair_of_vertices(Edge(*cit, i, j));
	    if(edges.find(pair) == edges.end())
	    {
	      *it++=Edge(*cit, i, j);
	      edges.insert(pair);
	    }
	  }
    }
  public:
    // Constructor
    Surface_mesher_edges_level_base (C2T3& co, 
				     const Surface& s, 
				     const Surface_mesh_traits& mesh_traits,
				     const Edges_criteria& c) :
      Triangulation_mesher_level_traits_3<Tr>(co.triangulation()),
      c2t3(co),
      tr(co.triangulation()),
      surf(s),
      meshtraits(mesh_traits),
      criteria(c)
    {
#ifdef CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
      std::cerr << "CONS: Surface_mesher_edges_level_base\n";
#endif // CGAL_SURFACE_MESHER_DEBUG_CONSTRUCTORS
    }

    Tr& triangulation_ref_impl()
    {
      return tr;
    }

    const Tr& triangulation_ref_impl() const
    {
      return tr;
    }

    Edge get_next_element_impl()
    {
      const typename Container::value_type container_value = 
	Container::get_next_element_impl();
      refinement_point_cache = container_value.second.second;
      const Edge& e = container_value.second.first;
      CGAL_assertion_code(const Vertex_handle& va = container_value.first.first);    
      CGAL_assertion_code(const Vertex_handle& vb = container_value.first.second);    
      CGAL_assertion_code(const Cell_handle& c = e.first);
      CGAL_assertion_code(const int i = e.second);
      CGAL_assertion_code(const int j = e.third);
      CGAL_assertion(va < vb);
      CGAL_assertion( (c->vertex(i) == va) ? c->vertex(j) == vb : (c->vertex(j) == va && c->vertex(i) == vb) );
      return e;
    }

    Point_3 refinement_point_impl(const Edge&) const
    {
      return refinement_point_cache;
    }

    Mesher_level_conflict_status private_test_point_conflict_impl(const Point_3& p,
								  Zone& zone)
    {
      if( zone.locate_type == Tr::VERTEX )
      {
	std::cerr << boost::format("Error: (%1%) is already inserted on edge\n") % p;
	return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
      }
      else
	return NO_CONFLICT;
    }

    Zone conflicts_zone_impl(const Point_3& p,
			     const Edge& e) const 
    {
      Zone zone;

      // TODO may be avoid the locate here
      zone.cell =
	tr.locate(p, zone.locate_type, zone.i, zone.j, e.first);
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      std::cerr << 
	boost::format("-> edge Edge(%1%, %2%, %3%)=(%4%, %5%)\n"
		      "     insertion point: %6% (locate_type=%7%)\n")
	% &*e.first % e.second % e.third
	% e.first->vertex(e.second)->point()
	% e.first->vertex(e.third)->point()
	% p
	% zone.locate_type;
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      if(zone.locate_type != Tr::VERTEX)
        tr.find_conflicts(p, zone.cell,
                          std::back_inserter(zone.boundary_facets),
                          std::back_inserter(zone.cells),
                          std::back_inserter(zone.internal_facets));
      return zone;
    }

    void scan_triangulation_impl() 
    {
#ifdef CGAL_SURFACE_MESHER_VERBOSE
      std::cout << "scanning edges (curves)..." << std::endl;
#endif // CGAL_SURFACE_MESHER_VERBOSE
      for(Finite_edges_iterator eit = tr.finite_edges_begin();
	  eit != tr.finite_edges_end();
	  ++eit)
	new_edge(*eit);
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
      std::cout << "number of edges: " << this->size() << std::endl;
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INTERSECTION
    }

    Mesher_level_conflict_status
    test_point_conflict_from_superior_impl(const Point_3& p,
					   Zone& zone)
    {
      std::set<Edge> edges;

      all_edges(zone.cells.begin(), zone.cells.end(), 
		inserter(edges));

      for(typename std::set<Edge>::iterator eit = edges.begin();
	  eit != edges.end();
	  ++eit)
	if( test_if_edge_is_encroached(*eit, p) )
	  return CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;
      
      return NO_CONFLICT;
    }

    void before_insertion_impl(const Edge& CGAL_assertion_code(e),
                               const Point_3&,
                               const Zone& zone)
    {
      CGAL_assertion_code(bool is_e_removed = false);
      CGAL_assertion_code(const Cell_handle& c = e.first);
      CGAL_assertion_code(Vertex_handle va = c->vertex(e.second));
      CGAL_assertion_code(Vertex_handle vb = c->vertex(e.third));
      CGAL_assertion_code((order_pair(va, vb)));

      std::set<Edge> edges;
      all_edges(zone.cells.begin(), zone.cells.end(), 
		CGAL::inserter(edges));

#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      int number_of_edges_removed = 0;
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS

      for(typename std::set<Edge>::iterator eit = edges.begin();
	  eit != edges.end();
	  ++eit)
	if( !tr.is_infinite(*eit) &&
	    c2t3.unmark(*eit) )
	{
	  CGAL_assertion_code(Vertex_handle eit_va = eit->first->vertex(eit->second));
	  CGAL_assertion_code(Vertex_handle eit_vb = eit->first->vertex(eit->third));
	  CGAL_assertion_code(order_pair(eit_va, eit_vb));
	  CGAL_assertion_code(if(va == eit_va && vb == eit_vb) is_e_removed = true);
// 	  c2t3.remove_from_complex(*eit);
	  this->remove_element(make_pair_of_vertices(*eit));
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
	  ++number_of_edges_removed;
#endif
	}
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      std::cerr << 
	boost::format("     before insertion: remove %1% edges\n")
	% number_of_edges_removed;
#endif
//       CGAL_assertion(is_e_removed == true);
    }

    // for visitors
    void remove_edges(const Point_3&, const Zone& zone)
    {
      std::set<Edge> edges;
      all_edges(zone.cells.begin(), zone.cells.end(), 
		CGAL::inserter(edges));

#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      int number_of_edges_removed = 0;
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS

      for(typename std::set<Edge>::iterator eit = edges.begin();
	  eit != edges.end();
	  ++eit)
	if( !tr.is_infinite(*eit) &&
	    c2t3.unmark(*eit) )
	{
	  this->remove_element(make_pair_of_vertices(*eit));
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
	  ++number_of_edges_removed;
#endif
	}
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      std::cerr << 
	boost::format("     before insertion: remove %1% edges\n")
	% number_of_edges_removed;
#endif
    }

    void after_insertion_impl(const Vertex_handle& v)
    {
      CGAL_MESHES_OUTPUT_STREAM << "-";
      std::vector<Cell_handle> cellules;
      tr.incident_cells (v, std::back_inserter(cellules));

      std::set<Edge> edges;

      all_edges(cellules.begin(), cellules.end(),
		CGAL::inserter(edges));

#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      std::cerr << 
	boost::format("     after insertion: %1% new edges \n")
	% edges.size();
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      for(typename std::set<Edge>::iterator eit = edges.begin();
	  eit != edges.end();
	  ++eit)
	new_edge(*eit);
#ifdef CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
      std::cerr << 
	boost::format("     c2t3.number_of_marked_edges=%1%\n"
		      "     number of bad edges=%2%\n")
	% c2t3.number_of_marked_edges()
	% this->size();
#endif // CGAL_SURFACE_MESHER_EDGES_DEBUG_INSERTIONS
    }

  private: 
    /** Test if the edge e is encroached by the point p, that is if the
	circumscribed sphere centered a the surfacic center of e contains
	the point p in its interior. */
    bool test_if_edge_is_encroached(const Edge e, const Point_3& p)
    {
      typename GT::Compute_squared_distance_3 sq_distance =
	tr.geom_traits().compute_squared_distance_3_object();

      if( !tr.is_infinite(e) &&
	  c2t3.is_marked(e) )
      {
	Edge_info& info = c2t3.get_info(e);
	if(!info.is_cached)
	{
	  info.lineic_center = lineic_center(tr, 
					     surf,
					     meshtraits,
					     e);
	  info.is_cached = true;
	}
	if(sq_distance(info.lineic_center, p) <
	   sq_distance(info.lineic_center,
		       e.first->vertex(e.second)->point()))
	{
	  this->add_bad_element(make_pair_of_vertices(e),
			    std::make_pair(e,
					   info.lineic_center));
	  return true;
	}
	else
	  return false;
      }
      else
	return false;
    }

    bool is_edge_in_restricted_triangulation(const Edge& e,
					     Point_3& p) const
    {
      typedef typename Surface_mesh_traits::Intersection_point Intersection_point;
      const Object obj = compute_edge_intersection_curve(tr,
                                                         surf,
                                                         meshtraits,
                                                         e);
      if(const Intersection_point* point = 
         object_cast<Intersection_point>(&obj))
      {
        p = static_cast<Point_3>(*point);
        return true;
      }
      else
        return false;
    }

    void new_edge(const Edge e)
    {
      if(tr.is_infinite(e)) return;

      Point_3 center;
      if(is_edge_in_restricted_triangulation(e, center))
      {
// 	CGAL_assertion(e == canonical_edge(tr, e));
	Edge_info info;
	info.is_cached = true;
	info.lineic_center = center;
	c2t3.mark(e, info);
	if(criteria.is_bad(e, center))
	{
	  this->add_bad_element(make_pair_of_vertices(e),
			    std::make_pair(e,
					   center));
	}
      }
    }

  protected: // --- public functions used in derived classed ---

    std::string debug_info() const
    {
      std::stringstream s;
      s << this->size();
      return s.str();
    }

    static std::string debug_info_header()
    {
      return "#edges";
    }

    bool check_restricted_delaunay () 
    {
      // to be done
      return true;
    }
  }; // end class Surface_mesher_edges_level_base

} // end namespace Surface_mesher
} // end namespace CGAL::

#endif // CGAL_SURFACE_MESHER_EDGES_LEVEL_H
