// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008,2013  GeometryFactory, Sophia Antipolis (France)
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
// Author(s)     : Steve Oudot, David Rey, Mariette Yvinec, Laurent Rineau, Andreas Fabri

#ifndef CGAL_MESH_3_REFINE_FACETS_MANIFOLD_BASE_H
#define CGAL_MESH_3_REFINE_FACETS_MANIFOLD_BASE_H

#include <CGAL/utility.h>
#include <set>
#include <vector>

namespace CGAL {

namespace Mesh_3 {

template <
  typename Base_,
  bool withBoundary = false
  >
class Refine_facets_manifold_base
  : public Base_
{
public:

  typedef Base_ Base ;

  typedef typename Base::C3T3 C3t3;
  typedef typename Base::Criteria Criteria;
  typedef typename Base::Mesh_domain Mesh_domain;
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Geom_traits GT;
  typedef typename GT::FT FT;
  typedef typename Tr::Point Point;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Facet_circulator Tr_facet_circulator;

  typedef std::pair<Vertex_handle, Vertex_handle> EdgeVV;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;

protected:
  Tr& r_tr_;
  C3t3& r_c3t3_;
  mutable std::set <std::pair <Vertex_handle, Vertex_handle> > m_bad_edges;
  mutable std::set<Vertex_handle> m_bad_vertices;

  mutable bool m_manifold_info_initialized;
  mutable bool m_bad_vertices_initialized;
  bool m_with_boundary;

protected:
  // computes and return an ordered pair of Vertex
  EdgeVV make_edgevv(const Vertex_handle vh1,
                     const Vertex_handle vh2) const {
    if (vh1 < vh2) {
      return std::make_pair(vh1, vh2);
    }
    else {
      return std::make_pair(vh2, vh1);
    }
  }

  // computes and returns the Edge type object from the EdgeVV object
  // use tr.is_edge(...)
  Edge edgevv_to_edge(const EdgeVV& arete) const {
    Vertex_handle sommet1 = arete.first;
    Vertex_handle sommet2 = arete.second;
    Cell_handle c;
    int index1, index2;

    CGAL_assertion_code(bool is_edge =)
    r_tr_.is_edge(sommet1, sommet2, c, index1, index2);
    CGAL_assertion(is_edge);

    return make_triple( c, index1, index2 );
  }

  // computes and returns the EdgeVV type object from the Edge object
  EdgeVV edge_to_edgevv(const Edge& arete) const {
    return make_edgevv(arete.first->vertex(arete.second),
                       arete.first->vertex(arete.third) );
  }

  FT compute_distance_to_facet_center(const Facet& f,
                                      const Vertex_handle v) const {
    const Point& fcenter = f.first->get_facet_surface_center(f.second);
    const Point& vpoint = v->point();

    return r_tr_.geom_traits().compute_squared_distance_3_object()(fcenter,
                                                                vpoint );
  }

  Facet
  biggest_incident_facet_in_complex(const Vertex_handle sommet) const
  {

    std::vector<Facet> facets;
    facets.reserve(64);
    r_tr_.incident_facets(sommet, std::back_inserter(facets));

    typename std::vector<Facet>::iterator fit = facets.begin();
    while(fit != facets.end() && !r_c3t3_.is_in_complex(*fit)) ++fit;
    CGAL_assertion(fit!=facets.end());

    Facet biggest_facet = *fit;

    for (++fit; fit != facets.end(); )
    {
      Facet current_facet = *fit;
      // is the current facet bigger than the current biggest one
      if ( compute_distance_to_facet_center(current_facet, sommet) >
           compute_distance_to_facet_center(biggest_facet, sommet) )
      {
        biggest_facet = current_facet;
      }
      ++fit;
      while(fit != facets.end() && !r_c3t3_.is_in_complex(*fit)) ++fit;
    }
    return biggest_facet;
  }

  Facet biggest_incident_facet_in_complex(const Edge& arete) const {
    // Find the first facet in the incident facets
    // of the edge which is in the Complex
    // use the list of incident facets in the complex
    Tr_facet_circulator fcirc = r_tr_.incident_facets(arete);
    while(!r_c3t3_.is_in_complex(*fcirc)) ++fcirc;
    Facet first_facet = *fcirc;
    Facet biggest_facet = *fcirc;

    for (++fcirc; *fcirc != first_facet; ++fcirc) 
    {
      while(!r_c3t3_.is_in_complex(*fcirc)) ++fcirc;
      if(*fcirc == first_facet) break;

      Vertex_handle fev = edge_to_edgevv(arete).first;
      // is the current facet bigger than the current biggest one
      if ( compute_distance_to_facet_center(*fcirc, fev) >
           compute_distance_to_facet_center(biggest_facet,
                                            fev) ) {
        biggest_facet = *fcirc;
      }
      else { // @TODO: il ne faut pas aller voir des deux cotes: c'est
        // le meme centre de facet!!!
        Facet autre_cote = r_tr_.mirror_facet(*fcirc);
        // is the current facet bigger than the current biggest one
        if ( compute_distance_to_facet_center(autre_cote, fev) >
             compute_distance_to_facet_center(biggest_facet, fev) ) {
          biggest_facet = autre_cote;
        }
      }
    }
    return biggest_facet;
  }

  ///////////////////////
  // For before_insertion

  // Actions to perform on a facet inside the conflict zone
  void
  before_insertion_handle_facet_inside_conflict_zone (const Facet& f) 
  {
    if ( r_c3t3_.is_in_complex(f) ) {
      // foreach edge of f
      const Cell_handle cell = f.first;
      const int i = f.second;
      for(int j = 0; j < 3; ++j)
      {
        const int edge_index_va = r_tr_.vertex_triple_index(i, j);
        const int edge_index_vb = r_tr_.vertex_triple_index(i, (j == 2) ? 0 : (j+1));
        const Vertex_handle edge_va = cell->vertex(edge_index_va);
        const Vertex_handle edge_vb = cell->vertex(edge_index_vb);
        m_bad_edges.erase( make_edgevv(edge_va, edge_vb));
      }
    }
  }

  // Action to perform on a facet on the boundary of the conflict zone
  void
  before_insertion_handle_facet_on_boundary_of_conflict_zone(const Facet& f)
  {
    // perform the same operations as for an internal facet
    before_insertion_handle_facet_inside_conflict_zone (f);

    if(m_bad_vertices_initialized) {
      const Cell_handle& c = f.first;
      const int i = f.second;

      // for each v of f
      for (int j = 0; j < 4; j++)
        if (i != j)
          m_bad_vertices.erase(c->vertex(j));
    }
  }

public:
  Refine_facets_manifold_base(Tr& triangulation,
                              const Criteria& criteria,
                              const Mesh_domain& oracle,
                              C3t3& c3t3,
                              bool with_boundary = withBoundary)
    : Base(triangulation,
           criteria,
           oracle,
           c3t3)
    , r_tr_(triangulation)
    , r_c3t3_(c3t3)
    , m_manifold_info_initialized(false)
    , m_bad_vertices_initialized(false)
    , m_with_boundary(with_boundary)
  {
#ifdef CGAL_MESH_3_DEBUG_CONSTRUCTORS
    std::cerr << "CONS: Refine_facets_manifold_base";
    if(m_with_boundary)
      std::cerr << " (with boundaries)\n";
    else
      std::cerr << " (without boundary)\n";
#endif
  }

private:
  // Initialization function
  void initialize_manifold_info() const {
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "\nscanning edges ";
    if(m_with_boundary)
      std::cerr << "(boundaries allowed)";
    std::cerr << "...\n";
    int n = 0;
#endif
    for (typename Tr::Finite_edges_iterator
           eit = r_tr_.finite_edges_begin(), end = r_tr_.finite_edges_end();
         eit != end; ++eit)
    {
      if ( (r_c3t3_.face_status(*eit) == C3t3::SINGULAR) ||
           ( (!m_with_boundary) &&
             (r_c3t3_.face_status(*eit) == C3t3::BOUNDARY) ) )
      {
        m_bad_edges.insert( edge_to_edgevv(*eit) );
#ifdef CGAL_MESH_3_VERBOSE
        ++n;
#endif
      }
    }
    m_manifold_info_initialized = true;
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "   -> found " << n << " bad edges\n";
#endif
  }

  void initialize_bad_vertices() const
  {
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "\nscanning vertices..." << std::endl;
    int n = 0;
#endif
    for (typename Tr::Finite_vertices_iterator
           vit = r_tr_.finite_vertices_begin(),
           end = r_tr_.finite_vertices_end();
         vit != end; ++vit)
    {
      if( r_c3t3_.face_status(vit) == C3t3::SINGULAR ) {
        m_bad_vertices.insert( vit );
#ifdef CGAL_MESH_3_VERBOSE
        ++n;
#endif
      }
    }
    m_bad_vertices_initialized = true;
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "   -> found " << n << " bad vertices\n";
#endif
  }

public:
  void scan_triangulation_impl() {
    Base::scan_triangulation_impl();
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "scanning edges (lazy)" << std::endl;
    std::cerr << "scanning vertices (lazy)" << std::endl;
#endif
  }

  // Tells whether there remain elements to refine
  bool no_longer_element_to_refine_impl() const {
    if(Base::no_longer_element_to_refine_impl())
    {
      if( ! m_manifold_info_initialized ) initialize_manifold_info();

      if(m_bad_edges.empty())
      {
        if( ! m_bad_vertices_initialized ) initialize_bad_vertices();

        return m_bad_vertices.empty();
      }
      else // m_bad_vertices is not empty
        return false;
    }
    else // Base::no_longer_element_to_refine_impl() returned false
      return false;
  }

  // Returns the next element to refine
  Facet get_next_element_impl() {

    if (!Base::no_longer_element_to_refine_impl()) {
      return Base::get_next_element_impl();
    }
    else {
      if(!m_bad_edges.empty()) {
        Edge first_bad_edge = edgevv_to_edge(*(m_bad_edges.begin()));
        return biggest_incident_facet_in_complex(first_bad_edge);
      } else {
        return biggest_incident_facet_in_complex(*m_bad_vertices.begin());
      }
    }
  }

  void before_insertion_impl(const Facet& f, const Point& s,
                             Zone& zone) {
    if( m_manifold_info_initialized )
    {
      for (typename Zone::Facets_iterator fit =
             zone.internal_facets.begin();
           fit != zone.internal_facets.end();
           ++fit)
        before_insertion_handle_facet_inside_conflict_zone (*fit);

      for (typename Zone::Facets_iterator fit =
             zone.boundary_facets.begin(); fit !=
             zone.boundary_facets.end(); ++fit)
        before_insertion_handle_facet_on_boundary_of_conflict_zone (*fit);
    }
    Base::before_insertion_impl(f, s, zone);
  }

  void after_insertion_impl(const Vertex_handle v) {
    Base::after_insertion_impl(v);

    if( ! m_manifold_info_initialized )
      return;

    // search for incident facets around v
    typedef std::vector<Facet> Facets;
    Facets facets;
    facets.reserve(64);
    r_tr_.incident_facets (v, std::back_inserter(facets));

    // foreach f in star(v)
    for (typename Facets::iterator fit = facets.begin();
         fit != facets.end();
         ++fit) 
    {
      // foreach edge of *fit
      const Cell_handle cell = fit->first;
      const int i = fit->second;
      for(int j = 0; j < 3; ++j)
      {
        const int edge_index_va = r_tr_.vertex_triple_index(i, j);
        const int edge_index_vb = r_tr_.vertex_triple_index(i, (j == 2) ? 0 : (j+1));
        Edge edge(cell, edge_index_va, edge_index_vb);
        // test if edge is in Complex
        if ( r_c3t3_.face_status(edge) != C3t3::NOT_IN_COMPLEX ) {
          // test if edge is not regular to store it as a "bad_edge"
          // e.g. more than or equal to 3 incident facets (SINGULAR)
          // or less than or equal to 1
          // (BOUNDARY only, because ISOLATED is NA)
          // This test is not efficient because
          // edges are tried to be inserted several times
          // TODO one day: test if the edge is still singular
          if ( (r_c3t3_.face_status(edge) == C3t3::SINGULAR) ||
               ( (!m_with_boundary) && 
                 (r_c3t3_.face_status(edge) == C3t3::BOUNDARY) )
               )
          {
            m_bad_edges.insert( edge_to_edgevv(edge) );
          }
          else {
            m_bad_edges.erase( edge_to_edgevv(edge) ); // @TODO: pourquoi?!
          }
        }
      }
    }

    // foreach v' in star of v
    std::vector<Vertex_handle> vertices;
    vertices.reserve(64);
    r_tr_.incident_vertices(v, std::back_inserter(vertices));

    // is_regular_or_boundary_for_vertices
    // is used here also incident edges are not known to be
    // REGULAR which may cause some singular vertices to be forgotten
    // This causes no problem because
    // those SINGULAR incident SINGULAR edges are going to be handled
    for (typename std::vector<Vertex_handle>::iterator
           vit = vertices.begin(), end = vertices.end();
         vit != end; ++vit)
    {
      if ( r_c3t3_.is_in_complex(*vit)  &&
           !r_c3t3_.is_regular_or_boundary_for_vertices(*vit))
      {
        m_bad_vertices.insert(*vit);
      }
    }

    if ( r_c3t3_.is_in_complex(v) &&
         !r_c3t3_.is_regular_or_boundary_for_vertices(v))
    {
      m_bad_vertices.insert(v);
    }
  }

  std::string debug_info() const
  {
    std::stringstream s;
    s << Base::debug_info()
      << "," << m_bad_edges.size()
      << "," << m_bad_vertices.size();
    return s.str();
  }

  std::string debug_info_header() const
  {
    std::stringstream s;
    s << Base::debug_info_header()
      << ",#bad edges,#bad vertices";
    return s.str();
  }
};  // end Refine_facets_manifold_base

}  // end namespace Mesh_3

}  // end namespace CGAL


#endif // CGAL_MESH_3_REFINE_FACETS_MANIFOLD_BASE_H
