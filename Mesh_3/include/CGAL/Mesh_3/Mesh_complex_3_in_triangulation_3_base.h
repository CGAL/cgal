// Copyright (c) 2003-2009  INRIA Sophia-Antipolis (France).
// Copyright (c) 2013       GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau, Stéphane Tayeb
//
//******************************************************************************
// File Description : Implements class Mesh_complex_3_in_triangulation_3.
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H
#define CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Mesh_3/utilities.h>
#include <CGAL/iterator.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/IO/File_maya.h>
#include <CGAL/Bbox_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/Mesh_3/io_signature.h>
#include <CGAL/Union_find.h>

#include <boost/functional/hash.hpp>

#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/atomic.h>
  #include <tbb/concurrent_hash_map.h>

namespace CGAL {
  template < class DSC, bool Const >
  std::size_t tbb_hasher(const CGAL::internal::CC_iterator<DSC, Const>& it)
  {
    return CGAL::internal::hash_value(it);
  }


  template < class DSC, bool Const >
  std::size_t tbb_hasher(const CGAL::CCC_internal::CCC_iterator<DSC, Const>& it)
  {
    return CGAL::CCC_internal::hash_value(it);
  }

  // As Marc Glisse pointed out the TBB hash of a std::pair is
  // simplistic and leads to the
  // TBB Warning: Performance is not optimal because the hash function
  //              produces bad randomness in lower bits in class
  //              tbb::interface5::concurrent_hash_map
  template < class DSC, bool Const >
  std::size_t tbb_hasher(const std::pair<CGAL::internal::CC_iterator<DSC, Const>,
                                         CGAL::internal::CC_iterator<DSC, Const> >& p)
  {
    return boost::hash<std::pair<CGAL::internal::CC_iterator<DSC, Const>,
                                 CGAL::internal::CC_iterator<DSC, Const> > >()(p);
  }


  template < class DSC, bool Const >
  std::size_t tbb_hasher(const std::pair<CGAL::CCC_internal::CCC_iterator<DSC, Const>,
                                         CGAL::CCC_internal::CCC_iterator<DSC, Const> >& p)
  {
    return boost::hash<std::pair<CGAL::CCC_internal::CCC_iterator<DSC, Const>,
                                 CGAL::CCC_internal::CCC_iterator<DSC, Const> > >()(p);
  }

}
#endif

namespace CGAL {
namespace Mesh_3 {

  namespace details {

    template <typename Tr>
    class C3t3_helper_class
    {
    protected:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle   Cell_handle;
      typedef typename Tr::Facet         Facet;
      typedef typename Tr::Edge          Edge;

      typedef std::pair<Vertex_handle, Vertex_handle> Pair_of_vertices;

      // computes and return an ordered pair of Vertex
      Pair_of_vertices
      make_ordered_pair(const Vertex_handle vh1, const Vertex_handle vh2) const {
	if (vh1 < vh2) {
	  return std::make_pair(vh1, vh2);
	}
	else {
	  return std::make_pair(vh2, vh1);
	}
      }

      // same from an Edge
      Pair_of_vertices
      make_ordered_pair(const Edge e) const {
        return make_ordered_pair(e.first->vertex(e.second),
                                 e.first->vertex(e.third));
      }

      Facet canonical_facet(Cell_handle c, int i) const {
        Cell_handle c2 = c->neighbor(i);
        return (c2 < c) ? std::make_pair(c2,c2->index(c)) : std::make_pair(c,i);
      }

    }; // end class template C3t3_helper_class

  } // end namespace Mesh_3::details

/**
 * @class Mesh_complex_3_in_triangulation_3_base
 * @brief A data-structure to represent and maintain a 3D complex embedded
 * in a 3D triangulation.
 */
template<typename Tr, typename Concurrency_tag>
class Mesh_complex_3_in_triangulation_3_base
  : public details::C3t3_helper_class<Tr>
{
  typedef Mesh_complex_3_in_triangulation_3_base<Tr, Concurrency_tag> Self;
  typedef details::C3t3_helper_class<Tr> Base;

public:
  // Triangulation types
  typedef Tr                            Triangulation;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Facet            Facet;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::size_type        size_type;

  // Indices types
  typedef typename Tr::Cell::Subdomain_index      Subdomain_index;
  typedef typename Tr::Cell::Surface_patch_index  Surface_patch_index;
  typedef typename Tr::Vertex::Index              Index;

  enum Face_status{ NOT_IN_COMPLEX = 0,
                    ISOLATED = 1, // - An ISOLATED edge is a marked edge,
                                  //   without any incident facets.
                    BOUNDARY,     // - An edge is on BOUNDARY if it has only
                                  //   one incident facet.
                                  // - A vertex is on BOUNDARY if all its
                                  //   incident edges are REGULAR or on
                                  //   BOUNDARY, at least one is on
                                  //   BOUNDARY, and the incident facets
                                  //   form only one connected component.
                    REGULAR,      // - A facet that is in the complex is
                                  //   REGULAR.
                                  // - An edge is REGULAR if it has
                                  //   exactly two incident facets.
                                  // - A vertex is REGULAR if all it
                                  //   incident edges are REGULAR, and the
                                  //   incident facets form only one
                                  //   connected component.
                    SINGULAR};    // - SINGULAR is for all other cases.

  //-------------------------------------------------------
  // Constructors / Destructors
  //-------------------------------------------------------
  /**
   * @brief Constructor
   * Builds an empty 3D complex.
   */
  Mesh_complex_3_in_triangulation_3_base()
    : Base()
    , tr_()
    , edge_facet_counter_() //TODO: parallel!
    , manifold_info_initialized_(false) //TODO: parallel!
  {
    // We don't put it in the initialization list because
    // tbb::atomic has no contructors
    number_of_facets_ = 0;
    number_of_cells_ = 0;
  }

  /// Copy constructor
  Mesh_complex_3_in_triangulation_3_base(const Self& rhs)
    : Base()
    , tr_(rhs.tr_)
    , edge_facet_counter_(rhs.edge_facet_counter_)
    , manifold_info_initialized_(rhs.manifold_info_initialized_)
  {
    number_of_facets_ = rhs.number_of_facets_;
    number_of_cells_ = rhs.number_of_cells_;
  }

  /// Destructor
  ~Mesh_complex_3_in_triangulation_3_base() {}

  void clear() {
    number_of_cells_ = 0;
    number_of_facets_ = 0;
    clear_manifold_info();
    tr_.clear();
  }

  /// Assignment operator
  Self& operator=(Self rhs)
  {
    swap(rhs);
    return *this;
  }

  /// Returns the reference to the triangulation
  Triangulation& triangulation() { return tr_; }
  /// Returns a const reference to the triangulation
  const Triangulation& triangulation() const { return tr_; }


  /// Adds facet \c facet to the 2D complex, with surface index \c index
  void add_to_complex(const Facet& facet, const Surface_patch_index& index)
  {
    add_to_complex(facet.first, facet.second, index);
  }

  /// Adds facet(\c cell, \c i) to the 2D complex, with surface index \c index
  void add_to_complex(const Cell_handle& cell,
                      const int i,
                      const Surface_patch_index& index);

  /// Removes facet \c facet from 2D complex
  void remove_from_complex(const Facet& facet);

  /// Removes facet(\c cell, \c i) from 2D complex
  void remove_from_complex(const Cell_handle& c, const int i) {
    remove_from_complex(Facet(c, i));
  }

  /// Sets surface index of facet \c facet to \c index
  void set_surface_patch_index(const Facet& f, const Surface_patch_index& index)
  {
    set_surface_patch_index(f.first, f.second, index);
  }

  /// Sets surface index of facet(\c cell, \c i) to \c index
  void set_surface_patch_index(const Cell_handle& cell,
                         const int i,
                         const Surface_patch_index& index) const
  {
    cell->set_surface_patch_index(i, index);
  }

  /// Returns `NOT_IN_COMPLEX`, `BOUNDARY`, `REGULAR`, or `SINGULAR`,
  /// depending on the number of incident facets in the complex, and the
  /// number of connected components of its link
  Face_status face_status(const Vertex_handle v) const
  {
    if(!manifold_info_initialized_) init_manifold_info();
    const std::size_t n = v->cached_number_of_incident_facets();

    if(n == 0) return NOT_IN_COMPLEX;

    //test incident edges for REGULARITY and count BOUNDARY edges
    typename std::vector<Edge> edges;
    edges.reserve(64);
    if(tr_.is_parallel()) {
      tr_.incident_edges_threadsafe(v, std::back_inserter(edges));
    } else {
      tr_.incident_edges(v, std::back_inserter(edges));
    }
    int number_of_boundary_incident_edges = 0; // could be a bool
    for (typename std::vector<Edge>::iterator
           eit=edges.begin(), end = edges.end();
	 eit != end; eit++)
    {
      switch( face_status(*eit) )
      {
      case NOT_IN_COMPLEX: case REGULAR: break;
      case BOUNDARY: ++number_of_boundary_incident_edges; break;
      default :
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
        std::cerr << "singular edge...\n";
        std::cerr << v->point() << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
        return SINGULAR;
      }
    }

    // From here all incident edges (in complex) are REGULAR or BOUNDARY.
    const std::size_t nb_components = union_find_of_incident_facets(v);
    if(nb_components > 1) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      std::cerr << "singular vertex: nb_components=" << nb_components << std::endl;
      std::cerr << v->point() << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      return SINGULAR;
    }
    else { // REGULAR OR BOUNDARY
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      std::cerr << "regular or boundary: " << v->point() << std::endl;
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
      if (number_of_boundary_incident_edges != 0)
        return BOUNDARY;
      else
        return REGULAR;
    }
  }

  /// This function should be called only when incident edges
  /// are known to be REGULAR OR BOUNDARY
  bool is_regular_or_boundary_for_vertices(Vertex_handle v) const {
    return union_find_of_incident_facets(v) == 1;
  }

  /// Returns `NOT_IN_COMPLEX`, `BOUNDARY`, `REGULAR`, or `SINGULAR`,
  /// depending on the number of incident facets in the complex
  Face_status face_status(const Edge& edge) const
  {
    if(!manifold_info_initialized_) init_manifold_info();

#ifdef CGAL_LINKED_WITH_TBB
    typename Edge_facet_counter::const_accessor accessor;
    if(!edge_facet_counter_.find(accessor,
				 this->make_ordered_pair(edge)))
      return NOT_IN_COMPLEX;
    switch(accessor->second)
#else // not CGAL_LINKED_WITH_TBB
    switch(edge_facet_counter_[this->make_ordered_pair(edge)])
#endif // not CGAL_LINKED_WITH_TBB
    {
    case 0: return NOT_IN_COMPLEX;
    case 1: return BOUNDARY;
    case 2: return REGULAR;
    default: return SINGULAR;
    }
  }

  /// Returns true if the vertex \c v has is incident to at least a facet
  /// of the complex
  bool has_incident_facets_in_complex(const Vertex_handle& v) const
  {
    if(!manifold_info_initialized_) init_manifold_info();
    return v->cached_number_of_incident_facets() > 0;
  }

  /// Returns true if facet \c facet is in complex
  bool is_in_complex(const Facet& facet) const
  {
    return is_in_complex(facet.first, facet.second);
  }

  /// Returns true if facet (\c cell, \c i) is in 2D complex
  bool is_in_complex(const Cell_handle& cell, const int i) const
  {
    return ( cell->is_facet_on_surface(i) );
  }

  /// Returns surface index of facet \c f
  Surface_patch_index surface_patch_index(const Facet& f) const
  {
    return surface_patch_index(f.first,f.second);
  }

  /// Returns surface index of facet(\c cell, \c i)
  Surface_patch_index surface_patch_index(const Cell_handle& cell,
                                          const int i) const
  {
    return cell->surface_patch_index(i);
  }

  /// Returns the number of surface facets of the mesh
  size_type number_of_facets_in_complex() const { return number_of_facets_; }

  /// Adds cell \c cell to the 3D complex, with subdomain index \c index
  void add_to_complex(const Cell_handle& cell, const Subdomain_index& index)
  {
    CGAL_precondition( !( index == Subdomain_index() ) );

    if ( ! is_in_complex(cell) )
    {
      set_subdomain_index(cell, index);
      ++number_of_cells_;
    }
  }

  /// Removes cell \c cell from the 3D complex
  void remove_from_complex(const Cell_handle& cell)
  {
    if ( is_in_complex(cell) )
    {
      set_subdomain_index(cell, Subdomain_index());
      --number_of_cells_;
    }
  }


  /// Sets subdomain index of cell \c cell to \c index
  void set_subdomain_index(const Cell_handle& cell,
                           const Subdomain_index& index) const
  {
    cell->set_subdomain_index(index);
  }

  /// Sets index of vertex \c vertex to \c index
  void set_index(const Vertex_handle& vertex, const Index& index) const
  {
    vertex->set_index(index);
  }

  /// Sets dimension of vertex \c vertex to \c dimension
  void set_dimension(const Vertex_handle& vertex, int dimension) const
  {
    vertex->set_dimension(dimension);
  }

  /// Returns the number of cells which belongs to the 3D complex
  size_type number_of_cells_in_complex() const
  {
    return number_of_cells_;
  }

  /// Returns \c true if cell \c cell belongs to the 3D complex
  bool is_in_complex(const Cell_handle& cell) const
  {
    return !( subdomain_index(cell) == Subdomain_index() );
  }

  /// Returns the subdomain index of cell \c cell
  Subdomain_index subdomain_index(const Cell_handle& cell) const
  {
    return cell->subdomain_index();
  }

  /// Returns the dimension of the lowest dimensional face of the input 3D
  /// complex that contains the vertex
  int in_dimension(const Vertex_handle& v) const { return v->in_dimension(); }

  /// Returns the index of vertex \c v
  Index index(const Vertex_handle& v) const { return v->index(); }

  /// Outputs the mesh to medit
  void output_to_medit(std::ostream& os,
                       bool rebind = true,
                       bool show_patches = false) const
  {
    // Call global function
    CGAL::output_to_medit(os,*this,rebind,show_patches);
  }
  
  /// Outputs the mesh to maya
  void output_to_maya(std::ofstream& os,
                      bool surfaceOnly = true) const
  {
    // Call global function
    CGAL::output_to_maya(os,*this,surfaceOnly);
  }

  //-------------------------------------------------------
  // Undocumented features
  //-------------------------------------------------------
  /**
   * @brief insert \c [first,last[ in the triangulation (with dimension 2)
   * @param first the iterator on the first point to insert
   * @param last the iterator past the last point to insert
   *
   * InputIterator value type must be \c std::pair<Tr::Point,Index>
   */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first, InputIterator last)
  {
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp =
      tr_.geom_traits().construct_weighted_point_3_object();

    while ( first != last )
    {
      Vertex_handle vertex = tr_.insert(cwp((*first).first));
      vertex->set_index((*first).second);
      vertex->set_dimension(2);
      ++first;
    }
  }

  /**
   * @brief insert \c [first,last[ in the triangulation (with dimension 2 and
   * index \c default_index)
   * @param first the iterator on the first point to insert
   * @param last the iterator past the last point to insert
   * @param default_index the index to be used to insert points
   *
   * InputIterator value type must be \c Tr::Point
   */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first,
                             InputIterator last,
                             const Index& default_index)
  {
    typename Tr::Geom_traits::Construct_weighted_point_3 cwp =
      tr_.geom_traits().construct_weighted_point_3_object();

    while ( first != last )
    {
      Vertex_handle vertex = tr_.insert(cwp(*first));
      vertex->set_index(default_index);
      vertex->set_dimension(2);
      ++first;
    }
  }

  /// Swaps this & rhs
  void swap(Self& rhs)
  {
    std::swap(rhs.number_of_facets_, number_of_facets_);
    tr_.swap(rhs.tr_);
    std::swap(rhs.number_of_cells_, number_of_cells_);
  }

  /// Returns bbox
  Bbox_3 bbox() const;

  void clear_cells_and_facets_from_c3t3() {
    for(typename Tr::Finite_cells_iterator
          cit = this->triangulation().finite_cells_begin(),
          end = this->triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      set_subdomain_index(cit, Subdomain_index());
    }
    this->number_of_cells_ = 0;
    for(typename Tr::Finite_facets_iterator
          fit = this->triangulation().finite_facets_begin(),
          end = this->triangulation().finite_facets_end();
        fit != end; ++fit)
    {
      Facet facet = *fit;
      Facet mirror = tr_.mirror_facet(facet);
      set_surface_patch_index(facet.first, facet.second, Surface_patch_index());
      set_surface_patch_index(mirror.first, mirror.second, Surface_patch_index());
    }
    this->number_of_facets_ = 0;
    clear_manifold_info();
  }

  void clear_manifold_info() {
    edge_facet_counter_.clear();
    manifold_info_initialized_ = false;
  }

private:
  void init_manifold_info() const {
    for(typename Tr::All_vertices_iterator
          vit = triangulation().finite_vertices_begin(),
          end = triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      vit->set_c2t3_cache(0, -1);
    }

    edge_facet_counter_.clear();

    for(typename Tr::Finite_facets_iterator
          fit = triangulation().finite_facets_begin(),
          end = triangulation().finite_facets_end();
        fit != end; ++fit)
    {
      if ( is_in_complex(*fit) ) {
        const Cell_handle cell = fit->first;
        const int i = fit->second;
        for(int j = 0; j < 3; ++j)
        {
          const int edge_index_va = tr_.vertex_triple_index(i, j);
          const int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j+1));
          const Vertex_handle edge_va = cell->vertex(edge_index_va);
          const Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifndef CGAL_LINKED_WITH_TBB
          ++edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#else // CGAL_LINKED_WITH_TBB
	  {
	    typename Edge_facet_counter::accessor accessor;
	    edge_facet_counter_.insert(accessor,
				       this->make_ordered_pair(edge_va, edge_vb));
	    ++accessor->second;
	  }
#endif // CGAL_LINKED_WITH_TBB

          const std::size_t n = edge_va->cached_number_of_incident_facets();
          edge_va->set_c2t3_cache(n+1, -1);
        }
      }
    }
    manifold_info_initialized_ = true;
  }

  /// Extract the subset `F` of facets of the complex incident to `v` and
  /// return the number of connected component of the adjacency graph of `F`.
  std::size_t union_find_of_incident_facets(const Vertex_handle v) const
  {
    if( v->is_c2t3_cache_valid() )
    {
      const std::size_t n = v->cached_number_of_components();
      if(n != std::size_t(-1)) return n;
    }

    Union_find<Facet> facets;
    { // fill the union find
      std::vector<Facet> non_filtered_facets;
      if(tr_.is_parallel()) {
	tr_.incident_facets_threadsafe(v, std::back_inserter(non_filtered_facets));
      } else {
	tr_.incident_facets(v, std::back_inserter(non_filtered_facets));
      }

      for(typename std::vector<Facet>::iterator
            fit = non_filtered_facets.begin(),
            end = non_filtered_facets.end();
          fit != end; ++fit)
      {
        if(is_in_complex(*fit)) facets.push_back(*fit);
      }
    }

    typedef std::map<Vertex_handle,
                     typename Union_find<Facet>::handle> Vertex_set_map;
    typedef typename Vertex_set_map::iterator Vertex_set_map_iterator;

    Vertex_set_map vsmap;

    for(typename Union_find<Facet>::iterator
          it = facets.begin(), end = facets.end();
        it != end; ++it)
    {
      const Cell_handle& ch = (*it).first;
      const int& i = (*it).second;
      for(int j=0; j < 3; ++j) {
	const Vertex_handle w = ch->vertex(tr_.vertex_triple_index(i,j));
	if(w != v){
	  Vertex_set_map_iterator vsm_it = vsmap.find(w);
	  if(vsm_it != vsmap.end()){
	    facets.unify_sets(vsm_it->second, it);
	  } else {
	    vsmap.insert(std::make_pair(w, it));
	  }
	}
      }
    }
    const std::size_t nb_components = facets.number_of_sets();

    const std::size_t n = v->cached_number_of_incident_facets();
    v->set_c2t3_cache(n, nb_components);
    return nb_components;
  }
  
  //-------------------------------------------------------
  // Traversal
  //-------------------------------------------------------
private:
  typedef Mesh_3::internal::Iterator_not_in_complex<Self> Iterator_not_in_complex;

  class Facet_iterator_not_in_complex
  {
    const Self* c3t3_;
    Surface_patch_index index_; //need by SWIG: should be const Surface_patch_index
  public:
    Facet_iterator_not_in_complex(){} //need by SWIG
    Facet_iterator_not_in_complex(const Self& c3t3,
                                  const Surface_patch_index& index = Surface_patch_index())
      : c3t3_(&c3t3)
      , index_(index) { }

    template <typename Iterator>
    bool operator()(Iterator it) const
    {
      if ( index_ == Surface_patch_index() ) { return ! c3t3_->is_in_complex(*it); }
      else { return !( c3t3_->surface_patch_index(*it) == index_ );  }
    }
  };

  /**
   * @class Cell_not_in_complex
   * @brief A class to filter cells which do not belong to the complex
   */
  class Cell_not_in_complex
  {
    const Self* r_self_;
    Subdomain_index index_;//needed by SWIG, should be const Subdomain_index
  public:
    Cell_not_in_complex(){}//needed by SWIG
    Cell_not_in_complex(const Self& self,
                        const Subdomain_index& index = Subdomain_index())
      : r_self_(&self)
      , index_(index) { }

    bool operator()(Cell_handle ch) const
    {
      if ( index_ == Subdomain_index() ) { return !r_self_->is_in_complex(ch); }
      else { return !( r_self_->subdomain_index(ch) == index_ ); }
    }
  }; // end class Cell_not_in_complex

public:
  /// Iterator type to visit the facets of the 2D complex.
  typedef Filter_iterator<
    typename Triangulation::Finite_facets_iterator,
    Facet_iterator_not_in_complex >               Facets_in_complex_iterator;

  /// Returns a Facets_in_complex_iterator to the first facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this),
                                 tr_.finite_facets_begin());
  }

  /// Returns a Facets_in_complex_iterator to the first facet of the 2D complex
  Facets_in_complex_iterator
  facets_in_complex_begin(const Surface_patch_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this,index),
                                 tr_.finite_facets_begin());
  }

  /// Returns past-the-end iterator on facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_end(const Surface_patch_index = Surface_patch_index()) const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this));
  }

  /**
   * @class Cells_in_complex_iterator
   * @brief Iterator type to visit the cells of triangulation belonging
   * to the 3D complex
   *
   * This class is usefull to ensure that Cells_in_complex_iterator is convertible
   * to Cell_handle
   */
  class Cells_in_complex_iterator :
    public Filter_iterator<typename Triangulation::Finite_cells_iterator,
                           Cell_not_in_complex>
  {
  private:
    typedef typename Triangulation::Finite_cells_iterator Tr_iterator;
    typedef Filter_iterator<typename Triangulation::Finite_cells_iterator,
                            Cell_not_in_complex> Base;
    typedef Cells_in_complex_iterator Self;

  public:
    Cells_in_complex_iterator() : Base() { }
    Cells_in_complex_iterator(Base i) : Base(i) { }

    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    operator Cell_handle() const { return Cell_handle(this->base()); }
  }; // end class Cells_in_complex_iterator


  /// Returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator cells_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this),
                                 tr_.finite_cells_begin());
  }

  /// Returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator
  cells_in_complex_begin(const Subdomain_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this,index),
                                 tr_.finite_cells_begin());
  }

  /// Returns the past-the-end iterator for the cells of the 3D complex
  Cells_in_complex_iterator cells_in_complex_end() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this));
  }

  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  void set_surface_index(const Facet& f, const Surface_index& index)
  { set_surface_patch_index(f, index); }

  void set_surface_index(const Cell_handle& c, const int i, const Surface_index& index)
  { set_surface_patch_index(c,i,index); }

  Surface_index surface_index(const Facet& f) const
  { return surface_patch_index(f); }

  Surface_index surface_index(const Cell_handle& c, const int i) const
  { return surface_patch_index(c,i); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

#ifndef CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
  typedef Facets_in_complex_iterator  Facet_iterator;
  typedef Cells_in_complex_iterator   Cell_iterator;

  Facet_iterator facets_begin() const
  { return facets_in_complex_begin(); }

  Facet_iterator facets_end() const
  { return facets_in_complex_end(); }

  Cell_iterator cells_begin() const
  { return cells_in_complex_begin(); }

  Cell_iterator cells_end() const
  { return cells_in_complex_end(); }
#endif // CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

  size_type number_of_facets() const
  { return number_of_facets_in_complex(); }

  size_type number_of_cells() const
  { return number_of_cells_in_complex(); }

public:
  template <typename Tr2, typename Ct2>
  friend
  std::istream &
  operator>> (std::istream& is,
              Mesh_complex_3_in_triangulation_3_base<Tr2,Ct2> &c3t3);

  void rescan_after_load_of_triangulation();

  static
  std::string io_signature()
  {
    return
      Get_io_signature<Tr>()();
  }
private:

  // Sequential: non-atomic
  // "dummy" is here to allow the specialization (see below)
  // See http://groups.google.com/group/comp.lang.c++.moderated/browse_thread/thread/285ab1eec49e1cb6
  template<typename Concurrency_tag2, typename dummy = void>
  struct Number_of_elements
  {
    typedef size_type type;
  };
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel: atomic
  template<typename dummy>
  struct Number_of_elements<Parallel_tag, dummy>
  {
    typedef tbb::atomic<size_type> type;
  };
#endif // CGAL_LINKED_WITH_TBB

  // Private date members
  Triangulation tr_;

  typedef typename Base::Pair_of_vertices Pair_of_vertices;
#ifdef CGAL_LINKED_WITH_TBB
  typedef tbb::concurrent_hash_map<Pair_of_vertices, int> Edge_facet_counter;
#else // not CGAL_LINKED_WITH_TBB
  typedef std::map<Pair_of_vertices, int> Edge_facet_counter;
#endif // not CGAL_LINKED_WITH_TBB

  mutable Edge_facet_counter edge_facet_counter_;

  typename Number_of_elements<Concurrency_tag>::type number_of_facets_;
  typename Number_of_elements<Concurrency_tag>::type number_of_cells_;

  mutable bool manifold_info_initialized_;
};  // end class Mesh_complex_3_in_triangulation_3_base


template <typename Tr, typename Ct>
void
Mesh_complex_3_in_triangulation_3_base<Tr,Ct>::add_to_complex(
    const Cell_handle& cell,
    const int i,
    const Surface_patch_index& index)
{
  CGAL_precondition( !( index == Surface_patch_index() ) );

  if ( ! is_in_complex(cell,i) )
  {
    Facet mirror = tr_.mirror_facet(std::make_pair(cell,i));
    set_surface_patch_index(cell, i, index);
    set_surface_patch_index(mirror.first, mirror.second, index);
    ++number_of_facets_;
    if(manifold_info_initialized_) {
      for(int j = 0; j < 3; ++j)
      {
        int edge_index_va = tr_.vertex_triple_index(i, j);
        int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j+1));
        Vertex_handle edge_va = cell->vertex(edge_index_va);
        Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifdef CGAL_LINKED_WITH_TBB
	{
	  typename Edge_facet_counter::accessor accessor;
	  edge_facet_counter_.insert(accessor,
				     this->make_ordered_pair(edge_va, edge_vb));
	  ++accessor->second;
	}
#else // not CGAL_LINKED_WITH_TBB
        ++edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#endif // not CGAL_LINKED_WITH_TBB

        const std::size_t n = edge_va->cached_number_of_incident_facets();
        const std::size_t m = edge_va->cached_number_of_components();
        edge_va->set_c2t3_cache(n+1, m);
      }
      const int dimension_plus_1 = tr_.dimension() + 1;
      // update c2t3 for vertices of f
      for (int j = 0; j < dimension_plus_1; j++) {
        if (j != i) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          if(cell->vertex(j)->is_c2t3_cache_valid())
            std::cerr << "(" << cell->vertex(j)->point() << ")->invalidate_c2t3_cache()\n";
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          cell->vertex(j)->invalidate_c2t3_cache();
        }
      }
    }
  }
}


template <typename Tr, typename Ct>
void
Mesh_complex_3_in_triangulation_3_base<Tr,Ct>::remove_from_complex(const Facet& facet)
{
  if ( is_in_complex(facet) )
  {
    Facet mirror = tr_.mirror_facet(facet);
    set_surface_patch_index(facet.first, facet.second, Surface_patch_index());
    set_surface_patch_index(mirror.first, mirror.second, Surface_patch_index());
    --number_of_facets_;
    if(manifold_info_initialized_) {
      const Cell_handle cell = facet.first;
      const int i = facet.second;
      for(int j = 0; j < 3; ++j)
      {
        const int edge_index_va = tr_.vertex_triple_index(i, j);
        const int edge_index_vb = tr_.vertex_triple_index(i, (j == 2) ? 0 : (j+1));
        const Vertex_handle edge_va = cell->vertex(edge_index_va);
        const Vertex_handle edge_vb = cell->vertex(edge_index_vb);
#ifdef CGAL_LINKED_WITH_TBB
	{
	  typename Edge_facet_counter::accessor accessor;
	  edge_facet_counter_.insert(accessor,
				     this->make_ordered_pair(edge_va, edge_vb));
	  --accessor->second;
	}
#else // not CGAL_LINKED_WITH_TBB
        --edge_facet_counter_[this->make_ordered_pair(edge_va, edge_vb)];
#endif // not CGAL_LINKED_WITH_TBB

        const std::size_t n = edge_va->cached_number_of_incident_facets();
        CGAL_assertion(n>0);
        const std::size_t m = edge_va->cached_number_of_components();
        edge_va->set_c2t3_cache(n-1, m);
      }
      const int dimension_plus_1 = tr_.dimension() + 1;
      // update c2t3 for vertices of f
      for (int j = 0; j < dimension_plus_1; j++) {
        if (j != facet.second) {
#ifdef CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          if(cell->vertex(j)->is_c2t3_cache_valid())
            std::cerr << "(" << cell->vertex(j)->point() << ")->invalidate_c2t3_cache()\n";
#endif // CGAL_MESHES_DEBUG_REFINEMENT_POINTS
          cell->vertex(j)->invalidate_c2t3_cache();
        }
      }
    }
  }
}


// -----------------------------------
// Undocumented
// -----------------------------------
template <typename Tr, typename Ct>
Bbox_3
Mesh_complex_3_in_triangulation_3_base<Tr,Ct>::
bbox() const
{
  if ( 0 == triangulation().number_of_vertices() )
  {
    return Bbox_3();
  }

  typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
  Bbox_3 result = (vit++)->point().bbox();

  for(typename Tr::Finite_vertices_iterator end = tr_.finite_vertices_end();
      vit != end ; ++vit)
  {
    result = result + vit->point().bbox();
  }

  return result;
}

template <typename Tr, typename Ct>
std::ostream &
operator<< (std::ostream& os,
            const Mesh_complex_3_in_triangulation_3_base<Tr,Ct> &c3t3)
{
  return os << c3t3.triangulation();
}


template <typename Tr, typename Ct>
std::istream &
operator>> (std::istream& is,
            Mesh_complex_3_in_triangulation_3_base<Tr,Ct> &c3t3)
{
  c3t3.clear();
  is >> c3t3.triangulation();

  if(!is) {
    c3t3.clear();
    return is;
  }

  c3t3.rescan_after_load_of_triangulation();
  return is;
}

template <typename Tr, typename Ct>
void
Mesh_complex_3_in_triangulation_3_base<Tr,Ct>::
rescan_after_load_of_triangulation() {
  this->number_of_facets_ = 0;
  for(typename Tr::Finite_facets_iterator
        fit = this->triangulation().finite_facets_begin(),
        end = this->triangulation().finite_facets_end();
      fit != end; ++fit)
  {
    if ( this->is_in_complex(*fit) ) {
      ++this->number_of_facets_;
    }
  }

  this->number_of_cells_ = 0;
  for(typename Tr::Finite_cells_iterator
        cit = this->triangulation().finite_cells_begin(),
        end = this->triangulation().finite_cells_end();
      cit != end; ++cit)
  {
    if ( this->is_in_complex(cit) ) {
      ++this->number_of_cells_;
    }
  }
}

}  // end namespace Mesh_3
}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H
