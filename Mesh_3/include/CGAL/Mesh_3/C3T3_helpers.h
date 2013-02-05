// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_C3T3_HELPERS_H
#define CGAL_MESH_3_C3T3_HELPERS_H

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/tuple.h>
#include <CGAL/iterator.h>

#include <boost/optional.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <functional>
#include <vector>
#include <set>

namespace CGAL {
namespace Mesh_3 {
  
template <typename C3T3, typename MeshDomain>
class C3T3_helpers
{
  // -----------------------------------
  // Private types
  // -----------------------------------
  typedef typename C3T3::Triangulation  Tr;
  typedef Tr                            Triangulation;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Gt::Vector_3         Vector_3;
  typedef typename Gt::Point_3          Point_3;
  typedef typename Gt::Plane_3          Plane_3;
  typedef typename Gt::FT               FT;
  
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Cell             Cell;
  typedef typename Tr::Facet            Facet;
  
  typedef typename C3T3::Surface_patch_index  Surface_patch_index;
  typedef typename C3T3::Subdomain_index      Subdomain_index;
  typedef typename C3T3::Index                Index;
  
  typedef std::vector<Cell_handle>      Cell_vector;
  typedef std::set<Cell_handle>         Cell_set;
  typedef std::vector<Facet>            Facet_vector;
  typedef std::vector<Vertex_handle>    Vertex_vector;
  typedef std::set<Vertex_handle>       Vertex_set;
  
  // Facet_boundary stores the boundary of surface facets
  typedef std::pair<Vertex_handle,Vertex_handle> Ordered_edge;
  typedef std::set<std::pair<Ordered_edge,Surface_patch_index> >  Facet_boundary;
  
  typedef Triangulation_helpers<Tr> Th;
  
public:
  // -----------------------------------
  // Public interface
  // -----------------------------------
  typedef boost::optional<Vertex_handle> Update_mesh;
  
  /**
   * Constructor
   */
  C3T3_helpers(C3T3& c3t3, const MeshDomain& domain)
    : c3t3_(c3t3)
    , tr_(c3t3.triangulation())
    , domain_(domain) { }
  
  /**
   * @brief tries to move \c old_vertex to \c new_point in the mesh
   * @param new_point the new location of \c old_vertex
   * @param old_vertex the old vertex
   * @param criterion the criterion which will be used to verify the new 
   *    location is ok. c3t3 minimal value of new criterion shall not decrease.
   * @param modified_vertices contains the vertices incident to cells which 
   *    may have been impacted by relocation
   * @return a pair which contains:
   *    - a bool which is \c true if the move has been done.
   *    - a Vertex_handle which is always filled and may be the new vertex (if
   *      the move is a success), or the vertex which lies at \c v's position in
   *      the updated c3t3.
   */
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh(const Point_3& new_point,
              const Vertex_handle& old_vertex,
              const SliverCriterion& criterion,
              OutputIterator modified_vertices);
  
  /**
   * Updates mesh moving vertex \c old_vertex to \c new_point. Returns the
   * new vertex of the triangulation.
   *
   * Insert into modified vertices the vertices which are impacted by to move.
   */
  template <typename OutputIterator>
  Vertex_handle update_mesh(const Point_3& new_point,
                            const Vertex_handle& old_vertex,
                            OutputIterator modified_vertices,
                            bool fill_modified_vertices = true);
  
  /**
   * Updates mesh moving vertex \c old_vertex to \c new_point. Returns the
   * new vertex of the triangulation.
   */
  Vertex_handle update_mesh(const Point_3& new_point,
                            const Vertex_handle& old_vertex)
  {
    std::vector<Vertex_handle> dummy;
    return update_mesh(new_point, old_vertex, std::back_inserter(dummy), false);
  }
  
  /**
   * Rebuilds restricted Delaunay
   */
  template <typename ForwardIterator>
  void rebuild_restricted_delaunay(ForwardIterator first_cell,
                                   ForwardIterator last_cell,
                                   Vertex_set& moving_vertices);
  
  /**
   * @brief Project \c p on surface, using incident facets of \c v
   * @param p The point to project
   * @param v The vertex from which p was moved
   * @return the projected point
   *
   * \c p is projected as follows using normal of least square fitting plane
   * on \c v incident surface points.
   */
  Point_3
  project_on_surface(const Point_3& p, const Vertex_handle& v) const;

  /**
   * Returns the minimum value for criterion for incident cells of \c vh 
   */
  template <typename SliverCriterion> 
  FT min_incident_value(const Vertex_handle& vh,
                        const SliverCriterion& criterion) const;
  
  /**
   * Moves \c old_vertex to \c new_location
   * Stores the cells which have to be updated in \c outdated_cells
   */
  Vertex_handle move_point(const Vertex_handle& old_vertex,
                           const Point_3& new_location,
                           Cell_set& outdated_cells);
  
  /**
   * Outputs to out the sliver (wrt \c criterion and \c sliver_bound) incident
   * to \c v
   */
  template <typename SliverCriterion, typename OutputIterator>
  OutputIterator
  incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound,
                   OutputIterator out) const;
  
  /**
   * Returns the sliver (wrt \c criterion and \c sliver_bound) incident to \c v
   */ 
  template <typename SliverCriterion> 
  Cell_vector
  incident_slivers(const Vertex_handle& v,
                   const SliverCriterion& criterion,
                   const FT& sliver_bound) const
  {
    Cell_vector slivers;
    incident_slivers(v, criterion, sliver_bound, std::back_inserter(slivers));
    return slivers;
  }
  
  /**
   * Returns the minimum criterion value of cells contained in \c cells
   * Precondition: cells of \c cells must not be infinite.
   * Warning: Here we don't check if cells are in c3t3
   */
  template <typename SliverCriterion>
  FT min_sliver_value(const Cell_vector& cells,
                      const SliverCriterion& criterion,
                      const bool use_cache = true) const;  
  
  /**
   * Reset cache validity of all cells of c3t3_
   */
  void reset_cache() const
  {
    namespace bl = boost::lambda;
    std::for_each(c3t3_.cells_in_complex_begin(),c3t3_.cells_in_complex_end(),
                  bl::bind(&Cell::reset_cache_validity, bl::_1) );
  }
  
private:  
  // -----------------------------------
  // Usefull Functors
  // -----------------------------------
  /**
   * @class Get_all_facets
   *
   * A functor which adds to an output iterator canonical facets of a cell
   */
  template <typename OutputIterator>
  class Get_all_facets
  {
  public:
    Get_all_facets(const Triangulation& tr, OutputIterator out)
      : tr_(tr)
      , out_(out) {}
    
    void operator()(const Cell_handle& cell)
    {
      for ( int i=0 ; i<4 ; ++i )
        if ( !tr_.is_infinite(cell,i) )
          *out_++ = canonical_facet(cell,i);
    }
    
  private:
    Facet canonical_facet(const Cell_handle& c, const int i) const
    {
      Facet facet(c,i);
      Facet mirror = tr_.mirror_facet(facet);
      return ( (mirror<facet)?mirror:facet );
    }
    
  private:
    const Triangulation& tr_;
    OutputIterator out_;
  };
  
  
  /**
   * @class Is_in_c3t3
   *
   * A functor which returns true if a given handle is in c3t3
   */
  template <typename Handle>
  class Is_in_c3t3 : public std::unary_function<Handle, bool>
  {
  public:
    Is_in_c3t3(const C3T3& c3t3) : c3t3_(c3t3) { }
    bool operator()(const Handle& h) const { return c3t3_.is_in_complex(h); }
    
  private:
    const C3T3& c3t3_;
  };
  

  /**
   * @class Is_sliver
   *
   * A functor which answers true if a Cell_handle is a sliver
   */
  template <typename SliverCriterion>
  struct Is_sliver : public std::unary_function<Cell_handle,bool>
  {
    Is_sliver(const C3T3& c3t3,
              const SliverCriterion& criterion,
              const FT& bound)
      : c3t3_(c3t3)
      , criterion_(criterion)
      , bound_(bound) { }
    
    bool operator()(const Cell_handle& c) const
    {
      if ( c3t3_.is_in_complex(c) )
      {
        CGAL_assertion(!c3t3_.triangulation().is_infinite(c));
        
        if ( ! c->is_cache_valid() )
        {
          FT sliver_value = criterion_(c3t3_.triangulation().tetrahedron(c));
          c->set_sliver_value(sliver_value);
        }
        return ( c->sliver_value() <= bound_ );
      }
      else
        return false;
    }
    
  private:
    const C3T3& c3t3_;
    const SliverCriterion& criterion_;
    const FT bound_;
  };
  
  
  /**
   * @class Update_c3t3
   *
   * A functor which updates c3t3 w.r.t the domain.
   */
  class Update_c3t3
  {
  public:
    Update_c3t3(const MeshDomain& domain, C3T3& c3t3)
      : domain_(domain)
      , c3t3_(c3t3) {}
    
    /**
     * @brief Updates facet \c facet in c3t3
     * @param facet the facet to update
     * @param update if set to \c false, checking only is done
     * @return true if \c facet is in c3t3
     */
    bool operator()(const Facet& facet, const bool update = true) const
    {
      typedef typename C3T3::Triangulation::Geom_traits Gt;
      typedef typename Gt::Segment_3 Segment_3;
      typedef typename Gt::Ray_3 Ray_3;
      typedef typename Gt::Line_3 Line_3;
      
      // Nothing to do for infinite facets
      if ( c3t3_.triangulation().is_infinite(facet) )
        return false;
      
      // Functors
      typename Gt::Is_degenerate_3 is_degenerate = 
        Gt().is_degenerate_3_object();
      
      // Get dual of facet
      Object dual = c3t3_.triangulation().dual(facet);

      // The dual is a segment, a ray or a line
      if ( const Segment_3* p_segment = object_cast<Segment_3>(&dual) )
      {
        if (is_degenerate(*p_segment)) 
          return false;
        
        return dual_intersect(*p_segment,facet,update);
      }
      else if ( const Ray_3* p_ray = object_cast<Ray_3>(&dual) )
      {
        if (is_degenerate(*p_ray))
          return false;
        
        return dual_intersect(*p_ray,facet,update);
      }
      else if ( const Line_3* p_line = object_cast<Line_3>(&dual) )
      {
        return dual_intersect(*p_line,facet,update);
      }
      
      // Should not happen
      CGAL_assertion(false);
      return false;
    }
    
    /**
     * @brief Updates cell \c ch in c3t3
     * @param ch the cell to update
     * @param update if set to \c false, checking only is done
     * @return true if \c ch is in c3t3
     */
    bool operator()(const Cell_handle& ch, const bool update = true) const
    {
      typedef boost::optional<typename MeshDomain::Subdomain_index> Subdomain;
      
      if ( c3t3_.triangulation().is_infinite(ch) )
        return false;
      
      // treat cell
      const Subdomain subdomain =
        domain_.is_in_domain_object()(c3t3_.triangulation().dual(ch));
      
      if ( subdomain && update )
      {
        c3t3_.add_to_complex(ch,*subdomain);
      }
      
      return subdomain;
    }
    
  private:
    
    // Returns true if query intersects the surface.
    template <typename Query>
    bool dual_intersect(const Query& dual,
                        const Facet& facet,
                        const bool update) const
    {
      typedef boost::optional<typename MeshDomain::Surface_patch_index> Surface_patch;
      
      typename MeshDomain::Do_intersect_surface do_intersect_surface =
        domain_.do_intersect_surface_object();
      
      Surface_patch surface = do_intersect_surface(dual);
      
      // Update if needed
      if ( surface && update )
      {
        // Update facet surface center
        typename MeshDomain::Construct_intersection intersection =
          domain_.construct_intersection_object();
        
        Point_3 surface_center = CGAL::cpp11::get<0>(intersection(dual));
        facet.first->set_facet_surface_center(facet.second,surface_center);
        
        // Update status in c3t3 
        c3t3_.add_to_complex(facet,*surface);          
      }
      
      return surface;
    }
    
    
  private:
    const MeshDomain& domain_;
    C3T3& c3t3_;
  };
  
  
  /**
   * @class Sliver_criterion_value
   *
   * A functor which returns sliver criterion value for a Cell_handle 
   */
  template <typename SliverCriterion>
  class Sliver_criterion_value
    : public std::unary_function<Cell_handle, FT>
  {
  public:
    Sliver_criterion_value(const Tr& tr,
                           const SliverCriterion& criterion)
      : p_tr_(&tr)
      , criterion_(criterion) {}
    
    FT operator()(const Cell_handle& ch) const
    {
      CGAL_precondition(!p_tr_->is_infinite(ch));
      
      if ( ! ch->is_cache_valid() )
      {
        FT sliver_value = criterion_(p_tr_->tetrahedron(ch));
        ch->set_sliver_value(sliver_value);
      }
      return ch->sliver_value();
    }
    
  private:
    // '=' is used, so p_tr_ must be a pointer ...
    const Tr* p_tr_;
    SliverCriterion criterion_;
  };
  
private:
  // -----------------------------------
  // Private methods
  // -----------------------------------
  /**
   * Returns the minimum criterion value of c3t3 cells contained in \c cells.
   */
  template <typename SliverCriterion>
  FT min_sliver_in_c3t3_value(const Cell_vector& cells,
                              const SliverCriterion& criterion,
                              const bool use_cache = true) const
  {
    // Get complex cells only
    Cell_vector c3t3_cells;
    std::remove_copy_if(cells.begin(),
                        cells.end(),
                        std::back_inserter(c3t3_cells),
                        std::not1(Is_in_c3t3<Cell_handle>(c3t3_)) );
    
    return min_sliver_value(c3t3_cells,criterion,use_cache);
  }
  
  /**
   * Removes objects of [begin,end[ range from \c c3t3_
   */
  template<typename ForwardIterator>
  void remove_from_c3t3(ForwardIterator begin, ForwardIterator end)
  {
    while ( begin != end )
      c3t3_.remove_from_complex(*begin++);
  }
  
  /**
   * Remove cells and facets of \c cells from c3t3
   */
  template < typename ForwardIterator >
  void remove_cells_and_facets_from_c3t3(ForwardIterator cells_begin,
                                         ForwardIterator cells_end)
  {
    Facet_vector facets = get_facets(cells_begin,cells_end);
    remove_from_c3t3(facets.begin(), facets.end());
    remove_from_c3t3(cells_begin, cells_end);    
  }
  
  /**
   * Insert into \c out the vertices of range [cells_begin,cells_end[
   */
  template <typename InputIterator, typename OutputIterator>
  void fill_modified_vertices(InputIterator cells_begin,
                              InputIterator cells_end,
                              const Vertex_handle& vertex,
                              OutputIterator out) const;
  
  
  /**
   * Update mesh iff sliver criterion value does not decrease.
   */
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh_no_topo_change(const Point_3& new_point,
                             const Vertex_handle& old_vertex,
                             const SliverCriterion& criterion,
                             OutputIterator modified_vertices);
  
  template <typename SliverCriterion, typename OutputIterator>
  std::pair<bool,Vertex_handle>
  update_mesh_topo_change(const Point_3& new_point,
                          const Vertex_handle& old_vertex,
                          const SliverCriterion& criterion,
                          OutputIterator modified_vertices);
  
  /**
   * Move point and returns the set of cells that are not valid anymore, and
   * the set of cells which have been deleted by the move process.
   */
  template < typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle move_point(const Vertex_handle& old_vertex,
                           const Point_3& new_location,
                           OutdatedCellsOutputIterator outdated_cells,
                           DeletedCellsOutputIterator deleted_cells);
  
  template < typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle
  move_point_topo_change(const Vertex_handle& old_vertex,
                         const Point_3& new_location,
                         OutdatedCellsOutputIterator outdated_cells,
                         DeletedCellsOutputIterator deleted_cells);
  
  template < typename ConflictCellsInputIterator,
             typename OutdatedCellsOutputIterator,
             typename DeletedCellsOutputIterator >
  Vertex_handle 
  move_point_topo_change_conflict_zone_known(
     const Vertex_handle& old_vertex,
     const Point_3& new_location,
     ConflictCellsInputIterator conflict_cells_begin,
     ConflictCellsInputIterator conflict_cells_end,
     OutdatedCellsOutputIterator outdated_cells,
     DeletedCellsOutputIterator deleted_cells);

  Vertex_handle move_point_topo_change(const Vertex_handle& old_vertex,
                                       const Point_3& new_location);
  
  template < typename OutdatedCellsOutputIterator >
  Vertex_handle
  move_point_no_topo_change(const Vertex_handle& old_vertex,
                            const Point_3& new_location,
                            OutdatedCellsOutputIterator outdated_cells);

  Vertex_handle
  move_point_no_topo_change(const Vertex_handle& old_vertex,
                            const Point_3& new_location);
  
  /**
   * Returns the least square plane from v, using adjacent surface points
   */
  Plane_3 get_least_square_surface_plane(const Vertex_handle& v,
                                         Point_3& ref_point) const;
  
  /**
   * @brief Returns the projection of \c p, using direction of 
   * \c projection_vector
   */
  Point_3
  project_on_surface_aux(const Point_3& p,
                         const Point_3& ref_point,
                         const Vector_3& projection_vector) const;
  
  /**
   * Reverts the move from \c old_point to \c new_vertex. Returns the inserted
   * vertex located at \c old_point.
   */
  Vertex_handle revert_move(const Vertex_handle& new_vertex,
                            const Point_3& old_point)
  {
    // Store vertex location
    Point_3 new_point = new_vertex->point();
    
    // Move vertex
    Vertex_handle revert_vertex = move_point_topo_change(new_vertex, old_point);
    CGAL_assertion(Vertex_handle() != revert_vertex);
    
    // Restore cell & facet attributes
    // TODO: optimize this to use knowledge on deleted elements ?
    Cell_vector conflict_cells;
    conflict_cells.reserve(64);
    get_conflict_zone_topo_change(revert_vertex, new_point,
                                  std::back_inserter(conflict_cells));
    
    restore_mesh(conflict_cells.begin(), conflict_cells.end());
    
    return revert_vertex;
  }
  
  /**
   * Returns the boundary of facets of \c facets
   */
  Facet_boundary get_surface_boundary(const Facet_vector& facets) const;
  
  /**
   * Returns the boundary of facets of \c cells
   */
  Facet_boundary get_surface_boundary(const Cell_vector& cells) const
  {
    return get_surface_boundary(get_facets(cells));
  }
  
  /**
   * Returns false if there is a vertex belonging to one facet of \c facets 
   * which has not his dimension < 3
   */
  bool check_no_inside_vertices(const Facet_vector& facets) const;
  
  /**
   * Returns the impacted cells when moving \c vertex to \c conflict_point
   */
  template <typename OutputIterator>
  OutputIterator
  get_conflict_zone_no_topo_change(const Vertex_handle& vertex,
                                   OutputIterator conflict_cells) const;
  
  template <typename OutputIterator>
  OutputIterator
  get_conflict_zone_topo_change(const Vertex_handle& vertex,
                                const Point_3& conflict_point,
                                OutputIterator conflict_cells) const;
  
  /**
   * Updates \c boundary wrt \c edge: if edge is already in boundary we remove
   * it, else we add it.
   */
  void update_boundary(Facet_boundary& boundary,
                       const Ordered_edge& edge,
                       const Surface_patch_index& surface_index) const
  {
    typename Facet_boundary::iterator boundary_it =
      boundary.find(std::make_pair(edge,surface_index));
    
    if ( boundary_it != boundary.end() )
      boundary.erase(boundary_it);
    else
      boundary.insert(std::make_pair(edge,surface_index));
  }
  
  /**
   * Returns the facets of \c cells (returns each facet only once i.e. use
   * canonical facet)
   */
  Facet_vector get_facets(const Cell_vector& cells) const
  {
    return get_facets(cells.begin(),cells.end());
  }
  
  /**
   * Returns the facets of \c cells (returns each facet only once i.e. use
   * canonical facet)
   */
  template <typename ForwardIterator>
  Facet_vector get_facets(ForwardIterator first_cell,
                          ForwardIterator last_cell) const
  {
    // Get all facets
    typedef Get_all_facets<std::back_insert_iterator<Facet_vector> > Get_facets;
    
    Facet_vector all_facets;
    all_facets.reserve(64);
    std::for_each(first_cell,
                  last_cell,
                  Get_facets(tr_,std::back_inserter(all_facets)));
    
    std::sort(all_facets.begin(), all_facets.end());
    
    // Keep one copy of each facet (maybe copy could be avoided)
    //    typename Facet_vector::iterator all_facets_end =
    //      std::unique(all_facets.begin(), all_facets.end());
    Facet_vector facets;
    facets.reserve(64);
    std::unique_copy(all_facets.begin(),
                     all_facets.end(),
                     std::back_inserter(facets));
    
    return facets;
  }
  
  /**
   * Returns true if all surface facets of cells are really in restricted
   * Delaunay.
   */
  bool verify_surface(const Cell_vector& cells) const
  {
    // Naive implementation.
    // Todo: improve this (maybe we don't have to check if no facet is on surface)
    Facet_vector facets = get_facets(cells);
    Facet_vector surface_facets;
    
    // Check that nothing changed
    Update_c3t3 checker(domain_,c3t3_);
    for ( typename Facet_vector::iterator fit = facets.begin() ;
          fit != facets.end() ;
          ++fit )
    {
      if ( c3t3_.is_in_complex(*fit) )
      {
        surface_facets.push_back(*fit);
      }
      
      if ( c3t3_.is_in_complex(*fit) != checker(*fit,false) )
        return false;
    }
    
    // Facet surface center must be updated if verify_surface is ok
    std::for_each(surface_facets.begin(),surface_facets.end(),checker);
    
    return true;
  }
  
  
  /**
   * Restore mesh for cells and facets of \c cells, using domain_
   */ 
  void restore_mesh(const Cell_vector& cells)
  {
    restore_mesh(cells.begin(), cells.end());
  }
  
  /**
   * Restore mesh for cells and facets of \c cells, using domain_
   */ 
  template <typename ForwardIterator>
  void restore_mesh(ForwardIterator first_cell, ForwardIterator last_cell)
  {
    Facet_vector facets = get_facets(first_cell, last_cell);
    restore_mesh(first_cell, last_cell, facets.begin(), facets.end());
  }
  
  /**
   * Restore mesh for cells of \c cells and facets of \c facets, using domain_
   */
  template <typename CellForwardIterator, typename FacetForwardIterator>
  void restore_mesh(CellForwardIterator first_cell,
                    CellForwardIterator last_cell,
                    FacetForwardIterator first_facet,
                    FacetForwardIterator last_facet)
  {
    // Update mesh
    Update_c3t3 updater(domain_,c3t3_);
    std::for_each(first_facet, last_facet, updater);
    std::for_each(first_cell, last_cell, updater);
  }
  
  /**
   * Returns true if facets of \c facets have the same boundary as 
   * \c old_boundary
   */
  bool check_surface_mesh(const Facet_vector& facets,
                          const Facet_boundary& old_boundary) const
  {
    Facet_boundary new_boundary = get_surface_boundary(facets);

    return ( old_boundary.size() == new_boundary.size()
            && std::equal(new_boundary.begin(),
                          new_boundary.end(),
                          old_boundary.begin()) );
  }
  
  /**
   * Restore mesh for cells and facets of \c cells, then check that the new
   * boundary of facets of \c cells is the same as \c old_boundary.
   */
  bool restore_and_check_mesh(const Cell_vector& cells,
                              const Facet_boundary& old_boundary)
  {
    Facet_vector facets = get_facets(cells);
    restore_mesh(cells.begin(), cells.end(), facets.begin(), facets.end());
    return check_mesh(facets, old_boundary);
  }
  
  void set_facet_visited(const Facet& facet)
  {
    facet.first->set_facet_visited(facet.second);
    const Facet mirror_facet = tr_.mirror_facet(facet);
    mirror_facet.first->set_facet_visited(mirror_facet.second);    
  }
  
  /**
   * Orders handles \c h1, \c h2 & \c h3
   */ 
  template <typename Handle>
  void order_handles(Handle& h1, Handle& h2, Handle& h3) const
  {
    if ( h2 < h1 )
      std::swap(h1,h2);
    
    if ( h3 < h2 )
    {
      std::swap(h2,h3);
      
      if ( h2 < h1 ) // don't need to compare h2 & h1 if h2 didn't change
        std::swap(h1,h2);
    }
  }

  template < typename ForwardIterator >
  void reset_cache_validity(ForwardIterator cells_begin,
                            ForwardIterator cells_end) const
  {
    namespace bl = boost::lambda;
    std::for_each(cells_begin, cells_end,
                  bl::bind(&Cell::reset_cache_validity, *bl::_1) );
  }
  
  
private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const MeshDomain& domain_;
};
  
  
template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh(const Point_3& new_location,
            const Vertex_handle& old_vertex,
            const SliverCriterion& criterion,
            OutputIterator modified_vertices)
{
  if ( Th().no_topological_change(tr_, old_vertex, new_location) )
  {
    return update_mesh_no_topo_change(new_location,
                                      old_vertex,
                                      criterion,
                                      modified_vertices);
  }
  else
  {
    return update_mesh_topo_change(new_location,
                                   old_vertex,
                                   criterion,
                                   modified_vertices);
  }
}

template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh_no_topo_change(const Point_3& new_location,
                           const Vertex_handle& vertex,
                           const SliverCriterion& criterion,
                           OutputIterator modified_vertices)
{
  // Get conflict zone
  Cell_vector conflict_cells;
  conflict_cells.reserve(64);
  get_conflict_zone_no_topo_change(vertex, std::back_inserter(conflict_cells));
  
  // Get old values
  FT old_sliver_value = min_sliver_in_c3t3_value(conflict_cells, criterion);
  Point_3 old_location = vertex->point();
  
  // Move point
  move_point_no_topo_change(vertex,new_location);
    
  // Get new criterion value (conflict_zone did not change)
  const FT new_sliver_value = 
    min_sliver_in_c3t3_value(conflict_cells, criterion, false);
  
  // Check that mesh is still valid
  if ( new_sliver_value > old_sliver_value && verify_surface(conflict_cells) )
  {
    fill_modified_vertices(conflict_cells.begin(), conflict_cells.end(),
                           vertex, modified_vertices);
    return std::make_pair(true,vertex);
  }
  else
  {
    // revert move
    move_point_no_topo_change(vertex,old_location);
    reset_cache_validity(conflict_cells.begin(), conflict_cells.end());
    return std::make_pair(false,vertex);
  }
}


template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
std::pair<bool,typename C3T3_helpers<C3T3,MD>::Vertex_handle>
C3T3_helpers<C3T3,MD>::  
update_mesh_topo_change(const Point_3& new_location,
                        const Vertex_handle& old_vertex,
                        const SliverCriterion& criterion,
                        OutputIterator modified_vertices)
{
  Cell_vector conflict_cells;
  conflict_cells.reserve(64);
  get_conflict_zone_topo_change(old_vertex, new_location,
                                std::back_inserter(conflict_cells));
  
  // May happen in case of new_point is located on a vertex
  if ( conflict_cells.empty() )
    return std::make_pair(false,old_vertex);
  
  FT old_sliver_value = min_sliver_in_c3t3_value(conflict_cells, criterion);
  Point_3 old_location = old_vertex->point();
  
  // Keep old boundary
  Facet_boundary old_surface_boundary = get_surface_boundary(conflict_cells);
  
  Cell_vector outdated_cells;
  Vertex_handle new_vertex = 
    move_point_topo_change_conflict_zone_known(old_vertex,
                                               new_location,
                                               conflict_cells.begin(),
                                               conflict_cells.end(),
                                               std::back_inserter(outdated_cells),
                                               CGAL::Emptyset_iterator());
  
  // If nothing changed, return
  if ( old_location == new_vertex->point() ) 
  {
    return std::make_pair(false,old_vertex);
  }
  
  restore_mesh(outdated_cells.begin(),outdated_cells.end());
  FT new_sliver_value = min_sliver_in_c3t3_value(outdated_cells, criterion);
  
  // Check that surface boundary do not change.
  // This check ensures that vertices which are inside c3t3 stay inside. 
  if ( new_sliver_value > old_sliver_value
      && check_surface_mesh(get_facets(outdated_cells), old_surface_boundary) )
  {
    fill_modified_vertices(outdated_cells.begin(), outdated_cells.end(),
                           new_vertex, modified_vertices);
    return std::make_pair(true,new_vertex);
  }
  else
  {
    // Remove from c3t3 cells which will be destroyed by revert_move
    remove_cells_and_facets_from_c3t3(outdated_cells.begin(),
                                      outdated_cells.end());
    
    // Revert move
    Vertex_handle revert_vertex = revert_move(new_vertex, old_location);
    return std::make_pair(false,revert_vertex);
  }
}
  
template <typename C3T3, typename MD>
template <typename OutputIterator>  
typename C3T3_helpers<C3T3,MD>::Vertex_handle
C3T3_helpers<C3T3,MD>::  
update_mesh(const Point_3& new_point,
            const Vertex_handle& old_vertex,
            OutputIterator modified_vertices,
            bool fill_vertices)
{
  Cell_vector outdated_cells;
  Vertex_handle new_vertex = move_point(old_vertex,
                                        new_point,
                                        std::back_inserter(outdated_cells),
                                        CGAL::Emptyset_iterator());
  
  restore_mesh(outdated_cells.begin(),outdated_cells.end());
  
  // Fill modified vertices
  if ( fill_vertices )
  {
    fill_modified_vertices(outdated_cells.begin(), outdated_cells.end(),
                           new_vertex, modified_vertices);        
  }
  
  return new_vertex;  
}
  
  
template <typename C3T3, typename MD>
template <typename ForwardIterator>
void
C3T3_helpers<C3T3,MD>:: 
rebuild_restricted_delaunay(ForwardIterator first_cell,
                            ForwardIterator last_cell,
                            Vertex_set& moving_vertices)
{
  Update_c3t3 updater(domain_,c3t3_);
  
  // Get facets (returns each canonical facet only once)
  Facet_vector facets = get_facets(first_cell, last_cell);
  
  // Updates cells
  while ( first_cell != last_cell )
  {
    const Cell_handle& cell = *first_cell++;
    c3t3_.remove_from_complex(cell);
    updater(cell);
    
    // Update moving_vertices
    if ( c3t3_.is_in_complex(cell) )
    {
      for ( int i=0 ; i<4 ; ++i )
      {
        moving_vertices.insert(cell->vertex(i)); 
      }
    }
  }
  
  // Updates facets
  std::set<Vertex_handle> vertex_to_proj;
  for ( typename Facet_vector::iterator fit = facets.begin() ;
       fit != facets.end() ;
       ++fit )
  {
    // Update facet
    c3t3_.remove_from_complex(*fit);
    updater(*fit);

    // Update vertex_to_proj
    if ( c3t3_.is_in_complex(*fit) )
    {
      // Iterate on vertices
      int k = fit->second;
      for ( int i=1 ; i<4 ; ++i )
      {
        const Vertex_handle& v = fit->first->vertex((k+i)&3);
        if ( c3t3_.in_dimension(v) > 2 )
        { 
          vertex_to_proj.insert(v);
        }
      }
    }
  }
  
  // Project interior vertices
  // TODO : iterate to be sure no interior vertice become on the surface
  // because of move ?
  for ( typename std::set<Vertex_handle>::iterator it = vertex_to_proj.begin() ;
       it != vertex_to_proj.end() ;
       ++it )
  {
    Point_3 new_pos = project_on_surface((*it)->point(),*it);

    if ( new_pos != Point_3() )
    {
      Vertex_handle new_vertex = update_mesh(new_pos,*it);
      c3t3_.set_dimension(new_vertex,2);
      
      // Update moving vertices (it has become new_vertex)
      moving_vertices.erase(*it);
      moving_vertices.insert(new_vertex);
    }
  }
}

  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point(const Vertex_handle& old_vertex,
           const Point_3& new_location,
           Cell_set& outdated_cells_set)
{
  Cell_vector outdated_cells;
  Cell_vector deleted_cells;
  
  Vertex_handle new_vertex =
    move_point(old_vertex,
               new_location,
               std::back_inserter(outdated_cells),
               std::back_inserter(deleted_cells));

  // Get Cell_set::erase and Cell_set::insert function pointers
  namespace bl = boost::lambda;
  typename Cell_set::size_type (Cell_set::*erase_cell)
    (const typename Cell_set::key_type&) = &Cell_set::erase;
  
  std::pair<typename Cell_set::iterator,bool> (Cell_set::*insert_cell)
    (const typename Cell_set::value_type&) = &Cell_set::insert;

  // Update outdated_cells_set
  std::for_each(deleted_cells.begin(),deleted_cells.end(),
                bl::bind(erase_cell, &outdated_cells_set, bl::_1) );
  
  std::for_each(outdated_cells.begin(),outdated_cells.end(),
                bl::bind(insert_cell, &outdated_cells_set, bl::_1) );
  
  return new_vertex;
}  


template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator,
          typename DeletedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point(const Vertex_handle& old_vertex,
           const Point_3& new_location,
           OutdatedCellsOutputIterator outdated_cells,
           DeletedCellsOutputIterator deleted_cells)
{
  if ( Th().no_topological_change(tr_, old_vertex, new_location) )
  {
    return move_point_no_topo_change(old_vertex,
                                     new_location,
                                     outdated_cells);
  }
  else
  {
    return move_point_topo_change(old_vertex,
                                  new_location,
                                  outdated_cells,
                                  deleted_cells);
  }
}
  
  
template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator,
          typename DeletedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change(const Vertex_handle& old_vertex,
                       const Point_3& new_location,
                       OutdatedCellsOutputIterator outdated_cells,
                       DeletedCellsOutputIterator deleted_cells)
{
  Cell_vector conflict_cells;
  conflict_cells.reserve(64);
  get_conflict_zone_topo_change(old_vertex, new_location,
                                std::back_inserter(conflict_cells));
  
  return move_point_topo_change_conflict_zone_known(old_vertex,
                                                    new_location,
                                                    conflict_cells.begin(),
                                                    conflict_cells.end(),
                                                    outdated_cells,
                                                    deleted_cells);
}
  

template <typename C3T3, typename MD>
template < typename ConflictCellsInputIterator,
           typename OutdatedCellsOutputIterator,
           typename DeletedCellsOutputIterator >
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_topo_change_conflict_zone_known(
    const Vertex_handle& old_vertex,
    const Point_3& new_location,
    ConflictCellsInputIterator conflict_cells_begin,
    ConflictCellsInputIterator conflict_cells_end,
    OutdatedCellsOutputIterator outdated_cells,
    DeletedCellsOutputIterator deleted_cells)
{
  Point_3 old_location = old_vertex->point();
  
  // Remove conflict zone cells from c3t3 (cells will be destroyed)  
  remove_cells_and_facets_from_c3t3(conflict_cells_begin, conflict_cells_end);
  
  // Move point
  Vertex_handle new_vertex = move_point_topo_change(old_vertex,new_location);
  
  // If nothing changed, return
  if ( Vertex_handle() == new_vertex )
  {
    std::copy(conflict_cells_begin,conflict_cells_end,outdated_cells);
    return old_vertex;
  }
  
  // Get conflict zone in new triangulation and set cells outdated
  Cell_vector new_conflict_cells;
  new_conflict_cells.reserve(64);
  get_conflict_zone_topo_change(new_vertex, old_location,
                                std::back_inserter(new_conflict_cells));
  
  std::copy(conflict_cells_begin,conflict_cells_end,deleted_cells);
  std::copy(new_conflict_cells.begin(),new_conflict_cells.end(),outdated_cells);

  return new_vertex;
}


template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>::
move_point_topo_change(const Vertex_handle& old_vertex,
                       const Point_3& new_location)
{
  // Insert new_vertex, remove old_vertex
  int dimension = c3t3_.in_dimension(old_vertex);
  Index vertice_index = c3t3_.index(old_vertex);
  FT meshing_info = old_vertex->meshing_info();
  
  // insert new point
  Vertex_handle new_vertex = tr_.insert(new_location,old_vertex->cell());
  // If new_point is hidden, return default constructed handle
  if ( Vertex_handle() == new_vertex ) { return Vertex_handle(); }
  // remove old point
  tr_.remove(old_vertex);
  
  c3t3_.set_dimension(new_vertex,dimension);
  c3t3_.set_index(new_vertex,vertice_index);
  new_vertex->set_meshing_info(meshing_info);
  
  return new_vertex;
}
  

template <typename C3T3, typename MD>
template <typename OutdatedCellsOutputIterator>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_no_topo_change(const Vertex_handle& old_vertex,
                          const Point_3& new_location,
                          OutdatedCellsOutputIterator outdated_cells)
{
  get_conflict_zone_no_topo_change(old_vertex, outdated_cells);
  return move_point_no_topo_change(old_vertex, new_location);
}


template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Vertex_handle 
C3T3_helpers<C3T3,MD>:: 
move_point_no_topo_change(const Vertex_handle& old_vertex,
                          const Point_3& new_location)
{  
  // Change vertex location
  old_vertex->set_point(new_location);
  return old_vertex;  
}


/**
 * @brief Returns the projection of \c p, using direction of 
 * \c projection_vector
 */
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Point_3
C3T3_helpers<C3T3,MD>:: 
project_on_surface_aux(const Point_3& p,
                       const Point_3& ref_point,
                       const Vector_3& projection_vector) const
{
  typedef typename Gt::Segment_3 Segment_3;
  
  // Build a segment directed as projection_direction,
  typename Gt::Compute_squared_distance_3 sq_distance =
    Gt().compute_squared_distance_3_object();
  
  typename Gt::Compute_squared_length_3 sq_length =
    Gt().compute_squared_length_3_object();
  
  typename Gt::Construct_scaled_vector_3 scale =
    Gt().construct_scaled_vector_3_object();
    
  typename Gt::Is_degenerate_3 is_degenerate =
    Gt().is_degenerate_3_object();
  
  typename MD::Do_intersect_surface do_intersect =
    domain_.do_intersect_surface_object();
  
  typename MD::Construct_intersection intersection =
    domain_.construct_intersection_object();
  
  const FT sq_dist = sq_distance(p,ref_point);
  const FT sq_proj_length = sq_length(projection_vector);
  
  if ( CGAL_NTS is_zero(sq_proj_length) )
    return ref_point;
  
  const Vector_3 projection_scaled_vector =
    scale(projection_vector, CGAL::sqrt(sq_dist/sq_proj_length));
  
  const Point_3 source = p + projection_scaled_vector;
  const Point_3 target = p - projection_scaled_vector;
  
  const Segment_3 proj_segment(source,target);
  
  if ( is_degenerate(proj_segment) )
    return ref_point;
  
  if ( do_intersect(proj_segment) )
    return CGAL::cpp11::get<0>(intersection(proj_segment));
  else
    return ref_point;  
}

  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Plane_3
C3T3_helpers<C3T3,MD>:: 
get_least_square_surface_plane(const Vertex_handle& v,
                               Point_3& reference_point) const
{
  // Get incident facets
  Facet_vector facets;
  tr_.finite_incident_facets(v,std::back_inserter(facets));

  // Get adjacent surface points
  std::vector<Point_3> surface_point_vector;
  for ( typename Facet_vector::iterator fit = facets.begin() ;
       fit != facets.end() ;
       ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      const Cell_handle& cell = fit->first;
      const int& i = fit->second;
      
      surface_point_vector.push_back(cell->get_facet_surface_center(i));
    }
  }

  // In some cases point is not a real surface point
  if ( surface_point_vector.empty() )
    return Plane_3();
    
  // Compute least square fitting plane
  Plane_3 plane;
  CGAL::linear_least_squares_fitting_3(surface_point_vector.begin(),
                                       surface_point_vector.end(),
                                       plane,
                                       Dimension_tag<0>());
  
  reference_point = surface_point_vector.front();

  return plane;
}
  
  
  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Point_3
C3T3_helpers<C3T3,MD>:: 
project_on_surface(const Point_3& p,
                   const Vertex_handle& v) const
{
  // Get plane
  Point_3 reference_point(CGAL::ORIGIN);
  Plane_3 plane = get_least_square_surface_plane(v,reference_point);
  
  if ( reference_point == CGAL::ORIGIN )
    return p;
  
  // Project
  if ( p != v->point() )
    return project_on_surface_aux(p,
                                  v->point(),
                                  plane.orthogonal_vector());
  else
    return project_on_surface_aux(p,
                                  reference_point,
                                  plane.orthogonal_vector());
}

  
  
template <typename C3T3, typename MD>
template <typename SliverCriterion> 
typename C3T3_helpers<C3T3,MD>::FT
C3T3_helpers<C3T3,MD>::
min_incident_value(const Vertex_handle& vh,
                   const SliverCriterion& criterion) const
{
  Cell_vector incident_cells;
  tr_.finite_incident_cells(vh,std::back_inserter(incident_cells));
  
  return min_sliver_in_c3t3_value(incident_cells, criterion);
}

  
template <typename C3T3, typename MD>
template <typename SliverCriterion, typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
incident_slivers(const Vertex_handle& v,
                 const SliverCriterion& criterion,
                 const FT& sliver_bound,
                 OutputIterator out) const
{
  typedef SliverCriterion Sc;
  
  std::vector<Cell_handle> incident_cells;
  tr_.incident_cells(v, std::back_inserter(incident_cells));
  
  std::remove_copy_if(incident_cells.begin(),
                      incident_cells.end(),
                      out,
                      std::not1(Is_sliver<Sc>(c3t3_,criterion,sliver_bound)));
  
  return out;
}

  
  
template <typename C3T3, typename MD>
template <typename SliverCriterion>
typename C3T3_helpers<C3T3,MD>::FT
C3T3_helpers<C3T3,MD>::
min_sliver_value(const Cell_vector& cells,
                 const SliverCriterion& criterion,
                 const bool use_cache) const
{
  using boost::make_transform_iterator;
  
  if ( cells.empty() )
    return SliverCriterion::max_value;
  
  if ( ! use_cache )
  { 
    reset_cache_validity(cells.begin(),cells.end());
  }
  
  // Return min dihedral angle
  Sliver_criterion_value<SliverCriterion> sc_value(tr_,criterion);
  
  return *(std::min_element(make_transform_iterator(cells.begin(),sc_value),
                            make_transform_iterator(cells.end(),sc_value)));
}
  
  
template <typename C3T3, typename MD>
template <typename InputIterator, typename OutputIterator>
void
C3T3_helpers<C3T3,MD>::
fill_modified_vertices(InputIterator cells_begin,
                       InputIterator cells_end,
                       const Vertex_handle& vertex,
                       OutputIterator out) const
{
  std::set<Vertex_handle> already_inserted_vertices;
  // Dont insert vertex in out
  already_inserted_vertices.insert(vertex);
  
  for ( InputIterator it = cells_begin ; it != cells_end ; ++it )
  {
    for ( int i=0 ; i<4 ; ++i )
    {
      // Insert vertices if not already inserted
      const Vertex_handle& current_vertex = (*it)->vertex(i);
      if ( !tr_.is_infinite(current_vertex)
          && already_inserted_vertices.find(current_vertex) ==
             already_inserted_vertices.end() )
      {
        *out++ = current_vertex;
        already_inserted_vertices.insert(current_vertex);
      }
    }
  }
}
  
  
template <typename C3T3, typename MD>
template <typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
get_conflict_zone_no_topo_change(const Vertex_handle& vertex,
                                 OutputIterator conflict_cells) const
{  
  return tr_.incident_cells(vertex,conflict_cells);
}
  
template <typename C3T3, typename MD>
template <typename OutputIterator>
OutputIterator
C3T3_helpers<C3T3,MD>::
get_conflict_zone_topo_change(const Vertex_handle& vertex,
                              const Point_3& conflict_point,
                              OutputIterator conflict_cells) const
{
  // Get triangulation_vertex incident cells
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(vertex, std::back_inserter(incident_cells));
  
  // Get conflict_point conflict zone
  Cell_vector deleted_cells;
  deleted_cells.reserve(64);
  
  // Vertex removal is forbidden 
  int li=0;
  int lj=0;
  typename Tr::Locate_type locate_type;
  Cell_handle cell = tr_.locate(conflict_point,
                                locate_type,
                                li,
                                lj,
                                vertex->cell());
  
  if ( Tr::VERTEX == locate_type )
    return conflict_cells;
  
  // Find conflict zone
  tr_.find_conflicts(conflict_point,
                     cell,
                     CGAL::Emptyset_iterator(),
                     std::back_inserter(deleted_cells),
                     CGAL::Emptyset_iterator());
  
  // Compute union of conflict_point conflict zone and triangulation_vertex
  // incident cells
  std::sort(deleted_cells.begin(),deleted_cells.end());
  std::sort(incident_cells.begin(),incident_cells.end());
  
  Cell_vector conflict_zone_cells;
  std::set_union(deleted_cells.begin(), deleted_cells.end(),
                 incident_cells.begin(), incident_cells.end(),
                 conflict_cells);
  
  return conflict_cells;
}
  
  
template <typename C3T3, typename MD>
typename C3T3_helpers<C3T3,MD>::Facet_boundary
C3T3_helpers<C3T3,MD>::
get_surface_boundary(const Facet_vector& facets) const
{
  Facet_boundary boundary;
  typename Facet_vector::const_iterator fit = facets.begin();
  for ( ; fit != facets.end() ; ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      const Surface_patch_index surface_index = c3t3_.surface_patch_index(*fit);
      const int k = fit->second;
      Vertex_handle v1 = fit->first->vertex((k+1)&3);
      Vertex_handle v2 = fit->first->vertex((k+2)&3);
      Vertex_handle v3 = fit->first->vertex((k+3)&3);
      
      // Check that each vertex is a surface one
      // This is a trick to ensure that in_domain vertices stay inside
      if ( c3t3_.in_dimension(v1) > 2
          || c3t3_.in_dimension(v2) > 2
          || c3t3_.in_dimension(v3) > 2 )
      {
        // Ordered_edge(tr_.infinite_vertex(),v1) can't be in boundary
        // So if there is a boundary facets which is not built on 3 boundary
        // vertices, check of boundary equality before and after the move will
        // fail (we know that this is not the case before)
        update_boundary(boundary,
                        Ordered_edge(Vertex_handle(),Vertex_handle()),
                        Surface_patch_index());
        return boundary;
      }
      
      order_handles(v1,v2,v3);
      
      CGAL_assertion(v1<v2);
      CGAL_assertion(v2<v3);
      
      update_boundary(boundary, Ordered_edge(v1,v2), surface_index);
      update_boundary(boundary, Ordered_edge(v1,v3), surface_index);
      update_boundary(boundary, Ordered_edge(v2,v3), surface_index);
    }
  }

  return boundary;
}
  
template <typename C3T3, typename MD>
bool
C3T3_helpers<C3T3,MD>::
check_no_inside_vertices(const Facet_vector& facets) const
{
  typename Facet_vector::const_iterator fit = facets.begin();
  for ( ; fit != facets.end() ; ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      const int k = fit->second;
      const Vertex_handle& v1 = fit->first->vertex((k+1)&3);
      const Vertex_handle& v2 = fit->first->vertex((k+2)&3);
      const Vertex_handle& v3 = fit->first->vertex((k+3)&3);
      
      // Check that each vertex is a surface one
      if ( c3t3_.in_dimension(v1) > 2
          || c3t3_.in_dimension(v2) > 2
          || c3t3_.in_dimension(v3) > 2 )
      {
        return false;
      }
    }
  }
  
  return true;
}
  
  
} // end namespace Mesh_3
} // end namespace CGAL

#endif // CGAL_MESH_3_C3T3_HELPERS_H
