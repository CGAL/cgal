// Copyright (c) 2004-2009  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_MESH_3_REFINE_CELLS_3_H
#define CGAL_MESH_3_REFINE_CELLS_3_H


#include <CGAL/Mesher_level.h>
#include <CGAL/Mesher_level_default_implementations.h>
#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>
#include <CGAL/Meshes/Double_map_container.h>

#include <boost/format.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <sstream>


namespace CGAL {
  
namespace Mesh_3 {

// Helper meta-programming functions, to allow backward compatibility.
//
//   - Has_Is_cell_bad and Had_Cell_badness are used to detect if a model
//     of the MeshCellCriteria_3 concept follows the specifications of
//     CGAL-3.7 (with Cell_badness) or later (with Is_cell_bad).
//
//   - Then the meta-function Get_Is_cell_bad is used to get the actual
//     type, either Cell_criteria::Cell_badness or
//      Cell_criteria::Is_cell_bad.

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_Is_cell_bad, Is_cell_bad, true)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_Cell_badness, Cell_badness, false)

// template class, used when use_cell_badness = false
template <typename Cell_criteria, 
          bool use_cell_badness = (!Has_Is_cell_bad<Cell_criteria>::value) &&
                                    Has_Cell_badness<Cell_criteria>::value >
struct Get_Is_cell_bad {
  typedef typename Cell_criteria::Is_cell_bad Type;
  typedef Type type; // compatibility with Boost
};

// partial specialization when use_cell_badness == true
template <typename Cell_criteria>
struct Get_Is_cell_bad<Cell_criteria, true> {
  typedef typename Cell_criteria::Cell_badness Type;
  typedef Type type;
};



// Class Refine_cells_3
//
// Template parameters should be models of
// Tr         : MeshTriangulation_3
// Criteria   : MeshCellsCriteria_3
// MeshDomain : MeshTraits_3
//
// Implements a Mesher_level for cells
template<class Tr,
         class Criteria,
         class MeshDomain,
         class Complex3InTriangulation3,
         class Previous_,
         class Container_ = Meshes::Double_map_container<
                                          typename Tr::Cell_handle,
                                          typename Criteria::Cell_quality> >
class Refine_cells_3
  : public Mesher_level<Tr,
                        Refine_cells_3<Tr,
                                       Criteria,
                                       MeshDomain,
                                       Complex3InTriangulation3,
                                       Previous_,
                                       Container_>,
                        typename Tr::Cell_handle,
                        Previous_,
                        Triangulation_mesher_level_traits_3<Tr> >
  , public Container_
  , public No_test_point_conflict
  , public No_after_no_insertion
  , public No_before_conflicts
{
private:
  // Internal types
  typedef typename Tr::Facet Facet;
  typedef typename MeshDomain::Subdomain_index  Subdomain_index;
  typedef typename MeshDomain::Index  Index;
  typedef typename Get_Is_cell_bad<Criteria>::Type Is_cell_bad;
  
  // Self
  typedef Refine_cells_3<Tr,
    Criteria,
    MeshDomain,
    Complex3InTriangulation3,
    Previous_,
    Container_>       Self;
  
public:  
  typedef typename Tr::Point Point;
  typedef typename Tr::Cell Cell;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Criteria::Cell_quality Cell_quality;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
  typedef Complex3InTriangulation3 C3T3;
  
  
  // Constructor
  Refine_cells_3(Tr& triangulation,
                 const Criteria& criteria,
                 const MeshDomain& oracle,
                 Previous_& previous,
                 C3T3& c3t3) ;
  
  // Destructor
  virtual ~Refine_cells_3() { }
  
  // Get a reference on triangulation
  Tr& triangulation_ref_impl() { return r_tr_; }
  const Tr& triangulation_ref_impl() const { return r_tr_; }
  
  // Initialization function
  void scan_triangulation_impl();
  
  // Gets the point to insert from the element to refine
  Point refinement_point_impl(const Cell_handle& cell) const
  {
    last_vertex_index_ = r_oracle_.index_from_subdomain_index(
        cell->subdomain_index());
    //    last_vertex_index_ = Index(cell->subdomain_index());
    // NB : dual() is optimized when the cell base class has circumcenter()
    return r_tr_.dual(cell);
  }
  
  // Returns the conflicts zone
  Zone conflicts_zone_impl(const Point& point, const Cell_handle& cell) const;
  
  // Job to do before insertion
  void before_insertion_impl(const Cell_handle&, const Point&, Zone& zone)
  {
    before_insertion_handle_cells_in_conflict_zone(zone);
  }
  
  // Job to do after insertion
  void after_insertion_impl(const Vertex_handle& v)
#ifndef CGAL_MESH_3_USE_OLD_SURFACE_RESTRICTED_DELAUNAY_UPDATE
  { update_star_self(v); }
#else
  { update_star(v); }
#endif

  // Insertion implementation ; returns the inserted vertex
  Vertex_handle insert_impl(const Point& p, const Zone& zone);
  
  // Updates cells incident to vertex, and add them to queue if needed
  void update_star(const Vertex_handle& vertex);
  
  /// Handle cells contained in \c zone (before their destruction by insertion)
  void before_insertion_handle_cells_in_conflict_zone(Zone& zone);
  
  std::string debug_info() const
  {
    std::stringstream s;
    s << this->previous().debug_info() << "," << this->size();
    return s.str();
  }
  
  std::string debug_info_header() const
  {
    std::stringstream s;
    s << this->previous().debug_info_header() <<  "," << "#tets to refine";
    return s.str();
  }
  
#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  std::size_t queue_size() const { return this->size(); }
#endif
  
private:
  /// Adds \c cell to the refinement queue if needed
  void treat_new_cell(const Cell_handle& cell);
  
  /// Computes badness and add to queue if needed
  void compute_badness(const Cell_handle& cell);
  
  // Updates cells incident to vertex, and add them to queue if needed
  void update_star_self(const Vertex_handle& vertex);
  
  /// Set \c cell to domain, with subdomain index \c index
  void set_cell_in_domain(const Cell_handle& cell,
                          const Subdomain_index& index)
  {
    r_c3t3_.add_to_complex(cell, index);
  }
  
  /// Removes \c cell from domain
  void remove_cell_from_domain(const Cell_handle& cell)
  {
    r_c3t3_.remove_from_complex(cell);
  }
  
  /// Sets index and dimension of vertex \c v
  void set_vertex_properties(Vertex_handle& v, const Index& index)
  {
    r_c3t3_.set_index(v, index);
    // Set dimension of v: v is inside volume by construction, so dimension=3
    v->set_dimension(3);
  }
  
  /// Get mirror facet
  Facet mirror_facet(const Facet& f) const { return r_tr_.mirror_facet(f); };
  Facet mirror_facet(const Cell_handle& c, const int i) const
  { return mirror_facet(std::make_pair(c,i)); }
  
private:
  /// The triangulation
  Tr& r_tr_;
  /// The cell criteria
  const Criteria& r_criteria_;
  /// The oracle
  const MeshDomain& r_oracle_;
  /// The mesh result
  C3T3& r_c3t3_;
  
  //-------------------------------------------------------
  // Cache objects
  //-------------------------------------------------------
  /// Stores index of vertex that may be inserted into triangulation
  mutable Index last_vertex_index_;
  
private:
  // Disabled copy constructor
  Refine_cells_3(const Self& src);
  // Disabled assignment operator
  Self& operator=(const Self& src);
  
};  // end class Refine_cells_3






template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
Refine_cells_3(Tr& triangulation,
               const Cr& criteria,
               const MD& oracle,
               P_& previous,
               C3T3& c3t3)
  : Mesher_level<Tr, Self, Cell_handle, P_,
      Triangulation_mesher_level_traits_3<Tr> >(previous)
  , C_()
  , No_test_point_conflict()
  , No_after_no_insertion()
  , No_before_conflicts()
  , r_tr_(triangulation)
  , r_criteria_(criteria)
  , r_oracle_(oracle)
  , r_c3t3_(c3t3)
  , last_vertex_index_()
{
}




template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
scan_triangulation_impl()
{
  typedef typename Tr::Finite_cells_iterator Finite_cell_iterator;
  
  for(Finite_cell_iterator cell_it = r_tr_.finite_cells_begin();
      cell_it != r_tr_.finite_cells_end();
      ++cell_it)
  {
    treat_new_cell(cell_it);
  }
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
typename Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::Zone
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
conflicts_zone_impl(const Point& point,
                    const Cell_handle& cell) const
{
  Zone zone;
  zone.cell = cell;
  zone.locate_type = Tr::CELL;
  
  r_tr_.find_conflicts(point,
                       zone.cell,
                       std::back_inserter(zone.boundary_facets),
                       std::back_inserter(zone.cells),
                       std::back_inserter(zone.internal_facets));
  return zone;
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
before_insertion_handle_cells_in_conflict_zone(Zone& zone)
{
  typename Zone::Cells_iterator cit = zone.cells.begin();
  for ( ; cit != zone.cells.end() ; ++cit )
  {
    // Remove cell from refinement queue
    this->remove_element(*cit);
    
    // Remove cell from complex
    remove_cell_from_domain(*cit);
  }
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
update_star(const Vertex_handle& vertex)
{
  typedef std::vector<Cell_handle> Cells;
  typedef typename Cells::iterator Cell_iterator;
  
  // Get the star of v
  Cells incident_cells;
  r_tr_.incident_cells(vertex, std::back_inserter(incident_cells));
  
  // Scan tets of the star of v
  for( Cell_iterator cell_it = incident_cells.begin();
      cell_it != incident_cells.end();
      ++cell_it )
  {
    if( ! r_tr_.is_infinite(*cell_it) )
    {
      // update queue with the new cell if needed
      treat_new_cell(*cell_it);
    }
  }
}

  
template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
update_star_self(const Vertex_handle& vertex)
{
  typedef std::vector<Cell_handle> Cells;
  typedef typename Cells::iterator Cell_iterator;
  
  // Get the star of v
  Cells incident_cells;
  r_tr_.incident_cells(vertex, std::back_inserter(incident_cells));
  
  // Get subdomain index
  Subdomain_index cells_subdomain = r_oracle_.subdomain_index(vertex->index());
  
  // Restore surface & domain
  for( Cell_iterator cell_it = incident_cells.begin();
      cell_it != incident_cells.end();
      ++cell_it )
  {
    CGAL_assertion(!r_tr_.is_infinite(*cell_it));
    
    // Restore surface
    const int& k = (*cell_it)->index(vertex);
    const Facet mirror_f = mirror_facet(*cell_it,k);
    const Cell_handle& neighbor_cell = mirror_f.first;
    const int& neighb_k = mirror_f.second;
    
    if ( neighbor_cell->is_facet_on_surface(neighb_k) )
    {
      // Facet(*cell_it,k) is on surface
      (*cell_it)->set_surface_patch_index(
        k,neighbor_cell->surface_patch_index(neighb_k));

      (*cell_it)->set_facet_surface_center(
        k,neighbor_cell->get_facet_surface_center(neighb_k));

      (*cell_it)->set_facet_surface_center_index(
        k,neighbor_cell->get_facet_surface_center_index(neighb_k));
    }
    
    // Set subdomain index
    set_cell_in_domain(*cell_it, cells_subdomain);
    
    // Add to queue
    compute_badness(*cell_it);
  }
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
treat_new_cell(const Cell_handle& cell)
{
  typedef boost::optional<typename MD::Subdomain_index> Subdomain;
  
  // treat cell
  const Subdomain subdomain = r_oracle_.is_in_domain_object()(r_tr_.dual(cell));
  if ( subdomain )
  {
    set_cell_in_domain(cell, *subdomain);
    
    // Add to refinement queue if needed
    compute_badness(cell);
  }
  else
  {
    remove_cell_from_domain(cell);
  }
}
  

template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
compute_badness(const Cell_handle& cell)
{
  const Is_cell_bad is_cell_bad = r_criteria_(cell);
  if( is_cell_bad )
  {
    this->add_bad_element(cell, *is_cell_bad);
  }
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
typename Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::Vertex_handle
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,C_>::
insert_impl(const Point& point,
            const Zone& zone)
{
  // TODO: look at this
  if( zone.locate_type == Tr::VERTEX )
  {
    return zone.cell->vertex(zone.i);
  }
  
  const Facet& facet = *(zone.boundary_facets.begin());
  
  Vertex_handle v = r_tr_.insert_in_hole(point,
                                         zone.cells.begin(),
                                         zone.cells.end(),
                                         facet.first,
                                         facet.second);
  
  // Set index and dimension of v
  set_vertex_properties(v, last_vertex_index_);
  
  return v;
}

}  // end namespace Mesh_3


}  // end namespace CGAL


#endif // CGAL_MESH_3_REFINE_CELLS_3_H
