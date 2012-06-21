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
#ifdef CGAL_LINKED_WITH_TBB
  #include <tbb/tbb.h>
#endif

#include <CGAL/Meshes/Filtered_deque_container.h>
#include <CGAL/Meshes/Filtered_multimap_container.h>
#include <CGAL/Meshes/Double_map_container.h>

#ifdef MESH_3_PROFILING
  #include <CGAL/Mesh_3/Profiling_tools.h>
#endif

#include <boost/format.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/if.hpp>
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

  // Predicate to know if a cell in a refinement queue is a zombie
  template<typename Cell_handle>
  class Cell_to_refine_is_not_zombie
  {
  public:
    Cell_to_refine_is_not_zombie() {}

    bool operator()(const std::pair<Cell_handle, unsigned int> &c) const
    {
      return (c.first->get_erase_counter() == c.second);
    }
  };
  
/************************************************
// Class Refine_cells_3_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename Index, typename Concurrency_tag>
class Refine_cells_3_base
{
protected:
  Refine_cells_3_base() : m_last_vertex_index() {}

  Index get_last_vertex_index() const
  {
    return m_last_vertex_index;
  }

  void set_last_vertex_index(Index i) const
  {
    m_last_vertex_index = i;
  }
  
#if defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE) \
 || defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
  template <typename Cell_handle>
  std::pair<Cell_handle, unsigned int> 
  from_cell_to_refinement_queue_element(Cell_handle ch) const
  {
    return std::make_pair(ch, ch->get_erase_counter());
  }
#else
  template <typename Cell_handle>
  Cell_handle
  from_cell_to_refinement_queue_element(Cell_handle ch) const
  {
    return ch;
  }
#endif

  /// Stores index of vertex that may be inserted into triangulation
  mutable Index m_last_vertex_index;
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename Index>
class Refine_cells_3_base<Index, Parallel_tag>
{
protected:
  Refine_cells_3_base() : m_last_vertex_index(Index()) {}

  Index get_last_vertex_index() const
  {
    return m_last_vertex_index.local();
  }

  void set_last_vertex_index(Index i) const
  {
    m_last_vertex_index.local() = i;
  }

  template <typename Cell_handle>
  std::pair<Cell_handle, unsigned int> 
  from_cell_to_refinement_queue_element(Cell_handle ch) const
  {
    return std::make_pair(ch, ch->get_erase_counter());
  }

  /// Stores index of vertex that may be inserted into triangulation
  mutable tbb::enumerable_thread_specific<Index> m_last_vertex_index;
};
#endif // CGAL_LINKED_WITH_TBB

/************************************************
// Class Refine_cells_3
//
// Template parameters should be models of
// Tr         : MeshTriangulation_3
// Criteria   : MeshCellsCriteria_3
// MeshDomain : MeshTraits_3
//
// Implements a Mesher_level for cells
************************************************/

template<class Tr,
         class Criteria,
         class MeshDomain,
         class Complex3InTriangulation3,
         class Previous_,
         class Concurrency_tag,
#ifdef CGAL_LINKED_WITH_TBB
         class Container_ = typename boost::mpl::if_c // (parallel/sequential?)
         <
          boost::is_base_of<Parallel_tag, Concurrency_tag>::value,

          // Parallel
# ifdef CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
          Meshes::Filtered_deque_container
# else
          Meshes::Filtered_multimap_container
# endif
          <
            std::pair<typename Tr::Cell_handle, unsigned int>,
            typename Criteria::Cell_quality,
            Cell_to_refine_is_not_zombie<typename Tr::Cell_handle>,
            Concurrency_tag
          >,
          // Sequential
# ifdef CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
          Meshes::Filtered_deque_container
          <
            std::pair<typename Tr::Cell_handle, unsigned int>,
            typename Criteria::Cell_quality,
            Cell_to_refine_is_not_zombie<typename Tr::Cell_handle>,
            Concurrency_tag
          >
# elif defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE)
          Meshes::Filtered_multimap_container
          <
            std::pair<typename Tr::Cell_handle, unsigned int>,
            typename Criteria::Cell_quality,
            Cell_to_refine_is_not_zombie<typename Tr::Cell_handle>,
            Concurrency_tag
          >
# else
          Meshes::Double_map_container<typename Tr::Cell_handle,
                                       typename Criteria::Cell_quality>
# endif
         >::type // boost::if (parallel/sequential)

#else // !CGAL_LINKED_WITH_TBB

        // Sequential
        class Container_ = 
# ifdef CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE
          Meshes::Filtered_deque_container
          <
            std::pair<typename Tr::Cell_handle, unsigned int>,
            typename Criteria::Cell_quality,
            Cell_to_refine_is_not_zombie<typename Tr::Cell_handle>,
            Concurrency_tag
          >
# elif defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE)
          Meshes::Filtered_multimap_container
          <
            std::pair<typename Tr::Cell_handle, unsigned int>,
            typename Criteria::Cell_quality,
            Cell_to_refine_is_not_zombie<typename Tr::Cell_handle>,
            Concurrency_tag
          >
# else
          Meshes::Double_map_container<typename Tr::Cell_handle,
                                       typename Criteria::Cell_quality>
# endif

#endif // CGAL_LINKED_WITH_TBB
>
class Refine_cells_3
: public Refine_cells_3_base<typename MeshDomain::Index, Concurrency_tag>
, public Mesher_level<Tr,
                      Refine_cells_3<Tr,
                                      Criteria,
                                      MeshDomain,
                                      Complex3InTriangulation3,
                                      Previous_,
                                      Concurrency_tag,
                                      Container_>,
                      typename Tr::Cell_handle,
                      Previous_,
                      Triangulation_mesher_level_traits_3<Tr>,
                      Concurrency_tag>
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
                         Concurrency_tag,
                         Container_>       Self;
  
public:  
  typedef Container_ Container; // Because we need it in Mesher_level
  typedef typename Container::Element Container_element;
  typedef typename Tr::Point Point;
  typedef typename Tr::Cell Cell;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Criteria::Cell_quality Cell_quality;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
  typedef Complex3InTriangulation3 C3T3;
  
  
  // Constructor
  // For sequential
  Refine_cells_3(Tr& triangulation,
                 const Criteria& criteria,
                 const MeshDomain& oracle,
                 Previous_& previous,
                 C3T3& c3t3);
  // For parallel
  Refine_cells_3(Tr& triangulation,
                 const Criteria& criteria,
                 const MeshDomain& oracle,
                 Previous_& previous,
                 C3T3& c3t3,
                 Mesh_3::LockDataStructureType *p_lock_ds,
                 Mesh_3::WorksharingDataStructureType *p_worksharing_ds);
  
  // Destructor
  virtual ~Refine_cells_3() { }
  
  // Get a reference on triangulation
  Tr& triangulation_ref_impl() { return r_tr_; }
  const Tr& triangulation_ref_impl() const { return r_tr_; }
  
  // Initialization function
  void scan_triangulation_impl();
  
  int get_number_of_bad_elements_impl();
  
  Point circumcenter_impl(const Cell_handle& cell) const
  {
    return r_tr_.dual(cell);
  }
  
  template <typename Mesh_visitor>
  void before_next_element_refinement_in_superior_impl(Mesh_visitor)
  {
  }

  void before_next_element_refinement_impl() 
  {
  }

  Cell_handle extract_element_from_container_value(const Container_element &e) const
  {
    // We get the Cell_handle from the pair
    return e.first;
  }
  
  Cell_handle get_next_element_impl() const
  {
    return extract_element_from_container_value(Container_::get_next_element_impl());
  }

  // Gets the point to insert from the element to refine
  Point refinement_point_impl(const Cell_handle& cell) const
  {
    set_last_vertex_index(
      r_oracle_.index_from_subdomain_index(cell->subdomain_index()) );
    
    //    last_vertex_index_ = Index(cell->subdomain_index());
    // NB : dual() is optimized when the cell base class has circumcenter()
    return r_tr_.dual(cell);
  }
  
  // Returns the conflicts zone
  Zone conflicts_zone_impl(const Point& point
                           , const Cell_handle& cell
                           , bool &facet_not_in_its_cz) const;
  Zone conflicts_zone_impl(const Point& point
                           , const Cell_handle& cell
                           , bool &facet_not_in_its_cz
                           , bool &could_lock_zone) const;
  
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
  
  // Sequential
  void remove_element_from_refinement_queue(Cell_handle c, Sequential_tag)
  {
    // If sequential AND NOT lazy, remove cell from refinement queue
  #if !defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE) \
   && !defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
    this->remove_element(c);
  #endif
  }
  // Parallel: it's always lazy, so do nothing
  void remove_element_from_refinement_queue(Cell_handle, Parallel_tag) {}

  /// Handle cells contained in \c zone (before their destruction by insertion)
  void before_insertion_handle_cells_in_conflict_zone(Zone& zone);
  
  /// debug info: class name
  std::string debug_info_class_name_impl() const
  {
    return "Refine_cells_3";
  }

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

  std::string debug_info_element_impl(const Cell_handle &ch) const
  {
    std::stringstream sstr;
    sstr << "Cell { " << std::endl
    << "  - " << *ch->vertex(0)  << std::endl
    << "  - " << *ch->vertex(1)  << std::endl
    << "  - " << *ch->vertex(2)  << std::endl
    << "  - " << *ch->vertex(3)  << std::endl
    << "}" << std::endl;

    return sstr.str();
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
  
private:
  // Disabled copy constructor
  Refine_cells_3(const Self& src);
  // Disabled assignment operator
  Self& operator=(const Self& src);
  
};  // end class Refine_cells_3



// For sequential
template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
Refine_cells_3(Tr& triangulation,
               const Cr& criteria,
               const MD& oracle,
               P_& previous,
               C3T3& c3t3)
  : Mesher_level<Tr, Self, Cell_handle, P_,
      Triangulation_mesher_level_traits_3<Tr>, Ct >(previous)
  , C_()
  , No_test_point_conflict()
  , No_after_no_insertion()
  , No_before_conflicts()
  , r_tr_(triangulation)
  , r_criteria_(criteria)
  , r_oracle_(oracle)
  , r_c3t3_(c3t3)
{
}


// For parallel
template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
Refine_cells_3(Tr& triangulation,
               const Cr& criteria,
               const MD& oracle,
               P_& previous,
               C3T3& c3t3,
               Mesh_3::LockDataStructureType *p_lock_ds,
               Mesh_3::WorksharingDataStructureType *p_worksharing_ds)
  : Mesher_level<Tr, Self, Cell_handle, P_,
      Triangulation_mesher_level_traits_3<Tr>, Ct >(previous, p_lock_ds, p_worksharing_ds)
  , C_()
  , No_test_point_conflict()
  , No_after_no_insertion()
  , No_before_conflicts()
  , r_tr_(triangulation)
  , r_criteria_(criteria)
  , r_oracle_(oracle)
  , r_c3t3_(c3t3)
{
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
scan_triangulation_impl()
{
  typedef typename Tr::Finite_cells_iterator Finite_cell_iterator;
  typedef typename Tr::All_cells_iterator All_cells_iterator;
  
#ifdef MESH_3_PROFILING
  // Refinement done
  std::cerr << "done in " << m_timer.elapsed() << " seconds." << std::endl;
  WallClockTimer t;
#endif

  
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_base_of<Parallel_tag, Ct>::value)
    {
    std::cerr << "Scanning triangulation for bad cells (in parallel)...";
    add_to_TLS_lists(true);
  
    // WITH PARALLEL_FOR
  
    //WallClockTimer t2;

    std::vector<Cell_handle> cells;
    for(All_cells_iterator cell_it = r_tr_.all_cells_begin();
        cell_it != r_tr_.all_cells_end();
        ++cell_it)
    {
      cells.push_back(cell_it);
    }

    //std::cerr << "Parallel_for - push_backs done: " << t2.elapsed() << " seconds." << std::endl;
    //t2.reset();

    std::cerr << "Num cells to scan: " << cells.size() << std::endl;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, cells.size(), 1000),
      [&]( const tbb::blocked_range<size_t>& r ) { // CJTODO: lambdas ok?
        for( size_t i = r.begin() ; i != r.end() ; ++i)
        {
          Cell_handle c = cells[i];
          if (!r_tr_.is_infinite(c))
            treat_new_cell(c);
        }
    });

    //std::cerr << "Parallel_for - iterations done: " << t2.elapsed() << " seconds." << std::endl;
    //t2.reset();  

    // WITH PARALLEL_DO
    /*tbb::parallel_do(r_tr_.finite_cells_begin(), r_tr_.finite_cells_end(),
      [=]( Cell &cell ) { // CJTODO: lambdas ok?
        // CJTODO: should use Compact_container::s_iterator_to, 
        // but we don't know the exact Compact_container type here
        Cell_handle c(&cell);
        treat_new_cell( c );
    });*/

    splice_local_lists();
    //std::cerr << "Parallel_for - splice done: " << t2.elapsed() << " seconds." << std::endl;
    add_to_TLS_lists(false);
  }
  // Sequential
  else
#endif // CGAL_LINKED_WITH_TBB
  {
    std::cerr << "Scanning triangulation for bad cells (sequential)... ";

    int count = 0;
    for(Finite_cell_iterator cell_it = r_tr_.finite_cells_begin();
        cell_it != r_tr_.finite_cells_end();
        ++cell_it)
    {
      treat_new_cell(cell_it);
      ++count;
    }
    std::cerr << count << " cells scanned, ";
  }

#ifdef MESH_3_PROFILING
  double cell_scan_time = t.elapsed();
  std::cerr << "done in " << cell_scan_time << " seconds." << std::endl;
# ifdef CGAL_MESH_3_EXPORT_PERFORMANCE_DATA
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Cells_scan_time", cell_scan_time);
# endif
  std::cerr << "Refining... ";
#else
  std::cerr << "done." << std::endl;
#endif
  
  std::cerr << "Number of bad cells: " << size() << std::endl;
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
int
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
get_number_of_bad_elements_impl()
{
  typedef typename MD::Subdomain Subdomain;
  typedef typename Tr::Finite_cells_iterator Finite_cell_iterator;
  
  int count = 0;
  for(Finite_cell_iterator cell_it = r_tr_.finite_cells_begin();
      cell_it != r_tr_.finite_cells_end();
      ++cell_it)
  {
    // treat cell
    const Subdomain subdomain = r_oracle_.is_in_domain_object()(r_tr_.dual(cell_it));
    if ( subdomain )
    {
      const Is_cell_bad is_cell_bad = r_criteria_(cell_it);
      if( is_cell_bad )
        ++count;
    }
  }

  return count;
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
typename Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::Zone
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
conflicts_zone_impl(const Point& point
                    , const Cell_handle& cell
                    , bool &facet_not_in_its_cz) const
{
  Zone zone;
  zone.cell = cell;
  zone.locate_type = Tr::CELL;
  
  r_tr_.find_conflicts(point,
                       zone.cell,
                       std::back_inserter(zone.boundary_facets),
                       std::back_inserter(zone.cells),
                       std::back_inserter(zone.internal_facets));

  facet_not_in_its_cz = false; // Always false

  return zone;
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
typename Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::Zone
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
conflicts_zone_impl(const Point& point
                    , const Cell_handle& cell
                    , bool &facet_not_in_its_cz
                    , bool &could_lock_zone
     ) const
{
  Zone zone;
  zone.cell = cell;
  zone.locate_type = Tr::CELL;
  
  r_tr_.find_conflicts(point,
                       zone.cell,
                       std::back_inserter(zone.boundary_facets),
                       std::back_inserter(zone.cells),
                       std::back_inserter(zone.internal_facets),
                       &could_lock_zone);

  facet_not_in_its_cz = false; // Always false

  return zone;
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
before_insertion_handle_cells_in_conflict_zone(Zone& zone)
{
  typename Zone::Cells_iterator cit = zone.cells.begin();
  for ( ; cit != zone.cells.end() ; ++cit )
  {
    // Remove element (if needed - see 
    // remove_element_from_refinement_queue implementation)
    this->remove_element_from_refinement_queue(*cit, Ct());
    
    // Remove cell from complex
    remove_cell_from_domain(*cit);
  }
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
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

  
template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
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


template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
treat_new_cell(const Cell_handle& cell)
{
  typedef typename MD::Subdomain Subdomain;
  
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

template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
void
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
compute_badness(const Cell_handle& cell)
{
  const Is_cell_bad is_cell_bad = r_criteria_(cell);
  if( is_cell_bad )
  {
    this->add_bad_element(from_cell_to_refinement_queue_element(cell), *is_cell_bad);
  }
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class Ct, class C_>
typename Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::Vertex_handle
Refine_cells_3<Tr,Cr,MD,C3T3_,P_,Ct,C_>::
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
  set_vertex_properties(v, get_last_vertex_index());  
  return v;
}

}  // end namespace Mesh_3


}  // end namespace CGAL


#endif // CGAL_MESH_3_REFINE_CELLS_3_H
