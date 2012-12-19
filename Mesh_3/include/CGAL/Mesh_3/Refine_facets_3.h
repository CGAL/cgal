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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implements a mesher level for facets.
//******************************************************************************

#ifndef CGAL_MESH_3_REFINE_FACETS_3_H
#define CGAL_MESH_3_REFINE_FACETS_3_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesher_level_default_implementations.h>
#include <CGAL/Meshes/Double_map_container.h>
#include <CGAL/Meshes/Triangulation_mesher_level_traits_3.h>

#include <CGAL/Object.h>

#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/tuple.h>
#include <sstream>

namespace CGAL {

namespace Mesh_3 {

// Helper meta-programming functions, to allow backward compatibility.
//
//   - Has_Is_facet_bad and Had_Facet_badness are used to detect if a model
//     of the MeshFacetCriteria_3 concept follows the specifications of
//     CGAL-3.7 (with Facet_badness) or later (with Is_facet_bad).
//
//   - Then the meta-function Get_Is_facet_bad is used to get the actual
//     type, either Facet_criteria::Facet_badness or
//      Facet_criteria::Is_facet_bad.

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_Is_facet_bad, Is_facet_bad, true)
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_Facet_badness, Facet_badness, false)

// template class, used when use_facet_badness = false
template <typename Facet_criteria, 
          bool use_facet_badness = (!Has_Is_facet_bad<Facet_criteria>::value) &&
                                    Has_Facet_badness<Facet_criteria>::value >
struct Get_Is_facet_bad {
  typedef typename Facet_criteria::Is_facet_bad Type;
  typedef Type type; // compatibility with Boost
};

// partial specialization when use_facet_badness == true
template <typename Facet_criteria>
struct Get_Is_facet_bad<Facet_criteria, true> {
  typedef typename Facet_criteria::Facet_badness Type;
  typedef Type type;
};


// Class Refine_facets_3
//
// Template parameters should be models of
// Tr         : MeshTriangulation_3
// Criteria   : SurfaceMeshFacetsCriteria_3
// MeshDomain : MeshTraits_3
//
// Implements a Mesher_level for facets

// TODO document Container_ requirements
template<class Tr,
         class Criteria,
         class MeshDomain,
         class Complex3InTriangulation3,
         class Previous_level_,
         class Container_ = Meshes::Double_map_container<
                                            typename Tr::Facet,
                                            typename Criteria::Facet_quality> >
class Refine_facets_3
: public Mesher_level<Tr,
                      Refine_facets_3<Tr,
                                      Criteria,
                                      MeshDomain,
                                      Complex3InTriangulation3,
                                      Previous_level_>,
                      typename Tr::Facet,
                      Previous_level_,
                      Triangulation_mesher_level_traits_3<Tr> >
, public Container_
, public No_after_no_insertion
, public No_before_conflicts
{
  // Self
  typedef Refine_facets_3<Tr,
                          Criteria,
                          MeshDomain,
                          Complex3InTriangulation3,
                          Previous_level_>                  Self;

public:
  typedef typename Tr::Point Point;
  typedef typename Tr::Facet Facet;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
  typedef Complex3InTriangulation3 C3T3;

  /// Constructor
  Refine_facets_3(Tr& triangulation,
                  const Criteria& criteria,
                  const MeshDomain& oracle,
                  Previous_level_& previous,
                  C3T3& c3t3);

  /// Destructor
  virtual ~Refine_facets_3() { }

  /// Get a reference on triangulation
  Tr& triangulation_ref_impl() { return r_tr_; }
  const Tr& triangulation_ref_impl() const { return r_tr_; }

  /// Initialization function
  void scan_triangulation_impl();

  /// Gets the point to insert from the element to refine
  Point refinement_point_impl(const Facet& facet) const
  {
    CGAL_assertion (is_facet_on_surface(facet));
    last_vertex_index_ = get_facet_surface_center_index(facet);
    return get_facet_surface_center(facet);
  };

  /// Tests if p encroaches facet from zone
  Mesher_level_conflict_status
  test_point_conflict_from_superior_impl(const Point& p, Zone& zone);

  /// Useless here
  Mesher_level_conflict_status private_test_point_conflict_impl(
      const Point& p,
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

  /// Returns the conflicts zone
  Zone conflicts_zone_impl(const Point& point, const Facet& facet) const;

  /// Job to do before insertion
  void before_insertion_impl(const Facet& facet,
                             const Point& point,
                             Zone& zone);

  /// Job to do after insertion
  void after_insertion_impl(const Vertex_handle& v)
  { restore_restricted_Delaunay(v); }

  /// Insert p into triangulation
  Vertex_handle insert_impl(const Point& p, const Zone& zone);

  /// Restore restricted Delaunay ; may be call by Cells_mesher visitor
  void restore_restricted_Delaunay(const Vertex_handle& v);

  /// debug info
  std::string debug_info() const
  {
    std::stringstream s;
    s << Container_::size();
    return s.str();
  }

  /// debug_info_header
  std::string debug_info_header() const
  {
    return "#facets to refine";
  }
  
#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  std::size_t queue_size() const { return this->size(); }
#endif


private:
  //-------------------------------------------------------
  // Private types
  //-------------------------------------------------------
  typedef typename Tr::Geom_traits Gt;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Criteria::Facet_quality Quality;
  typedef typename Get_Is_facet_bad<Criteria>::Type Is_facet_bad;
  typedef typename MeshDomain::Surface_patch_index Surface_patch_index;
  typedef typename MeshDomain::Index Index;
  typedef typename Gt::Segment_3 Segment_3;
  typedef typename Gt::Ray_3 Ray_3;
  typedef typename Gt::Line_3 Line_3;

  typedef typename boost::optional<CGAL::cpp11::tuple<Surface_patch_index, Index, Point> >
                                                               Facet_properties;

private:
  /// Get mirror facet
  Facet mirror_facet(const Facet& f) const { return r_tr_.mirror_facet(f); };

  /// Returns canonical facet of facet
  Facet canonical_facet(const Facet& facet) const
  {
    const Facet mirror = mirror_facet(facet);
    return ( (facet < mirror)?facet:mirror );
  }

  /// Returns true if f has already been visited
  bool is_facet_visited(const Facet& f) const
  {
    return f.first->is_facet_visited(f.second);
  }

  /// Sets facet f to visited
  void set_facet_visited(Facet& f) const
  {
    f.first->set_facet_visited(f.second);
  }

  /// Sets facet f to not visited
  void reset_facet_visited(Facet& f) const
  {
    f.first->reset_visited(f.second);
  }

  /// Sets facet f and it's mirror one surface center to point p
  void set_facet_surface_center(const Facet& f,
                                const Point& p,
                                const Index& index) const
  {
    const Facet mirror = mirror_facet(f);

    f.first->set_facet_surface_center(f.second, p);
    mirror.first->set_facet_surface_center(mirror.second, p);

    f.first->set_facet_surface_center_index(f.second,index);
    mirror.first->set_facet_surface_center_index(mirror.second,index);
  }

  /// Returns facet surface center of \c f
  Point get_facet_surface_center(const Facet& f) const
  {
    return f.first->get_facet_surface_center(f.second);
  }

  /// Returns index of surface center of facet \c f
  Index get_facet_surface_center_index(const Facet& f) const
  {
    return f.first->get_facet_surface_center_index(f.second);
  }

  /// Sets \c f to surface facets, with index \c index
  void set_facet_on_surface(const Facet& f,
                            const Surface_patch_index& index)
  {
    r_c3t3_.add_to_complex(f, index);
  }

  /// Returns index of facet \c f
  Surface_patch_index get_facet_surface_index(const Facet& f) const
  {
    return r_c3t3_.surface_patch_index(f.first, f.second);
  }

  /// Removes \c f from surface facets
  void remove_facet_from_surface(const Facet& f)
  {
    r_c3t3_.remove_from_complex(f);
  }

  /// Returns to if f is on surface
  bool is_facet_on_surface(const Facet& f) const
  {
    return r_c3t3_.is_in_complex(f) ;
  }

  /// Sets index and dimension of vertex \v
  void set_vertex_properties(Vertex_handle& v, const Index& index)
  {
    r_c3t3_.set_index(v, index);
    // Set dimension of v: v is on surface by construction, so dimension=2
    v->set_dimension(2);
  }

  /// Computes facet properties and add facet to the refinement queue if needed
  void treat_new_facet(Facet& facet);

  /**
   * Computes at once is_facet_on_surface and facet_surface_center.
   * @param facet The input facet
   * @return \c true if \c facet is on surface, \c false otherwise
   */
  Facet_properties compute_facet_properties(const Facet& facet) const;

  /// Returns true if point encroaches facet
  bool is_facet_encroached(const Facet& facet, const Point& point) const;

  /// Returns whethere an encroached facet is refinable or not
  bool is_encroached_facet_refinable(Facet& facet) const;

  /// Insert facet into refinement queue
  void insert_bad_facet(Facet& facet, const Quality& quality)
  {
    // Insert canonical facet
    this->add_bad_element(this->canonical_facet(facet), quality);
  }
  
  /// Insert encroached facet into refinement queue
  void insert_encroached_facet_in_queue(Facet& facet)
  {
    insert_bad_facet(facet,Quality());
  }

  /// Removes facet from refinement queue
  void remove_bad_facet(Facet& facet)
  {
    // Remove canonical facet
    Facet canonical_facet = this->canonical_facet(facet);
    this->remove_element(canonical_facet);
  }

  /**
   * Action to perform on a facet inside the conflict zone before insertion
   * @return true if source_facet is the same as facet or mirror(facet)
   */
  bool before_insertion_handle_facet_in_conflict_zone(Facet& facet,
                                                     const Facet& source_facet);

  /**
   * Action to perform on a facet on the boundary of the conflict zone
   * before insertion
   * @return true if source_facet is the same as facet or mirror(facet)
   */
  bool before_insertion_handle_facet_on_conflict_zone(Facet& facet,
                                                      const Facet& source_facet)
  {
    // perform the same operations as for an internal facet
    return before_insertion_handle_facet_in_conflict_zone(facet, source_facet);
  }

  /// Action to perform on a facet incident to the new vertex
  void after_insertion_handle_incident_facet(Facet& facet);

  /// Action to perform on a facet opposite to the new vertex
  void after_insertion_handle_opposite_facet(Facet& facet)
  {
    // perform the same operations as for a facet incident to the new vertex
    after_insertion_handle_incident_facet(facet);
  }

private:
  /// The triangulation
  Tr& r_tr_;
  /// The facet criteria
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
  Refine_facets_3(const Self& src);
  // Disabled assignment operator
  Self& operator=(const Self& src);
};  // end class Refine_facets_3





template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
Refine_facets_3(Tr& triangulation,
                const Cr& criteria,
                const MD& oracle,
                P_& previous,
                C3T3& c3t3)
  : Mesher_level<Tr, Self, Facet, P_,
                               Triangulation_mesher_level_traits_3<Tr> >(previous)
  , C_()
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
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
scan_triangulation_impl()
{
  typedef typename Tr::Finite_facets_iterator Finite_facet_iterator;

  for(Finite_facet_iterator facet_it = r_tr_.finite_facets_begin();
      facet_it != r_tr_.finite_facets_end();
      ++facet_it)
  {
    // Cannot be const, see treat_new_facet signature
    Facet facet = *facet_it;
    treat_new_facet(facet);
  }
}




template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
Mesher_level_conflict_status
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
test_point_conflict_from_superior_impl(const Point& point,
                                       Zone& zone)
{
  typedef typename Zone::Facets_iterator Facet_iterator;

  for (Facet_iterator facet_it = zone.internal_facets.begin();
       facet_it != zone.internal_facets.end();
       ++facet_it)
  {
    // surface facets which are internal facets of the conflict zone are
    // encroached
    if( is_facet_on_surface(*facet_it) )
    {
      if ( is_encroached_facet_refinable(*facet_it) )
      {
        insert_encroached_facet_in_queue(*facet_it);
        return CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;
      }
      else
        return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
    }
  }

  for (Facet_iterator facet_it = zone.boundary_facets.begin();
       facet_it != zone.boundary_facets.end();
       ++facet_it)
  {
    if( is_facet_encroached(*facet_it, point) )
    {
      // Insert already existing surface facet into refinement queue
      if ( is_encroached_facet_refinable(*facet_it) )
      {
        insert_encroached_facet_in_queue(*facet_it);
        return CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;
      }
      else
        return CONFLICT_AND_ELEMENT_SHOULD_BE_DROPPED;
    }
  }

  return NO_CONFLICT;
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
typename Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::Zone
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
conflicts_zone_impl(const Point& point,
                    const Facet& facet) const
{
  Zone zone;

  // TODO may be avoid the locate here
  zone.cell = r_tr_.locate(point,
                           zone.locate_type,
                           zone.i,
                           zone.j,
                           facet.first);

  if(zone.locate_type != Tr::VERTEX)
  {
    r_tr_.find_conflicts(point,
                         zone.cell,
                         std::back_inserter(zone.boundary_facets),
                         std::back_inserter(zone.cells),
                         std::back_inserter(zone.internal_facets));
  }

  return zone;
}




template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
before_insertion_impl(const Facet& facet,
                      const Point& point,
                      Zone& zone)
{
  typedef typename Zone::Facets_iterator Facets_iterator;

  bool source_facet_is_in_conflict = false;

  // Iterate on conflict zone facets
  for (Facets_iterator facet_it = zone.internal_facets.begin();
       facet_it != zone.internal_facets.end();
       ++facet_it)
  {
    if (before_insertion_handle_facet_in_conflict_zone(*facet_it, facet) )
    {
      source_facet_is_in_conflict = true;
    }
  }

  for (Facets_iterator facet_it = zone.boundary_facets.begin() ;
       facet_it != zone.boundary_facets.end() ;
       ++facet_it)
  {
    if (before_insertion_handle_facet_on_conflict_zone(*facet_it, facet))
    {
      source_facet_is_in_conflict = true;
    }
  }

  // source_facet == Facet() when this->before_insertion_impl is
  // called from a Mesh_3 visitor.
  if ( !source_facet_is_in_conflict && facet != Facet()  )
  {
    const Facet source_other_side = mirror_facet(facet);
    std::stringstream error_msg;
    error_msg <<
      boost::format("Mesh_3 ERROR: "
                    "A facet is not in conflict with its refinement point!\n"
                    "Debugging informations:\n"
                    "  Facet: (%1%, %2%) = (%6%, %7%, %8%)\n"
                    "  Dual: (%3%, %4%)\n"
                    "  Refinement point: %5%\n")
      % (&*facet.first)
      % facet.second
      % triangulation_ref_impl().dual(facet.first)
      % triangulation_ref_impl().dual(source_other_side.first)
      % point
      % facet.first->vertex((facet.second + 1)&3)->point()
      % facet.first->vertex((facet.second + 2)&3)->point()
      % facet.first->vertex((facet.second + 3)&3)->point();

    CGAL_error_msg(error_msg.str().c_str());
  }
}



template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
typename Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::Vertex_handle
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
insert_impl(const Point& point,
            const Zone& zone)
{
  if( zone.locate_type == Tr::VERTEX )
  {
    // TODO look at this
    std::cerr<<"VERTEX\n";
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



template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
restore_restricted_Delaunay(const Vertex_handle& vertex)
{
  typedef std::vector<Cell_handle> Cell_handle_vector;
  typedef typename Cell_handle_vector::iterator Cell_handle_vector_iterator;

  // Update incident facets
  Cell_handle_vector cells;
  r_tr_.incident_cells(vertex, std::back_inserter(cells));

  for(Cell_handle_vector_iterator cell_it = cells.begin();
      cell_it != cells.end();
      ++cell_it)
  {
    // Look at all four facets of the cell, starting with the
    // facet opposite to the new vertex
    int index = (*cell_it)->index(vertex);

    Facet opposite_facet(*cell_it, index);
    after_insertion_handle_opposite_facet(opposite_facet);

    for (int i = 1; i <= 3; ++i)
    {
      Facet incident_facet(*cell_it, (index+i)&3);
      after_insertion_handle_incident_facet(incident_facet);
    }
  }
}

//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
treat_new_facet(Facet& facet)
{
  // Treat facet
  Facet_properties properties = compute_facet_properties(facet);
  if ( properties )
  {
    const Surface_patch_index& surface_index = CGAL::cpp11::get<0>(*properties);
    const Index& surface_center_index = CGAL::cpp11::get<1>(*properties);
    const Point& surface_center = CGAL::cpp11::get<2>(*properties);

    // Facet is on surface: set facet properties
    set_facet_surface_center(facet, surface_center, surface_center_index);
    set_facet_on_surface(facet, surface_index);

    // Insert facet into refinement queue if needed
    const Is_facet_bad is_facet_bad = r_criteria_(facet);
    if ( is_facet_bad )
    {
      insert_bad_facet(facet, *is_facet_bad);
    }
  }
  else
  {
    // Facet is not on surface
    remove_facet_from_surface(facet);
  }

  // Set facet visited
  Facet mirror = mirror_facet(facet);
  set_facet_visited(facet);
  set_facet_visited(mirror);
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
typename Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::Facet_properties
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
compute_facet_properties(const Facet& facet) const
{
  //-------------------------------------------------------
  // Facet must be finite
  //-------------------------------------------------------
  CGAL_assertion( ! r_tr_.is_infinite(facet) );

  // types
  typedef boost::optional<typename MD::Surface_patch_index> Surface_patch;
  typedef typename MD::Intersection Intersection;

  // Functor
  typename Gt::Is_degenerate_3 is_degenerate = Gt().is_degenerate_3_object();
  typename Gt::Compare_xyz_3 compare_xyz = Gt().compare_xyz_3_object();
  typename MD::Do_intersect_surface do_intersect_surface =
      r_oracle_.do_intersect_surface_object();


  // Get dual of facet
  Object dual = r_tr_.dual(facet);

  // If the dual is a segment
  if ( const Segment_3* p_segment = object_cast<Segment_3>(&dual) )
  {
    if (is_degenerate(*p_segment)) { return Facet_properties(); }

    // If facet is on surface, compute intersection point and return true
    Surface_patch surface = do_intersect_surface(*p_segment);
    if ( surface )
    {
      typename MD::Construct_intersection construct_intersection =
          r_oracle_.construct_intersection_object();

      // Trick to have canonical vector : thus, we compute alwais the same
      // intersection
      Segment_3 segment = *p_segment;
      if ( compare_xyz(p_segment->source(),p_segment->target())
              == CGAL::LARGER )
      {
        typename Gt::Construct_opposite_segment_3 opposite =
            Gt().construct_opposite_segment_3_object();

        segment = opposite(*p_segment);
      }

      Intersection intersect = construct_intersection(segment);
      return Facet_properties(CGAL::cpp11::make_tuple(*surface,
                                                CGAL::cpp11::get<1>(intersect),
                                                CGAL::cpp11::get<0>(intersect)));
    }
  }
  // If the dual is a ray
  else if ( const Ray_3* p_ray = object_cast<Ray_3>(&dual) )
  {
    // If a facet is on the convex hull, and if its finite incident
    // cell has a very bid Delaunay ball, then the dual of the facet is
    // a ray constructed with a point with very big coordinates, and a
    // vector with small coordinates. Its can happen than the
    // constructed ray is degenerate (the point(1) of the ray is
    // point(0) plus a vector whose coordinates are espilon).
    if (is_degenerate(*p_ray)) { return Facet_properties(); }

    Surface_patch surface = do_intersect_surface(*p_ray);
    if ( surface )
    {
      typename MD::Construct_intersection construct_intersection =
          r_oracle_.construct_intersection_object();

      Intersection intersect = construct_intersection(*p_ray);
      return Facet_properties(CGAL::cpp11::make_tuple(*surface,
                                                CGAL::cpp11::get<1>(intersect),
                                                CGAL::cpp11::get<0>(intersect)));
    }
  }
  // If the dual is a line
  else if ( const Line_3* p_line = object_cast<Line_3>(&dual) )
  {
    Surface_patch surface = do_intersect_surface(*p_line);
    if ( surface )
    {
      typename MD::Construct_intersection construct_intersection =
          r_oracle_.construct_intersection_object();

      // Trick to have canonical vector : thus, we compute alwais the same
      // intersection
      Line_3 line = *p_line;
      typename Gt::Compare_xyz_3 compare_xyz = Gt().compare_xyz_3_object();
      if ( compare_xyz(p_line->point(0),p_line->point(1))
              == CGAL::LARGER )
      {
        typename Gt::Construct_opposite_line_3 opposite =
            Gt().construct_opposite_line_3_object();

        line = opposite(*p_line);
      }

      Intersection intersect = construct_intersection(line);
      return Facet_properties(CGAL::cpp11::make_tuple(*surface,
                                                CGAL::cpp11::get<1>(intersect),
                                                CGAL::cpp11::get<0>(intersect)));
    }
  }
  else
  {
    // Else there is a problem with the dual
    std::cerr << "In is_facet_on_surface(const Facet& f, Point& center)\n"
    << "file " << __FILE__ << ", line " << __LINE__ << "\n";
    std::cerr << "Incorrect object type: " << dual.type().name() << "\n";
    CGAL_error();
  }

  return Facet_properties();
}


template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
bool
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
is_facet_encroached(const Facet& facet,
                    const Point& point) const
{
  if ( r_tr_.is_infinite(facet) || ! is_facet_on_surface(facet) )
  {
    return false;
  }
  
  typename Gt::Compare_power_distance_3 compare_distance =
    r_tr_.geom_traits().compare_power_distance_3_object();
  
  const Cell_handle& cell = facet.first;
  const int& facet_index = facet.second;
  const Point& center = get_facet_surface_center(facet);
  const Point& reference_point = cell->vertex((facet_index+1)&3)->point();
  
  // facet is encroached if the new point is near from center than
  // one vertex of the facet
  return ( compare_distance(center, reference_point, point) != CGAL::SMALLER );
}

template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
bool
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
is_encroached_facet_refinable(Facet& facet) const
{
  typedef typename Gt::Point_3 Point_3;
  typedef typename Gt::FT      FT;
  
  typename Gt::Compute_squared_radius_smallest_orthogonal_sphere_3 sq_radius =
    Gt().compute_squared_radius_smallest_orthogonal_sphere_3_object();
  
  typename Gt::Compare_weighted_squared_radius_3 compare =
    Gt().compare_weighted_squared_radius_3_object();
  
  const Cell_handle& c = facet.first;
  const int& k = facet.second;
  
  int k1 = (k+1)&3;
  int k2 = (k+2)&3;
  int k3 = (k+3)&3;
  
  // Get number of weighted points, and ensure that they will be accessible
  // using k1...ki, if i is the number of weighted points.
  int wp_nb = 0;
  if(c->vertex(k1)->point().weight() > FT(0))
  { 
    ++wp_nb;
  }
  
  if(c->vertex(k2)->point().weight() > FT(0))
  { 
    if ( 0 == wp_nb ) { std::swap(k1,k2); }
    ++wp_nb;
  }
  
  if(c->vertex(k3)->point().weight() > FT(0))
  { 
    if ( 0 == wp_nb ) { std::swap(k1,k3); }
    if ( 1 == wp_nb ) { std::swap(k2,k3); }
    ++wp_nb;
  }
  
  const Point_3& p1 = c->vertex(k1)->point();
  const Point_3& p2 = c->vertex(k2)->point();
  const Point_3& p3 = c->vertex(k3)->point();
  
  const FT min_ratio (0.16); // (0.2*2)^2
  
  // Check ratio
  switch ( wp_nb )
  {
    case 1:
    {
      FT r = (std::max)(sq_radius(p1,p2),sq_radius(p1,p3));
      if ( r < min_ratio*p1.weight() ) { return false; }
      break;
    }
      
    case 2:
    {
      bool do_spheres_intersect = (compare(p1,p2,FT(0)) != CGAL::LARGER);
      if ( do_spheres_intersect )
      {
        FT r13 = sq_radius(p1,p3) / p1.weight();
        FT r23 = sq_radius(p2,p3) / p2.weight();
        FT r = (std::max)(r13,r23);
        
        if ( r < min_ratio ) { return false; }
      }
      break;
    }
      
    case 3:
    {
      bool do_spheres_intersect = (compare(p1,p2,p3,FT(0)) != CGAL::LARGER);
      if ( do_spheres_intersect ) { return false; }
      break;
    }  
      
    default: break;
  }

  return true;
}

 
template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
bool
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
before_insertion_handle_facet_in_conflict_zone(Facet& facet,
                                               const Facet& source_facet)
{
  Facet other_side = mirror_facet(facet);

  if ( is_facet_on_surface(facet) )
  {
    // Remove facet from refinement queue
    remove_bad_facet(facet);

    // Remove facet from complex
    remove_facet_from_surface(facet);

    // Reset visited
    reset_facet_visited(facet);
    reset_facet_visited(other_side);
  }

  return ( (facet == source_facet) || (other_side == source_facet) );
}



template<class Tr, class Cr, class MD, class C3T3_, class P_, class C_>
void
Refine_facets_3<Tr,Cr,MD,C3T3_,P_,C_>::
after_insertion_handle_incident_facet(Facet& facet)
{
  // If the facet is infinite or has been already visited,
  // then there is nothing to do as for it or its edges
  // Facet other_side = mirror_facet(facet);
  if ( r_tr_.is_infinite(facet) || (is_facet_visited(facet)) ) // && is_facet_visited(other_side)) )
  {
    return;
  }

  treat_new_facet(facet);
}


}  // end namespace Mesh_3


}  // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_FACETS_3_H
