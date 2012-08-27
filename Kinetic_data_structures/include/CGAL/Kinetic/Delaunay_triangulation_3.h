// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_KINETIC_DELAUNAY_3_H
#define CGAL_KINETIC_KINETIC_DELAUNAY_3_H

#include <CGAL/Kinetic/basic.h>

#include <CGAL/Kinetic/internal/Delaunay_triangulation_base_3.h>

#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

// Triangulations
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

// Local helpers
#include <CGAL/Kinetic/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_3.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4355) // complaint about using 'this' to
#endif                          // initialize a member

namespace CGAL { namespace Kinetic { namespace internal {

template <class Traits>
struct Delaunay_triangulation_3_types
{
  typedef typename Traits::Active_points_3_table MPT;
  typedef typename Traits::Kinetic_kernel KK;
  typedef CGAL::Kinetic::Delaunay_triangulation_cell_base_3<Traits> CFBI;
  /*typedef CGAL::Triangulation_cell_base_with_info_3<Delaunay_cache_3<MPT, KK>,
    typename Traits::Instantaneous_kernel, CFB> CFBI;*/
  typedef CGAL::Triangulation_vertex_base_3<typename Traits::Instantaneous_kernel> CVB;
  typedef CGAL::Triangulation_data_structure_3<CVB, CFBI> TDS;

  typedef CGAL::Delaunay_triangulation_3<typename Traits::Instantaneous_kernel, TDS> Default_triangulation;

  //friend class CGAL::Delaunay_triangulation_3<typename P::Instantaneous_kernel, TDS>;

};

} } } //namespace CGAL::Kinetic::internal

namespace CGAL { namespace Kinetic {

//! A 3D kinetic Delaunay triangulation.
template <class TraitsT,
	  class Visitor= Delaunay_triangulation_visitor_base_3,
	  class TriangulationT=typename internal::Delaunay_triangulation_3_types<TraitsT>::Default_triangulation>
class Delaunay_triangulation_3: public Ref_counted<Delaunay_triangulation_3<TraitsT, Visitor, TriangulationT> > {
private:
  typedef Delaunay_triangulation_3<TraitsT, Visitor, TriangulationT> This_DT3;
  typedef Delaunay_triangulation_3<TraitsT, Visitor, TriangulationT> This;
public:
  typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3::result_type Root_stack;
  typedef typename TriangulationT::Facet Facet;
  typedef typename TriangulationT::Edge Edge;
  typedef typename TraitsT::Simulator::Event_key Event_key;
  typedef typename TraitsT::Active_points_3_table::Key Point_key;
private:
  struct Base_traits: public TraitsT {
    typedef TriangulationT Triangulation;
    typedef typename TraitsT::Kinetic_kernel::Side_of_oriented_sphere_3 Side_of_oriented_sphere_3;
    typedef typename TraitsT::Kinetic_kernel::Orientation_3 Orientation_3;
    typedef internal::Delaunay_3_edge_flip_event<This_DT3, Root_stack> Edge_flip;
    typedef typename internal::Delaunay_3_facet_flip_event<This_DT3, Root_stack> Facet_flip;

    Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object() const
    {
      return TraitsT::kinetic_kernel_object().side_of_oriented_sphere_3_object();
    }

    Orientation_3 orientation_3_object() const
    {
      return TraitsT::kinetic_kernel_object().orientation_3_object();
    }

    Base_traits(This_DT3 *t, const TraitsT &tr): TraitsT(tr), wr_(t) {}

    This_DT3* wrapper_handle() {
      return wr_;
    }
    const This_DT3* wrapper_handle() const
    {
      return wr_;
    }

    This_DT3 *wr_;
  };

 
  friend class internal::Delaunay_event_base_3<This, Root_stack>;  

  friend class internal::Delaunay_3_edge_flip_event<This, Root_stack>;

  friend class internal::Delaunay_3_facet_flip_event<This, Root_stack>;

  

  typedef internal::Delaunay_triangulation_base_3<Base_traits, Visitor> KDel;

  CGAL_KINETIC_DECLARE_LISTENERS(typename TraitsT::Simulator,
				 typename TraitsT::Active_points_3_table)

public:
  //! Initialize it.
  Delaunay_triangulation_3(TraitsT tr, Visitor v= Visitor()): kdel_(Base_traits(this, tr), v) {
    CGAL_KINETIC_INITIALIZE_LISTENERS(tr.simulator_handle(),
				      tr.active_points_3_table_handle());
  }

  //! The type of the underlying triangulation
  typedef TriangulationT Triangulation;
  //! access the underlying triangulation
  const Triangulation& triangulation() const
  {
    return kdel_.triangulation();
  }

  Visitor& visitor() {
    return kdel_.visitor();
  }

  const Visitor& visitor() const
  {
    return kdel_.visitor();
  }



  void write(std::ostream &out) const
  {
    kdel_.write(out);
  }


  //! make the structure have or not have certificates
  void set_has_certificates(bool tf) {
    kdel_.set_has_certificates(tf);
  }
  
  void audit() const
  {
    kdel_.audit();
  }

  //! true if the structure has certificates
  bool has_certificates() const
  {
    return kdel_.has_certificates();
  }

  void erase(Point_key k) {
    kdel_.delete_vertex(k);
    on_geometry_changed();
  }

  void set(Point_key k) {
    kdel_.change_vertex(k);
  }

  void insert(Point_key k) {

    kdel_.insert(k);
    /*if (kdel_.triangulation()->dimension() ==3){
      kdel_.set_has_certificates(true);
      }*/
    on_geometry_changed();
  }

  void flip(const typename KDel::Edge &edge) {
    kdel_.flip(edge);
    on_geometry_changed();
  }

  void flip(const typename KDel::Facet &flip_facet) {
    kdel_.flip(flip_facet);
    on_geometry_changed();
  }

 

  CGAL_KINETIC_LISTENER1(TRIANGULATION)

 void on_geometry_changed() {
    CGAL_KINETIC_NOTIFY(TRIANGULATION);
  }

  KDel kdel_;
};

} } //namespace CGAL::Kinetic

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif
