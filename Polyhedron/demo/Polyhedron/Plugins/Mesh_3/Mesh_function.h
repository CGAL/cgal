// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_DEMO_MESH_3_MESH_FUNCTION_H
#define CGAL_DEMO_MESH_3_MESH_FUNCTION_H

#define CGAL_MESH_3_MESHER_STATUS_ACTIVATED 1

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

#include <QStringList>
#include <QString>

#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/Mesh_3/initialize_triangulation_from_labeled_image.h>

#include "C3t3_type.h"
#include "Meshing_thread.h"
#include <CGAL/make_mesh_3.h> // for C3t3_initializer
#include <CGAL/use.h>

#include <boost/any.hpp>
#include <CGAL/Mesh_3/experimental/Get_facet_patch_id.h>

namespace CGAL {
  class Image_3;
}
namespace internal{
//general case for polyhedron
template<class Domain>
struct Get_facet_patch_id_selector {
  typedef CGAL::Default type;
};
//specialization for surface_mesh
template<>
struct Get_facet_patch_id_selector<Polyhedral_mesh_domain_sm> {
  typedef CGAL::Mesh_3::Get_facet_patch_id_sm<Polyhedral_mesh_domain_sm> type;
};
}//end internal
struct Mesh_parameters
{
  double facet_angle;
  double facet_sizing;
  double facet_approx;
  
  double tet_shape;
  double tet_sizing;
  double edge_sizing;
  bool protect_features;
  bool detect_connected_components;
  int manifold;
  const CGAL::Image_3* image_3_ptr;
  bool use_sizing_field_with_aabb_tree;
  
  inline QStringList log() const;
};

namespace Mesh_fnt {
  struct Domain_tag {};
  struct Polyhedral_domain_tag : public Domain_tag {};
  struct Labeled_domain_tag : public Domain_tag {};
  struct Image_domain_tag : Labeled_domain_tag {};
  struct Implicit_domain_tag : Labeled_domain_tag {};
  struct Labeled_image_domain_tag : Image_domain_tag {};
  struct Gray_image_domain_tag : Image_domain_tag {};
} // end Mesh_fnt

template < typename Domain_, typename Domain_tag >
class Mesh_function
  : public Mesh_function_interface
{
  typedef Domain_ Domain;
  
public:
  Mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p);
  
  ~Mesh_function();
  
  // Launch
  virtual void launch();
  
  // Stop
  virtual void stop();
  
  // Logs
  virtual QStringList parameters_log() const;
  virtual QString status(double time_period) const;

private:
  typedef typename Domain::Point_3                  Point_3;
  typedef typename Domain::Index                    Index;
  typedef std::vector<std::pair<Point_3, Index> >   Initial_points_vector;
  typedef typename Initial_points_vector::iterator  Ipv_iterator;
  typedef C3t3::Vertex_handle                       Vertex_handle;
  
  typedef C3t3::Triangulation                       Tr;
  typedef CGAL::Mesh_criteria_3<Tr>                 Mesh_criteria;
  typedef Mesh_criteria::Edge_criteria              Edge_criteria;
  typedef Mesh_criteria::Facet_criteria             Facet_criteria;
  typedef Mesh_criteria::Cell_criteria              Cell_criteria;
  
  typedef CGAL::Mesh_3::Mesher_3<C3t3, Mesh_criteria, Domain>   Mesher;

  void initialize(const Mesh_criteria& criteria, Mesh_fnt::Domain_tag);
  void initialize(const Mesh_criteria& criteria, Mesh_fnt::Labeled_image_domain_tag);

  Edge_criteria edge_criteria(double b, Mesh_fnt::Domain_tag);
  Edge_criteria edge_criteria(double b, Mesh_fnt::Polyhedral_domain_tag);

  void tweak_criteria(Mesh_criteria&, Mesh_fnt::Domain_tag) {}
  void tweak_criteria(Mesh_criteria&, Mesh_fnt::Polyhedral_domain_tag);
private:
  boost::any object_to_destroy;
  C3t3& c3t3_;
  Domain* domain_;
  Mesh_parameters p_;
  bool continue_;
  Mesher* mesher_;
#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  mutable typename Mesher::Mesher_status last_report_;
#endif
};



// -----------------------------------
// Class Mesh_parameters
// -----------------------------------
inline
QStringList
Mesh_parameters::
log() const
{
  return QStringList()
  << QString("edge max size: %1").arg(edge_sizing)
  << QString("facet min angle: %1").arg(facet_angle)
  << QString("facet max size: %1").arg(facet_sizing)
  << QString("facet approx error: %1").arg(facet_approx)
  << QString("tet shape (radius-edge): %1").arg(tet_shape)
  << QString("tet max size: %1").arg(tet_sizing)
  << QString("detect connected components: %1")
    .arg(detect_connected_components)
  << QString("protect features: %1").arg(protect_features);
}


// -----------------------------------
// Class Mesh_function
// -----------------------------------
template < typename D_, typename Tag >
Mesh_function<D_,Tag>::
Mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p)
: c3t3_(c3t3)
, domain_(domain)
, p_(p)
, continue_(true)
, mesher_(NULL)
#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
, last_report_(0,0,0)
#endif
{
#ifdef CGAL_CONCURRENT_MESH_3
  Concurrent_mesher_config::load_config_file(CONFIG_FILENAME, false);
#endif
}


template < typename D_, typename Tag >
Mesh_function<D_,Tag>::
~Mesh_function()
{
  delete domain_;
  delete mesher_;
}


CGAL::Mesh_facet_topology topology(int manifold) {
  switch(manifold) {
  case 1: return static_cast<CGAL::Mesh_facet_topology>
      (CGAL::MANIFOLD |
       CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH);
  case 2: return CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK;
  case 3: return static_cast<CGAL::Mesh_facet_topology>
      (CGAL::MANIFOLD |
       CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK);
  default: return CGAL::FACET_VERTICES_ON_SURFACE;
  }
}

template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
initialize(const Mesh_criteria& criteria, Mesh_fnt::Labeled_image_domain_tag)
// for a labeled image
{
  if(p_.detect_connected_components) {
    initialize_triangulation_from_labeled_image(c3t3_
                                                , *domain_
                                                , *p_.image_3_ptr
                                                , criteria
                                                , typename D_::Image_word_type()
                                                , p_.protect_features);
  } else {
    initialize(criteria, Mesh_fnt::Domain_tag());
  }
}

template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
initialize(const Mesh_criteria& criteria, Mesh_fnt::Domain_tag)
// for the other domain types
{
  // Initialization of the mesh, either with the protection of sharp
  // features, or with the initial points (or both).
  // If `detect_connected_components==true`, the initialization is
  // already done.
  CGAL::internal::Mesh_3::C3t3_initializer<
    C3t3,
    Domain,
    Mesh_criteria,
    CGAL::internal::Mesh_3::has_Has_features<Domain>::value >()
    (c3t3_,
     *domain_,
     criteria,
     p_.protect_features,
     p_.use_sizing_field_with_aabb_tree);
}

template < typename D_, typename Tag >
typename Mesh_function<D_,Tag>::Edge_criteria
Mesh_function<D_,Tag>::
edge_criteria(double b, Mesh_fnt::Domain_tag)
{
  return Edge_criteria(b);
}

#include <CGAL/Mesh_3/experimental/Sizing_field_with_aabb_tree.h>
#include <CGAL/Mesh_3/experimental/Facet_topological_criterion_with_adjacency.h>

template < typename D_, typename Tag >
typename Mesh_function<D_,Tag>::Edge_criteria
Mesh_function<D_,Tag>::
edge_criteria(double edge_size, Mesh_fnt::Polyhedral_domain_tag)
{
  if(p_.use_sizing_field_with_aabb_tree) {
    typedef typename Domain::Surface_patch_index_set Set_of_patch_ids;
    typedef Sizing_field_with_aabb_tree
      <
      Kernel
      , Domain
      , typename Domain::AABB_tree
      , CGAL::Default
      , typename internal::Get_facet_patch_id_selector<Domain>::type
      > Mesh_sizing_field; // type of sizing field for 0D and 1D features
    typedef std::vector<Set_of_patch_ids> Patches_ids_vector;
    typedef typename Domain::Curve_index Curve_index;
    const Curve_index max_index = domain_->maximal_curve_index();
    Patches_ids_vector* patches_ids_vector_p =
      new Patches_ids_vector(max_index+1);
    for(Curve_index curve_id = 1; curve_id <= max_index; ++curve_id)
    {
      (*patches_ids_vector_p)[curve_id] = domain_->get_incidences(curve_id);
    }
    Mesh_sizing_field* sizing_field_ptr =
      new Mesh_sizing_field(edge_size,
                            domain_->aabb_tree(),
                            *domain_);
    // The sizing field object, as well as the `patch_ids_vector` are
    // allocated on the heap, and the following `boost::any` object,
    // containing two shared pointers, is used to make the allocated
    // objects be destroyed at the destruction of the thread object, using
    // type erasure (`boost::any`).
    object_to_destroy =
      std::make_pair(QSharedPointer<Mesh_sizing_field>(sizing_field_ptr),
                     QSharedPointer<Patches_ids_vector>(patches_ids_vector_p));

    std::cerr << "USE SIZING FIELD!\n";
    return Edge_criteria(*sizing_field_ptr);
  } else {
    return Edge_criteria(edge_size);
  }
}

template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
launch()
{
#ifdef CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
  CGAL::default_random = CGAL::Random(0);
#endif

  // Create mesh criteria
  Mesh_criteria criteria(edge_criteria(p_.edge_sizing, Tag()),
                         Facet_criteria(p_.facet_angle,
                                        p_.facet_sizing,
                                        p_.facet_approx),
                         Cell_criteria(p_.tet_shape,
                                       p_.tet_sizing));

  tweak_criteria(criteria, Tag());
  initialize(criteria, Tag());

  // Build mesher and launch refinement process
  mesher_ = new Mesher(c3t3_, *domain_, criteria,
                       topology(p_.manifold) & CGAL::MANIFOLD);

#ifdef CGAL_MESH_3_PROFILING
  CGAL::Real_timer t;
  t.start();
#endif

#if CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  mesher_->initialize();
  while ( ! mesher_->is_algorithm_done() && continue_ )
  {
    mesher_->one_step();
  }
#else // not CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  mesher_->refine_mesh();
#endif

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Full refinement time (without fix_c3t3): " << t.time() << " seconds." << std::endl;
#endif

  // Ensure c3t3 is ok (usefull if process has been stop by the user)
  mesher_->fix_c3t3();
}


template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
tweak_criteria(Mesh_criteria& c, Mesh_fnt::Polyhedral_domain_tag) {
  typedef CGAL::Mesh_3::Facet_topological_criterion_with_adjacency<Tr,
       Domain, typename Facet_criteria::Visitor> New_topo_adj_crit;

  if((int(c.facet_criteria().topology()) &
      CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK) != 0)
  {
    c.add_facet_criterion(new New_topo_adj_crit(this->domain_));
  }
}


template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
stop()
{
  continue_ = false;
}


template < typename D_, typename Tag >
QStringList
Mesh_function<D_,Tag>::
parameters_log() const
{
  return p_.log();
}


template < typename D_, typename Tag >
QString
Mesh_function<D_,Tag>::
status(double time_period) const
{
  QString result;

  CGAL_USE(time_period); // to avoid a warning when the macro
                         // CGAL_MESH_3_MESHER_STATUS_ACTIVATED is not
                         // defined
#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  // If mesher_ is not yet created, it means that either launch() has not
  // been called or that initial points have not been founded
  if ( NULL == mesher_ ) /// @TODO: there is a race-condition, here.
  {
    typename Mesher::Mesher_status s(c3t3_.triangulation().number_of_vertices(),
                                     0,
                                     0);
    result = QString("Initialization in progress...\n\n"
                     "Vertices: %1 \n"
                     "Vertices inserted last %2s: %3 \n")
      .arg(s.vertices)
      .arg(time_period)
      .arg(s.vertices - last_report_.vertices);
    last_report_ = s;
  } else {
    // Get status and return a string corresponding to it
    typename Mesher::Mesher_status s = mesher_->status();
  
    result = QString("Vertices: %1 \n"
                     "Vertices inserted last %2s: %3 \n\n"
                     "Bad facets: %4 \n"
                     "Bad cells: %5")
      .arg(s.vertices)
      .arg(time_period)
      .arg(s.vertices - last_report_.vertices)
      .arg(s.facet_queue)
      .arg(s.cells_queue);
    last_report_ = s;
  }
  
#endif
  return result;
}

#endif // CGAL_DEMO_MESH_3_MESH_FUNCTION_H
