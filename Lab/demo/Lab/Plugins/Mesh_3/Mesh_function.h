// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <any>

namespace CGAL {
  class Image_3;
}

struct Mesh_parameters
{
  double facet_angle;
  double facet_sizing;
  double facet_min_sizing;
  double facet_approx;

  double tet_shape;
  double tet_sizing;
  double tet_min_sizing;
  double edge_sizing;
  double edge_min_sizing;
  double edge_distance;
  bool protect_features;
  bool detect_connected_components;
  int manifold;
  const CGAL::Image_3* image_3_ptr;
  const CGAL::Image_3* weights_ptr;
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

  Edge_criteria edge_criteria(double b, double minb, double d, Mesh_fnt::Domain_tag);
  Edge_criteria edge_criteria(double b, double minb, double d, Mesh_fnt::Polyhedral_domain_tag);

  void tweak_criteria(Mesh_criteria&, Mesh_fnt::Domain_tag) {}
  void tweak_criteria(Mesh_criteria&, Mesh_fnt::Polyhedral_domain_tag);
private:
  std::any object_to_destroy;
  C3t3& c3t3_;
  Domain* const domain_;
  Mesh_parameters const p_;
  std::atomic<bool> stop_;
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
  QStringList res("Mesh criteria");

  // doubles
  if(edge_sizing > 0)
    res << QString("edge max size: %1").arg(edge_sizing);
  if(edge_min_sizing > 0)
    res << QString("edge min size: %1").arg(edge_min_sizing);
  if(edge_distance > 0)
    res << QString("edge max distance: %1").arg(edge_distance);
  if(facet_angle > 0)
    res << QString("facet min angle: %1").arg(facet_angle);
  if(facet_sizing > 0)
    res << QString("facet max size: %1").arg(facet_sizing);
  if(facet_min_sizing > 0)
    res << QString("facet min size: %1").arg(facet_min_sizing);
  if(facet_approx > 0)
    res << QString("facet approx error: %1").arg(facet_approx);
  if(tet_shape > 0)
    res << QString("tet shape (radius-edge): %1").arg(tet_shape);
  if(tet_sizing > 0)
    res << QString("tet max size: %1").arg(tet_sizing);
  if(tet_min_sizing > 0)
    res << QString("tet min size: %1").arg(tet_min_sizing);

  // booleans
  res << QString("protect features: %1").arg(protect_features);
  if(image_3_ptr != nullptr)
  {
    res << QString("detect connected components: %1")
             .arg(detect_connected_components);
    res << QString("use weights: %1").arg(weights_ptr != nullptr);
  }
  res << QString("use aabb tree: %1").arg(use_sizing_field_with_aabb_tree);
  res << QString("manifold: %1").arg(manifold);

  return res;
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
, stop_()
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
    CGAL_IMAGE_IO_CASE(p_.image_3_ptr->image(),
            initialize_triangulation_from_labeled_image(c3t3_
                                                        , *domain_
                                                        , *p_.image_3_ptr
                                                        , criteria
                                                        , Word()
                                                        , p_.protect_features);
                       );
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
  namespace p = CGAL::parameters;
  // Initialization of the mesh, either with the protection of sharp
  // features, or with the initial points (or both).
  // If `detect_connected_components==true`, the initialization is
  // already done.
  CGAL::Mesh_3::internal::C3t3_initializer<
    C3t3,
    Domain,
    Mesh_criteria,
    CGAL::internal::has_Has_features<Domain>::value >()
    (c3t3_,
     *domain_,
     criteria,
     p_.protect_features,
     p::mesh_3_options(p::pointer_to_stop_atomic_boolean = &stop_,
                       p::nonlinear_growth_of_balls = true).v);
}

template < typename D_, typename Tag >
typename Mesh_function<D_,Tag>::Edge_criteria
Mesh_function<D_,Tag>::
edge_criteria(double edge_size, double minb, double edge_dist, Mesh_fnt::Domain_tag)
{
  return Edge_criteria(edge_size, minb, edge_dist);
}

#include <CGAL/Sizing_field_with_aabb_tree.h>
#include <CGAL/Mesh_3/experimental/Facet_topological_criterion_with_adjacency.h>

template < typename D_, typename Tag >
typename Mesh_function<D_,Tag>::Edge_criteria
Mesh_function<D_,Tag>::
edge_criteria(double edge_size, double minb, double edge_dist, Mesh_fnt::Polyhedral_domain_tag)
{
  if(p_.use_sizing_field_with_aabb_tree) {
    typedef typename Domain::Surface_patch_index_set Set_of_patch_ids;
    typedef CGAL::Sizing_field_with_aabb_tree<Kernel, Domain> Mesh_sizing_field; // type of sizing field for 0D and 1D features
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
      new Mesh_sizing_field(edge_size, *domain_, domain_->aabb_tree());
    // The sizing field object, as well as the `patch_ids_vector` are
    // allocated on the heap, and the following `std::any` object,
    // containing two shared pointers, is used to make the allocated
    // objects be destroyed at the destruction of the thread object, using
    // type erasure (`std::any`).
    object_to_destroy =
      std::make_pair(QSharedPointer<Mesh_sizing_field>(sizing_field_ptr),
                     QSharedPointer<Patches_ids_vector>(patches_ids_vector_p));

    std::cerr << "Note: Mesh_3 is using a sizing field based on AABB tree.\n";
    return Edge_criteria(*sizing_field_ptr, minb, edge_dist);
  } else {
    return Edge_criteria(edge_size, minb, edge_dist);
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
  Mesh_criteria criteria(edge_criteria(p_.edge_sizing,
                                       p_.edge_min_sizing,
                                       p_.edge_distance,
                                       Tag()),
                         Facet_criteria(p_.facet_angle,
                                        p_.facet_sizing,
                                        p_.facet_approx,
                                        CGAL::FACET_VERTICES_ON_SURFACE,
                                        p_.facet_min_sizing),
                         Cell_criteria(p_.tet_shape,
                                       p_.tet_sizing,
                                       p_.tet_min_sizing));

  tweak_criteria(criteria, Tag());
  initialize(criteria, Tag());

  // Build mesher and launch refinement process
  mesher_ = new Mesher(c3t3_, *domain_, criteria,
                       topology(p_.manifold) & CGAL::MANIFOLD,
                       0,
                       0,
                       &stop_);

#ifdef CGAL_MESH_3_PROFILING
  CGAL::Real_timer t;
  t.start();
#endif

#if CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  mesher_->initialize();
  while ( ! mesher_->is_algorithm_done() && ! stop_ )
  {
    mesher_->one_step();
  }
#else // not CGAL_MESH_3_MESHER_STATUS_ACTIVATED
  mesher_->refine_mesh();
#endif

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << "Full refinement time (without fix_c3t3): " << t.time() << " seconds." << std::endl;
#endif

  // Ensure c3t3 is ok (useful if process has been stop by the user)
  mesher_->fix_c3t3();
  std::cerr<<"Done."<<std::endl;
}


template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
tweak_criteria(Mesh_criteria& c, Mesh_fnt::Polyhedral_domain_tag) {
  typedef CGAL::Mesh_3::Facet_topological_criterion_with_adjacency<Tr,
       Domain, typename Facet_criteria::Visitor> New_topo_adj_crit;

  if(((int(c.facet_criteria_object().topology()) &
      CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK) != 0)
    && c.edge_criteria_object().min_length_bound() == 0)
  {
    c.add_facet_criterion(new New_topo_adj_crit(this->domain_));
  }
}


template < typename D_, typename Tag >
void
Mesh_function<D_,Tag>::
stop()
{
  stop_.store(true, std::memory_order_release);
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