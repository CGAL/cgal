// Copyright (c) 2009, 2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stéphane Tayeb,
//                 Mikhail Bogdanov,
//                 Mael Rouxel-Labbé
//
//******************************************************************************
// File Description : make_periodic_3_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H
#define CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/Periodic_3_mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/refine_periodic_3_mesh_3.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/parameter.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <boost/parameter/preprocessor.hpp>

namespace CGAL {
namespace Periodic_3_mesh_3 {
namespace internal {

template<typename C3T3>
void mark_dummy_points(C3T3& c3t3)
{
  CGAL_precondition(c3t3.triangulation().is_1_cover());

  typedef typename C3T3::Triangulation::Vertex_iterator       Vertex_iterator;

  for(Vertex_iterator vit = c3t3.triangulation().vertices_begin();
                      vit != c3t3.triangulation().vertices_end(); ++vit)
  {
    c3t3.set_index(vit, 0);
  }
}

template <typename C3T3, typename MeshDomain, typename MeshCriteria>
void init_c3t3_with_features(C3T3& c3t3,
                             const MeshDomain& domain,
                             const MeshCriteria& criteria,
                             bool nonlinear = false)
{
  typedef typename MeshCriteria::Edge_criteria                                Edge_criteria;
  typedef Mesh_3::internal::Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;

  CGAL::Periodic_3_mesh_3::Protect_edges_sizing_field<C3T3, MeshDomain, Sizing_field>
    protect_edges(c3t3, domain, Sizing_field(criteria.edge_criteria_object()));
  protect_edges.set_nonlinear_growth_of_balls(nonlinear);

  protect_edges(true);
}

// C3t3_initializer: initialize c3t3
template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          bool MeshDomainHasHasFeatures,
          typename HasFeatures = int>
struct C3t3_initializer_base
  : public CGAL::Mesh_3::internal::C3t3_initializer<
      C3T3, MeshDomain, MeshCriteria, MeshDomainHasHasFeatures, HasFeatures>
{
  typedef CGAL::Mesh_3::internal::C3t3_initializer<
            C3T3, MeshDomain, MeshCriteria,
            MeshDomainHasHasFeatures, HasFeatures>              Base;

  void operator()(C3T3& c3t3,
                  const MeshDomain& domain,
                  const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    c3t3.triangulation().set_domain(domain.bounding_box());
    c3t3.triangulation().insert_dummy_points();
    mark_dummy_points(c3t3);

    // Call the basic initialization from c3t3, which handles features and
    // adds a bunch of points on the surface
    Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          bool MeshDomainHasHasFeatures,
          typename HasFeatures = int>
struct C3t3_initializer
  : public C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, MeshDomainHasHasFeatures, HasFeatures>
{
  typedef C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria,
                                MeshDomainHasHasFeatures, HasFeatures> Base;

  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

// Specialization when the mesh domain has 'Has_features'
template <typename C3T3,
          typename MeshDomain,
          typename MeshCriteria,
          typename HasFeatures>
struct C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, HasFeatures>
{
  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, typename MeshDomain::Has_features>()
        (c3t3, domain, criteria, with_features, mesh_options);
  }
};

// Specialization when the mesh domain has 'Has_features' and it's set to CGAL::Tag_true
template < typename C3T3,
           typename MeshDomain,
           typename MeshCriteria>
struct C3t3_initializer<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true>
  : public C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true>
{
  typedef C3t3_initializer_base<C3T3, MeshDomain, MeshCriteria, true, CGAL::Tag_true> Base;

  virtual ~C3t3_initializer() { }

  // this override will be used when initialize_features() is called, in make_mesh_3.h
  virtual void
  initialize_features(C3T3& c3t3,
                      const MeshDomain& domain,
                      const MeshCriteria& criteria,
                      const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Periodic_3_mesh_3::internal::init_c3t3_with_features
      (c3t3, domain, criteria, mesh_options.nonlinear_growth_of_balls);
  }

  void operator()(C3T3& c3t3, const MeshDomain& domain, const MeshCriteria& criteria,
                  bool with_features,
                  const parameters::internal::Mesh_3_options& mesh_options)
  {
    return Base::operator()(c3t3, domain, criteria, with_features, mesh_options);
  }
};

} // namespace internal
} // namespace Periodic_3_mesh_3

// -----------------------------------
// make_periodic_3_mesh_3 stuff
// -----------------------------------

// Manual redirections
// boost::parameter can't handle make_periodic_3_mesh_3 return_type alone...
template <typename C3T3, typename MD, typename MC, typename ... T>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const T& ...t)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,t...);
  return c3t3;
}

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/boost/parameter.h>
CGAL_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

BOOST_PARAMETER_FUNCTION(
  (void),
  make_periodic_3_mesh_3_bp,
  parameters::tag,
  (required (in_out(c3t3),*) (domain,*) (criteria,*) ) // nondeduced
  (deduced
   (optional
    (features_param, (parameters::internal::Features_options), parameters::features(domain))
    (exude_param, (parameters::internal::Exude_options), parameters::exude())
    (perturb_param, (parameters::internal::Perturb_options), parameters::perturb())
    (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
    (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
    (mesh_options_param, (parameters::internal::Mesh_3_options),
                         parameters::internal::Mesh_3_options())
    (manifold_options_param, (parameters::internal::Manifold_options),
                             parameters::internal::Manifold_options())
    )
   )
  )
{
  make_periodic_3_mesh_3_impl(c3t3, domain, criteria,
                              exude_param, perturb_param, odt_param, lloyd_param,
                              features_param.features(), mesh_options_param,
                              manifold_options_param);
}
CGAL_PRAGMA_DIAG_POP

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

/**
 * @brief This function meshes the domain defined by mesh_traits
 * (respecting criteria), and outputs the mesh to c3t3
 *
 * @param domain the domain to be discretized
 * @param criteria the criteria
 * @param exude if it is set to \c true, an exudation step will be done at
 *   the end of the Delaunay refinement process
 *
 * @return The mesh as a C3T3 object
 */
template<class C3T3, class MeshDomain, class MeshCriteria>
void make_periodic_3_mesh_3_impl(C3T3& c3t3,
                                 const MeshDomain& domain,
                                 const MeshCriteria& criteria,
                                 const parameters::internal::Exude_options& exude,
                                 const parameters::internal::Perturb_options& perturb,
                                 const parameters::internal::Odt_options& odt,
                                 const parameters::internal::Lloyd_options& lloyd,
                                 const bool with_features,
                                 const parameters::internal::Mesh_3_options&
                                   mesh_options = parameters::internal::Mesh_3_options(),
                                 const parameters::internal::Manifold_options&
                                   manifold_options = parameters::internal::Manifold_options())
{
  // Initialize c3t3
  Periodic_3_mesh_3::internal::C3t3_initializer<
    C3T3, MeshDomain, MeshCriteria,
    Mesh_3::internal::has_Has_features<MeshDomain>::value>()(c3t3,
                                                             domain,
                                                             criteria,
                                                             with_features,
                                                             mesh_options);

  // Build mesher and launch refinement process
  refine_periodic_3_mesh_3(c3t3, domain, criteria,
                           exude, perturb, odt, lloyd,
                           parameters::no_reset_c3t3(), // do not reset c3t3 as we just created it
                           mesh_options, manifold_options);
}

} // end namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_MAKE_PERIODIC_3_MESH_3_H
