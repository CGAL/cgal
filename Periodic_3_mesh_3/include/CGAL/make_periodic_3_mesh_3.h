// Copyright (c) 2009, 2017 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb, Mikhail Bogdanov
//
//******************************************************************************
// File Description : make_periodic_3_mesh_3 function definition.
//******************************************************************************

#ifndef CGAL_MAKE_PERIODIC_3_MESH_3_H
#define CGAL_MAKE_PERIODIC_3_MESH_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>
#include <CGAL/Mesh_3/global_parameters.h>

#include <CGAL/refine_periodic_3_mesh_3.h>

#include <CGAL/assertions.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <boost/parameter/preprocessor.hpp>

#include <iterator>

namespace CGAL {

namespace internal {
namespace Periodic_3_mesh_3 {

template<typename C3T3>
void give_dummy_points_artificial_index(C3T3& c3t3)
{
  CGAL_precondition(c3t3.triangulation().is_1_cover());

  typedef typename C3T3::Triangulation::Vertex_iterator       Vertex_iterator;

  for(Vertex_iterator vit = c3t3.triangulation().vertices_begin();
                      vit != c3t3.triangulation().vertices_end(); ++vit)
  {
    std::cout << "dummy point: " << vit->point() << std::endl;
    c3t3.set_index(vit, 0);
  }
}

// C3t3_initializer: initialize c3t3
template < typename C3T3,
           typename MeshDomain,
           typename MeshCriteria,
           bool MeshDomainHasHasFeatures,
           typename HasFeatures = int>
class C3t3_initializer
  : public CGAL::internal::Mesh_3::C3t3_initializer<
      C3T3, MeshDomain, MeshCriteria, MeshDomainHasHasFeatures, HasFeatures>
{
  typedef CGAL::internal::Mesh_3::C3t3_initializer<
            C3T3, MeshDomain, MeshCriteria,
            MeshDomainHasHasFeatures, HasFeatures>              Base;

public:
  void operator()(C3T3& c3t3,
                  const MeshDomain& domain,
                  const MeshCriteria& criteria,
                  bool with_features,
                  bool nonlinear = false,
                  const int nb_initial_points = -1)
  {
    c3t3.triangulation().set_domain(domain.periodic_bounding_box());
    c3t3.triangulation().insert_dummy_points();

    give_dummy_points_artificial_index(c3t3);

    // Call the basic initialization from c3t3, which handles features and
    // adds a bunch of points on the surface
    Base::operator()(c3t3, domain, criteria, with_features,
                     nonlinear, nb_initial_points);
  }
};

} // namespace Periodic_3_mesh_3
} // namespace internal

// -----------------------------------
// make_periodic_3_mesh_3 stuff
// -----------------------------------

// Manual redirections
// boost::parameter can't handle make_periodic_3_mesh_3 return_type alone...
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template <typename C3T3, typename MD, typename MC, typename ... T>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const T& ...t)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,t...);
  return c3t3;
}
#else

template <typename C3T3, typename MD, typename MC>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC, typename Arg1>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,a1);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC, typename Arg1, typename Arg2>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,a1,a2);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                            const Arg3& a3)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,a1,a2,a3);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3, typename Arg4>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                            const Arg3& a3, const Arg4& a4)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,a1,a2,a3,a4);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3,
          typename Arg4, typename Arg5>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc,
                            const Arg1& a1, const Arg2& a2, const Arg3& a3,
                            const Arg4& a4, const Arg5& a5)
{
  C3T3 c3t3;
  make_periodic_3_mesh_3_bp(c3t3,md,mc,a1,a2,a3,a4,a5);
  return c3t3;
}
#endif

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/Mesh_3/config.h>
CGAL_MESH_3_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

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
  internal::Periodic_3_mesh_3::C3t3_initializer<
    C3T3, MeshDomain, MeshCriteria,
    internal::Mesh_3::has_Has_features<MeshDomain>::value>()(c3t3,
                                                             domain,
                                                             criteria,
                                                             with_features,
                                                             mesh_options.nonlinear_growth_of_balls,
                                                             mesh_options.number_of_initial_points);

  // Build mesher and launch refinement process
  refine_periodic_3_mesh_3(c3t3, domain, criteria,
                           exude, perturb, odt, lloyd,
                           parameters::no_reset_c3t3(), // do not reset c3t3 as we just created it
                           mesh_options, manifold_options);
}

}  // end namespace CGAL

#endif // CGAL_MAKE_PERIODIC_3_MESH_3_H
