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

#include <CGAL/Mesh_domain_holder_with_corners_3.h>
#include <CGAL/refine_periodic_3_mesh_3.h>

#include <CGAL/make_mesh_3.h>

#include <boost/parameter/preprocessor.hpp>

namespace CGAL {

// -----------------------------------
// make_periodic_3_mesh_3 stuff
// -----------------------------------

// Manual redirections
// boost::parameter can't handle make_periodic_3_mesh_3 return_type alone...
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template <typename C3T3, typename MD, typename MC, typename ... T>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const T& ...t)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,t...);
  return c3t3;
}
#else

template <typename C3T3, typename MD, typename MC>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC, typename Arg1>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,a1);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC, typename Arg1, typename Arg2>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,a1,a2);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                            const Arg3& a3)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,a1,a2,a3);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3, typename Arg4>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc, const Arg1& a1, const Arg2& a2,
                            const Arg3& a3, const Arg4& a4)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,a1,a2,a3,a4);
  return c3t3;
}

template <typename C3T3, typename MD, typename MC,
          typename Arg1, typename Arg2, typename Arg3,
          typename Arg4, typename Arg5>
C3T3 make_periodic_3_mesh_3(const MD& md, const MC& mc,
                            const Arg1& a1, const Arg2& a2, const Arg3& a3,
                            const Arg4& a4, const Arg5& a5)
{
  typedef CGAL::Mesh_domain_holder_with_corners_3<MD> MDH;

  C3T3 c3t3;
  MDH mdh(md);

  make_periodic_3_mesh_3_bp(c3t3,mdh,mc,a1,a2,a3,a4,a5);
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
    (exude_param, (parameters::internal::Exude_options), parameters::no_exude()) // another default parameter distinct from Mesh_3
    (perturb_param, (parameters::internal::Perturb_options), parameters::no_perturb()) // another default parameter distinct from Mesh_3
    (odt_param, (parameters::internal::Odt_options), parameters::no_odt())
    (lloyd_param, (parameters::internal::Lloyd_options), parameters::no_lloyd())
    (mesh_options_param, (parameters::internal::Mesh_3_options),
                         parameters::internal::Mesh_3_options())
    )
   )
  )
{
  make_periodic_3_mesh_3_impl(c3t3, domain, criteria,
                              exude_param, perturb_param, odt_param, lloyd_param,
                              features_param.features());
}
CGAL_PRAGMA_DIAG_POP

template<class C3T3, class MeshDomain>
void init_triangulation(C3T3& c3t3, MeshDomain& domain)
{
  typedef typename C3T3::Triangulation Tr;

  Tr& tr = c3t3.triangulation();
  tr.set_domain(domain.periodic_bounding_box());

  init_default_triangulation_vertices(c3t3, domain);

  init_domain(c3t3, domain);
}

/**
 * @brief If the triangulation of a complex contains any vertices, then
 *  we initialize them with the proper information
 */
template<class C3T3, class MD>
void init_default_triangulation_vertices(C3T3& c3t3, const MD& oracle)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename MD::Subdomain Subdomain;
  typedef typename MD::Is_in_domain Is_in_domain;

  typename Tr::Geom_traits::Construct_point_3 wp2p =
      c3t3.triangulation().geom_traits().construct_point_3_object();

  // test whether a point is in domain
  Is_in_domain is_in_domain = oracle.is_in_domain_object();
  Tr& tr = c3t3.triangulation();

  // go over the vertices
  for(Finite_vertices_iterator vertex_it = tr.finite_vertices_begin();
      vertex_it != tr.finite_vertices_end();
      ++vertex_it) {

    // test
    const Subdomain subdomain = is_in_domain(wp2p(vertex_it->point()));

    // if the point is inside
    if(subdomain) {
      // set the index that is equal to the index of the subdomain
      c3t3.set_index(vertex_it, oracle.index_from_subdomain_index(*subdomain));

      // mark the dimension, by construction it's 3
      vertex_it->set_dimension(3);
    }
  }
}

template<class C3T3, class MD>
void init_domain(C3T3& c3t3, MD& oracle)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  Tr& tr = c3t3.triangulation();

  typename Tr::Geom_traits::Construct_point_3 wp2p =
    tr.geom_traits().construct_point_3_object();

  // At this point, the triangulation contains the dummy points.
  // Mark them all as corners so they are kept throughout the refinement.
  // [ Does that really matter ? @fixme ]

  for(Finite_vertices_iterator vertex_it = tr.finite_vertices_begin();
      vertex_it != tr.finite_vertices_end(); ++vertex_it) {
    oracle.add_corner(wp2p(vertex_it->point()));
  }
}

template<class C3T3, class MeshDomain>
void projection_of_external_points_of_surface(C3T3& c3t3, const MeshDomain& domain)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef std::set<Vertex_handle> Vertex_set;

  Vertex_set vertex_set;

  find_points_to_project(c3t3, std::insert_iterator<Vertex_set>(vertex_set, vertex_set.begin()));

  std::cout << "nb of points to project" << vertex_set.size() << std::endl;
  projection_of_points(c3t3, domain, vertex_set.begin(), vertex_set.end());
}

template<class C3T3, class OutputIterator>
void find_points_to_project(C3T3& c3t3, OutputIterator vertices)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;

  for(typename C3T3::Facets_in_complex_iterator face_it = c3t3.facets_in_complex_begin();
      face_it != c3t3.facets_in_complex_end();
      ++face_it) {

    int ind = face_it->second;
    Cell_handle c = face_it->first;

    for(int i = 1; i < 4; i++) {
      Vertex_handle v = c->vertex((ind+i)&3);

      if(0 == c3t3.index(v)) {
        *vertices++ = v;
      }
    }
  }
}

template<class C3T3, class MeshDomain, class InputIterator>
void projection_of_points(C3T3& c3t3,
                          const MeshDomain& domain,
                          InputIterator vertex_begin,
                          InputIterator vertex_end)
{
  typedef typename C3T3::Triangulation::Point Point;
  typedef typename C3T3::Vertex_handle Vertex_handle;

  Mesh_3::C3T3_helpers<C3T3, MeshDomain> helper(c3t3, domain);

  for(InputIterator it = vertex_begin; it != vertex_end; ++it) {
    Point new_point = helper.project_on_surface((*it)->point(),*it);

    std::cout << "squared distance from dummy to surface: " << CGAL::squared_distance(new_point, (*it)->point()) << std::endl;

    if( new_point != Point()) {
      //Vertex_handle new_vertex = c3t3.triangulation().insert(new_point);
      Vertex_handle new_vertex = helper.update_mesh(new_point, *it);

      c3t3.set_dimension(new_vertex, 2);
    }
  }
}

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
  // Initialize the periodic triangulation
  init_triangulation(c3t3, const_cast<MeshDomain&>(domain));

  return make_mesh_3_impl(c3t3, domain, criteria,
                          exude, perturb, odt, lloyd, with_features,
                          mesh_options, manifold_options);
}

}  // end namespace CGAL

#endif // CGAL_MAKE_PERIODIC_3_MESH_3_H
