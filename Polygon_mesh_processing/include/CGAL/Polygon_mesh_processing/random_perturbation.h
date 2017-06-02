// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H
#define CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Random.h>

#include <boost/foreach.hpp>

#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

  template<typename GT, typename PM, typename VCMap, typename VPMap>
  void random_perturbation_impl(PM& pmesh,
                                const double& max_size,
                                VCMap vcmap,
                                VPMap vpmap,
                                const int seed)
  {
    typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
    typedef typename GT::Point_3 Point_3;
    typedef typename GT::FT      FT;

    CGAL::Random rng(seed);

    BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
    {
      if (!get(vcmap, v) && !is_border(v, pmesh))
      {
        const Point_3& p = get(vpmap, v);
        const Point_3 np(p.x() + FT(rng.get_double(0., max_size)),
                         p.y() + FT(rng.get_double(0., max_size)),
                         p.z() + FT(rng.get_double(0., max_size)));
        put(vpmap, v, np);
      }
    }
  }

} //end namespace internal

/*!
* \ingroup PMP_meshing_grp
* @brief perturbs the locations of vertices of a polygon mesh.

* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during remeshing
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
*/
template<typename PolygonMesh, typename NamedParameters>
void random_perturbation(PolygonMesh& pmesh
                       , const double& perturbation_max_size
                       , const NamedParameters& np)
{
  typedef PolygonMesh PM;
  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
  std::cout << std::endl;
  CGAL::Timer t;
  std::cout << "Random perturbation (max size = "<< perturbation_max_size<<")...";
  std::cout.flush();
  t.start();
#endif

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pmesh));

  typedef typename boost::lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      internal::No_constraint_pmap<vertex_descriptor>//default
    > ::type VCMap;
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             internal::No_constraint_pmap<vertex_descriptor>());

  internal::random_perturbation_impl<GT>(pmesh,
                                         perturbation_max_size,
                                         vcmap,
                                         vpmap,
                                         0/*seed*/);

#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
  t.stop();
  std::cout << "Perturbation done (;
  std::cout << t.time() << " sec )." << std::endl;
#endif
}

template<typename PolygonMesh>
void random_perturbation(PolygonMesh& pmesh
                       , const double& perturbation_max_size)
{
  random_perturbation(pmesh,
                      perturbation_max_size,
                      parameters::all_default());
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H

