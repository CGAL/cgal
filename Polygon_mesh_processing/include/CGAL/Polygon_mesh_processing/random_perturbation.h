// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H
#define CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Random.h>


#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

  template<typename GT, typename RNG>
  typename GT::Vector_3 construct_random_vector_3(const double& max_size,
                                                  RNG& rng,
                                                  const GT& gt)
  {
    typedef typename GT::FT FT;
    typename GT::Construct_vector_3 vec = gt.construct_vector_3_object();

    return vec(FT(rng.get_double(-max_size, max_size)),
               FT(rng.get_double(-max_size, max_size)),
               FT(rng.get_double(-max_size, max_size)));
  }

  template<typename GT, typename VertexRange,
           typename PM, typename VCMap, typename VPMap, typename RNG>
  void random_perturbation_impl(VertexRange vrange,
                                PM& tmesh,
                                const double& max_size,
                                VCMap vcmap,
                                VPMap vpmap,
                                bool do_project,
                                RNG& rng,
                                const GT& gt)
  {
    typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
    typedef typename GT::Point_3    Point_3;

    typedef CGAL::AABB_face_graph_triangle_primitive<PM> Primitive;
    typedef CGAL::AABB_traits<GT, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    Tree tree;
    if(do_project)
    {
      tree.rebuild(faces(tmesh).first, faces(tmesh).second, tmesh);
    }
    typename GT::Construct_translated_point_3 translate
      = gt.construct_translated_point_3_object();

    for(vertex_descriptor v : vrange)
    {
      if (!get(vcmap, v) && !is_border(v, tmesh))
      {
        const Point_3& p = get(vpmap, v);
        const Point_3 np = translate(p, construct_random_vector_3<GT>(max_size, rng, gt));

        if (do_project)
          put(vpmap, v, tree.closest_point(np)); //project on input surface
        else
          put(vpmap, v, np);
      }
    }
  }

} //end namespace internal

/*!
* \ingroup PMP_meshing_grp
* @brief randomly perturbs the locations of vertices of a triangulated surface mesh.
* By default, the vertices are re-projected onto the input surface after perturbation.
* Note that no geometric checks are done after the perturbation
* (face orientation might become incorrect and self-intersections might be introduced).
*
* @tparam VertexRange model of `Range`, holding
*         vertices of type `boost::graph_traits<TriangleMesh>::%vertex_descriptor`.
*         Its iterator type is `ForwardIterator`.
* @tparam TriangleMesh model of `MutableFaceGraph`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param vertices the range of vertices to be perturbed
* @param tmesh the triangulated surface mesh
* @param perturbation_max_size the maximal length of moves that can be applied to
*        vertices of `tmesh`.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `tmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertex is constrained}
*     \cgalParamExtra{A constrained vertex cannot be modified at all during perturbation}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{do_project}
*     \cgalParamDescription{indicates whether vertices are reprojected on the input surface
*                           after their coordinates random perturbation}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{random_seed}
*     \cgalParamDescription{a value to seed the random number generator, and make the perturbation deterministic}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`unsigned int(-1)`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template<typename VertexRange, typename TriangleMesh, typename NamedParameters>
void random_perturbation(VertexRange vertices
                       , TriangleMesh& tmesh
                       , const double& perturbation_max_size
                       , const NamedParameters& np)
{
  typedef TriangleMesh PM;
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
  std::cout << std::endl;
  CGAL::Timer t;
  std::cout << "Random perturbation (max size = "<< perturbation_max_size<<")...";
  std::cout.flush();
  t.start();
#endif

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_property_map(vertex_point, tmesh));

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Static_boolean_property_map<vertex_descriptor, false> // default
    > ::type VCMap;
  VCMap vcmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                 Static_boolean_property_map<vertex_descriptor, false>());

  unsigned int seed = choose_parameter(get_parameter(np, internal_np::random_seed), -1);
  bool do_project = choose_parameter(get_parameter(np, internal_np::do_project), true);

  CGAL::Random rng = (seed == unsigned(-1)) ? CGAL::Random() : CGAL::Random(seed);

  internal::random_perturbation_impl(vertices,
          tmesh,
          perturbation_max_size,
          vcmap,
          vpmap,
          do_project,
          rng,
          gt);

#ifdef CGAL_PMP_RANDOM_PERTURBATION_VERBOSE
  t.stop();
  std::cout << "Perturbation done (";
  std::cout << t.time() << " sec )." << std::endl;
#endif
}

/*!
* \ingroup PMP_meshing_grp
* @brief same as above, but all non-border vertices of `tmesh` are perturbed.
*/
template<typename TriangleMesh, typename NamedParameters>
void random_perturbation(TriangleMesh& tmesh
                       , const double& perturbation_max_size
                       , const NamedParameters& np)
{
  random_perturbation(vertices(tmesh), tmesh, perturbation_max_size, np);
}

template<typename VertexRange, typename TriangleMesh>
void random_perturbation(VertexRange vertices
                       , TriangleMesh& tmesh
                       , const double& perturbation_max_size)
{
  random_perturbation(vertices, tmesh, perturbation_max_size,
                      parameters::all_default());
}

template<typename TriangleMesh>
void random_perturbation(TriangleMesh& tmesh
                       , const double& perturbation_max_size)
{
  random_perturbation(tmesh,
                      perturbation_max_size,
                      parameters::all_default());
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_RANDOM_PERTURBATION_H
