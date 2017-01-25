#ifndef CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
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


#define CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H

#include <CGAL/license/Surface_mesh_segmentation.h>


/**
 * @file mesh_segmentation.h
 * @brief The API which contains free template functions for SDF computation and mesh segmentation.
 */
#include <CGAL/internal/Surface_mesh_segmentation/Surface_mesh_segmentation.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/config.hpp>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL
{


/// @cond SKIP_IN_MANUAL
template <bool Fast_sdf_calculation_mode, class TriangleMesh,
         class SDFPropertyMap,
         class PointPropertyMap
#ifdef DOXYGEN_RUNNING
         = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type
#endif
         , class GeomTraits
#ifdef DOXYGEN_RUNNING
         = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
         >
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true,
            PointPropertyMap ppmap = PointPropertyMap(),
            GeomTraits traits = GeomTraits())
{
  typedef PointPropertyMap VPMap;
  internal::Surface_mesh_segmentation<TriangleMesh, GeomTraits, VPMap, Fast_sdf_calculation_mode>
    algorithm(triangle_mesh, traits, ppmap);
  return algorithm.calculate_sdf_values(cone_angle, number_of_rays,
                                        sdf_values_map, postprocess);
}
/// @endcond

/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the Shape Diameter Function over a surface mesh.
 *
 * This function implements the Shape Diameter Function (SDF) as described in \cgalCite{shapira2008consistent}.
 * It is possible to compute raw SDF values (without post-processing). In such a case,
 * -1 is used to indicate when no SDF value could be computed for a facet.
 *
 * @pre `is_triangle_mesh(triangle_mesh)`
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam SDFPropertyMap  a `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `double` as value type
 * @tparam GeomTraits a model of `SegmentationGeomTraits`
 * @tparam PointPropertyMap a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and `GeomTraits::Point_3` as value type.
 *
 * @param triangle_mesh surface mesh on which SDF values are computed
 * @param[out] sdf_values_map the SDF value of each facet
 * @param cone_angle opening angle in radians for the cone of each facet
 * @param number_of_rays number of rays picked in the cone of each facet. In our experiments, we observe that increasing the number of rays beyond the default has little effect on the quality of the segmentation result
 * @param postprocess if `true`, `CGAL::sdf_values_postprocessing()` is called on raw SDF value computed.
 * @param traits traits class
 * @param ppmap point property map. An overload is provided with `get(boost::vertex_point,triangle_mesh)` as default.
 *
 * @return minimum and maximum raw SDF values if  `postprocess` is `true`, otherwise minimum and maximum SDF values (before linear normalization)
 */
template <class TriangleMesh, class SDFPropertyMap, class PointPropertyMap
#ifdef DOXYGEN_RUNNING
         = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type
#endif
, class GeomTraits
#ifdef DOXYGEN_RUNNING
=  typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
>
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true,
            PointPropertyMap ppmap = PointPropertyMap(),
            GeomTraits traits = GeomTraits())
{
  return sdf_values<true, TriangleMesh, SDFPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, cone_angle, number_of_rays, postprocess, ppmap, traits);
}


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function post-processing raw SDF values computed per facet.
 *
 * Post-processing steps applied :
 *   - Facets with -1 SDF values are assigned the average SDF value of their edge-adjacent neighbors.
 *     If there is still a facet having -1 SDF value, the minimum valid SDF value assigned to it. Note that this step is not inherited from the paper.
 *     The main reason for not assigning 0 to facets with no SDF values (i.e. -1) is that it can obstruct log-normalization process which takes place at the beginning of `CGAL::segmentation_from_sdf_values()`.
 *   - SDF values are smoothed with bilateral filtering.
 *   - SDF values are linearly normalized between [0,1].
 *
 * See the section \ref Surface_mesh_segmentationPostprocessing for more details.
 *
 * @pre `is_triangle_mesh(triangle_mesh)`
 * @pre Raw values should be greater or equal to 0. -1 indicates when no value could be computed
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam SDFPropertyMap  a `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `double` as value type
 *
 * @param triangle_mesh surface mesh on which SDF values are computed
 * @param[in, out] sdf_values_map the SDF value of each facet
 *
 * @return minimum and maximum SDF values before linear normalization
 */
template<class TriangleMesh, class SDFPropertyMap>
std::pair<double, double>
sdf_values_postprocessing(const TriangleMesh& triangle_mesh,
                          SDFPropertyMap sdf_values_map)
{
  CGAL_precondition(is_triangle_mesh(triangle_mesh));
  return internal::Postprocess_sdf_values<TriangleMesh>().postprocess(triangle_mesh,
         sdf_values_map);
}


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh given an SDF value per facet.
 *
 * This function fills a property map which associates a segment-id (in [0, number of segments -1])
 * or a cluster-id (in [0, `number_of_clusters` -1]) to each facet.
 * A segment is a set of connected facets which are placed under the same cluster (see \cgalFigureRef{Cluster_vs_segment}).
 *
 * \note Log-normalization is applied on `sdf_values_map` before segmentation.
 *       As described in the original paper \cgalCite{shapira2008consistent},
 *       this normalization is done to preserve thin parts of the mesh
 *       by increasing the distance between smaller SDF values and reducing
 *       it between larger ones.
 * \note There is no direct relation between the parameter `number_of_clusters`
 * and the final number of segments after segmentation. However, setting a large number of clusters will result in a detailed segmentation of the mesh with a large number of segments.
 *
 * @pre `is_triangle_mesh(triangle_mesh)`
 * @pre `number_of_clusters > 0`
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam SDFPropertyMap  a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `double` as value type
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type
 * @tparam GeomTraits a model of `SegmentationGeomTraits`
 * @tparam PointPropertyMap a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and `GeomTraits::Point_3` as value type.
 *
 * @param triangle_mesh surface mesh corresponding to the SDF values
 * @param sdf_values_map the SDF value of each facet between [0-1]
 * @param[out] segment_ids the segment or cluster id of each facet
 * @param number_of_clusters number of clusters for the soft clustering
 * @param smoothing_lambda factor which indicates the importance of the surface features for the energy minimization. It is recommended to choose a value in the interval [0,1]. See the section \ref Surface_mesh_segmentationGraphCut for more details.
 * @param output_cluster_ids if `false` fill `segment_ids` with segment-ids, and with cluster-ids otherwise (see \cgalFigureRef{Cluster_vs_segment})
 * @param traits traits class
 * @param ppmap point property map. An overload is provided with `get(boost::vertex_point,triangle_mesh)` as default.
 *
 * @return number of segments if `output_cluster_ids` is set to `false` and `number_of_clusters` otherwise
 */
template <class TriangleMesh, class SDFPropertyMap, class SegmentPropertyMap,
          class PointPropertyMap
#ifdef DOXYGEN_RUNNING
         = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type
#endif
         , class GeomTraits
#ifdef DOXYGEN_RUNNING
         = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
         >
std::size_t
segmentation_from_sdf_values( const TriangleMesh& triangle_mesh,
                              SDFPropertyMap sdf_values_map,
                              SegmentPropertyMap segment_ids,
                              std::size_t number_of_clusters = 5,
                              double smoothing_lambda = 0.26,
                              bool output_cluster_ids = false,
                              PointPropertyMap ppmap=PointPropertyMap(),
                              GeomTraits traits=GeomTraits())
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VPMap;
  internal::Surface_mesh_segmentation<TriangleMesh, GeomTraits, VPMap> algorithm(triangle_mesh, traits, ppmap);
  return algorithm.partition(number_of_clusters, smoothing_lambda, sdf_values_map,
                             segment_ids, !output_cluster_ids);
}

///\cond SKIP_IN_MANUAL
template < bool Fast_sdf_calculation_mode, class TriangleMesh,
         class SegmentPropertyMap, class PointPropertyMap
#ifdef DOXYGEN_RUNNING
         = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type
#endif
        , class GeomTraits
#ifdef DOXYGEN_RUNNING
         = typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
        >
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false,
                            PointPropertyMap ppmap=PointPropertyMap(),
                            GeomTraits traits=GeomTraits() )
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef std::map<face_descriptor, double>
  Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(
    internal_sdf_map);

  sdf_values<Fast_sdf_calculation_mode, TriangleMesh, boost::associative_property_map<Facet_double_map>, PointPropertyMap, GeomTraits>
  (triangle_mesh, sdf_property_map, cone_angle, number_of_rays, true, ppmap, traits);
  return segmentation_from_sdf_values<TriangleMesh, boost::associative_property_map<Facet_double_map>, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_property_map, segment_ids, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}
/// \endcond


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh.
 *
 * This function is equivalent to calling the functions `CGAL::sdf_values()` and
 * `CGAL::segmentation_from_sdf_values()` with the same parameters.
 *
 * \note There is no direct relation between the parameter `number_of_clusters`
 * and the final number of segments after segmentation. However, setting a large number of clusters will result in a detailed segmentation of the mesh with a large number of segments.
 * \note For computing segmentations of the mesh with different parameters (i.e. number of levels, and smoothing lambda),
 * it is more efficient to first compute the SDF values using `CGAL::sdf_values()` and use them in different calls to
 * `CGAL::segmentation_from_sdf_values()`.
 *
 * @pre `is_triangle_mesh(triangle_mesh)`
 * @pre `number_of_clusters > 0`
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type
 * @tparam GeomTraits a model of `SegmentationGeomTraits`
 * @tparam PointPropertyMap a `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and `GeomTraits::Point_3` as value type.
 *
 * @param triangle_mesh surface mesh on which SDF values are computed
 * @param[out] segment_ids the segment or cluster id of each facet
 * @param cone_angle opening angle in radians for the cone of each facet
 * @param number_of_rays number of rays picked in the cone of each facet. In our experiments, we observe that increasing the number of rays beyond the default has a little effect on the quality of the segmentation result
 * @param number_of_clusters number of clusters for the soft clustering
 * @param smoothing_lambda factor which indicates the importance of the surface features for the energy minimization. It is recommended to choose a value in the interval [0,1]. See the section \ref Surface_mesh_segmentationGraphCut for more details.
 * @param output_cluster_ids if `false` fill `segment_ids` with segment-ids, and with cluster-ids otherwise (see \cgalFigureRef{Cluster_vs_segment})
 * @param traits traits class
 * @param ppmap point property map. An overload is provided with `get(boost::vertex_point,triangle_mesh)` as default.
 *
 * @return number of segments if `output_cluster_ids` is set to `false` and `number_of_clusters` otherwise
 */
template < class TriangleMesh, class SegmentPropertyMap, class PointPropertyMap
#ifdef DOXYGEN_RUNNING
         = typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type
#endif
, class GeomTraits
#ifdef DOXYGEN_RUNNING
= typename Kernel_traits<typename boost::property_traits<PointPropertyMap>::value_type>::Kernel
#endif
>
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false,
                            PointPropertyMap ppmap=PointPropertyMap(),
                            GeomTraits traits=GeomTraits())
{
  return segmentation_via_sdf_values<true, TriangleMesh, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, segment_ids, cone_angle, number_of_rays, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}

#ifndef DOXYGEN_RUNNING
// we need these overloads for the default of the point property map

/// sdf_values ///
template < bool Fast_sdf_calculation_mode, class TriangleMesh, class SDFPropertyMap, class PointPropertyMap>
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true,
            PointPropertyMap ppmap = PointPropertyMap())
{
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return sdf_values<Fast_sdf_calculation_mode, TriangleMesh, SDFPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, cone_angle, number_of_rays, postprocess, ppmap, traits);
}

template < bool Fast_sdf_calculation_mode, class TriangleMesh, class SDFPropertyMap>
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<TriangleMesh&>(triangle_mesh));
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return sdf_values<Fast_sdf_calculation_mode, TriangleMesh, SDFPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, cone_angle, number_of_rays, postprocess, ppmap, traits);
}

template < class TriangleMesh, class SDFPropertyMap, class PointPropertyMap>
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true,
            PointPropertyMap ppmap = PointPropertyMap())
{
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return sdf_values<true, TriangleMesh, SDFPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, cone_angle, number_of_rays, postprocess, ppmap, traits);
}

template < class TriangleMesh, class SDFPropertyMap>
std::pair<double, double>
sdf_values( const TriangleMesh& triangle_mesh,
            SDFPropertyMap sdf_values_map,
            double cone_angle = 2.0 / 3.0 * CGAL_PI,
            std::size_t number_of_rays = 25,
            bool postprocess = true)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<TriangleMesh&>(triangle_mesh));
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return sdf_values<true, TriangleMesh, SDFPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, cone_angle, number_of_rays, postprocess, ppmap, traits);
}

/// segmentation_from_sdf_values ///
template <class TriangleMesh, class SDFPropertyMap, class SegmentPropertyMap, class PointPropertyMap>
std::size_t
segmentation_from_sdf_values(const TriangleMesh& triangle_mesh,
                             SDFPropertyMap sdf_values_map,
                             SegmentPropertyMap segment_ids,
                             std::size_t number_of_clusters = 5,
                             double smoothing_lambda = 0.26,
                             bool output_cluster_ids = false,
                             PointPropertyMap ppmap = PointPropertyMap() )
{
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_from_sdf_values<TriangleMesh, SDFPropertyMap, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, segment_ids, number_of_clusters, smoothing_lambda,
          output_cluster_ids, ppmap, traits);
}

template <class TriangleMesh, class SDFPropertyMap, class SegmentPropertyMap>
std::size_t
segmentation_from_sdf_values(const TriangleMesh& triangle_mesh,
                             SDFPropertyMap sdf_values_map,
                             SegmentPropertyMap segment_ids,
                             std::size_t number_of_clusters = 5,
                             double smoothing_lambda = 0.26,
                             bool output_cluster_ids = false)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<TriangleMesh&>(triangle_mesh));
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_from_sdf_values<TriangleMesh, SDFPropertyMap, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, sdf_values_map, segment_ids, number_of_clusters, smoothing_lambda,
          output_cluster_ids, ppmap, traits);
}

/// segmentation_via_sdf_values ///
template <bool Fast_sdf_calculation_mode, class TriangleMesh, class SegmentPropertyMap, class PointPropertyMap>
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false,
                            PointPropertyMap ppmap = PointPropertyMap() )
{
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_via_sdf_values<Fast_sdf_calculation_mode, TriangleMesh, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, segment_ids, cone_angle, number_of_rays, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}

template <bool Fast_sdf_calculation_mode, class TriangleMesh, class SegmentPropertyMap>
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<TriangleMesh&>(triangle_mesh));
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_via_sdf_values<Fast_sdf_calculation_mode, TriangleMesh, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, segment_ids, cone_angle, number_of_rays, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}

template <class TriangleMesh, class SegmentPropertyMap, class PointPropertyMap>
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false,
                            PointPropertyMap ppmap = PointPropertyMap() )
{
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_via_sdf_values<true, TriangleMesh, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, segment_ids, cone_angle, number_of_rays, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}

template <class TriangleMesh, class SegmentPropertyMap>
std::size_t
segmentation_via_sdf_values(const TriangleMesh& triangle_mesh,
                            SegmentPropertyMap segment_ids,
                            double cone_angle = 2.0 / 3.0 * CGAL_PI,
                            std::size_t number_of_rays = 25,
                            std::size_t number_of_clusters = 5,
                            double smoothing_lambda = 0.26,
                            bool output_cluster_ids = false)
{
  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<TriangleMesh&>(triangle_mesh));
  typedef typename boost::property_traits<PointPropertyMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel GeomTraits;
  GeomTraits traits;
  return segmentation_via_sdf_values<true, TriangleMesh, SegmentPropertyMap, PointPropertyMap, GeomTraits>
         (triangle_mesh, segment_ids, cone_angle, number_of_rays, number_of_clusters,
          smoothing_lambda, output_cluster_ids, ppmap, traits);
}
#endif


}//namespace CGAL

#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //
