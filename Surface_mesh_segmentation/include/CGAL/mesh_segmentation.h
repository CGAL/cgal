#ifndef CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H

/**
 * @file mesh_segmentation.h
 * @brief The API which contains free template functions for SDF computation and mesh segmentation.
 */
#include <CGAL/internal/Surface_mesh_segmentation/Surface_mesh_segmentation.h>
#include <boost/config.hpp>

namespace CGAL
{


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the Shape Diameter Function over a surface mesh.
 *
 * This function implements the Shape Diameter Function (SDF) as described in \cite shapira2008consistent.
 * After the computation of SDF values for each facet, the following post-processing steps are applied:
 *   - Facets with no SDF values (i.e. zero) are assigned the average SDF value of their one-ring edge-adjacent neighbors.
 *     If there is still a facet having a zero SDF value, the minimum SDF value greater than zero is assigned to it. Note that this step is not inherited from the paper.
 *     The main reason for avoiding zero SDF values is that it can obstruct log-normalization process which takes place at the beginning of segment_from_sdf_values()`.
 *   - SDF values are smoothed with bilateral filtering.
 *   - SDF values are linearly normalized between [0,1].
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Fast_sdf_calculation_mode regardless of `GeomTraits`, use inexact predicates while traversing AABB tree nodes.
 *          It is set by default to true, and can be omitted.
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SDFPropertyMap  a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param[out] sdf_values the SDF value of each facet
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet. In general, increasing the number of rays beyond the default value has little influence upon the resulting segmentation.
 * @param traits traits object
 * @return minimum and maximum SDF values before linear normalization
 */
template <bool Fast_sdf_calculation_mode, class Polyhedron,
         class SDFPropertyMap, class GeomTraits
#ifndef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
         = typename Polyhedron::Traits
#endif
         >
std::pair<double, double>
compute_sdf_values(const Polyhedron& polyhedron,
                   SDFPropertyMap sdf_values,
                   double cone_angle = 2.0 / 3.0 * CGAL_PI,
                   int number_of_rays = 25,
                   GeomTraits traits = GeomTraits())
{
  internal::Surface_mesh_segmentation<Polyhedron, GeomTraits, Fast_sdf_calculation_mode>
  algorithm(polyhedron, traits);
  return algorithm.calculate_sdf_values(cone_angle, number_of_rays, sdf_values);
}

/// @cond SKIP_IN_MANUAL
template <class Polyhedron, class SDFPropertyMap, class GeomTraits
#ifndef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
= typename Polyhedron::Traits
#endif
>
std::pair<double, double>
compute_sdf_values(const Polyhedron& polyhedron,
                   SDFPropertyMap sdf_values,
                   double cone_angle = 2.0 / 3.0 * CGAL_PI,
                   int number_of_rays = 25,
                   GeomTraits traits = GeomTraits())
{
  return compute_sdf_values<true, Polyhedron, SDFPropertyMap, GeomTraits>
         (polyhedron, sdf_values, cone_angle, number_of_rays, traits);
}
/// @endcond

/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh given an SDF value per facet.
 *
 * This function fills a property map which associates a segment-id (between [0, number of segments -1])
 * or a cluster-id (between [0, @a number_of_levels -1]) to each facet.
 * A segment is a set of connected facets which are placed under the same cluster \cgalFigureRef{Cluster_vs_segment}.
 *
 * \note Log-normalization is applied on @a sdf_values before segmentation.
 * \note There is no direct relation between the parameter @a number_of_levels
 * and the final number of segments after segmentation. However, large number of clusters likely to result in detailed segmentation of the mesh with a large number of segments.
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SDFPropertyMap  a `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param sdf_values the SDF value of each facet between [0-1]
 * @param[out] segment_ids the segment id of each facet
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor which indicates the importance of the surface features for the energy minimization. It is recommended to choose a value in the interval [0,1]. See the section \ref Surface_mesh_segmentationGraphCut for more details.
 * @param extract_segments if true fill @a segment_ids with segment-ids, otherwise fill with cluster-ids \cgalFigureRef{Cluster_vs_segment}
 * @param traits traits object
 * @return number of segments if @a extract_segments true, @a number_of_levels otherwise
 */
template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap,
         class GeomTraits
#ifndef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
         = typename Polyhedron::Traits
#endif
         >
int
segment_from_sdf_values(const Polyhedron& polyhedron,
                        SDFPropertyMap sdf_values,
                        SegmentPropertyMap segment_ids,
                        int number_of_levels = 5,
                        double smoothing_lambda = 0.26,
                        bool extract_segments = true,
                        GeomTraits traits = GeomTraits())
{
  internal::Surface_mesh_segmentation<Polyhedron, GeomTraits> algorithm(
    polyhedron, traits);
  return algorithm.partition(number_of_levels, smoothing_lambda, sdf_values,
                             segment_ids, extract_segments);
}


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh.
 *
 * This function combines `CGAL::compute_sdf_values` and
 * `CGAL::segment_from_sdf_values` functions by computing SDF values and segmenting the mesh in one go.
 *
 * \note For computing segmentations of the mesh with different parameters (i.e. number of levels, and smoothing lambda),
 * it is more efficient to first compute the SDF values using `CGAL::compute_sdf_values` and use them for each call to
 * `CGAL::segment_from_sdf_values`.
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Fast_sdf_calculation_mode regardless of `GeomTraits`, use inexact predicates while traversing AABB tree nodes.
 *          It is set by default to true, and can be omitted.
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param[out] segment_ids the segment id of each facet
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet. In general, increasing the number of rays has a little effect on the quality of the result
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor which indicates the importance of the surface features for the energy minimization. It is recommended to choose a value in the interval [0,1]. See the section \ref Surface_mesh_segmentationGraphCut for more details.
 * @param extract_segments if true fill @a segment_ids with segment-ids, otherwise fill with cluster-ids \cgalFigureRef{Cluster_vs_segment}
 * @param traits traits object
 * @return number of segments if @a extract_segments true, @a number_of_levels otherwise
 */
template < bool Fast_sdf_calculation_mode, class Polyhedron,
         class SegmentPropertyMap, class GeomTraits
#ifndef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
         = typename Polyhedron::Traits
#endif
         >
int
compute_sdf_values_and_segment(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = 2.0 / 3.0 * CGAL_PI,
                               int number_of_rays = 25,
                               int number_of_levels = 5,
                               double smoothing_lambda = 0.26,
                               bool extract_segments = true,
                               GeomTraits traits = GeomTraits())
{
  typedef std::map< typename Polyhedron::Facet_const_handle, double>
  Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(
    internal_sdf_map);

  compute_sdf_values<Fast_sdf_calculation_mode, Polyhedron, boost::associative_property_map<Facet_double_map>, GeomTraits>
  (polyhedron, sdf_property_map, cone_angle, number_of_rays, traits);
  return segment_from_sdf_values<Polyhedron, boost::associative_property_map<Facet_double_map>, SegmentPropertyMap, GeomTraits>
         (polyhedron, sdf_property_map, segment_ids, number_of_levels, smoothing_lambda,
          extract_segments, traits);
}

/// @cond SKIP_IN_MANUAL
template < class Polyhedron, class SegmentPropertyMap, class GeomTraits
#ifndef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
= typename Polyhedron::Traits
#endif
>
int
compute_sdf_values_and_segment(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = 2.0 / 3.0 * CGAL_PI,
                               int number_of_rays = 25,
                               int number_of_levels = 5,
                               double smoothing_lambda = 0.26,
                               bool extract_segments = true,
                               GeomTraits traits = GeomTraits())
{
  return compute_sdf_values_and_segment<true, Polyhedron, SegmentPropertyMap, GeomTraits>
         (polyhedron, segment_ids, cone_angle, number_of_rays, number_of_levels,
          smoothing_lambda, extract_segments, traits);
}
/// @endcond

#ifdef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
template <bool Fast_sdf_calculation_mode, class Polyhedron, class SDFPropertyMap>
std::pair<double, double>
compute_sdf_values(const Polyhedron& polyhedron,
                   SDFPropertyMap sdf_values,
                   double cone_angle = 2.0 / 3.0 * CGAL_PI,
                   int number_of_rays = 25,
                   typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values<Fast_sdf_calculation_mode, Polyhedron, SDFPropertyMap, typename Polyhedron::Traits>
         (polyhedron, sdf_values, cone_angle, number_of_rays, traits);
}

template < class Polyhedron, class SDFPropertyMap>
std::pair<double, double>
compute_sdf_values(const Polyhedron& polyhedron,
                   SDFPropertyMap sdf_values,
                   double cone_angle = 2.0 / 3.0 * CGAL_PI,
                   int number_of_rays = 25,
                   typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values<true, Polyhedron, SDFPropertyMap, typename Polyhedron::Traits>
         (polyhedron, sdf_values, cone_angle, number_of_rays, traits);
}

template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap>
int
segment_from_sdf_values(const Polyhedron& polyhedron,
                        SDFPropertyMap sdf_values,
                        SegmentPropertyMap segment_ids,
                        int number_of_levels = 5,
                        double smoothing_lambda = 0.26,
                        bool extract_segments = true,
                        typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return segment_from_sdf_values<Polyhedron, SDFPropertyMap, SegmentPropertyMap, typename Polyhedron::Traits>
         (polyhedron, sdf_values, segment_ids, number_of_levels, smoothing_lambda,
          extract_segments, traits);
}

template <bool Fast_sdf_calculation_mode, class Polyhedron, class SegmentPropertyMap>
int
compute_sdf_values_and_segment(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = 2.0 / 3.0 * CGAL_PI,
                               int number_of_rays = 25,
                               int number_of_levels = 5,
                               double smoothing_lambda = 0.26,
                               bool extract_segments = true,
                               typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values_and_segment< Fast_sdf_calculation_mode, Polyhedron, SegmentPropertyMap, typename Polyhedron::Traits>
         (polyhedron, segment_ids, cone_angle, number_of_rays, number_of_levels,
          smoothing_lambda, extract_segments, traits);
}

template <class Polyhedron, class SegmentPropertyMap>
int
compute_sdf_values_and_segment(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = 2.0 / 3.0 * CGAL_PI,
                               int number_of_rays = 25,
                               int number_of_levels = 5,
                               double smoothing_lambda = 0.26,
                               bool extract_segments = true,
                               typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values_and_segment<true, Polyhedron, SegmentPropertyMap, typename Polyhedron::Traits>
         (polyhedron, segment_ids, cone_angle, number_of_rays, number_of_levels,
          smoothing_lambda, extract_segments, traits);
}


#endif

}//namespace CGAL

#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //