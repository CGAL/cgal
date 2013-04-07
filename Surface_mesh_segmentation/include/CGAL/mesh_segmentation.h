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
 *  - Facets with no SDF values (i.e. zero) are assigned to average SDF value of its neighbors.
 * If still there is any facet which has no SDF value, minimum SDF value greater than zero is assigned to it.
 *  - Smoothed with bilateral filtering.
 *  - Linearly normalized between [0,1].
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SDFPropertyMap  a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param[out] sdf_values the sdf value of each facet
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet. In general, increasing the number of rays has a little effect on the quality of the result.
 * @param traits traits object
 * @return minimum and maximum SDF values before linear normalization
 */
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
  internal::Surface_mesh_segmentation<Polyhedron, GeomTraits> algorithm(
    polyhedron, traits);
  return algorithm.calculate_sdf_values(cone_angle, number_of_rays, sdf_values);
}


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh given an SDF value per facet.
 *
 * This function fills a property map which associates a segment-id (between [0, number of segments -1]) to each facet.
 * Formally, a segment is a set of connected facets which are placed under same cluster.
 *
 * \note There is no direct relation between the parameter @a number_of_levels
 * and number of segments. However, large number of clusters likely to result in detailed segmentation of the mesh with large number of segments.
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SDFPropertyMap  a `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param sdf_values sdf_values the sdf value of each facet
 * @param[out] segment_ids the segment id of each facet
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor in the interval [0,1] (sugggested but not forced) which indicates the importance of surface features in energy minimization
 * @param traits traits object
 * @return number of segments
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
                        GeomTraits traits = GeomTraits())
{
  internal::Surface_mesh_segmentation<Polyhedron, GeomTraits> algorithm(
    polyhedron, traits);
  return algorithm.partition(number_of_levels, smoothing_lambda, sdf_values,
                             segment_ids);
}


/*!
 * \ingroup PkgSurfaceSegmentation
 * @brief Function computing the segmentation of a surface mesh.
 *
 * Basically this function combines `CGAL::sdf_values_computation` and
 * `CGAL::surface_mesh_segmentation_from_sdf_values` functions by computing SDF values and segmenting the mesh in one go.
 *
 * \note For computing several segmentation of the mesh with different parameters (i.e. number of levels, and smoothing lambda),
 * it is advised to first compute the SDF values using `CGAL::compute_sdf_values` and use them each time you want to
 * call `CGAL::segment_from_sdf_values`.
 *
 * @pre @a polyhedron.is_pure_triangle()
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam SegmentPropertyMap a `ReadWritePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @param polyhedron surface mesh on which SDF values are computed
 * @param[out] segment_ids the segment id of each facet
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet. In general, increasing the number of rays has a little effect on the quality of the result.
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor in the interval [0,1] (sugggested but not forced) which indicates the importance of surface features in energy minimization
 * @param traits traits object
 * @return number of segments
 */
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
                               GeomTraits traits = GeomTraits())
{
  typedef std::map< typename Polyhedron::Facet_const_handle, double>
  Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(
    internal_sdf_map);

  compute_sdf_values<Polyhedron, boost::associative_property_map<Facet_double_map>, GeomTraits>
  (polyhedron, sdf_property_map, cone_angle, number_of_rays, traits);
  return segment_from_sdf_values<Polyhedron, boost::associative_property_map<Facet_double_map>, SegmentPropertyMap, GeomTraits>
         (polyhedron, sdf_property_map, segment_ids, number_of_levels, smoothing_lambda,
          traits);
}


#ifdef BOOST_NO_FUNCTION_TEMPLATE_DEFAULT_ARGS
template <class Polyhedron, class SDFPropertyMap>
std::pair<double, double>
compute_sdf_values(const Polyhedron& polyhedron,
                   SDFPropertyMap sdf_values,
                   double cone_angle = 2.0 / 3.0 * CGAL_PI,
                   int number_of_rays = 25,
                   typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values<Polyhedron, SDFPropertyMap, typename Polyhedron::Traits>
         (polyhedron, sdf_values, cone_angle, number_of_rays, traits);
}

template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap>
int
segment_from_sdf_values(const Polyhedron& polyhedron,
                        SDFPropertyMap sdf_values,
                        SegmentPropertyMap segment_ids,
                        int number_of_levels = 5,
                        double smoothing_lambda = 0.26,
                        typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return segment_from_sdf_values<Polyhedron, SDFPropertyMap, SegmentPropertyMap, typename Polyhedron::Traits>
         (polyhedron, sdf_values, segment_ids, number_of_levels, smoothing_lambda,
          traits);
}

template <class Polyhedron, class SegmentPropertyMap>
int
compute_sdf_values_and_segment(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = 2.0 / 3.0 * CGAL_PI,
                               int number_of_rays = 25,
                               int number_of_levels = 5,
                               double smoothing_lambda = 0.26,
                               typename Polyhedron::Traits traits = typename Polyhedron::Traits())
{
  return compute_sdf_values_and_segment<Polyhedron, SegmentPropertyMap, typename Polyhedron::Traits>
         (polyhedron, segment_ids, cone_angle, number_of_rays, number_of_levels,
          smoothing_lambda, traits);
}
#endif

}//namespace CGAL

#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //