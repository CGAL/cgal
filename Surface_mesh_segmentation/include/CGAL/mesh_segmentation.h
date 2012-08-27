#ifndef CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H

/**
 * @file mesh_segmentation.h
 * @brief The API which contains free template functions for SDF computation and mesh segmentation.
 */
#include <CGAL/internal/Surface_mesh_segmentation/Surface_mesh_segmentation.h>

/** CGAL */
namespace CGAL
{

/*!
 * @brief Function computing the Shape Diameter Function over a surface mesh.
 *
 * This function implements the Shape Diameter Function (SDF) as described in \cite shapira2008consistent.
 * After the computation of SDF values for each facet, the following post-processing steps are applied:
 *  - Facets with no SDF values (i.e. zero) are assigned to average SDF value of its neighbors.
 * If still there is any facet which has no SDF value, minimum SDF value greater than zero is assigned to it.
 *  - Smoothed with bilateral filtering.
 *  - Linearly normalized between [0,1].
 *
 * @pre parameter @a polyhedron should consist of triangles.
 * @param polyhedron `CGAL Polyhedron` on which SDF values are computed
 * @param[out] sdf_values <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a>  with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet
 * @return minimum and maximum SDF values before linear normalization
 */
template <class Polyhedron, class SDFPropertyMap>
std::pair<double, double>
sdf_values_computation(const Polyhedron& polyhedron,
                       SDFPropertyMap sdf_values,
                       double cone_angle = 2.0 / 3.0 * CGAL_PI,
                       int number_of_rays = 25)
{
  internal::Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  return algorithm.calculate_sdf_values(cone_angle, number_of_rays, sdf_values);
}

/*!
 * @brief Function computing the segmentation of a surface mesh given an SDF value per facet.
 *
 * This function fills a property map which associates a segment-id (between [0, number of segments -1]) to each facet.
 * Formally, a segment is a set of connected facets which are placed under same cluster.
 *
 * Note that there is no direct relation between the parameter @a number_of_levels
 * and number of segments. However, large number of clusters likely to result in detailed segmentation of the mesh with large number of segments.
 *
 * @pre parameter @a polyhedron should consist of triangles.
 * @param polyhedron `CGAL Polyhedron` on which segmentation is applied
 * @param sdf_values <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadablePropertyMap.html">`ReadablePropertyMap`</a>  with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @param[out] segment_ids <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a> with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor in the interval [0,1] which indicates the importance of surface features in energy minimization
 * @return number of segments
 */
template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap>
int
surface_mesh_segmentation_from_sdf_values(const Polyhedron& polyhedron,
    SDFPropertyMap sdf_values,
    SegmentPropertyMap segment_ids,
    int number_of_levels = 5,
    double smoothing_lambda = 0.26)
{
  smoothing_lambda = (std::max)(0.0, (std::min)(1.0,
                                smoothing_lambda)); // clip into [0-1]

  internal::Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  return algorithm.partition(number_of_levels, smoothing_lambda, sdf_values,
                             segment_ids);
}

/*!
 * @brief Function computing the segmentation of a surface mesh.
 *
 * Basically this function combines CGAL::sdf_values_computation and
 * CGAL::surface_mesh_segmentation_from_sdf_values functions by computing SDF values and segmenting the mesh in one go.
 *
 * Note that for segmenting the mesh several times with different parameters (i.e. number of levels, and smoothing lambda),
 * it is wise to first compute SDF values using CGAL::sdf_values_computation,
 * and then call CGAL::surface_mesh_segmentation_from_sdf_values with the same SDF values.
 *
 * @pre parameter @a polyhedron should consist of triangles.
 * @param polyhedron `CGAL Polyhedron` on which segmentation is applied
 * @param[out] segment_ids <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/ReadWritePropertyMap.html">`ReadWritePropertyMap`</a> with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @param number_of_rays number of rays picked from cone for each facet
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor in the interval [0,1] which indicates the importance of surface features in energy minimization
 * @return number of segments
 */
template <class Polyhedron, class SegmentPropertyMap>
int
surface_mesh_segmentation(const Polyhedron& polyhedron,
                          SegmentPropertyMap segment_ids,
                          double cone_angle = 2.0 / 3.0 * CGAL_PI,
                          int number_of_rays = 25,
                          int number_of_levels = 5,
                          double smoothing_lambda = 0.26)
{
  smoothing_lambda = (std::max)(0.0, (std::min)(1.0,
                                smoothing_lambda)); // clip into [0-1]

  typedef std::map< typename Polyhedron::Facet_const_handle, double>
  Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(
    internal_sdf_map);

  sdf_values_computation(polyhedron, sdf_property_map, cone_angle,
                         number_of_rays);
  return surface_mesh_segmentation_from_sdf_values(
           polyhedron, sdf_property_map, segment_ids, number_of_levels, smoothing_lambda);
}

}//namespace CGAL

#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //