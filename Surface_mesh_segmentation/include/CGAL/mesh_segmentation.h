#ifndef CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
/**
 * @file mesh_segmentation.h
 * The API which contains free template functions for SDF computation and mesh segmentation.
 */
#include <CGAL/Surface_mesh_segmentation.h>

#define CGAL_DEFAULT_CONE_ANGLE (2.0 / 3.0) * CGAL_PI /**< Default opening angle for cone */
#define CGAL_DEFAULT_NUMBER_OF_RAYS 25                /**< Default number of rays picked from cone for each facet */
#define CGAL_DEFAULT_NUMBER_OF_LEVELS 5               /**< Default number of clusters for soft clustering */
#define CGAL_DEFAULT_SMOOTHING_LAMBDA 23.0            /**< Default factor which indicates importance of surface features in energy minization*/

/** CGAL */
namespace CGAL
{

/*!
 * @brief Free template function for SDF computation.
 * After computation of SDF, applied postprocesses are as follows:
 *  - Facets with no SDF values (i.e. zero) are assigned to average SDF value of its neighbors.
 * If still there is any facet which has no SDF value, minimum SDF value greater than zero is assigned to it.
 *  - Smoothed with bilateral filtering.
 *  - Linearly normalized between [0-1].
 * @param polyhedron `CGAL Polyhedron` on which SDF is computed
 * @param[out] sdf_values `WritablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet
 */
template <class Polyhedron, class SDFPropertyMap>
void sdf_values_computation(const Polyhedron& polyhedron,
                            SDFPropertyMap sdf_values,
                            double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                            int number_of_rays = CGAL_DEFAULT_NUMBER_OF_RAYS)
{
  Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  algorithm.calculate_sdf_values(sdf_values, cone_angle, number_of_rays);
}

/*!
 * @brief Free template function for surface mesh segmentation.
 * By taking SDF values and segmentation parameters as input, a map which contains a segment-id
 * (between [0, number of segments -1]) for each facet is provided as output.
 * Formally, a segment is set of connected facets which are placed under same cluster after graph-cut.
 *
 * Note that there is no direct relation between parameter @a number_of_levels
 * and number of segments. However, large number of levels likely to result in detailed segmentation of the mesh with large number of segments.
 *
 * @param polyhedron `CGAL Polyhedron` on which segmentation is applied
 * @param sdf_values `ReadablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @param[out] segment_ids `WritablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor which indicates importance of surface features in energy minimization
 */
template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap>
void surface_mesh_segmentation_from_sdf_values(const Polyhedron& polyhedron,
    SDFPropertyMap sdf_values,
    SegmentPropertyMap segment_ids,
    int number_of_levels = CGAL_DEFAULT_NUMBER_OF_LEVELS,
    double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA)
{
  Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  algorithm.partition(sdf_values, segment_ids, number_of_levels,
                      smoothing_lambda);
}

/*!
 * @brief Free template function for SDF computation and mesh segmentation.
 * Basically this function combines CGAL::sdf_values_computation and
 * CGAL::surface_mesh_segmentation_from_sdf_values functions by computing SDF values and segmenting the mesh in one go.
 *
 * Note that for segmenting the mesh several times with different parameters (i.e. number of levels, and smoothing lambda),
 * it is wise to first compute SDF values using CGAL::sdf_values_computation,
 * and then call CGAL::surface_mesh_segmentation_from_sdf_values with same SDF values.
 *
 * @param polyhedron `CGAL Polyhedron` on which segmentation is applied
 * @param[out] segment_ids `WritablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `int` as value type
 * @param cone_angle opening angle for cone, expressed in radians
 * @param number_of_rays number of rays picked from cone for each facet
 * @param number_of_levels number of clusters for soft clustering
 * @param smoothing_lambda factor which indicates importance of surface features in energy minimization
 */
template <class Polyhedron, class SegmentPropertyMap>
void surface_mesh_segmentation(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                               int number_of_rays = CGAL_DEFAULT_NUMBER_OF_RAYS,
                               int number_of_levels = CGAL_DEFAULT_NUMBER_OF_LEVELS,
                               double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA)
{

  Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  algorithm.calculate_sdf_values(cone_angle, number_of_rays);
  algorithm.partition(segment_ids, number_of_levels, smoothing_lambda);
}

}//namespace CGAL

#undef CGAL_DEFAULT_CONE_ANGLE
#undef CGAL_DEFAULT_NUMBER_OF_RAYS
#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //