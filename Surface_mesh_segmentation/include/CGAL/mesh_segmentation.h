#ifndef CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H
/**
 * @file mesh_segmentation.h
 * The API for mesh segmentation which contains free template functions.
 */
#include <CGAL/Surface_mesh_segmentation.h>

#define CGAL_DEFAULT_CONE_ANGLE (2.0 / 3.0) * CGAL_PI /**< Default cone opening angle for ray picking */
#define CGAL_DEFAULT_NUMBER_OF_RAYS 25                /**< Default number of rays for ray picking */
#define CGAL_DEFAULT_NUMBER_OF_CLUSTERS 5
#define CGAL_DEFAULT_SMOOTHING_LAMBDA 23.0

/** CGAL */
namespace CGAL
{

/*!
 * @brief Free template function for SDF computation.
 * Returned SDF values are smoothed by bilateral filtering and linearly normalized between [0-1].
 * @param polyhedron CGAL Polyhedron on which SDF is computed
 * @param[out] sdf_values Writable Property-Map with `Polyhedron::Facet_const_handle` as key and `double` as value type
 * @param cone_angle opening angle for cone
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

template <class Polyhedron, class SDFPropertyMap, class SegmentPropertyMap>
void surface_mesh_segmentation_from_sdf_values(const Polyhedron& polyhedron,
    SDFPropertyMap sdf_values,
    SegmentPropertyMap segment_ids,
    int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
    double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA)
{
  Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  algorithm.partition(sdf_values, segment_ids, number_of_centers,
                      smoothing_lambda);
}

template <class Polyhedron, class SegmentPropertyMap>
void surface_mesh_segmentation(const Polyhedron& polyhedron,
                               SegmentPropertyMap segment_ids,
                               double cone_angle = CGAL_DEFAULT_CONE_ANGLE,
                               int number_of_rays = CGAL_DEFAULT_NUMBER_OF_RAYS,
                               int number_of_centers = CGAL_DEFAULT_NUMBER_OF_CLUSTERS,
                               double smoothing_lambda = CGAL_DEFAULT_SMOOTHING_LAMBDA)
{

  Surface_mesh_segmentation<Polyhedron> algorithm(polyhedron);
  algorithm.calculate_sdf_values(cone_angle, number_of_rays);
  //algorithm.partition(segment_ids, number_of_centers, smoothing_lambda);
}

}//namespace CGAL

#undef CGAL_DEFAULT_CONE_ANGLE
#undef CGAL_DEFAULT_NUMBER_OF_RAYS
#endif // CGAL_SURFACE_MESH_SEGMENTATION_MESH_SEGMENTATION_H //