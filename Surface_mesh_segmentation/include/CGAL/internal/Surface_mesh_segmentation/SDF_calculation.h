#ifndef CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H
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


#define CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/internal/Surface_mesh_segmentation/AABB_traversal_traits.h>
#include <CGAL/internal/Surface_mesh_segmentation/AABB_traits.h>
#include <CGAL/internal/Surface_mesh_segmentation/Disk_samplers.h>
#include <CGAL/constructions/kernel_ftC3.h>
#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>

#define CGAL_NUMBER_OF_MAD 1.5

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{

template<class ID>
struct SkipPrimitiveFunctor {

  bool operator()(const ID& skip_me) const {
    return skip_me == skip;
  }

  SkipPrimitiveFunctor(const ID& skip) : skip(skip) { }

  ID skip;
};

template<class ID>
struct FirstIntersectionVisitor {

  void operator()(const ID& /*closest*/, double /*distance*/) const {
  }
};

/**
 * @brief Responsible for calculating Shape Diameter Function over surface of the mesh.
 *
 * @tparam Polyhedron a CGAL polyhedron
 * @tparam GeomTraits a model of SegmentationGeomTraits
 */
template <
      class Polyhedron,
      class VertexPointPmap,
      class GeomTraits = typename Polyhedron::Traits,
      bool fast_bbox_intersection = true
      >
class SDF_calculation
{
//type definitions
private:

  typedef typename GeomTraits::Vector_3   Vector;
  typedef typename GeomTraits::Point_3    Point;
  typedef typename GeomTraits::Ray_3      Ray;
  typedef typename GeomTraits::Plane_3    Plane;
  typedef typename GeomTraits::Segment_3  Segment;
  typedef typename GeomTraits::FT         FT;

  typedef typename boost::graph_traits<Polyhedron>::face_iterator face_iterator;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor   face_handle;

  typedef AABB_face_graph_triangle_primitive<Polyhedron, VertexPointPmap> Primitive;
  typedef AABB_traits_SDF<GeomTraits, Primitive, fast_bbox_intersection>
  AABB_traits_internal;
  typedef typename CGAL::AABB_tree<AABB_traits_internal>                 Tree;

  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;
  typedef typename Tree::Primitive_id
  Primitive_id;

  // Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = weight.
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;

  // DiskSampling class responsible for the sampling points in a disk. It is used for generating rays in the cones. For different example see Disk_samplers.h
  typedef Vogel_disk_sampling<boost::tuple<double, double, double> >
  Default_sampler;

// member variables
private:
  GeomTraits traits;
  const Polyhedron& mesh;
  VertexPointPmap vertex_point_map;

  typename GeomTraits::Angle_3                         angle_functor;
  typename GeomTraits::Construct_scaled_vector_3       scale_functor;
  typename GeomTraits::Construct_sum_of_vectors_3      sum_functor;
  typename GeomTraits::Construct_normal_3              normal_functor;
  typename GeomTraits::Construct_translated_point_3    translated_point_functor;
  typename GeomTraits::Construct_centroid_3            centroid_functor;

  Tree tree;

  double max_diagonal;
  bool   use_diagonal;
public:
  /**
   * Construct AABB tree with a mesh.
   * @param mesh `CGAL Polyhedron` on which AABB tree constructed
   * @param build_kd_tree requirement on internal kd-tree (it is only required if find_closest_with_AABB_distance is planned to use)
   * @param use_diagonal if true: calculates diagonal of AABB tree and cast segments instead of rays using diagonal length
   * @param traits trait object
   */
  SDF_calculation(const Polyhedron& mesh,
                  VertexPointPmap vertex_point_map,
                  bool build_kd_tree = false,
                  bool use_diagonal = true, GeomTraits traits = GeomTraits())
    :
    traits(traits),
    mesh(mesh),
    vertex_point_map(vertex_point_map),
    angle_functor(traits.angle_3_object()),
    scale_functor(traits.construct_scaled_vector_3_object()),
    sum_functor(traits.construct_sum_of_vectors_3_object()),
    normal_functor(traits.construct_normal_3_object()),
    translated_point_functor(traits.construct_translated_point_3_object()),
    centroid_functor(traits.construct_centroid_3_object()),
    use_diagonal(use_diagonal) 
  {
    tree.insert(faces(mesh).first, faces(mesh).second, mesh, vertex_point_map);
    tree.build();

    if(build_kd_tree) {
      tree.accelerate_distance_queries();
    }

    if(use_diagonal) {
      CGAL::Bbox_3 bbox = tree.bbox();
      max_diagonal =
        std::sqrt(
          CGAL::squared_distanceC3( bbox.xmin(), bbox.ymin(), bbox.zmin(), bbox.xmax(),
                                    bbox.ymax(), bbox.zmax() )
        );
    }
  }

  /**
   * Calculates SDF values for each facet in a range, and stores them in @a sdf_values. Note that sdf values are neither smoothed nor normalized.
   * @tparam FacetValueMap `WritablePropertyMap` with `boost::graph_traits<Polyhedron>::face_handle` as key and `double` as value type
   * @tparam InputIterator Iterator over polyhedrons. Its value type is `pointer to polyhedron`.
   * @param facet_begin range begin
   * @param facet_end range past-the-end
   * @param cone_angle opening angle for cone, expressed in radians
   * @param number_of_rays number of rays picked from cone for each facet
   * @param[out] sdf_values
   */
  template <class FacetValueMap, class InputIterator, class DiskSampling>
  void calculate_sdf_values(
    InputIterator facet_begin,
    InputIterator facet_end,
    double cone_angle,
    std::size_t number_of_rays,
    FacetValueMap sdf_values,
    DiskSampling disk_sampler) const {
    Disk_samples_list disk_samples;
    disk_sampler(number_of_rays, std::back_inserter(disk_samples));

    for( ; facet_begin != facet_end; ++facet_begin) {
      boost::optional<double> sdf_value = calculate_sdf_value_of_facet(*facet_begin,
                                          cone_angle, true, disk_samples);

      if(sdf_value) {
        sdf_values[*facet_begin] = *sdf_value;
      } else          {
        sdf_values[*facet_begin] = -1.0;
      }
    }
  }

  /**
   * Overload for default sampling parameter
   */
  template <class FacetValueMap, class InputIterator>
  void calculate_sdf_values(
    InputIterator facet_begin,
    InputIterator facet_end,
    double cone_angle,
    std::size_t number_of_rays,
    FacetValueMap sdf_values) const {
    calculate_sdf_values(facet_begin, facet_end, cone_angle, number_of_rays,
                         sdf_values, Default_sampler());
  }

  /**
   * Cast rays inside code located around normal with apex at center.
   * Any encountered primitives are tested with `skip`, and accepted if `skip` returns false.
   * For each ray, closest encountered primitive send to `visitor`.
   *
   * \note: normal should have unit length
   */
  template<class SkipPrimitiveFunctor, class FirstIntersectionVisitor>
  boost::optional<double> calculate_sdf_value_of_point(
    Point center,
    Vector normal,
    SkipPrimitiveFunctor skip,
    FirstIntersectionVisitor visitor,
    double cone_angle,
    std::size_t number_of_rays,
    bool accept_if_acute) const {
    return calculate_sdf_value_of_point(center, normal, skip, visitor, cone_angle,
                                        number_of_rays, accept_if_acute,
                                        Default_sampler() );
  }

  /**
   * Overload for taking DiskSampling as template parameter
   */
  template<class SkipPrimitiveFunctor, class FirstIntersectionVisitor, class DiskSampling>
  boost::optional<double> calculate_sdf_value_of_point(
    Point center,
    Vector normal,
    SkipPrimitiveFunctor skip,
    FirstIntersectionVisitor visitor,
    double cone_angle,
    std::size_t number_of_rays,
    bool accept_if_acute,
    DiskSampling disk_sampler) const {
    Disk_samples_list disk_samples;
    disk_sampler(number_of_rays, std::back_inserter(disk_samples));

    return calculate_sdf_value_of_point(center, normal, skip, visitor, cone_angle,
                                        accept_if_acute, disk_samples);
  }

  /**
   * Overload for directly taking sampled points from disk as parameter
   */
  template<class SkipPrimitiveFunctor, class FirstIntersectionVisitor>
  boost::optional<double> calculate_sdf_value_of_point(
    const Point& center,
    const Vector& normal,
    SkipPrimitiveFunctor skip,
    FirstIntersectionVisitor visitor,
    double cone_angle,
    bool accept_if_acute,
    const Disk_samples_list& disk_samples) const {
    if(cone_angle < 0.0 || cone_angle > CGAL_PI) {
      CGAL_warning(false && "Cone angle is clamped between [0, CGAL_PI].");
      cone_angle = (std::min)(CGAL_PI, (std::max)(0.0, cone_angle));
    }

    Plane plane(center, normal);
    Vector v1 = plane.base1(), v2 = plane.base2();
    v1 = scale_functor(v1, FT(1.0 / CGAL::sqrt( to_double(v1.squared_length()) )));
    v2 = scale_functor(v2, FT(1.0 / CGAL::sqrt( to_double(v2.squared_length()) )));

    std::vector<std::pair<double, double> > ray_distances;
    ray_distances.reserve(disk_samples.size());

    const FT normal_multiplier( cos(cone_angle / 2.0) );
    const FT disk_multiplier  ( sin(cone_angle / 2.0) );

    const Vector& scaled_normal = scale_functor(normal, normal_multiplier);

    for(Disk_samples_list::const_iterator sample_it = disk_samples.begin();
        sample_it != disk_samples.end(); ++sample_it) {
      bool is_intersected, intersection_is_acute;
      double min_distance;
      Primitive_id closest_id;

      Vector disk_vector = sum_functor(
                             scale_functor(v1, FT(disk_multiplier * sample_it->get<0>())),
                             scale_functor(v2, FT(disk_multiplier * sample_it->get<1>())) );
      Vector ray_direction = sum_functor(scaled_normal, disk_vector);

      if(use_diagonal) {
        FT max_distance( max_diagonal / std::sqrt(to_double(
                           ray_direction.squared_length())));
        const Vector scaled_direction = scale_functor(ray_direction, max_distance);
        const Vector target_vector = sum_functor( Vector(Point(ORIGIN), center),
                                     scaled_direction);
        const Point  target_point = translated_point_functor(Point(ORIGIN),
                                    target_vector);
        Segment segment(center, target_point);

        if(traits.is_degenerate_3_object()(segment)) {
          CGAL_warning(false &&
                       "A degenerate segment is constructed. Most probable reason is using CGAL_PI as cone_angle parameter and also picking center of disk as a sample.");
        }

        boost::tie(is_intersected, intersection_is_acute, min_distance, closest_id)
          = cast_and_return_minimum(segment, skip, accept_if_acute);
      } else {
        Ray ray(center, ray_direction);

        if(traits.is_degenerate_3_object()(ray)) {
          CGAL_warning(false &&
                       "A degenerate ray is constructed. Most probable reason is using CGAL_PI as cone_angle parameter and also picking center of disk as a sample.");
        }

        boost::tie(is_intersected, intersection_is_acute, min_distance, closest_id)
          = cast_and_return_minimum(ray, skip, accept_if_acute);
      }

      if(!intersection_is_acute) {
        continue;
      }

      visitor(closest_id, min_distance);

      ray_distances.push_back(std::make_pair(min_distance, sample_it->get<2>()));
    }

    if(ray_distances.empty()) {
      return boost::optional<double>();
    }

    return boost::optional<double>(remove_outliers_and_calculate_sdf_value(
                                     ray_distances));
  }

  /**
   * Finds closest distance to center using `closest_point` query on AABB.
   * @return squared distance
   */
  double find_closest_with_AABB_distance(Point center) const {
    Point closest = tree.closest_point(center);
    return (closest - center).squared_length();
  }
private:
  /**
   * Calculates SDF value for parameter @a facet.
   * @param facet
   * @param tree AABB tree which is built from polyhedron
   * @param samples sampled points from a unit-disk which are corresponds to rays picked from cone
   * @return calculated SDF value
   */
  boost::optional<double> calculate_sdf_value_of_facet(
    face_handle facet,
    double cone_angle,
    bool accept_if_acute,
    const Disk_samples_list& disk_samples) const {
    
    const Point p1 = get(vertex_point_map,target(halfedge(facet,mesh),mesh));
    const Point p2 = get(vertex_point_map,target(next(halfedge(facet,mesh),mesh),mesh));
    const Point p3 = get(vertex_point_map,target(prev(halfedge(facet,mesh),mesh),mesh));
    const Point center  = centroid_functor(p1, p2, p3);
    Vector normal = normal_functor(p2, p1, p3);
    normal=scale_functor(normal,
                         FT(1.0/std::sqrt(to_double(normal.squared_length()))));

    CGAL::internal::SkipPrimitiveFunctor<face_handle>
    skip(facet);
    CGAL::internal::FirstIntersectionVisitor<face_handle>
    visitor;

    return calculate_sdf_value_of_point(center, normal, skip, visitor, cone_angle,
                                        accept_if_acute, disk_samples);
  }

  /**
   * Finds closest intersection for parameter @a query.
   * @param query `Segment` or `Ray` type query
   * @param tree AABB tree which includes polyhedron
   * @param facet parent facet of @a query
   * (since numerical limitations on both center calculation and intersection test, query might intersect with related facet, should be skipped in such case)
   * @return tuple of:
   *   - get<0> bool   : true if any intersection is found
   *   - get<1> bool   : true if intersection is found and is acceptable (i.e. accute angle with surface normal)
   *   - get<2> double : distance between ray/segment origin and intersection point (0.0 if get<0> is false)
   *   - get<3> Primitive_id : closest intersected primitive if get<0> is true, else Primitive_id()
   */
  template <class Query, class SkipPrimitiveFunctor> // Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double, Primitive_id> cast_and_return_minimum(
    const Query& query, SkipPrimitiveFunctor skip, bool accept_if_acute) const {
    boost::tuple<bool, bool, double, Primitive_id>
    min_distance(false, false, 0.0, Primitive_id());
    std::list<Object_and_primitive_id> intersections;

    //SL: the difference with all_intersections is that in the traversal traits, we do do_intersect before calling intersection.
    typedef  std::back_insert_iterator< std::list<Object_and_primitive_id> >
    Output_iterator;
    Listing_intersection_traits_ray_or_segment_triangle<typename Tree::AABB_traits,Query,Output_iterator>
    traversal_traits(std::back_inserter(intersections), tree.traits());
    tree.traversal(query,traversal_traits);

    Vector min_i_ray(NULL_VECTOR);
    Primitive_id min_id;
    for(typename std::list<Object_and_primitive_id>::iterator op_it =
          intersections.begin();
        op_it != intersections.end() ; ++op_it) {
      Object object = op_it->first;
      Primitive_id id = op_it->second;
      if( skip(id) ) {
        continue;
      }

      const Point* i_point;
      if(!(i_point = object_cast<Point>(&object))) {
        continue;  // continue in case of segment.
      }

      Vector i_ray(*i_point, query.source());
      double new_distance = to_double( i_ray.squared_length() );
      if(!min_distance.template get<0>()
          || new_distance < min_distance.template get<2>()) {
        min_distance.template get<3>() = id;
        min_distance.template get<2>() = new_distance;
        min_distance.template get<0>() = true;
        min_id = id;
        min_i_ray = i_ray;
      }
    }
    if(!min_distance.template get<0>()) {
      return min_distance;
    }

    if(accept_if_acute) {
      // check whether the ray makes acute angle with intersected facet
      const Point& min_v1 = get(vertex_point_map,target(halfedge(min_id,mesh),mesh));
      const Point& min_v2 = get(vertex_point_map,target(next(halfedge(min_id,mesh),mesh),mesh));
      const Point& min_v3 = get(vertex_point_map,target(prev(halfedge(min_id,mesh),mesh),mesh));
      Vector min_normal = scale_functor(normal_functor(min_v1, min_v2, min_v3), -1.0);

      if(angle_functor(translated_point_functor(Point(ORIGIN), min_i_ray),
                       Point(ORIGIN),
                       translated_point_functor(Point(ORIGIN), min_normal)) != ACUTE) {
        return min_distance;
      }
    }

    min_distance.template get<1>() = true; // founded intersection is acceptable.
    min_distance.template get<2>() = std::sqrt(min_distance.template get<2>());
    return min_distance;
  }

  /**
   * Uses Median Absolute Deviation and removes rays which don't fall into `CGAL_NUMBER_OF_MAD` * MAD.
   * Also takes weighted average of accepted rays and calculate final sdf value.
   * @param ray_distances contains distance & weight pairs for each ray
   * @return outlier removed and averaged sdf value
   */
  double remove_outliers_and_calculate_sdf_value(
    std::vector<std::pair<double, double> >& ray_distances) const {
    // pair first -> distance, second -> weight

    const std::size_t accepted_ray_count = ray_distances.size();
    if(accepted_ray_count == 0)      {
      return 0.0;
    } else if(accepted_ray_count == 1) {
      return ray_distances[0].first;
    }

    /* Calculate median sdf */
    const std::size_t half_ray_count = accepted_ray_count / 2;
    std::nth_element(ray_distances.begin(), ray_distances.begin() + half_ray_count,
                     ray_distances.end());
    double median_sdf = ray_distances[half_ray_count].first;
    if(accepted_ray_count % 2 == 0) {
      median_sdf += std::max_element(ray_distances.begin(),
                                     ray_distances.begin() + half_ray_count)->first;
      median_sdf /= 2.0;
    }

    /* Calculate median absolute deviation */
    std::vector<double> absolute_deviation;
    absolute_deviation.reserve(accepted_ray_count);
    for(std::vector<std::pair<double, double> >::iterator it =
          ray_distances.begin(); it != ray_distances.end(); ++it) {
      absolute_deviation.push_back(std::abs(it->first - median_sdf));
    }

    std::nth_element(absolute_deviation.begin(),
                     absolute_deviation.begin() + half_ray_count, absolute_deviation.end());
    double median_deviation = absolute_deviation[half_ray_count];
    if(accepted_ray_count % 2 == 0) {
      median_deviation += *std::max_element(absolute_deviation.begin(),
                                            absolute_deviation.begin() + half_ray_count);
      median_deviation /= 2.0;
    }

    /* Calculate sdf, accept rays if ray_dist - median < (median_deviation * Factor) */
    double total_weights = 0.0, total_distance = 0.0;
    for(std::vector<std::pair<double, double> >::iterator it =
          ray_distances.begin(); it != ray_distances.end(); ++it) {
      if(std::abs(it->first - median_sdf) > (median_deviation * CGAL_NUMBER_OF_MAD)) {
        continue;
      }
      total_distance += it->first * it->second;
      total_weights += it->second;
    }

    if(total_distance == 0.0) {
      return median_sdf;  // no ray is accepted, return median.
    } else                      {
      return total_distance / total_weights;
    }
  }
};
}//namespace internal
/// @endcond
}//namespace CGAL
#undef CGAL_NUMBER_OF_MAD

#endif //CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H
