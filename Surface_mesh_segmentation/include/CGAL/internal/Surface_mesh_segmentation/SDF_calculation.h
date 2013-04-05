#ifndef CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/internal/Surface_mesh_segmentation/AABB_traversal_traits.h>
#include <CGAL/internal/Surface_mesh_segmentation/Disk_samplers.h>
#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>

#define CGAL_ST_DEV_MULTIPLIER 1 //0.75

namespace CGAL
{
/// @cond CGAL_DOCUMENT_INTERNAL
namespace internal
{

/**
 * @brief Responsible for calculating Shape Diameter Function over surface of the mesh.
 *
 * @tparam Polyhedron a CGAL polyhedron
 * @tparam GeomTraits a model of SegmentationGeomTraits
 * @tparam DiskSampling class responsible for the sampling points in a disk. It is used for generating rays in the cones. For different example see Disk_samplers.h
 */
template <
class Polyhedron,
      class GeomTraits,
      class DiskSampling = Vogel_disk_sampling<boost::tuple<double, double, double> >
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

  typedef typename Polyhedron::Facet  Facet;

  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle   Facet_const_handle;

  typedef AABB_const_polyhedron_triangle_primitive<GeomTraits, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<AABB_traits<GeomTraits, Primitive> >    Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;
  typedef typename Tree::Primitive_id
  Primitive_id;

  // Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = weight.
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;

// member variables
private:
  double cone_angle;
  int    number_of_rays;

  Disk_samples_list disk_samples;

  bool use_minimum_segment;
  double multiplier_for_segment;
  GeomTraits traits;

  typename GeomTraits::Angle_3                         angle_functor;
  typename GeomTraits::Construct_scaled_vector_3       scale_functor;
  typename GeomTraits::Construct_sum_of_vectors_3      sum_functor;
  typename GeomTraits::Construct_normal_3              normal_functor;
  typename GeomTraits::Construct_unit_normal_3         unit_normal_functor;
  typename GeomTraits::Construct_translated_point_3    translated_point_functor;
  typename GeomTraits::Construct_centroid_3            centroid_functor;
public:
  /**
   * Assign default values to member variables.
   */
  SDF_calculation(GeomTraits traits)
    : traits(traits), use_minimum_segment(false), multiplier_for_segment(1),
      angle_functor(traits.angle_3_object()),
      scale_functor(traits.construct_scaled_vector_3_object()),
      sum_functor(traits.construct_sum_of_vectors_3_object()),
      normal_functor(traits.construct_normal_3_object()),
      unit_normal_functor(traits.construct_unit_normal_3_object()),
      translated_point_functor(traits.construct_translated_point_3_object()),
      centroid_functor(traits.construct_centroid_3_object()) {
  }

  /**
   * Calculates SDF values for each facet, and stores them in @a sdf_values. Note that sdf values are neither smoothed nor normalized.
   * @param mesh `CGAL Polyhedron` on which SDF values are computed
   * @param cone_angle opening angle for cone, expressed in radians
   * @param number_of_rays number of rays picked from cone for each facet
   * @param[out] sdf_values `WritablePropertyMap` with `Polyhedron::Facet_const_handle` as key and `double` as value type
   */
  template <class FacetValueMap>
  void calculate_sdf_values(const Polyhedron& mesh, double cone_angle,
                            int number_of_rays, FacetValueMap sdf_values) {
    this->cone_angle = cone_angle;
    this->number_of_rays = number_of_rays;

    disk_samples.clear();
    DiskSampling()(number_of_rays, cone_angle, std::back_inserter(disk_samples));

    Tree tree(mesh.facets_begin(), mesh.facets_end());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      sdf_values[facet_it] = calculate_sdf_value_of_facet(facet_it, tree);
    }
  }

private:
  /**
   * Calculates SDF value for parameter @a facet.
   * @param facet
   * @param tree AABB tree which is built from polyhedron
   * @param samples sampled points from a unit-disk which are corresponds to rays picked from cone
   * @return calculated SDF value
   */
  double calculate_sdf_value_of_facet(Facet_const_handle facet,
                                      const Tree& tree) const {
    const Point& p1 = facet->halfedge()->vertex()->point();
    const Point& p2 = facet->halfedge()->next()->vertex()->point();
    const Point& p3 = facet->halfedge()->prev()->vertex()->point();
    Point center  = centroid_functor(p1, p2, p3);
    Vector normal = unit_normal_functor(p2, p1, p3);

    Plane plane(center, normal);
    Vector v1 = plane.base1(), v2 = plane.base2();
    v1 = scale_functor(v1, static_cast<GeomTraits::FT>(1.0 / CGAL::sqrt(
                         v1.squared_length())));
    v2 = scale_functor(v2, static_cast<GeomTraits::FT>(1.0 / CGAL::sqrt(
                         v2.squared_length())));

    std::vector<std::pair<double, double> > ray_distances;
    ray_distances.reserve(disk_samples.size());

    const GeomTraits::FT length_of_normal = static_cast<GeomTraits::FT>( 1.0 / tan(
        cone_angle / 2.0) );
    normal = scale_functor(normal, length_of_normal);

    for(Disk_samples_list::const_iterator sample_it = disk_samples.begin();
        sample_it != disk_samples.end(); ++sample_it) {
      bool is_intersected, intersection_is_acute;
      double min_distance;
      Vector disk_vector = sum_functor(
                             scale_functor(v1, static_cast<GeomTraits::FT>(sample_it->get<0>())),
                             scale_functor(v2, static_cast<GeomTraits::FT>(sample_it->get<1>())) );
      Vector ray_direction = sum_functor(normal, disk_vector);

      Ray ray(center, ray_direction);
      boost::tie(is_intersected, intersection_is_acute,
                 min_distance) = cast_and_return_minimum(ray, tree, facet);
      if(!intersection_is_acute) {
        continue;
      }

      ray_distances.push_back(std::make_pair(min_distance, sample_it->get<2>()));
    }
    return remove_outliers_and_calculate_sdf_value(ray_distances);
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
   *   - get<2> double : distance between ray/segment origin and intersection point (0.0 if get<0> or get<1> is false)
   */
  template <class Query> // Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double> cast_and_return_minimum(
    const Query& query, const Tree& tree, Facet_const_handle facet) const {
    boost::tuple<bool, bool, double> min_distance(false, false, 0.0);
    std::list<Object_and_primitive_id> intersections;

    //SL: the difference with all_intersections is that in the traversal traits, we do do_intersect before calling intersection.
    typedef  std::back_insert_iterator< std::list<Object_and_primitive_id> >
    Output_iterator;
    Listing_intersection_traits_ray_or_segment_triangle<typename Tree::AABB_traits,Query,Output_iterator>
    traversal_traits(std::back_inserter(intersections));
    tree.traversal(query,traversal_traits);

    Vector min_i_ray;
    Primitive_id min_id;
    for(typename std::list<Object_and_primitive_id>::iterator op_it =
          intersections.begin();
        op_it != intersections.end() ; ++op_it) {
      Object object = op_it->first;
      Primitive_id id = op_it->second;
      if(id == facet) {
        continue;  // since center is located on related facet, we should skip it if there is an intersection with it.
      }

      const Point* i_point;
      if(!(i_point = object_cast<Point>(&object))) {
        continue;  // continue in case of segment.
      }

      Vector i_ray(*i_point, query.source());
      double new_distance = i_ray.squared_length();
      if(!min_distance.get<0>() || new_distance < min_distance.get<2>()) {
        min_distance.get<2>() = new_distance;
        min_distance.get<0>() = true;
        min_id = id;
        min_i_ray = i_ray;
      }
    }
    if(!min_distance.get<0>()) {
      return min_distance;
    }

    // check whether the ray makes acute angle with intersected facet
    const Point& min_v1 = min_id->halfedge()->vertex()->point();
    const Point& min_v2 = min_id->halfedge()->next()->vertex()->point();
    const Point& min_v3 = min_id->halfedge()->prev()->vertex()->point();
    Vector min_normal = scale_functor(normal_functor(min_v1, min_v2, min_v3), -1.0);

    if(angle_functor(translated_point_functor(Point(ORIGIN), min_i_ray),
                     Point(ORIGIN),
                     translated_point_functor(Point(ORIGIN), min_normal)) != ACUTE) {
      return min_distance;
    }
    min_distance.get<1>() = true; // founded intersection is acceptable.
    min_distance.get<2>() = std::sqrt(min_distance.get<2>());
    return min_distance;
  }

  /**
   * Uses Median Absolute Deviation and removes rays which don't fall into `CGAL_ST_DEV_MULTIPLIER` * MAD.
   * Also takes weighted average of accepted rays and calculate final sdf value.
   * @param ray_distances contains distance & weight pairs for each ray
   * @return outlier removed and averaged sdf value
   */
  double remove_outliers_and_calculate_sdf_value(
    std::vector<std::pair<double, double> >& ray_distances) const {
    // pair first -> distance, second -> weight

    const int accepted_ray_count = ray_distances.size();
    if(accepted_ray_count == 0)      {
      return 0.0;
    } else if(accepted_ray_count == 1) {
      return ray_distances[0].first * ray_distances[0].second;
    }

    /* Calculate median sdf */
    const int half_ray_count = accepted_ray_count / 2;
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
      if(std::abs(it->first - median_sdf) > (median_deviation *
                                             CGAL_ST_DEV_MULTIPLIER)) {
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
#undef CGAL_ST_DEV_MULTIPLIER
#undef CGAL_ACCEPTANCE_RATE_THRESHOLD

#endif //CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H