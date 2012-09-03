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

#define CGAL_ACCEPTANCE_RATE_THRESHOLD 0.5
#define CGAL_ST_DEV_MULTIPLIER 0.75

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
 * @tparam DiskSampling chosen sampling method from Disk_samplers.h
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

  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Polyhedron::Facet  Vertex;

  typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle   Facet_const_handle;


  typedef AABB_const_polyhedron_triangle_primitive<GeomTraits, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<AABB_traits<GeomTraits, Primitive> >    Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;
  typedef typename Tree::Primitive_id
  Primitive_id;

  // Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = weight (angle with cone-normal).
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;

// member variables
private:
  double cone_angle;

  Disk_samples_list disk_samples_sparse;
  Disk_samples_list disk_samples_dense;

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

    const int sparse_ray_count = number_of_rays;
    const int dense_ray_count = sparse_ray_count * 2;

    DiskSampling()(sparse_ray_count, cone_angle,
                   std::back_inserter(disk_samples_sparse));
    DiskSampling()(dense_ray_count, cone_angle,
                   std::back_inserter(disk_samples_dense));

    Tree tree(mesh.facets_begin(), mesh.facets_end());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      double sdf = calculate_sdf_value_of_facet(facet_it, tree, disk_samples_sparse);
      // Note that calculate_sdf_value_of_facet call itself again with disk_samples_dense if
      // number of accepted rays after outlier removal is below CGAL_ACCEPTANCE_RATE_THRESHOLD
      sdf_values[facet_it] = sdf;
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
  double calculate_sdf_value_of_facet(Facet_const_handle& facet, const Tree& tree,
                                      const Disk_samples_list& samples) const {
    const Point& p1 = facet->halfedge()->vertex()->point();
    const Point& p2 = facet->halfedge()->next()->vertex()->point();
    const Point& p3 = facet->halfedge()->prev()->vertex()->point();
    Point center  = centroid_functor(p1, p2, p3);
    Vector normal = unit_normal_functor(p2, p1,
                                        p3); //Assuming triangles are CCW oriented.

    Plane plane(center, normal);
    Vector v1 = plane.base1(), v2 = plane.base2();
    v1 = scale_functor(v1, 1.0 / CGAL::sqrt(v1.squared_length()));
    v2 = scale_functor(v2, 1.0 / CGAL::sqrt(v2.squared_length()));

    //arrange_center_orientation(plane, normal, center);

    std::vector<double> ray_distances, ray_weights;
    ray_distances.reserve(samples.size());
    ray_weights.reserve(samples.size());

    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    normal = scale_functor(normal, length_of_normal);
    // stores segment length,
    // making it too large might cause a non-filtered bboxes in traversal,
    // making it too small might cause a miss and consecutive ray casting.
    // for now storing maximum found distance so far.

#define SHOOT_ONLY_RAYS
#ifndef SHOOT_ONLY_RAYS
    boost::optional<double> segment_distance;
#endif
    for(Disk_samples_list::const_iterator sample_it = samples.begin();
        sample_it != samples.end(); ++sample_it) {
      bool is_intersected, intersection_is_acute;
      double min_distance;
      Vector disk_vector = sum_functor(scale_functor(v1, sample_it->get<0>()),
                                       scale_functor(v2, sample_it->get<1>()));
      Vector ray_direction = sum_functor(normal, disk_vector);

#ifdef SHOOT_ONLY_RAYS
      Ray ray(center, ray_direction);
      boost::tie(is_intersected, intersection_is_acute,
                 min_distance) = cast_and_return_minimum(ray, tree, facet);
      if(!intersection_is_acute) {
        continue;
      }
#else
      // at first cast ray
      if(!segment_distance) {
        Ray ray(center, ray_direction);
        boost::tie(is_intersected, intersection_is_acute,
                   min_distance) = cast_and_return_minimum(ray, tree, facet);
        if(!intersection_is_acute) {
          continue;
        }
        segment_distance = min_distance; //first assignment of the segment_distance
      } else { // use segment_distance to limit rays as segments
        ray_direction =  scale_functor(ray_direction,
                                       1.0 / CGAL::sqrt(ray_direction.squared_length()));
        ray_direction = scale_functor(ray_direction,
                                      (*segment_distance * multiplier_for_segment));

        Segment segment(center, translated_point_functor(center, ray_direction));
        boost::tie(is_intersected, intersection_is_acute,
                   min_distance) = cast_and_return_minimum(segment, tree, facet);
        if(!is_intersected) { //no intersection is found
          //++miss_counter;
          Ray ray(segment.target(), ray_direction);
          boost::tie(is_intersected, intersection_is_acute,
                     min_distance) = cast_and_return_minimum(ray, tree, facet);
          if(!intersection_is_acute) {
            continue;
          }
          min_distance += std::sqrt(
                            segment.squared_length()); // since our source is segment.target()
        } else if(
          !intersection_is_acute) { // intersection is found, but it is not acceptable (so, do not continue ray-segment casting)
          continue;
        }


        if(use_minimum_segment) {
          if(min_distance <
              *segment_distance) { // update segment_distance (minimum / maximum)
            *segment_distance = min_distance;
          }
        } else {
          if(min_distance >
              *segment_distance) { // update segment_distance (minimum / maximum)
            *segment_distance = min_distance;
          }
        }
      }
#endif
      ray_weights.push_back(sample_it->get<2>());
      ray_distances.push_back(min_distance);
    }
    double sdf_result, acceptance_rate;
    boost::tie(sdf_result, acceptance_rate) =
      remove_outliers_and_calculate_sdf_value(ray_distances, ray_weights);
    // If after outlier removal accepted rays are below threshold, we cast more rays.
    if(acceptance_rate > CGAL_ACCEPTANCE_RATE_THRESHOLD
        || samples.size() == disk_samples_dense.size()) {
      return sdf_result;
    }
    return calculate_sdf_value_of_facet(facet, tree, disk_samples_dense);
  }

  /**
   * Finds closest intersection for parameter @a query.
   * @param query `Segment` or `Ray` type query
   * @param tree AABB tree which includes polyhedron
   * @param facet parent facet of @a query
   * (since numerical limitations on both center calculation and intersection test, query might intersect with related facet, should be skipped in such case)
   * @return tuple of:
   *   - get<0> bool   : true if any intersection is found
   *   - get<1> bool   : true if intersection is acceptable (i.e. accute angle with surface normal)
   *   - get<2> double : distance between ray/segment origin and intersection point (0.0 if get<0> or get<1> is false)
   */
  template <class Query> // Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double> cast_and_return_minimum(
    const Query& query, const Tree& tree, Facet_const_handle& facet) const {
    boost::tuple<bool, bool, double> min_distance(false, false, 0.0);
    std::list<Object_and_primitive_id> intersections;
#if 1
    //SL: the difference with all_intersections is that in the traversal traits, we do do_intersect before calling intersection.
    typedef  std::back_insert_iterator< std::list<Object_and_primitive_id> >
    Output_iterator;
    Listing_intersection_traits_ray_or_segment_triangle<typename Tree::AABB_traits,Query,Output_iterator>
    traversal_traits(std::back_inserter(intersections));
    tree.traversal(query,traversal_traits);
#else
    tree.all_intersections(query, std::back_inserter(intersections));
#endif
    Vector min_i_ray;
    Primitive_id min_id;
    for(typename std::list<Object_and_primitive_id>::iterator op_it =
          intersections.begin();
        op_it != intersections.end() ; ++op_it) {
      Object object = op_it->first;
      Primitive_id id = op_it->second;
      if(id == facet) {
        continue;  //Since center is located on related facet, we should skip it if there is an intersection with it.
      }

      const Point* i_point;
      if(!(i_point = object_cast<Point>(&object))) {
        continue;  // Continue in case of segment.
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

  template <class Query> //Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double> cast_and_return_minimum_use_closest (
    const Query& ray, const Tree& tree,
    Facet_const_handle& facet) const {
    // get<0> : if any intersection is found then true
    // get<1> : if found intersection is acceptable (i.e. accute angle with surface normal) then true
    // get<2> : distance between ray/segment origin and intersection point.
    boost::tuple<bool, bool, double> min_distance(false, false, 0);
#if 1
    Closest_intersection_traits<typename Tree::AABB_traits, Query> traversal_traits;
    tree.traversal(ray, traversal_traits);
    boost::optional<Object_and_primitive_id> intersection =
      traversal_traits.result();
#else
    boost::optional<Object_and_primitive_id> intersection = tree.any_intersection(
          ray);
#endif
    if(!intersection) {
      return min_distance;
    }
    min_distance.get<0>() = true; // intersection is found

    Object object = intersection->first;
    Primitive_id min_id = intersection->second;

    if(min_id == facet) {
      CGAL_error();  // There should be no such case, after center-facet arrangments.
    }

    const Point* i_point;
    if(!(i_point = object_cast<Point>(&object))) {
      return min_distance;
    }

    Vector min_i_ray(*i_point, ray.source());
    const Point& min_v1 = min_id->halfedge()->vertex()->point();
    const Point& min_v2 = min_id->halfedge()->next()->vertex()->point();
    const Point& min_v3 = min_id->halfedge()->prev()->vertex()->point();

    Vector min_normal = scale_functor(normal_functor(min_v1, min_v2, min_v3), -1.0);

    if(angle_functor(translated_point_functor(Point(ORIGIN), min_i_ray),
                     Point(ORIGIN),
                     translated_point_functor(Point(ORIGIN), min_normal)) != ACUTE) {
      return min_distance;
    }
    min_distance.get<1>() = true;
    min_distance.get<2>() = CGAL::sqrt(min_i_ray.squared_length());
    return min_distance;
  }

  /**
   * Removes outliers by using median as mean and filtering rays which don't fall into `CGAL_ST_DEV_MULTIPLIER` * standard deviation.
   * @param ray_distances distances associated to casted rays
   * @param ray_weights weights associated to casted rays
   * @return tuple of:
   *   - get<0> double : outlier removed and averaged sdf value
   *   - get<1> double : ratio of (not outlier ray count) / (total ray count)
   */
  boost::tuple<double, double>
  remove_outliers_and_calculate_sdf_value(std::vector<double>& ray_distances,
                                          std::vector<double>& ray_weights) const {
    int accepted_ray_count = ray_distances.size();
    if(accepted_ray_count == 0)      {
      return 0.0;
    } else if(accepted_ray_count == 1) {
      return ray_distances[0];
    }
    /* local variables */
    double total_weights = 0.0, total_distance = 0.0;
    double median_sdf = 0.0, mean_sdf = 0.0, st_dev = 0.0;
    /* Calculate mean sdf */
    std::vector<double>::iterator w_it = ray_weights.begin();
    for(std::vector<double>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end();
        ++dist_it, ++w_it) {
      total_distance += (*dist_it) * (*w_it);
      total_weights += (*w_it);
    }
    mean_sdf = total_distance / total_weights;
    total_weights = total_distance = 0.0;
    /* Calculate median sdf */
    int half_ray_count = accepted_ray_count / 2;
    std::nth_element(ray_distances.begin(), ray_distances.begin() + half_ray_count,
                     ray_distances.end());
    median_sdf = ray_distances[half_ray_count];
    if( accepted_ray_count % 2 == 0) {
      median_sdf += *std::max_element(ray_distances.begin(),
                                      ray_distances.begin() + half_ray_count);
      median_sdf /= 2.0;
    }
    /* Calculate st dev using mean_sdf as mean */
    for(std::vector<double>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end(); ++dist_it) {
      double dif = (*dist_it) - mean_sdf;
      st_dev += dif * dif;
    }
    st_dev = std::sqrt(st_dev / accepted_ray_count);
    /* Calculate sdf, accept rays if ray_dist - median < st dev */
    int not_outlier_count = 0;
    w_it = ray_weights.begin();
    for(std::vector<double>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end();
        ++dist_it, ++w_it) {
      if(std::abs((*dist_it) - median_sdf) > (st_dev * CGAL_ST_DEV_MULTIPLIER)) {
        continue;
      }
      total_distance += (*dist_it) * (*w_it);
      total_weights += (*w_it);
      not_outlier_count++;
    }
    double sdf_res;
    if(total_distance == 0.0) {
      sdf_res = median_sdf;  // no ray is accepted, return median.
    } else                      {
      sdf_res = total_distance / total_weights;
    }
    double acceptance_rate = not_outlier_count / static_cast<double>
                             (accepted_ray_count);
    return boost::tuple<double, double>(sdf_res, acceptance_rate);
  }

  /**
   * Going to be removed.
   */
  void arrange_center_orientation(const Plane& plane, const Vector& unit_normal,
                                  Point& center) const {
    /*
    *  Not sure how to decide on how much I should move center ?
    */
    double epsilon = 1e-8;
    // Option-1
    //if(plane.has_on_positive_side(center)) { return; }
    //Vector epsilon_normal = unit_normal * epsilon;
    //do {
    //    center = operator+(center, epsilon_normal);
    //} while(!plane.has_on_positive_side(center));

    //Option-2
    //if(plane.has_on_positive_side(center)) { return; }
    //double distance = sqrt(squared_distance(plane, center));
    //distance = distance > epsilon ? (distance + epsilon) : epsilon;
    //Vector distance_normal = unit_normal * distance;
    //
    //do {
    //    center = operator+(center, distance_normal);
    //} while(!plane.has_on_positive_side(center));

    //Option-3
    Ray ray(center, unit_normal);
    typename Kernel::Do_intersect_3 intersector = Kernel().do_intersect_3_object();
    if(!intersector(ray, plane)) {
      return;
    }

    Vector epsilon_normal = unit_normal * epsilon;
    do {
      center = operator+(center, epsilon_normal);
      ray = Ray(center, unit_normal);
    } while(intersector(ray, plane));
  }
};
}//namespace internal
/// @endcond
}//namespace CGAL
#undef CGAL_ST_DEV_MULTIPLIER
#undef CGAL_ACCEPTANCE_RATE_THRESHOLD

#endif //CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H