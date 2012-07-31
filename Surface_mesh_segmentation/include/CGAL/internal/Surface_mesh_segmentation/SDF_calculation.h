#ifndef CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H
#define CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/internal/Surface_mesh_segmentation/AABB_traversal_traits.h>
#include <CGAL/internal/Surface_mesh_segmentation/AABB_const_polyhedron_triangle_primitive.h>

#include <map>
#include <vector>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/property_map/property_map.hpp>

#define CGAL_ANGLE_ST_DEV_DIVIDER 2.0
#define CGAL_ACCEPTANCE_RATE_THRESHOLD 0.5
#define CGAL_ST_DEV_MULTIPLIER 0.75

namespace CGAL
{
namespace internal
{

template <
class Polyhedron
>
class SDF_calculation
{
public:

//type definitions
protected:
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Facet  Facet;
  typedef typename Polyhedron::Facet  Vertex;
  typedef typename Kernel::Vector_3   Vector;
  typedef typename Kernel::Point_3    Point;

  typedef typename Polyhedron::Facet_const_iterator    Facet_const_iterator;

  typedef typename Kernel::Ray_3     Ray;
  typedef typename Kernel::Plane_3   Plane;
  typedef typename Kernel::Segment_3 Segment;

  typedef typename AABB_const_polyhedron_triangle_primitive<Kernel, Polyhedron>
  Primitive;
  typedef typename CGAL::AABB_tree<AABB_traits<Kernel, Primitive> >
  Tree;
  typedef typename Tree::Object_and_primitive_id
  Object_and_primitive_id;
  typedef typename Tree::Primitive_id
  Primitive_id;

  /*Sampled points from disk, t1 = coordinate-x, t2 = coordinate-y, t3 = weight (angle with cone-normal). */
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;

//Member variables
protected:
  double cone_angle;
  int    number_of_rays;

  Disk_samples_list disk_samples_sparse;
  Disk_samples_list disk_samples_dense;

  bool use_minimum_segment;
  double multiplier_for_segment;

public:
  SDF_calculation(double cone_angle, int number_of_rays)
    : cone_angle(cone_angle), number_of_rays(number_of_rays),
      use_minimum_segment(false), multiplier_for_segment(1) {
  }

  template < class FacetValueMap >
  void calculate_sdf_values(const Polyhedron& mesh, FacetValueMap sdf_values) {
    int sparse_ray_count = number_of_rays;
    int dense_ray_count = sparse_ray_count * 2;
    disk_sampling_vogel_method(disk_samples_sparse, sparse_ray_count);
    disk_sampling_vogel_method(disk_samples_dense, dense_ray_count);

    Tree tree(mesh.facets_begin(), mesh.facets_end());
    for(Facet_const_iterator facet_it = mesh.facets_begin();
        facet_it != mesh.facets_end(); ++facet_it) {
      CGAL_precondition(facet_it->is_triangle()); //Mesh should contain triangles.

      double sdf = calculate_sdf_value_of_facet(facet_it, tree, disk_samples_sparse);
      boost::put(sdf_values, facet_it, sdf);
    }
  }

protected:
  double calculate_sdf_value_of_facet(Facet_const_iterator& facet,
                                      const Tree& tree, const Disk_samples_list& samples) const {
    const Point& p1 = facet->halfedge()->vertex()->point();
    const Point& p2 = facet->halfedge()->next()->vertex()->point();
    const Point& p3 = facet->halfedge()->prev()->vertex()->point();
    Point center  = centroid(p1, p2, p3);
    Vector normal = unit_normal(p2, p1, p3); //Assuming triangles are CCW oriented.

    Plane plane(center, normal);
    Vector v1 = plane.base1(), v2 = plane.base2();
    v1 = v1 / std::sqrt(v1.squared_length());
    v2 = v2 / std::sqrt(v2.squared_length());

    //arrange_center_orientation(plane, normal, center);

    std::vector<double> ray_distances, ray_weights;
    ray_distances.reserve(number_of_rays);
    ray_weights.reserve(number_of_rays);

    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    normal = normal * length_of_normal;
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
      Vector disk_vector = v1 * sample_it->get<0>() + v2 * sample_it->get<1>();
      Vector ray_direction = normal + disk_vector;

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
        ray_direction =  ray_direction / std::sqrt(ray_direction.squared_length());
        ray_direction = ray_direction * (*segment_distance * multiplier_for_segment);

        Segment segment(center, operator+(center, ray_direction));
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

    if(acceptance_rate > CGAL_ACCEPTANCE_RATE_THRESHOLD
        || samples.size() == disk_samples_dense.size()) {
      return sdf_result;
    }
    return calculate_sdf_value_of_facet(facet, tree, disk_samples_dense);
  }

  template <class Query> //Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double> cast_and_return_minimum(
    const Query& query, const Tree& tree, Facet_const_iterator& facet) const {
    // get<0> : if any intersection is found then true
    // get<1> : if found intersection is acceptable (i.e. accute angle with surface normal) then true
    // get<2> : distance between ray/segment origin and intersection point.
    boost::tuple<bool, bool, double> min_distance(false, false, 0);
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
      Primitive_id id     = op_it->second;
      if(id == facet) {
        continue;  //Since center is located on related facet, we should skip it if there is an intersection with it.
      }

      const Point* i_point;
      if(!(i_point = object_cast<Point>(&object))) {
        continue;  // Continue in case of segment.
      }

      Vector i_ray = (query.source() - *i_point);
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
    Vector min_normal = normal(min_v1, min_v2, min_v3) * -1.0;

    if(angle(ORIGIN + min_i_ray, Point(ORIGIN), ORIGIN + min_normal) != ACUTE) {
      return min_distance;
    }
    min_distance.get<1>() = true; // founded intersection is acceptable.
    min_distance.get<2>() = std::sqrt(min_distance.get<2>());
    return min_distance;
  }

  template <class Query> //Query can be templated for just Ray and Segment types.
  boost::tuple<bool, bool, double> cast_and_return_minimum_use_closest (
    const Query& ray, const Tree& tree,
    Facet_const_iterator& facet) const {
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

    Vector min_i_ray = ray.source() - *i_point;
    const Point& min_v1 = min_id->halfedge()->vertex()->point();
    const Point& min_v2 = min_id->halfedge()->next()->vertex()->point();
    const Point& min_v3 = min_id->halfedge()->prev()->vertex()->point();
    Vector min_normal = normal(min_v1, min_v2, min_v3) * -1.0;

    if(angle(ORIGIN + min_i_ray, Point(ORIGIN), ORIGIN + min_normal) != ACUTE) {
      return min_distance;
    }
    min_distance.get<1>() = true;
    min_distance.get<2>() = std::sqrt(min_i_ray.squared_length());
    return min_distance;
  }

  boost::tuple<double, double> remove_outliers_and_calculate_sdf_value(
    std::vector<double>& ray_distances,
    std::vector<double>& ray_weights) const {
    // get<0> : sdf value
    // get<1> : not outlier ray count / total ray count
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
    st_dev = std::sqrt(st_dev / ray_distances.size());
    /* Calculate sdf, accept rays : ray_dist - median < st dev */
    int not_outlier_count = 0;
    w_it = ray_weights.begin();
    for(std::vector<double>::iterator dist_it = ray_distances.begin();
        dist_it != ray_distances.end();
        ++dist_it, ++w_it) {

      if(abs((*dist_it) - median_sdf) > (st_dev * CGAL_ST_DEV_MULTIPLIER)) {
        continue;
      }
      total_distance += (*dist_it) * (*w_it);
      total_weights += (*w_it);
      not_outlier_count++;
    }

    if(total_distance == 0.0) {
      return median_sdf;  // no ray is accepted, return median.
    }
    double sdf_res = total_distance / total_weights;
    double acceptance_rate = not_outlier_count / static_cast<double>
                             (accepted_ray_count);
    return boost::tuple<double, double>(sdf_res, acceptance_rate);
  }

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

  void disk_sampling_vogel_method(Disk_samples_list& samples, int ray_count) {
    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    const double angle_st_dev = cone_angle / CGAL_ANGLE_ST_DEV_DIVIDER;
    const double golden_ratio = 3.0 - std::sqrt(5.0);

#if 0
    for(int i = 0; i < number_of_points; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::sqrt(static_cast<double>(i) / ray_count);
      double angle = atan(R / length_of_normal);
      angle =  exp(-0.5 * (std::pow(angle / angle_st_dev, 2))); // weight
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), angle));
    }
#else
    double custom_power = 8.0 / 8.0;
    for(int i = 0; i < ray_count; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::pow(static_cast<double>(i) / ray_count, custom_power);
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), 1.0));
    }
#endif
  }
};
}//namespace internal
}//namespace CGAL
#undef CGAL_ANGLE_ST_DEV_DIVIDER
#undef CGAL_ST_DEV_MULTIPLIER
#undef CGAL_ACCEPTANCE_RATE_THRESHOLD

#endif //CGAL_SURFACE_MESH_SEGMENTATION_SDF_CALCULATION_H