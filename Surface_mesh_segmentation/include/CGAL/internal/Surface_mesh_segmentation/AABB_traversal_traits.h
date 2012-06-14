#ifndef CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H
#define CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H

namespace CGAL
{

/**
 * @class Closest_intersection_traits
 * Not a generic implementation !!
 * Assumes Query is Ray, and intersection with primitive is point
 * Also assumes intersection with Bbox is Segment.
 * It is returning correct minimum, but not reducing running time.
 */
template<typename AABBTraits, typename Query>
class Closest_intersection_traits
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef typename ::CGAL::AABB_tree<AABBTraits>::size_type size_type;

public:
  typedef typename boost::optional<Object_and_primitive_id> Result;

public:
  Closest_intersection_traits()
    : min_result(), min_distance((std::numeric_limits<double>::max)()) {
  }

  bool go_further() const {
    return true;
  }

  void intersection(const Query& query, const Primitive& primitive) {
    Result m_result = AABBTraits().intersection_object()(query, primitive);
    if(m_result) {
      Point i_point;
      double distance;
      if(CGAL::assign(i_point, m_result->first)) {
        distance = (query.source() - i_point).squared_length();
      }

      if(!min_result || distance < min_distance) {
        min_distance = distance;
        min_result = m_result;
      }
    }
  }

  bool do_intersect(const Query& query, const Node& node) const {
    bool intersected = AABBTraits().do_intersect_object()(query, node.bbox());
    if(!intersected || !min_result) {
      return intersected;
    }

    CGAL::Object intersection = CGAL::intersection(query, node.bbox());
    if(intersection.empty()) {
      return intersected;
    }
    //SL: warning this is another requirements for the traits!!
    typedef typename Kernel_traits<Query>::Kernel GeomTraits;
    typename GeomTraits::Segment i_segment;
    if(CGAL::assign(i_segment, intersection)) {
      double distance_1 = (i_segment.source() -
                           query.source()).squared_length(); // source returns closest intersection ?
      double distance_2 = (i_segment.target() - query.source()).squared_length();
      double distance = (CGAL::min)(distance_1, distance_2);
      return distance < min_distance;
    }
    return true;
  }

  Result result() const {
    return min_result;
  }
  bool is_intersection_found() const {
    return min_result;
  }

private:
  Result min_result;
  double min_distance;
};

/**
 * @class Special case for ray/segment-triangle
 * the only difference with the offical one (Listing_intersection_traits) is that
 * is the do_intersect which is made prior to the intersection call.
 */
#define TRAITS_USE_COUNTER
template<typename AABBTraits, typename Query, typename Output_iterator>
class Listing_intersection_traits_ray_or_segment_triangle
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point Point;
  typedef typename AABBTraits::Primitive Primitive;
  typedef typename AABBTraits::Bounding_box Bounding_box;
  typedef typename AABBTraits::Primitive::Id Primitive_id;
  typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
  typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
  typedef ::CGAL::AABB_node<AABBTraits> Node;
  typedef typename ::CGAL::AABB_tree<AABBTraits>::size_type size_type;

public:
  Listing_intersection_traits_ray_or_segment_triangle(Output_iterator out_it)
    : m_out_it(out_it) {}

  bool go_further() const {
    return true;
  }

  void intersection(const Query& query, const Primitive& primitive) {
#ifdef TRAITS_USE_COUNTER
    ++inter_counter;
#endif
    //SL: using Kernel_traits is not bad in this context cause we expect a Ray/Segment from a CGAL Kernel here
    typedef typename Kernel_traits<Query>::Kernel GeomTraits;
    if ( GeomTraits().do_intersect_3_object()(query,primitive.datum()) ) {
#ifdef TRAITS_USE_COUNTER
      ++true_inter_counter;
#endif
      boost::optional<Object_and_primitive_id> intersection;
      intersection = AABBTraits().intersection_object()(query, primitive);
      if(intersection) {
        *m_out_it++ = *intersection;
      }
    }
  }

  bool do_intersect(const Query& query, const Node& node) const {
#ifdef TRAITS_USE_COUNTER
    ++do_inter_counter;
#endif
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

#ifdef TRAITS_USE_COUNTER
  static long true_inter_counter;
  static long inter_counter;
  static long do_inter_counter;
#endif
private:
  Output_iterator m_out_it;

};
#ifdef TRAITS_USE_COUNTER
template<typename AABBTraits, typename Query, typename Output_iterator>
long Listing_intersection_traits_ray_or_segment_triangle<AABBTraits, Query, Output_iterator>::inter_counter(
  0);
template<typename AABBTraits, typename Query, typename Output_iterator>
long Listing_intersection_traits_ray_or_segment_triangle<AABBTraits, Query, Output_iterator>::do_inter_counter(
  0);
template<typename AABBTraits, typename Query, typename Output_iterator>
long Listing_intersection_traits_ray_or_segment_triangle<AABBTraits, Query, Output_iterator>::true_inter_counter(
  0);
#endif
} //namespace CGAL
#undef TRAITS_USE_COUNTER
#endif


