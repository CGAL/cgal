#ifndef CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H
#define CGAL_INTERNAL_SURFACE_MESH_SEGMENTATION_AABB_TRAVERSAL_TRAITS_H

namespace CGAL
{

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @class Special case for ray/segment-triangle
 * the only difference with the offical one (Listing_intersection_traits) is that
 * is the do_intersect which is made prior to the intersection call.
 */
template<typename AABBTraits, typename Query, typename Output_iterator>
class Listing_intersection_traits_ray_or_segment_triangle
{
  typedef typename AABBTraits::FT FT;
  typedef typename AABBTraits::Point_3 Point;
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
    //SL: using Kernel_traits is not bad in this context cause we expect a Ray/Segment from a CGAL Kernel here
    typedef typename Kernel_traits<Query>::Kernel GeomTraits;
    if ( GeomTraits().do_intersect_3_object()(query,primitive.datum()) ) {
      boost::optional<Object_and_primitive_id> intersection;
      intersection = AABBTraits().intersection_object()(query, primitive);
      if(intersection) {
        *m_out_it++ = *intersection;
      }
    }
  }

  bool do_intersect(const Query& query, const Node& node) const {
    return AABBTraits().do_intersect_object()(query, node.bbox());
  }

private:
  Output_iterator m_out_it;

};

/// @endcond

} //namespace CGAL
#endif


