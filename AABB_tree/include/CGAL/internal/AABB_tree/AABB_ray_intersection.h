#ifndef CGAL_AABB_RAY_INTERSECTION_H
#define CGAL_AABB_RAY_INTERSECTION_H

namespace CGAL {

template<typename AABBTraits>
template<typename Ray>
boost::optional< typename typename AABB_tree<AABBTraits>::template Intersection_and_primitive_id<Query>::Type >
AABB_tree<AABBTraits>::ray_intersection(const Ray& query) const {
  switch(size()) // copy-paste from AABB_tree::traversal
  {
  case 0: // Tree empty
    break;
  case 1: // Tree has 1 node
    /// TODO
    break;
  default: // Tree has >= 2 nodes
    if(traits().do_intersect_object()(query, root_node()->bvolume())) {
      // TODO
    } else {
      break;
    }
  }
  return boost::none;
}

}


#endif /* CGAL_AABB_RAY_INTERSECTION_H */
