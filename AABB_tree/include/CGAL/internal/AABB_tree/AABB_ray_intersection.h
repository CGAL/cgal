#ifndef CGAL_AABB_RAY_INTERSECTION_H
#define CGAL_AABB_RAY_INTERSECTION_H

namespace CGAL {

template<typename AABBTree>
class AABB_ray_intersection {
  AABB_ray_intersection(const AABBTree& tree) : tree_(tree) {}
public:
  template<typename Ray>
  boost::optional< typename typename AABBtree<AABBTraits>::template Intersection_and_primitive_id<Query>::Type >
  ray_intersection(const Ray& query) const {
    // We hit the root, now continue on the children. Keep track of
    // nb_primitives through a variable in each Node on the stack. In
    // BVH_node::traversal this is done through the function parameter
    // nb_primitives in the recursion.
    typedef std::greater< Node_ptr_with_ft> Node_ptr_comparison;
    typedef std::priority_queue< Node_ptr_with_ft, std::vector<Node_ptr_with_ft>, Node_ptr_comparison > Heap_type;

    Heap_type pq = Heap_type(Node_ptr_comparison(std::greater<Node_ptr_with_ft>()));
    boost::optional< typename Intersection_and_primitive_id<Ray>::Type > p;
    typename AABB_traits::Intersection
      intersection_obj = traits().intersection_object();
    typename AABB_traits::Intersection_distance
      intersection_distance_obj = traits().intersection_distance_object();

    // this is not the right way to do it, but using
    // numeric_limits<FT>::{max,infinity} will not work with Epeck.
    FT t = std::numeric_limits<double>::max();
    // Start with the root node.
    pq.push(Node_ptr_with_ft(root_node(), 0, m_primitives.size()));

    while(!pq.empty() && pq.top().value < t) {
      Node_ptr_with_ft current = pq.top();
      pq.pop();

      switch(current.nb_primitives) { // almost copy-paste from BVH_node::traversal
      case 2: // Left & right child both leaves
      {

      }
      case 3: // Left child leaf, right child inner node
      {

      }
      default: // Children both inner nodes
      {
        const Node* child = &(current.node->left_child());
        boost::optional<FT> dist = intersection_distance_object(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, current.nb_primitives/2));

        child = &(current.node->right_child());
        dist = intersection_distance_object(query, child->bbox());
        if(dist)
          pq.push(Node_ptr_with_ft(child, *dist, current.nb_primitives - current.nb_primitives/2));

        break;
      }
      }
    }

    return p;
  }
private:
  const AABBTree& tree_;
  typedef typename AABBTree::FT FT;
  typedef typename AABBTree::Node Node;
  typedef typename AABBTree::FT size_type;

  struct Node_ptr_with_ft {
    Node_ptr_with_ft(const Node* node, const AABBTree::FT& value, size_type nb_primitives)
      : node(const_cast<Node*>(node)), value(value), nb_primitives(nb_primitives) {}
    Node* node;
    FT value;
    size_type nb_primitives;
    bool operator<(const Node_ptr_with_ft& other) const { return value < other.value; }
    bool operator>(const Node_ptr_with_ft& other) const { return value > other.value; }
  };

};

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
