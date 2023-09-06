#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Dimension.h>

#include <CGAL/Orthtree.h>
#include <CGAL/Orthtree_traits_2_base.h>

using Kernel = CGAL::Simple_cartesian<double>;

namespace CGAL {

struct empty_type {
};

template <typename K>
struct Orthtree_traits_empty_2 : public Orthtree_traits_2_base<K> {

  using Self = Orthtree_traits_empty_2<K>;
  using Tree = Orthtree<Self>;

  using Node_data = std::array<empty_type, 0>;

  Orthtree_traits_empty_2(typename Self::Bbox_d bbox) : m_bbox(bbox) {};

  auto construct_point_d_from_array_object() const {
    return [](const typename Self::Array& array) -> typename Self::Point_d { return {array[0], array[1]}; };
  }
  using Construct_point_d_from_array = std::invoke_result_t<decltype(&Self::construct_point_d_from_array_object), Self>;

  auto root_node_bbox_object() const {
    return [&]() -> typename Self::Bbox_d { return m_bbox; };
  }

  auto root_node_contents_object() const { return [&]() -> Node_data { return {}; }; }

  auto distribute_node_contents_object() {
    return [&](typename Tree::Node_index n, Tree& tree, const typename Self::Point_d& center) -> void {};
  }

private:

  typename Self::Bbox_d m_bbox;

};
}

using EmptyQuadtree = CGAL::Orthtree<CGAL::Orthtree_traits_empty_2<Kernel>>;

int main() {

  // Build an empty quadtree which covers the domain (-1, -1) to (1, 1)
  EmptyQuadtree quadtree{{{-1, -1, 1, 1}}};

  // Split several nodes of the tree
  quadtree.split(quadtree.root());
  quadtree.split(quadtree.node(0));
  quadtree.split(quadtree.node(3));
  quadtree.split(quadtree.node(3, 0));

  std::cout << quadtree.bbox(quadtree.root()) << std::endl;
  std::cout << quadtree << std::endl;
  return EXIT_SUCCESS;
}
