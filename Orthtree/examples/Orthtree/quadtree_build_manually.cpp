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
  using Node_data_element = empty_type;

  Orthtree_traits_empty_2(typename Self::Bbox_d bbox) : m_bbox(bbox) {};

  auto construct_point_d_from_array_object() const {
    return [](const typename Self::Array& array) -> typename Self::Point_d { return {array[0], array[1]}; };
  }
  using Construct_point_d_from_array = std::invoke_result_t<decltype(&Self::construct_point_d_from_array_object), Self>;

  auto construct_bbox_d_object() const {
    return [](const typename Self::Array& min, const typename Self::Array& max) -> typename Self::Bbox_d {
      return {min[0], min[1], max[0], max[1]};
    };
  }
  using Construct_bbox_d = std::invoke_result_t<decltype(&Self::construct_bbox_d_object), Self>;

  auto root_node_bbox_object() const {
    return [&]() -> std::pair<typename Self::Array, typename Self::Array> {
      return {{m_bbox.xmax(), m_bbox.ymax()},
              {m_bbox.xmax(), m_bbox.ymax()}};
    };
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

  std::cout << quadtree << std::endl;
  return EXIT_SUCCESS;
}
