// Copyright (c) 2007-2020  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro, Simon Giraudot, Cédric Portaneri, Tong Zhao

#ifndef CGAL_ORTHTREE_H
#define CGAL_ORTHTREE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Orthtree/Split_predicates.h>
#include <CGAL/Orthtree/Traversals.h>
#include <CGAL/Orthtree/Traversal_iterator.h>
#include <CGAL/Orthtree/IO.h>

#include <CGAL/NT_converter.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Property_container.h>
#include <CGAL/property_map.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h>

#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>
#include <fstream>
#include <ostream>
#include <functional>

#include <bitset>
#include <stack>
#include <queue>
#include <vector>
#include <math.h>
#include <utility>

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

namespace Orthtree_impl {

BOOST_MPL_HAS_XXX_TRAIT_DEF(Node_data)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Squared_distance_of_element)

template <class GT, bool has_data>
struct Node_data_wrapper;

template <class GT>
struct Node_data_wrapper<GT, true>
{
  using Node_index = typename GT::Node_index;
  using Node_data = typename GT::Node_data;
  typename CGAL::Properties::Experimental::Property_container<Node_index>::template Array<Node_data>& m_node_contents;

  template <class Property_container>
  Node_data_wrapper(Property_container& node_properties)
    : m_node_contents(node_properties.template get_or_add_property<Node_data>("contents").first)
  {}

  const Node_data& operator[](Node_index n) const
  {
    return m_node_contents[n];
  }

  Node_data& operator[](Node_index n)
  {
    return m_node_contents[n];
  }
};

template <class GT>
struct Node_data_wrapper<GT, false>
{
  using Node_index = typename GT::Node_index;
  using Node_data = void*;

  template <class Property_container>
  Node_data_wrapper(Property_container&) {}

  void* operator[](Node_index) const
  {
    return nullptr;
  }
};

} // end of Orthtree_impl namespace

/*!
  \ingroup PkgOrthtreeRef

  \brief A data structure using an axis-aligned hyperrectangle
  decomposition of dD space for efficient access and
  computation.

  \details It builds a hierarchy of nodes which subdivides the space.
  Each node represents an axis-aligned hyperrectangle region of space.
  The contents of nodes depend on the traits class, non-leaf nodes also
  contain \f$2^{dim}\f$ other nodes which further subdivide the
  region.

  \sa `CGAL::Quadtree`
  \sa `CGAL::Octree`

  \tparam GeomTraits must be a model of `OrthtreeTraits` or `OrthtreeTraitsWithData`.
 */
template <typename GeomTraits>
class Orthtree {
public:
  /// \name Template Types
  /// @{
  using Traits = GeomTraits; ///< Geometry traits
  /// @}

  /// \name Traits Types
  /// @{
#ifndef DOXYGEN_RUNNING
  static inline constexpr bool has_data = Orthtree_impl::has_Node_data<GeomTraits>::value;
  static inline constexpr bool supports_neighbor_search = true;// Orthtree_impl::has_Squared_distance_of_element<GeomTraits>::value;
#else
  static inline constexpr bool has_data = bool_value; ///< `true` if `GeomTraits` is a model of `OrthtreeTraitsWithData` and `false` otherwise.
  static inline constexpr bool supports_neighbor_search = bool_value; ///< `true` if `GeomTraits` is a model of `CollectionPartitioningOrthtreeTraits` and `false` otherwise.
#endif
  static constexpr int dimension = Traits::dimension; ///< Dimension of the tree
  using Kernel = typename Traits::Kernel; ///< Kernel type.
  using FT = typename Traits::FT; ///< Number type.
  using Point = typename Traits::Point_d; ///< Point type.
  using Bbox = typename Traits::Bbox_d; ///< Bounding box type.
  using Sphere = typename Traits::Sphere_d; ///< Sphere type.
  using Adjacency = typename Traits::Adjacency; ///< Adjacency type.

  using Node_index = typename Traits::Node_index; ///< Index of a given node in the tree; the root always has index 0.
#ifndef DOXYGEN_RUNNING
  using Node_data = typename Orthtree_impl::Node_data_wrapper<Traits, has_data>::Node_data;
#else
  using Node_data = std::conditional_t<has_data,typename GeomTraits::Node_data,void*>;
#endif

  /// @}

  /// \name Public Types
  /// @{

  /*!
   * \brief Self alias for convenience.
   */
  using Self = Orthtree<Traits>;

  /*!
   * \brief Degree of the tree (number of children of non-leaf nodes).
   */
  static constexpr int degree = (2 << (dimension - 1));

  /*!
    \brief Set of bits representing this node's relationship to its parent.

    Equivalent to an array of Booleans, where index[0] is whether `x`
    is greater, index[1] is whether `y` is greater, index[2] is whether
    `z` is greater, and so on for higher dimensions if needed.
    Used to represent a node's relationship to the center of its parent.
   */
  using Local_coordinates = std::bitset<dimension>;

  /*!
    \brief Coordinates representing this node's relationship
    with the rest of the tree.

    Each value `(x, y, z, ...)` of global coordinates is calculated by doubling
    the parent's global coordinates and adding the local coordinates.
   */
  using Global_coordinates = std::array<std::uint32_t, dimension>;

  /*!
   * \brief A predicate that determines whether a node must be split when refining a tree.
   */
  using Split_predicate = std::function<bool(Node_index, const Self&)>;

  /*!
   * \brief A model of `ForwardRange` whose value type is `Node_index`.
   */
#ifdef DOXYGEN_RUNNING
  using Node_index_range = unspecified_type;
#else
  using Node_index_range = boost::iterator_range<Index_traversal_iterator<Self>>;
#endif

  /*!
   * \brief A model of `LvaluePropertyMap` with `Node_index` as key type and `T` as value type.
   */
#ifdef DOXYGEN_RUNNING
  template <class T>
  using Property_map = unspecified_type;
#else
  template <class T>
  using Property_map = Properties::Experimental::Property_array_handle<Node_index, T>;
#endif

  /// @}

private: // data members :

  using Cartesian_ranges = Orthtrees::internal::Cartesian_ranges<Traits>;
  using Node_property_container = Properties::Experimental::Property_container<Node_index>;

  template <typename T>
  using Property_array = typename Properties::Experimental::Property_container<Node_index>::template Array<T>;

  Traits m_traits; /* the tree traits */

  Node_property_container m_node_properties;
  Orthtree_impl::Node_data_wrapper<Traits, has_data> m_node_contents;
  Property_array<std::uint8_t>& m_node_depths;
  Property_array<Global_coordinates>& m_node_coordinates;
  Property_array<std::optional<Node_index>>& m_node_parents;
  Property_array<std::optional<Node_index>>& m_node_children;

  using Bbox_dimensions = std::array<FT, dimension>;
  Bbox m_bbox;
  std::vector<Bbox_dimensions> m_side_per_depth;      /* precomputed (potentially approximated) side length per node's depth */

  Cartesian_ranges cartesian_range; /* a helper to easily iterate over coordinates of points */

public:

  /// \name Constructor
  /// @{

  /*!
    \brief constructs an orthtree for a traits instance.

    The constructed orthtree has a root node with no children,
    containing the contents determined by `Construct_root_node_contents` from the traits class.
    That root node has a bounding box determined by `Construct_root_node_bbox` from the traits class,
    which typically encloses its contents.

    This single-node orthtree is valid and compatible
    with all orthtree functionality, but any performance benefits are
    unlikely to be realized until `refine()` is called.

    \param traits the traits object.
  */
  explicit Orthtree(Traits traits) :
    m_traits(traits),
    m_node_contents(m_node_properties),
    m_node_depths(m_node_properties.template get_or_add_property<std::uint8_t>("depths", 0).first),
    m_node_coordinates(m_node_properties.template get_or_add_property<Global_coordinates>("coordinates").first),
    m_node_parents(m_node_properties.template get_or_add_property<std::optional<Node_index>>("parents").first),
    m_node_children(m_node_properties.template get_or_add_property<std::optional<Node_index>>("children").first) {

    m_node_properties.emplace();

    // init bbox with first values found
    m_bbox = m_traits.construct_root_node_bbox_object()();

    // Determine dimensions of the root bbox

    Bbox_dimensions size;
    for (int i = 0; i < dimension; ++i)
    {
      size[i] = (m_bbox.max)()[i] - (m_bbox.min)()[i];
    }
    // save orthtree attributes
    m_side_per_depth.push_back(size);

    if constexpr (has_data)
      data(root()) = m_traits.construct_root_node_contents_object()();
  }

  /*!
    constructs an orthtree from a set of arguments provided to the traits constructor
   */
  template <class ... Args, class = std::enable_if_t<sizeof...(Args)>= 2>>
  explicit Orthtree(Args&& ... args)
    : Orthtree(Traits(std::forward<Args>(args)...))
  {}

  /// copy constructor
  explicit Orthtree(const Orthtree& other) :
    m_traits(other.m_traits),
    m_node_properties(other.m_node_properties),
    m_node_contents(m_node_properties),
    m_node_depths(m_node_properties.template get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.template get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.template get_property<std::optional<Node_index>>("parents")),
    m_node_children(m_node_properties.template get_property<std::optional<Node_index>>("children")),
    m_bbox(other.m_bbox), m_side_per_depth(other.m_side_per_depth) {}

  /// move constructor
  explicit Orthtree(Orthtree&& other) :
    m_traits(other.m_traits),
    m_node_properties(std::move(other.m_node_properties)),
    m_node_contents(m_node_properties),
    m_node_depths(m_node_properties.template get_property<std::uint8_t>("depths")),
    m_node_coordinates(m_node_properties.template get_property<Global_coordinates>("coordinates")),
    m_node_parents(m_node_properties.template get_property<std::optional<Node_index>>("parents")),
    m_node_children(m_node_properties.template get_property<std::optional<Node_index>>("children")),
    m_bbox(other.m_bbox), m_side_per_depth(other.m_side_per_depth)
  {
    other.m_node_properties.emplace();
  }

  /// @}


  // Non-necessary but just to be clear on the rule of 5:

  // assignment operators deleted
  Orthtree& operator=(const Orthtree& other) = delete;

  Orthtree& operator=(Orthtree&& other) = delete;

  /// \name Tree Building
  /// @{

  /*!
    \brief recursively subdivides the orthtree until it meets the given criteria.

    The split predicate should return `true` if a leaf node should be split and `false` otherwise.

    This function may be called several times with different
    predicates: in that case, nodes already split are left unaltered,
    while nodes that were not split and for which `split_predicate`
    returns `true` are split.

    \param split_predicate determines whether or not a leaf node needs to be subdivided.
   */
  void refine(const Split_predicate& split_predicate) {

    // Initialize a queue of nodes that need to be refined
    std::queue<Node_index> todo;
    todo.push(0);

    // Process items in the queue until it's consumed fully
    while (!todo.empty()) {

      // Get the next element
      auto current = todo.front();
      todo.pop();

      // Check if this node needs to be processed
      if (split_predicate(current, *this)) {

        // Split the node, redistributing its contents to its children
        split(current);

      }

      // Check if the node has children which need to be processed
      if (!is_leaf(current)) {

        // Process each of its children
        for (int i = 0; i < degree; ++i)
          todo.push(child(current, i));
      }
    }
  }

  /*!
    \brief convenience overload that refines an orthtree using a
    maximum depth and maximum number of contained elements in a node as split
    predicate.

    This is equivalent to calling
    `refine(Orthtrees::Maximum_depth_and_maximum_contained_elements(max_depth,
    bucket_size))`.

    The refinement is stopped as soon as one of the conditions is
    violated: if a node contains more elements than `bucket_size` but is
    already at `max_depth`, it is not split. Similarly, a node that is
    at a depth smaller than `max_depth` but already contains fewer elements
    than `bucket_size`, it is not split.

    \warning This convenience method is only appropriate for trees with traits classes where
    `Node_data` is a model of `Range`. `RandomAccessRange` is suggested for performance.

    \param max_depth deepest a tree is allowed to be (nodes at this depth will not be split).
    \param bucket_size maximum number of items a node is allowed to contain.
   */
  void refine(size_t max_depth = 10, size_t bucket_size = 20) {
    refine(Orthtrees::Maximum_depth_and_maximum_contained_elements(max_depth, bucket_size));
  }

  /*!
    \brief refines the orthtree such that the difference of depth
    between two immediate neighbor leaves is never more than 1.

    This is done only by adding nodes, nodes are never removed.
   */
  void grade() {

    // Collect all the leaf nodes
    std::queue<Node_index> leaf_nodes;
    for (Node_index leaf: traverse(Orthtrees::Leaves_traversal<Self>(*this))) {
      leaf_nodes.push(leaf);
    }

    // Iterate over the nodes
    while (!leaf_nodes.empty()) {

      // Get the next node
      Node_index node = leaf_nodes.front();
      leaf_nodes.pop();

      // Skip this node if it isn't a leaf anymore
      if (!is_leaf(node))
        continue;

      // Iterate over each of the neighbors
      for (int direction = 0; direction < 6; ++direction) {

        // Get the neighbor
        auto neighbor = adjacent_node(node, direction);

        // If it doesn't exist, skip it
        if (!neighbor)
          continue;

        // Skip if this neighbor is a direct sibling (it's guaranteed to be the same depth)
        // TODO: This check might be redundant, if it doesn't affect performance maybe I could remove it
        if (parent(*neighbor) == parent(node))
          continue;

        // If it's already been split, skip it
        if (!is_leaf(*neighbor))
          continue;

        // Check if the neighbor breaks our grading rule
        // TODO: could the rule be parametrized?
        if ((depth(node) - depth(*neighbor)) > 1) {

          // Split the neighbor
          split(*neighbor);

          // Add newly created children to the queue
          for (int i = 0; i < degree; ++i) {
            leaf_nodes.push(child(*neighbor, i));
          }
        }
      }
    }
  }

  /// @}

  /// \name Accessors
  /// @{

  /*!
   * \brief provides direct read-only access to the tree traits.
   */
  const Traits& traits() const { return m_traits; }

  /*!
    \brief provides access to the root node, and by
    extension the rest of the tree.
   */
  Node_index root() const { return 0; }

  /*!
    \brief returns the deepest level reached by a leaf node in this tree (root being level 0).
   */
  std::size_t depth() const { return m_side_per_depth.size() - 1; }

  /*!
    \brief constructs a node index range using a tree-traversal function.

    This method allows iteration over the nodes of the tree with a
    user-selected order (preorder, postorder, leaves-only, etc.).

    \tparam Traversal a model of `OrthtreeTraversal`

    \param traversal class defining the traversal strategy

    \return a `ForwardRange` over the node indices of the tree
   */
  template <typename Traversal>
  Node_index_range traverse(Traversal traversal) const {

    Node_index first = traversal.first_index();

    auto next = [=](const Self&, Node_index index) -> std::optional<Node_index> {
      return traversal.next_index(index);
    };

    return boost::make_iterator_range(Index_traversal_iterator<Self>(*this, first, next),
                                      Index_traversal_iterator<Self>());
  }


  /*!
    \brief convenience method for using a traversal without constructing it yourself

    \tparam Traversal a model of `OrthtreeTraversal`

    \param args Arguments to to pass to the traversal's constructor, excluding the first (always an orthtree reference)

    \return a `ForwardRange` over the node indices of the tree
   */
  template <typename Traversal, typename ...Args>
  Node_index_range traverse(Args&& ...args) const {
    return traverse(Traversal{*this, std::forward<Args>(args)...});
  }

  // TODO shall we document it?
  FT
  compute_cartesian_coordinate(std::uint32_t gc, std::size_t depth, int ci) const
  {
    CGAL_assertion(depth <= m_side_per_depth.size());
    // an odd coordinate will be first compute at the current depth,
    // while an even coordinate has already been computed at a previous depth.
    // So while the coordinate is even, we decrease the depth to end up of the first
    // non-even coordinate to compute it (with particular case for bbox limits).
    // Note that is depth becomes too large, we might end up with incorrect coordinates
    // due to rounding errors.
    if (gc == (1u << depth)) return (m_bbox.max)()[ci]; // gc == 2^node_depth
    if (gc == 0) return (m_bbox.min)()[ci];
    if (gc % 2 !=0)
    {
      FT size = depth < m_side_per_depth.size()
              ? m_side_per_depth[depth][ci]
              : m_side_per_depth[depth-1][ci]/FT(2);
      return (m_bbox.min)()[ci] + int(gc) * size;
    }
    std::size_t nd = depth;
    do{
      --nd;
      gc = gc >> 1;
    }
    while((gc&1)==0); // while even, shift
    return (m_bbox.min)()[ci] + int(gc) * m_side_per_depth[nd][ci];
  }

  /*!
    \brief constructs the bounding box of a node.

    \note The object constructed is not the bounding box of the node's contents,
    but the bounding box of the node itself.

    \param n node to generate a bounding box for

    \return the bounding box of the node n
   */
  Bbox bbox(Node_index n) const {
    using Cartesian_coordinate = std::array<FT, dimension>;
    Cartesian_coordinate min_corner, max_corner;
    std::size_t node_depth = depth(n);

    for (int i = 0; i < dimension; i++)
    {
      min_corner[i]=compute_cartesian_coordinate(global_coordinates(n)[i], node_depth, i);
      max_corner[i]=compute_cartesian_coordinate(global_coordinates(n)[i]+1, node_depth, i);
    }
    return {std::apply(m_traits.construct_point_d_object(), min_corner),
            std::apply(m_traits.construct_point_d_object(), max_corner)};
  }

  /// @}

  /// \name Custom Properties
  /// @{

  /*!
    \brief gets a property for nodes, adding it if it does not already exist.

    \tparam T the type of the property to add

    \param name the name of the new property
    \param default_value the default value assigned to nodes for this property

    \return pair of the property map and a Boolean which is `true` if the property needed to be created
   */
  template <typename T>
  std::pair<Property_map<T>, bool> add_property(const std::string& name, const T default_value = T()) {
    auto p = m_node_properties.get_or_add_property(name, default_value);
    return std::pair<Property_map<T>, bool>(Property_map<T>(p.first), p.second);
  }

  /*!
    \brief gets a property of the nodes if it exists.

    \tparam T the type of the property to retrieve

    \param name the name of the property

    \return an optional containing the property map if it exists
   */
  template <typename T>
  std::optional<Property_map<T>> property(const std::string& name) const {
    auto p = m_node_properties.template get_property_if_exists<T>(name);
    if (p)
      return std::optional<Property_map<T> >(Property_map<T>(*p));
    else
      return std::nullopt;
  }

  /*!
    \brief returns a vector of all property names.
   */
  std::vector<std::string> properties() const {
    return m_node_properties.properties();
  }

  /*!
    \brief removes the node property from the tree.

    \tparam T the type of the property to remove

    \param property the property to be removed from the tree.

    \return true if property was a valid property of the tree.
   */
  template <typename T>
  bool remove_property(Property_map<T> property) {
    return m_node_properties.remove_property(property.array());
  }

  /// @}

  /// \name Queries
  /// @{

  /*!
    \brief finds the leaf node which contains a particular point in space.

    Traverses the orthtree and finds the leaf cell that has a
    domain enclosing the point passed. The point passed must be within
    the region enclosed by the orthtree (bbox of the root node). The point is contained in the
    lower cell of each direction if its coordinate is lower than the center.

    \param point query point.

    \return the index of the node which contains the point.
   */
  Node_index locate(const Point& point) const {

    // Make sure the point is enclosed by the orthtree
    CGAL_precondition (CGAL::do_intersect(point, bbox(root())));

    // Start at the root node
    Node_index node_for_point = root();

    // Descend the tree until reaching a leaf node
    while (!is_leaf(node_for_point)) {

      // Find the point to split around
      Point center = barycenter(node_for_point);

      // Find the index of the correct sub-node
      Local_coordinates local_coords;
      std::size_t dim = 0;
      for (const auto& r : cartesian_range(center, point))
        local_coords[dim++] = (get<0>(r) <= get<1>(r));

      // Find the correct sub-node of the current node
      node_for_point = child(node_for_point, local_coords.to_ulong());
    }

    // Return the result
    return node_for_point;
  }

  /*!
  \brief finds the `k` nearest neighbors of the point `query`.

  Nearest neighbors are outputted in order of increasing distance to `query`.

  \tparam OutputIterator a model of `OutputIterator` that accepts `GeomTraits::Node_data_element` objects.

  \param query query point
  \param k number of neighbors to find
  \param output output iterator

  \warning Nearest neighbor searches requires `GeomTraits` to be a model of `CollectionPartitioningOrthtreeTraits`.
 */
  template<typename OutputIterator>
  auto nearest_k_neighbors(const Point& query,
    std::size_t k,
    OutputIterator output) const -> std::enable_if_t<supports_neighbor_search, OutputIterator> {
    Sphere query_sphere(query, (std::numeric_limits<FT>::max)());
    CGAL_precondition(k > 0);

    return nearest_k_neighbors_within_radius(query_sphere, k, output);
  }

  /*!
  \brief finds the elements in the sphere `query`.

  Elements are outputted in order of increasing distance to
  the center of the sphere.

  \tparam OutputIterator a model of `OutputIterator` that accepts `GeomTraits::Node_data_element` objects.

  \param query query sphere
  \param output output iterator

  \warning Nearest neighbor searches requires `GeomTraits` to be a model of `CollectionPartitioningOrthtreeTraits`.
 */
  template<typename OutputIterator>
  auto neighbors_within_radius(const Sphere& query, OutputIterator output) const -> std::enable_if_t<supports_neighbor_search, OutputIterator> {
    return nearest_k_neighbors_within_radius(query, (std::numeric_limits<std::size_t>::max)(), output);
  }

  /*!
  \brief finds at most `k` elements within a specific radius that are
  nearest to the center of the sphere `query`: if `query` does not contain
  at least `k` elements, only contained elements will be returned.

  This function is useful when the user already knows how sparse the elements are,
  or if they do not care about elements that are too far away.
  Setting a small radius may have performance benefits.

  \tparam OutputIterator must be a model of `OutputIterator` that accepts `GeomTraits::Node_data_element`

  \param query the region to search within
  \param k the number of elements to find
  \param output the output iterator to add the found elements to (in order of increasing distance)

  \warning Nearest neighbor searches requires `GeomTraits` to be a model of `CollectionPartitioningOrthtreeTraits`.
 */
  template <typename OutputIterator>
  auto nearest_k_neighbors_within_radius(
    const Sphere& query,
    std::size_t k,
    OutputIterator output
  ) const -> std::enable_if_t<supports_neighbor_search, OutputIterator> {
    CGAL_precondition(k > 0);
    Sphere query_sphere = query;

    // todo: this type is over-constrained, this must be made more generic
    struct Node_element_with_distance {
      typename Traits::Node_data_element element;
      FT distance;
    };

    // Create an empty list of elements
    std::vector<Node_element_with_distance> element_list;
    if (k != (std::numeric_limits<std::size_t>::max)())
      element_list.reserve(k);

    // Invoking the recursive function adds those elements to the vector (passed by reference)
    nearest_k_neighbors_recursive(query_sphere, root(), element_list, k);

    // Add all the points found to the output
    for (auto& item : element_list)
      *output++ = item.element;

    return output;
  }

  /*!
    \brief finds the leaf nodes that intersect with any primitive.

    \note this function requires the function
    `bool CGAL::do_intersect(QueryType, Traits::Bbox_d)` to be defined.

    This function finds all the intersecting leaf nodes and writes their indices to the output iterator.

    \tparam Query the primitive class (e.g., sphere, ray)
    \tparam OutputIterator a model of `OutputIterator` that accepts `Node_index` types

    \param query the intersecting primitive.
    \param output output iterator.

    \return the output iterator after writing
   */
  template <typename Query, typename OutputIterator>
  OutputIterator intersected_nodes(const Query& query, OutputIterator output) const {
    return intersected_nodes_recursive(query, root(), output);
  }

  /// @}

  /// \name Operators
  /// @{

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.

    Trees may be considered equivalent even if they have different contents.
    Equivalent trees must have the same root bounding box and the same node structure.

    \param rhs the other orthtree

    \return `true` if the trees have the same topology, and `false` otherwise
   */
  bool operator==(const Self& rhs) const {

    // Identical trees should have the same bounding box
    if (rhs.m_bbox != m_bbox || rhs.m_side_per_depth[0] != m_side_per_depth[0])
      return false;

    // Identical trees should have the same depth
    if (rhs.depth() != depth())
      return false;

    // If all else is equal, recursively compare the trees themselves
    return is_topology_equal(*this, rhs);
  }

  /*!
    \brief compares the topology of the orthtree with that of `rhs`.

    \param rhs the other orthtree

    \return `false` if the trees have the same topology, and `true` otherwise
   */
  bool operator!=(const Self& rhs) const {
    return !operator==(rhs);
  }

  /// @}

  /// \name Node Access
  /// @{

  /*!
    \brief determines whether the node specified by index `n` is a leaf node.
   */
  bool is_leaf(Node_index n) const {
    return !m_node_children[n].has_value();
  }

  /*!
    \brief determines whether the node specified by index `n` is the root node.
   */
  bool is_root(Node_index n) const {
    return n == 0;
  }

  /*!
    \brief determines the depth of the node specified.

    The root node has depth 0, its children have depth 1, and so on.

    \param n index of the node to check.

    \return the depth of the node n within its tree.
   */
  std::size_t depth(Node_index n) const {
    return m_node_depths[n];
  }

  /*!
    \brief retrieves a reference to the `Node_data` associated with the node specified by `n` if
    `GeomTraits` is a model of `OrthtreeTraitswithData`, and `nullptr` otherwise.
   */
  std::conditional_t<has_data,Node_data&,void*>& data(Node_index n){
    return m_node_contents[n];
  }

  /*!
    \brief retrieves a const reference to the `Node_data` associated with the node specified by `n` if
    `GeomTraits` is a model of `OrthtreeTraitswithData`, and `nullptr` otherwise.
   */
  std::conditional_t<has_data,const Node_data&,void*> data(Node_index n) const{
    return m_node_contents[n];
  }

  /*!
    \brief retrieves the global coordinates of the node.
   */
  Global_coordinates global_coordinates(Node_index n) const {
    return m_node_coordinates[n];
  }

  /*!
    \brief retrieves the local coordinates of the node.
   */
  Local_coordinates local_coordinates(Node_index n) const {
    Local_coordinates result;
    for (std::size_t i = 0; i < dimension; ++i)
      result[i] = global_coordinates(n)[i] & 1;
    return result;
  }

  /*!
    \brief returns this n's parent.

    \pre `!is_root()`

    \param n index of the node to retrieve the parent of

    \return the index of the parent of node n
   */
  Node_index parent(Node_index n) const {
    CGAL_precondition (!is_root(n));
    return *m_node_parents[n];
  }

  /*!
    \brief returns a node's `i`th child.

    \pre `!is_leaf()`

    \param n index of the node to retrieve the child of
    \param i in [0, 2^D) specifying the child to retrieve

    \return the index of the `i`th child of node n
   */
  Node_index child(Node_index n, std::size_t i) const {
    CGAL_precondition (!is_leaf(n));
    return *m_node_children[n] + i;
  }

  /*!
    \brief retrieves an arbitrary descendant of the node specified by `node`.

    Convenience function to avoid the need to call `orthtree.child(orthtree.child(node, 0), 1)`.

    Each index in `indices` specifies which child to enter as descending the tree from `node` down.
    Indices are evaluated in the order they appear as parameters, so
    `descendant(root, 0, 1)` returns the second child of the first child of the root.

    \param node the node to descend
    \param indices the integer indices specifying the descent to perform

    \return the index of the specified descendant node
   */
  template <typename... Indices>
  Node_index descendant(Node_index node, Indices... indices) {
    return recursive_descendant(node, indices...);
  }

  /*!
    \brief convenience function for retrieving arbitrary nodes.

    Equivalent to `tree.descendant(tree.root(), indices...)`.

    \param indices the integer indices specifying the descent to perform, starting from the root

    \return the index of the specified node
   */
  template <typename... Indices>
  Node_index node(Indices... indices) {
    return descendant(root(), indices...);
  }

  /*!
    \brief finds the next sibling in the parent of the node specified by the index `n`.

    Traverses the tree in increasing order of local index (e.g., 000, 001, 010, etc.)

    \param n the index of the node to find the sibling of

    \return the index of the next sibling of n
    if n is not the last node in its parent, otherwise `std::nullopt`.
   */
  const std::optional<Node_index> next_sibling(Node_index n) const {

    // Root node has no siblings
    if (is_root(n)) return {};

    // Find out which child this is
    std::size_t local_coords = local_coordinates(n).to_ulong();

    // The last child has no more siblings
    if (int(local_coords) == degree - 1)
      return {};

    // The next sibling is the child of the parent with the following local coordinates
    return child(parent(n), local_coords + 1);
  }

  /*!
    \brief finds the next sibling of the parent of the node specified by `n` if it exists.

    \param n the index node to find the sibling up of.

    \return The index of the next sibling of the parent of n
    if n is not the root and its parent has a sibling, otherwise nothing.
   */
  const std::optional<Node_index> next_sibling_up(Node_index n) const {

    // the root node has no next sibling up
    if (n == 0) return {};

    auto up = std::optional<Node_index>{parent(n)};
    while (up) {

      if (next_sibling(*up)) return {next_sibling(*up)};

      up = is_root(*up) ? std::optional<Node_index>{} : std::optional<Node_index>{parent(*up)};
    }

    return {};
  }

  /*!
    \brief finds the leaf node reached when descending the tree and always choosing child 0.

    This is the starting point of a depth-first traversal.

    \param n the index of the node to find the deepest first child of.

    \return the index of the deepest first child of node n.
   */
  Node_index deepest_first_child(Node_index n) const {

    auto first = n;
    while (!is_leaf(first))
      first = child(first, 0);

    return first;
  }

  /*!
    \brief finds node reached when descending the tree to a depth `d` and always choosing child 0.

    Similar to `deepest_first_child()`, but does go to a fixed depth.

    \param n the index of the node to find the `d`th first child of.
    \param d the depth to descend to.

    \return the index of the `d`th first child, nothing if the tree is not deep enough.
   */
  std::optional<Node_index> first_child_at_depth(Node_index n, std::size_t d) const {

    std::queue<Node_index> todo;
    todo.push(n);

    while (!todo.empty()) {
      Node_index node = todo.front();
      todo.pop();

      if (depth(node) == d)
        return node;

      if (!is_leaf(node))
        for (int i = 0; i < degree; ++i)
          todo.push(child(node, i));
    }

    return {};
  }

  /*!
    \brief splits a node into subnodes.

    Only leaf nodes should be split.
    When a node is split it is no longer a leaf node.
    The full set of `degree` children are constructed automatically, and their values are set.
    Contents of this node are _not_ propagated automatically, this is responsibility of the
    `distribute_node_contents_object` in the traits class.

    \param n index of the node to split
 */
  void split(Node_index n) {

    // Make sure the node hasn't already been split
    CGAL_precondition (is_leaf(n));

    // Split the node to create children
    using Local_coordinates = Local_coordinates;
    m_node_children[n] = m_node_properties.emplace_group(degree);
    for (std::size_t i = 0; i < degree; i++) {

      Node_index c = *m_node_children[n] + i;

      // Make sure the node isn't one of its own children
      CGAL_assertion(n != *m_node_children[n] + i);

      Local_coordinates local_coordinates{i};
      for (int i = 0; i < dimension; i++)
        m_node_coordinates[c][i] = (2 * m_node_coordinates[n][i]) + local_coordinates[i];
      m_node_depths[c] = m_node_depths[n] + 1;
      m_node_parents[c] = n;
    }

    // Check if we've reached a new max depth
    if (depth(n) + 1 == m_side_per_depth.size()) {
      // Update the side length map with the dimensions of the children
      Bbox_dimensions size = m_side_per_depth.back();
      Bbox_dimensions child_size;
      for (int i = 0; i < dimension; ++i)
        child_size[i] = size[i] / FT(2);
      m_side_per_depth.push_back(child_size);
    }

    // Find the point around which the node is split
    Point center = barycenter(n);

    // Add the node's contents to its children
    if constexpr (has_data)
      m_traits.distribute_node_contents_object()(n, *this, center);
  }

  /*!
   * \brief returns the center point of a node.
   *
   * @param n index of the node to find the center point for
   *
   * @return the center point of node n
   */
  Point barycenter(Node_index n) const {
    std::size_t node_depth = depth(n);
    // the barycenter is computed as the lower corner of the lexicographically top child node
    std::array<FT, dimension> bary;
    for (int i = 0; i < dimension; i++)
      bary[i] = compute_cartesian_coordinate(2 * global_coordinates(n)[i]+1, node_depth+1, i);

    return std::apply(m_traits.construct_point_d_object(), bary);
  }

  /*!
    \brief determines whether a pair of subtrees have the same topology.

    \param lhsNode index of a node in lhsTree
    \param lhsTree an Orthtree
    \param rhsNode index of a node in rhsTree
    \param rhsTree another Orthtree

    @return true if lhsNode and rhsNode have the same topology, false otherwise
   */
  static bool is_topology_equal(Node_index lhsNode, const Self& lhsTree, Node_index rhsNode, const Self& rhsTree) {

    // If one node is a leaf, and the other isn't, they're not the same
    if (lhsTree.is_leaf(lhsNode) != rhsTree.is_leaf(rhsNode))
      return false;

    // If both nodes are non-leaf
    if (!lhsTree.is_leaf(lhsNode)) {

      // Check all the children
      for (int i = 0; i < degree; ++i) {
        // If any child cell is different, they're not the same
        if (!is_topology_equal(lhsTree.child(lhsNode, i), lhsTree,
                               rhsTree.child(rhsNode, i), rhsTree))
          return false;
      }
    }

    return (lhsTree.global_coordinates(lhsNode) == rhsTree.global_coordinates(rhsNode));
  }

  /*!
    \brief helper function for calling `is_topology_equal()` on the root nodes of two trees.

    \param lhs an Orthtree
    \param rhs another Orthtree

    \return `true` if `lhs` and `rhs` have the same topology, and `false` otherwise
   */
  static bool is_topology_equal(const Self& lhs, const Self& rhs) {
    return is_topology_equal(lhs.root(), lhs, rhs.root(), rhs);
  }

  /*!
    \brief finds the directly adjacent node in a specific direction

    \pre `direction.to_ulong < 2 * dimension`

    Adjacent nodes are found according to several properties:
    - adjacent nodes may be larger than the seek node, but never smaller
    - a node has at most `2 * dimension` different adjacent nodes (in 3D: left, right, up, down, front, back)
    - adjacent nodes are not required to be leaf nodes

    Here's a diagram demonstrating the concept for a quadtree:

    ```
    +---------------+---------------+
    |               |               |
    |               |               |
    |               |               |
    |       A       |               |
    |               |               |
    |               |               |
    |               |               |
    +-------+-------+---+---+-------+
    |       |       |   |   |       |
    |   A   |  (S)  +---A---+       |
    |       |       |   |   |       |
    +---+---+-------+---+---+-------+
    |   |   |       |       |       |
    +---+---+   A   |       |       |
    |   |   |       |       |       |
    +---+---+-------+-------+-------+
    ```

    - (S) : Seek node
    - A  : Adjacent node

    Note how the top adjacent node is larger than the seek node.  The
    right adjacent node is the same size, even though it contains
    further subdivisions.

    This implementation returns the adjacent node if it's found.  If
    there is no adjacent node in that direction, it returns a null
    node.

    \param n index of the node to find a neighbor of
    \param direction which way to find the adjacent node relative to
    this one. Each successive bit selects the direction for the
    corresponding dimension: for an octree in 3D, 010 means: negative
    direction in X, position direction in Y, negative direction in Z.

    \return the index of the adjacent node if it exists, nothing otherwise.
  */
  std::optional<Node_index> adjacent_node(Node_index n, const Local_coordinates& direction) const {

    // Direction:   LEFT  RIGHT  DOWN    UP  BACK FRONT
    // direction:    000    001   010   011   100   101

    // Nodes only have up to 2*dim different adjacent nodes (since boxes have 6 sides)
    CGAL_precondition(direction.to_ulong() < dimension * 2);

    // The root node has no adjacent nodes!
    if (is_root(n)) return {};

    // The least significant bit indicates the sign (which side of the node)
    bool sign = direction[0];

    // The first two bits indicate the dimension/axis (x, y, z)
    uint8_t dim = uint8_t((direction >> 1).to_ulong());

    // Create an offset so that the bit-significance lines up with the dimension (e.g., 1, 2, 4 --> 001, 010, 100)
    int8_t offset = (uint8_t) 1 << dim;

    // Finally, apply the sign to the offset
    offset = (sign ? offset : -offset);

    // Check if this child has the opposite sign along the direction's axis
    if (local_coordinates(n)[dim] != sign) {
      // This means the adjacent node is a direct sibling, the offset can be applied easily!
      return {child(parent(n), local_coordinates(n).to_ulong() + offset)};
    }

    // Find the parent's neighbor in that direction, if it exists
    auto adjacent_node_of_parent = adjacent_node(parent(n), direction);

    // If the parent has no neighbor, then this node doesn't have one
    if (!adjacent_node_of_parent) return {};

    // If the parent's adjacent node has no children, then it's this node's adjacent node
    if (is_leaf(*adjacent_node_of_parent))
      return adjacent_node_of_parent;

    // Return the nearest node of the parent by subtracting the offset instead of adding
    return {child(*adjacent_node_of_parent, local_coordinates(n).to_ulong() - offset)};
  }

  /*!
    \brief equivalent to `adjacent_node()`, with an adjacency direction rather than a bitset.

    \param n index of the node to find a neighbor of
    \param adjacency which way to find the adjacent node relative to this one
   */
  std::optional<Node_index> adjacent_node(Node_index n, Adjacency adjacency) const {
    return adjacent_node(n, std::bitset<dimension>(static_cast<int>(adjacency)));
  }

  /// @}

private: // functions :

  Node_index recursive_descendant(Node_index node, std::size_t i) { return child(node, i); }

  template <typename... Indices>
  Node_index recursive_descendant(Node_index node, std::size_t i, Indices... remaining_indices) {
    return recursive_descendant(child(node, i), remaining_indices...);
  }

  bool do_intersect(Node_index n, const Sphere& sphere) const {

    // Create a bounding box from the node
    Bbox node_box = bbox(n);

    // Check for intersection between the node and the sphere
    return CGAL::do_intersect(node_box, sphere);
  }

  template <typename Query, typename Node_output_iterator>
  Node_output_iterator intersected_nodes_recursive(const Query& query, Node_index node,
                                                   Node_output_iterator output) const {

    // Check if the current node intersects with the query
    if (CGAL::do_intersect(query, bbox(node))) {

      // if this node is a leaf, then it's considered an intersecting node
      if (is_leaf(node)) {
        *output++ = node;
        return output;
      }

      // Otherwise, each of the children need to be checked
      for (int i = 0; i < degree; ++i) {
        intersected_nodes_recursive(query, child(node, i), output);
      }
    }
    return output;
  }

  template <typename Result>
  auto nearest_k_neighbors_recursive(
    Sphere& search_bounds,
    Node_index node,
    std::vector<Result>& results,
    std::size_t k,
    FT epsilon = 0) const -> std::enable_if_t<supports_neighbor_search> {

    // Check whether the node has children
    if (is_leaf(node)) {

      // Base case: the node has no children

      // Loop through each of the elements contained by the node
      // Note: there might be none, and that should be fine!
      for (auto& e : data(node)) {

        // Pair that element with its distance from the search point
        Result current_element_with_distance =
        { e, m_traits.squared_distance_of_element_object()(e, m_traits.construct_center_d_object()(search_bounds)) };

        // Check if the new element is within the bounds
        if (current_element_with_distance.distance < m_traits.compute_squared_radius_d_object()(search_bounds)) {

          // Check if the results list is full
          if (results.size() == k) {
            // Delete a element if we need to make room
            results.pop_back();
          }

          // Add the new element
          results.push_back(current_element_with_distance);

          // Sort the list
          std::sort(results.begin(), results.end(), [=](auto& left, auto& right) {
            return left.distance < right.distance;
            });

          // Check if the results list is full
          if (results.size() == k) {

            // Set the search radius
            search_bounds = m_traits.construct_sphere_d_object()(m_traits.construct_center_d_object()(search_bounds), results.back().distance + epsilon);
          }
        }
      }
    }
    else {

      struct Node_index_with_distance {
        Node_index index;
        FT distance;

        Node_index_with_distance(const Node_index& index, FT distance) :
          index(index), distance(distance) {}
      };

      // Recursive case: the node has children

      // Create a list to map children to their distances
      std::vector<Node_index_with_distance> children_with_distances;
      children_with_distances.reserve(Self::degree);

      // Fill the list with child nodes
      for (int i = 0; i < Self::degree; ++i) {
        auto child_node = child(node, i);

        FT squared_distance = 0;
        Point c = m_traits.construct_center_d_object()(search_bounds);
        Point b = barycenter(child_node);
        for (const auto r : cartesian_range(c, b)) {
          FT d = (get<0>(r) - get<1>(r));
          squared_distance += d * d;
        }

        // Add a child to the list, with its distance
        children_with_distances.emplace_back(
          child_node,
          squared_distance
        );
      }

      // Sort the children by their distance from the search point
      std::sort(children_with_distances.begin(), children_with_distances.end(), [=](auto& left, auto& right) {
        return left.distance < right.distance;
        });

      // Loop over the children
      for (auto child_with_distance : children_with_distances) {

        // Check whether the bounding box of the child intersects with the search bounds
        if (CGAL::do_intersect(bbox(child_with_distance.index), search_bounds)) {

          // Recursively invoke this function
          nearest_k_neighbors_recursive(search_bounds, child_with_distance.index, results, k);
        }
      }
    }
  }

public:

  /// \cond SKIP_IN_MANUAL
  void dump_to_polylines(std::ostream& os) const {
    for (const auto n: traverse<Orthtrees::Preorder_traversal>())
      if (is_leaf(n)) {
        Bbox box = bbox(n);
        dump_box_to_polylines(box, os);
      }
  }

  void dump_box_to_polylines(const Bbox_2& box, std::ostream& os) const {
    // dump in 3D for visualization
    os << "5 "
       << box.xmin() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymax() << " 0 "
       << box.xmax() << " " << box.ymin() << " 0 "
       << box.xmin() << " " << box.ymin() << " 0" << std::endl;
  }

  void dump_box_to_polylines(const Bbox_3& box, std::ostream& os) const {
    // Back face
    os << "5 "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << std::endl;

    // Front face
    os << "5 "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmax() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmax() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmax() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << std::endl;

    // Traversal edges
    os << "2 "
       << box.xmin() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymin() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmin() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmin() << " " << box.ymax() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmax() << " " << box.ymin() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymin() << " " << box.zmax() << std::endl;
    os << "2 "
       << box.xmax() << " " << box.ymax() << " " << box.zmin() << " "
       << box.xmax() << " " << box.ymax() << " " << box.zmax() << std::endl;
  }

  std::string to_string(Node_index node) {
    std::stringstream stream;
    internal::print_orthtree_node(stream, node, *this);
    return stream.str();
  }

  friend std::ostream& operator<<(std::ostream& os, const Self& orthtree) {
    // Iterate over all nodes
    for (auto n: orthtree.traverse(Orthtrees::Preorder_traversal<Self>(orthtree))) {
      // Show the depth
      for (std::size_t i = 0; i < orthtree.depth(n); ++i)
        os << ". ";
      // Print the node
      internal::print_orthtree_node(os, n, orthtree);
      os << std::endl;
    }
    return os;
  }

  /// \endcond

}; // end class Orthtree

} // namespace CGAL

#endif // CGAL_ORTHTREE_H
