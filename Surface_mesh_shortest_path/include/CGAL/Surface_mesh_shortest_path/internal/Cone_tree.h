// Copyright (c) 2014 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_TREE_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_TREE_H

#include <CGAL/license/Surface_mesh_shortest_path.h>

#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Surface_mesh_shortest_path/internal/Cone_expansion_event.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Surface_mesh_shortest_paths_3 {
namespace internal {

template<class Traits>
class Cone_tree_node
{
public:
  enum Node_type
  {
    ROOT = 0,
    FACE_SOURCE = 1,
    EDGE_SOURCE = 2,
    VERTEX_SOURCE = 3,
    INTERVAL = 4
  };

private:
  typedef typename Traits::Triangle_mesh Triangle_mesh;
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Triangle_2 Triangle_2;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Ray_2 Ray_2;
  typedef typename boost::graph_traits<Triangle_mesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Surface_mesh_shortest_paths_3::internal::Cone_expansion_event<Traits> Cone_expansion_event;

private:
  // These could be pulled back into a 'context' class to save space
  const Traits& m_traits;
  const Triangle_mesh& m_graph;

  const halfedge_descriptor m_entryEdge;

  const Point_2 m_sourceImage;
  const Triangle_2 m_layoutFace;
  const FT m_pseudoSourceDistance;

  const Point_2 m_windowLeft;
  const Point_2 m_windowRight;

  std::size_t m_level;
  std::size_t m_treeId;

  const Node_type m_nodeType;

  Cone_tree_node* m_leftChild;
  std::vector<Cone_tree_node*> m_middleChildren;
  Cone_tree_node* m_rightChild;

  Cone_tree_node* m_parent;

public:
  Cone_expansion_event* m_pendingLeftSubtree;
  Cone_expansion_event* m_pendingRightSubtree;
  Cone_expansion_event* m_pendingMiddleSubtree;

private:
  void on_child_link(Cone_tree_node* child)
  {
    child->m_parent = this;
    child->m_level = m_level + 1;
    child->m_treeId = m_treeId;
  }

public:
  Cone_tree_node(const Traits& traits,
                 const Triangle_mesh& g,
                 const std::size_t treeId)
    : m_traits(traits)
    , m_graph(g)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_pseudoSourceDistance(0.0)
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_leftChild(nullptr)
    , m_rightChild(nullptr)
    , m_pendingLeftSubtree(nullptr)
    , m_pendingRightSubtree(nullptr)
    , m_pendingMiddleSubtree(nullptr)
  {
  }

  Cone_tree_node(const Traits& traits,
                 const Triangle_mesh& g,
                 const std::size_t treeId,
                 const halfedge_descriptor entryEdge)
    : m_traits(traits)
    , m_graph(g)
    , m_entryEdge(entryEdge)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_pseudoSourceDistance(0.0)
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_leftChild(nullptr)
    , m_rightChild(nullptr)
    , m_pendingLeftSubtree(nullptr)
    , m_pendingRightSubtree(nullptr)
    , m_pendingMiddleSubtree(nullptr)
  {
  }

  Cone_tree_node(const Traits& traits,
                 const Triangle_mesh& g,
                 const halfedge_descriptor entryEdge,
                 const Triangle_2& layoutFace,
                 const Point_2& sourceImage,
                 const FT& pseudoSourceDistance,
                 const Point_2& windowLeft,
                 const Point_2& windowRight,
                 const Node_type nodeType = INTERVAL)
    : m_traits(traits)
    , m_graph(g)
    , m_entryEdge(entryEdge)
    , m_sourceImage(sourceImage)
    , m_layoutFace(layoutFace)
    , m_pseudoSourceDistance(pseudoSourceDistance)
    , m_windowLeft(windowLeft)
    , m_windowRight(windowRight)
    , m_nodeType(nodeType)
    , m_leftChild(nullptr)
    , m_rightChild(nullptr)
    , m_pendingLeftSubtree(nullptr)
    , m_pendingRightSubtree(nullptr)
    , m_pendingMiddleSubtree(nullptr)
  {
  }

  std::size_t tree_id() const
  {
    return m_treeId;
  }

  std::size_t level() const
  {
    return m_level;
  }

  bool is_source_node() const
  {
    return m_nodeType == FACE_SOURCE || m_nodeType == EDGE_SOURCE || m_nodeType == VERTEX_SOURCE;
  }

  bool is_vertex_node() const
  {
    return m_nodeType == VERTEX_SOURCE;
  }

  bool is_root_node() const
  {
    return m_nodeType == ROOT;
  }

  const Triangle_2& layout_face() const
  {
    return m_layoutFace;
  }

  face_descriptor current_face() const
  {
    return face(m_entryEdge, m_graph);
  }

  bool is_null_face() const
  {
    return current_face() == Graph_traits::null_face();
  }

  std::size_t edge_face_index() const
  {
    return edge_index(entry_edge(), m_graph);
  }

  halfedge_descriptor entry_edge() const
  {
    return m_entryEdge;
  }

  halfedge_descriptor left_child_edge() const
  {
    return opposite(prev(m_entryEdge, m_graph), m_graph);
  }

  halfedge_descriptor right_child_edge() const
  {
    return opposite(next(m_entryEdge, m_graph), m_graph);
  }

  vertex_descriptor target_vertex() const
  {
    return target(next(m_entryEdge, m_graph), m_graph);
  }

  const Point_2& source_image() const
  {
    return m_sourceImage;
  }

  Node_type node_type() const
  {
    return m_nodeType;
  }

  FT distance_to_root(const Point_2& point) const
  {
    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());
    return CGAL::approximate_sqrt(csd2(point, m_sourceImage)) + m_pseudoSourceDistance;
  }

  FT distance_from_source_to_root() const
  {
    return m_pseudoSourceDistance;
  }

  FT distance_from_target_to_root() const
  {
    return distance_to_root(target_point());
  }

  Ray_2 left_boundary() const
  {
    return Ray_2(source_image(), m_windowLeft);
  }

  Ray_2 right_boundary() const
  {
    return Ray_2(source_image(), m_windowRight);
  }

  const Point_2& window_left() const
  {
    return m_windowLeft;
  }

  const Point_2& window_right() const
  {
    return m_windowRight;
  }

  Ray_2 ray_to_target_vertex() const
  {
    return Ray_2(source_image(), target_point());
  }

  bool inside_window(const Point_2& point) const
  {
    typename Traits::Orientation_2 orientation_2(m_traits.orientation_2_object());

    Point_2 sourceImagePoint(source_image());
    CGAL::Orientation leftOrientation = orientation_2(sourceImagePoint, m_windowLeft, point);
    CGAL::Orientation rightOrientation = orientation_2(sourceImagePoint, m_windowRight, point);

    return (leftOrientation == CGAL::RIGHT_TURN || leftOrientation == CGAL::COLLINEAR) &&
           (rightOrientation == CGAL::LEFT_TURN || rightOrientation == CGAL::COLLINEAR);
  }

  Point_2 target_point() const
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    return cv2(m_layoutFace, 2);
  }

  bool is_target_vertex_inside_window() const
  {
    return inside_window(target_point());
  }

  bool has_left_side() const
  {
    if (is_source_node())
    {
      return true;
    }

    return (m_traits.orientation_2_object()(source_image(), m_windowLeft, target_point()) != CGAL::LEFT_TURN);
  }

  bool has_right_side() const
  {
    return (m_traits.orientation_2_object()(source_image(), m_windowRight, target_point()) != CGAL::RIGHT_TURN);
  }

  Segment_2 left_child_base_segment() const
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(cv2(m_layoutFace, 0), cv2(m_layoutFace, 2));
  }

  Segment_2 right_child_base_segment() const
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(cv2(m_layoutFace, 2), cv2(m_layoutFace, 1));
  }

  Segment_2 entry_segment() const
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    return Segment_2(cv2(m_layoutFace, 0), cv2(m_layoutFace, 1));
  }

  bool has_middle_children() const
  {
    return m_middleChildren.size() > 0;
  }

  std::size_t num_middle_children() const
  {
    return m_middleChildren.size();
  }

  Cone_tree_node* get_middle_child(const std::size_t i) const
  {
    return m_middleChildren.at(i);
  }

  void push_middle_child(Cone_tree_node* child)
  {
    if (m_pendingMiddleSubtree != nullptr)
    {
      m_pendingMiddleSubtree->m_cancelled = true;
      m_pendingMiddleSubtree = nullptr;
    }

    m_middleChildren.push_back(child);
    on_child_link(child);
  }

  Cone_tree_node* pop_middle_child()
  {
    Cone_tree_node* temp = m_middleChildren.back();
    m_middleChildren.pop_back();
    return temp;
  }

  void set_left_child(Cone_tree_node* child)
  {
    if (m_pendingLeftSubtree != nullptr)
    {
      m_pendingLeftSubtree->m_cancelled = true;
      m_pendingLeftSubtree = nullptr;
    }

    m_leftChild = child;
    on_child_link(child);
  }

  Cone_tree_node* get_left_child() const
  {
    return m_leftChild;
  }

  Cone_tree_node* remove_left_child()
  {
    Cone_tree_node* temp = m_leftChild;
    m_leftChild = nullptr;
    return temp;
  }

  void set_right_child(Cone_tree_node* child)
  {
    if (m_pendingRightSubtree != nullptr)
    {
      m_pendingRightSubtree->m_cancelled = true;
      m_pendingRightSubtree = nullptr;
    }

    m_rightChild = child;
    on_child_link(child);
  }

  Cone_tree_node* get_right_child() const
  {
    return m_rightChild;
  }

  Cone_tree_node* remove_right_child()
  {
    Cone_tree_node* temp = m_rightChild;
    m_rightChild = nullptr;
    return temp;
  }

  Cone_tree_node* parent() const
  {
    return m_parent;
  }

  bool is_left_child() const
  {
    return m_parent != nullptr && m_parent->m_leftChild == this;
  }

  bool is_right_child() const
  {
    return m_parent != nullptr && m_parent->m_rightChild == this;
  }
};

} // namespace internal
} // namespace Surface_mesh_shortest_paths_3
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_INTERNAL_CONE_TREE_H
