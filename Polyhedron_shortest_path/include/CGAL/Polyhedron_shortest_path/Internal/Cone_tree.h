// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_shortest_path/Internal/Cone_expansion_event.h>

#include <vector>

namespace CGAL
{

namespace internal
{

// Computes the conditions of Theorem 3.2 of "Improving Chen and Han's Algorithm on the Discrete Geodesic Problem" 
// Shi-Qing Xin and Guo-Jin Wang
// ACM Transactions on Graph Algorithms, Volume 28, No. 4, Article 104 (August 2009)
// http://doi.acm.org/10.1145/1559755.1559761
template <class Traits>
class Xin_wang_window_distance_filter
{
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::FT FT;
  typedef typename Traits::Compute_squared_distance_2 Compute_squared_distance_2;
  
  Compute_squared_distance_2 compute_squared_distance_2;
  
public:
  Xin_wang_window_distance_filter()
  {
  }
  
  Xin_wang_window_distance_filter(const Compute_squared_distance_2& cds)
    : compute_squared_distance_2(cds)
  {
  }

  bool operator()(
    Point_2 v1,                         // Near vertex
    FT d1,                      // Near vertex distance 
    Point_2 v2,                         // Far vertex
    FT d2,                      // Far vertex distance 
    Point_2 v3,                         // Opposite vertex
    FT d3,                      // Opposite vertex distance 
    Point_2 A,                          // Near intersection
    Point_2 B,                          // Far intersection
    Point_2 I,                          // Source image
    FT d                        // Distance to source image 
  )    
  {
    if (d + CGAL::sqrt(compute_squared_distance_2(I, B)) > d1 + CGAL::sqrt(compute_squared_distance_2(v1, B)))
    {
      return false;
    }
    
    if (d + CGAL::sqrt(compute_squared_distance_2(I, A)) > d2 + CGAL::sqrt(compute_squared_distance_2(v2, A)))
    {
      return false;
    }
    
    if (d + CGAL::sqrt(compute_squared_distance_2(I, A)) > d3 + CGAL::sqrt(compute_squared_distance_2(v3, A)))
    {
      return false;
    }
    
    return true;
  }
};

template<class Traits>
class Cone_tree_node
{
public:
  enum Node_type
  {
    ROOT,
    FACE_SOURCE,
    EDGE_SOURCE,
    VERTEX_SOURCE,
    INTERVAL,
  };

private:
  typedef typename Traits::Polyhedron Polyhedron;
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Triangle_2 Triangle_2;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Ray_2 Ray_2;
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::Orientation_2 Orientation_2;
  typedef typename Traits::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef typename Traits::Construct_triangle_location_2 Construct_triangle_location_2;
  typedef typename CGAL::internal::Cone_expansion_event<Traits> Cone_expansion_event;

private:
  Orientation_2 m_orientation_2;
  Compute_squared_distance_2 m_compute_squared_distance_2;
  Construct_triangle_location_2 m_construct_triangle_location_2;
  
  Polyhedron* m_polyhedron;
  halfedge_descriptor m_entryEdge;
  
  Triangle_2 m_layoutFace;
  Cone_tree_node* m_leftChild;
  Cone_tree_node* m_rightChild;
  std::vector<Cone_tree_node*> m_middleChildren;
  Cone_tree_node* m_parent;
  
  Point_2 m_sourceImage;
  FT m_pseudoSourceDistance;

  Point_2 m_windowLeft;
  Point_2 m_windowRight;

  Node_type m_nodeType;
  size_t m_treeId;
  size_t m_level;
  
private:
  void on_child_link(Cone_tree_node* child)
  {
    child->m_parent = this;
    child->m_level = m_level + 1;
    child->m_treeId = m_treeId;
  }
  
public:
  Cone_tree_node(Polyhedron* polyhedron, size_t treeId)
    : m_polyhedron(polyhedron)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
  {
  }
  
  Cone_tree_node(Polyhedron* polyhedron, size_t treeId, halfedge_descriptor entryEdge)
    : m_polyhedron(polyhedron)
    , m_entryEdge(entryEdge)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
  {
  }

  Cone_tree_node(Polyhedron* polyhedron, halfedge_descriptor entryEdge, const Triangle_2& layoutFace, const Point_2& sourceImage, const FT& pseudoSourceDistance, const Point_2& windowLeft, const Point_2& windowRight, Node_type nodeType = INTERVAL)
    : m_polyhedron(polyhedron)
    , m_entryEdge(entryEdge)
    , m_layoutFace(layoutFace)
    , m_sourceImage(sourceImage)
    , m_pseudoSourceDistance(pseudoSourceDistance)
    , m_windowLeft(windowLeft)
    , m_windowRight(windowRight)
    , m_nodeType(nodeType)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
  {
  }

  size_t tree_id()
  {
    return m_treeId;
  }
  
  size_t level()
  {
    return m_level;
  }
  
  bool is_source_node()
  {
    return m_nodeType == FACE_SOURCE || m_nodeType == EDGE_SOURCE || m_nodeType == VERTEX_SOURCE;
  }
  
  bool is_vertex_node()
  {
    return m_nodeType == VERTEX_SOURCE;
  }
  
  bool is_root_node()
  {
    return m_nodeType == ROOT;
  }
  
  Triangle_2 layout_face()
  {
    return m_layoutFace;
  }
  
  face_descriptor current_face()
  {
    return CGAL::face(m_entryEdge, m_polyhedron);
  }
  
  halfedge_descriptor entry_edge()
  {
    return m_entryEdge;
  }
  
  halfedge_descriptor left_child_edge()
  {
    return CGAL::opposite(CGAL::prev(m_entryEdge, *m_polyhedron), *m_polyhedron);
  }
  
  halfedge_descriptor right_child_edge()
  {
    return CGAL::opposite(CGAL::next(m_entryEdge, *m_polyhedron), *m_polyhedron);
  }
  
  vertex_descriptor target_vertex()
  {
    return CGAL::target(CGAL::next(m_entryEdge, *m_polyhedron), *m_polyhedron);
  }
  
  Point_2 source_image()
  {
    return m_sourceImage;
  }
  
  Node_type node_type()
  {
    return m_nodeType;
  }
  
  FT distance_to_root(const Point_2& point)
  {
    return CGAL::sqrt(m_compute_squared_distance_2(point, m_sourceImage)) + m_pseudoSourceDistance;
  }
  
  FT distance_from_source_to_root()
  {
    return m_pseudoSourceDistance;
  }
  
  FT distance_from_target_to_root()
  {
    return distance_to_root(target_vertex_location());
  }
  
  Ray_2 left_boundary()
  {
    return Ray_2(source_image(), m_windowLeft);
  }
  
  Ray_2 right_boundary()
  {
    return Ray_2(source_image(), m_windowRight);
  }
  
  Point_2 window_left()
  {
    return m_windowLeft;
  }
  
  Point_2 window_right()
  {
    return m_windowRight;
  }
  
  Ray_2 ray_to_target_vertex()
  {
    return Ray_2(source_image(), target_vertex_location());
  }
  
  bool inside_window(const Point_2& point)
  {
    Point_2 sourceImagePoint(source_image());
    return m_orientation_2(sourceImagePoint, m_windowLeft, point) == CGAL::RIGHT_TURN && m_orientation_2(sourceImagePoint, m_windowRight, point) == CGAL::LEFT_TURN;
  }

  Point_2 target_vertex_location()
  {
    return m_layoutFace[2];
  }
  
  bool is_target_vertex_inside_window()
  {
    return inside_window(target_vertex_location());
  }
  
  bool has_left_side()
  {
    return m_orientation_2(source_image(), m_windowLeft, target_vertex_location()) == CGAL::RIGHT_TURN;
  }
  
  bool has_right_side()
  {
    return m_orientation_2(source_image(), m_windowRight, target_vertex_location()) == CGAL::LEFT_TURN;
  }
  
  Segment_2 left_child_base_segment()
  {
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(m_layoutFace[0], m_layoutFace[2]);
  }
  
  Segment_2 right_child_base_segment()
  {
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(m_layoutFace[2], m_layoutFace[1]);
  }
  
  Segment_2 entry_segment()
  {
    return Segment_2(m_layoutFace[0], m_layoutFace[1]);
  }
  
  bool has_middle_children()
  {
    return m_middleChildren.size() > 0;
  }
  
  size_t num_middle_children()
  {
    return m_middleChildren.size();
  }
  
  Cone_tree_node* get_middle_child(size_t i)
  {
    return m_middleChildren.at(i);
  }

  void push_middle_child(Cone_tree_node* child)
  {
    if (m_pendingMiddleSubtree != NULL)
    {
      m_pendingMiddleSubtree->m_cancelled = true;
      m_pendingMiddleSubtree = NULL;
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
    if (m_pendingLeftSubtree != NULL)
    {
      m_pendingLeftSubtree->m_cancelled = true;
      m_pendingLeftSubtree = NULL;
    }
    
    m_leftChild = child;
    on_child_link(child);
  }
  
  Cone_tree_node* get_left_child()
  {
    return m_leftChild;
  }
  
  Cone_tree_node* remove_left_child()
  {
    Cone_tree_node* temp = m_leftChild;
    m_leftChild = NULL;
    return temp;
  }
  
  void set_right_child(Cone_tree_node* child)
  {
    if (m_pendingRightSubtree != NULL)
    {
      m_pendingRightSubtree->m_cancelled = true;
      m_pendingRightSubtree = NULL;
    }
  
    m_rightChild = child;
    on_child_link(child);
  }
  
  Cone_tree_node* get_right_child()
  {
    return m_rightChild;
  }
  
  Cone_tree_node* remove_right_child()
  {
    Cone_tree_node* temp = m_rightChild;
    m_rightChild = NULL;
    return temp;
  }
  
public:
  Cone_expansion_event* m_pendingLeftSubtree;
  Cone_expansion_event* m_pendingRightSubtree;
  Cone_expansion_event* m_pendingMiddleSubtree;
  
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H
