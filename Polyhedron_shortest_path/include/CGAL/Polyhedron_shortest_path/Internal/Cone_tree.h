// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H

#include <vector>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Polyhedron_shortest_path/Internal/Cone_expansion_event.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>

#if !defined(NDEBUG)

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Interval_nt.h>

#endif 

namespace CGAL
{

namespace internal
{

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
    INTERVAL = 4,
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
  
#if !defined(NDEBUG)
  typedef CGAL::Simple_cartesian<Interval_nt<true> > IntervalKernel;
#endif

private:
  // These could be pulled back into a 'context' class to save space
  Orientation_2 m_orientation_2;
  Compute_squared_distance_2 m_compute_squared_distance_2;
  Construct_triangle_location_2 m_construct_triangle_location_2;
  Polyhedron& m_polyhedron;
  
  halfedge_descriptor m_entryEdge;
  
  Point_2 m_sourceImage;
  Triangle_2 m_layoutFace;
  FT m_pseudoSourceDistance;

  Point_2 m_windowLeft;
  Point_2 m_windowRight;

  size_t m_level;
  size_t m_treeId;

  Node_type m_nodeType;
  
  Cone_tree_node* m_leftChild;
  Cone_tree_node* m_rightChild;
  
  std::vector<Cone_tree_node*> m_middleChildren;
  Cone_tree_node* m_parent;

private:
  void on_child_link(Cone_tree_node* child)
  {
    child->m_parent = this;
    child->m_level = m_level + 1;
    child->m_treeId = m_treeId;
  }
  
public:
  Cone_tree_node(Traits& traits, Polyhedron& polyhedron, size_t treeId)
    : m_orientation_2(traits.orientation_2_object())
    , m_compute_squared_distance_2(traits.compute_squared_distance_2_object())
    , m_construct_triangle_location_2(traits.construct_triangle_location_2_object())
    , m_polyhedron(polyhedron)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_pseudoSourceDistance(0.0)
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
  {
  }
  
  Cone_tree_node(Traits& traits, Polyhedron& polyhedron, size_t treeId, halfedge_descriptor entryEdge)
    : m_orientation_2(traits.orientation_2_object())
    , m_compute_squared_distance_2(traits.compute_squared_distance_2_object())
    , m_construct_triangle_location_2(traits.construct_triangle_location_2_object())
    , m_polyhedron(polyhedron)
    , m_entryEdge(entryEdge)
    , m_sourceImage(Point_2(CGAL::ORIGIN))
    , m_layoutFace(Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN),Point_2(CGAL::ORIGIN))
    , m_pseudoSourceDistance(0.0)
    , m_level(0)
    , m_treeId(treeId)
    , m_nodeType(ROOT)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
  {
  }

  Cone_tree_node(Traits& traits, Polyhedron& polyhedron, halfedge_descriptor entryEdge, const Triangle_2& layoutFace, const Point_2& sourceImage, const FT& pseudoSourceDistance, const Point_2& windowLeft, const Point_2& windowRight, Node_type nodeType = INTERVAL)
    : m_orientation_2(traits.orientation_2_object())
    , m_compute_squared_distance_2(traits.compute_squared_distance_2_object())
    , m_construct_triangle_location_2(traits.construct_triangle_location_2_object())
    , m_polyhedron(polyhedron)
    , m_entryEdge(entryEdge)
    , m_sourceImage(sourceImage)
    , m_layoutFace(layoutFace)
    , m_pseudoSourceDistance(pseudoSourceDistance)
    , m_windowLeft(windowLeft)
    , m_windowRight(windowRight)
    , m_nodeType(nodeType)
    , m_leftChild(NULL)
    , m_rightChild(NULL)
    , m_pendingLeftSubtree(NULL)
    , m_pendingRightSubtree(NULL)
    , m_pendingMiddleSubtree(NULL)
  {
  }

  size_t tree_id() const
  {
    return m_treeId;
  }
  
  size_t level() const
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
  
  Triangle_2 layout_face() const
  {
    return m_layoutFace;
  }
  
  face_descriptor current_face() const
  {
    return CGAL::face(m_entryEdge, m_polyhedron);
  }
  
  bool is_null_face() const
  {
    return current_face() == GraphTraits::null_face();
  }
  
  size_t edge_face_index() const
  {
    return CGAL::internal::edge_index(entry_edge(), m_polyhedron);
  }
  
  halfedge_descriptor entry_edge() const
  {
    return m_entryEdge;
  }
  
  halfedge_descriptor left_child_edge() const
  {
    return CGAL::opposite(CGAL::prev(m_entryEdge, m_polyhedron), m_polyhedron);
  }
  
  halfedge_descriptor right_child_edge() const
  {
    return CGAL::opposite(CGAL::next(m_entryEdge, m_polyhedron), m_polyhedron);
  }
  
  vertex_descriptor target_vertex() const
  {
    return CGAL::target(CGAL::next(m_entryEdge, m_polyhedron), m_polyhedron);
  }
  
  Point_2 source_image() const
  {
    return m_sourceImage;
  }
  
  Node_type node_type() const
  {
    return m_nodeType;
  }
  
  FT distance_to_root(const Point_2& point) const
  {
    return CGAL::sqrt(m_compute_squared_distance_2(point, m_sourceImage)) + m_pseudoSourceDistance;
  }
  
  FT distance_from_source_to_root() const
  {
    return m_pseudoSourceDistance;
  }
  
  FT distance_from_target_to_root() const
  {
    return distance_to_root(target_vertex_location());
  }
  
#if !defined(NDEBUG)

  Interval_nt<true> distance_to_root_interval(const Point_2& point) const
  {
    if (is_root_node())
    {
      return Interval_nt<true>(0.0);
    }
    
    IntervalKernel::Compute_squared_distance_2 compute_squared_distance_2;
    IntervalKernel::Point_2 iPoint(point.x(), point.y());
    IntervalKernel::Point_2 iSource(m_sourceImage.x(), m_sourceImage.y());
    
    return CGAL::sqrt(compute_squared_distance_2(iPoint, iSource)) + distance_from_source_to_root_interval();
  }
  
  Interval_nt<true> distance_from_source_to_root_interval() const
  {
    if (is_root_node())
    {
      return Interval_nt<true>(0.0);
    }
    else if (is_source_node())
    {
      return m_parent->distance_from_target_to_root_interval();
    }
    else
    {
      return m_parent->distance_from_source_to_root_interval();
    }
  }
  
  Interval_nt<true> distance_from_target_to_root_interval() const
  {
    if (is_root_node())
    {
      return Interval_nt<true>(0.0);
    }
    else
    {
      return distance_to_root_interval(target_vertex_location());
    }
  }

#endif
  
  Ray_2 left_boundary() const
  {
    return Ray_2(source_image(), m_windowLeft);
  }
  
  Ray_2 right_boundary() const
  {
    return Ray_2(source_image(), m_windowRight);
  }
  
  Point_2 window_left() const
  {
    return m_windowLeft;
  }
  
  Point_2 window_right() const
  {
    return m_windowRight;
  }
  
  Ray_2 ray_to_target_vertex() const
  {
    return Ray_2(source_image(), target_vertex_location());
  }
  
  bool inside_window(const Point_2& point) const
  {
    Point_2 sourceImagePoint(source_image());
    return m_orientation_2(sourceImagePoint, m_windowLeft, point) == CGAL::RIGHT_TURN && m_orientation_2(sourceImagePoint, m_windowRight, point) == CGAL::LEFT_TURN;
  }

  Point_2 target_vertex_location() const
  {
    return m_layoutFace[2];
  }
  
  bool is_target_vertex_inside_window() const
  {
    return inside_window(target_vertex_location());
  }
  
  bool has_left_side() const
  {
    if (is_source_node())
    {
      return true;
    }
    //Ray_2 leftBoundary = left_boundary();
    //std::cout << "Left side check: " << m_orientation_2(source_image(), m_windowLeft, target_vertex_location()) << " vs. " << m_orientation_2(leftBoundary.source(), leftBoundary.point(1), target_vertex_location()) << std::endl;
    
    return m_orientation_2(source_image(), m_windowLeft, target_vertex_location()) == CGAL::RIGHT_TURN;
  }
  
  bool has_right_side() const
  {
    //Ray_2 rightBoundary = right_boundary();
    //std::cout << "Right side check: " << m_orientation_2(source_image(), m_windowRight, target_vertex_location()) << " vs. " << m_orientation_2(rightBoundary.source(), rightBoundary.point(1), target_vertex_location()) << std::endl;
    
    return m_orientation_2(source_image(), m_windowRight, target_vertex_location()) == CGAL::LEFT_TURN;
  }
  
  Segment_2 left_child_base_segment() const
  {
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(m_layoutFace[0], m_layoutFace[2]);
  }
  
  Segment_2 right_child_base_segment() const
  {
    // reversed to maintain consistent triangle winding on the child
    return Segment_2(m_layoutFace[2], m_layoutFace[1]);
  }
  
  Segment_2 entry_segment() const
  {
    return Segment_2(m_layoutFace[0], m_layoutFace[1]);
  }
   
  bool has_middle_children() const
  {
    return m_middleChildren.size() > 0;
  }
  
  size_t num_middle_children() const
  {
    return m_middleChildren.size();
  }
  
  Cone_tree_node* get_middle_child(size_t i) const
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
  
  Cone_tree_node* get_left_child() const
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
  
  Cone_tree_node* get_right_child() const
  {
    return m_rightChild;
  }
  
  Cone_tree_node* remove_right_child()
  {
    Cone_tree_node* temp = m_rightChild;
    m_rightChild = NULL;
    return temp;
  }
  
  Cone_tree_node* parent() const
  {
    return m_parent;
  }
  
  bool is_left_child() const
  {
    return m_parent != NULL && m_parent->m_leftChild == this;
  }
  
  bool is_right_child() const
  {
    return m_parent != NULL && m_parent->m_rightChild == this;
  }
  
public:
  Cone_expansion_event* m_pendingLeftSubtree;
  Cone_expansion_event* m_pendingRightSubtree;
  Cone_expansion_event* m_pendingMiddleSubtree;
  
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_CONE_TREE_H
