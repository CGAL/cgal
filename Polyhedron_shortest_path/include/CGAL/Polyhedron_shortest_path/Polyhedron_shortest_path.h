// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <map>
#include <utility>
#include <queue>
#include <boost/array.hpp>
#include <CGAL/Polyhedron_shortest_path/Internal/Cone_tree.h>
#include <CGAL/Polyhedron_shortest_path/Internal/Barycentric.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {


template<class Traits, 
  class VIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_external_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type,
  class VPM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_point_t>::type>
class Polyhedron_shortest_path
{
public:
  typedef VIM VertexIndexMap;
  typedef HIM HalfedgeIndexMap;
  typedef FIM FaceIndexMap;
  typedef VPM VertexPointMap;

  typedef typename Traits::Polyhedron Polyhedron;
  typedef typename Traits::Triangle_3 Triangle_3;
  typedef typename Traits::Triangle_2 Triangle_2;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Ray_2 Ray_2;
  typedef typename Traits::Point_3 Point_3;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Vector_2 Vector_2;

  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::halfedge_iterator halfedge_iterator;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;
  typedef typename Traits::FT FT;
  
  typedef typename internal::Cone_tree_node<Traits> Cone_tree_node;
  typedef typename internal::Cone_expansion_event<Traits> Cone_expansion_event;
  typedef typename Traits::Intersect_2 Intersect_2;

  typedef typename std::priority_queue<Cone_expansion_event, std::vector<Cone_expansion_event*>, internal::Cone_expansion_event_min_priority_queue_comparator<Traits> > Expansion_priqueue;
  
  typedef typename std::pair<face_descriptor, Barycentric_coordinate> FaceLocationPair;
  typedef typename std::pair<Cone_tree_node*, FT> NodeDistancePair;

private:

  template <class Visitor>
  struct Point_path_visitor_wrapper
  {
    Visitor& m_visitor;
    Polyhedron& m_polyhedron;
    VertexPointMap& m_vertexPointMap;
    
    Point_path_visitor_wrapper(Visitor& visitor, Polyhedron& polyhedron, VertexPointMap& vertexPointMap)
      : m_visitor(visitor)
      , m_polyhedron(polyhedron)
      , m_vertexPointMap(vertexPointMap)
    {
    }
    
    void edge(halfedge_descriptor he, FT alpha)
    {
      return m_visitor.point(CGAL::internal::interpolate_points<Point_3, Vector_3, FT>(m_vertexPointMap[CGAL::source(he, m_polyhedron)], m_vertexPointMap[CGAL::target(he, m_polyhedron)], alpha));
    }
    
    void vertex(vertex_descriptor v)
    {
      return m_vertexPointMap[v];
    }
    
    void face(face_descriptor face, Barycentric_coordinate alpha)
    {
      return m_traits.construct_barycentric_coordinate_3_object()(triangle_from_halfedge(CGAL::halfedge(face, m_polyhedron)), alpha);
    }
  };

private:
  Traits m_traits;

  std::vector<bool> m_vertexIsPseudoSource;
  
  std::vector<NodeDistancePair> m_vertexOccupiers;
  std::vector<NodeDistancePair> m_closestToVertices;
  
  std::vector<Cone_tree_node*> m_rootNodes;
  std::vector<FaceLocationPair> m_faceLocations;
  
  std::vector<std::vector<Cone_tree_node*> > m_faceOccupiers;
  
  Expansion_priqueue m_expansionPriqueue;

  VertexIndexMap m_vertexIndexMap;
  FaceIndexMap m_faceIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  VertexPointMap m_vertexPointMap;
  
  Polyhedron& m_polyhedron;
  
public:

  bool m_debugOutput;
  
private:

  Triangle_3 triangle_from_halfedge(halfedge_descriptor edge)
  {
    halfedge_descriptor start = edge;
    halfedge_descriptor current = edge;
    
    Point_3 points[3];
    size_t currentPoint = 0;
    
    do
    {
      points[currentPoint] = m_vertexPointMap[boost::source(current, m_polyhedron)];
      current = CGAL::next(current, m_polyhedron);
      ++currentPoint;
    }
    while (current != start);
    
    return Triangle_3(points[0], points[1], points[2]);
  }

  bool window_distance_filter(Cone_tree_node* cone, Segment_2 windowSegment, bool reversed)
  {
    
    Segment_2 parentEntrySegment = cone->entry_segment();
    Point_2 v2 = cone->target_vertex_location();
    Point_2 I = cone->source_image();
    FT d = cone->distance_from_source_to_root();
    FT d1;
    FT d2;
    FT d3;
    Point_2 A;
    Point_2 B;
    Point_2 v1;
    Point_2 v3;
    
    size_t v1Index = m_vertexIndexMap[CGAL::source(cone->entry_edge(), m_polyhedron)];
    size_t v2Index = m_vertexIndexMap[cone->target_vertex()];
    size_t v3Index = m_vertexIndexMap[CGAL::target(cone->entry_edge(), m_polyhedron)];
    
    NodeDistancePair v1Distance = m_closestToVertices[v1Index];
    NodeDistancePair v2Distance = m_closestToVertices[v2Index];
    NodeDistancePair v3Distance = m_closestToVertices[v3Index];
    
    if (reversed)
    {
      std::swap(v1Distance, v3Distance);
      A = windowSegment[1];
      B = windowSegment[0];
      v1 = parentEntrySegment[1];
      v3 = parentEntrySegment[0];
    }
    else
    {
      A = windowSegment[0];
      B = windowSegment[1];
      v1 = parentEntrySegment[0];
      v3 = parentEntrySegment[1];
    }
    
    d1 = v1Distance.second;
    d2 = v2Distance.second;
    d3 = v3Distance.second;
    
    bool hasD1 = v1Distance.first != NULL;
    bool hasD2 = v2Distance.first != NULL;
    bool hasD3 = v3Distance.first != NULL;
    
    if (hasD1 && (d + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(I, B)) > d1 + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(v1, B))))
    {
      return false;
    }
    
    if (hasD2 && (d + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(I, A)) > d2 + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(v2, A))))
    {
      return false;
    }
    
    if (hasD3 && (d + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(I, A)) > d3 + CGAL::sqrt(m_traits.compute_squared_distance_2_object()(v3, A))))
    {
      return false;
    }
    
    return true;
  }
    
  void expand_left_child(Cone_tree_node* cone, Segment_2 windowSegment)
  {
    assert(cone->m_pendingLeftSubtree != NULL);
    
    cone->m_pendingLeftSubtree = NULL;
    
    if (window_distance_filter(cone, windowSegment, false))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->left_child_edge());
      Triangle_2 layoutFace = m_traits.flatten_triangle_3_along_segment_2_object()(adjacentFace, 0, cone->left_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, cone->left_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), windowSegment[0], windowSegment[1], Cone_tree_node::INTERVAL);
      cone->set_left_child(child);
      process_node(child);
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNode was filtered." << std::endl;
    }
  }
  
  void expand_right_child(Cone_tree_node* cone, Segment_2 windowSegment)
  {
    assert(cone->m_pendingRightSubtree != NULL);
    
    cone->m_pendingRightSubtree = NULL;
    
    if (window_distance_filter(cone, windowSegment, true))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->right_child_edge());
      Triangle_2 layoutFace = m_traits.flatten_triangle_3_along_segment_2_object()(adjacentFace, 0, cone->right_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, cone->right_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), windowSegment[0], windowSegment[1], Cone_tree_node::INTERVAL);
      cone->set_right_child(child);
      process_node(child);
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNode was filtered." << std::endl;
    }
  }
  
  void expand_root(face_descriptor face, Barycentric_coordinate location)
  {
    size_t associatedEdge;
    CGAL::internal::Barycentric_coordinate_type type = classify_barycentric_coordinate(location, associatedEdge);
    
    switch (type)
    {
      case CGAL::internal::BARYCENTRIC_COORDINATE_INTERNAL:
        expand_face_root(face, location);
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_EDGE:
        {
          halfedge_descriptor he = CGAL::halfedge(face, m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            he = CGAL::next(he, m_polyhedron);
          }
          expand_edge_root(he, location[associatedEdge], location[(associatedEdge + 1) % 3]);
        }
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX:
        {
          halfedge_descriptor he = CGAL::halfedge(face, m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            he = CGAL::next(he, m_polyhedron);
          }
          expand_vertex_root(CGAL::source(he, m_polyhedron));
        }
        break;
      default:
        assert(false && "Invalid face location");
        // Perhaps hit an assertion that the type must not be external or invalid?
    }
  }
  
  void expand_face_root(face_descriptor faceId, Barycentric_coordinate faceLocation)
  {
    halfedge_descriptor start = CGAL::halfedge(faceId, m_polyhedron);
    halfedge_descriptor current = start;
    
    Cone_tree_node* faceRoot = new Cone_tree_node(m_polyhedron, m_rootNodes.size());
    m_rootNodes.push_back(faceRoot);
    
    if (m_debugOutput)
    {
      std::cout << "\tFace Root Expansion: face = " << m_faceIndexMap[faceId] << " , Location = " << faceLocation << std::endl;
    }
    
    for (size_t currentVertex = 0; currentVertex < 3; ++currentVertex)
    {
      Triangle_3 face3d(triangle_from_halfedge(current));
      Triangle_2 layoutFace(m_traits.project_triangle_3_to_triangle_2_object()(face3d));
      Barycentric_coordinate rotatedFaceLocation(faceLocation[currentVertex], faceLocation[(currentVertex + 1) % 3], faceLocation[(currentVertex + 2) % 3]);
      Point_2 sourcePoint(m_traits.construct_triangle_location_2_object()(layoutFace, rotatedFaceLocation));
      
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, current, layoutFace, sourcePoint, FT(0.0), layoutFace[0], layoutFace[2], Cone_tree_node::FACE_SOURCE);
      faceRoot->push_middle_child(child);
      
      if (m_debugOutput)
      {
        std::cout << "\tExpanding face root #" << currentVertex << " : " << std::endl;;
        std::cout << "\t\tFace = " << layoutFace << std::endl;
        std::cout << "\t\tLocation = " << sourcePoint << std::endl;
      }
      
      process_node(child);

      current = CGAL::next(current, m_polyhedron);
    }
  }

  void expand_edge_root(halfedge_descriptor baseEdge, FT t0, FT t1)
  {
    if (m_debugOutput)
    {
      std::cout << "\tEdge Root Expansion: faceA = " << m_faceIndexMap[CGAL::face(baseEdge, m_polyhedron)] << " , faceB = " << m_faceIndexMap[CGAL::face(CGAL::opposite(baseEdge, m_polyhedron), m_polyhedron)] << " , t0 = " << t0 << " , t1 = " << t1 << std::endl;
    }
    
    halfedge_descriptor baseEdges[2];
    baseEdges[0] = baseEdge;
    baseEdges[1] = CGAL::opposite(baseEdge, m_polyhedron);
    
    Triangle_3 faces3d[2];
    Triangle_2 layoutFaces[2];

    for (size_t i = 0; i < 2; ++i)
    {
       faces3d[i] = triangle_from_halfedge(baseEdges[i]);
       layoutFaces[i] = m_traits.project_triangle_3_to_triangle_2_object()(faces3d[i]);
    }
    
    Point_2 sourcePoints[2];
    sourcePoints[0] = Point_2(layoutFaces[0][0][0] * t0 + layoutFaces[0][1][0] * t1, layoutFaces[0][0][1] * t0 + layoutFaces[0][1][1] * t1); 
    sourcePoints[1] = Point_2(layoutFaces[1][0][0] * t0 + layoutFaces[1][1][0] * t1, layoutFaces[1][0][1] * t0 + layoutFaces[1][1][1] * t1); 
    
    Cone_tree_node* edgeRoot = new Cone_tree_node(m_polyhedron, m_rootNodes.size());
    m_rootNodes.push_back(edgeRoot);
    
    for (size_t side = 0; side < 2; ++side)
    {
      if (m_debugOutput)
      {
        std::cout << "\tExpanding edge root #" << side << " : " << std::endl;;
        std::cout << "\t\tFace = " << layoutFaces[side] << std::endl;
        std::cout << "\t\tLocation = " << sourcePoints[side] << std::endl;
      }
      
      Cone_tree_node* mainChild = new Cone_tree_node(m_polyhedron, baseEdges[side], layoutFaces[side], sourcePoints[side], FT(0.0), layoutFaces[side][0], layoutFaces[side][2], Cone_tree_node::EDGE_SOURCE);
      edgeRoot->push_middle_child(mainChild);
      process_node(mainChild);

      Cone_tree_node* oppositeChild = new Cone_tree_node(m_polyhedron, baseEdges[side], Triangle_2(layoutFaces[side][2], layoutFaces[side][1], layoutFaces[side][2]), sourcePoints[side], FT(0.0), layoutFaces[side][1], layoutFaces[side][2], Cone_tree_node::EDGE_SOURCE);
      edgeRoot->push_middle_child(oppositeChild);
      process_node(oppositeChild);
    }
  }

  void expand_vertex_root(vertex_descriptor vertex)
  {
    if (m_debugOutput)
    {
      std::cout << "\tVertex Root Expansion: Vertex = " << m_vertexIndexMap[vertex] << std::endl;
    }
    
    Cone_tree_node* vertexRoot = new Cone_tree_node(m_polyhedron, m_rootNodes.size(), CGAL::prev(CGAL::halfedge(vertex, m_polyhedron), m_polyhedron));
    m_rootNodes.push_back(vertexRoot);
    
    m_closestToVertices[m_vertexIndexMap[vertex]] = NodeDistancePair(vertexRoot, FT(0.0));
    
    expand_psuedo_source(vertexRoot);
  }

  void expand_psuedo_source(Cone_tree_node* parent)
  {
    parent->m_pendingMiddleSubtree = NULL;
    
    vertex_descriptor expansionVertex = parent->target_vertex();
  
    halfedge_descriptor startEdge = CGAL::halfedge(expansionVertex, m_polyhedron);
    halfedge_descriptor currentEdge = CGAL::halfedge(expansionVertex, m_polyhedron);
        
    FT distanceFromTargetToRoot = parent->distance_from_target_to_root();
      
    if (m_debugOutput)
    {
      std::cout << "Distance from target to root: " << distanceFromTargetToRoot << std::endl;
    }
    
    do
    {
      Triangle_3 face3d(triangle_from_halfedge(currentEdge));
      Triangle_2 layoutFace(m_traits.project_triangle_3_to_triangle_2_object()(face3d));
      
      if (m_debugOutput)
      {
        std::cout << "Expanding PsuedoSource: id = " << m_faceIndexMap[CGAL::face(currentEdge, m_polyhedron)] << " , face = " << layoutFace << std::endl;
      }
      
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, currentEdge, layoutFace, layoutFace[1], distanceFromTargetToRoot, layoutFace[0], layoutFace[2], Cone_tree_node::VERTEX_SOURCE);
      parent->push_middle_child(child);
      process_node(child);
      
      currentEdge = CGAL::opposite(CGAL::next(currentEdge, m_polyhedron), m_polyhedron);
    }
    while (currentEdge != startEdge);

  }

  Segment_2 clip_to_bounds(Segment_2 segment, Ray_2 leftBoundary, Ray_2 rightBoundary)
  {
    typedef typename cpp11::result_of<Intersect_2(Segment_2, Ray_2)>::type SegmentRayIntersectResult;

    SegmentRayIntersectResult leftIntersection = m_traits.intersect_2_object()(segment, leftBoundary);
    Point_2 leftPoint;
    bool containsLeft = true;
    
    if (leftIntersection)
    {
      Point_2* result = boost::get<Point_2>(&*leftIntersection);
      
      if (result)
      {
        leftPoint = *result;
        containsLeft = false;
      }
    }
    
    if (containsLeft)
    {
      leftPoint = segment[0];
    }
    
    SegmentRayIntersectResult rightIntersection = m_traits.intersect_2_object()(segment, rightBoundary);
    Point_2 rightPoint;
    bool containsRight = true;
    
    if (rightIntersection)
    {
      Point_2* result = boost::get<Point_2>(&*rightIntersection);
      
      if (result)
      {
        rightPoint = *result;
        containsRight = false;
      }
    }
    
    if (containsRight)
    {
      rightPoint = segment[1];
    }
    
    return Segment_2(leftPoint, rightPoint);
  }

  void process_node(Cone_tree_node* node)
  {
    bool leftSide = node->has_left_side();
    bool rightSide = node->has_right_side();
    
    bool propagateLeft = false;
    bool propagateRight = false;
    bool propagateMiddle = false;
  
    if (m_debugOutput)
    {
      std::cout << " Processing node " << node << " , level = " << node->level() << std::endl;
      std::cout << "\tFace = " << node->layout_face() << std::endl;
      std::cout << "\tSource Image = " << node->source_image() << std::endl;
      std::cout << "\tWindow Left = " << node->window_left() << std::endl;
      std::cout << "\tWindow Right = " << node->window_right() << std::endl;
      std::cout << "\t Has Left : " << (leftSide ? "yes" : "no") << " , Has Right : " << (rightSide ? "yes" : "no") << std::endl;
    }
    
    if (node->is_source_node() || (leftSide && rightSide))
    {
      if (m_debugOutput)
      {
        std::cout << "\tContains target vertex" << std::endl;
      }
      
      size_t entryEdgeIndex = m_halfedgeIndexMap[node->entry_edge()];
      
      NodeDistancePair currentOccupier = m_vertexOccupiers[entryEdgeIndex];
      FT currentNodeDistance = node->distance_from_target_to_root();

      bool isLeftOfCurrent = false;
      
      if (m_debugOutput)
      {
        std::cout << "\t Target vertex = " << m_vertexIndexMap[node->target_vertex()] << std::endl;
      }
      
      if (currentOccupier.first != NULL)
      {
        if (node->is_vertex_node())
        {
          isLeftOfCurrent = false;
        }
        else if (currentOccupier.first->is_vertex_node())
        {
          isLeftOfCurrent = true;
        }
        else
        {
          CGAL::Comparison_result comparison = m_traits.compare_relative_intersection_along_segment_2_object()(
            node->entry_segment(), 
            node->ray_to_target_vertex(), 
            currentOccupier.first->entry_segment(),
            currentOccupier.first->ray_to_target_vertex()
          );
          
          if (comparison == CGAL::SMALLER)
          {
            isLeftOfCurrent = true;
          }
        }
        
        if (m_debugOutput)
        {
          std::cout << "\t Current occupier = " << currentOccupier.first << std::endl;
          std::cout << "\t Current Occupier Distance = " << currentOccupier.second << std::endl;
          std::cout << "\t " << (isLeftOfCurrent ? "Left" : "Right") << " of current" << std::endl;
        }
      }
      
      if (m_debugOutput)
      {
        std::cout << "\t New Distance = " << currentNodeDistance << std::endl;
      }

      if (currentOccupier.first == NULL || currentOccupier.second > currentNodeDistance)
      {
        if (m_debugOutput)
        {
          std::cout << "\t Current node is now the occupier" << std::endl;
        }
        
        m_vertexOccupiers[entryEdgeIndex] = std::make_pair(node, currentNodeDistance);
        
        propagateLeft = true;
        propagateRight = true;
        
        if (node->node_type() != Cone_tree_node::INTERVAL)
        {
          propagateRight = false;
        }
        
        if (currentOccupier.first != NULL)
        {
          if (isLeftOfCurrent)
          {
            if (currentOccupier.first->get_left_child())
            {
              delete_node(currentOccupier.first->remove_left_child());
            }
            else if (currentOccupier.first->m_pendingLeftSubtree != NULL)
            {
              currentOccupier.first->m_pendingLeftSubtree->m_cancelled = true;
              currentOccupier.first->m_pendingLeftSubtree = NULL;
            }
          }
          else
          {
            if (currentOccupier.first->get_right_child())
            {
              delete_node(currentOccupier.first->remove_right_child());
            }
            else if (currentOccupier.first->m_pendingRightSubtree != NULL)
            {
              currentOccupier.first->m_pendingRightSubtree->m_cancelled = true;
              currentOccupier.first->m_pendingRightSubtree = NULL;
            }
          }
        }
        
        size_t targetVertexIndex = m_vertexIndexMap[node->target_vertex()];
        
        // Check if this is now the absolute closest node, and replace the current closest as appropriate
        NodeDistancePair currentClosest = m_closestToVertices[targetVertexIndex];
        
        if (m_debugOutput && currentClosest.first != NULL)
        {
          std::cout << "\t Current Closest Distance = " << currentClosest.second << std::endl;
        }
        
        if (currentClosest.first == NULL || currentClosest.second > currentNodeDistance)
        {
          if (m_debugOutput)
          {
            std::cout << "\t Current node is now the closest" << std::endl;
          }

          // if this is a saddle vertex, then evict previous closest vertex
          if (m_vertexIsPseudoSource[targetVertexIndex])
          {
            if (currentClosest.first != NULL)
            {
              if (m_debugOutput)
              {
                std::cout << "\tEvicting old pseudo-source: " << currentClosest.first << std::endl;
              }
              
              if (currentClosest.first->m_pendingMiddleSubtree != NULL)
              {
                currentClosest.first->m_pendingMiddleSubtree->m_cancelled = true;
                currentClosest.first->m_pendingMiddleSubtree = NULL;
              }

              while (currentClosest.first->has_middle_children())
              {
                delete_node(currentClosest.first->pop_middle_child());
              }
              
              if (m_debugOutput)
              {
                std::cout << "\tFinished Evicting" << std::endl;
              }
            }

            propagateMiddle = true;
          }
          
          m_closestToVertices[targetVertexIndex] = NodeDistancePair(node, currentNodeDistance);
        }
      }
      else
      {
        if (isLeftOfCurrent)
        {
          propagateLeft = true;
        }
        else
        {
          propagateRight = true;
        }
      }
    }
    else
    {
      propagateLeft = leftSide;
      propagateRight = rightSide;
    }
    
    if (node->level() < num_faces(m_polyhedron))
    {
      if (propagateLeft)
      {
        push_left_child(node);
      }
      
      if (propagateRight)
      {
        push_right_child(node);
      }
      
      if (propagateMiddle)
      {
        push_middle_child(node);
      }
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNo expansion since level limit reached" << std::endl;
    }
  }

  void push_left_child(Cone_tree_node* parent)
  {
    Segment_2 leftWindow(clip_to_bounds(parent->left_child_base_segment(), parent->left_boundary(), parent->right_boundary()));
    FT distanceEstimate = std::min(parent->distance_to_root(leftWindow[0]), parent->distance_to_root(leftWindow[1]));
    
    if (m_debugOutput)
    {
      std::cout << "\tPushing Left Child, Segment = " << parent->left_child_base_segment() << " , clipped = " << leftWindow << " , Estimate = " << distanceEstimate << std::endl;
    }

    Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::LEFT_CHILD, leftWindow);
    parent->m_pendingLeftSubtree = event;
    
    m_expansionPriqueue.push(event);
  }

  void push_right_child(Cone_tree_node* parent)
  {
    Segment_2 rightWindow(clip_to_bounds(parent->right_child_base_segment(), parent->left_boundary(), parent->right_boundary()));
    FT distanceEstimate = std::min(parent->distance_to_root(rightWindow[0]), parent->distance_to_root(rightWindow[1]));
    
    if (m_debugOutput)
    {
      std::cout << "\tPushing Right Child, Segment = " << parent->right_child_base_segment() << " , clipped = " << rightWindow << " , Estimate = " << distanceEstimate << std::endl;
    }
    
    Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::RIGHT_CHILD, rightWindow);
    parent->m_pendingRightSubtree = event;

    m_expansionPriqueue.push(event);
  }

  void push_middle_child(Cone_tree_node* parent)
  {
    if (m_debugOutput)
    {
      std::cout << "\tPushing Middle Child, Estimate = " << parent->distance_from_target_to_root() << std::endl;
    }
    
    Cone_expansion_event* event = new Cone_expansion_event(parent, parent->distance_from_target_to_root(), Cone_expansion_event::PSEUDO_SOURCE);
    parent->m_pendingMiddleSubtree = event;
    
    m_expansionPriqueue.push(event);
  }
  
  void delete_node(Cone_tree_node* node)
  {
    if (node != NULL)
    {
      if (m_debugOutput)
      {
        std::cout << "Deleting node " << node << std::endl;
      }
      
      if (node->m_pendingLeftSubtree != NULL)
      {
        node->m_pendingLeftSubtree->m_cancelled = true;
        node->m_pendingLeftSubtree = NULL;
      }

      if (node->get_left_child() != NULL)
      {
        if (m_debugOutput)
        {
          std::cout << "\t"  << node << " Descending left." << std::endl;
        }
      
        delete_node(node->remove_left_child());
      }
      
      
      if (node->m_pendingRightSubtree != NULL)
      {
        node->m_pendingRightSubtree->m_cancelled = true;
        node->m_pendingRightSubtree = NULL;
      }

      if (node->get_right_child() != NULL)
      {
        if (m_debugOutput)
        {
          std::cout << "\t"  << node << " Descending right." << std::endl;
        }
        
        delete_node(node->remove_right_child());
      }
      
      if (node->m_pendingMiddleSubtree != NULL)
      {
        node->m_pendingMiddleSubtree->m_cancelled = true;
        node->m_pendingMiddleSubtree = NULL;
      }
      
      if (node->has_middle_children() && m_debugOutput)
      {
        std::cout << "\t"  << node << " Descending middle." << std::endl;
      }
      
      while (node->has_middle_children())
      {
        delete_node(node->pop_middle_child());
      }
      
      size_t entryEdgeIndex = m_halfedgeIndexMap[node->entry_edge()];
      
      if (m_vertexOccupiers[entryEdgeIndex].first == node)
      {
        m_vertexOccupiers[entryEdgeIndex].first = NULL;
        
        size_t targetVertexIndex = m_vertexIndexMap[node->target_vertex()];
        
        if (m_closestToVertices[targetVertexIndex].first == node)
        {
          m_closestToVertices[targetVertexIndex].first = NULL;
        }
      }
    }
  }

  void set_vertex_types()
  {
    vertex_iterator current, end;
    
    for (boost::tie(current, end) = boost::vertices(m_polyhedron); current != end; ++current)
    {
      size_t vertexIndex = m_vertexIndexMap[*current];
    
      if (is_saddle_vertex(*current) || is_boundary_vertex(*current))
      {
        m_vertexIsPseudoSource[vertexIndex] = true;
      }
      else
      {
        m_vertexIsPseudoSource[vertexIndex] = false;
      }
      
      m_closestToVertices[vertexIndex] = NodeDistancePair(NULL, FT(0.0));
    }
    
    halfedge_iterator currHe, endHe;
    
    m_vertexOccupiers.clear();
    
    for (boost::tie(currHe, endHe) = CGAL::halfedges(m_polyhedron); currHe != endHe; ++currHe)
    {
      m_vertexOccupiers[m_halfedgeIndexMap[*currHe]] = NodeDistancePair(NULL, FT(0.0));
    }
  }
  
  bool is_saddle_vertex(vertex_descriptor v)
  {
    return m_traits.is_saddle_vertex_object()(v, m_polyhedron);
  }
  
  bool is_boundary_vertex(vertex_descriptor v) // TODO: confirm that this actually works
  {
    halfedge_descriptor h = CGAL::halfedge(v, m_polyhedron);
    halfedge_descriptor first = h;
    
    do
    {
      if (CGAL::face(h, m_polyhedron) == GraphTraits::null_face() || CGAL::face(CGAL::opposite(h, m_polyhedron), m_polyhedron) == GraphTraits::null_face())
      {
        return true;
      }
      
      h = CGAL::opposite(CGAL::next(h, m_polyhedron), m_polyhedron);
    }
    while(h != first);
    
    return false;
  }
  
  void reset_containers()
  {
    m_closestToVertices.resize(boost::num_vertices(m_polyhedron));
    m_vertexOccupiers.resize(CGAL::num_halfedges(m_polyhedron));
    
    while (!m_expansionPriqueue.empty())
    {
      delete m_expansionPriqueue.top();
      m_expansionPriqueue.pop();
    }
    
    m_faceLocations.clear();
    m_rootNodes.clear();
    m_vertexIsPseudoSource.resize(boost::num_vertices(m_polyhedron));
  }
  
  template <class Visitor>
  void visit_shortest_path(Cone_tree_node* startNode, const Point_2& startLocation, Visitor& visitor)
  {
    typedef typename Traits::Intersect_2 Intersect_2;
    typedef typename cpp11::result_of<Intersect_2(Segment_2, Ray_2)>::type SegmentRayIntersectResult;
    
    Cone_tree_node* current = startNode;
    Point_2 currentLocation(startLocation);
    
    while (!current->is_root_node())
    {
      switch (current->node_type())
      {
        case Cone_tree_node::INTERVAL:
        case Cone_tree_node::EDGE_SOURCE:
        {
          Segment_2 entrySegment = current->entry_segment();
          Ray_2 rayToLocation(current->source_image(), currentLocation);
          SegmentRayIntersectResult intersection = m_traits.intersect_2_object()(entrySegment, rayToLocation);
          assert(intersection && "Line from source did not cross entry segment");
          Point_2* result = boost::get<Point_2>(&*intersection);
          assert(result != NULL && "Intersection with entry segment was not a single point");
          FT parametricLocation = m_traits.parameteric_distance_along_segment_2_object()(entrySegment[0], entrySegment[1], *result);
          visitor.edge(current->entry_edge(), parametricLocation);

          if (current->is_left_child())
          {
            Segment_2 baseSegment = current->parent()->left_child_base_segment();
            currentLocation = CGAL::internal::interpolate_points<Point_2, Vector_2, FT>(baseSegment[0], baseSegment[1], parametricLocation);
          }
          else if (current->is_right_child())
          {
            Segment_2 baseSegment = current->parent()->right_child_base_segment();
            currentLocation = CGAL::internal::interpolate_points<Point_2, Vector_2, FT>(baseSegment[0], baseSegment[1], parametricLocation);
          }
          
          current = current->parent();
        }
          break;
        case Cone_tree_node::VERTEX_SOURCE:
          visitor.vertex(CGAL::target(current->entry_edge(), m_polyhedron));
          currentLocation = current->parent()->target_vertex_location();
          current = current->parent();
          break;
        case Cone_tree_node::FACE_SOURCE:
          // This is guaranteed to be the final node in any sequence
          visitor.face(m_faceLocations[current->tree_id()].first, m_faceLocations[current->tree_id()].second);
          current = current->parent();
          break;
      }
    }
  }
  
  void add_to_face_list(Cone_tree_node* node)
  {
    if (!node->is_root_node())
    {
      size_t faceIndex = m_faceIndexMap[node->current_face()];
      m_faceOccupiers[faceIndex].push_back(node);
    }
    
    if (node->get_left_child() != NULL)
    {
      add_to_face_list(node->get_left_child());
    }
    
    if (node->get_right_child() != NULL)
    {
      add_to_face_list(node->get_right_child());
    }
    
    for (size_t i = 0; node->num_middle_children(); ++i)
    {
      add_to_face_list(node->get_middle_child(i));
    }
  }
  
public:
  
  Polyhedron_shortest_path(const Traits& traits, Polyhedron& p)
    : m_traits(traits)
    , m_polyhedron(p)
    , m_vertexIndexMap(CGAL::get(boost::vertex_external_index, p))
    , m_halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, p))
    , m_faceIndexMap(CGAL::get(CGAL::face_external_index, p))
    , m_vertexPointMap(CGAL::get(CGAL::vertex_point, p))
    , m_debugOutput(false)
  {
  }
  
  Polyhedron_shortest_path(const Traits& traits, Polyhedron& p, VertexIndexMap& vertexIndexMap, HalfedgeIndexMap& halfedgeIndexMap, FaceIndexMap& faceIndexMap, VertexPointMap& vertexPointMap)
    : m_traits(traits)
    , m_polyhedron(p)
    , m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
    , m_vertexPointMap(vertexPointMap)
    , m_debugOutput(false)
  {
  }
  
  void compute_shortest_paths(face_descriptor face, Barycentric_coordinate location)
  {
    typedef FaceLocationPair* FaceLocationPairIterator;

    FaceLocationPair faceLocation(std::make_pair(face, location));
    compute_shortest_paths<FaceLocationPairIterator>(&faceLocation, (&faceLocation) + 1);
  }
  
  template<class InputIterator>
  void compute_shortest_paths(InputIterator faceLocationsBegin, InputIterator faceLocationsEnd)
  {
    reset_containers();
    set_vertex_types();
    
    m_vertexOccupiers.resize(CGAL::num_halfedges(m_polyhedron));
    m_closestToVertices.resize(CGAL::num_vertices(m_polyhedron));

    if (m_debugOutput)
    {
      vertex_iterator current, end;
      
      size_t numVertices = 0;

      for (boost::tie(current,end) = boost::vertices(m_polyhedron); current != end; ++current)
      {
        std::cout << "Vertex#" << numVertices << ": p = " << m_vertexPointMap[*current] << " , Concave: " << (m_vertexIsPseudoSource[numVertices] ? "yes" : "no") << std::endl;
        ++numVertices;
      }
    }
    
    face_iterator facesCurrent;
    face_iterator facesEnd;
    
    if (m_debugOutput)
    {
      size_t numFaces = 0;
      
      for (boost::tie(facesCurrent, facesEnd) = CGAL::faces(m_polyhedron); facesCurrent != facesEnd; ++facesCurrent)
      {
        std::cout << "Face#" << numFaces << ": Vertices = (";
        ++numFaces;
        halfedge_iterator faceEdgesStart = CGAL::halfedge(*facesCurrent, m_polyhedron);
        halfedge_iterator faceEdgesCurrent = faceEdgesStart;
        
        do
        {
          std::cout << m_vertexIndexMap[CGAL::source(*faceEdgesCurrent, m_polyhedron)];
            
          faceEdgesCurrent = CGAL::next(*faceEdgesCurrent, m_polyhedron);
          
          if (faceEdgesCurrent != faceEdgesStart)
          {
            std::cout << ", ";
          }
          else
          {
            std::cout << ")";
          }
        }
        while (faceEdgesCurrent != faceEdgesStart);
        
        std::cout << std::endl;
      }
    
    }
    
    for (InputIterator it = faceLocationsBegin; it != faceLocationsEnd; ++it)
    {
      m_faceLocations.push_back(*it);
      
      if (m_debugOutput)
      {
        std::cout << "Root: " << m_faceIndexMap[it->first] << " , " << it->second << std::endl;
      }
      
      expand_root(it->first, it->second);
    }
    
    if (m_debugOutput)
    {
      std::cout << "PriQ start size = " << m_expansionPriqueue.size() << std::endl;

      std::cout << "Num face locations: " << m_faceLocations.size() << std::endl;
      std::cout << "Num root nodes: " << m_rootNodes.size() << " (Hint: these should be the same size)" << std::endl;
   
    }
    
    while (m_expansionPriqueue.size() > 0)
    {
      Cone_expansion_event* event = m_expansionPriqueue.top();
      m_expansionPriqueue.pop();
      
      if (!event->m_cancelled)
      {
        typename Cone_expansion_event::Expansion_type type = event->m_type;
        Cone_tree_node* parent = event->m_parent;

        switch (type)
        {
          case Cone_expansion_event::PSEUDO_SOURCE:
            if (m_debugOutput)
            {
              std::cout << "PseudoSource Expansion: Parent = " << parent << " , Vertex = " << m_vertexIndexMap[event->m_parent->target_vertex()] << " , Distance = " << event->m_distanceEstimate << " , Level = " << event->m_parent->level() + 1 << std::endl;
            }
            
            expand_psuedo_source(parent);
            break;
          case Cone_expansion_event::LEFT_CHILD:
            if (m_debugOutput)
            {
              std::cout << "Left Expansion: Parent = " << parent << " Edge = (" << m_vertexIndexMap[CGAL::source(event->m_parent->left_child_edge(), m_polyhedron)] << "," << m_vertexIndexMap[CGAL::target(event->m_parent->left_child_edge(), m_polyhedron)] << ") , Distance = " << event->m_distanceEstimate << " , Level = " << event->m_parent->level() + 1 << std::endl;
            }
            
            expand_left_child(parent, event->m_windowSegment);
            break;
          case Cone_expansion_event::RIGHT_CHILD:
            if (m_debugOutput)
            {
              std::cout << "Right Expansion: Parent = " << parent << " , Edge = (" << m_vertexIndexMap[CGAL::source(event->m_parent->right_child_edge(), m_polyhedron)] << "," << m_vertexIndexMap[CGAL::target(event->m_parent->right_child_edge(), m_polyhedron)] << ") , Distance = " << event->m_distanceEstimate << " , Level = " << event->m_parent->level() + 1 << std::endl;
            }
            
            expand_right_child(parent, event->m_windowSegment);
            break;
        }
      }
      else if (m_debugOutput)
      {
        std::cout << "Found cancelled event for node: " << event->m_parent << std::endl;
      }
      
      delete event;
    }
    
    if (m_debugOutput)
    {   
      std::cout << "Closest distances: " << std::endl;
      
      for (size_t i = 0; i < m_closestToVertices.size(); ++i)
      {
        std::cout << "\tVertex = " << i << std::endl;
        std::cout << "\tDistance = " << m_closestToVertices[i].second << std::endl;
      }
      
      std::cout << std::endl << "Done!" << std::endl;
    }
    
    m_faceOccupiers.clear();
    m_faceOccupiers.resize(CGAL::num_faces(m_polyhedron));
    
    for (size_t i = 0; i < m_rootNodes.size(); ++i)
    {
      add_to_face_list(m_rootNodes[i]);
    }
  }

  FT shortest_distance_to_vertex(vertex_descriptor v)
  {
    return m_closestToVertices[m_vertexIndexMap[v]].second;
  }
  
  template <class Visitor>
  void shortest_edge_sequence(vertex_descriptor v, Visitor& visitor)
  {
    Cone_tree_node* current = m_closestToVertices[m_vertexIndexMap[v]].first;
    
    visit_shortest_path(current, current->target_vertex_location(), visitor);
  }
  
  template <class Visitor>
  void shortest_path_points(vertex_descriptor v, Visitor& visitor)
  {
    shortest_edge_sequence(v, Point_path_visitor_wrapper<Visitor>(visitor, m_polyhedron, m_vertexPointMap));
  }

};

} // namespace CGAL
