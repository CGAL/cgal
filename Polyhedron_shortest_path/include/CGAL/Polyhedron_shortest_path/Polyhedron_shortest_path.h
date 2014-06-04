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
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

template<class Traits>
class Polyhedron_shortest_path
{
public:
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

  typedef typename std::priority_queue<Cone_expansion_event, std::vector<Cone_expansion_event>, internal::Cone_expansion_event_min_priority_queue_comparator<Traits> > Expansion_priqueue;
  
  typedef typename std::pair<face_descriptor, Barycentric_coordinate> FaceLocationPair;
  
  typedef typename std::pair<Cone_tree_node*, FT> NodeDistancePair;

private:
  Traits m_traits;

private:
  typedef typename std::map<vertex_descriptor, bool> PsuedoSourceMap;

  PsuedoSourceMap m_vertexIsPsuedoSource;
  
  std::map<halfedge_descriptor, NodeDistancePair> m_vertexOccupiers;
  std::map<vertex_descriptor, NodeDistancePair> m_closestToVertices;
  
  std::vector<Cone_tree_node*> m_rootNodes;
  
  Expansion_priqueue m_expansionPriqueue;
  
  std::vector<FaceLocationPair> m_faceLocations;
  
  Polyhedron* m_polyhedron;
  
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
      points[currentPoint] = boost::source(current, *m_polyhedron)->point();
      current = CGAL::next(current, *m_polyhedron);
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
    
    NodeDistancePair v1Distance = m_closestToVertices[CGAL::source(cone->entry_edge(), *m_polyhedron)];
    NodeDistancePair v2Distance = m_closestToVertices[cone->target_vertex()];
    NodeDistancePair v3Distance = m_closestToVertices[CGAL::target(cone->entry_edge(), *m_polyhedron)];
    
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
    if (cone->m_hasLeftSubtree && window_distance_filter(cone, windowSegment, false))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->left_child_edge());
      Triangle_2 layoutFace = m_traits.flatten_triangle_3_along_segment_2_object()(adjacentFace, 0, cone->left_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, cone->left_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), windowSegment[0], windowSegment[1], false);
      cone->set_left_child(child);
      process_node(child);
    }
    else
    {
      if (!cone->m_hasLeftSubtree)
      {
        std::cout << "\tNode was evicted." << std::endl;
      }
      else
      {
        std::cout << "\tNode was filtered." << std::endl;
      }
    }
  }
  
  void expand_right_child(Cone_tree_node* cone, Segment_2 windowSegment)
  {
    if (cone->m_hasRightSubtree && window_distance_filter(cone, windowSegment, true))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->right_child_edge());
      Triangle_2 layoutFace = m_traits.flatten_triangle_3_along_segment_2_object()(adjacentFace, 0, cone->right_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, cone->right_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), windowSegment[0], windowSegment[1], false);
      cone->set_right_child(child);
      process_node(child);
    }
    else
    {
      if (!cone->m_hasRightSubtree)
      {
        std::cout << "\tNode was evicted." << std::endl;
      }
      else
      {
        std::cout << "\tNode was filtered." << std::endl;
      }
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
          halfedge_descriptor he = CGAL::halfedge(face, *m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            he = CGAL::next(he, *m_polyhedron);
          }
          expand_edge_root(he, location[associatedEdge], location[(associatedEdge + 1) % 3]);
        }
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX:
        {
          halfedge_descriptor he = CGAL::halfedge(face, *m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            he = CGAL::next(he, *m_polyhedron);
          }
          expand_vertex_root(CGAL::source(he, *m_polyhedron));
        }
        break;
      default:
        assert(false && "Invalid face location");
        // Perhaps hit an assertion that the type must not be external or invalid?
    }
  }
  
  void expand_face_root(face_descriptor faceId, Barycentric_coordinate faceLocation)
  {
    halfedge_descriptor start = CGAL::halfedge(faceId, *m_polyhedron);
    halfedge_descriptor current = start;
    
    Cone_tree_node* faceRoot = new Cone_tree_node(m_polyhedron, m_rootNodes.size());
    m_rootNodes.push_back(faceRoot);
    
    if (m_debugOutput)
    {
      std::cout << "\tFace Root Expansion: face = " << m_facesMap[faceId] << " , Location = " << faceLocation << std::endl;
    }
    
    for (size_t currentVertex = 0; currentVertex < 3; ++currentVertex)
    {
      Triangle_3 face3d(triangle_from_halfedge(current));
      Triangle_2 layoutFace(m_traits.project_triangle_3_to_triangle_2_object()(face3d));
      Barycentric_coordinate rotatedFaceLocation(faceLocation[currentVertex], faceLocation[(currentVertex + 1) % 3], faceLocation[(currentVertex + 2) % 3]);
      Point_2 sourcePoint(m_traits.construct_triangle_location_2_object()(layoutFace, rotatedFaceLocation));
      
      Cone_tree_node* child = new Cone_tree_node(m_polyhedron, current, layoutFace, sourcePoint, FT(0.0), layoutFace[0], layoutFace[2], true);
      faceRoot->push_middle_child(child);
      
      if (m_debugOutput)
      {
        std::cout << "\tExpanding face root #" << currentVertex << " : " << std::endl;;
        std::cout << "\t\tFace = " << layoutFace << std::endl;
        std::cout << "\t\tLocation = " << sourcePoint << std::endl;
      }
      process_node(child);

      current = CGAL::next(current, *m_polyhedron);
    }
  }

  void expand_edge_root(halfedge_descriptor baseEdge, FT t0, FT t1)
  {
    std::cout << "\tEdge Root Expansion: faceA = " << m_facesMap[CGAL::face(baseEdge, *m_polyhedron)] << " , faceB = " << m_facesMap[CGAL::face(CGAL::opposite(baseEdge, *m_polyhedron), *m_polyhedron)] << " , t0 = " << t0 << " , t1 = " << t1 << std::endl;
  
    halfedge_descriptor baseEdges[2];
    baseEdges[0] = baseEdge;
    baseEdges[1] = CGAL::opposite(baseEdge, *m_polyhedron);
    
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
      std::cout << "\tExpanding edge root #" << side << " : " << std::endl;;
      std::cout << "\t\tFace = " << layoutFaces[side] << std::endl;
      std::cout << "\t\tLocation = " << sourcePoints[side] << std::endl;
    
      Cone_tree_node* mainChild = new Cone_tree_node(m_polyhedron, baseEdges[side], layoutFaces[side], sourcePoints[side], FT(0.0), layoutFaces[side][0], layoutFaces[side][2], true);
      edgeRoot->push_middle_child(mainChild);
      process_node(mainChild);

      Cone_tree_node* oppositeChild = new Cone_tree_node(m_polyhedron, baseEdges[side], Triangle_2(layoutFaces[side][2], layoutFaces[side][1], layoutFaces[side][2]), sourcePoints[side], FT(0.0), layoutFaces[side][1], layoutFaces[side][2], true);
      edgeRoot->push_middle_child(oppositeChild);
      process_node(oppositeChild);
    }
  }

  void expand_vertex_root(vertex_descriptor vertex)
  {
    std::cout << "\tVertex Root Expansion: Vertex = " << m_vertexMap[vertex] << std::endl;

    Cone_tree_node* vertexRoot = new Cone_tree_node(m_polyhedron, m_rootNodes.size(), CGAL::halfedge(vertex, *m_polyhedron));
    m_rootNodes.push_back(vertexRoot);
    
    vertexRoot->m_hasMiddleSubtree = true;
    m_closestToVertices[vertex] = NodeDistancePair(vertexRoot, FT(0.0));
    
    expand_psuedo_source(vertexRoot);
  }

  void expand_psuedo_source(Cone_tree_node* parent)
  {
    if (parent->m_hasMiddleSubtree)
    {
      vertex_descriptor expansionVertex = parent->target_vertex();
    
      halfedge_descriptor startEdge = CGAL::halfedge(expansionVertex, *m_polyhedron);
      halfedge_descriptor currentEdge = CGAL::halfedge(expansionVertex, *m_polyhedron);
          
      do
      {
        Triangle_3 face3d(triangle_from_halfedge(currentEdge));
        Triangle_2 layoutFace(m_traits.project_triangle_3_to_triangle_2_object()(face3d));
        
        if (m_debugOutput)
        {
          std::cout << "\tExpanding PsuedoSource: id = " << m_facesMap[CGAL::face(currentEdge, *m_polyhedron)] << " , face = " << layoutFace << std::endl;
        }
        
        Cone_tree_node* child = new Cone_tree_node(m_polyhedron, currentEdge, layoutFace, layoutFace[1], FT(0.0), layoutFace[0], layoutFace[2], true);
        parent->push_middle_child(child);
        process_node(child);
        
        currentEdge = CGAL::opposite(CGAL::next(currentEdge, *m_polyhedron), *m_polyhedron);
      }
      while (currentEdge != startEdge);
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNode was evicted." << std::endl;
    }
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
      
      NodeDistancePair currentOccupier = m_vertexOccupiers[node->entry_edge()];
      FT currentNodeDistance = node->distance_from_target_to_root();

      bool isLeftOfCurrent = false;
      
      if (m_debugOutput)
      {
        std::cout << "\t Target vertex = " << m_vertexMap[node->target_vertex()] << std::endl;
      }
      
      if (currentOccupier.first != NULL)
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
        
        m_vertexOccupiers[node->entry_edge()] = std::make_pair(node, currentNodeDistance);
        
        node->m_hasLeftSubtree = true;
        node->m_hasRightSubtree = true;
        
        if (node->is_source_node())
        {
          node->m_hasRightSubtree = false;
        }
        
        if (currentOccupier.first != NULL)
        {
          if (isLeftOfCurrent)
          {
            currentOccupier.first->m_hasLeftSubtree = false;
            
            if (currentOccupier.first->get_left_child())
            {
              delete_node(currentOccupier.first->remove_left_child());
            }
          }
          else
          {
            currentOccupier.first->m_hasRightSubtree = false;
            
            if (currentOccupier.first->get_right_child())
            {
              delete_node(currentOccupier.first->remove_right_child());
            }
          }
        }
        
        // Check if this is now the absolute closest node, and replace the current closest as appropriate
        NodeDistancePair currentClosest = m_closestToVertices[node->target_vertex()];
        
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
          if (m_vertexIsPsuedoSource[node->target_vertex()])
          {
            if (currentClosest.first != NULL)
            {
              currentClosest.first->m_hasMiddleSubtree = false;
              
              while (currentClosest.first->has_middle_children())
              {
                delete_node(currentClosest.first->pop_middle_child());
              }
            }

            node->m_hasMiddleSubtree = true;
          }
          
          m_closestToVertices[node->target_vertex()] = NodeDistancePair(node, currentNodeDistance);
        }
      }
      else
      {
        if (isLeftOfCurrent)
        {
          node->m_hasLeftSubtree = true;
        }
        else
        {
          node->m_hasRightSubtree = true;
        }
      }
    }
    else
    {
      node->m_hasLeftSubtree = leftSide;
      node->m_hasRightSubtree = rightSide;
    }
    
    if (node->level() < num_faces(*m_polyhedron))
    {
      if (node->m_hasLeftSubtree)
      {
        push_left_child(node);
      }
      
      if (node->m_hasRightSubtree)
      {
        push_right_child(node);
      }
      
      if (node->m_hasMiddleSubtree)
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

    m_expansionPriqueue.push(Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::LEFT_CHILD, leftWindow));
  }

  void push_right_child(Cone_tree_node* parent)
  {
    Segment_2 rightWindow(clip_to_bounds(parent->right_child_base_segment(), parent->left_boundary(), parent->right_boundary()));
    FT distanceEstimate = std::min(parent->distance_to_root(rightWindow[0]), parent->distance_to_root(rightWindow[1]));
    
    if (m_debugOutput)
    {
      std::cout << "\tPushing Right Child, Segment = " << parent->right_child_base_segment() << " , clipped = " << rightWindow << " , Estimate = " << distanceEstimate << std::endl;
    }

    m_expansionPriqueue.push(Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::RIGHT_CHILD, rightWindow));
  }

  void push_middle_child(Cone_tree_node* parent)
  {
    if (m_debugOutput)
    {
      std::cout << "\tPushing Middle Child, Estimate = " << parent->distance_from_target_to_root() << std::endl;
    }
  
    m_expansionPriqueue.push(Cone_expansion_event(parent, parent->distance_from_target_to_root(), Cone_expansion_event::PSEUDO_SOURCE));
  }
  
  void delete_node(Cone_tree_node* node)
  {
    if (node != NULL)
    {
      if (m_debugOutput)
      {
        std::cout << "Deleting node " << node << std::endl;
      }
      
      delete_node(node->remove_left_child());
      delete_node(node->remove_right_child());
      
      while (node->has_middle_children())
      {
        delete_node(node->pop_middle_child());
      }
      
      if (m_vertexOccupiers[node->entry_edge()].first == node)
      {
        m_vertexOccupiers[node->entry_edge()].first = NULL;
        
        if (m_closestToVertices[node->target_vertex()].first == node)
        {
          m_closestToVertices[node->target_vertex()].first = NULL;
        }
      }
    }
  }

  void set_vertex_types()
  {
    vertex_iterator current, end;
    
    m_vertexIsPsuedoSource.clear();
    
    for (boost::tie(current, end) = boost::vertices(*m_polyhedron); current != end; ++current)
    {
      if (is_saddle_vertex(*current) || is_boundary_vertex(*current))
      {
        m_vertexIsPsuedoSource[*current] = true;
      }
      else
      {
        m_vertexIsPsuedoSource[*current] = false;
      }
    }
  }
  
  bool is_saddle_vertex(vertex_descriptor v)
  {
    return m_traits.is_saddle_vertex_object()(v);
  }
  
  bool is_boundary_vertex(vertex_descriptor v) // TODO: confirm that this actually works
  {
    halfedge_descriptor h = CGAL::halfedge(v, *m_polyhedron);
    halfedge_descriptor first = h;
    
    do
    {
      if (h->is_border_edge())
      {
        return true;
      }
      
      h = CGAL::opposite(CGAL::next(h, *m_polyhedron), *m_polyhedron);
    }
    while(h != first);
    
    return false;
  }
  
  void reset_containers()
  {
    m_closestToVertices.clear();
    m_vertexOccupiers.clear();
    
    while (!m_expansionPriqueue.empty())
    {
      m_expansionPriqueue.pop();
    }
    
    m_faceLocations.clear();
    m_rootNodes.clear();
    m_vertexMap.clear();
    m_vertexIsPsuedoSource.clear();
  }
  
  std::map<vertex_descriptor, size_t> m_vertexMap;
  std::map<face_descriptor, size_t> m_facesMap;
  
public:

  Polyhedron_shortest_path(const Traits& traits)
    : m_traits(traits)
    , m_polyhedron(NULL)
    , m_debugOutput(false)
  {
  }
  
  void compute_shortest_paths(Polyhedron& p, face_descriptor face, Barycentric_coordinate location)
  {
    typedef FaceLocationPair* FaceLocationPairIterator;

    FaceLocationPair faceLocation(std::make_pair(face, location));
    compute_shortest_paths<FaceLocationPairIterator>(p, &faceLocation, (&faceLocation) + 1);
  }
  
  template<class InputIterator>
  void compute_shortest_paths(Polyhedron& p, InputIterator faceLocationsBegin, InputIterator faceLocationsEnd)
  {
    m_polyhedron = &p;
    
    reset_containers();
    set_vertex_types();

    size_t vertexCount = 0;
    
    for (typename PsuedoSourceMap::iterator it = m_vertexIsPsuedoSource.begin(); it != m_vertexIsPsuedoSource.end(); ++it)
    {
      m_vertexMap[it->first] = vertexCount;
      ++vertexCount;
    
      if (m_debugOutput)
      {
        std::cout << "Vertex#" << vertexCount << ": p = " << it->first->point() << " , Concave: " << (it->second ? "yes" : "no") << std::endl;
      }
    }
    
    face_iterator facesCurrent;
    face_iterator facesEnd;
    
    if (m_debugOutput)
    {
      size_t faceCount = 0;
      
      for (boost::tie(facesCurrent, facesEnd) = CGAL::faces(*m_polyhedron); facesCurrent != facesEnd; ++facesCurrent)
      {
        m_facesMap[*facesCurrent] = faceCount;
        ++faceCount;

        std::cout << "Face#" << faceCount << ": Vertices = (";

        halfedge_iterator faceEdgesStart = CGAL::halfedge(*facesCurrent, *m_polyhedron);
        halfedge_iterator faceEdgesCurrent = faceEdgesStart;
        
        do
        {
          std::cout << m_vertexMap[CGAL::source(*faceEdgesCurrent, *m_polyhedron)];
            
          faceEdgesCurrent = CGAL::next(*faceEdgesCurrent, *m_polyhedron);
          
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
        std::cout << "Root: " << m_facesMap[it->first] << " , " << it->second << std::endl;
      }
      
      expand_root(it->first, it->second);
    }
    
    if (m_debugOutput)
    {
      std::cout << "PriQ start size = " << m_expansionPriqueue.size() << std::endl;

      std::cout << "Num face locations: " << m_faceLocations.size() << std::endl;
      std::cout << "Num root nodes: " << m_rootNodes.size() << " (Hint: these should be the same size)" << std::endl;
    
      for (size_t i = 0; i < m_rootNodes.size(); ++i)
      {
        std::cout << "Root Node #" << i << ": " << std::endl;
      }
    }
    
    while (m_expansionPriqueue.size() > 0)
    {
      Cone_expansion_event event = m_expansionPriqueue.top();
      m_expansionPriqueue.pop();
      typename Cone_expansion_event::Expansion_type type = event.m_type;
      Cone_tree_node* parent = event.m_parent;

      switch (type)
      {
        case Cone_expansion_event::PSEUDO_SOURCE:
          if (m_debugOutput)
          {
            std::cout << "PseudoSource Expansion: Parent = " << parent << " , Vertex = " << m_vertexMap[event.m_parent->target_vertex()] << " , Distance = " << event.m_distanceEstimate << " , Level = " << event.m_parent->level() + 1 << std::endl;
          }
          
          expand_psuedo_source(parent);
          break;
        case Cone_expansion_event::LEFT_CHILD:
          if (m_debugOutput)
          {
            std::cout << "Left Expansion: Parent = " << parent << " Edge = (" << m_vertexMap[CGAL::source(event.m_parent->left_child_edge(), *m_polyhedron)] << "," << m_vertexMap[CGAL::target(event.m_parent->left_child_edge(), *m_polyhedron)] << ") , Distance = " << event.m_distanceEstimate << " , Level = " << event.m_parent->level() + 1 << std::endl;
          }
          
          expand_left_child(parent, event.m_windowSegment);
          break;
        case Cone_expansion_event::RIGHT_CHILD:
          if (m_debugOutput)
          {
            std::cout << "Right Expansion: Parent = " << parent << " , Edge = (" << m_vertexMap[CGAL::source(event.m_parent->right_child_edge(), *m_polyhedron)] << "," << m_vertexMap[CGAL::target(event.m_parent->right_child_edge(), *m_polyhedron)] << ") , Distance = " << event.m_distanceEstimate << " , Level = " << event.m_parent->level() + 1 << std::endl;
          }
          
          expand_right_child(parent, event.m_windowSegment);
          break;
      }
    }
    
    if (m_debugOutput)
    {   
      std::cout << "Closest distances: " << std::endl;
      
      for (typename std::map<vertex_descriptor, NodeDistancePair>::iterator it = m_closestToVertices.begin(); it != m_closestToVertices.end(); ++it)
      {
        std::cout << "\tVertex = " << m_vertexMap[it->first] << std::endl;
        std::cout << "\tDistance = " << it->second.second << std::endl;
      }
      
      std::cout << std::endl << "Done!" << std::endl;
    }
  }
  
  FT shortest_distance_to_vertex(vertex_descriptor v)
  {
    if (m_polyhedron != NULL)
    {
      return m_closestToVertices[v].second;
    }
    
    return FT(-1.0);
  }

};

} // namespace CGAL
