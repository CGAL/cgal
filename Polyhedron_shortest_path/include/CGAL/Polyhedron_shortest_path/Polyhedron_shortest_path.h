// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_POLYHEDRON_SHORTEST_PATH_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_POLYHEDRON_SHORTEST_PATH_H

#include <iterator>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <iterator>

#include <CGAL/Polyhedron_shortest_path/internal/Cone_tree.h>
#include <CGAL/Polyhedron_shortest_path/internal/misc_functions.h>
#include <CGAL/Polyhedron_shortest_path/internal/Barycentric.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

namespace CGAL {

/*!
\ingroup PkgPolyhedronShortestPath

\brief Computes shortest surface paths from one or more source points on a polyhedral surface

\details Uses an optimized variation of Chen and Han's O(n^2) algorithm by Xin and Wang. 
Refer to those respective papers for the details of the implementation.
 
\tparam Traits The geometric traits for this algorithm, a model of PolyhedronShortestPathTraits concept.

\tparam VIM A model of the boost ReadablePropertyMap concept, provides a vertex index property map.

\tparam HIM A model of the boost ReadablePropertyMap concept, provides a halfedges index property map.

\tparam FIM A model of the boost ReadablePropertyMap concept, provides a face index property map.

\tparam VPM A model of the boost ReadablePropertyMap concept, provides a vertex point property map.

 */
 
template<class Traits, 
  class VIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_external_index_t>::type,
  class HIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::halfedge_external_index_t>::type,
  class FIM = typename boost::property_map<typename Traits::Polyhedron, CGAL::face_external_index_t>::type,
  class VPM = typename boost::property_map<typename Traits::Polyhedron, CGAL::vertex_point_t>::type>
class Polyhedron_shortest_path
{
public:
/// \name Types
/// @{

  /// The vertex index property map class
  typedef VIM VertexIndexMap;
  
  /// The halfedge index property map class
  typedef HIM HalfedgeIndexMap;
  
  /// The face index property map class
  typedef FIM FaceIndexMap;
  
  /// The vertex point property map class
  typedef VPM VertexPointMap;

  /// The polyhedron type which this algorithm acts on.
  typedef typename Traits::Polyhedron Polyhedron;

  /// The BGL graph traits for this polyhedron
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;

  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  
  /// The numeric type used by this algorithm.
  typedef typename Traits::FT FT;
  
  /// The 3-dimensional point type of the polyhedron.
  typedef typename Traits::Point_3 Point_3;
  
  /// An ordered triple which specifies the location within a triangle as
  /// a convex combination of its three vertices.
  typedef typename Traits::Barycentric_coordinate Barycentric_coordinate;

  /// \brief An ordered pair specifying a location on the surface of the polyhedron.
  /// \detail Assuming you are given the pair (`face`, `location`), the weights of 
  /// `location` are applied to the vertices of `face` in the following way
  /// the following way:
  /// 0 - CGAL::source(CGAL::halfedge(`face`))
  /// 1 - CGAL::target(CGAL::halfedge(`face`))
  /// 2 - CGAL::target(CGAL::next(CGAL::halfedge(`face`)))
  typedef typename std::pair<face_descriptor, Barycentric_coordinate> Face_location;
  
  /// AABB primitive type for face graph polyhedron type
  typedef AABB_face_graph_triangle_primitive<Polyhedron, VertexPointMap> AABB_polyhedron_primitive;
  
  /// Traits class for point location AABB tree
  typedef CGAL::AABB_traits<typename Traits::Kernel, AABB_polyhedron_primitive> AABB_polyhedron_traits;
  
  /// An AABB tree to perform point location queries to get surface locations
  typedef AABB_tree<AABB_polyhedron_traits> AABB_polyhedron_tree;
  
/// @}
  
private:
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::halfedge_iterator halfedge_iterator;
  typedef typename GraphTraits::face_iterator face_iterator;

  typedef typename Traits::Triangle_3 Triangle_3;
  typedef typename Traits::Triangle_2 Triangle_2;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Ray_3 Ray_3;
  typedef typename Traits::Ray_2 Ray_2;
  typedef typename Traits::Line_2 Line_2;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Vector_2 Vector_2;

  typedef typename internal::Cone_tree_node<Traits> Cone_tree_node;
  typedef typename internal::Cone_expansion_event<Traits> Cone_expansion_event;

  typedef typename std::priority_queue<Cone_expansion_event, std::vector<Cone_expansion_event*>, internal::Cone_expansion_event_min_priority_queue_comparator<Traits> > Expansion_priqueue;
  typedef typename std::pair<Cone_tree_node*, FT> Node_distance_pair;
  
private:

  template <class OutputIterator>
  struct Point_path_visitor_wrapper
  {
    Polyhedron_shortest_path& m_owner;
    OutputIterator m_output;
    
    Point_path_visitor_wrapper(Polyhedron_shortest_path& owner, OutputIterator output)
      : m_owner(owner)
      , m_output(output)
    {
    }
    
    void edge(halfedge_descriptor edge, FT t)
    {
      *m_output = m_owner.get_edge_location(edge, t);
      ++m_output;
    }
    
    void vertex(vertex_descriptor vertex)
    {
      *m_output = m_owner.get_vertex_location(vertex);
      ++m_output;
    }
    
    void face(face_descriptor face, Barycentric_coordinate location)
    {
      *m_output = m_owner.get_face_location(face, location);
      ++m_output;
    }
  };

private:
  Traits m_traits;
  Polyhedron& m_polyhedron;
  VertexIndexMap m_vertexIndexMap;
  HalfedgeIndexMap m_halfedgeIndexMap;
  FaceIndexMap m_faceIndexMap;
  VertexPointMap m_vertexPointMap;

  std::vector<bool> m_vertexIsPseudoSource;
  
  std::vector<Node_distance_pair> m_vertexOccupiers;
  std::vector<Node_distance_pair> m_closestToVertices;
  
  std::vector<Cone_tree_node*> m_rootNodes;
  std::vector<Face_location> m_faceLocations;
  
  std::vector<std::vector<Cone_tree_node*> > m_faceOccupiers;
  
  Expansion_priqueue m_expansionPriqueue;
  
#if !defined(NDEBUG)
  
  size_t m_currentNodeCount;
  size_t m_peakNodeCount;
  size_t m_queueAtPeakNodes;
  size_t m_peakQueueSize;
  size_t m_nodesAtPeakQueue;
  
#endif

#if !defined(NDEBUG)
public:

  size_t peak_node_count()
  {
    return m_peakNodeCount;
  }
  
  size_t current_node_count()
  {
    return m_currentNodeCount;
  }
  
  size_t peak_queue_size()
  {
    return m_peakQueueSize;
  }
  
  size_t current_memory_usage()
  {
    size_t baseUsage = m_rootNodes.size() * sizeof(Cone_tree_node*) + m_closestToVertices.size() * sizeof(Node_distance_pair);

    size_t finalUsage = baseUsage + sizeof(Cone_tree_node) * m_currentNodeCount;
    
    for (size_t i = 0; i < m_faceOccupiers.size(); ++i)
    {
      finalUsage += (m_faceOccupiers[i].size() * sizeof(Cone_tree_node*)) + sizeof(std::vector<Cone_tree_node*>);
    }
    
    return finalUsage;
  }
  
  size_t peak_memory_usage()
  {
    size_t baseUsage = m_rootNodes.size() * sizeof(Cone_tree_node*) + m_vertexOccupiers.size() * sizeof(Node_distance_pair) + m_closestToVertices.size() * sizeof(Node_distance_pair);

    size_t peakNodeUsage = baseUsage + (sizeof(Cone_tree_node) * m_peakNodeCount) + ((sizeof(Cone_expansion_event) + sizeof(Cone_expansion_event*)) * m_queueAtPeakNodes);
    
    size_t peakQueueUsage = baseUsage + (sizeof(Cone_expansion_event) + (sizeof(Cone_expansion_event*)) * m_peakQueueSize) + (sizeof(Cone_tree_node) * m_nodesAtPeakQueue);
    
    return std::max(peakNodeUsage, peakQueueUsage);
  }

#endif


public:

  /// This is just a placeholder for a proper debug output verbosity switch method
  bool m_debugOutput;
  
private:

  void node_created()
  {
#if !defined(NDEBUG)
    ++m_currentNodeCount;
    if (m_currentNodeCount > m_peakNodeCount)
    {
      m_peakNodeCount = m_currentNodeCount;
      m_queueAtPeakNodes = m_expansionPriqueue.size();
    }
#endif
  }
  
  void queue_pushed()
  {
#if !defined(NDEBUG)
    if (m_expansionPriqueue.size() > m_peakQueueSize)
    {
      m_peakQueueSize = m_expansionPriqueue.size();
      m_nodesAtPeakQueue = m_currentNodeCount;
    }
#endif
  }
  
  void node_deleted()
  {
#if !defined(NDEBUG)
    --m_currentNodeCount;
#endif
  }

  Point_2 construct_barycenter_in_triangle_2(const Triangle_2& t, const Barycentric_coordinate& b) const 
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
    typename Traits::Construct_barycenter_2 cb2(m_traits.construct_barycenter_2_object());
    
    return cb2(cv2(t, 0), cbcw(b, 0), cv2(t, 1), cbcw(b, 1), cv2(t, 2), cbcw(b, 2));
  }
  
  Point_3 construct_barycenter_in_triangle_3(const Triangle_3& t, const Barycentric_coordinate& b) const 
  {
    typename Traits::Construct_vertex_3 cv3(m_traits.construct_vertex_3_object());
    typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
    typename Traits::Construct_barycenter_3 cb3(m_traits.construct_barycenter_3_object());
    
    return cb3(cv3(t, 0), cbcw(b, 0), cv3(t, 1), cbcw(b, 1), cv3(t, 2), cbcw(b, 2));
  }
  
  Triangle_3 triangle_from_halfedge(halfedge_descriptor edge) const
  {
    return CGAL::internal::triangle_from_halfedge<Triangle_3, Polyhedron, VertexPointMap>(edge, m_polyhedron, m_vertexPointMap);
  }
  
  Triangle_3 triangle_from_face(face_descriptor face) const
  {
    return triangle_from_halfedge(CGAL::halfedge(face, m_polyhedron));
  }

  bool window_distance_filter(Cone_tree_node* cone, Segment_2 windowSegment, bool reversed)
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());
  
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
    
    Node_distance_pair v1Distance = m_closestToVertices[v1Index];
    Node_distance_pair v2Distance = m_closestToVertices[v2Index];
    Node_distance_pair v3Distance = m_closestToVertices[v3Index];
    
    if (reversed)
    {
      std::swap(v1Distance, v3Distance);
      std::swap(v1Index, v3Index);
      A = cv2(windowSegment, 1);
      B = cv2(windowSegment, 0);
      v1 = cv2(parentEntrySegment, 1);
      v3 = cv2(parentEntrySegment, 0);
    }
    else
    {
      A = cv2(windowSegment, 0);
      B = cv2(windowSegment, 1);
      v1 = cv2(parentEntrySegment, 0);
      v3 = cv2(parentEntrySegment, 1);
    }
    
    d1 = v1Distance.second;
    d2 = v2Distance.second;
    d3 = v3Distance.second;
    
    bool hasD1 = v1Distance.first != NULL;
    bool hasD2 = v2Distance.first != NULL;
    bool hasD3 = v3Distance.first != NULL;
    
    if (hasD1 && (d + CGAL::internal::my_sqrt(csd2(I, B)) > d1 + CGAL::internal::my_sqrt(csd2(v1, B))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,B| > d1 + |v1,B|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v1,B = " << CGAL::internal::my_sqrt(csd2(v1, B)) << std::endl;
        std::cout << "I,B = " << CGAL::internal::my_sqrt(csd2(I, B)) << std::endl;
        std::cout << "I,A = " << CGAL::internal::my_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::internal::my_sqrt(csd2(I, B))) << " vs. " << (d1 + CGAL::internal::my_sqrt(csd2(v1, B))) << std::endl;
      }
      
      return false;
    }
    
    if (hasD2 && (d + CGAL::internal::my_sqrt(csd2(I, A)) > d2 + CGAL::internal::my_sqrt(csd2(v2, A))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,A| > d1 + |v2,A|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v2,A = " << CGAL::internal::my_sqrt(csd2(v2, A)) << std::endl;
        std::cout << "I,A = " << CGAL::internal::my_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::internal::my_sqrt(csd2(I, A))) << " vs. " << (d2 + CGAL::internal::my_sqrt(csd2(v2, A))) << std::endl;
      }
      
      return false;
    }
    
    if (hasD3 && (d + CGAL::internal::my_sqrt(csd2(I, A)) > d3 + CGAL::internal::my_sqrt(csd2(v3, A))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,A| > d1 + |v3,A|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v3,A = " << CGAL::internal::my_sqrt(csd2(v3, A)) << std::endl;
        std::cout << "I,A = " << CGAL::internal::my_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::internal::my_sqrt(csd2(I, A))) << " vs. " << (d3 + CGAL::internal::my_sqrt(csd2(v3, A))) << std::endl;
      }
      
      return false;
    }

    return true;
  }
    
  void expand_left_child(Cone_tree_node* cone, Segment_2 windowSegment)
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Flatten_triangle_3_along_segment_2 ft3as2(m_traits.flatten_triangle_3_along_segment_2_object());
    
    assert(cone->m_pendingLeftSubtree != NULL);
    
    cone->m_pendingLeftSubtree = NULL;
    
    if (window_distance_filter(cone, windowSegment, false))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->left_child_edge());
      Triangle_2 layoutFace = ft3as2(adjacentFace, 0, cone->left_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_traits, m_polyhedron, cone->left_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), cv2(windowSegment, 0), cv2(windowSegment, 1), Cone_tree_node::INTERVAL);
      node_created();
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
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Flatten_triangle_3_along_segment_2 ft3as2(m_traits.flatten_triangle_3_along_segment_2_object());
    
    assert(cone->m_pendingRightSubtree != NULL);
    
    cone->m_pendingRightSubtree = NULL;
    
    if (window_distance_filter(cone, windowSegment, true))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->right_child_edge());
      Triangle_2 layoutFace = ft3as2(adjacentFace, 0, cone->right_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_traits, m_polyhedron, cone->right_child_edge(), layoutFace, cone->source_image(), cone->distance_from_source_to_root(), cv2(windowSegment, 0), cv2(windowSegment, 1), Cone_tree_node::INTERVAL);
      node_created();
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
    typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
    typename Traits::Classify_barycentric_coordinate classify_barycentric_coordinate(m_traits.classify_barycentric_coordinate_object());
  
    size_t associatedEdge;
    CGAL::internal::Barycentric_coordinate_type type;
    boost::tie(type, associatedEdge) = classify_barycentric_coordinate(location);
    
    switch (type)
    {
      case CGAL::internal::BARYCENTRIC_COORDINATE_INTERNAL:
        expand_face_root(face, location);
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_EDGE:
        {
          halfedge_descriptor halfedge = CGAL::halfedge(face, m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            halfedge = CGAL::next(halfedge, m_polyhedron);
          }
          expand_edge_root(halfedge, cbcw(location, associatedEdge), cbcw(location, (associatedEdge + 1) % 3));
        }
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX:
        {
          halfedge_descriptor halfedge = CGAL::halfedge(face, m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            halfedge = CGAL::next(halfedge, m_polyhedron);
          }
          expand_vertex_root(CGAL::source(halfedge, m_polyhedron));
        }
        break;
      default:
        assert(false && "Invalid face location");
        // Perhaps hit an assertion that the type must not be external or invalid?
    }
  }
  
  void expand_face_root(face_descriptor faceId, Barycentric_coordinate faceLocation)
  {
    typename Traits::Project_triangle_3_to_triangle_2 pt3t2(m_traits.project_triangle_3_to_triangle_2_object());
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
  
    halfedge_descriptor start = CGAL::halfedge(faceId, m_polyhedron);
    halfedge_descriptor current = start;
    
    Cone_tree_node* faceRoot = new Cone_tree_node(m_traits, m_polyhedron, m_rootNodes.size());
    node_created();
    m_rootNodes.push_back(faceRoot);
    
    if (m_debugOutput)
    {
      typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
      std::cout << "\tFace Root Expansion: id = " << m_faceIndexMap[faceId] << " , Location = " << cbcw(faceLocation, 0) << " " << cbcw(faceLocation, 1) << " " << cbcw(faceLocation, 2) << " " << std::endl;
    }
    
    for (size_t currentVertex = 0; currentVertex < 3; ++currentVertex)
    {
      Triangle_3 face3d(triangle_from_halfedge(current));
      Triangle_2 layoutFace(pt3t2(face3d));
      Barycentric_coordinate rotatedFaceLocation(shifted_coordiate(faceLocation, currentVertex));
      Point_2 sourcePoint(construct_barycenter_in_triangle_2(layoutFace, rotatedFaceLocation));
      
      Cone_tree_node* child = new Cone_tree_node(m_traits, m_polyhedron, current, layoutFace, sourcePoint, FT(0.0), cv2(layoutFace, 0), cv2(layoutFace, 2), Cone_tree_node::FACE_SOURCE);
      node_created();
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
    typename Traits::Construct_barycenter_2 cb2(m_traits.construct_barycenter_2_object());
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Project_triangle_3_to_triangle_2 pt3t2(m_traits.project_triangle_3_to_triangle_2_object());
    typename Traits::Construct_triangle_2 ct2(m_traits.construct_triangle_2_object());
    
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
       layoutFaces[i] = pt3t2(faces3d[i]);
    }
    
    Point_2 sourcePoints[2];
    sourcePoints[0] = cb2(cv2(layoutFaces[0], 0), t0, cv2(layoutFaces[0], 1), t1); 
    sourcePoints[1] = cb2(cv2(layoutFaces[1], 0), t0, cv2(layoutFaces[1], 1), t1); 
    
    Cone_tree_node* edgeRoot = new Cone_tree_node(m_traits, m_polyhedron, m_rootNodes.size());
    node_created();
    m_rootNodes.push_back(edgeRoot);
    
    for (size_t side = 0; side < 2; ++side)
    {
      if (m_debugOutput)
      {
        std::cout << "\tExpanding edge root #" << side << " : " << std::endl;;
        std::cout << "\t\tFace = " << layoutFaces[side] << std::endl;
        std::cout << "\t\tLocation = " << sourcePoints[side] << std::endl;
      }
      
      Cone_tree_node* mainChild = new Cone_tree_node(m_traits, m_polyhedron, baseEdges[side], layoutFaces[side], sourcePoints[side], FT(0.0), cv2(layoutFaces[side], 0), cv2(layoutFaces[side], 2), Cone_tree_node::EDGE_SOURCE);
      node_created();
      edgeRoot->push_middle_child(mainChild);
      process_node(mainChild);

      Cone_tree_node* oppositeChild = new Cone_tree_node(m_traits, m_polyhedron, CGAL::prev(baseEdges[side], m_polyhedron), ct2(cv2(layoutFaces[side], 2), cv2(layoutFaces[side], 0), cv2(layoutFaces[side], 1)), sourcePoints[side], FT(0.0), cv2(layoutFaces[side], 2), cv2(layoutFaces[side], 1), Cone_tree_node::EDGE_SOURCE);
      node_created();
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

    Cone_tree_node* vertexRoot = new Cone_tree_node(m_traits, m_polyhedron, m_rootNodes.size(), CGAL::prev(CGAL::halfedge(vertex, m_polyhedron), m_polyhedron));

    node_created();
    m_rootNodes.push_back(vertexRoot);

    m_closestToVertices[m_vertexIndexMap[vertex]] = Node_distance_pair(vertexRoot, FT(0.0));

    expand_pseudo_source(vertexRoot);
  }

  void expand_pseudo_source(Cone_tree_node* parent)
  {
    typename Traits::Project_triangle_3_to_triangle_2 pt3t2(m_traits.project_triangle_3_to_triangle_2_object());
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    
    parent->m_pendingMiddleSubtree = NULL;
    
    vertex_descriptor expansionVertex = parent->target_vertex();
  
    halfedge_descriptor startEdge = CGAL::halfedge(expansionVertex, m_polyhedron);
    halfedge_descriptor currentEdge = CGAL::halfedge(expansionVertex, m_polyhedron);

    FT distanceFromTargetToRoot = parent->distance_from_target_to_root();

    if (m_debugOutput)
    {
      std::cout << "Distance from target to root: " << distanceFromTargetToRoot << std::endl;
    }
    
    // A potential optimization could be made by only expanding in the 'necessary' range (i.e. the range outside of geodesic visibility), but the
    // benefits may be small, since the node filtering would prevent more than one-level propagation.
    do
    {
      Triangle_3 face3d(triangle_from_halfedge(currentEdge));
      Triangle_2 layoutFace(pt3t2(face3d));

      if (m_debugOutput)
      {
        std::cout << "Expanding PsuedoSource: id = ";
        if (CGAL::face(currentEdge, m_polyhedron) != GraphTraits::null_face())
        {
          std::cout << m_faceIndexMap[CGAL::face(currentEdge, m_polyhedron)];
        }
        else
        {
          std::cout << "EXTERNAL";
        }
        std::cout << std::endl;
      }

      Cone_tree_node* child = new Cone_tree_node(m_traits, m_polyhedron, currentEdge, layoutFace, cv2(layoutFace, 1), distanceFromTargetToRoot, cv2(layoutFace, 0), cv2(layoutFace, 2), Cone_tree_node::VERTEX_SOURCE);

      node_created();
      parent->push_middle_child(child);
      process_node(child);
      
      currentEdge = CGAL::opposite(CGAL::next(currentEdge, m_polyhedron), m_polyhedron);
    }
    while (currentEdge != startEdge);

  }

  bool clip_to_bounds(const Segment_2& segment, const Ray_2& leftBoundary, const Ray_2& rightBoundary, Segment_2& outSegment)
  {
    typename Traits::Construct_source_2 cs2(m_traits.construct_source_2_object());
    typename Traits::Construct_segment_2 cseg2(m_traits.construct_segment_2_object());
    typename Traits::Construct_target_2 ct2(m_traits.construct_target_2_object());
    typename Traits::Construct_line_2 cl2(m_traits.construct_line_2_object());
    typename Traits::Intersect_2 i2(m_traits.intersect_2_object());
    typename Traits::Orientation_2 o2(m_traits.orientation_2_object());
    typename Traits::Construct_point_on_2 cpo2(m_traits.construct_point_on_2_object());
    typename Traits::Parametric_distance_along_segment_2 pdas2(m_traits.parametric_distance_along_segment_2_object());
  
    typedef typename cpp11::result_of<typename Traits::Intersect_2(Line_2, Line_2)>::type LineLineIntersectResult;

    Point_2 leftPoint;
    Point_2 rightPoint;
    
    if (m_debugOutput)
    {
      std::cout << "Clipping Segment " << segment << " with left = " << leftBoundary << " and right = " << rightBoundary << std::endl;
    }
    
    FT leftT;
    FT rightT;

    CGAL::Orientation leftOrientation = o2(cs2(leftBoundary), cpo2(leftBoundary, 1), cs2(segment));

    if (leftOrientation == CGAL::RIGHT_TURN || leftOrientation == CGAL::COLLINEAR)
    {
      if (m_debugOutput)
      {
        std::cout << "\tLeft is completely covered." << std::endl;
      }
      leftPoint = cs2(segment);
      leftT = FT(0.0);
    }
    else
    {
      LineLineIntersectResult cgalIntersection = i2(cl2(segment), cl2(leftBoundary));

      if (!cgalIntersection || !boost::get<Point_2>(&*cgalIntersection))
      {
        if (m_debugOutput)
        {
          std::cout << "Dropping left due to co-linearity of boundary. " << bool(cgalIntersection) << std::endl;
        }
        return false;
      }
      else
      {
        Point_2* result = boost::get<Point_2>(&*cgalIntersection);
        FT t0 = pdas2(cs2(segment), ct2(segment), *result);

        if (t0 >= FT(1.00000))
        {
          if (m_debugOutput)
          {
            std::cout << "Dropping due to missing left intersect. " << t0 << std::endl;
          }

          return false;
        }
        else if (t0 <= FT(0.00000))
        {
          if (m_debugOutput)
          {
            std::cout << "\tLeft is completely covered (secondary check). " << t0 << std::endl;
          }

          leftPoint = cs2(segment);
          leftT = FT(0.0);
        }
        else
        {
          if (m_debugOutput)
          {
            std::cout << "\tLeft intersects at t = " << t0 << std::endl;
          }
         
          leftPoint = *result;
          leftT = t0;
        }
      }
    }
    
    CGAL::Orientation rightOrientation = o2(cs2(rightBoundary), cpo2(rightBoundary, 1), ct2(segment));
    
    if (rightOrientation == CGAL::LEFT_TURN || rightOrientation == CGAL::COLLINEAR)
    {
      if (m_debugOutput)
      {
        std::cout << "Right is completely covered." << std::endl;
      }
      rightPoint = ct2(segment);
      rightT = FT(1.0);
    }
    else
    {
      LineLineIntersectResult cgalIntersection = i2(cl2(segment), cl2(rightBoundary));

      if (!cgalIntersection || !boost::get<Point_2>(&*cgalIntersection))
      {
        if (m_debugOutput)
        {
          std::cout << "Dropping due to co-linearity of right boundary." << std::endl;
        }
        return false;
      }
      else
      {
        Point_2* result = boost::get<Point_2>(&*cgalIntersection);
        FT t0 = pdas2(cs2(segment), ct2(segment), *result);
      
        if (t0 <= FT(0.00000))
        {
          if (m_debugOutput)
          {
            std::cout << "Dropping due to missing right intersect. " << t0 << std::endl;
          }
          return false;
        }
        else if (t0 >= FT(1.00000))
        {
          if (m_debugOutput)
          {
            std::cout << "\tRight is completely covered (secondary check). " << t0 << std::endl;
          }
          rightPoint = ct2(segment);
          rightT = FT(1.0);
        }
        else
        {
          if (m_debugOutput)
          {
            std::cout << "\tRight intersects at t = " << t0 << std::endl;
          }
          rightPoint = *result;
          rightT = t0;
        }
      }
    }
    
    if (leftT >= rightT)
    {
      if (m_debugOutput)
      {
        std::cout << "Dropping due to overlap. " << leftT << " : " << rightT << std::endl;
      }
      return false;
    }
    
    outSegment = cseg2(leftPoint, rightPoint);
    
    return true;
  }

  void process_node(Cone_tree_node* node)
  {
    typename Traits::Compare_relative_intersection_along_segment_2 crias2(m_traits.compare_relative_intersection_along_segment_2_object());
  
    bool leftSide = false;
    bool rightSide = false;
    
    if (!node->is_source_node())
    {
      leftSide = node->has_left_side();
      rightSide = node->has_right_side();
    }
    else
    {
      leftSide = true;
      rightSide = false;
    }

    bool propagateLeft = false;
    bool propagateRight = false;
    bool propagateMiddle = false;
  
    if (m_debugOutput)
    {
      std::cout << " Processing node " << node << " , level = " << node->level() << std::endl;
      std::cout << "\tFace = " << node->layout_face() << std::endl;
      std::cout << "\tVertices = ";
      halfedge_descriptor current = node->entry_edge();
      for (size_t i = 0; i < 3; ++i)
      {
        std::cout << m_vertexIndexMap[CGAL::source(current, m_polyhedron)] << " ";
        current = CGAL::next(current, m_polyhedron);
      }
      std::cout << std::endl;
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

      Node_distance_pair currentOccupier = m_vertexOccupiers[entryEdgeIndex];
      FT currentNodeDistance = node->distance_from_target_to_root();

      bool isLeftOfCurrent = false;
      
      if (m_debugOutput)
      {
        std::cout << "\t Entry Edge = " << entryEdgeIndex << std::endl;
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
          CGAL::Comparison_result comparison = crias2(
            node->entry_segment(), 
            node->ray_to_target_vertex().supporting_line(), 
            currentOccupier.first->entry_segment(),
            currentOccupier.first->ray_to_target_vertex().supporting_line()
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
        
        // This is a consequence of using the same basic node type for source and interval nodes
        // If this is a source node, it is only pointing to one of the two opposite edges (the left one by convention)
        if (node->node_type() != Cone_tree_node::INTERVAL)
        {
          propagateRight = false;
          
          // Propagating a pseudo-source on a boundary vertex can result in a cone on a null face
          // In such a case, we only care about the part of the cone pointing at the vertex (i.e. the middle child),
          // so we can avoid propagating over the (non-existant) left opposite edge
          if (node->is_null_face())
          {
            propagateLeft = false;
          }
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
        Node_distance_pair currentClosest = m_closestToVertices[targetVertexIndex];
        
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

          m_closestToVertices[targetVertexIndex] = Node_distance_pair(node, currentNodeDistance);
        }
      }
      else
      {
        if (isLeftOfCurrent)
        {
          propagateLeft = true;
        }
        else if (!node->is_source_node())
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
    
    if (node->level() <= num_faces(m_polyhedron))
    {
      if (propagateLeft)
      {
        push_left_child(node);
      }
      
      if (propagateRight && !node->is_source_node())
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
    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());
  
    if (CGAL::face(parent->left_child_edge(), m_polyhedron) != GraphTraits::null_face())
    {
      Segment_2 leftWindow;
      
      if (parent->is_source_node())
      {
        leftWindow = parent->left_child_base_segment();
      }
      else
      {
        bool result = clip_to_bounds(parent->left_child_base_segment(), parent->left_boundary(), parent->right_boundary(), leftWindow);
        if (!result)
        {
          if (m_debugOutput)
          {
            std::cout << "Left child clip failed, killing node." << std::endl;
          }
          return;
        }
      }
      
      FT distanceEstimate = parent->distance_from_source_to_root() + CGAL::internal::my_sqrt(csd2(parent->source_image(), leftWindow));

      if (m_debugOutput)
      {
        std::cout << "\tPushing Left Child, Segment = " << parent->left_child_base_segment() << " , clipped = " << leftWindow << " , Estimate = " << distanceEstimate << std::endl;
      }

      Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::LEFT_CHILD, leftWindow);
      parent->m_pendingLeftSubtree = event;
      
      m_expansionPriqueue.push(event);
      
      queue_pushed();
    }
  }

  void push_right_child(Cone_tree_node* parent)
  {
    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());
    
    if (CGAL::face(parent->right_child_edge(), m_polyhedron) != GraphTraits::null_face())
    {
      Segment_2 rightWindow;
      bool result = clip_to_bounds(parent->right_child_base_segment(), parent->left_boundary(), parent->right_boundary(), rightWindow);
      
      if (!result)
      {
        if (m_debugOutput)
        {
          std::cout << "Right child clip failed, killing node." << std::endl;
        }
        return;
      }

      FT distanceEstimate = parent->distance_from_source_to_root() + CGAL::internal::my_sqrt(csd2(parent->source_image(), rightWindow));
      
      if (m_debugOutput)
      {
        std::cout << "\tPushing Right Child, Segment = " << parent->right_child_base_segment() << " , clipped = " << rightWindow << " , Estimate = " << distanceEstimate << std::endl;
      }
      
      Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate, Cone_expansion_event::RIGHT_CHILD, rightWindow);
      queue_pushed();
      parent->m_pendingRightSubtree = event;

      m_expansionPriqueue.push(event);
    }
  }

  void push_middle_child(Cone_tree_node* parent)
  {
    if (m_debugOutput)
    {
      std::cout << "\tPushing Middle Child, Estimate = " << parent->distance_from_target_to_root() << std::endl;
    }

    Cone_expansion_event* event = new Cone_expansion_event(parent, parent->distance_from_target_to_root(), Cone_expansion_event::PSEUDO_SOURCE);

    queue_pushed();

    parent->m_pendingMiddleSubtree = event;
    
    m_expansionPriqueue.push(event);
  }
  
  void delete_node(Cone_tree_node* node, bool destruction = false)
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
      
        delete_node(node->remove_left_child(), destruction);
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
        
        delete_node(node->remove_right_child(), destruction);
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
        delete_node(node->pop_middle_child(), destruction);
      }
      
      // At the point of destruction, the polyhedron referenced may have gone out of scope, we wish to distinguish between deletion with an assumed reference
      // to the original polyhedron, and deletion without
      if (!node->is_root_node() && !destruction)
      {
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
      
      delete node;
    }
    
    node_deleted();
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
    }
  }
  
  bool is_saddle_vertex(vertex_descriptor v)
  {
    return m_traits.is_saddle_vertex_object()(v, m_polyhedron, m_vertexPointMap);
  }
  
  bool is_boundary_vertex(vertex_descriptor v)
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
  
  void delete_all_nodes()
  {
    for (size_t i = 0; i < m_rootNodes.size(); ++i)
    {
      delete_node(m_rootNodes[i], true);
    }
  }
  
  void reset_algorithm(bool clearFaceLocations = true)
  {
    m_closestToVertices.resize(boost::num_vertices(m_polyhedron));
    std::fill(m_closestToVertices.begin(), m_closestToVertices.end(), Node_distance_pair(NULL, FT(0.0)));
    m_vertexOccupiers.resize(CGAL::num_halfedges(m_polyhedron));
    std::fill(m_vertexOccupiers.begin(), m_vertexOccupiers.end(), Node_distance_pair(NULL, FT(0.0)));

    while (!m_expansionPriqueue.empty())
    {
      delete m_expansionPriqueue.top();
      m_expansionPriqueue.pop();
    }

    if (clearFaceLocations)
    {
      m_faceLocations.clear();
    }
    
    delete_all_nodes();
    m_rootNodes.clear();
    m_vertexIsPseudoSource.resize(boost::num_vertices(m_polyhedron));

#if !defined(NDEBUG)
    m_currentNodeCount = 0;
    m_peakNodeCount = 0;
    m_queueAtPeakNodes = 0;
    m_peakQueueSize = 0;
    m_nodesAtPeakQueue = 0;
#endif

  }
  
  template <class Visitor>
  void visit_shortest_path(Cone_tree_node* startNode, const Point_2& startLocation, Visitor& visitor)
  {
    typename Traits::Parametric_distance_along_segment_2 parametric_distance_along_segment_2(m_traits.parametric_distance_along_segment_2_object());
    typename Traits::Construct_ray_2 construct_ray_2(m_traits.construct_ray_2_object());
    typename Traits::Construct_segment_2 construct_segment_2(m_traits.construct_segment_2_object());
    typename Traits::Construct_line_2 construct_line_2(m_traits.construct_line_2_object());
    typename Traits::Construct_source_2 construct_source_2(m_traits.construct_source_2_object());
    typename Traits::Construct_target_2 construct_target_2(m_traits.construct_target_2_object());
    typename Traits::Intersect_2 intersect_2(m_traits.intersect_2_object());
    
    typedef typename cpp11::result_of<typename Traits::Intersect_2 (Line_2, Line_2)>::type LineLineIntersectResult;
    
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
          Ray_2 rayToLocation(construct_ray_2(current->source_image(), currentLocation));
          
          LineLineIntersectResult cgalIntersection = intersect_2(construct_line_2(entrySegment), construct_line_2(rayToLocation));

          assert(cgalIntersection);
          
          Point_2* result = boost::get<Point_2>(&*cgalIntersection);
          
          assert(result && "Error, did not get point intersection on path walk to source");
          
          FT t0 = parametric_distance_along_segment_2(construct_source_2(entrySegment), construct_target_2(entrySegment), *result);
          
          if (m_debugOutput)
          {
            std::cout << "Current Node: " << current << " , Face = " << current->layout_face() << std::endl;
            halfedge_descriptor halfedge = current->entry_edge();
            std::cout << "Face vertices: ";
            for (size_t i = 0; i < 3; ++i)
            {
              std::cout << m_vertexIndexMap[CGAL::source(halfedge, m_polyhedron)] << ",";
              halfedge = CGAL::next(halfedge, m_polyhedron);
            }
            std::cout << std::endl;
            std::cout << "Current Location: " << currentLocation << std::endl;
            std::cout << "Distance: " << current->distance_to_root(currentLocation) << std::endl;
            std::cout << "Inside cone: " << (current->inside_window(currentLocation) ? "Yes" : "No") << std::endl;
            std::cout << "Current Source: " << current->source_image() << std::endl;
            std::cout << "Current Segment: " << entrySegment << std::endl;
            std::cout << "Current Left Window: " << current->window_left() << "  ,  " << m_traits.parametric_distance_along_segment_2_object()(entrySegment.start(), entrySegment.end(), current->window_left()) << std::endl;
            std::cout << "Current Right Window: " << current->window_right() << "  ,  " << m_traits.parametric_distance_along_segment_2_object()(entrySegment.start(), entrySegment.end(), current->window_right()) << std::endl;
            std::cout << "Current Segment Intersection: " << *result << std::endl;
            std::cout << "Edge: (" << m_vertexIndexMap[CGAL::source(current->entry_edge(), m_polyhedron)] << "," << m_vertexIndexMap[CGAL::target(current->entry_edge(), m_polyhedron)] << ")  :  " << t0 << std::endl;
          }
          
          visitor.edge(current->entry_edge(), t0);

          if (current->is_left_child())
          {
            Segment_2 baseSegment = current->parent()->left_child_base_segment();
            currentLocation = *result;
          }
          else if (current->is_right_child())
          {
            Segment_2 baseSegment = current->parent()->right_child_base_segment();
            currentLocation = *result;
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
        default:
          assert(false && "Unhandled node type found in tree");
      }
    }
  }
  
  void add_to_face_list(Cone_tree_node* node)
  {
    if (!node->is_root_node() && !node->is_null_face())
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
    
    for (size_t i = 0; i < node->num_middle_children(); ++i)
    {
      add_to_face_list(node->get_middle_child(i));
    }
  }
  
  Point_2 face_location_with_normalized_coordinate(Cone_tree_node* node, Barycentric_coordinate location)
  {
    return construct_barycenter_in_triangle_2(node->layout_face(), localized_coordiate(node, location));
  }
  
  Barycentric_coordinate localized_coordiate(Cone_tree_node* node, Barycentric_coordinate location)
  {
    return shifted_coordiate(location, node->edge_face_index());
  }
  
  Barycentric_coordinate shifted_coordiate(Barycentric_coordinate location, size_t shift)
  {
    typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
    typename Traits::Construct_barycentric_coordinate cbc(m_traits.construct_barycentric_coordinate_object());
    return cbc(cbcw(location, shift), cbcw(location, (shift + 1) % 3), cbcw(location, (shift + 2) % 3));
  }
  
  std::pair<Node_distance_pair, Barycentric_coordinate> nearest_on_face(face_descriptor face, Barycentric_coordinate location)
  {
    size_t faceIndex = m_faceIndexMap[face];
    
    halfedge_descriptor halfedge = CGAL::halfedge(face, m_polyhedron);

    Cone_tree_node* closest = NULL;
    FT closestDistance;
    
    std::vector<Cone_tree_node*>& currentFaceList = m_faceOccupiers[faceIndex];
    
    for (size_t i = 0; i < currentFaceList.size(); ++i)
    {
      Cone_tree_node* current = currentFaceList[i];
      
      if (closest != NULL && current->distance_from_source_to_root() >= closestDistance)
      {
        continue;
      }
      
      Point_2 locationInContext = face_location_with_normalized_coordinate(current, location);

      if (current->inside_window(locationInContext))
      {
        FT currentDistance = current->distance_to_root(locationInContext);
        
        if (closest == NULL || currentDistance < closestDistance)
        {
          closest = current;
          closestDistance = currentDistance;
        }
      }
    }
  
    if (closest)
    {
      return std::make_pair(Node_distance_pair(closest, closestDistance), localized_coordiate(closest, location));
    }
    else
    {
      return std::make_pair(Node_distance_pair(NULL, FT(0.0)), Barycentric_coordinate(FT(0.0), FT(0.0), FT(0.0)));
    }
  }
  
  std::pair<Node_distance_pair, Barycentric_coordinate> nearest_to_location(face_descriptor face, Barycentric_coordinate location)
  {
    typename Traits::Construct_barycentric_coordinate_weight cbcw(m_traits.construct_barycentric_coordinate_weight_object());
    typename Traits::Construct_barycentric_coordinate cbc(m_traits.construct_barycentric_coordinate_object());
    typename Traits::Classify_barycentric_coordinate classify_barycentric_coordinate(m_traits.classify_barycentric_coordinate_object());
    
    size_t associatedEdge;
    CGAL::internal::Barycentric_coordinate_type type;
    boost::tie(type, associatedEdge) = classify_barycentric_coordinate(location);
    
    switch (type)
    {
      case CGAL::internal::BARYCENTRIC_COORDINATE_INTERNAL:
        return nearest_on_face(face, location);
      case CGAL::internal::BARYCENTRIC_COORDINATE_EDGE:
        {
          halfedge_descriptor halfedge = CGAL::halfedge(face, m_polyhedron);
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            halfedge = CGAL::next(halfedge, m_polyhedron);
          }
          expand_edge_root(halfedge, cbcw(location, associatedEdge), cbcw(location, (associatedEdge + 1) % 3));
          
          halfedge_descriptor oppositeHalfedge = CGAL::opposite(halfedge, m_polyhedron);
          
          size_t oppositeIndex = internal::edge_index(oppositeHalfedge, m_polyhedron);
          
          FT oppositeLocationCoords[3] = { FT(0.0), FT(0.0), FT(0.0) };
          
          oppositeLocationCoords[oppositeIndex] = cbcw(location, (associatedEdge + 1) % 3);
          oppositeLocationCoords[(oppositeIndex + 1) % 3] = cbcw(location, associatedEdge);

          std::pair<Node_distance_pair,Barycentric_coordinate> mainFace = nearest_on_face(face, location);
          Barycentric_coordinate oppositeLocation(cbc(oppositeLocationCoords[0], oppositeLocationCoords[1], oppositeLocationCoords[2]));
          std::pair<Node_distance_pair,Barycentric_coordinate> otherFace = nearest_on_face(CGAL::face(oppositeHalfedge, m_polyhedron), oppositeLocation);
          
          if (mainFace.first.first == NULL)
          {
            return otherFace;
          }
          else if (otherFace.first.first == NULL)
          {
            return mainFace;
          }
          else
          {
            return mainFace.first.second < otherFace.first.second ? mainFace : otherFace;
          }
        }
        break;
      case CGAL::internal::BARYCENTRIC_COORDINATE_VERTEX:
        {
          halfedge_descriptor halfedge = CGAL::halfedge(face, m_polyhedron);
          
          for (size_t i = 0; i < associatedEdge; ++i)
          {
            halfedge = CGAL::next(halfedge, m_polyhedron);
          }
          
          vertex_descriptor vertex = CGAL::source(halfedge, m_polyhedron);

          return std::make_pair(m_closestToVertices[m_vertexIndexMap[vertex]], Barycentric_coordinate(FT(0.0), FT(0.0), FT(1.0)));
        }
        break;
        
      default:
        assert(false && "Invalid face location");
    }
  }
  
  static bool cone_comparator(const Cone_tree_node* lhs, const Cone_tree_node* rhs)
  {
    return lhs->distance_from_source_to_root() < rhs->distance_from_source_to_root();
  }
  
  void construct_sequence_tree_internal()
  {
    reset_algorithm(false);
    set_vertex_types();
    
    m_vertexOccupiers.resize(CGAL::num_halfedges(m_polyhedron));
    m_closestToVertices.resize(CGAL::num_vertices(m_polyhedron));

    if (m_debugOutput)
    {
      vertex_iterator current, end;
      
      size_t numVertices = 0;

      for (boost::tie(current,end) = boost::vertices(m_polyhedron); current != end; ++current)
      {
        std::cout << "Vertex#" << numVertices << ": p = " << m_vertexPointMap[*current] << " , Saddle Vertex: " << (is_saddle_vertex(*current) ? "yes" : "no") << " , Boundary Vertex: " << (is_boundary_vertex(*current) ? "yes" : "no") << std::endl;
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
    
    for (size_t i = 0; i < m_faceLocations.size(); ++i)
    {
      if (m_debugOutput)
      {
        std::cout << "Root: " << m_faceIndexMap[m_faceLocations[i].first] << " , " << m_faceLocations[i].second[0] << " " << m_faceLocations[i].second[1] << " " << m_faceLocations[i].second[2] << " " << std::endl;
      }
      
      expand_root(m_faceLocations[i].first, m_faceLocations[i].second);
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
            
            expand_pseudo_source(parent);
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
    
    m_faceOccupiers.clear();
    m_faceOccupiers.resize(CGAL::num_faces(m_polyhedron));
    
    for (size_t i = 0; i < m_rootNodes.size(); ++i)
    {
      add_to_face_list(m_rootNodes[i]);
    }
    
    for (size_t i = 0; i < m_faceOccupiers.size(); ++i)
    {
      std::vector<Cone_tree_node*>& currentFaceList = m_faceOccupiers[i];
      std::sort(currentFaceList.begin(), currentFaceList.end(), cone_comparator);
    }
    
    if (m_debugOutput)
    {   
      std::cout << "Closest distances: " << std::endl;
      
      for (size_t i = 0; i < m_closestToVertices.size(); ++i)
      {
        std::cout << "\tVertex = " << i << std::endl;
        std::cout << "\tDistance = " << m_closestToVertices[i].second << std::endl;
      }
      
      std::cout << std::endl;
      
      for (size_t i = 0; i < m_faceOccupiers.size(); ++i)
      {
        std::cout << "\tFace = " << i << std::endl;
        std::cout << "\t#Occupiers = " << m_faceOccupiers[i].size() << std::endl;
      }
      
      std::cout << std::endl << "Done!" << std::endl;
    }
  }
  
public:
  
  /// \name Constructors
  /// @{
  
  /*!
  \brief Creates a shortest paths object associated with a specific polyhedron.
  
  \details No copy of the polyhedron is made, only a reference to the polyhedron is held.
  Default versions of the necessary polyhedron property maps are created and
  used with this constructor.
  
  \param polyhedron The polyhedral surface to use.  Note that it must be triangulated.
  
  \param traits An optional instance of the traits class to use.
  
  */
  Polyhedron_shortest_path(Polyhedron& polyhedron, const Traits& traits = Traits())
    : m_traits(traits)
    , m_polyhedron(polyhedron)
    , m_vertexIndexMap(CGAL::get(boost::vertex_external_index, polyhedron))
    , m_halfedgeIndexMap(CGAL::get(CGAL::halfedge_external_index, polyhedron))
    , m_faceIndexMap(CGAL::get(CGAL::face_external_index, polyhedron))
    , m_vertexPointMap(CGAL::get(CGAL::vertex_point, polyhedron))
    , m_debugOutput(false)
  {
    reset_algorithm();
  }
  
  /*!
  \brief Creates a shortest paths object associated with a specific polyhedron.
  
  \details No copy of the polyhedron is made, only a reference to the polyhedron is held.
  
  \param polyhedron The polyhedral surface to use.  Note that it must be triangulated.
  
  \param vertexIndexMap Maps between vertices and their index.
  
  \param halfedgeIndexMap Maps between halfedges and their index.
  
  \param faceIndexMap Maps between faces and their index.
  
  \param vertexPointMap Maps between vertices and their 3-dimensional coordinates.
  
  \param traits An optional instance of the traits class to use.
  */
  Polyhedron_shortest_path(Polyhedron& polyhedron, VertexIndexMap vertexIndexMap, HalfedgeIndexMap halfedgeIndexMap, FaceIndexMap faceIndexMap, VertexPointMap vertexPointMap, const Traits& traits = Traits())
    : m_traits(traits)
    , m_polyhedron(polyhedron)
    , m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
    , m_vertexPointMap(vertexPointMap)
    , m_debugOutput(false)
  {
    reset_algorithm();
  }
  
  /// @}

  ~Polyhedron_shortest_path()
  {
    delete_all_nodes();
    
#if !defined(NDEBUG)
    if (m_debugOutput)
    {
      std::cout << "Final node count: " << m_currentNodeCount << std::endl;
    }
    return;
    assert(m_currentNodeCount == 0);
#endif
  }
  
  /// \name Methods
  /// @{
  
  /*!
  \brief Compute shortest paths sequence tree from a single vertex
  
  \details Constructs a shortest paths sequence tree that covers shortest surface paths
  to all locations on the polyhedron from the given source vertex.
  
  \param face Handle to the face on which the source originates.
  
  \param location Barycentric coordinate on face specifying the source location.
  */
  void construct_sequence_tree(vertex_descriptor vertex)
  {
    m_faceLocations.clear();
    m_faceLocations.push_back(get_vertex_as_face_location(vertex));
    construct_sequence_tree_internal();
  }
  
  /*!
  \brief Compute shortest paths from a single source location
  
  \details Constructs a shortest paths sequence tree that covers shortest surface paths
  to all locations on the polyhedron reachable from the given source point.
  
  \param face Handle to the face on which the source originates.
  
  \param location Barycentric coordinate on face specifying the source location.
  */
  void construct_sequence_tree(face_descriptor face, Barycentric_coordinate location)
  {
    m_faceLocations.clear();
    m_faceLocations.push_back(std::make_pair(face, location));
    construct_sequence_tree_internal();
  }
  
  /*!
  \brief Compute a shortest path sequence tree from multiple source vertices
  
  \details Constructs a shortest paths sequence tree that covers shortest paths
  to all locations on the polyhedron reachable from the supplied source locations.
  
  \tparam InputIterator a ForwardIterator type which dereferences to Face_location.
  
  \param faceLocationsBegin iterator to the first in the list of face location pairs.
  
  \param faceLocationsEnd iterator to one past the end of the list of face location pairs.
  */
  template <class InputIterator>
  typename boost::enable_if<typename boost::is_same<typename std::iterator_traits<InputIterator>::value_type, vertex_descriptor>::value, void>::type construct_sequence_tree(InputIterator begin, InputIterator end)
  {
    m_faceLocations.clear();
    for (InputIterator it = begin; it != end; ++it)
    {
      m_faceLocations.push_back(get_vertex_as_face_location(*it));
    }
    construct_sequence_tree_internal();
  }
  
    /*!
  \brief Compute a shortest path sequence tree from multiple source locations
  
  \details Constructs a shortest paths sequence tree that covers shortest surface paths
  to all locations on the polyhedron reachable from the supplied source locations.
  
  \tparam InputIterator a ForwardIterator type which dereferences to Face_location.
  
  \param faceLocationsBegin iterator to the first in the list of face location pairs.
  
  \param faceLocationsEnd iterator to one past the end of the list of face location pairs.
  */
  template <class InputIterator>
  typename boost::enable_if<typename boost::is_same<typename std::iterator_traits<InputIterator>::value_type, Face_location>::type, void>::type construct_sequence_tree(InputIterator begin, InputIterator end)
  {
    m_faceLocations.clear();
    for (InputIterator it = begin; it != end; ++it)
    {
      m_faceLocations.push_back(*it);
    }
    construct_sequence_tree_internal();
  }
  
  /*!
  \brief Gets the face location of the `i`th source point
    given to this algorithm.
    
  \param i Index of the source point to get.  Precondition: 0 <= i < num_source_locations()
  */
  const Face_location& get_source_location(size_t i) const
  {
    return m_faceLocations[i];
  }
  
  /*!
  \brief The total number of source points in the current
    sequence tree, or 0 if no sequence tree is computed.
  */
  size_t num_source_locations() const
  {
    return m_faceLocations.size();
  }
  
  /*!
  Computes the shortest surface distance from a vertex to any source point
  
  \param v The vertex to act as the query point
  
  \return A pair, containing the distance to the source location, and the
    index of the source location itself.  If no source location was 
    reachable, the distance will be a negative value and the source 
    location will be an index greater than the number of source points.
  */
  std::pair<FT, size_t> shortest_distance_to_vertex(vertex_descriptor v)
  {
    Node_distance_pair result = m_closestToVertices[m_vertexIndexMap[v]];
    
    Cone_tree_node* current = result.first;
    
    if (current)
    {
      return std::make_pair(result.second, current->tree_id());
    }
    else
    {
      return std::make_pair(FT(-1.0), num_source_locations());
    }
  }
  
  /*!
  \brief Computes the shortest surface distance from any surface location to any source point
  
  \param face Face of the polyhedron of the query point
  
  \param location Barycentric coordinate on face of the query point
  
  \return A pair, containing the distance to the source location, and the
    index of the source location itself.  If no source location was 
    reachable, the distance will be a negative value and the source 
    location will be an index greater than the number of source points.
  */
  std::pair<FT, size_t> shortest_distance_to_location(face_descriptor face, Barycentric_coordinate location)
  {
    std::pair<Node_distance_pair, Barycentric_coordinate> result = nearest_to_location(face, location);
    
    Cone_tree_node* current = result.first.first;
    
    if (current)
    {
      return std::make_pair(result.first.second, current->tree_id());
    }
    else
    {
      return std::make_pair(FT(-1.0), num_source_locations());
    }
  }
  
  /*!
  \brief Visits the sequence of edges, vertices and faces traversed by the shortest path
  from a vertex to any source point.
  \param v The vertex to act as the query point
  \param visitor A model of PolyhedronShortestPathVisitor to receive the shortest path
  \return true if there exists a shortest path to v, false otherwise (may occur if the face graph is disconnected)
  */
  template <class Visitor>
  bool shortest_path_sequence(vertex_descriptor v, Visitor& visitor)
  {
    Cone_tree_node* current = m_closestToVertices[m_vertexIndexMap[v]].first;
    
    if (current)
    {
      visit_shortest_path(current, current->target_vertex_location(), visitor);
      return true;
    }
    else
    {
      return false;
    }
  }
  
  /*!
  \brief Visits the sequence of edges, vertices and faces traversed by the shortest path
  from any surface location to any source point.
  
  \param face Face of the polyhedron of the query point
  
  \param location Barycentric coordinate on face of the query point
  
  \param visitor A model of PolyhedronShortestPathVisitor to receive the shortest path
  */
  template <class Visitor>
  bool shortest_path_sequence(face_descriptor face, Barycentric_coordinate location, Visitor& visitor)
  {
    std::pair<Node_distance_pair, Barycentric_coordinate> result = nearest_to_location(face, location);
    Cone_tree_node* current = result.first.first;
    
    if (current)
    {
      Point_2 locationInContext = construct_barycenter_in_triangle_2(current->layout_face(), result.second);
      visit_shortest_path(current, locationInContext, visitor);
      return true;
    }
    else
    {
      return false;
    }
  }

  /*!
  \brief Visits the sequence of points in the surface-restricted polyline from a vertex
  to any source point (used for visualization of the shortest path).
  
  \param v The vertex to act as the query point
  
  \param output An OutputIterator to receive the shortest path points as Point_3
  */
  template <class OutputIterator>
  void shortest_path_points(vertex_descriptor v, OutputIterator output)
  {
    *output = get_vertex_location(v);
    ++output;
    Point_path_visitor_wrapper<OutputIterator> wrapper(*this, output);
    shortest_path_sequence(v, wrapper);
  }
  
  /*!
  \brief Visits the sequence of points in the surface-restricted polyline from any surface location
  to any source point (used for visualization of the shortest path).
 
  \param face Face of the polyhedron of the query point
  
  \param location Barycentric coordinate on face of the query point
  
  \param output An OutputIterator to receive the shortest path points as Point_3
  */
  template <class OutputIterator>
  void shortest_path_points(face_descriptor face, Barycentric_coordinate location, OutputIterator output)
  {
    *output = get_face_location(face, location);
    ++output;
    Point_path_visitor_wrapper<OutputIterator> wrapper(*this, output);
    shortest_path_sequence(face, location, wrapper);
  }
  
  /*!
  \brief Returns the 3-dimensional coordinate of the given face and face location on the polyhedron.
  
  \param face Face of the polyhedron of the query point
  
  \param location Barycentric coordinate on face of the query point
  */
  Point_3 get_face_location(face_descriptor face, Barycentric_coordinate location) const
  {
    return construct_barycenter_in_triangle_3(triangle_from_face(face), location);
  }
  
  /*!
  \brief Returns the 3-dimensional coordinate of the given edge and a parametric location along that edge.
  
  \param edge Edge of the polyhedron to use
  
  \param t Parametric distance along edge
  */
  Point_3 get_edge_location(halfedge_descriptor edge, FT t) const
  {
    typename Traits::Construct_barycenter_3 construct_barycenter_3(m_traits.construct_barycenter_3_object());
    
    // Note: the parameter t is meant to be the weighted coordinate on the _endpoint_ (i.e. target) of the segment
    return construct_barycenter_3(m_vertexPointMap[CGAL::target(edge, m_polyhedron)], t, m_vertexPointMap[CGAL::source(edge, m_polyhedron)]);
  }
  
  /*!
  \brief Returns the 3-dimensional of the given vertex.
  
  \param vertex Vertex of the polyhedron
  */
  Point_3 get_vertex_location(vertex_descriptor vertex) const
  {
    return m_vertexPointMap[vertex];
  }
  
  /*!
  \brief Returns a vertex location as a face location object
  
  \param vertex Vertex of the polyhedron
  */
  Face_location get_vertex_as_face_location(vertex_descriptor vertex) const
  {
    typename Traits::Construct_barycentric_coordinate construct_barycentric_coordinate(m_traits.construct_barycentric_coordinate_object());
    halfedge_descriptor he = CGAL::next(CGAL::halfedge(vertex, m_polyhedron), m_polyhedron);
    face_descriptor locationFace = CGAL::face(he, m_polyhedron);
    size_t edgeIndex = CGAL::internal::edge_index(he, m_polyhedron);
    
    FT coords[3] = { FT(0.0), FT(0.0), FT(0.0) };
    
    coords[edgeIndex] = FT(1.0);
    
    return Face_location(locationFace, construct_barycentric_coordinate(coords[0], coords[1], coords[2]));
  }
  
  /*!
  \brief Returns an edge location as a face location pair
  
  \param he halfedge of the polyhedron
  \param t parametric distance along he
  */
  Face_location get_edge_as_face_location(halfedge_descriptor he, FT t) const
  {
    typename Traits::Construct_barycentric_coordinate cbc(m_traits.construct_barycentric_coordinate_object());
    face_descriptor locationFace = CGAL::face(he, m_polyhedron);
    size_t edgeIndex = CGAL::internal::edge_index(he, m_polyhedron);
    
    const FT oneMinusT(FT(1.0) - t);
    
    FT coords[3] = { FT(0.0), FT(0.0), FT(0.0) };
    
    coords[edgeIndex] = oneMinusT;
    coords[(edgeIndex + 1) % 3] = t;

    return Face_location(locationFace, cbc(coords[0], coords[1], coords[2]));
  }
  
  /*!
  \brief Return the nearest face location to the given point.
  
  \param point Point to locate on the polyhedron
  \param tree A cached AABB to perform the point location
  */
  Face_location get_nearest_face_location(const Point_3& location, const AABB_polyhedron_tree& tree) const
  {
    typename Traits::Construct_barycentric_coordinate_in_triangle_3 cbcit3(m_traits.construct_barycentric_coordinate_in_triangle_3_object());
    typename AABB_polyhedron_tree::Point_and_primitive_id result = tree.closest_point_and_primitive(location);
    
    face_descriptor face = result.second;
    Barycentric_coordinate b = cbcit3(triangle_from_face(face), result.first);
    return Face_location(face, b);
  }
  
  /*!
  \brief Return the nearest face location to the given point.
    Note that this will fully build an AABB on each call, use the
    other version in conjunction with `construct_aabb_tree' 
    if you need to call this method more than once.
  
  \param point Point to locate on the polyhedron
  */
  Face_location get_nearest_face_location(const Point_3& location) const
  {
    AABB_polyhedron_tree tree;
    construct_aabb_tree(tree);
    return get_nearest_face_location(location, tree);
  }
  
  /*!
  \brief Return the face location along `ray` nearest to
    its source point.
  
  \param ray Ray to intersect with the polyhedron
  \param tree A cached AABB to perform the intersection
  */
  Face_location get_nearest_face_location(const Ray_3& ray, const AABB_polyhedron_tree& tree) const
  {
    typename Traits::Construct_barycentric_coordinate_in_triangle_3 cbcit3(m_traits.construct_barycentric_coordinate_in_triangle_3_object());
    typename Traits::Compute_squared_distance_3 csd3(m_traits.compute_squared_distance_3_object());
    typedef typename AABB_polyhedron_traits::template Intersection_and_primitive_id<Ray_3>::Type Intersection_type;
    typedef boost::optional<Intersection_type> Ray_intersection;
    
    std::vector<Ray_intersection> intersections;
    
    tree.all_intersections(ray, std::back_inserter(intersections));
    
    bool foundOne = false;
    FT nearestDistance;
    Point_3 nearestPoint;
    face_descriptor nearestFace;
    
    for (size_t i = 0; i < intersections.size(); ++i)
    {
      if (intersections[i])
      {
        Point_3* intersectionPoint = boost::get<Point_3>(&(intersections[i]->first));
        
        if (intersectionPoint)
        {
          FT distance = csd3(*intersectionPoint, ray.source());
          
          if (!foundOne || distance < nearestDistance)
          {
            foundOne = true;
            nearestPoint = *intersectionPoint;
            nearestDistance = distance;
            nearestFace = intersections[i]->second;
          }
        }
      }
    }
    
    if (foundOne)
    {
      Barycentric_coordinate b = cbcit3(triangle_from_face(nearestFace), nearestPoint);
      return Face_location(nearestFace, b);
    }
    else
    {
      return Face_location(GraphTraits::null_face(), Barycentric_coordinate());
    }
  }
  
  /*!
  \brief Return the face location along `ray` nearest to
    its source point.
    Note that this will fully build an AABB on each call, use the
    other version in conjunction with `construct_aabb_tree' 
    if you need to call this method more than once.
  
  \param ray Ray to intersect with the polyhedron
  */
  Face_location get_nearest_face_location(const Ray_3& ray) const
  {
    AABB_polyhedron_tree tree;
    construct_aabb_tree(tree);
    return get_nearest_face_location(ray, tree);
  }
  
  /*!
  \brief Creates an AABB tree suitable for use with `get_nearest_face_location`.
  
  \param outTree Output parameter to hold the created AABB tree
  */
  void construct_aabb_tree(AABB_polyhedron_tree& outTree) const
  {
    face_iterator facesStart, facesEnd;
    boost::tie(facesStart, facesEnd) = CGAL::faces(m_polyhedron);
    outTree.rebuild(facesStart, facesEnd, m_polyhedron);
    outTree.build();
  }
  
/// @}

};

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_POLYHEDRON_SHORTEST_PATH_H
