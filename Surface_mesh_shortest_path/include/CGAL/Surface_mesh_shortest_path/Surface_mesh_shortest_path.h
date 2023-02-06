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

#ifndef CGAL_SURFACE_MESH_SHORTEST_PATH_SURFACE_MESH_SHORTEST_PATH_H
#define CGAL_SURFACE_MESH_SHORTEST_PATH_SURFACE_MESH_SHORTEST_PATH_H

#include <CGAL/license/Surface_mesh_shortest_path.h>

#include <CGAL/Surface_mesh_shortest_path/barycentric.h>
#include <CGAL/Surface_mesh_shortest_path/internal/Cone_tree.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/assertions.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Default.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <boost/lexical_cast.hpp>
#include <boost/variant/get.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <list>
#include <queue>
#include <utility>
#include <vector>
#include <type_traits>

namespace CGAL {

/*!
\ingroup PkgSurfaceMeshShortestPathRef

\brief Computes shortest surface paths from one or more source points on a surface mesh.

\details Uses an optimized variation of Chen and Han's \f$ O(n^2) \f$ algorithm by Xin and Wang.
Refer to those respective papers for the details of the implementation.

\tparam Traits a model of `SurfaceMeshShortestPathTraits`.
\tparam VIM a model of `ReadablePropertyMap` with `vertex_descriptor` as key and `unsigned int` as value type.
            The default is `boost::property_map<HG, boost::%vertex_index_t>::%const_type`.
\tparam HIM a model of `ReadablePropertyMap` with `halfedge_descriptor` as key and `unsigned int` as value type.
            The default is `boost::property_map<HG, boost::%halfedge_index_t>::%const_type`.
\tparam FIM a model of `ReadablePropertyMap` with `face_descriptor` as key and `unsigned int` as value type.
            The default is `boost::property_map<HG, boost::%face_index_t>::%const_type`.
\tparam VPM a model of `ReadablePropertyMap` with `vertex_descriptor` as key and `Traits::Point_3` as value type.
            The default is `boost::property_map<HG, CGAL::vertex_point_t>::%const_type`.

If index property maps are not provided through the constructor of the class, internal property maps must
be available and initialized.

\sa \link BGLGraphExternalIndices `CGAL::set_halfedgeds_items_id()`\endlink
*/

template<class Traits,
  class VIM = Default,
  class HIM = Default,
  class FIM = Default,
  class VPM = Default>
class Surface_mesh_shortest_path
{
public:
/// \name Types
/// @{

  /// The triangle mesh type which this algorithm acts on.
  typedef typename Traits::Triangle_mesh Triangle_mesh;

  /// The BGL graph traits for this triangle mesh
  typedef boost::graph_traits<Triangle_mesh> Graph_traits;

  /// Descriptors for the vertices of `Triangle_mesh`
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;

  /// Descriptors for the halfedges of `Triangle_mesh`
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  /// Descriptors of the faces of `Triangle_mesh`
  typedef typename Graph_traits::face_descriptor face_descriptor;

#ifndef DOXYGEN_RUNNING

  typedef typename Default::Get<
    VIM,
    typename boost::property_map<Triangle_mesh, boost::vertex_index_t>::const_type
      >::type Vertex_index_map;

  typedef typename Default::Get<
    HIM,
    typename boost::property_map<Triangle_mesh, boost::halfedge_index_t>::const_type
      >::type Halfedge_index_map;

  typedef typename Default::Get<
    FIM,
    typename boost::property_map<Triangle_mesh, boost::face_index_t>::const_type
      >::type Face_index_map;

  typedef typename Default::Get<
    VPM,
    typename boost::property_map<Triangle_mesh, CGAL::vertex_point_t>::const_type
      >::type Vertex_point_map;

#else
  /// The vertex index property map class
  typedef VIM Vertex_index_map;

  /// The halfedge index property map class
  typedef HIM Halfedge_index_map;

  /// The face index property map class
  typedef FIM Face_index_map;

  /// The vertex point property map class
  typedef VPM Vertex_point_map;
#endif

  /// The numeric type used by this algorithm.
  typedef typename Traits::FT FT;

  /// The 3-dimensional point type, which must coincide with the value type of `Vertex_point_map`.
  typedef typename Traits::Point_3 Point_3;

  /// An ordered triple which specifies a location inside a triangle as
  /// a convex combination of its three vertices.
  typedef typename Traits::Barycentric_coordinates Barycentric_coordinates;

#ifndef CGAL_NO_DEPRECATED_CODE
  // deprecated in CGAL 4.10
  /// \deprecated
  typedef Barycentric_coordinates Barycentric_coordinate;
#endif

  /// \brief An ordered pair specifying a location on the surface of the `Triangle_mesh`.
  /// \details If `tm` is the input graph and given the pair (`f`, `bc`) such that `bc` is `(w0, w1, w2)`,
  ///  the correspondence with the weights in `bc` and the vertices of the face `f` is the following:
  /// - `w0 = source(halfedge(f,tm),tm)`
  /// - `w1 = target(halfedge(f,tm),tm)`
  /// - `w2 = target(next(halfedge(f,tm),tm),tm)`
  typedef std::pair<face_descriptor, Barycentric_coordinates> Face_location;

private:

  typedef std::list<Face_location> Source_point_list;
  typedef typename Source_point_list::iterator Source_point_underlying_iterator;

public:

  /*!
  \brief A model of `BidirectionalIterator` to access the source points

  \details An iterator becomes invalid if:
   - the corresponding point is removed (either with `Surface_mesh_shortest_path::remove_source_point()`
     or `Surface_mesh_shortest_path::remove_all_source_points()`).
   - the structure is re-built (triggered by a shortest path query or a call to `Surface_mesh_shortest_path::build_sequence_tree()`).
   - the structure is cleared (`Surface_mesh_shortest_path::clear()`).

  Dereferencing this iterator yields a `const Surface_mesh_shortest_path::Face_location&`.

  This iterator supports equality comparison operations.
  */
  class Source_point_iterator
  {
  public:
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef Face_location value_type;
    typedef Face_location* pointer;
    typedef Face_location& reference;

  private:
    Source_point_underlying_iterator m_iterator;

    friend class Surface_mesh_shortest_path;

    Source_point_iterator(Source_point_underlying_iterator it)
      : m_iterator(it)
    {
    }

  public:
    /*!
    \brief %Default constructor
    */
    Source_point_iterator()
    {
    }

    /*!
    \brief Copy constructor
    */
    Source_point_iterator(const Source_point_iterator& other)
      : m_iterator(other.m_iterator)
    {
    }

    /*
    \brief Copy the contents of another `Source_point_iterator`

    \param other The iterator to be copied
    */
    Source_point_iterator& operator=(const Source_point_iterator& other)
    {
      m_iterator = other.m_iterator;
      return *this;
    }

    const Face_location& operator*() const
    {
      return *m_iterator;
    }

    const Face_location* operator->() const
    {
      return &(*m_iterator);
    }

    Source_point_iterator& operator++()
    {
      ++m_iterator;
      return *this;
    }

    Source_point_iterator operator++(int)
    {
      Source_point_iterator temp(*this);
      ++m_iterator;
      return temp;
    }

    Source_point_iterator& operator--()
    {
      --m_iterator;
      return *this;
    }

    Source_point_iterator operator--(int)
    {
      Source_point_iterator temp(*this);
      --m_iterator;
      return temp;
    }

    bool operator==(const Source_point_iterator& other) const
    {
      return m_iterator == other.m_iterator;
    }

    bool operator!=(const Source_point_iterator& other) const
    {
      return m_iterator != other.m_iterator;
    }
  };

  /// The return type from shortest path distance queries. Stores the distance
  /// to the nearest source point, and a `Source_point_iterator` to the
  /// source point itself.
  typedef std::pair<FT, Source_point_iterator> Shortest_path_result;

/// @}

private:
  typedef typename Graph_traits::face_iterator face_iterator;

  typedef typename Traits::Triangle_3 Triangle_3;
  typedef typename Traits::Triangle_2 Triangle_2;
  typedef typename Traits::Segment_2 Segment_2;
  typedef typename Traits::Ray_3 Ray_3;
  typedef typename Traits::Ray_2 Ray_2;
  typedef typename Traits::Line_2 Line_2;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Traits::Vector_2 Vector_2;

  typedef Surface_mesh_shortest_paths_3::internal::Cone_tree_node<Traits> Cone_tree_node;
  typedef Surface_mesh_shortest_paths_3::internal::Cone_expansion_event<Traits> Cone_expansion_event;

  typedef std::priority_queue<Cone_expansion_event*,
                              std::vector<Cone_expansion_event*>,
                              Surface_mesh_shortest_paths_3::internal::Cone_expansion_event_min_priority_queue_comparator<Traits> > Expansion_priqueue;
  typedef std::pair<Cone_tree_node*, FT> Node_distance_pair;

private:

  template <class OutputIterator>
  struct Point_path_visitor_wrapper
  {
    Surface_mesh_shortest_path& m_owner;
    OutputIterator m_output;

    Point_path_visitor_wrapper(Surface_mesh_shortest_path& owner, OutputIterator output)
      : m_owner(owner)
      , m_output(output)
    {
    }

    void operator()(halfedge_descriptor edge, FT t)
    {
      *m_output = m_owner.point(edge, t);
      ++m_output;
    }

    void operator()(vertex_descriptor vertex)
    {
      *m_output = m_owner.point(vertex);
      ++m_output;
    }

    void operator()(face_descriptor f, Barycentric_coordinates location)
    {
      *m_output = m_owner.point(f, location);
      ++m_output;
    }
  };

private:
  const Traits m_traits;
  const Triangle_mesh& m_graph;

  Vertex_index_map m_vertexIndexMap;
  Halfedge_index_map m_halfedgeIndexMap;
  Face_index_map m_faceIndexMap;
  Vertex_point_map m_vertexPointMap;

  std::vector<bool> m_vertexIsPseudoSource;

  std::vector<Node_distance_pair> m_vertexOccupiers;
  std::vector<Node_distance_pair> m_closestToVertices;

  std::vector<std::pair<Cone_tree_node*, Source_point_iterator> > m_rootNodes;
  Source_point_list m_faceLocations;
  Source_point_underlying_iterator m_firstNewSourcePoint;
  Source_point_list m_deletedSourceLocations;

  std::vector<std::vector<Cone_tree_node*> > m_faceOccupiers;

  Expansion_priqueue m_expansionPriqueue;

#if !defined(NDEBUG)
  std::size_t m_currentNodeCount;
  std::size_t m_peakNodeCount;
  std::size_t m_queueAtPeakNodes;
  std::size_t m_peakQueueSize;
  std::size_t m_nodesAtPeakQueue;
#endif

#if !defined(NDEBUG)
public:

  /// \cond

  std::size_t peak_node_count() const
  {
    return m_peakNodeCount;
  }

  std::size_t current_node_count() const
  {
    return m_currentNodeCount;
  }

  std::size_t peak_queue_size() const
  {
    return m_peakQueueSize;
  }

  std::size_t current_memory_usage() const
  {
    std::size_t baseUsage = m_rootNodes.size() * sizeof(Cone_tree_node*)
                            + m_closestToVertices.size() * sizeof(Node_distance_pair);

    std::size_t finalUsage = baseUsage + sizeof(Cone_tree_node) * m_currentNodeCount;

    for (std::size_t i = 0; i < m_faceOccupiers.size(); ++i)
    {
      finalUsage += (m_faceOccupiers[i].size() * sizeof(Cone_tree_node*))
                    + sizeof(std::vector<Cone_tree_node*>);
    }

    return finalUsage;
  }

  std::size_t peak_memory_usage() const
  {
    std::size_t baseUsage = m_rootNodes.size() * sizeof(Cone_tree_node*)
                            + m_vertexOccupiers.size() * sizeof(Node_distance_pair)
                            + m_closestToVertices.size() * sizeof(Node_distance_pair);

    std::size_t peakNodeUsage = baseUsage + (sizeof(Cone_tree_node) * m_peakNodeCount)
                                + ((sizeof(Cone_expansion_event)
                                    + sizeof(Cone_expansion_event*)) * m_queueAtPeakNodes);

    std::size_t peakQueueUsage = baseUsage
                                 + (sizeof(Cone_expansion_event) + (sizeof(Cone_expansion_event*)) * m_peakQueueSize)
                                 + (sizeof(Cone_tree_node) * m_nodesAtPeakQueue);

    return (std::max)(peakNodeUsage, peakQueueUsage);
  }

  /// \endcond

#endif

public:

  /// \cond

  /// This is just a placeholder for a proper debug output verbosity switch method
  bool m_debugOutput;

  /// \endcond

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

  Point_2 construct_barycenter_in_triangle_2(const Triangle_2& t,
                                             const Barycentric_coordinates& b) const
  {
    return construct_barycenter_in_triangle_2(t, b, m_traits);
  }

  static Point_2 construct_barycenter_in_triangle_2(const Triangle_2& t,
                                                    const Barycentric_coordinates& b,
                                                    const Traits& traits)
  {
    typename Traits::Construct_vertex_2 cv2(traits.construct_vertex_2_object());
    typename Traits::Construct_barycentric_coordinates_weight cbcw(traits.construct_barycentric_coordinates_weight_object());
    typename Traits::Construct_barycenter_2 cb2(traits.construct_barycenter_2_object());

    return cb2(cv2(t, 0), cbcw(b, 0), cv2(t, 1), cbcw(b, 1), cv2(t, 2), cbcw(b, 2));
  }

  Point_3 construct_barycenter_in_triangle_3(const Triangle_3& t,
                                             const Barycentric_coordinates& b) const
  {
    return construct_barycenter_in_triangle_3(t, b, m_traits);
  }

  static Point_3 construct_barycenter_in_triangle_3(const Triangle_3& t,
                                                    const Barycentric_coordinates& b,
                                                    const Traits& traits)
  {
    typename Traits::Construct_vertex_3 cv3(traits.construct_vertex_3_object());
    typename Traits::Construct_barycentric_coordinates_weight cbcw(traits.construct_barycentric_coordinates_weight_object());
    typename Traits::Construct_barycenter_3 cb3(traits.construct_barycenter_3_object());

    return cb3(cv3(t, 0), cbcw(b, 0), cv3(t, 1), cbcw(b, 1), cv3(t, 2), cbcw(b, 2));
  }

  Triangle_3 triangle_from_halfedge(const halfedge_descriptor edge) const
  {
    return triangle_from_halfedge(edge, m_graph, m_vertexPointMap);
  }

  static Triangle_3 triangle_from_halfedge(const halfedge_descriptor edge,
                                           const Triangle_mesh& tm)
  {
    return triangle_from_halfedge(edge, tm, get(vertex_point, tm));
  }

  static Triangle_3 triangle_from_halfedge(const halfedge_descriptor edge,
                                           const Triangle_mesh& tm,
                                           const Vertex_point_map vertexPointMap)
  {
    return Surface_mesh_shortest_paths_3::internal::triangle_from_halfedge<Triangle_3, Triangle_mesh, Vertex_point_map>(edge, tm, vertexPointMap);
  }

  Triangle_3 triangle_from_face(const face_descriptor f) const
  {
    return triangle_from_face(f, m_graph, m_vertexPointMap);
  }

  static Triangle_3 triangle_from_face(const face_descriptor f,
                                       const Triangle_mesh& tm)
  {
    return triangle_from_halfedge(halfedge(f, tm), tm, get(vertex_point, tm));
  }

  static Triangle_3 triangle_from_face(const face_descriptor f,
                                       const Triangle_mesh& tm,
                                       const Vertex_point_map vertexPointMap)
  {
    return triangle_from_halfedge(halfedge(f, tm), tm, vertexPointMap);
  }

  /*
    Filtering algorithm described in Xin and Wang (2009)
    "Improving chen and han's algorithm on the discrete geodesic problem."
    https://dl.acm.org/citation.cfm?doid=1559755.1559761
  */
  bool window_distance_filter(Cone_tree_node* cone,
                              const Segment_2& windowSegment,
                              const bool reversed)
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());

    const Segment_2& parentEntrySegment = cone->entry_segment();
    const Point_2& v2 = cone->target_point();
    const Point_2& I = cone->source_image();
    const FT d = cone->distance_from_source_to_root();

    FT d1;
    FT d2;
    FT d3;
    Point_2 A;
    Point_2 B;
    Point_2 v1;
    Point_2 v3;

    std::size_t v1Index = get(m_vertexIndexMap, source(cone->entry_edge(), m_graph));
    std::size_t v2Index = get(m_vertexIndexMap, cone->target_vertex());
    std::size_t v3Index = get(m_vertexIndexMap, target(cone->entry_edge(), m_graph));

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

    const bool hasD1 = v1Distance.first != nullptr && v1Distance.first != cone->parent();
    const bool hasD2 = v2Distance.first != nullptr && v2Distance.first != cone->parent();
    const bool hasD3 = v3Distance.first != nullptr && v3Distance.first != cone->parent();

    if (hasD1 && (d + CGAL::approximate_sqrt(csd2(I, B)) > d1 + CGAL::approximate_sqrt(csd2(v1, B))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,B| > d1 + |v1,B|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index
                  << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v1,B = " << CGAL::approximate_sqrt(csd2(v1, B)) << std::endl;
        std::cout << "I,B = " << CGAL::approximate_sqrt(csd2(I, B)) << std::endl;
        std::cout << "I,A = " << CGAL::approximate_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::approximate_sqrt(csd2(I, B)))
                  << " vs. " << (d1 + CGAL::approximate_sqrt(csd2(v1, B))) << std::endl;
      }

      return false;
    }

    if (hasD2 && (d + CGAL::approximate_sqrt(csd2(I, A)) > d2 + CGAL::approximate_sqrt(csd2(v2, A))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,A| > d2 + |v2,A|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index
                  << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v2,A = " << CGAL::approximate_sqrt(csd2(v2, A)) << std::endl;
        std::cout << "I,A = " << CGAL::approximate_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::approximate_sqrt(csd2(I, A)))
                  << " vs. " << (d2 + CGAL::approximate_sqrt(csd2(v2, A))) << std::endl;
      }

      return false;
    }

    if (hasD3 && (d + CGAL::approximate_sqrt(csd2(I, A)) > d3 + CGAL::approximate_sqrt(csd2(v3, A))))
    {
      if (m_debugOutput)
      {
        std::cout << "Filter: d + |I,A| > d3 + |v3,A|: " << std::endl;
        std::cout << "v1 = " << v1Index << " , " << d1 << " , v2 = " << v2Index
                  << " , " << d2 << " , v3 = " << v3Index << " , " << d3 << std::endl;
        std::cout << "d = " << d << std::endl;
        std::cout << "v3,A = " << CGAL::approximate_sqrt(csd2(v3, A)) << std::endl;
        std::cout << "I,A = " << CGAL::approximate_sqrt(csd2(I, A)) << std::endl;
        std::cout << (d + CGAL::approximate_sqrt(csd2(I, A)))
                  << " vs. " << (d3 + CGAL::approximate_sqrt(csd2(v3, A))) << std::endl;
      }

      return false;
    }

    return true;
  }

  /*
    Push a new node representing crossing the edge to the left of `cone`'s target vertex
  */
  void expand_left_child(Cone_tree_node* cone,
                         const Segment_2& windowSegment)
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Construct_triangle_3_along_segment_2_flattening ft3as2(m_traits.construct_triangle_3_along_segment_2_flattening_object());

    if (m_debugOutput)
    {
      std::cout << std::endl << " >>>>>>>>>>>>>>>>>>> Expanding LEFT CHILD <<<<<<<<<<<<<<<<<<<" <<std::endl;
    }

    CGAL_assertion(cone->m_pendingLeftSubtree != nullptr);
    cone->m_pendingLeftSubtree = nullptr;

    if (window_distance_filter(cone, windowSegment, false))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->left_child_edge());
      Triangle_2 layoutFace = ft3as2(adjacentFace, 0, cone->left_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_traits, m_graph, cone->left_child_edge(), layoutFace,
                                                 cone->source_image(), cone->distance_from_source_to_root(),
                                                 cv2(windowSegment, 0), cv2(windowSegment, 1),
                                                 Cone_tree_node::INTERVAL);
      node_created();
      cone->set_left_child(child);
      process_node(child);
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNode was filtered." << std::endl;
    }
  }

  /*
    Push a new node representing crossing the edge to the right of `cone`'s target vertex
  */
  void expand_right_child(Cone_tree_node* cone,
                          const Segment_2& windowSegment)
  {
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Construct_triangle_3_along_segment_2_flattening ft3as2(m_traits.construct_triangle_3_along_segment_2_flattening_object());

    if (m_debugOutput)
    {
      std::cout << std::endl << " >>>>>>>>>>>>>>>>>>> Expanding RIGHT CHILD <<<<<<<<<<<<<<<<<<<" <<std::endl;
    }

    CGAL_assertion(cone->m_pendingRightSubtree != nullptr);
    cone->m_pendingRightSubtree = nullptr;

    if (window_distance_filter(cone, windowSegment, true))
    {
      Triangle_3 adjacentFace = triangle_from_halfedge(cone->right_child_edge());
      Triangle_2 layoutFace = ft3as2(adjacentFace, 0, cone->right_child_base_segment());
      Cone_tree_node* child = new Cone_tree_node(m_traits, m_graph, cone->right_child_edge(), layoutFace,
                                                 cone->source_image(), cone->distance_from_source_to_root(),
                                                 cv2(windowSegment, 0), cv2(windowSegment, 1),
                                                 Cone_tree_node::INTERVAL);
      node_created();
      cone->set_right_child(child);
      process_node(child);
    }
    else if (m_debugOutput)
    {
      std::cout << "\tNode was filtered." << std::endl;
    }
  }

  /*
    Determines whether to expand `location` as a face, edge, or vertex root, depending on
    whether it is near to a given edge or vertex, or is an internal face location
  */
  void expand_root(const face_descriptor f,
                   const Barycentric_coordinates& location,
                   Source_point_iterator sourcePointIt)
  {
    typename Traits::Construct_barycentric_coordinates_weight cbcw(m_traits.construct_barycentric_coordinates_weight_object());
    typename Traits::Classify_barycentric_coordinates classify_barycentric_coordinates(m_traits.classify_barycentric_coordinates_object());

    std::size_t associatedEdge;
    CGAL::Surface_mesh_shortest_paths_3::Barycentric_coordinates_type type;
    std::tie(type, associatedEdge) = classify_barycentric_coordinates(location);

    switch (type)
    {
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_BOUNDED_SIDE:
        expand_face_root(f, location, sourcePointIt);
        break;
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_BOUNDARY:
        {
          halfedge_descriptor he = halfedge(f, m_graph);
          for (std::size_t i = 0; i < associatedEdge; ++i)
          {
            he = next(he, m_graph);
          }
          expand_edge_root(he, cbcw(location, associatedEdge), cbcw(location, (associatedEdge + 1) % 3), sourcePointIt);
        }
        break;
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_VERTEX:
        {
          halfedge_descriptor he = halfedge(f, m_graph);
          for (std::size_t i = 0; i < associatedEdge; ++i)
          {
            he = next(he, m_graph);
          }
          expand_vertex_root(source(he, m_graph), sourcePointIt);
        }
        break;
      default:
        CGAL_assertion(false && "Invalid face location");
        // Perhaps hit an assertion that the type must not be external or invalid?
    }
  }

  /*
    Create source nodes facing each edge of `f`, rooted at the given `faceLocation`
  */
  void expand_face_root(const face_descriptor f,
                        const Barycentric_coordinates& faceLocation,
                        Source_point_iterator sourcePointIt)
  {
    typename Traits::Construct_triangle_3_to_triangle_2_projection pt3t2(m_traits.construct_triangle_3_to_triangle_2_projection_object());
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());

    const halfedge_descriptor start = halfedge(f, m_graph);
    halfedge_descriptor current = start;

    Cone_tree_node* faceRoot = new Cone_tree_node(m_traits, m_graph, m_rootNodes.size());
    node_created();
    m_rootNodes.emplace_back(faceRoot, sourcePointIt);

    if (m_debugOutput)
    {
      typename Traits::Construct_barycentric_coordinates_weight cbcw(m_traits.construct_barycentric_coordinates_weight_object());

      std::cout << std::endl << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "\tFace Root Expansion: id = " << get(m_faceIndexMap, f)
                << " , Location = " << cbcw(faceLocation, 0) << " " << cbcw(faceLocation, 1)
                << " " << cbcw(faceLocation, 2) << " " << std::endl;
    }

    for (std::size_t currentVertex = 0; currentVertex < 3; ++currentVertex)
    {
      const Triangle_3 face3d(triangle_from_halfedge(current));
      const Triangle_2 layoutFace(pt3t2(face3d));
      const Barycentric_coordinates rotatedFaceLocation(shifted_coordinates(faceLocation, currentVertex));
      const Point_2 sourcePoint(construct_barycenter_in_triangle_2(layoutFace, rotatedFaceLocation));

      Cone_tree_node* child = new Cone_tree_node(m_traits, m_graph,
                                                 current /*entryEdge*/,
                                                 layoutFace, sourcePoint,
                                                 FT(0) /*pseudoSourceDistance*/,
                                                 cv2(layoutFace, 0) /*windowLeft*/,
                                                 cv2(layoutFace, 2) /*windowRight*/,
                                                 Cone_tree_node::FACE_SOURCE);
      node_created();
      faceRoot->push_middle_child(child);

      if (m_debugOutput)
      {
        std::cout << "\tExpanding face root #" << currentVertex << " : " << std::endl;;
        std::cout << "\t\t3D Face = " << face3d << std::endl;
        std::cout << "\t\t2D Face = " << layoutFace << std::endl;
        std::cout << "\t\tLocation = " << sourcePoint << std::endl;
      }

      process_node(child);

      current = next(current, m_graph);
    }
  }

  /*
    Create 'source' nodes to each size of the given edge, rooted at the specified parametric location
  */
  void expand_edge_root(const halfedge_descriptor baseEdge,
                        const FT t0, const FT t1,
                        Source_point_iterator sourcePointIt)
  {
    CGAL_precondition(!is_border(baseEdge, m_graph));
    CGAL_precondition(t0 + t1 == FT(1));

    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());
    typename Traits::Construct_triangle_3_to_triangle_2_projection pt3t2(m_traits.construct_triangle_3_to_triangle_2_projection_object());

    if (m_debugOutput)
    {
      std::cout << std::endl << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "\tEdge Root Expansion: faceA = " << get(m_faceIndexMap, face(baseEdge, m_graph))
                << " , faceB = " << get(m_faceIndexMap, face(opposite(baseEdge, m_graph), m_graph))
                << " , t0 = " << t0 << " , t1 = " << t1 << std::endl;
      std::cout << "\t\tBoundary: " << is_border_edge(baseEdge, m_graph) << std::endl;
    }

    Cone_tree_node* edgeRoot = new Cone_tree_node(m_traits, m_graph, m_rootNodes.size());
    node_created();
    m_rootNodes.emplace_back(edgeRoot, sourcePointIt);

    /* If v0v1 is not a border edge:
     *
     *      v2
     *     /  \
     *    /    \
     *   /      \
     * v0 - S - v1
     *   \      /
     *    \    /
     *     \  /
     *      v3
     * The source S must reach all Vi, so for each side of the edge, there are two windwows being spawned:
     * - v0v1 targeting v2 propagating only on the left (v0v2)
     * - v2v0 targeting v1 propagating only on the left (v2v1)
     * - v1v0 targeting v3 propagating only on the left (v1v3)
     * - v3v1 targeting v0 propagating only on the left (v3v0)
     *
     * If v0v1 is a border edge, spawn 3 children in the face, and none on the other side
     */

    if(is_border_edge(baseEdge, m_graph))
    {
      const Face_location edgeSourceLocation = face_location(baseEdge, t0);
      return expand_face_root(face(baseEdge, m_graph), edgeSourceLocation.second, sourcePointIt);
    }

    // From here on, it is not a border edge --> spawn 2 children on each side

    halfedge_descriptor baseEdges[2];
    baseEdges[0] = baseEdge;
    baseEdges[1] = opposite(baseEdge, m_graph);

    // shift is because the entry halfedge is not necessarily equal to halfedge(face(entry_h, g), g)
    Barycentric_coordinates edgeSourceLocations[2];
    edgeSourceLocations[0] = shifted_coordinates(face_location(baseEdges[0], t0).second,
                                                 Surface_mesh_shortest_paths_3::internal::edge_index(baseEdges[0], m_graph));
    edgeSourceLocations[1] = shifted_coordinates(face_location(baseEdges[1], t1).second,
                                                 Surface_mesh_shortest_paths_3::internal::edge_index(baseEdges[1], m_graph));

    for (std::size_t side = 0; side < 2; ++side)
    {
      Triangle_3 face3d(triangle_from_halfedge(baseEdges[side]));
      Triangle_2 layoutFace(pt3t2(face3d));
      Point_2 sourcePoint(construct_barycenter_in_triangle_2(layoutFace, edgeSourceLocations[side]));

      // v0v1 targeting v2
      if (m_debugOutput)
      {
        std::cout << std::endl << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "\tExpanding edge root, side #" << side << ", targeting LOCAL 'v2'" << std::endl;
        std::cout << "\t\t3D Face = " << face3d << std::endl;
        std::cout << "\t\t2D Face = " << layoutFace << std::endl;
        std::cout << "\t\tBarycentric coordinates: " << edgeSourceLocations[side][0]
                                              << " " << edgeSourceLocations[side][1]
                                              << " " << edgeSourceLocations[side][2] << std::endl;
        std::cout << "\t\tLocation = " << sourcePoint << std::endl;
      }

      Cone_tree_node* v2_Child = new Cone_tree_node(m_traits, m_graph,
                                                    baseEdges[side] /*entryEdge*/,
                                                    layoutFace,
                                                    sourcePoint /*sourceImage*/,
                                                    FT(0) /*pseudoSourceDistance*/,
                                                    cv2(layoutFace, 0) /*windowLeft*/,
                                                    cv2(layoutFace, 2) /*windowRight*/,
                                                    Cone_tree_node::EDGE_SOURCE);
      node_created();
      edgeRoot->push_middle_child(v2_Child);
      process_node(v2_Child);

      // v2v0 targeting v1
      face3d = triangle_from_halfedge(prev(baseEdges[side], m_graph));
      layoutFace = pt3t2(face3d);

      // shift the barycentric coordinates to correspond to the new layout
      std::swap(edgeSourceLocations[side][1], edgeSourceLocations[side][2]);
      std::swap(edgeSourceLocations[side][0], edgeSourceLocations[side][1]);
      sourcePoint = Point_2(construct_barycenter_in_triangle_2(layoutFace, edgeSourceLocations[side]));

      if (m_debugOutput)
      {
        std::cout << std::endl << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "\tExpanding edge root, side #" << side << ", targeting LOCAL 'v1'" << std::endl;
        std::cout << "\t\t3D Face = " << face3d << std::endl;
        std::cout << "\t\t2D Face = " << layoutFace << std::endl;
        std::cout << "\t\tBarycentric coordinates: " << edgeSourceLocations[side][0]
                                              << " " << edgeSourceLocations[side][1]
                                              << " " << edgeSourceLocations[side][2] << std::endl;
        std::cout << "\t\tLocation = " << sourcePoint << std::endl;
      }

      Cone_tree_node* v1_Child = new Cone_tree_node(m_traits, m_graph,
                                                    prev(baseEdges[side], m_graph) /*entryEdge*/,
                                                    layoutFace,
                                                    sourcePoint /*sourceImage*/,
                                                    FT(0) /*pseudoSourceDistance*/,
                                                    cv2(layoutFace, 0) /*windowLeft*/,
                                                    cv2(layoutFace, 2) /*windowRight*/,
                                                    Cone_tree_node::EDGE_SOURCE);
      node_created();
      edgeRoot->push_middle_child(v1_Child);
      process_node(v1_Child);
    }
  }

  /*
    Create a 'source' node for each face surrounding the given vertex.
  */
  void expand_vertex_root(const vertex_descriptor vertex,
                          Source_point_iterator sourcePointIt)
  {
    if (m_debugOutput)
    {
      std::cout << std::endl << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "\tVertex Root Expansion: Vertex = " << get(m_vertexIndexMap, vertex) << std::endl;
    }

    Cone_tree_node* vertexRoot = new Cone_tree_node(m_traits, m_graph, m_rootNodes.size(),
                                                    prev(halfedge(vertex, m_graph), m_graph));

    node_created();
    m_rootNodes.emplace_back(vertexRoot, sourcePointIt);

    m_closestToVertices[get(m_vertexIndexMap, vertex)] = Node_distance_pair(vertexRoot, FT(0));

    expand_pseudo_source(vertexRoot);
  }

  /*
    Create child nodes for each face surrounding the vertex occupied by `parent`, and push them to the queue

    By convention, source windows are left & middle (no right) as to not create overlaps.
    A child is also created for the border halfedge (if any) as to propagate distance
    to the next vertex on the border (aka `target(next(border_h, g), g)`).
    This creates a nonsensical triangle made of `source(border_h, g)`, `target(border_h, g)`,
    and `target(next(border_h, g), g)` but propagation is only done to the vertex (vertex source:
    no propagation on the right by convention, and left is a border halfedge so no propagation either).
  */
  void expand_pseudo_source(Cone_tree_node* parent)
  {
    typename Traits::Construct_triangle_3_to_triangle_2_projection pt3t2(m_traits.construct_triangle_3_to_triangle_2_projection_object());
    typename Traits::Construct_vertex_2 cv2(m_traits.construct_vertex_2_object());

    parent->m_pendingMiddleSubtree = nullptr;

    vertex_descriptor expansionVertex = parent->target_vertex();

    halfedge_descriptor startEdge = halfedge(expansionVertex, m_graph);
    halfedge_descriptor currentEdge = halfedge(expansionVertex, m_graph);

    FT distanceFromTargetToRoot = parent->distance_from_target_to_root();

    if (m_debugOutput)
    {
      std::cout << "Pseudo source: V" << get(m_vertexIndexMap, expansionVertex) << std::endl;
      std::cout << "Distance from pseudo source to root: " << distanceFromTargetToRoot << std::endl;
    }

    // A potential optimization could be made by only expanding in the 'necessary' range (i.e. the range outside of geodesic visibility), but the
    // benefits may be small, since the node filter would prevent more than one-level propagation.
    // It would also be necessary to distinguish expanding a root vertex node from a pseudo-source node

    do
    {
      const Triangle_3 face3d(triangle_from_halfedge(currentEdge));
      const Triangle_2 layoutFace(pt3t2(face3d));

      if (m_debugOutput)
      {
        std::cout << std::endl << " >>>>>>>>>>>>>>>>>>> Expanding PseudoSource <<<<<<<<<<<<<<<<<<<" <<std::endl;
        std::cout << "currentEdge: "
                  << get(m_vertexIndexMap, source(currentEdge, m_graph)) << " "
                  << get(m_vertexIndexMap, target(currentEdge, m_graph)) << std::endl;
        std::cout << "face id = ";
        if (!is_border(currentEdge, m_graph))
        {
          std::cout << get(m_faceIndexMap, face(currentEdge, m_graph)) << std::endl;

          std::cout << "3D face:" << std::endl << face3d[0] << std::endl << face3d[1] << std::endl << face3d[2] << std::endl;
          std::cout << "current face: " << get(m_faceIndexMap, face(currentEdge, m_graph)) << " gives 2D layout: " << layoutFace << std::endl;
          std::cout << "source: " << face3d[1] << std::endl;
        }
        else
        {
          std::cout << "EXTERNAL" << std::endl;
        }
      }

      Cone_tree_node* child = new Cone_tree_node(m_traits, m_graph, currentEdge /*entryEdge*/,
                                                 layoutFace, cv2(layoutFace, 1) /*sourceImage*/,
                                                 distanceFromTargetToRoot,
                                                 cv2(layoutFace, 0) /*windowLeft*/,
                                                 cv2(layoutFace, 2) /*windowRight*/,
                                                 Cone_tree_node::VERTEX_SOURCE);

      node_created();
      parent->push_middle_child(child);
      process_node(child);

      currentEdge = opposite(next(currentEdge, m_graph), m_graph);
    }
    while (currentEdge != startEdge);

  }

  /*
    Returns the intersection of `segment` and the cone defined by the region to the left `leftBoundary` and right of `rightBoundary`
  */
  bool clip_to_bounds(const Segment_2& segment,
                      const Ray_2& leftBoundary,
                      const Ray_2& rightBoundary,
                      Segment_2& outSegment) const
  {
    typename Traits::Construct_source_2 cs2(m_traits.construct_source_2_object());
    typename Traits::Construct_segment_2 cseg2(m_traits.construct_segment_2_object());
    typename Traits::Construct_target_2 ct2(m_traits.construct_target_2_object());
    typename Traits::Construct_line_2 cl2(m_traits.construct_line_2_object());
    typename Traits::Intersect_2 i2(m_traits.intersect_2_object());
    typename Traits::Orientation_2 o2(m_traits.orientation_2_object());
    typename Traits::Construct_point_on_2 cpo2(m_traits.construct_point_on_2_object());
    typename Traits::Compute_parametric_distance_along_segment_2 pdas2(m_traits.compute_parametric_distance_along_segment_2_object());

    Point_2 leftPoint;
    Point_2 rightPoint;

    if (m_debugOutput)
    {
      std::cout << "Clipping Segment " << segment << std::endl
                << "\t with left = " << leftBoundary << std::endl
                << "\t with right = " << rightBoundary << std::endl;
    }

    FT leftT;
    FT rightT;

    if (o2(cs2(leftBoundary), cpo2(leftBoundary, 1), cs2(segment)) != CGAL::LEFT_TURN)
    {
      if (m_debugOutput)
      {
        std::cout << "\tLeft is completely covered." << std::endl;
      }
      leftPoint = cs2(segment);
      leftT = FT(0);
    }
    else
    {
      const auto cgalIntersection = i2(cl2(segment), cl2(leftBoundary));

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
        const Point_2* result = boost::get<Point_2>(&*cgalIntersection);
        FT t0 = pdas2(cs2(segment), ct2(segment), *result);

        if (t0 >= FT(1))
        {
          if (m_debugOutput)
          {
            std::cout << "Dropping due to missing left intersect. " << t0 << std::endl;
          }

          return false;
        }
        else if (t0 <= FT(0))
        {
          if (m_debugOutput)
          {
            std::cout << "\tLeft is completely covered (secondary check). " << t0 << std::endl;
          }

          leftPoint = cs2(segment);
          leftT = FT(0);
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

    if (o2(cs2(rightBoundary), cpo2(rightBoundary, 1), ct2(segment)) != CGAL::RIGHT_TURN)
    {
      if (m_debugOutput)
      {
        std::cout << "\tRight is completely covered." << std::endl;
      }
      rightPoint = ct2(segment);
      rightT = FT(1);
    }
    else
    {
      const auto cgalIntersection = i2(cl2(segment), cl2(rightBoundary));

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
        const Point_2* result = boost::get<Point_2>(&*cgalIntersection);
        FT t0 = pdas2(cs2(segment), ct2(segment), *result);

        if (t0 <= FT(0))
        {
          if (m_debugOutput)
          {
            std::cout << "Dropping due to missing right intersect. " << t0 << std::endl;
          }
          return false;
        }
        else if (t0 >= FT(1))
        {
          if (m_debugOutput)
          {
            std::cout << "\tRight is completely covered (secondary check). " << t0 << std::endl;
          }
          rightPoint = ct2(segment);
          rightT = FT(1);
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

  /*
    Take a node and compute whether it is an occupier/evicts older nodes, then push any children it may have.
  */
  void process_node(Cone_tree_node* node)
  {
    if (m_debugOutput)
    {
      std::cout << std::endl << " ---------------- Processing node ---------------" << std::endl;
      std::cout << "Node: " << node << std::endl;
      std::cout << "Node type: " << node->node_type() << std::endl;
      std::cout << "Tree ID: " << node->tree_id() << " at level = " << node->level() << std::endl;
      std::cout << "\tParent node: " << node->parent() << std::endl;
      std::cout << "\tParent node type: " << node->parent()->node_type() << std::endl;
      std::cout << "\tFace = " << node->layout_face() << std::endl;
      std::cout << "\tVertices =";
      halfedge_descriptor current = node->entry_edge();
      for (std::size_t i = 0; i<3; ++i)
      {
        std::cout << " " << get(m_vertexIndexMap, source(current, m_graph));
        current = next(current, m_graph);
      }
      std::cout << std::endl;
      std::cout << "\tSource Image = " << node->source_image() << std::endl;
      std::cout << "\tEntry Halfedge = (V" << get(m_vertexIndexMap, source(node->entry_edge(), m_graph)) << " V"
                                           << get(m_vertexIndexMap, target(node->entry_edge(), m_graph)) << ")" << std::endl;
      std::cout << "\tTarget vertex = V" << get(m_vertexIndexMap, node->target_vertex()) << std::endl;

      std::cout << "\tWindow Left = " << node->window_left() << std::endl;
      std::cout << "\tWindow Right = " << node->window_right() << std::endl;
    }

    bool leftSide = false;
    bool rightSide = false;

    if (!node->is_source_node())
    {
      leftSide = node->has_left_side();
      rightSide = node->has_right_side();
    }
    else // source nodes only have left sides
    {
      leftSide = true;
      rightSide = false;
    }

    if (m_debugOutput)
    {
      std::cout << "\t Has Left : " << (leftSide ? "yes" : "no")
                << " , Has Right : " << (rightSide ? "yes" : "no") << std::endl;
    }

    bool propagateLeft = false;
    bool propagateRight = false;
    bool propagateMiddle = false;

    if (node->is_source_node() || (leftSide && rightSide))
    {
      if (m_debugOutput)
      {
        std::cout << "\tContains target vertex" << std::endl;
      }

      std::size_t entryHalfEdgeIndex = get(m_halfedgeIndexMap, node->entry_edge());

      const Node_distance_pair& currentOccupier = m_vertexOccupiers[entryHalfEdgeIndex];
      FT currentNodeDistance = node->distance_from_target_to_root();

      if (m_debugOutput)
      {
        std::cout << "\t Distance to target: " << currentNodeDistance << std::endl;
        std::cout << "\t Is there a current occupier? " << (currentOccupier.first != nullptr) << std::endl;
      }

      // the relative position of the ray between node.source() and node.target_vertex() and the ray
      // from occupier.source() (-1 left, 0 collinear, 1 right)
      CGAL::Comparison_result c = CGAL::EQUAL; // initializing to please weak compilers

      if (currentOccupier.first != nullptr)
      {
        CGAL_assertion(node->entry_edge() == currentOccupier.first->entry_edge());
        CGAL_assertion(node->target_vertex() == currentOccupier.first->target_vertex());

        // for a vertex source, the ray is along the halfedge pointing towards the target
        if (node->is_vertex_node())
        {
          if (currentOccupier.first->is_vertex_node())
            c = CGAL::EQUAL;
          else
            c = CGAL::LARGER;
        }
        else if (currentOccupier.first->is_vertex_node()) // node is not a vertex source
        {
          c = CGAL::SMALLER;
        }
        else // generic case
        {
          // must compute intersections because although entry edges are identical, their 2D representation is not
          c = m_traits.compare_relative_intersection_along_segment_2_object()(
                node->entry_segment(),
                node->ray_to_target_vertex().supporting_line(),
                currentOccupier.first->entry_segment(),
                currentOccupier.first->ray_to_target_vertex().supporting_line());
        }

        if (m_debugOutput)
        {
          std::cout << "\t Current occupier, EH (V"
                    << get(m_vertexIndexMap, source(currentOccupier.first->entry_edge(), m_graph)) << " V"
                    << get(m_vertexIndexMap, target(currentOccupier.first->entry_edge(), m_graph)) << ")" << std::endl;
          std::cout << "\t Current occupier, Source = " << currentOccupier.first->source_image() << std::endl;
          std::cout << "\t Current occupier Distance = " << currentOccupier.second << std::endl;
          std::cout << "\t smaller (-1)/equal (0)/larger (1) comparison? " << c << std::endl;
        }
      }

      bool is_node_new_occupier = false;
      if (currentOccupier.first == nullptr)
      {
        m_vertexOccupiers[entryHalfEdgeIndex] = std::make_pair(node, currentNodeDistance);
        is_node_new_occupier = true;
      }
      else
      {
        // Only replace the current occupier if the time is _strictly_ larger
        // and yield the way to vertex sources (cleaner than manipulating 0-length intervals)
        if (currentOccupier.second > currentNodeDistance ||
           (currentOccupier.second == currentNodeDistance && node->node_type() == Cone_tree_node::VERTEX_SOURCE))
        {
          m_vertexOccupiers[entryHalfEdgeIndex] = std::make_pair(node, currentNodeDistance);
          is_node_new_occupier = true;
        }
      }

      if (is_node_new_occupier)
      {
        if (m_debugOutput)
        {
          std::cout << "\t Current node is now the occupier of target vertex "
                    << get(m_vertexIndexMap, node->target_vertex()) << std::endl;
        }

        propagateLeft = true;
        propagateRight = true;

        // This is a consequence of using the same basic node type for source and interval nodes
        // If this is a source node, it is only pointing to one of the two opposite edges (the left one by convention)
        if (node->is_source_node())
        {
          propagateRight = false;

          // Propagating a pseudo-source on a boundary vertex can result in a cone on a null face
          // In such a case, we only care about the part of the cone pointing at the vertex (i.e. the middle child),
          // so we can avoid propagating over the (non-existent) left opposite edge
          if (node->is_null_face())
          {
            propagateLeft = false;
          }
        }

        // Some branches from the old occupier that has been superseded can now be pruned
        if (currentOccupier.first != nullptr)
        {
          if (c == CGAL::SMALLER) // node's ray is left of occupier's ray
          {
            if (currentOccupier.first->get_left_child())
            {
              delete_node(currentOccupier.first->remove_left_child());
            }
            else if (currentOccupier.first->m_pendingLeftSubtree != nullptr)
            {
              currentOccupier.first->m_pendingLeftSubtree->m_cancelled = true;
              currentOccupier.first->m_pendingLeftSubtree = nullptr;
            }
          }
          else if (c == CGAL::LARGER) // node's ray is right of occupier's ray
          {
            if (currentOccupier.first->get_right_child())
            {
              delete_node(currentOccupier.first->remove_right_child());
            }
            else if (currentOccupier.first->m_pendingRightSubtree != nullptr)
            {
              currentOccupier.first->m_pendingRightSubtree->m_cancelled = true;
              currentOccupier.first->m_pendingRightSubtree = nullptr;
            }
          }
        }

        // Check if `node` is now the absolute closest node, and replace the current closest as appropriate
        std::size_t targetVertexIndex = get(m_vertexIndexMap, node->target_vertex());
        const Node_distance_pair& currentClosest = m_closestToVertices[targetVertexIndex];

        if (m_debugOutput && currentClosest.first != nullptr)
        {
          std::cout << "\t Current Closest Distance = " << currentClosest.second << std::endl;
        }

        // If equal times, give priority to vertex sources since it's cleaner and simpler to handle than interval windows
        if (currentClosest.first == nullptr ||
            currentClosest.second > currentNodeDistance ||
            (currentClosest.second == currentNodeDistance &&
             node->node_type() == Cone_tree_node::VERTEX_SOURCE))
        {
          if (m_debugOutput)
          {
            std::cout << "\t Current node is now the closest at target vertex V"
                      << get(m_vertexIndexMap, node->target_vertex()) << std::endl;
          }

          // if this is a saddle vertex, then evict previous closest vertex
          if (m_vertexIsPseudoSource[targetVertexIndex])
          {
            if (m_debugOutput)
            {
              std::cout << "\t Vertex V" << targetVertexIndex << " is a pseudo-source" << std::endl;
            }

            if (currentClosest.first != nullptr)
            {
              if (m_debugOutput)
              {
                std::cout << "\tEvicting old pseudo-source: " << currentClosest.first << std::endl;
              }

              if (currentClosest.first->m_pendingMiddleSubtree != nullptr)
              {
                currentClosest.first->m_pendingMiddleSubtree->m_cancelled = true;
                currentClosest.first->m_pendingMiddleSubtree = nullptr;
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
      else // there is already an occupier, at a strictly smaller distance
      {
        // this is an application of "one angle one split"
        if (c != CGAL::LARGER) // propagate on the left if the node's ray is left of the occupier's
          propagateLeft = true;
        if (c != CGAL::SMALLER && !node->is_source_node()) // by convention a source node only points at the left edge
          propagateRight = true;

        // No point propagating middle because we know the current occupier has a better time
        // at the target and if the target is a saddle, middle children have already been spawned
        // when we were in the current occupier.
      }
    }
    else // interval node that does not contain the target vertex
    {
      propagateLeft = leftSide;
      propagateRight = rightSide;
    }

    if (m_debugOutput)
    {
      std::cout << "Propagate (L/M/R): " << propagateLeft << " " << propagateMiddle << " " << propagateRight << std::endl;
    }

    if (node->level() <= static_cast<std::size_t>(num_faces(m_graph)))
    {
      if (propagateLeft)
      {
        CGAL_assertion(!node->is_null_face());
        push_left_child(node);
      }

      if (propagateRight)
      {
        CGAL_assertion(!node->is_source_node());
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
    if (m_debugOutput)
    {
      std::cout << "Tentative push of left child edge "
                << " (V" << get(m_vertexIndexMap, source(parent->left_child_edge(), m_graph))
                << " V" << get(m_vertexIndexMap, target(parent->left_child_edge(), m_graph)) << ")" << std::endl;
      std::cout << "Boundary? " << is_border(parent->left_child_edge(), m_graph) << std::endl;
    }

    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());

    if (!is_border(parent->left_child_edge(), m_graph))
    {
      Segment_2 leftWindow;

      if (parent->is_source_node())
      {
        leftWindow = parent->left_child_base_segment();
      }
      else
      {
        bool result = clip_to_bounds(parent->left_child_base_segment(), parent->left_boundary(),
                                     parent->right_boundary(), leftWindow);
        if (!result)
        {
          if (m_debugOutput)
          {
            std::cout << "Left child clip failed, killing node." << std::endl;
          }
          return;
        }
      }

      FT distanceEstimate = parent->distance_from_source_to_root()
                            + CGAL::approximate_sqrt(csd2(parent->source_image(), leftWindow));

      if (m_debugOutput)
      {
        std::cout << ">>> Pushing Left Child, Segment = " << parent->left_child_base_segment()
                  << " , clipped = " << leftWindow << " , Estimate = " << distanceEstimate << std::endl;
      }

      Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate,
                                                             Cone_expansion_event::LEFT_CHILD, leftWindow);
      parent->m_pendingLeftSubtree = event;

      m_expansionPriqueue.push(event);
      queue_pushed();
    }
  }

  void push_right_child(Cone_tree_node* parent)
  {
    if (m_debugOutput)
    {
      std::cout << "Tentative push of right child edge"
                << " (V" << get(m_vertexIndexMap, source(parent->right_child_edge(), m_graph))
                << " V" << get(m_vertexIndexMap, target(parent->right_child_edge(), m_graph)) << ")" << std::endl;
      std::cout << "Boundary? " << is_border(parent->right_child_edge(), m_graph) << std::endl;
    }

    typename Traits::Compute_squared_distance_2 csd2(m_traits.compute_squared_distance_2_object());

    if (!is_border(parent->right_child_edge(), m_graph))
    {
      Segment_2 rightWindow;
      bool result = clip_to_bounds(parent->right_child_base_segment(), parent->left_boundary(),
                                   parent->right_boundary(), rightWindow);

      if (!result)
      {
        if (m_debugOutput)
        {
          std::cout << "Right child clip failed, killing node." << std::endl;
        }
        return;
      }

      FT distanceEstimate = parent->distance_from_source_to_root()
                            + CGAL::approximate_sqrt(csd2(parent->source_image(), rightWindow));

      if (m_debugOutput)
      {
        std::cout << ">>> Pushing Right Child, Segment = " << parent->right_child_base_segment()
                  << " , clipped = " << rightWindow << " , Estimate = " << distanceEstimate << std::endl;
      }

      Cone_expansion_event* event = new Cone_expansion_event(parent, distanceEstimate,
                                                             Cone_expansion_event::RIGHT_CHILD, rightWindow);
      parent->m_pendingRightSubtree = event;

      m_expansionPriqueue.push(event);
      queue_pushed();
    }
  }

  void push_middle_child(Cone_tree_node* parent)
  {
    if (m_debugOutput)
    {
      std::cout << ">>> Pushing Middle Child, Estimate = " << parent->distance_from_target_to_root() << std::endl;
    }

    Cone_expansion_event* event = new Cone_expansion_event(parent, parent->distance_from_target_to_root(), Cone_expansion_event::PSEUDO_SOURCE);
    parent->m_pendingMiddleSubtree = event;

    m_expansionPriqueue.push(event);
    queue_pushed();
  }

  void delete_node(Cone_tree_node* node,
                   const bool destruction = false)
  {
    if (node != nullptr)
    {
      if (m_debugOutput)
      {
        std::cout << "Deleting node " << node << std::endl;
      }

      if (node->m_pendingLeftSubtree != nullptr)
      {
        node->m_pendingLeftSubtree->m_cancelled = true;
        node->m_pendingLeftSubtree = nullptr;
      }

      if (node->get_left_child() != nullptr)
      {
        if (m_debugOutput)
        {
          std::cout << "\t"  << node << " Descending left." << std::endl;
        }

        delete_node(node->remove_left_child(), destruction);
      }

      if (node->m_pendingRightSubtree != nullptr)
      {
        node->m_pendingRightSubtree->m_cancelled = true;
        node->m_pendingRightSubtree = nullptr;
      }

      if (node->get_right_child() != nullptr)
      {
        if (m_debugOutput)
        {
          std::cout << "\t"  << node << " Descending right." << std::endl;
        }

        delete_node(node->remove_right_child(), destruction);
      }

      if (node->m_pendingMiddleSubtree != nullptr)
      {
        node->m_pendingMiddleSubtree->m_cancelled = true;
        node->m_pendingMiddleSubtree = nullptr;
      }

      if (node->has_middle_children() && m_debugOutput)
      {
        std::cout << "\t"  << node << " Descending middle." << std::endl;
      }

      while (node->has_middle_children())
      {
        delete_node(node->pop_middle_child(), destruction);
      }

      // At the point of destruction, the `Triangle_mesh` referenced may have gone out of scope,
      // we wish to distinguish between deletion with an assumed reference
      // to the original `Triangle_mesh`, and deletion without
      if (!node->is_root_node() && !destruction)
      {
        std::size_t entryHalfEdgeIndex = get(m_halfedgeIndexMap, node->entry_edge());

        if (m_vertexOccupiers[entryHalfEdgeIndex].first == node)
        {
          m_vertexOccupiers[entryHalfEdgeIndex].first = nullptr;

          std::size_t targetVertexIndex = get(m_vertexIndexMap, node->target_vertex());

          if (m_closestToVertices[targetVertexIndex].first == node)
          {
            m_closestToVertices[targetVertexIndex].first = nullptr;
          }
        }
      }

      delete node;
    }

    node_deleted();
  }

  void set_vertex_types()
  {
    for(vertex_descriptor v : vertices(m_graph))
    {
      std::size_t vertexIndex = get(m_vertexIndexMap, v);
      m_vertexIsPseudoSource[vertexIndex] = !internal::is_isolated(v, m_graph) &&
                                            (is_saddle_vertex(v) || is_boundary_vertex(v));
    }
  }

  bool is_saddle_vertex(const vertex_descriptor v) const
  {
    return m_traits.is_saddle_vertex_object()(v, m_graph, m_vertexPointMap);
  }

  bool is_boundary_vertex(const vertex_descriptor v) const
  {
    return bool(is_border(v, m_graph));
  }

  void delete_all_nodes()
  {
    for (std::size_t i = 0; i < m_rootNodes.size(); ++i)
    {
      delete_node(m_rootNodes[i].first, true);
    }
  }

  void reset_algorithm(const bool clearFaceLocations = true)
  {
    m_closestToVertices.assign(num_vertices(m_graph), Node_distance_pair(nullptr, FT(-1)));
    m_vertexOccupiers.assign(num_halfedges(m_graph), Node_distance_pair(nullptr, FT(-1)));

    while (!m_expansionPriqueue.empty())
    {
      delete m_expansionPriqueue.top();
      m_expansionPriqueue.pop();
    }

    if (clearFaceLocations)
    {
      m_faceLocations.clear();
      m_firstNewSourcePoint = m_faceLocations.end();
      m_deletedSourceLocations.clear();
    }

    delete_all_nodes();
    m_rootNodes.clear();
    m_vertexIsPseudoSource.assign(num_vertices(m_graph), false);

#if !defined(NDEBUG)
    m_currentNodeCount = 0;
    m_peakNodeCount = 0;
    m_queueAtPeakNodes = 0;
    m_peakQueueSize = 0;
    m_nodesAtPeakQueue = 0;
#endif

  }

  template <class Visitor>
  void visit_shortest_path(const Cone_tree_node* startNode,
                           const Point_2& startLocation,
                           Visitor& visitor)
  {
    typename Traits::Compute_parametric_distance_along_segment_2 parametric_distance_along_segment_2(m_traits.compute_parametric_distance_along_segment_2_object());
    typename Traits::Construct_ray_2 construct_ray_2(m_traits.construct_ray_2_object());
    typename Traits::Construct_line_2 construct_line_2(m_traits.construct_line_2_object());
    typename Traits::Construct_source_2 construct_source_2(m_traits.construct_source_2_object());
    typename Traits::Construct_target_2 construct_target_2(m_traits.construct_target_2_object());
    typename Traits::Intersect_2 intersect_2(m_traits.intersect_2_object());

    const Cone_tree_node* current = startNode;
    Point_2 currentLocation(startLocation);

    while (!current->is_root_node())
    {
      switch (current->node_type())
      {
        case Cone_tree_node::INTERVAL:
        {
          const Segment_2& entrySegment = current->entry_segment();
          const Point_2& currentSourceImage = current->source_image();
          Ray_2 rayToLocation(construct_ray_2(currentSourceImage, currentLocation));

          const auto cgalIntersection = intersect_2(construct_line_2(entrySegment),
                                                    construct_line_2(rayToLocation));

          CGAL_assertion(bool(cgalIntersection));

          const Point_2* result = boost::get<Point_2>(&*cgalIntersection);
          if (!result)
            result = &currentSourceImage;

          FT t0 = parametric_distance_along_segment_2(construct_source_2(entrySegment),
                                                      construct_target_2(entrySegment), *result);

          if (m_debugOutput)
          {
            std::cout << "Current Node: " << current << " , Face = " << current->layout_face() << std::endl;
            halfedge_descriptor he = current->entry_edge();
            std::cout << "Face vertices: ";
            for (std::size_t i = 0; i < 3; ++i)
            {
              std::cout << get(m_vertexIndexMap, source(he, m_graph)) << ",";
              he = next(he, m_graph);
            }
            std::cout << std::endl;
            std::cout << "Current Location: " << currentLocation << std::endl;
            std::cout << "Distance: " << current->distance_to_root(currentLocation) << std::endl;
            std::cout << "Inside cone: " << (current->inside_window(currentLocation) ? "Yes" : "No") << std::endl;
            std::cout << "Current Source: " << current->source_image() << std::endl;
            std::cout << "Current Segment: " << entrySegment << std::endl;
            std::cout << "Current Left Window: " << current->window_left() << "  ,  "
                      << m_traits.compute_parametric_distance_along_segment_2_object()(entrySegment.start(), entrySegment.end(), current->window_left()) << std::endl;
            std::cout << "Current Right Window: " << current->window_right() << "  ,  "
                      << m_traits.compute_parametric_distance_along_segment_2_object()(entrySegment.start(), entrySegment.end(), current->window_right()) << std::endl;
            std::cout << "Current Segment Intersection: " << *result << std::endl;
            std::cout << "Edge: (V" << get(m_vertexIndexMap, source(current->entry_edge(), m_graph))
                      << ", V" << get(m_vertexIndexMap, target(current->entry_edge(), m_graph)) << ")  :  " << t0 << std::endl;
          }

          visitor(current->entry_edge(), t0);

          currentLocation = *result;

          current = current->parent();

        }
          break;
        case Cone_tree_node::VERTEX_SOURCE:
          // This might be a pseudo source
          visitor(target(current->entry_edge(), m_graph));
          currentLocation = current->parent()->target_point();
          current = current->parent();
          break;
        case Cone_tree_node::EDGE_SOURCE:
        case Cone_tree_node::FACE_SOURCE:
          // This is guaranteed to be the final node in any sequence
          visitor(m_rootNodes[current->tree_id()].second->first,
                  m_rootNodes[current->tree_id()].second->second);
          current = current->parent();
          break;
        default:
          CGAL_assertion(false && "Unhandled node type found in tree");
      }
    }
  }

  void add_to_face_list(Cone_tree_node* node)
  {
    if (!node->is_root_node() && !node->is_null_face())
    {
      std::size_t faceIndex = get(m_faceIndexMap, node->current_face());
      m_faceOccupiers[faceIndex].push_back(node);
    }

    if (node->get_left_child() != nullptr)
    {
      add_to_face_list(node->get_left_child());
    }

    if (node->get_right_child() != nullptr)
    {
      add_to_face_list(node->get_right_child());
    }

    for (std::size_t i = 0; i < node->num_middle_children(); ++i)
    {
      add_to_face_list(node->get_middle_child(i));
    }
  }

  Point_2 face_location_with_normalized_coordinates(const Cone_tree_node* node,
                                                    const Barycentric_coordinates& location) const
  {
    return construct_barycenter_in_triangle_2(node->layout_face(), localized_coordiate(node, location));
  }

  Barycentric_coordinates localized_coordiate(const Cone_tree_node* node,
                                              const Barycentric_coordinates& location) const
  {
    return shifted_coordinates(location, node->edge_face_index());
  }

  Barycentric_coordinates shifted_coordinates(const Barycentric_coordinates& location,
                                              const std::size_t shift) const
  {
    typename Traits::Construct_barycentric_coordinates_weight cbcw(m_traits.construct_barycentric_coordinates_weight_object());
    typename Traits::Construct_barycentric_coordinates cbc(m_traits.construct_barycentric_coordinates_object());
    return cbc(cbcw(location, shift), cbcw(location, (shift + 1) % 3), cbcw(location, (shift + 2) % 3));
  }

  std::pair<Node_distance_pair, Barycentric_coordinates> nearest_on_face(const face_descriptor f,
                                                                         const Barycentric_coordinates& location) const
  {
    typename Traits::Construct_barycentric_coordinates cbc(m_traits.construct_barycentric_coordinates_object());

    const std::size_t faceIndex = get(m_faceIndexMap, f);

    Cone_tree_node* closest = nullptr;
    FT closestDistance = 0;

    const std::vector<Cone_tree_node*>& currentFaceList = m_faceOccupiers[faceIndex];

    for (std::size_t i = 0; i < currentFaceList.size(); ++i)
    {
      Cone_tree_node* current = currentFaceList[i];

      if (closest != nullptr && current->distance_from_source_to_root() >= closestDistance)
      {
        continue;
      }

      const Point_2 locationInContext = face_location_with_normalized_coordinates(current, location);

      if (current->inside_window(locationInContext))
      {
        FT currentDistance = current->distance_to_root(locationInContext);

        if (closest == nullptr || currentDistance < closestDistance)
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
      return std::make_pair(Node_distance_pair((Cone_tree_node*)nullptr, FT(0)), cbc(FT(0), FT(0), FT(0)));
    }
  }

  std::pair<Node_distance_pair, Barycentric_coordinates> nearest_to_location(const face_descriptor f,
                                                                             const Barycentric_coordinates& location) const
  {
    typename Traits::Construct_barycentric_coordinates_weight cbcw(m_traits.construct_barycentric_coordinates_weight_object());
    typename Traits::Construct_barycentric_coordinates cbc(m_traits.construct_barycentric_coordinates_object());
    typename Traits::Classify_barycentric_coordinates classify_barycentric_coordinates(m_traits.classify_barycentric_coordinates_object());

    std::size_t associatedEdge;
    CGAL::Surface_mesh_shortest_paths_3::Barycentric_coordinates_type type;
    std::tie(type, associatedEdge) = classify_barycentric_coordinates(location);

    switch (type)
    {
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_BOUNDED_SIDE:
        return nearest_on_face(f, location);
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_BOUNDARY:
        {
          halfedge_descriptor he = halfedge(f, m_graph);
          for (std::size_t i = 0; i < associatedEdge; ++i)
          {
            he = next(he, m_graph);
          }

          std::pair<Node_distance_pair,Barycentric_coordinates> mainFace = nearest_on_face(f, location);

          halfedge_descriptor oppositeHalfedge = opposite(he, m_graph);
          if (!CGAL::is_border(oppositeHalfedge, m_graph))
          {
              std::size_t oppositeIndex = Surface_mesh_shortest_paths_3::internal::edge_index(oppositeHalfedge, m_graph);

              FT oppositeLocationCoords[3] = { FT(0), FT(0), FT(0) };
              oppositeLocationCoords[oppositeIndex] = cbcw(location, (associatedEdge + 1) % 3);
              oppositeLocationCoords[(oppositeIndex + 1) % 3] = cbcw(location, associatedEdge);
              Barycentric_coordinates oppositeLocation(cbc(oppositeLocationCoords[0], oppositeLocationCoords[1], oppositeLocationCoords[2]));
              std::pair<Node_distance_pair,Barycentric_coordinates> otherFace = nearest_on_face(face(oppositeHalfedge, m_graph), oppositeLocation);

              if (mainFace.first.first == nullptr)
              {
                return otherFace;
              }
              else if (otherFace.first.first == nullptr)
              {
                return mainFace;
              }
              else
              {
                return mainFace.first.second < otherFace.first.second ? mainFace : otherFace;
              }
          }

          return mainFace;
        }
        break;
      case CGAL::Surface_mesh_shortest_paths_3::BARYCENTRIC_COORDINATES_ON_VERTEX:
        {
          halfedge_descriptor he = halfedge(f, m_graph);

          for (std::size_t i = 0; i < associatedEdge; ++i)
          {
            he = next(he, m_graph);
          }

          vertex_descriptor vertex = source(he, m_graph);

          return std::make_pair(m_closestToVertices[get(m_vertexIndexMap, vertex)], cbc(FT(0), FT(0), FT(1)));
        }
        break;

      default:
        CGAL_assertion(false && "Invalid face location");
        return std::pair<Node_distance_pair, Barycentric_coordinates>();
    }
  }

  static bool cone_comparator(const Cone_tree_node* lhs, const Cone_tree_node* rhs)
  {
    return lhs->distance_from_source_to_root() < rhs->distance_from_source_to_root();
  }

  template <class InputIterator>
  Source_point_iterator add_source_points_internal(InputIterator begin, InputIterator end,
                                                   const vertex_descriptor)
  {
    Source_point_iterator firstAdded;

    for (InputIterator it = begin; it != end; ++it)
    {
      Source_point_iterator added = add_source_point(face_location(*it));

      if (it == begin)
      {
        firstAdded = added;
      }
    }

    return firstAdded;
  }

  template <class InputIterator>
  Source_point_iterator add_source_points_internal(InputIterator begin, InputIterator end,
                                                   const Face_location&)
  {
    Source_point_iterator firstAdded;

    for (InputIterator it = begin; it != end; ++it)
    {
      Source_point_iterator added = add_source_point(it->first, it->second);

      if (it == begin)
      {
        firstAdded = added;
      }
    }

    return firstAdded;
  }

  void construct_sequence_tree_internal()
  {
    reset_algorithm(false);
    set_vertex_types();

    m_vertexOccupiers.assign(num_halfedges(m_graph), Node_distance_pair(nullptr, FT(-1)));
    m_closestToVertices.assign(num_vertices(m_graph), Node_distance_pair(nullptr, FT(-1)));

    if (m_debugOutput)
    {
      std::size_t numVertices = 0;
      for (vertex_descriptor v : vertices(m_graph))
      {
        std::cout << "Vertex#" << numVertices
                  << ": p = " << get(m_vertexPointMap, v)
                  << " , Saddle Vertex: " << (is_saddle_vertex(v) ? "yes" : "no")
                  << " , Boundary Vertex: " << (is_boundary_vertex(v) ? "yes" : "no") << std::endl;
        ++numVertices;
      }
    }

    if (m_debugOutput)
    {
      std::size_t numFaces = 0;
      for (face_descriptor f : faces(m_graph))
      {
        std::cout << "Face#" << numFaces << ": Vertices = (";
        ++numFaces;
        halfedge_descriptor faceEdgesStart = halfedge(f, m_graph);
        halfedge_descriptor faceEdgesCurrent = faceEdgesStart;

        do
        {
          std::cout << "V" << get(m_vertexIndexMap, source(faceEdgesCurrent, m_graph));

          faceEdgesCurrent = next(faceEdgesCurrent, m_graph);

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

    for (typename Source_point_list::iterator it = m_faceLocations.begin(); it != m_faceLocations.end(); ++it)
    {
      if (m_debugOutput)
      {
        std::cout << "Root: F" << get(m_faceIndexMap, it->first)
                  << " , bar " << it->second[0] << " " << it->second[1] << " " << it->second[2] << " "
                  << ", pos " << point(it->first, it->second) << std::endl;
      }

      expand_root(it->first, it->second, Source_point_iterator(it));
    }

    if (m_debugOutput)
    {
      std::cout << "PriQ start size = " << m_expansionPriqueue.size() << std::endl;

      std::cout << "Num face locations: " << m_faceLocations.size() << std::endl;
      std::cout << "Num root nodes: " << m_rootNodes.size() << " (Hint: these should be the same size)" << std::endl;

    }

    while (!m_expansionPriqueue.empty())
    {
      if (m_debugOutput)
      {
        std::cout << " -----------------------------------------------------------------------" << std::endl;
        std::cout << " -----------------------------------------------------------------------" << std::endl;
        std::cout << "Num face locations: " << m_faceLocations.size() << std::endl;
        std::cout << "Num root nodes: " << m_rootNodes.size() << " (Hint: these should be the same size)" << std::endl;

        std::cout << "Prio Queue size = " << m_expansionPriqueue.size() << std::endl;
        std::cout << "Queue:" << std::endl;
        auto duplicate_queue = m_expansionPriqueue;
        while(duplicate_queue.size() > 0)
        {
          Cone_expansion_event* event = duplicate_queue.top();

          std::cout << "event type: " << event->m_type << " "
                    << " time: " << event->m_distanceEstimate << " ";
          std::cout << "cancelled? " << event->m_cancelled << " " ;

          if (!event->m_cancelled)
          {
            std::cout << " ------ Parent (" << event->m_parent << ") INFO: ";
            std::cout << "EH = (V" << get(m_vertexIndexMap, source(event->m_parent->entry_edge(), m_graph)) << " V"
                                  << get(m_vertexIndexMap, target(event->m_parent->entry_edge(), m_graph)) << ") ";
            std::cout << "Src = (" << event->m_parent->source_image() << ") ";
            std::cout << "Tar = V" << get(m_vertexIndexMap, target(next(event->m_parent->entry_edge(), m_graph), m_graph));
          }

          std::cout << std::endl;

          duplicate_queue.pop();
        }
      }

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
              std::cout << "PseudoSource Expansion: Parent = " << parent
                        << " , Vertex = " << get(m_vertexIndexMap, event->m_parent->target_vertex())
                        << " , Distance = " << event->m_distanceEstimate
                        << " , Level = " << event->m_parent->level() + 1 << std::endl;
            }

            expand_pseudo_source(parent);
            break;
          case Cone_expansion_event::LEFT_CHILD:
            if (m_debugOutput)
            {
              std::cout << "Left Expansion: Parent = " << parent
                        << " Edge = (V" << get(m_vertexIndexMap, source(event->m_parent->left_child_edge(), m_graph))
                        << ", V" << get(m_vertexIndexMap, target(event->m_parent->left_child_edge(), m_graph))
                        << ") , Distance = " << event->m_distanceEstimate
                        << " , Level = " << event->m_parent->level() + 1 << std::endl;
            }

            expand_left_child(parent, event->m_windowSegment);
            break;
          case Cone_expansion_event::RIGHT_CHILD:
            if (m_debugOutput)
            {
              std::cout << "Right Expansion: Parent = " << parent
                        << " , Edge = (V" << get(m_vertexIndexMap, source(event->m_parent->right_child_edge(), m_graph))
                        << ", V" << get(m_vertexIndexMap, target(event->m_parent->right_child_edge(), m_graph))
                        << ") , Distance = " << event->m_distanceEstimate
                        << " , Level = " << event->m_parent->level() + 1 << std::endl;
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
    m_faceOccupiers.resize(num_faces(m_graph));

    for (std::size_t i = 0; i < m_rootNodes.size(); ++i)
    {
      add_to_face_list(m_rootNodes[i].first);
    }

    for (std::size_t i = 0; i < m_faceOccupiers.size(); ++i)
    {
      std::vector<Cone_tree_node*>& currentFaceList = m_faceOccupiers[i];
      std::sort(currentFaceList.begin(), currentFaceList.end(), cone_comparator);
    }

    if (m_debugOutput)
    {
      std::cout << "Closest distances: " << std::endl;

      for (std::size_t i = 0; i < m_closestToVertices.size(); ++i)
      {
        std::cout << "\tVertex = " << i << std::endl;
        std::cout << "\tDistance = " << m_closestToVertices[i].second << " to " << m_closestToVertices[i].first << std::endl;
      }

      std::cout << std::endl;

      for (std::size_t i = 0; i < m_faceOccupiers.size(); ++i)
      {
        std::cout << "\tFace = " << i << std::endl;
        std::cout << "\t#Occupiers = " << m_faceOccupiers[i].size() << std::endl;
      }

      std::cout << std::endl << "Done!" << std::endl;
    }

    m_firstNewSourcePoint = m_faceLocations.end();
    m_deletedSourceLocations.clear();

  }

public:

  /// \name Constructors
  /// @{

  /*!
  \brief creates a shortest paths object using `tm` as input.

  Equivalent to `Surface_mesh_shortest_path(tm, get(boost::vertex_index, tm), get(boost::halfedge_index, tm),
                                            get(boost::face_index, tm), get(CGAL::vertex_point, tm), traits)`.

  Internal property maps must be available and initialized.

  \sa \link BGLGraphExternalIndices `CGAL::set_halfedgeds_items_id()`\endlink
  */
  Surface_mesh_shortest_path(const Triangle_mesh& tm,
                             const Traits& traits = Traits())
    : m_traits(traits)
    , m_graph(tm)
    , m_vertexIndexMap(get(boost::vertex_index, tm))
    , m_halfedgeIndexMap(get(boost::halfedge_index, tm))
    , m_faceIndexMap(get(boost::face_index, tm))
    , m_vertexPointMap(get(CGAL::vertex_point, m_graph))
    , m_debugOutput(false)
  {
    reset_algorithm();
  }

  /*!
  \brief creates a shortest paths object using `tm` as input.

  \details No copy of the `Triangle_mesh` is made, only a reference to the `tm` is held.

  \param tm The surface mesh to compute shortest paths on.  Note that it must be triangulated.

  \param vertexIndexMap Property map associating an id to each vertex, from 0 to `num_vertices(tm) - 1`.

  \param halfedgeIndexMap Property map associating an id to each halfedge, from 0 to `num_halfedges(tm) - 1`.

  \param faceIndexMap Property map associating an id to each face, from 0 to `num_faces(tm) - 1`.

  \param vertexPointMap Property map used to access the points associated to each vertex of the graph.

  \param traits Optional instance of the traits class to use.
  */
  Surface_mesh_shortest_path(const Triangle_mesh& tm,
                             Vertex_index_map vertexIndexMap,
                             Halfedge_index_map halfedgeIndexMap,
                             Face_index_map faceIndexMap,
                             Vertex_point_map vertexPointMap,
                             const Traits& traits = Traits())
    : m_traits(traits)
    , m_graph(tm)
    , m_vertexIndexMap(vertexIndexMap)
    , m_halfedgeIndexMap(halfedgeIndexMap)
    , m_faceIndexMap(faceIndexMap)
    , m_vertexPointMap(vertexPointMap)
    , m_debugOutput(false)
  {
    reset_algorithm();
  }

  /// @}

  /// \cond

  ~Surface_mesh_shortest_path()
  {
    delete_all_nodes();

#if !defined(NDEBUG)
    if (m_debugOutput)
    {
      std::cout << "Final node count: " << m_currentNodeCount << std::endl;
      std::cout << "Peak node count: " << m_peakNodeCount << std::endl;
    }
#endif
  }

  /// \endcond

  /// \name Addition and Removal of Source Points
  /// @{

  /*!
  \brief adds `v` as a source for the shortest path queries.

  \details No change to the internal shortest paths data structure occurs
  until either `Surface_mesh_shortest_path::build_sequence_tree()` or
  the first shortest path query is done.

  \return An iterator to the source point added
  */
  Source_point_iterator add_source_point(vertex_descriptor v)
  {
    Face_location location = face_location(v);

    if (m_debugOutput)
    {
      std::cout << "Face location from V" << get(m_vertexIndexMap, v) << " is F" << get(m_faceIndexMap, location.first) << " "
                << location.second[0] << " " << location.second[1] << " " << location.second[2] << std::endl;
    }

    return add_source_point(location);
  }

  /*!
  \brief adds a point inside the face `f` as a source for the shortest path queries.

  \details No change to the internal shortest paths data structure occurs
  until either `Surface_mesh_shortest_path::build_sequence_tree()` or
  the first shortest path query is done.

  \param f A face of the input face graph
  \param location Barycentric coordinates in face `f` specifying the source point.
  \return An iterator to the source point added
  */
  Source_point_iterator add_source_point(const face_descriptor f,
                                         const Barycentric_coordinates& location)
  {
    return add_source_point(std::make_pair(f, location));
  }

  /*!
  \brief adds a point inside a face as a source for the shortest path queries,
  equivalent to `Surface_mesh_shortest_path::add_source_point(location.first, location.second);`
  */
  Source_point_iterator add_source_point(const Face_location& location)
  {
    if (m_debugOutput)
    {
      std::cout << "Add source point at position " << point(location.first, location.second) << std::endl;
    }

    Source_point_underlying_iterator added = m_faceLocations.insert(m_faceLocations.end(), location);

    if (m_firstNewSourcePoint == m_faceLocations.end())
    {
      m_firstNewSourcePoint = added;
    }

    return Source_point_iterator(added);
  }

  /*!
  \brief adds a range of points as sources for the shortest path queries.

  \details No change to the internal shortest paths data structure occurs
  until either `Surface_mesh_shortest_path::build_sequence_tree()` or
  the first shortest path query is done.

  \tparam InputIterator A `ForwardIterator` which dereferences to either `Surface_mesh_shortest_path::Face_location`,
                        or `Surface_mesh_shortest_path::vertex_descriptor`.

  \param begin iterator to the first in the list of source point locations.
  \param end iterator to one past the end of the list of source point locations.
  \return An iterator to the first source point added.
  */
  template <class InputIterator>
  Source_point_iterator add_source_points(InputIterator begin, InputIterator end)
  {
    return add_source_points_internal(begin, end, typename std::iterator_traits<InputIterator>::value_type());
  }

  /*!
  \brief removes a source point for the shortest path queries.

  \details No change to the internal shortest paths data structure occurs
  until either `Surface_mesh_shortest_path::build_sequence_tree()` or
  the first shortest path query is done.
  Behaviour is undefined if the source point `it` was already removed.

  \param it iterator to the source point to be removed
  */
  void remove_source_point(Source_point_iterator it)
  {
    if (it == m_firstNewSourcePoint)
    {
      ++m_firstNewSourcePoint;
    }

    m_deletedSourceLocations.splice(m_deletedSourceLocations.begin(), m_faceLocations, it.m_iterator);
  }

  /*!
  \brief removes all source points for the shortest path queries.

  \details No change to the internal shortest paths data structure occurs
  until either `Surface_mesh_shortest_path::build_sequence_tree()` or
  the first shortest path query is done.
  For a version which deletes all data immediately, use `clear()` instead.
  */
  void remove_all_source_points()
  {
    m_deletedSourceLocations.splice(m_deletedSourceLocations.begin(), m_faceLocations,
                                    m_faceLocations.begin(), m_faceLocations.end());
    m_firstNewSourcePoint = m_faceLocations.end();
  }

  /// @}

  /// \name Creation and Destruction of the Shortest Paths Sequence Tree
  /// @{

  /*!
  \brief Computes all pending changes to the internal sequence tree

  \details A call to this method will only trigger a computation only if some
  change to the set of source points occurred since the last time
  the sequence tree was computed.
  */
  void build_sequence_tree()
  {
    if (changed_since_last_build())
    {
      construct_sequence_tree_internal();
    }
  }

  /*!
  \brief removes all data, the class is as if it was constructed.

  \details All internal containers are cleared  and the internal
  sequence tree is also cleared.  For a version which defers deletion until
  it is necessary, use `Surface_mesh_shortest_path::remove_all_source_points()`.
  */
  void clear()
  {
    reset_algorithm();
  }

  /// @}

  /// \name Accessors
  /// @{

  /*!
  \brief returns an iterator to the first source point location

  \details The elements will appear in the order they were inserted to the
  structure by calls to `add_source_point()` or `add_source_points()`.  Deleted
  points will not appear in the sequence.

  \return An iterator to the first of the stored source points.
  */
  Source_point_iterator source_points_begin() const
  {
    // It is a feature of C++11 that `const_iterator` may be used in calls to `erase()`, however
    // in order to support C++98, we must use `iterator`.  Semantically, this is correct, but
    // I must cast away the const-ness to hide the internal ugliness
    return Source_point_iterator(const_cast<std::list<Face_location>&>(m_faceLocations).begin());
  }

  /*!
  \brief returns an iterator to one past the last source point location

  \return An iterator to one past-the-end in the list of stored source points.
  */
  Source_point_iterator source_points_end() const
  {
    return Source_point_iterator(const_cast<std::list<Face_location>&>(m_faceLocations).end());
  }

  /*!
  \brief returns the total number of source points used for the shortest path queries.
  */
  std::size_t number_of_source_points() const
  {
    return m_faceLocations.size();
  }

  /*!
  \brief determines if the internal sequence tree is valid (already built and no new source point has been added).

  \return true if the structure needs to be rebuilt, false otherwise
  */
  bool changed_since_last_build() const
  {
    return m_firstNewSourcePoint != m_faceLocations.end() || !m_deletedSourceLocations.empty();
  }

  /// @}

  /// \name Shortest Distance Queries
  /// @{

  /*!
  \brief Computes the shortest surface distance from a vertex to any source point

  \param v A vertex of the input face graph
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  Shortest_path_result shortest_distance_to_source_points(const vertex_descriptor v)
  {
    build_sequence_tree();

    const Node_distance_pair& result = m_closestToVertices[get(m_vertexIndexMap, v)];
    const Cone_tree_node* current = result.first;

    if (current)
    {
      return std::make_pair(result.second, m_rootNodes[current->tree_id()].second);
    }
    else
    {
      return std::make_pair(FT(-1), source_points_end());
    }
  }

  /*!
  \brief Computes the shortest surface distance from any surface location to any source point

  \param f A face of the input face graph
  \param location Barycentric coordinates of the query point on face `f`
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  Shortest_path_result shortest_distance_to_source_points(const face_descriptor f,
                                                          const Barycentric_coordinates& location)
  {
    build_sequence_tree();

    const std::pair<Node_distance_pair, Barycentric_coordinates>& result = nearest_to_location(f, location);
    const Cone_tree_node* current = result.first.first;

    if (current)
    {
      return std::make_pair(result.first.second, m_rootNodes[current->tree_id()].second);
    }
    else
    {
      return std::make_pair(FT(-1), source_points_end());
    }
  }

  /// @}

  /// \name Shortest Path Sequence Queries
  /// @{

  /*!
  \brief visits the sequence of edges, vertices and faces traversed by the shortest path
  from a vertex to any source point.

  \details Visits simplices, starting from the query vertex, back to
  the nearest source point. If no shortest path could be found (for example,
  the surface is disconnected), then no calls to the visitor will be made
  (not even for the query vertex).

  \param v A vertex of the input face graph
  \param visitor A model of `SurfaceMeshShortestPathVisitor` to receive the shortest path
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  template <class Visitor>
  Shortest_path_result
  shortest_path_sequence_to_source_points(const vertex_descriptor v,
                                          Visitor& visitor)
  {
    build_sequence_tree();

    const Node_distance_pair& result = m_closestToVertices[get(m_vertexIndexMap, v)];
    const Cone_tree_node* current = result.first;

    if (current)
    {
      visitor(v);
      visit_shortest_path(current, current->target_point(), visitor);
      return std::make_pair(result.second, m_rootNodes[current->tree_id()].second);
    }
    else
    {
      return std::make_pair(FT(-1), source_points_end());
    }
  }

  /*!
  \brief visits the sequence of edges, vertices and faces traversed by the shortest path
  from any surface location to any source point.

  \details Visits simplices, starting from the query point, back to
  the nearest source point. If no shortest path could be found (for example,
  the surface is disconnected), then no calls to the visitor will be made
  (not even for the query point).

  \param f A face of the input face graph
  \param location Barycentric coordinates of the query point on face `f`
  \param visitor A model of `SurfaceMeshShortestPathVisitor` to receive the shortest path
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  template <class Visitor>
  Shortest_path_result
  shortest_path_sequence_to_source_points(const face_descriptor f,
                                          const Barycentric_coordinates& location,
                                          Visitor& visitor)
  {
    build_sequence_tree();

    std::pair<Node_distance_pair, Barycentric_coordinates> result = nearest_to_location(f, location);
    Cone_tree_node* current = result.first.first;

    if (current)
    {
      Point_2 locationInContext = construct_barycenter_in_triangle_2(current->layout_face(), result.second);
      visitor(f, location);
      visit_shortest_path(current, locationInContext, visitor);
      return std::make_pair(result.first.second, m_rootNodes[current->tree_id()].second);
    }
    else
    {
      return std::make_pair(FT(-1), source_points_end());
    }
  }

  /// @}

  /// \name Shortest Path Point Queries
  /// @{

  /*!
  \brief Computes the sequence of points in the shortest path along the
    surface of the input face graph from the given vertex to the closest
    source point.

  \param v A vertex of the input face graph
  \param output An OutputIterator to receive the shortest path points as `Point_3` objects
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  template <class OutputIterator>
  Shortest_path_result
  shortest_path_points_to_source_points(const vertex_descriptor v, OutputIterator output)
  {
    build_sequence_tree();

    Point_path_visitor_wrapper<OutputIterator> wrapper(*this, output);
    return shortest_path_sequence_to_source_points(v, wrapper);
  }

  /*!
  \brief Computes the sequence of points in the shortest path along the
    surface of the input face graph from the given query location to the closest
    source point.

  \param f A face of on the input face graph
  \param location The barycentric coordinates of the query point on face `f`
  \param output An OutputIterator to receive the shortest path points as `Point_3` objects
  \return A pair, containing the distance to the source point, and an
    iterator to the source point.  If no source point was reachable (can
    occur when the graph is disconnected), the distance will be a negative
    value and the source point iterator will be equal to `source_points_end()`.
  */
  template <class OutputIterator>
  Shortest_path_result
  shortest_path_points_to_source_points(const face_descriptor f,
                                        const Barycentric_coordinates& location, OutputIterator output)
  {
    build_sequence_tree();

    Point_path_visitor_wrapper<OutputIterator> wrapper(*this, output);
    return shortest_path_sequence_to_source_points(f, location, wrapper);
  }

  /// @}

  /// \name Surface Point Constructions
  /// @{

  /*!
  \brief returns the 3-dimensional coordinates at the barycentric coordinates
    of the given face.

  \details The following static overloads are also available:
    - `static Point_3 point(face_descriptor f, Barycentric_coordinates location, const Triangle_mesh& tm,
                            const Traits& traits = Traits())`
    - `static Point_3 point(face_descriptor f, Barycentric_coordinates location, const Triangle_mesh& tm,
                            Vertex_point_map vertexPointMap, const Traits& traits = Traits())`

  \param f A face of on the input face graph
  \param location The barycentric coordinates of the query point on face `f`
  */
  Point_3 point(const face_descriptor f,
                const Barycentric_coordinates& location) const
  {
    return point(f, location, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  static Point_3 point(const face_descriptor f,
                       const Barycentric_coordinates& location,
                       const Triangle_mesh& tm,
                       const Traits& traits = Traits())
  {
    using boost::get;
    return point(f, location, tm, get(CGAL::vertex_point, tm), traits);
  }

  static Point_3 point(const face_descriptor f,
                       const Barycentric_coordinates& location,
                       const Triangle_mesh& tm,
                       Vertex_point_map vertexPointMap,
                       const Traits& traits = Traits())
  {
    return construct_barycenter_in_triangle_3(triangle_from_face(f, tm, vertexPointMap), location, traits);
  }

  /// \endcond

  /*!
  \brief returns the 3-dimensional coordinates at the parametric location
    along the given edge.

  \details The following static overloads are also available:
    - `static Point_3 point(halfedge_descriptor edge, FT t, const Triangle_mesh& tm, const Traits& traits = Traits())`
    - `static Point_3 point(halfedge_descriptor edge, FT t, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())`

  \param edge An edge of the input face graph
  \param t The parametric distance along edge of the desired point
  */
  Point_3 point(const halfedge_descriptor edge, const FT t) const
  {
    return point(edge, t, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  static Point_3 point(const halfedge_descriptor edge, const FT t,
                       const Triangle_mesh& tm,
                       const Traits& traits = Traits())
  {
    using boost::get;
    return point(edge, t, tm, get(CGAL::vertex_point, tm), traits);
  }

  static Point_3 point(const halfedge_descriptor edge, const FT t,
                       const Triangle_mesh& tm,
                       Vertex_point_map vertexPointMap,
                       const Traits& traits = Traits())
  {
    typename Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());

    // Note: the parameter t is meant to be the weighted coordinates on the _endpoint_ (i.e. target) of the segment
    return construct_barycenter_3(get(vertexPointMap, target(edge, tm)), t, get(vertexPointMap, source(edge, tm)));
  }

  /// \endcond

  /*!
  \brief returns the 3-dimensional coordinates of the given vertex.

  \param v A vertex of the input face graph
  */
  decltype(auto) point(const vertex_descriptor v) const
  {
    return get(m_vertexPointMap, v);
  }

  /// \cond

  static decltype(auto) point(const vertex_descriptor v,
                              const Triangle_mesh& tm)
  {
    return get(CGAL::vertex_point, tm, v);
  }

  /// \endcond

  /// @}

  /// \name Surface Face Location Constructions
  /// @{

  /*!
  \brief returns the location of the given vertex as a `Face_location`

  \details The following static overload is also available:
    - `static Face_location face_location(vertex_descriptor vertex, const Triangle_mesh& tm, const Traits& traits = Traits())`

  \param vertex A vertex of the input face graph
  */
  Face_location face_location(const vertex_descriptor vertex) const
  {
    return face_location(vertex, m_graph, m_traits);
  }

  /// \cond

  static Face_location face_location(const vertex_descriptor vertex,
                                     const Triangle_mesh& tm,
                                     const Traits& traits = Traits())
  {
    typename Traits::Construct_barycentric_coordinates construct_barycentric_coordinates(traits.construct_barycentric_coordinates_object());
    halfedge_descriptor hinit=halfedge(vertex, tm);
    while (is_border(hinit, tm))
      hinit = opposite(next(hinit, tm), tm);

    halfedge_descriptor he = next(hinit, tm);
    face_descriptor locationFace = face(he, tm);
    std::size_t edgeIndex = Surface_mesh_shortest_paths_3::internal::edge_index(he, tm);

    FT coords[3] = { FT(0), FT(0), FT(0) };

    coords[edgeIndex] = FT(1);

    return Face_location(locationFace, construct_barycentric_coordinates(coords[0], coords[1], coords[2]));
  }

  /// \endcond

  /*!
  \brief returns a location along the given edge as a `Face_location`.

  \details The following static overload is also available:
    - `static Face_location face_location(halfedge_descriptor he, FT t, const Triangle_mesh& tm, const Traits& traits = Traits())`

  \param he A halfedge of the input face graph
  \param t Parametric distance of the desired point along `he`
  */
  Face_location face_location(const halfedge_descriptor he, const FT t) const
  {
    return face_location(he, t, m_graph, m_traits);
  }

  /// \cond

  static Face_location face_location(const halfedge_descriptor he, FT t,
                                     const Triangle_mesh& tm,
                                     const Traits& traits = Traits())
  {
    typename Traits::Construct_barycentric_coordinates cbc(traits.construct_barycentric_coordinates_object());
    face_descriptor locationFace = face(he, tm);
    std::size_t edgeIndex = Surface_mesh_shortest_paths_3::internal::edge_index(he, tm);

    const FT oneMinusT(FT(1) - t);

    FT coords[3] = { FT(0), FT(0), FT(0) };

    coords[edgeIndex] = oneMinusT;
    coords[(edgeIndex + 1) % 3] = t;

    return Face_location(locationFace, cbc(coords[0], coords[1], coords[2]));
  }

  /// \endcond

  /// @}

  /// \name Nearest Face Location Queries
  /// @{

  /*!
  \brief returns the nearest face location to the given point.
    Note that this will (re-)build an `AABB_tree` on each call. If you need
    to  call this function more than once, use `build_aabb_tree()` to cache a
    copy of the `AABB_tree`, and use the overloads of this function
    that accept a reference to an `AABB_tree` as input.

  \details The following static overload is also available:
    - `static Face_location locate(const %Point_3& p, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())`

  \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.

  \param p Point to locate on the input face graph
  */
  template <class AABBTraits>
  Face_location locate(const Point_3& p) const
  {
    return locate<AABBTraits>(p, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  template <class AABBTraits>
  static Face_location locate(const Point_3& location,
                              const Triangle_mesh& tm,
                              Vertex_point_map vertexPointMap,
                              const Traits& traits = Traits())
  {
    AABB_tree<AABBTraits> tree;
    build_aabb_tree(tm, tree, vertexPointMap);
    return locate(location, tree, tm, vertexPointMap, traits);
  }

  /// \endcond

  /*!
  \brief returns the face location nearest to the given point.

  \details The following static overload is also available:
    - static Face_location locate(const %Point_3& p, const AABB_tree<AABBTraits>& tree, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())

  \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.

  \param p Point to locate on the input face graph
  \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
  */
  template <class AABBTraits>
  Face_location locate(const Point_3& p,
                       const AABB_tree<AABBTraits>& tree) const
  {
    return locate(p, tree, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  template <class AABBTraits>
  static Face_location locate(const Point_3& location,
                              const AABB_tree<AABBTraits>& tree,
                              const Triangle_mesh& tm,
                              Vertex_point_map vertexPointMap,
                              const Traits& traits = Traits())
  {
    typename Traits::Construct_barycentric_coordinates_in_triangle_3 cbcit3(traits.construct_barycentric_coordinates_in_triangle_3_object());
    typename AABB_tree<AABBTraits>::Point_and_primitive_id result = tree.closest_point_and_primitive(location);

    face_descriptor f = result.second;
    Barycentric_coordinates b = cbcit3(triangle_from_face(f, tm, vertexPointMap), result.first);
    return Face_location(f, b);
  }

  /// \endcond

  /*!
  \brief returns the face location along `ray` nearest to its source point.
    Note that this will (re-)build an `AABB_tree` on each call. If you need
    to  call this function more than once, use `build_aabb_tree()` to cache a
    copy of the `AABB_tree`, and use the overloads of this function
    that accept a reference to an `AABB_tree` as input.

  \details The following static overload is also available:
    - `static Face_location locate(const %Ray_3& ray, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())`

  \tparam AABBTraits A model of `AABBTraits` used to define an `AABB_tree`.

  \param ray Ray to intersect with the input face graph
  */
  template <class AABBTraits>
  Face_location locate(const Ray_3& ray) const
  {
    return locate<AABBTraits>(ray, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  template <class AABBTraits>
  static Face_location locate(const Ray_3& ray,
                              const Triangle_mesh& tm,
                              Vertex_point_map vertexPointMap,
                              const Traits& traits = Traits())
  {
    AABB_tree<AABBTraits> tree;
    build_aabb_tree(tm, tree, vertexPointMap);
    return locate(ray, tree, tm, vertexPointMap, traits);
  }

  /// \endcond

  /*!
  \brief returns the face location along `ray` nearest to
    its source point.

  \details The following static overload is also available:
    - static Face_location locate(const %Ray_3& ray, const AABB_tree<AABBTraits>& tree, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())

  \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.

  \param ray Ray to intersect with the input face graph
  \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
  */
  template <class AABBTraits>
  Face_location locate(const Ray_3& ray,
                       const AABB_tree<AABBTraits>& tree) const
  {
    return locate(ray, tree, m_graph, m_vertexPointMap, m_traits);
  }

  /// \cond

  template <class AABBTraits>
  static Face_location locate(const Ray_3& ray,
                              const AABB_tree<AABBTraits>& tree,
                              const Triangle_mesh& tm,
                              Vertex_point_map vertexPointMap,
                              const Traits& traits = Traits())
  {
    typedef AABB_tree<AABBTraits> AABB_face_graph_tree;
    typename Traits::Construct_barycentric_coordinates_in_triangle_3 cbcit3(traits.construct_barycentric_coordinates_in_triangle_3_object());
    typename Traits::Construct_barycentric_coordinates cbc(traits.construct_barycentric_coordinates_object());
    typename Traits::Compute_squared_distance_3 csd3(traits.compute_squared_distance_3_object());
    typedef typename AABB_face_graph_tree::template Intersection_and_primitive_id<Ray_3>::Type Intersection_type;
    typedef boost::optional<Intersection_type> Ray_intersection;

    std::vector<Ray_intersection> intersections;

    tree.all_intersections(ray, std::back_inserter(intersections));

    bool foundOne = false;
    FT nearestDistance = 0;
    Point_3 nearestPoint = CGAL::ORIGIN;
    face_descriptor nearestFace;

    for (std::size_t i = 0; i < intersections.size(); ++i)
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
      Barycentric_coordinates b = cbcit3(triangle_from_face(nearestFace, tm, vertexPointMap), nearestPoint);
      return Face_location(nearestFace, b);
    }
    else
    {
      return Face_location(Graph_traits::null_face(), cbc(FT(0), FT(0), FT(0)));
    }
  }

  /// \endcond

  /*!
  \brief creates an `AABB_tree` suitable for use with `locate`.

  \details The following static overload is also available:
    - `static void build_aabb_tree(const Triangle_mesh& tm, AABB_tree<AABBTraits>& outTree)`

  \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.

  \param outTree Output parameter to store the computed `AABB_tree`
  */
  template <class AABBTraits>
  void build_aabb_tree(AABB_tree<AABBTraits>& outTree) const
  {
    build_aabb_tree(m_graph, outTree, m_vertexPointMap);
  }

  /// \cond

  template <class AABBTraits>
  void build_aabb_tree(AABB_tree<AABBTraits>& outTree,
                       Vertex_point_map vertexPointMap) const
  {
    build_aabb_tree(m_graph, outTree, vertexPointMap);
  }

  template <class AABBTraits>
  static void build_aabb_tree(const Triangle_mesh& tm,
                              AABB_tree<AABBTraits>& outTree,
                              Vertex_point_map vertexPointMap)
  {
    face_iterator facesStart, facesEnd;
    std::tie(facesStart, facesEnd) = faces(tm);
    outTree.rebuild(facesStart, facesEnd, tm, vertexPointMap);
    outTree.build();
  }
  /// \endcond

  /// @}
};

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SHORTEST_PATH_SURFACE_MESH_SHORTEST_PATH_H
