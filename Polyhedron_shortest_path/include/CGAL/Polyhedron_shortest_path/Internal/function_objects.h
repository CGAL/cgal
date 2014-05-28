// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
//#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
//#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

namespace internal {

template<class K>
class Project_triangle_3_to_triangle_2
{
public:
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;

private:
  Compute_squared_distance_3 m_compute_squared_distance_3;
  
public:
  Project_triangle_3_to_triangle_2()
  {
  }
  
  Project_triangle_3_to_triangle_2(const Compute_squared_distance_3& cds)
    : m_compute_squared_distance_3(cds)
  {
  }

  Triangle_2 operator() (const Triangle_3& t3) const
  {
    Vector_3 v01 = t3[1] - t3[0];
    Vector_3 v02 = t3[2] - t3[0];
    
    FT scalePoint = (v01 * v02) / (v01 * v01);
    Point_3 projectedLocation3d = t3[0] + (scalePoint * v01);
    FT triangleHeight = CGAL::sqrt(m_compute_squared_distance_3(projectedLocation3d, t3[2]));
    FT v01Len = CGAL::sqrt(m_compute_squared_distance_3(t3[1], t3[0])); 
    
    Point_2 A(0.0, 0.0);
    Point_2 B(v01Len, 0.0);
    Point_2 C(v01Len * scalePoint, triangleHeight);
    
    return Triangle_2(A, B, C);
  }
};

template<class K>
class Flatten_triangle_3_along_segment_2
{
public:
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::FT FT;
  typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;
  typedef typename K::Segment_2 Segment_2;
  
private:
  Compute_squared_distance_3 m_compute_squared_distance_3;

public:
  Flatten_triangle_3_along_segment_2()
  {
  }
  
  Flatten_triangle_3_along_segment_2(const Compute_squared_distance_3& cds)
    : m_compute_squared_distance_3(cds)
  {
  }

  Triangle_2 operator() (const Triangle_3& t3, size_t edgeIndex, const Segment_2& segment) const
  {
    Vector_3 v01 = t3.vertex(edgeIndex + 1) - t3.vertex(edgeIndex + 0);
    Vector_3 v02 = t3.vertex(edgeIndex + 2) - t3.vertex(edgeIndex + 0);
    
    FT scalePoint = (v01 * v02) / (v01 * v01);
    Point_3 projectedLocation3d = t3.vertex(edgeIndex) + (scalePoint * v01);
    FT triangleHeight = CGAL::sqrt(m_compute_squared_distance_3(projectedLocation3d, t3.vertex(edgeIndex + 2)));

    Vector_2 edgeVector = segment.to_vector();
    Point_2 projectionPoint = segment.start() + (segment.to_vector() * scalePoint);
    
    Vector_2 perpendicularEdgeVector(-edgeVector[1], edgeVector[0]);
    perpendicularEdgeVector = perpendicularEdgeVector / CGAL::sqrt(perpendicularEdgeVector.squared_length());
    
    Point_2 points[3];
    points[edgeIndex] = segment.start();
    points[(edgeIndex + 1) % 3] = segment.end();
    points[(edgeIndex + 2) % 3] = segment.start() + (edgeVector * scalePoint) + (perpendicularEdgeVector * triangleHeight);
    return Triangle_2(points[0], points[1], points[2]);
  }
};

/*
// TODO: test these in isoloation at least
// This is BGL based, but I'm not really convinced that makes a lot of sense...
template <class P>
class Is_vertex_convex_BGL
{
public:
  typedef P Polyhedron;
  
  typedef typename Polyhedron::FT FT;
  typedef typename Polyhedron::Point_3 Point_3;
  typedef typename Polyhedron::Vector_3 Vector_3;
  
  // typedef typename boost::graph_traits<Polyhedron>::face_descriptor Face_descriptor;
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor Vertex_descriptor;
  typedef typename GraphTraits::vertex_iterator Vertex_iterator;
  typedef typename GraphTraits::adjacency_iterator Adjacency_iterator;
  typedef typename GraphTraits::vertex_descriptor Edge_descriptor;
  
  bool operator() (Vertex_descriptor v, Polyhedron& polyhedron)
  {
    Adjacency_iterator begin, end;
    boost::tie(begin, end) = boost::adjacent_vertices(v, polyhedron);
    
    Adjacency_iterator beforeEnd = end;
    --beforeEnd;
    
    Vector_3 previousEdge = (*beforeEnd)->point() - v->point();
    Vector_3 currentEdge = (*begin)->point() - v->point();
    
    for (Adjacency_iterator current = begin; current != end; ++current)
    {
      Adjacency_iterator next = current;
      ++next;
      
      if (next == end)
      {
        next = begin;
      }
      
      Vector_3 nextEdge = (*next)->point() - v->point();
      
      Vector_3 currentPlane = CGAL::cross_product(previousEdge, currentEdge);
      
      if (CGAL::is_positive(currentPlane * nextEdge))
      {
        return false;
      }
      
      previousEdge = currentEdge;
      nextEdge = currentEdge;
    }
    
    return true;
  }
};
*/

template <class Kernel, class Polyhedron>
class Is_saddle_vertex
{
public:
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_2 Point_2;
  
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  
  typedef typename CGAL::internal::Project_triangle_3_to_triangle_2<Kernel> Project_triangle_3_to_triangle_2;
  typedef typename CGAL::internal::Flatten_triangle_3_along_segment_2<Kernel> Flatten_triangle_3_along_segment_2;
  typedef typename Kernel::Orientation_2 Orientation_2;
  
private:
  Project_triangle_3_to_triangle_2 m_project_triangle_3_to_triangle_2;
  Flatten_triangle_3_along_segment_2 m_flatten_triangle_3_along_segment_2;
  Orientation_2 m_orientation_2;
  
public:

  Is_saddle_vertex()
  {
  }
  
  Is_saddle_vertex(Project_triangle_3_to_triangle_2 pt3tt2, Flatten_triangle_3_along_segment_2 ft3as2, Orientation_2 o2)
    : m_project_triangle_3_to_triangle_2(pt3tt2)
    , m_flatten_triangle_3_along_segment_2(ft3as2)
    , m_orientation_2(o2)
  {
  }
  
  bool operator() (Vertex_handle v)
  {
    Halfedge_handle startEdge = v->halfedge();
    
    Halfedge_handle currentEdge = startEdge;
    
    Point_3 rootPoint(v->point());
    Point_3 nextPoint(currentEdge->next()->vertex()->point());
    Point_3 prevPoint(currentEdge->prev()->vertex()->point());
    Triangle_3 baseFace3(rootPoint, nextPoint, prevPoint);
    
    currentEdge = currentEdge->next()->opposite();
    
    Triangle_2 baseFace2(m_project_triangle_3_to_triangle_2(baseFace3));
    
    Segment_2 baseSegment(baseFace2[0], baseFace2[2]);

    Segment_2 nextSegment(baseFace2[1], baseFace2[0]);
    
    CGAL::Orientation baseOrientation = m_orientation_2(baseFace2[0], baseFace2[2], baseFace2[1]);
    
    if (baseOrientation == CGAL::COLLINEAR)
    {
      // I would say this violates a precondition
    }
    
    do
    {
      prevPoint = nextPoint;
      currentEdge = currentEdge->next();
      nextPoint = currentEdge->vertex()->point();
      currentEdge = currentEdge->opposite();
      
      Triangle_3 currentFace3(rootPoint, nextPoint, prevPoint);
      Triangle_2 currentFace2(m_flatten_triangle_3_along_segment_2(currentFace3, 2, nextSegment));

      if (m_orientation_2(baseSegment[0], baseSegment[1], currentFace2[2]) != baseOrientation && m_orientation_2(baseSegment[0], baseSegment[1], currentFace2[1]) == baseOrientation)
      {
        return true;
      }
      
      nextSegment = Segment_2(currentFace2[1], currentFace2[0]);
    }
    while (currentEdge != startEdge);
    
    return false;
  }
};

} // namespace internal

} // namespace CGAL