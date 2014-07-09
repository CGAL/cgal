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
#include <CGAL/result_of.h>
#include <CGAL/Polyhedron_shortest_path/Internal/misc_functions.h>

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H

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
  typedef typename K::Segment_2 Segment_2;
  
  typedef typename K::Compute_squared_distance_3 Compute_squared_distance_3;

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

template <class Kernel>
class Parametric_distance_along_segment_2
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Segment_2 Segment_2;
  
  typedef typename Kernel::Intersect_2 Intersect_2;
  
private:
  Intersect_2 m_intersect_2;
  
public:
  FT operator () (const Point_2& x0, const Point_2& x1, const Point_2& point)
  {
    Vector_2 lineDiff = x1 - x0;
    Vector_2 pointDiff = point - x0;
    
    if (CGAL::abs(lineDiff[0]) > CGAL::abs(lineDiff[1]))
    {
      return pointDiff[0] / lineDiff[0];
    }
    else 
    {
      return pointDiff[1] / lineDiff[1];
    }
  }
};

template <class K>
class Compare_relative_intersection_along_segment_2
{
public:
  typedef typename K::FT FT;
  typedef typename K::Ray_2 Ray_2;
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::Triangle_2 Triangle_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Line_2 Line_2;

  typedef typename K::Intersect_2 Intersect_2;
 
  
private:
  Parametric_distance_along_segment_2<K> m_parametric_distance_along_segment_2;
  Intersect_2 m_intersect_2;
  
  bool in_range(FT x, FT end1, FT end2)
  {
    return x >= std::min(end1, end2) && x <= std::max(end1, end2);
  }
  
public:
  Compare_relative_intersection_along_segment_2()
  {
  }
  
  Compare_relative_intersection_along_segment_2(const Parametric_distance_along_segment_2<K>& pds, const Intersect_2& i2)
    : m_parametric_distance_along_segment_2(pds)
    , m_intersect_2(i2)
  {
  }

  CGAL::Comparison_result operator () (const Segment_2& s1, const Line_2& l1, const Segment_2& s2, const Line_2& l2)
  {
    typedef typename cpp11::result_of<Intersect_2(Line_2, Line_2)>::type LineLineIntersectResult;
  
    Line_2 s1Line(s1.source(), s1.target());
    
    Line_2 s2Line(s2.source(), s2.target());
    
    LineLineIntersectResult intersectResult1 = m_intersect_2(s1Line, l1);
    
    assert(intersectResult1);
    
    FT t1(0.0);
    
    if (intersectResult1)
    {
      Point_2* result = boost::get<Point_2>(&*intersectResult1);
      
      assert(result && "Intersection should have been a point");
      
      if (result)
      {
        t1 = m_parametric_distance_along_segment_2(s1.source(), s1.target(), *result);
        assert(t1 >= FT(-0.00001) && t1 <= FT(1.00001));
      }
    }
    
    LineLineIntersectResult intersectResult2 = m_intersect_2(s2Line, l2);
    
    assert(intersectResult2);
    
    FT t2(0.0);
    
    if (intersectResult2)
    {
      Point_2* result = boost::get<Point_2>(&*intersectResult2);
      
      assert(result && "Intersection should have been a point");
      
      if (result)
      {
        t2 = m_parametric_distance_along_segment_2(s2.source(), s2.target(), *result);
        assert(t2 >= FT(-0.00001) && t2 <= FT(1.00001));
      }
    }
    
    if (t1 == t2)
    {
      return CGAL::EQUAL;
    }
    else if (t1 < t2)
    {
      return CGAL::SMALLER;
    }
    else
    {
      return CGAL::LARGER;
    }
  }
};

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
  
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  
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
  
  bool operator() (vertex_descriptor v, Polyhedron& P)
  {
    return (*this)(v, P, CGAL::get(CGAL::vertex_point, P));
  }
  
  template<class VertexPointMap>
  bool operator() (vertex_descriptor v, Polyhedron& P, VertexPointMap const& pointMap)
  {
    halfedge_descriptor startEdge = CGAL::halfedge(v, P);
    
    Point_3 rootPoint(pointMap[v]);
    Point_3 prevPoint(pointMap[CGAL::source(startEdge, P)]);
    
    halfedge_descriptor currentEdge = CGAL::next(startEdge, P);
    
    Point_3 nextPoint(pointMap[CGAL::target(currentEdge, P)]);
    
    Triangle_3 baseFace3(rootPoint, nextPoint, prevPoint);
    
    currentEdge = CGAL::opposite(currentEdge, P);
    
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
      //std::cout << "Here:" << __LINE__ << std::endl;
      
      prevPoint = nextPoint;
      currentEdge = CGAL::next(currentEdge, P);
      nextPoint = pointMap[CGAL::target(currentEdge, P)];
      currentEdge = CGAL::opposite(currentEdge, P);
      
      Triangle_3 currentFace3(rootPoint, nextPoint, prevPoint);
      //std::cout << "Here:" << __LINE__ << std::endl;
      Triangle_2 currentFace2(m_flatten_triangle_3_along_segment_2(currentFace3, 2, nextSegment));
      //std::cout << "Here:" << __LINE__ << std::endl;

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

#endif /* CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_FUNCTION_OBJECTS_H */
