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
  
  typedef typename K::Compute_squared_distance_2 Compute_squared_distance_2;
  typedef typename K::Intersect_2 Intersect_2;
 
  
private:
  Compute_squared_distance_2 m_compute_squared_distance_2;
  Intersect_2 m_intersect_2;
  
public:
  Compare_relative_intersection_along_segment_2()
  {
  }
  
  Compare_relative_intersection_along_segment_2(const Compute_squared_distance_2& cds, const Intersect_2& i2)
    : m_compute_squared_distance_2(cds)
    , m_intersect_2(i2)
  {
  }

  CGAL::Comparison_result operator () (const Segment_2& s1, const Ray_2& r1, const Segment_2& s2, const Ray_2& r2)
  {
    typedef typename cpp11::result_of<Intersect_2(Segment_2, Line_2)>::type SegmentLineIntersectResult;

    SegmentLineIntersectResult s1r1Intersection = m_intersect_2(s1, r1.supporting_line());
    Point_2 p1;

    if (s1r1Intersection)
    {
      Point_2* result = boost::get<Point_2>(&*s1r1Intersection);
      
      if (result)
      {
        p1 = *result;
       
      }
      else
      {
        assert(false && "Ray entering triangle must not be parallel to entry segment.");
      }
    }
    else
    {
      // TODO: figure out what is causing this, i.e. is it just out of range, or is the algorithm incorrect
      std::cout << "Segment = " << s1 << std::endl;
      std::cout << "Ray = " << r1 << std::endl;
      Point_2 projection = s1.supporting_line().projection(r1.source());
      std::cout << "Proj = " << projection << std::endl;
      std::cout << "Dist_0^2 = " << m_compute_squared_distance_2(projection, s1[0]) << std::endl;
      std::cout << "Dist_1^2 = " << m_compute_squared_distance_2(projection, s1[1]) << std::endl;
      
      typedef typename K::Line_2 Line_2;
      typedef typename cpp11::result_of<Intersect_2(Line_2, Ray_2)>::type LineRayIntersectResult;

      LineRayIntersectResult lri = m_intersect_2(s1.supporting_line(), r1);
      
      if (lri)
      {
        Point_2* result = boost::get<Point_2>(&*lri);
      
        if (result)
        {
          std::cout << "Line Intersection = " << *result << std::endl;
         
        }
      }
      
      assert(s1r1Intersection && "Ray must enter triangle via entry segment.");
    }
    
    SegmentLineIntersectResult s2r2Intersection = m_intersect_2(s2, r2.supporting_line());
    Point_2 p2;

    if (s2r2Intersection)
    {
      Point_2* result = boost::get<Point_2>(&*s2r2Intersection);
      
      if (result)
      {
        p2 = *result;
       
      }
      else
      {
        assert(false && "Ray entering triangle must not be parallel to entry segment.");
      }
    }
    else
    {
      // TODO: same as above
      std::cout << "Segment = " << s2 << std::endl;
      std::cout << "Ray = " << r2 << std::endl;
      Point_2 projection = s2.supporting_line().projection(r1.source());
      std::cout << "Proj = " << projection << std::endl;
      std::cout << "Dist_0^2 = " << m_compute_squared_distance_2(projection, s2[0]) << std::endl;
      std::cout << "Dist_1^2 = " << m_compute_squared_distance_2(projection, s2[1]) << std::endl;
      
      typedef typename K::Line_2 Line_2;
      typedef typename cpp11::result_of<Intersect_2(Line_2, Ray_2)>::type LineRayIntersectResult;

      LineRayIntersectResult lri = m_intersect_2(s2.supporting_line(), r2);
      
      if (lri)
      {
        Point_2* result = boost::get<Point_2>(&*lri);
      
        if (result)
        {
          std::cout << "Line Intersection = " << *result << std::endl;
         
        }
      }
      
      assert(s2r2Intersection && "Ray must enter triangle via entry segment.");
    }
    
    FT d1 = m_compute_squared_distance_2(s1[0], p1);
    FT d2 = m_compute_squared_distance_2(s2[0], p2);
    
    if (d1 == d2)
    {
      return CGAL::EQUAL;
    }
    else if (d1 < d2)
    {
      return CGAL::SMALLER;
    }
    else if (d1 > d2)
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
    
    halfedge_descriptor currentEdge = currentEdge = CGAL::next(startEdge, P);
    
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
      prevPoint = nextPoint;
      currentEdge = CGAL::next(currentEdge, P);
      nextPoint = pointMap[CGAL::target(currentEdge, P)];
      currentEdge = CGAL::opposite(currentEdge, P);
      
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