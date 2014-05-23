

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

} // namespace internal

} // namespace CGAL