#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;

//A class to recover Voronoi diagram from stream.
//Rays, lines and segments are cropped to a rectangle
//so that only segments are stored
struct Cropped_voronoi_from_delaunay{
  std::list<Segment_2> m_cropped_vd;
  Iso_rectangle_2 m_bbox;
  
  Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}
  
  template <class RSL>
  void crop_and_extract_segment(const RSL& rsl){
    CGAL::Object obj = CGAL::intersection(rsl,m_bbox);
    const Segment_2* s=CGAL::object_cast<Segment_2>(&obj);
    if (s) m_cropped_vd.push_back(*s);
  }
  
  void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
  void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
  void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
};

int main(){
  //consider some points
  std::vector<Point_2> points;
  points.push_back(Point_2(0,0));
  points.push_back(Point_2(1,1));
  points.push_back(Point_2(0,1));
  
  Delaunay_triangulation_2 dt2;
  //insert points into the triangulation
  dt2.insert(points.begin(),points.end());
  //construct a rectangle
  Iso_rectangle_2 bbox(-1,-1,2,2);
  Cropped_voronoi_from_delaunay vor(bbox);
  //extract the cropped Voronoi diagram
  dt2.draw_dual(vor);
  //print the cropped Voronoi diagram as segments
  std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(),
    std::ostream_iterator<Segment_2>(std::cout,"\n"));
}

