#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;

struct Box
  : public CGAL::Box_intersection_d::Box_d< double, 3, CGAL::Box_intersection_d::ID_EXPLICIT>
{
  Segment_3 segment;

  typedef CGAL::Box_intersection_d::Box_d< double, 3, CGAL::Box_intersection_d::ID_EXPLICIT> Base;
  Box(const Point_3& p1,
      const Point_3& p2)
    :Base(p1.bbox()+p2.bbox()), segment(p1,p2)
  {}
};

struct Report_inters{
  void operator() ( const Box& a, const Box& b) {
    const Segment_3& s1=a.segment;
    const Segment_3& s2=b.segment;
    if (s1.is_degenerate() || s2.is_degenerate() ) return;
    if( do_intersect(s1, s2) )
      if ( !(s1[0]==s2[0] || s1[0]==s2[1] || s1[1]==s2[0] || s1[1]==s2[1] ) )
        std::cout << "2 " << s1 << "\n2 " << s2 << "\n";
  }
};

int main( int argc, char** argv)
{
  if (argc!=2)
  {
    std::cerr << "Please provide a polyline file\n";
    return 1;
  }
  std::ifstream input (argv[1]);
  std::vector< Box> boxes;
  int n;
  while(input >> n)
  {
    Point_3 p1, p2;
    input >> p1;
    for(int k=0; k<n;++k )
    {
      input >> p2;
      boxes.push_back( Box(p1,p2) );
      p1=p2;
    }

    CGAL::box_self_intersection_d( boxes.begin(), boxes.end(), Report_inters());
  }
}