#define CGAL_INTERSECT_WITH_ITERATORS_2 1
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_kernel.h>
#include <list>
#include "Cartesian_I.h"




typedef CGAL::Cartesian_I<CGAL::Gmpq> SC;
typedef CGAL::Lazy_kernel<SC, CGAL::Cartesian_I<CGAL::Interval_nt_advanced > > K;



typedef K::FT FT;

typedef K::Point_2 Point_2;

typedef K::Vector_2 Vector_2;
typedef K::Segment_2 Segment_2;
typedef K::Circle_2 Circle_2;

typedef CGAL::Bbox_2 Bbox_2;
typedef CGAL::Object Object;





int main()
{
  CGAL::Lazy_exact_nt<CGAL::Gmpq> nt;
  nt = nt + nt * nt;
  
  K::Intersect_with_iterators_2 iwi;

  CGAL::set_pretty_mode(std::cout);
  K::Intersect_with_iterators_2 intersect;

  Segment_2 s1(Point_2(0,1), Point_2(2,1));
  Segment_2 s2(Point_2(1,0), Point_2(1,2));

  Point_2 ip;
  std::list<Object> intersections;
  intersect(s1, s2, std::back_inserter(intersections));
  for(std::list<CGAL::Object>::iterator it = intersections.begin(); it != intersections.end(); it++){
    if(CGAL::assign(ip, *it)){
      std::cout << "intersection at " << ip << std::endl;
    }
  }


  

  FT ft = 3.1415;
  std::cout << "ft = " << ft << std::endl;
  
  ft *= ft;
  std::cout << "ft^2 = " << ft << std::endl;
  std::cout << "ft^2.depth() = " << ft.depth() << std::endl;
  std::cout << "ft^2.exact() = " << ft.exact() << std::endl;
  Point_2 p(ft, 2.22);
  Point_2 q(9,9);
  
  CGAL::Bbox_2 bb = p.bbox();
  
  Segment_2 s(p,q);
  
  Segment_2 s3(Point_2(0,1), Point_2(2,1));
  Segment_2 s4(Point_2(1,0), Point_2(1,2));
  
  CGAL::Object o = intersection(s3,s4);
  Point_2 rp;

  if(CGAL::assign(rp, o)){
    std::cout << "Intersection is a point:" << std::endl;
    std::cout << rp;
  }
    

  Point_2 r = K::Construct_vertex_2()(s,0);
  assert(r == s.source());
  std::cout << r << std::endl;


  Point_2 mp = midpoint(p,q);


  FT rx = r.x();

  std::cout << rx << std::endl;
  

    Vector_2 v1(1,1), v2(1,1);

    v1 = p - q;

    v1 = mp - CGAL::ORIGIN;

    q = CGAL::ORIGIN + v1;

    std::cout << q << std::endl;

    if(v1 == v2){}
  
    if(K::Compare_distance_2()(p,q,r)== CGAL::SMALLER)
      {
	std::cout << "smaller" << std::endl;
      }

  Circle_2 circ(CGAL::ORIGIN, p, q);
  std::cout << "\nCircle:\n " << circ << std::endl;

  Point_2 center = circ.center();
  FT sr = circ.squared_radius();
  std::cout << "\nCenter = " << center << "\nSquared radius = " << sr << std::endl;

  sr += 7812;
  std::cout << "squared radius + 7812 = " << sr << std::endl;


  circ.exact();
  std::cout << "\nCircle after circ.exaxt():\n " << circ << std::endl;
  std::cout << "\nCenter = " << center << "\nSquared radius = " << sr << std::endl;
  
  std::cout << "Done"  << std::endl;

  return 0;
}










