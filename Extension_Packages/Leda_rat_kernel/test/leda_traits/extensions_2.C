// test 2d extensions ...

#include <CGAL/basic.h>
#include <CGAL/convex_hull_2.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <LEDA/random_rat_point.h>
#include <iostream>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;
typedef K::Do_intersect_to_right_2                             Do_intersect_to_right_2;
typedef K::Do_intersect_to_left_2                              Do_intersect_to_left_2;


bool cgal_do_intersect_to_right(const leda_rat_segment& c1, const leda_rat_segment& c2,
                                const leda_rat_point& pt) 
{
    leda_rat_segment seg;
    bool res = c1.intersection(c2, seg);
    
    if (! res) return false;
    
    if (seg.start() == seg.end()) { // point ...
      return (leda_rat_point::cmp_xy(seg.start(), pt) == 1);
    }
    
    // intersection result is a segment ...
    return ( (leda_rat_point::cmp_xy(seg.start(), pt) == 1) || (leda_rat_point::cmp_xy(seg.end(), pt) == 1) );
}


bool cgal_do_intersect_to_left(const leda_rat_segment& c1, const leda_rat_segment& c2,
                               const leda_rat_point& pt) 
{
    leda_rat_segment seg;
    bool res = c1.intersection(c2, seg);
    
    if (! res) return false;
    
    if (seg.start() == seg.end()) { // point ...
      return (leda_rat_point::cmp_xy(seg.start(), pt) == -1);
    }
    
    // intersection result is a segment ...
    return ( (leda_rat_point::cmp_xy(seg.start(), pt) == -1) || (leda_rat_point::cmp_xy(seg.end(), pt) == -1) );
}

void my_random_point_in_square(rat_point& p, int maxc)
{
  random_point_in_square(p, maxc);
  
  // generate w value ...
  leda_random_source S(1,10);
  
  int w;
  S >> w;  
  p = leda_rat_point(p.X()*w,p.Y()*w,p.W()*w);
}


int main()
{
  int number_of_tests = 100000, i;
  
  std::cout << "test intersection to right ...\n";
  
  Do_intersect_to_right_2   do_intersect_to_right; 
  Do_intersect_to_left_2    do_intersect_to_left;
  
  for(i=0; i<number_of_tests;i++){
      rat_point p1,p2,p3,p4, pt;
  
      my_random_point_in_square(p1, 10000);
      my_random_point_in_square(p2, 10000);
      my_random_point_in_square(p3, 10000);
      my_random_point_in_square(p4, 10000);
      my_random_point_in_square(pt, 10000);
      
      leda_rat_segment s1(p1,p2), s2(p3,p4);
      
      bool res1 = cgal_do_intersect_to_right(s1,s2,pt);
      bool res2 = do_intersect_to_right(s1,s2,pt);
      
      if (res1 != res2){
         std::cout << "error !\n";
	 std::cout << "res1/res2:" << res1 << " " << res2 << "\n";
	 std::cout << "arguments: s1-  " << s1 << "  s2-  " << s2 << "  pt-  " << pt << "\n";
	 return 1;
      }
  }
  
  std::cout << "test intersection to left ...\n";  
  
  for(i=0; i<number_of_tests;i++){
      rat_point p1,p2,p3,p4, pt;
  
      my_random_point_in_square(p1, 10000);
      my_random_point_in_square(p2, 10000);
      my_random_point_in_square(p3, 10000);
      my_random_point_in_square(p4, 10000);
      my_random_point_in_square(pt, 10000);
      
      leda_rat_segment s1(p1,p2), s2(p3,p4);
      
      bool res1 = cgal_do_intersect_to_left(s1,s2,pt);
      bool res2 = do_intersect_to_left(s1,s2,pt);
      
      if (res1 != res2){
         std::cout << "error !\n";
	 std::cout << "res1/res2:" << res1 << " " << res2 << "\n";
	 std::cout << "arguments: s1-  " << s1 << "  s2-  " << s2 << "  pt-  " << pt << "\n";
	 return 1;
      }
  }  
  return 0;
}
