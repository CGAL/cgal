#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_S2/Sphere_point.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/test_macros.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float NT;
#endif

typedef CGAL::Homogeneous<NT>          Kernel;
typedef CGAL::Sphere_point<Kernel>     SPoint;
typedef CGAL::Sphere_segment<Kernel>   SSegment;
typedef CGAL::Sphere_circle<Kernel>    SCircle;
typedef CGAL::Sphere_direction<Kernel> SDirection;
typedef CGAL::Plane_3<Kernel>          Plane;

int main() 
{
  CGAL_TEST_START;    
  CGAL::set_pretty_mode ( std::cout );
  SPoint p(0,0,1), q(1,1,0), r(1,-1,0), s(1,1,1);
  SSegment s1(p,q), s2(p,r,false), s3(SPoint(0,-1,0),SPoint(-1,0,0));
  SCircle c0, c1(p,q), c2(1,1,1), c3(Plane(1,1,1,0));
  CGAL_TEST(p.x() == NT(0)){}
  CGAL_TEST(p.y() == NT(0)){}
  CGAL_TEST(p.z() == NT(1)){}
  CGAL_TEST(p.antipode().antipode()==p){}
  CGAL_TEST(p.antipode()!=p){}
  
  CGAL_TEST(c1.opposite().opposite()==c1){}
  CGAL_TEST(c1.has_on(p)&&c1.has_on(q)){}
  CGAL_TEST(c3.plane()==Plane(1,1,1,0)){}
  c1.split_at(p);
  c1.split_at_xy_plane();

  CGAL_TEST(s1.is_short()){}
  CGAL_TEST(s2.is_long()){}
  CGAL_TEST(s1.source()==p){}
  CGAL_TEST(s1.target()==q){}
  CGAL_TEST(s1.sphere_circle()==c1){}
  CGAL_TEST(s1.opposite().opposite()==s1){}
  CGAL_TEST(s1.complement().complement()==s1){}
  CGAL_TEST(SSegment(p,p,c1).is_degenerate()){}
  CGAL_TEST(SSegment(p,p.antipode(),c1).is_halfcircle()){}
  CGAL_TEST(s1.has_on(p)){}
  CGAL_TEST(SSegment(p,p.antipode(),c1).has_in_relative_interior(q)){}

  std::list<SSegment> L,Lp;
  std::list<SSegment>::iterator it;
  L.push_back(s1);
  L.push_back(s2);
  L.push_back(s3);
  SCircle pos_xy(0,0,1);
  CGAL::partition( pos_xy, L.begin(), L.end(), Lp);
  for (it = Lp.begin(); it != Lp.end(); ++it) {
    std::cout << *it << std::endl;
  }
  CGAL_TEST_END;     
}


