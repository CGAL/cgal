#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>

typedef leda_integer                   RT;
typedef CGAL::Homogeneous<RT>          Kernel;
typedef CGAL::Sphere_point<Kernel>     SPoint;
typedef CGAL::Sphere_segment<Kernel>   SSegment;
typedef CGAL::Sphere_circle<Kernel>    SCircle;
typedef CGAL::Sphere_direction<Kernel> SDirection;
typedef CGAL::Plane_3<Kernel>          Plane;

int main() 
{
  CGAL::set_pretty_mode ( std::cout );
  SPoint p(0,0,1), q(1,1,0), r(1,-1,0), s(1,1,1);
  SSegment s1(p,q), s2(p,r,false), s3(SPoint(0,-1,0),SPoint(-1,0,0));
  SCircle c0, c1(p,q), c2(1,1,1), c3(Plane(1,1,1,0));
  CGAL_assertion(p.opposite().opposite()==p);
  CGAL_assertion(c1.opposite().opposite()==c1);
  CGAL_assertion(c1.has_on(p)&&c1.has_on(q));
  CGAL_assertion(c3.plane()==Plane(1,1,1,0));
  c1.split_at(p);
  c1.split_at_xy_plane();
  CGAL_assertion(s1.is_short());
  CGAL_assertion(s2.is_long());

  std::list<SSegment> L,Lp;
  std::list<SSegment>::iterator it;
  L.push_back(s1);
  L.push_back(s2);
  L.push_back(s3);
  SCircle pos_xy(0,0,1);
  partition( pos_xy, L.begin(), L.end(), Lp);
  for (it = Lp.begin(); it != Lp.end(); ++it) {
    std::cout << *it << std::endl;
  }
  return 0;
}


