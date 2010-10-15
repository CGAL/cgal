#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <cassert>

typedef CGAL::Quotient< CGAL::MP_Float >                    FT;
typedef CGAL::Cartesian<FT>                                 K;

// The has_on test has to be working
// The construction test has to be working

template <class K>
void _test_intersection_construct(K k)
{
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Plane_3                          Plane_3;
  typedef typename K::Sphere_3                         Sphere_3;
  typedef typename K::Circle_3                         Circle_3;
  typedef typename K::Has_on_3                         Has_on_3;
  typedef typename K::Intersect_3                      Intersect_3;
  typedef typename K::Construct_sphere_3               Construct_sphere_3;
  typedef typename K::Construct_plane_3                Construct_plane_3;

  Has_on_3 theHas_on_3 = k.has_on_3_object();
  Intersect_3 theIntersect_3 = k.intersect_3_object();
  Construct_sphere_3 theConstruct_sphere_3 = k.construct_sphere_3_object();
  Construct_plane_3 theConstruct_plane_3 = k.construct_plane_3_object();

  std::cout << "Testing intersection(Sphere,Sphere)..." << std::endl;
  Sphere_3 s = theConstruct_sphere_3(Point_3(0,0,0), 1);
  Sphere_3 s_t10 = theConstruct_sphere_3(Point_3(10,10,10), 1);
  for(int vx=-3;vx<4;vx++) {
    for(int vy=-3;vy<4;vy++) {
      for(int vz=-3;vz<4;vz++) {
        for(int vr=1;vr<6;vr++) {
          const FT x = FT(vx);
          const FT y = FT(vy);
          const FT z = FT(vz);
          const FT r = FT(vr)/FT(2);
          Sphere_3 sl_1 = theConstruct_sphere_3(Point_3(x,y,z), r*r);
          Sphere_3 sl_2 = theConstruct_sphere_3(Point_3(x+10,y+10,z+10),r*r);
          int d2 = (vx*vx + vy*vy + vz*vz);

          CGAL::Object obj1 = theIntersect_3(s, sl_1);
          CGAL::Object obj2 = theIntersect_3(s_t10, sl_2);
          CGAL::Object obj3 = CGAL::intersection(s, sl_1);
          CGAL::Object obj4 = CGAL::intersection(s_t10, sl_2);

          // No intersection
          if((d2 > (r+1)*(r+1)) || (d2 < (r-1)*(r-1))) {
            assert(!CGAL::do_intersect(s, sl_1));
            assert(!CGAL::do_intersect(s_t10, sl_2));
            assert(obj1.is_empty());
            assert(obj2.is_empty());
            assert(obj3.is_empty());
            assert(obj4.is_empty());
          }

          // All the sphere intersect
          else if(x == 0 && y == 0 && z == 0 && r == 1) {
            assert(CGAL::do_intersect(s, sl_1));
            assert(CGAL::do_intersect(s_t10, sl_2));
            Sphere_3 sphere1, sphere2, sphere3, sphere4;
            assert(assign(sphere1, obj1));
            assert(assign(sphere2, obj2));
            assert(assign(sphere3, obj3));
            assert(assign(sphere4, obj4));
          } 

          // Tangent, 1 Intersection
          else if((d2 == (r+1)*(r+1)) || (d2 == (r-1)*(r-1))) {
            assert(CGAL::do_intersect(s, sl_1));
            assert(CGAL::do_intersect(s_t10, sl_2));
            Point_3 interp1, interp2, interp3, interp4;
            assert(assign(interp1, obj1));
            assert(assign(interp2, obj2));
            assert(theHas_on_3(s, interp1));
            assert(theHas_on_3(sl_1, interp1));
            assert(theHas_on_3(s_t10, interp2));
            assert(theHas_on_3(sl_2, interp2));
            assert(assign(interp3, obj3));
            assert(assign(interp4, obj4));
            assert(theHas_on_3(s, interp3));
            assert(theHas_on_3(sl_1, interp3));
            assert(theHas_on_3(s_t10, interp4));
            assert(theHas_on_3(sl_2, interp4));
          }
          // 1 Intersection Circle
          else {
            assert(CGAL::do_intersect(s, sl_1));
            assert(CGAL::do_intersect(s_t10, sl_2));
            Circle_3 circle1, circle2, circle3, circle4;
            assert(assign(circle1, obj1));
            assert(assign(circle2, obj2));
            assert(assign(circle3, obj3));
            assert(assign(circle4, obj4));
            assert(theHas_on_3(s, circle1));
            assert(theHas_on_3(sl_1, circle1));
            assert(theHas_on_3(s_t10, circle2));
            assert(theHas_on_3(sl_2, circle2));
            assert(theHas_on_3(s, circle3));
            assert(theHas_on_3(sl_1, circle3));
            assert(theHas_on_3(s_t10, circle4));
            assert(theHas_on_3(sl_2, circle4));
          }
        }
      }
    }
  }

  std::cout << "Testing intersection(Sphere,Plane)..." << std::endl;
  Plane_3 p = theConstruct_plane_3(1,1,1,-1);
  for(int vx=-2;vx<3;vx++) {
    for(int vy=-2;vy<3;vy++) {
      for(int vz=-2;vz<3;vz++) {
        for(int vr=1;vr<27;vr++) {
          const FT x = FT(vx);
          const FT y = FT(vy);
          const FT z = FT(vz);
          const FT sq_r = FT(vr)/FT(3);
          Sphere_3 sl = theConstruct_sphere_3(Point_3(x,y,z), sq_r);
          const FT d2 = ((x+y+z-1)*(x+y+z-1)/3);
          CGAL::Object obj = theIntersect_3(p, sl);
          CGAL::Object objl = CGAL::intersection(p, sl);

          // No intersection
          if(d2 > sq_r) {
            assert(!CGAL::do_intersect(p, sl));
            assert(obj.is_empty());
            assert(objl.is_empty());
          } 
          // Tangent, 1 Intersection
          else if(d2 == sq_r) {
            assert(CGAL::do_intersect(p, sl));
            Point_3 interp, interpl;
            assert(assign(interp, obj));
            assert(theHas_on_3(sl, interp));
            assert(theHas_on_3(p, interp));
            assert(assign(interpl, objl));
            assert(theHas_on_3(sl, interpl));
            assert(theHas_on_3(p, interpl));
          }
          // 1 Intersection Circle
          else {
            assert(CGAL::do_intersect(p, sl));
            Circle_3 circle1, circle2;
            assert(assign(circle1, obj));
            assert(theHas_on_3(sl, circle1));
            assert(theHas_on_3(p, circle1));
            assert(assign(circle2, objl));
            assert(theHas_on_3(sl, circle2));
            assert(theHas_on_3(p, circle2));
          }
        }
      }
    }
  }
}

int main()
{
  K k;
  _test_intersection_construct(k);
  return 0;
}
