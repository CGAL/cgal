
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//#include <CGAL/Cartesian.h>
//typedef CGAL::Cartesian<double> K;


#if (SEL == 1)
#include <CGAL/Intersections_3/Bbox_3_Bbox_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Bbox_3 T2;
#elif (SEL == 2)
#include <CGAL/Intersections_3/Bbox_3_Iso_cuboid_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Iso_cuboid_3<K> T2;
#elif (SEL == 3)
#include <CGAL/Intersections_3/Bbox_3_Line_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Line_3<K> T2;
#elif (SEL == 4)
#include <CGAL/Intersections_3/Bbox_3_Plane_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Plane_3<K> T2;
#elif (SEL == 5)
#include <CGAL/Intersections_3/Bbox_3_Ray_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 6)
#include <CGAL/Intersections_3/Bbox_3_Segment_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 7)
#include <CGAL/Intersections_3/Bbox_3_Sphere_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 8)
#include <CGAL/Intersections_3/Bbox_3_Tetrahedron_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 9)
#include <CGAL/Intersections_3/Bbox_3_Triangle_3.h>
typedef CGAL::Bbox_3 T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 10)
#include <CGAL/Intersections_3/Iso_cuboid_3_Iso_cuboid_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Iso_cuboid_3<K> T2;
#elif (SEL == 11)
#include <CGAL/Intersections_3/Iso_cuboid_3_Line_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Line_3<K> T2;
#elif (SEL == 12)
#include <CGAL/Intersections_3/Iso_cuboid_3_Plane_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Plane_3<K> T2;
#elif (SEL == 13)
#include <CGAL/Intersections_3/Iso_cuboid_3_Point_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Point_3<K> T2;
#elif (SEL == 14)
#include <CGAL/Intersections_3/Iso_cuboid_3_Ray_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 15)
#include <CGAL/Intersections_3/Iso_cuboid_3_Segment_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 16)
#include <CGAL/Intersections_3/Iso_cuboid_3_Sphere_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 17)
#include <CGAL/Intersections_3/Iso_cuboid_3_Tetrahedron_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 18)
#include <CGAL/Intersections_3/Iso_cuboid_3_Triangle_3.h>
typedef CGAL::Iso_cuboid_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 19)
#include <CGAL/Intersections_3/Line_3_Line_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Line_3<K> T2;
#elif (SEL == 20)
#include <CGAL/Intersections_3/Line_3_Plane_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Plane_3<K> T2;
#elif (SEL == 21)
#include <CGAL/Intersections_3/Line_3_Point_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Point_3<K> T2;
#elif (SEL == 22)
#include <CGAL/Intersections_3/Line_3_Ray_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 23)
#include <CGAL/Intersections_3/Line_3_Segment_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 24)
#include <CGAL/Intersections_3/Line_3_Sphere_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 25)
#include <CGAL/Intersections_3/Line_3_Tetrahedron_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 26)
#include <CGAL/Intersections_3/Line_3_Triangle_3.h>
typedef CGAL::Line_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 27)
#include <CGAL/Intersections_3/Plane_3_Plane_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Plane_3<K> T2;
#elif (SEL == 28)
#include <CGAL/Intersections_3/Plane_3_Point_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Point_3<K> T2;
#elif (SEL == 29)
#include <CGAL/Intersections_3/Plane_3_Ray_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 30)
#include <CGAL/Intersections_3/Plane_3_Segment_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 31)
#include <CGAL/Intersections_3/Plane_3_Sphere_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 32)
#include <CGAL/Intersections_3/Plane_3_Tetrahedron_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 33)
#include <CGAL/Intersections_3/Plane_3_Triangle_3.h>
typedef CGAL::Plane_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 34)
// duplicate
#elif (SEL == 35)
#include <CGAL/Intersections_3/Point_3_Point_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Point_3<K> T2;
#elif (SEL == 36)
#include <CGAL/Intersections_3/Point_3_Ray_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 37)
#include <CGAL/Intersections_3/Point_3_Segment_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 38)
#include <CGAL/Intersections_3/Point_3_Sphere_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 39)
#include <CGAL/Intersections_3/Point_3_Tetrahedron_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 40)
#include <CGAL/Intersections_3/Point_3_Triangle_3.h>
typedef CGAL::Point_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 41)
#include <CGAL/Intersections_3/Ray_3_Ray_3.h>
typedef CGAL::Ray_3<K> T1;
typedef CGAL::Ray_3<K> T2;
#elif (SEL == 42)
#include <CGAL/Intersections_3/Ray_3_Segment_3.h>
typedef CGAL::Ray_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 43)
#include <CGAL/Intersections_3/Ray_3_Sphere_3.h>
typedef CGAL::Ray_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 44)
#include <CGAL/Intersections_3/Ray_3_Tetrahedron_3.h>
typedef CGAL::Ray_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 45)
#include <CGAL/Intersections_3/Ray_3_Triangle_3.h>
typedef CGAL::Ray_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 46)
#include <CGAL/Intersections_3/Segment_3_Segment_3.h>
typedef CGAL::Segment_3<K> T1;
typedef CGAL::Segment_3<K> T2;
#elif (SEL == 47)
#include <CGAL/Intersections_3/Segment_3_Sphere_3.h>
typedef CGAL::Segment_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 48)
#include <CGAL/Intersections_3/Segment_3_Tetrahedron_3.h>
typedef CGAL::Segment_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 49)
#include <CGAL/Intersections_3/Segment_3_Triangle_3.h>
typedef CGAL::Segment_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 50)
#include <CGAL/Intersections_3/Sphere_3_Sphere_3.h>
typedef CGAL::Sphere_3<K> T1;
typedef CGAL::Sphere_3<K> T2;
#elif (SEL == 51)
#include <CGAL/Intersections_3/Sphere_3_Tetrahedron_3.h>
typedef CGAL::Sphere_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 52)
#include <CGAL/Intersections_3/Sphere_3_Triangle_3.h>
typedef CGAL::Sphere_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#elif (SEL == 53)
#include <CGAL/Intersections_3/Tetrahedron_3_Tetrahedron_3.h>
typedef CGAL::Tetrahedron_3<K> T1;
typedef CGAL::Tetrahedron_3<K> T2;
#elif (SEL == 54)
#include <CGAL/Intersections_3/Tetrahedron_3_Triangle_3.h>
typedef CGAL::Tetrahedron_3<K> T1;
typedef CGAL::Triangle_3<K> T2;

#else
#include <CGAL/Intersections_3/Triangle_3_Triangle_3.h>
typedef CGAL::Triangle_3<K> T1;
typedef CGAL::Triangle_3<K> T2;
#endif

int main()
{
  T1 t1;
  T2 t2;
  try {
    do_intersect(t1,t2);
    do_intersect(t2,t1);
    intersection(t1,t2);
    intersection(t2,t1);
  } catch(...){}
  return 0;
}
 
