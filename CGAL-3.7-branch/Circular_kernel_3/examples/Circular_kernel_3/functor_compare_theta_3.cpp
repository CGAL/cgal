#include <CGAL/Exact_spherical_kernel_3.h>

typedef CGAL::Exact_spherical_kernel_3               SK;

int main(){
  //construction of 3 spheres from their centers and squared radii
  SK::Sphere_3 s1(SK::Point_3(0,0,0),2);
  SK::Sphere_3 s2(SK::Point_3(0,1,0),1);
  SK::Sphere_3 s3(SK::Point_3(1,0,0),3);

  //construct two circles lying on sphere s1
  SK::Circle_3 C1(s1,s2);
  SK::Circle_3 C2(s1,s3);
  
  SK::Intersect_3 inter;
  //create a functor to compare theta-coordinates on sphere s1
  SK::Compare_theta_z_3 cmp(s1);
  std::vector< CGAL::Object > intersections;
  inter(C1,C2,std::back_inserter(intersections));

  //unsigned integer indicates multiplicity of intersection point 
  std::pair<SK::Circular_arc_point_3,unsigned> p1=
    CGAL::object_cast< std::pair<SK::Circular_arc_point_3,unsigned> >(intersections[0]);
  std::pair<SK::Circular_arc_point_3,unsigned> p2=
    CGAL::object_cast< std::pair<SK::Circular_arc_point_3,unsigned> >(intersections[1]);


  SK::Circular_arc_point_3 t_extreme[2];
  //Compute theta extremal points of circle C1 on sphere s1
  CGAL::theta_extremal_points(C1,s1,t_extreme);
  
  //The theta coordinates of theta extremal points of C1 enclose that of each intersection point.
  assert(cmp(t_extreme[0],p1.first)==CGAL::SMALLER);
  assert(cmp(t_extreme[0],p2.first)==CGAL::SMALLER);
  assert(cmp(t_extreme[1],p1.first)==CGAL::LARGER);
  assert(cmp(t_extreme[1],p2.first)==CGAL::LARGER);
  
  return 0;
}
