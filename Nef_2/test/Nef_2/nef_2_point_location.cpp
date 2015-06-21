
#include <CGAL/Exact_integer.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <cassert>

template<typename Kernel>
void test() {

  typedef CGAL::Nef_polyhedron_2<Kernel> Nef_polyhedron;
  typedef typename Nef_polyhedron::Point Point;
  typedef typename Nef_polyhedron::Direction Direction;
  typedef typename Nef_polyhedron::Explorer Explorer;
  typedef typename Nef_polyhedron::Object_handle    Object_handle;
  typedef typename Explorer::Vertex_const_handle    Vertex_const_handle;

  Point p11(0,0), p12(-1,1), p13(0,2);
  Point p31(0,0), p32(0,1), p33(1,2);
  Point p41(0,0), p42(1,1), p43(0,2);
  Point line1[3] = { p11, p12, p13};  
  Point line3[3] = { p31, p32, p33};  
  Point line4[3] = { p41, p42, p43};  
  std::pair<Point*, Point*> pr1(line1, line1+3);
  std::pair<Point*, Point*> pr3(line3, line3+3);
  std::pair<Point*, Point*> pr4(line4, line4+3);
  std::list<std::pair<Point*, Point*> > poly;
  poly.push_back(pr1);
  poly.push_back(pr3);
  poly.push_back(pr4);
  Nef_polyhedron N(poly.begin(), poly.end(), Nef_polyhedron::POLYLINES);
  
  //  Point q1(0,0), q2(1,0), q3(1,1), q4(0,1);
  //  Point x[4] = { q1,q2,q3,q4};
  //  Nef_polyhedron M(x,x+4, Nef_polyhedron::INCLUDED); 
  //  std::cerr << M;

  //  CGAL_NEF_SETDTHREAD(17);
  Object_handle o = N.locate(Point(2,1));
  o = N.locate(Point(2,1));
  o = N.ray_shoot(Point(4,3), Direction(-2,-1));
  o = N.ray_shoot_to_boundary(Point(2,1), Direction(-1,0));

  Vertex_const_handle v;
  Explorer E = N.explorer();
  assert(CGAL::assign(v,o));
  assert(E.point(v) == Point(1,1));
}
 
int main() {

  typedef CGAL::Exact_integer RT;
  typedef CGAL::Filtered_extended_homogeneous<RT> EKernel;
  typedef CGAL::Homogeneous<RT> HOM;
  typedef CGAL::Bounded_kernel<HOM> BKernel;

  test<EKernel>();
  test<BKernel>();

  return 0;
}


