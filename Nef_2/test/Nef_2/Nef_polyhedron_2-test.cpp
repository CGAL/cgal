#include <CGAL/basic.h>
#include <CGAL/test_macros.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Exact_integer.h>

#if defined (CGAL_USE_LEDA) || defined (CGAL_USE_GMP)

typedef CGAL::Exact_integer Integer;
typedef CGAL::Exact_rational Rational;
#else
typedef long Integer;
typedef double Rational;
#endif



int main()
{
#ifdef CGAL_CFG_ISTREAM_INT_BUG
  std::locale::global(std::locale("C")); 
#endif

  CGAL_NEF_SETDTHREAD(911); // 911
  CGAL::set_pretty_mode ( std::cerr );
  std::cerr << "using " << CGAL::pointlocationversion << std::endl;
  std::cerr << "using " << CGAL::sweepversion << std::endl;
  CGAL_TEST_START;
  
{
  typedef  CGAL::Extended_homogeneous<Integer> EKernel;
  typedef  CGAL::Nef_polyhedron_2<EKernel> Nef_polyhedron;
  typedef  Nef_polyhedron::Point     Point;
  typedef  Nef_polyhedron::Direction Direction;
  typedef  Nef_polyhedron::Line      Line;

  typedef  Nef_polyhedron::Object_handle Object_handle;
  typedef  Nef_polyhedron::Explorer Explorer;
  typedef  Explorer::Vertex_const_handle Vertex_const_handle;
  typedef  Explorer::Halfedge_const_handle Halfedge_const_handle;
  typedef  Explorer::Face_const_handle Face_const_handle;
  typedef  Explorer::Vertex_const_iterator Vertex_const_iterator;
  typedef  Explorer::Ray Ray;

  Point p1(0,0), p2(0,1), p3(1,0), p4(-1,-1), p5(0,-1), p6(-1,0), p7(1,1);
  Line l1(p2,p1); // neg y-axis
  Line l2(p1,p3); // pos x-axis
  Nef_polyhedron N1(l1), N2(l2, Nef_polyhedron::EXCLUDED),
                 EMPTY(Nef_polyhedron::EMPTY),PLANE(Nef_polyhedron::COMPLETE);
  CGAL_TEST((N1*N1) == N1);
  CGAL_TEST((N1*!N1) == EMPTY);
  CGAL_TEST((N1+!N1) == PLANE);
  CGAL_TEST((N1^N2) == ((N1-N2)+(N2-N1)));
  CGAL_TEST((!(N1*N2)) == (!N1+!N2));

  Nef_polyhedron N3 = N1.intersection(N2);
  // N3 is the first quadrant including the positive y-axis
   //  but excluding the origin and the positive x-axis 

  CGAL_TEST(N3 < N1 && N3 < N2);
  CGAL_TEST(N3 <= N1 && N3 <= N2);
  CGAL_TEST(N1 > N3 && N2 > N3);
  CGAL_TEST(N1 >= N3 && N2 >= N3);

  Explorer E = N3.explorer();
  Vertex_const_iterator v = E.vertices_begin(); 
  CGAL_TEST( !E.is_standard(v) && E.ray(v) == Ray(p1,p4) );
  Halfedge_const_handle e = E.first_out_edge(v);
  CGAL_TEST( E.is_frame_edge(e) );
  ++(++v); // third vertex
  CGAL_TEST( E.is_standard(v) && E.point(v) == p1 );

  Vertex_const_handle v1,v2;
  Halfedge_const_handle e1,e2;
  Face_const_handle f1,f2;
  Object_handle h1,h2,h3;
  h1 = N3.locate(p1);
  h2 = N3.locate(p1,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  CGAL_TEST( E.is_standard(v1) && E.point(v1) == p1 );
  h1 = N3.locate(p2);
  h2 = N3.locate(p2,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) );
  CGAL_TEST( (e1==e2 || e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.locate(p4);
  h2 = N3.locate(p4,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
               f1 == f2 && !E.mark(f1) );
  // shooting along angular bisector:
  h1 = N3.ray_shoot(p4,Direction(1,1));
  h2 = N3.ray_shoot(p4,Direction(1,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
             f1 == f2 && E.mark(f1) );
  // shooting along x-axis:
  h1 = N3.ray_shoot(p6,Direction(1,0));
  h2 = N3.ray_shoot(p6,Direction(1,0),Nef_polyhedron::NAIVE);
  CGAL_TEST( h1.empty() && h2.empty() );
  // shooting along y-axis:
  h1 = N3.ray_shoot(p5,Direction(0,1));
  h2 = N3.ray_shoot(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  e = e1;
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2||e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.ray_shoot_to_boundary(p5,Direction(0,1));
  h2 = N3.ray_shoot_to_boundary(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  h1 = N3.ray_shoot_to_boundary(p7,Direction(0,-1));
  h2 = N3.ray_shoot_to_boundary(p7,Direction(0,-1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2 || e1==E.twin(e2)) );

  std::list<Point> L;
  L.push_back(p1);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(p2);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(f1,h2) && !E.mark(f1) );
   
  L.push_back(p2);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );
    
  L.push_back(p3);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && E.mark(f1) );
  h3 = N3.locate(Point(1,1,3));
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );

  CGAL_IO_TEST(N1,N2,CGAL::IO::ASCII);
  //CGAL_IO_TEST(N1,N2,CGAL::IO::PRETTY);



}


{
  typedef  CGAL::Filtered_extended_homogeneous<Integer> EKernel;
  typedef  CGAL::Nef_polyhedron_2<EKernel> Nef_polyhedron;
  typedef  Nef_polyhedron::Point     Point;
  typedef  Nef_polyhedron::Direction Direction;
  typedef  Nef_polyhedron::Line      Line;

  typedef  Nef_polyhedron::Object_handle Object_handle;
  typedef  Nef_polyhedron::Explorer Explorer;
  typedef  Explorer::Vertex_const_handle Vertex_const_handle;
  typedef  Explorer::Halfedge_const_handle Halfedge_const_handle;
  typedef  Explorer::Face_const_handle Face_const_handle;
  typedef  Explorer::Vertex_const_iterator Vertex_const_iterator;
  typedef  Explorer::Ray Ray;

  Point p1(0,0), p2(0,1), p3(1,0), p4(-1,-1), p5(0,-1), p6(-1,0), p7(1,1);
  Line l1(p2,p1); // neg y-axis
  Line l2(p1,p3); // pos x-axis
  Nef_polyhedron N1(l1), N2(l2, Nef_polyhedron::EXCLUDED),
                 EMPTY(Nef_polyhedron::EMPTY),PLANE(Nef_polyhedron::COMPLETE);
  CGAL_TEST((N1*N1) == N1);
  CGAL_TEST((N1*!N1) == EMPTY);
  CGAL_TEST((N1+!N1) == PLANE);
  CGAL_TEST((N1^N2) == ((N1-N2)+(N2-N1)));
  CGAL_TEST((!(N1*N2)) == (!N1+!N2));

  Nef_polyhedron N3 = N1.intersection(N2);
  // N3 is the first quadrant including the positive y-axis
   //  but excluding the origin and the positive x-axis 

  CGAL_TEST(N3 < N1 && N3 < N2);
  CGAL_TEST(N3 <= N1 && N3 <= N2);
  CGAL_TEST(N1 > N3 && N2 > N3);
  CGAL_TEST(N1 >= N3 && N2 >= N3);

  Explorer E = N3.explorer();
  Vertex_const_iterator v = E.vertices_begin(); 
  CGAL_TEST( !E.is_standard(v) && E.ray(v) == Ray(p1,p4) );
  Halfedge_const_handle e = E.first_out_edge(v);
  CGAL_TEST( E.is_frame_edge(e) );
  ++(++v); // third vertex
  CGAL_TEST( E.is_standard(v) && E.point(v) == p1 );

  Vertex_const_handle v1,v2;
  Halfedge_const_handle e1,e2;
  Face_const_handle f1,f2;
  Object_handle h1,h2,h3;
  h1 = N3.locate(p1);
  h2 = N3.locate(p1,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  CGAL_TEST( E.is_standard(v1) && E.point(v1) == p1 );
  h1 = N3.locate(p2);
  h2 = N3.locate(p2,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) );
  CGAL_TEST( (e1==e2 || e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.locate(p4);
  h2 = N3.locate(p4,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
               f1 == f2 && !E.mark(f1) );
  // shooting along angular bisector:
  h1 = N3.ray_shoot(p4,Direction(1,1));
  h2 = N3.ray_shoot(p4,Direction(1,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
             f1 == f2 && E.mark(f1) );
  // shooting along x-axis:
  h1 = N3.ray_shoot(p6,Direction(1,0));
  h2 = N3.ray_shoot(p6,Direction(1,0),Nef_polyhedron::NAIVE);
  CGAL_TEST( h1.empty() && h2.empty() );
  // shooting along y-axis:
  h1 = N3.ray_shoot(p5,Direction(0,1));
  h2 = N3.ray_shoot(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  e = e1;
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2||e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.ray_shoot_to_boundary(p5,Direction(0,1));
  h2 = N3.ray_shoot_to_boundary(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  h1 = N3.ray_shoot_to_boundary(p7,Direction(0,-1));
  h2 = N3.ray_shoot_to_boundary(p7,Direction(0,-1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2 || e1==E.twin(e2)) );

  std::list<Point> L;
  L.push_back(p1);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(p2);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(f1,h2) && !E.mark(f1) );
   
  L.push_back(p2);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );
    
  L.push_back(p3);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && E.mark(f1) );
  h3 = N3.locate(Point(1,1,3));
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );

  CGAL_IO_TEST(N1,N2,CGAL::IO::ASCII);
  //CGAL_IO_TEST(N1,N2,CGAL::IO::PRETTY);



  Nef_polyhedron::EK.print_statistics();
}


{
  typedef  CGAL::Extended_cartesian<Rational> EKernel;
  typedef  CGAL::Nef_polyhedron_2<EKernel> Nef_polyhedron;
  typedef  Nef_polyhedron::Point     Point;
  typedef  Nef_polyhedron::Direction Direction;
  typedef  Nef_polyhedron::Line      Line;

  typedef  Nef_polyhedron::Object_handle Object_handle;
  typedef  Nef_polyhedron::Explorer Explorer;
  typedef  Explorer::Vertex_const_handle Vertex_const_handle;
  typedef  Explorer::Halfedge_const_handle Halfedge_const_handle;
  typedef  Explorer::Face_const_handle Face_const_handle;
  typedef  Explorer::Vertex_const_iterator Vertex_const_iterator;
  typedef  Explorer::Ray Ray;

  Point p1(0,0), p2(0,1), p3(1,0), p4(-1,-1), p5(0,-1), p6(-1,0), p7(1,1);
  Line l1(p2,p1); // neg y-axis
  Line l2(p1,p3); // pos x-axis
  Nef_polyhedron N1(l1), N2(l2, Nef_polyhedron::EXCLUDED),
                 EMPTY(Nef_polyhedron::EMPTY),PLANE(Nef_polyhedron::COMPLETE);
  CGAL_TEST((N1*N1) == N1);
  CGAL_TEST((N1*!N1) == EMPTY);
  Nef_polyhedron  negN1 = ! N1;
  Nef_polyhedron  N1pnegN1 = N1 + negN1; 
  CGAL_TEST(N1pnegN1 == PLANE);
  CGAL_TEST((N1^N2) == ((N1-N2)+(N2-N1)));
  Nef_polyhedron N1tN2 = N1 * N2;
  Nef_polyhedron negN1tN2 = !N1tN2;
  Nef_polyhedron negN1pnegP2 = !N1+!N2;
  CGAL_TEST((!(N1*N2)) == (!N1+!N2));

  Nef_polyhedron N3 = N1.intersection(N2);
  // N3 is the first quadrant including the positive y-axis
   //  but excluding the origin and the positive x-axis 

  CGAL_TEST(N3 < N1 && N3 < N2);
  CGAL_TEST(N3 <= N1 && N3 <= N2);
  CGAL_TEST(N1 > N3 && N2 > N3);
  CGAL_TEST(N1 >= N3 && N2 >= N3);

  Explorer E = N3.explorer();
  Vertex_const_iterator v = E.vertices_begin(); 
  CGAL_TEST( !E.is_standard(v) && E.ray(v) == Ray(p1,p4) );
  Halfedge_const_handle e = E.first_out_edge(v);
  CGAL_TEST( E.is_frame_edge(e) );
  ++(++v); // third vertex
  CGAL_TEST( E.is_standard(v) && E.point(v) == p1 );

  Vertex_const_handle v1,v2;
  Halfedge_const_handle e1,e2;
  Face_const_handle f1,f2;
  Object_handle h1,h2,h3;
  h1 = N3.locate(p1);
  h2 = N3.locate(p1,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  CGAL_TEST( E.is_standard(v1) && E.point(v1) == p1 );
  h1 = N3.locate(p2);
  h2 = N3.locate(p2,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) );
  CGAL_TEST( (e1==e2 || e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.locate(p4);
  h2 = N3.locate(p4,Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
               f1 == f2 && !E.mark(f1) );
  // shooting along angular bisector:
  h1 = N3.ray_shoot(p4,Direction(1,1));
  h2 = N3.ray_shoot(p4,Direction(1,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( CGAL::assign(f1,h1) && CGAL::assign(f2,h2) && 
             f1 == f2 && E.mark(f1) );
  // shooting along x-axis:
  h1 = N3.ray_shoot(p6,Direction(1,0));
  h2 = N3.ray_shoot(p6,Direction(1,0),Nef_polyhedron::NAIVE);
  CGAL_TEST( h1.empty() && h2.empty() );
  // shooting along y-axis:
  h1 = N3.ray_shoot(p5,Direction(0,1));
  h2 = N3.ray_shoot(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  e = e1;
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2||e1==E.twin(e2)) && E.mark(e1) );
  h1 = N3.ray_shoot_to_boundary(p5,Direction(0,1));
  h2 = N3.ray_shoot_to_boundary(p5,Direction(0,1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(v1,h1) && CGAL::assign(v2,h2) && v1 == v2 );
  h1 = N3.ray_shoot_to_boundary(p7,Direction(0,-1));
  h2 = N3.ray_shoot_to_boundary(p7,Direction(0,-1),Nef_polyhedron::NAIVE);
  CGAL_TEST( N3.contained_in_boundary(h1) && N3.contained_in_boundary(h2) );
  CGAL_TEST( CGAL::assign(e1,h1) && CGAL::assign(e2,h2) && 
             (e1==e2 || e1==E.twin(e2)) );

  std::list<Point> L;
  L.push_back(p1);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(p2);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(f1,h2) && !E.mark(f1) );
   
  L.push_back(p2);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );
    
  L.push_back(p3);
  N3 = Nef_polyhedron(L.begin(), L.end(), Nef_polyhedron::INCLUDED);
  E = N3.explorer();
  h1 = N3.locate(p1);
  h2 = N3.locate(CGAL::midpoint(p1,p2));
  h3 = N3.locate(p6);
  CGAL_TEST( CGAL::assign(v1,h1) && E.point(v1)==p1 && E.mark(v1) );
  CGAL_TEST( CGAL::assign(e1,h2) && E.mark(e1) );
  CGAL_TEST( CGAL::assign(f1,h3) && E.mark(f1) );
  h3 = N3.locate(Point(1,1,3));
  CGAL_TEST( CGAL::assign(f1,h3) && !E.mark(f1) );

  CGAL_IO_TEST(N1,N2,CGAL::IO::ASCII);
  //CGAL_IO_TEST(N1,N2,CGAL::IO::PRETTY);



}


  CGAL_TEST_END;
}

