#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/test_macros.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#else
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif

typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Sphere_point   Sphere_point;
typedef Nef_polyhedron::Sphere_segment Sphere_segment;
typedef Nef_polyhedron::Sphere_circle  Sphere_circle;
typedef Nef_polyhedron::Explorer Explorer;
typedef Explorer::SVertex_const_handle SVertex_const_handle;
typedef Explorer::SHalfedge_const_handle SHalfedge_const_handle;
typedef Explorer::SFace_const_handle SFace_const_handle;
typedef Nef_polyhedron::Object_handle Object_handle;

int main(int, char**)
{
  CGAL::set_pretty_mode ( std::cerr );
  CGAL_NEF_SETDTHREAD(911);
  std::cerr << CGAL::sweepversion << std::endl;

  // Sphere_map 109
  // Sphere_geometry 113
  // Nef_polyhedron_S2 121
  // SM_decoration 127
  // SM_overlayer 131
  // SM_triangulator 137
  // SM_constrained_triang_traits 139
  // SM_point_locator 143

  CGAL_TEST_START;    
  Sphere_point p1(0,1,0), p2(1,1,0), p3(1,1,1);
  Sphere_point p4(0,-1,0);
  Sphere_circle c1(1,0,0), c2(0,1,0), c3(0,0,1); 
  Sphere_segment s1(p4,p1,c3), s2(p1,p4,c1);
  Sphere_segment S[2] = { s1, s2 };
  Sphere_circle C[3] = { c1, c2, c3 };
  Nef_polyhedron N1(c1), N2(c3, Nef_polyhedron::EXCLUDED),
                 EMPTY(Nef_polyhedron::EMPTY),
                 SPHERE(Nef_polyhedron::COMPLETE), N11(N1), N22(N2);

  // missing construction from segments and random init

  CGAL_TEST((N1*N1) == N1){}
  CGAL_TEST((N1*!N1) == EMPTY){}
  CGAL_NEF_TRACEV(N1);
  CGAL_NEF_TRACEV(!N1);
  CGAL_NEF_TRACEV(N1+!N1);
  CGAL_TEST((N1+!N1) == SPHERE){}
  CGAL_TEST((N1^N2) == ((N1-N2)+(N2-N1))){} // xor reformulation

  CGAL_NEF_TRACEV((N1*N2)); CGAL_NEF_TRACEV(!(N1*N2)); CGAL_NEF_TRACEV((!N1+!N2)); 
  CGAL_NEF_TRACEV(!(N1*N2) ^ (!N1+!N2));
  CGAL_TEST( (!(N1*N2)) == (!N1+!N2) ){} // deMorgan
#if 1
  Nef_polyhedron N3 = N1.intersection(N2);
  /* N3 is the two octants +++ and +-+ including the z-positive yz-halfplane
     but excluding the y-axis and the x-positive xy-halfplane */

  Nef_polyhedron N4(S,S+2,Nef_polyhedron::INCLUDED);
  Nef_polyhedron N5(C,C+3,0.5);

  CGAL_TEST(N3 < N1 && N3 < N2){}
  CGAL_TEST(N3 <= N1 && N3 <= N2){}
  CGAL_TEST(N1 > N3 && N2 > N3){}
  CGAL_TEST(N1 >= N3 && N2 >= N3){}

  SVertex_const_handle v;
  SHalfedge_const_handle e;
  SFace_const_handle f;
  Object_handle h;
  h = N3.locate(p1); // on y-axis
  CGAL_TEST( CGAL::assign(v,h) ){}
  CGAL_TEST( v->point() == p1 && !v->mark() ){}
  h = N3.locate(p2);
  CGAL_TEST( CGAL::assign(e,h) ){}
  CGAL_TEST( e->circle() == c3 && !e->mark() ){}
  h = N3.locate(p3);
  CGAL_TEST( CGAL::assign(f,h) ){}
  CGAL_TEST( f->mark() ){}
#endif
  CGAL_TEST_END;     
  return 0;
}


