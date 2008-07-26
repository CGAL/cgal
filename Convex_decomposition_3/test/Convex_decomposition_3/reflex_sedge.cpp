#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Nef_S2/Sphere_point.h>
#include <CGAL/Nef_S2/Sphere_circle.h>
#include <CGAL/Nef_S2/Sphere_segment.h>
#include <CGAL/Convex_decomposition_3/is_reflex_sedge.h>

template<typename Kernel>
class Test_SNC {

public:
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef CGAL::Sphere_point<Kernel> Sphere_point;
  typedef CGAL::Sphere_circle<Kernel> Sphere_circle;
  typedef CGAL::Sphere_segment<Kernel> Sphere_segment;

  struct Vertex {
    Point_3 pnt;
    
    Vertex(Point_3 p) : pnt(p) {}
    Point_3 point() {return pnt; }
  };

  struct Halfedge {
    Vertex* src;
    Sphere_point spt;

    Halfedge(Vertex* v, Sphere_point sp) : src(v), spt(sp) {}
    Sphere_point point() { return spt; }
    Vertex* source() { return src; }
    Halfedge* twin() { assert(false); return this; }
  };

  struct SHalfedge {
    Halfedge* src;
    SHalfedge* twn;
    SHalfedge* prv;
    Sphere_circle crc;

    SHalfedge(Halfedge* s) : src(s) {}
    Halfedge* source() { return src; }
    SHalfedge* twin() { return twn; }
    void set_twin(SHalfedge* t) {
      twn = t;
      crc = Sphere_circle(src->point(), twn->source()->point());
    }
    SHalfedge*& sprev() { return prv; }
    Sphere_circle circle() { return crc; }
  };

  typedef Halfedge* Halfedge_handle;
  typedef SHalfedge* SHalfedge_handle;
};


int main() {

  typedef CGAL::Homogeneous<CGAL::Gmpz> Kernel;
  typedef Test_SNC<Kernel> SNC;
  typedef SNC::Point_3 Point_3;
  typedef SNC::Sphere_point Sphere_point;
  typedef SNC::Vertex Vertex;
  typedef SNC::Halfedge Halfedge;
  typedef SNC::SHalfedge SHalfedge;

  Sphere_point dir(1,0,0);

  Vertex v0(Point_3(0,0,0));
  
  Halfedge e0(&v0, Sphere_point( 0,  0,  1));
  Halfedge e1(&v0, Sphere_point( 1,  0,  5));
  Halfedge e2(&v0, Sphere_point( 1, -1,  5));
  Halfedge e3(&v0, Sphere_point( 0, -1,  5));
  Halfedge e4(&v0, Sphere_point(-1, -1,  5));
  Halfedge e5(&v0, Sphere_point(-1,  0,  5));
  Halfedge e6(&v0, Sphere_point(-1,  1,  5));
  Halfedge e7(&v0, Sphere_point( 0,  1,  5));
  Halfedge e8(&v0, Sphere_point( 1,  1,  5));

  SHalfedge se1(&e0);
  SHalfedge st1(&e1);
  SHalfedge se2(&e0);
  SHalfedge st2(&e2);
  SHalfedge se3(&e0);
  SHalfedge st3(&e3);
  SHalfedge se4(&e0);
  SHalfedge st4(&e4);
  SHalfedge se5(&e0);
  SHalfedge st5(&e5);
  SHalfedge se6(&e0);
  SHalfedge st6(&e6);
  SHalfedge se7(&e0);
  SHalfedge st7(&e7);
  SHalfedge se8(&e0);
  SHalfedge st8(&e8);
  
  se1.set_twin(&st1);
  st1.set_twin(&se1);
  se2.set_twin(&st2);
  st2.set_twin(&se2);
  se3.set_twin(&st3);
  st3.set_twin(&se3);
  se4.set_twin(&st4);
  st4.set_twin(&se4);
  se5.set_twin(&st5);
  st5.set_twin(&se5);
  se6.set_twin(&st6);
  st6.set_twin(&se6);
  se7.set_twin(&st7);
  st7.set_twin(&se7);
  se8.set_twin(&st8);
  st8.set_twin(&se8);

  se1.sprev() = &se1;

  assert(CGAL::is_reflex_sedge<SNC>(&se1, dir, false) == 0);
}
