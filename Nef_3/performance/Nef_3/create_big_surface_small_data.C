#include <CGAL/leda_integer.h>
typedef leda_integer RT;

#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;

#include <CGAL/basic.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include "tetrahedron_generator2.h"
#include <CGAL/Random.h>

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;
typedef tetrahedron_generator<Kernel> tgen;

bool cgal_nef3_timer_on = false;

RT x_min,x_max,y_min,y_max,z_min,z_max;

void transform_big(Nef_polyhedron& N,int n, int s) {
  Vertex_const_iterator vi(N.vertices_begin());
  x_min = vi->point().hx()/vi->point().hw();
  x_max = vi->point().hx()/vi->point().hw();
  y_min = vi->point().hy()/vi->point().hw();
  y_max = vi->point().hy()/vi->point().hw();
  z_min = vi->point().hz()/vi->point().hw();
  z_max = vi->point().hz()/vi->point().hw();
  for(;vi != N.vertices_end();++vi) {
    if(vi->point().hx()/vi->point().hw() < x_min) 
      x_min = vi->point().hx()/vi->point().hw();
    if(vi->point().hx()/vi->point().hw() > x_max) 
      x_max = vi->point().hx()/vi->point().hw();
    if(vi->point().hy()/vi->point().hw() < y_min) 
      y_min = vi->point().hy()/vi->point().hw();
    if(vi->point().hy()/vi->point().hw() > y_max) 
      y_max = vi->point().hy()/vi->point().hw();
    if(vi->point().hz()/vi->point().hw() < z_min) 
      z_min = vi->point().hz()/vi->point().hw();
    if(vi->point().hz()/vi->point().hw() > z_max) 
      z_max = vi->point().hz()/vi->point().hw();
  }

  N.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(-x_min,-y_min,-z_max)));
  N.transform(Aff_transformation_3(CGAL::SCALING, n*s+2,x_max-x_min));
  N.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0,0,s/2)));
}

void transform_small(Nef_polyhedron& N, int s, int l, int x, int y) {
  N.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(-x_min,-y_min,-z_min))); 
  N.transform(Aff_transformation_3(CGAL::SCALING, s*l, x_max-x_min));
  N.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3((2*x+1)*s/2,(2*y+1)*s/2,s/4)));
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc == 6);

  int n = argc>2 ? std::atoi(argv[2]) : 10;
  int s = argc>3 ? std::atoi(argv[3]) : 100;
  int l = argc>4 ? std::atoi(argv[4]) : 2;

  std::ostringstream out1;
  if(argc>5) {
    tgen t1(out1,s,std::atoi(argv[5]));
    t1.create_tetrahedra(n,n,1);
  } else {
    tgen t1(out1,s);
    t1.create_tetrahedra(n,n,1);
  }
  std::istringstream in1(out1.str());
  Nef_polyhedron N1;
  in1 >> N1;
  CGAL_assertion(N1.is_valid());

  Nef_polyhedron N2;
  std::ifstream in2(argv[1]);
  in2 >> N2;
  Nef_polyhedron N3=N2;
  transform_big(N2,n,s);
  
  CGAL::Random r;
  int x=r.get_int(0,n-l-1);
  int y=r.get_int(0,n-1-1);
  transform_small(N3,s,l,x,y);


};
