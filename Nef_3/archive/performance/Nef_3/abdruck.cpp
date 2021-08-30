#ifdef CGAL_NEF3_USE_LEDA_INTEGER
#include <CGAL/leda_integer.h>
typedef leda_integer NT;
#endif

#ifdef CGAL_NEF3_USE_GMPZ
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz NT;
#endif

#ifdef CGAL_NEF3_USE_LAZY_EXACT_NT
#include <CGAL/Lazy_nt>
typedef CGAL::Lazy_nt<NT> RT;
#else
typedef NT RT;
#endif

#ifdef CGAL_NEF3_USE_SIMPLE_HOMOGENEOUS
#include <CGAL/Simple_homogeneous.h>
typedef CGAL::Simple_homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_HOMOGENEOUS
#include <CGAL/Homogeneous.h>
typedef CGAL::Homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_FILTERED_HOMOGENEOUS_3
#include <CGAL/Filtered_homogeneous_3.h>
typedef CGAL::Filtered_homogeneous_3<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_EXTENDED_HOMOGENEOUS
#include <CGAL/Extended_homogeneous.h>
typedef CGAL::Extended_homogeneous<RT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_LEDA_RATIONAL
#include <CGAL/leda_rational.h>
typedef leda_rational NT;
#endif

#ifdef CGAL_NEF3_USE_GMPQ
#include <CGAL/Gmpq.h>
typedef CGAL::Gmpq NT;
#endif

#ifdef CGAL_NEF3_USE_SIMPLE_CARTESIAN
#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<NT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_CARTESIAN
#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<NT> Kernel;
#endif

#ifdef CGAL_NEF3_USE_EXTENDED_CARTESIAN
#include <CGAL/Extended_cartesian.h>
typdef CGAL::Extended_cartesian<NT> Kernel;
#endif

#include <CGAL/basic.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vector_3 Vector_3;
typedef Nef_polyhedron::Aff_transformation_3 Aff_transformation_3;
typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;

bool cgal_nef3_timer_on = false;

RT x_max,y_max,z_max;

void transform_form(const Nef_polyhedron& N, Nef_polyhedron& F) {

  Vertex_const_iterator vi(N.vertices_begin());
  x_max = CGAL_NTS abs(vi->point().hx()/vi->point().hw());
  y_max = CGAL_NTS abs(vi->point().hy()/vi->point().hw());
  z_max = CGAL_NTS abs(vi->point().hz()/vi->point().hw());
  for(;vi != N.vertices_end();++vi) {
    if(CGAL_NTS abs(vi->point().hx()/vi->point().hw()) > x_max)
      x_max = CGAL_NTS abs(vi->point().hx()/vi->point().hw());
    if(CGAL_NTS abs(vi->point().hy()/vi->point().hw()) > y_max)
      y_max = CGAL_NTS abs(vi->point().hy()/vi->point().hw());
    if(CGAL_NTS abs(vi->point().hz()/vi->point().hw()) > z_max)
      z_max = CGAL_NTS abs(vi->point().hz()/vi->point().hw());
  }

  x_max*105/100;
  y_max*105/100;
  z_max*105/100;

  F.transform(Aff_transformation_3(x_max,0,0, 0,y_max,0, 0,0,z_max, 1));
  F.transform(Aff_transformation_3(CGAL::TRANSLATION,Vector_3(0,0,z_max)));
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);

  Nef_polyhedron N,F;

  std::ifstream Fin("centered_cube.nef3");
  Fin >> F;

  std::ifstream Nin(argv[1]);
  Nin >> N;

  transform_form(N,F);

  cgal_nef3_timer_on = true;
  F.difference(N);
};
