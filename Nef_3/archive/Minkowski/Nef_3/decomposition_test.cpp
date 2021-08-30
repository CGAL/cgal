#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/convex_decomposition_3.h>
#include <CGAL/Nef_3/is_reflex_sedge.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convexity_check_3.h>

#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>

// #define CGAL_WITH_LAZY_KERNEL
#ifdef CGAL_WITH_LAZY_KERNEL
  #include<CGAL/Gmpq.h>
  #include<CGAL/Lazy_kernel.h>
  typedef CGAL::Gmpq NT;
  //typedef leda_rational NT;
  typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
  #ifdef CGAL_USE_LEDA
    #include <CGAL/leda_integer.h>
    typedef leda_integer NT;
  #else
    #include<CGAL/Gmpz.h>
    typedef CGAL::Gmpz NT;
  #endif
  typedef CGAL::Homogeneous<NT> Kernel;
#endif

#ifdef CGAL_NEF_INDEXED_ITEMS
#include<CGAL/Nef_3/SNC_indexed_items.h>
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron_3;
#else
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron_3;
#endif

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==2);
  std::ifstream in(argv[1]);
  Nef_polyhedron_3 N;
  in >> N;

  convex_decomposition_3<Nef_polyhedron_3>(N);

  N.is_valid(0,0);

  Nef_polyhedron_3::SHalfedge_const_iterator cse;
  CGAL_forall_shalfedges(cse, N)
    if(cse->incident_sface()->mark())
      CGAL_assertion(!CGAL::is_reflex_sedge_in_any_direction<Nef_polyhedron_3>(cse));

  Nef_polyhedron_3::Volume_const_iterator ci;
  CGAL_forall_volumes(ci, N) {
    if(!ci->mark()) continue;
    ci != N.volumes_begin();
    CGAL_assertion(++ci->shells_begin() == ci->shells_end());
    // TODO: test whether shell is outer shell
    Nef_polyhedron_3::SFace_const_handle sf(ci->shells_begin());
    CGAL::Polyhedron_3<Kernel> P;
    N.convert_inner_shell_to_polyhedron(sf, P);
    CGAL::is_strongly_convex_3(P);
  }

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron_3>(N);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
