#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_2/Object_index.h>
#include <fstream>
#include <map>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
//#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_3/convex_decomposition_3.h>
#include <CGAL/Nef_3/bipartite_nary_union_sequential.h>
#include <CGAL/Nef_3/bipartite_nary_union_sorted_separately.h>
#include <CGAL/Nef_3/bipartite_nary_union_sorted_combined.h>
#include <CGAL/Nef_3/bipartite_nary_union_sorted_within_grid.h>

#define CGAL_NEF3_SPHERE_SWEEP_OPTIMIZATION_OFF

#ifdef CGAL_WITH_LAZY_KERNEL
#include <CGAL/Lazy_kernel.h>
typedef CGAL::Gmpq NT;
//typedef leda_rational NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
#endif
#ifdef CGAL_NEF_INDEXED_ITEMS
#include <CGAL/Nef_3/SNC_indexed_items.h>
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron;
#else
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
#endif
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Nef_polyhedron::Halfedge_iterator  Halfedge_iterator;
//typedef Nef_polyhedron::Halfedge_handle Halfedge_handle;
typedef Nef_polyhedron::SHalfedge_iterator SHalfedge_iterator;
typedef Nef_polyhedron::SHalfedge_handle SHalfedge_handle;
typedef Nef_polyhedron::Sphere_segment Sphere_segment;
typedef Nef_polyhedron::Sphere_point   Sphere_point;
typedef Nef_polyhedron::Volume_const_handle  Volume_const_handle;
typedef Nef_polyhedron::Vertex_const_handle  Vertex_const_handle;
typedef Nef_polyhedron::Halfedge_const_handle  Halfedge_const_handle;
typedef Nef_polyhedron::Halffacet_const_handle  Halffacet_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle  SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfloop_const_handle  SHalfloop_const_handle;
typedef Nef_polyhedron::SFace_const_handle  SFace_const_handle;
typedef Nef_polyhedron::Halffacet_cycle_const_iterator      Halffacet_cycle_const_iterator;
typedef Nef_polyhedron::SHalfedge_around_facet_const_circulator      SHalfedge_around_facet_const_circulator;

/*
template <class HDS, class Const_decorator>
class Build_polyhedron : public CGAL::Modifier_base<HDS> {

  class Statistic_visitor {
  public:
    int vs, es, fs;
    Statistic_visitor() : vs(0), es(0), fs(0) {}
    void visit(SFace_const_handle s) {}
    void visit(Halfedge_const_handle e) { if(e->is_twin()) ++es; }
    void visit(Halffacet_const_handle f) { ++fs; }
    void visit(Vertex_const_handle v) { ++vs; }
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
  };

  class Vertex_creator {
    int vertex_index;
    CGAL::Polyhedron_incremental_builder_3<HDS>* B;
    CGAL::Object_index<Vertex_const_handle>& VI;
  public:
    Vertex_creator(CGAL::Polyhedron_incremental_builder_3<HDS>* Bin,
                   CGAL::Object_index<Vertex_const_handle>& VIin)
      : vertex_index(0), B(Bin), VI(VIin) {}
    void visit(SFace_const_handle s) {}
    void visit(Halfedge_const_handle e) {}
    void visit(Halffacet_const_handle f) {}
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
    void visit(Vertex_const_handle v) {
      //      std::cerr << CGAL::to_double(v->point().x()) << " "
      //                << CGAL::to_double(v->point().y()) << " "
      //                << CGAL::to_double(v->point().z()) << std::endl;
      VI[v]=vertex_index++;
      B->add_vertex(v->point());
    }
  };

  class Facet_creator {
    CGAL::Polyhedron_incremental_builder_3<HDS>* B;
    CGAL::Object_index<Vertex_const_handle>& VI;
  public:
    Facet_creator(CGAL::Polyhedron_incremental_builder_3<HDS>* Bin,
                   CGAL::Object_index<Vertex_const_handle>& VIin)
    : B(Bin), VI(VIin) {}
    void visit(Halffacet_const_handle opposite_facet) {
      //      std::cerr << "visit facet " << opposite_facet->plane() << std::endl;
      SHalfedge_const_handle se;
      Halffacet_cycle_const_iterator fc;
      Halffacet_const_handle f = opposite_facet->twin();

      B->begin_facet();
      fc = f->facet_cycles_begin();
      se = SHalfedge_const_handle(fc);
      CGAL_assertion(se!=0);
      SHalfedge_around_facet_const_circulator hc_start(se);
      SHalfedge_around_facet_const_circulator hc_end(hc_start), hcx(hc_start);
      //      std::cerr << std::distance(++hcx,hc_end)+1;
      CGAL_For_all(hc_start,hc_end) {
        //        std::cerr << " " << VI[hc_start->source()->center_vertex()];
        //        std::cerr << "  add vertex " << hc_start->source()->center_vertex()->point()
        //                  << " with index " << VI[hc_start->source()->center_vertex()] << std::endl;
        B->add_vertex_to_facet(VI[hc_start->source()->center_vertex()]);
      }
      //      std::cerr << std::endl;
      B->end_facet();
    }

    void visit(SFace_const_handle s) {}
    void visit(Halfedge_const_handle e) {}
    void visit(Vertex_const_handle v) {}
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
  };

public:

  const Const_decorator& scd;
  Volume_const_iterator c;
  CGAL::Object_index<Vertex_const_handle> VI;
  CGAL::Polyhedron_incremental_builder_3<HDS>* B;

  Build_polyhedron(const Const_decorator& s, Volume_const_iterator cin) :
    scd(s), c(cin) {}

  void operator()( HDS& hds) {

    B= new CGAL::Polyhedron_incremental_builder_3<HDS>(hds, true);

    Statistic_visitor SV;
    scd.visit_shell_objects(SFace_const_handle(c->shells_begin()),SV);

    B->begin_surface(SV.vs,SV.fs,SV.es);
    //    std::cerr << "OFF " << std::endl
    //              << SV.vs << " " << SV.fs << " 0" << std::endl;

    Vertex_creator VC(B,VI);
    scd.visit_shell_objects(SFace_const_handle(c->shells_begin()),VC);

    Facet_creator FC(B,VI);
    scd.visit_shell_objects(SFace_const_handle(c->shells_begin()),FC);

    B->end_surface();
    delete B;
  }

};

void convert_volume2polyhedron(const Nef_polyhedron& N,
                               Volume_const_iterator c, Polyhedron& P) {

  typedef Polyhedron::HalfedgeDS HalfedgeDS;
  Build_polyhedron<HalfedgeDS,Nef_polyhedron> bp(N,c);
  P.delegate(bp);

}
*/

bool loadFile(char* filename, Nef_polyhedron& N) {
  std::ifstream in(filename);
  std::ifstream test(filename);
  char c;
  test >> c;
  if(c!='S' && c!='O') return false;
  if(c == 'S')
    in >> N;
  else {
    Polyhedron P;
    in >> P;
    N = Nef_polyhedron(P);
  }
  return true;
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==3);
  /*
  std::ifstream in(argv[1]);
  Nef_polyhedron N0;
  in >> N0;

  std::ifstream inc(argv[2]);
  Nef_polyhedron N1;
  inc >> N1;
  */

  Nef_polyhedron N0, N1;
  if(!loadFile(argv[1], N0)) {
    std::cerr << "parameter 1 is not a valid input file" << std::endl;
    return 0;
  }
  if(!loadFile(argv[2], N1)) {
    std::cerr << "parameter 2 is not a valid input file" << std::endl;
    return 0;
  }

  CGAL::Timer t1, t2, t3;
  t1.start();
  t2.start();
  //  CGAL_NEF_SETDTHREAD(503*509);
  convex_decomposition_3<Nef_polyhedron>(N0);
  t2.stop();
  std::cerr << "Decomposition of Polyhedron 1: " << t2.time() << std::endl;
  t3.start();
  convex_decomposition_3<Nef_polyhedron>(N1);
  t3.stop();
  std::cerr << "Decomposition of Polyhedron 2: " << t3.time() << std::endl;

  CGAL_assertion(N0.is_valid());
  CGAL_assertion(N1.is_valid());

  Nef_polyhedron result =
#ifdef CGAL_MINKOWSKI_BIPARTITE_NARY_UNION_SEQUENTIAL
    CGAL::bipartite_nary_union_sequential(N0, N1);
#elif defined CGAL_MINKOWSKI_BIPARTITE_NARY_UNION_SORTED_SEPARATELY
    CGAL::bipartite_nary_union_sorted_separately(N0, N1);
#elif defined CGAL_MINKOWSKI_BIPARTITE_NARY_UNION_SORTED_WITH_GRID
    CGAL::bipartite_nary_union_sorted_within_grid(N0, N1);
#else
    CGAL::bipartite_nary_union_sorted_combined(N0, N1);
#endif

    std::cerr << "Decomposition: " << t2.time()+t3.time() << std::endl;
    std::cerr << "Total runtime: " << t1.time() << std::endl;

    std::cout << result;

    /*
  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w =
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
    */
}
