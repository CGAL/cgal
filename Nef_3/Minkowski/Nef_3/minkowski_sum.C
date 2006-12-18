#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_2/Object_index.h>
#include <fstream>
#include <map>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Nef_S2/gausian_map_to_polyhedron_3.h>
#include <CGAL/Nef_S2/gausian_map_to_nef_3.h>
#include <CGAL/Nef_3/convex_decomposition_3.h> 
#include <CGAL/convexity_check_3.h>

#include <CGAL/Nef_3/Nary_union_by_queue.h>
#include <CGAL/Nef_3/Nary_union_by_pq.h>

#ifdef CGAL_WITH_LAZY_KERNEL
typedef CGAL::Gmpq NT;
//typedef leda_rational NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
#endif
#ifdef CGAL_NEF_INDEXED_ITEMS
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron;
#else
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
#endif
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Nef_polyhedron::SNC_structure  SNC_structure;
typedef CGAL::SNC_decorator<SNC_structure>  SNC_decorator;
typedef Nef_polyhedron::Halfedge_iterator  Halfedge_iterator;
typedef SNC_structure::Halfedge_handle Halfedge_handle;
typedef Nef_polyhedron::SHalfedge_iterator SHalfedge_iterator;
typedef Nef_polyhedron::SHalfedge_handle SHalfedge_handle;
typedef Nef_polyhedron::Sphere_segment Sphere_segment;
typedef Nef_polyhedron::Sphere_point   Sphere_point;
typedef Nef_polyhedron::Volume_const_iterator  Volume_const_iterator;
typedef Nef_polyhedron::Vertex_const_handle  Vertex_const_handle;
typedef Nef_polyhedron::Halfedge_const_handle  Halfedge_const_handle;
typedef Nef_polyhedron::Halffacet_const_handle  Halffacet_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle  SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfloop_const_handle  SHalfloop_const_handle;
typedef Nef_polyhedron::SFace_const_handle  SFace_const_handle;
typedef Nef_polyhedron::Halffacet_cycle_const_iterator      Halffacet_cycle_const_iterator;
typedef Nef_polyhedron::SHalfedge_around_facet_const_circulator      SHalfedge_around_facet_const_circulator;
typedef Nef_polyhedron::Vector_3           Vector_3;
typedef Kernel::Plane_3            Plane_3;
typedef Kernel::Line_3             Line_3;
typedef Kernel::Point_3            Point_3;

typedef CGAL::Nary_union_by_queue<Nef_polyhedron> NUBQ;
typedef CGAL::Nary_union_by_pq<Nef_polyhedron> NUBPQ;
typedef CGAL::Nary_union_by_summup<Nef_polyhedron> NUBPS;

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
      //		<< CGAL::to_double(v->point().y()) << " "
      //		<< CGAL::to_double(v->point().z()) << std::endl;
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
	//	std::cerr << " " << VI[hc_start->source()->center_vertex()];
	//	std::cerr << "  add vertex " << hc_start->source()->center_vertex()->point()
	//		  << " with index " << VI[hc_start->source()->center_vertex()] << std::endl;
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
    //	      << SV.vs << " " << SV.fs << " 0" << std::endl;
      
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

int main(int argc, char* argv[]) {

  CGAL_assertion(argc==3);
  std::ifstream in(argv[1]);
  Nef_polyhedron N;
  in >> N;
  CGAL::Timer t1, t2, t3, t4;
  
  t1.start();
  SNC_decorator D(*const_cast<SNC_structure*>(N.sncp()));

  t2.start();
  //  CGAL_NEF_SETDTHREAD(227*229*233);
  convex_decomposition_3<Nef_polyhedron>(N);

  t2.stop();

  CGAL_assertion(N.is_valid());

  typedef CGAL::Gausian_map<Kernel> Gausian_map;

  std::ifstream inc(argv[2]);
  Nef_polyhedron NC;
  inc >> NC;
  Gausian_map GC(NC, --NC.volumes_end());

#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
  NUBS nary_union;
#elif defined CGAL_NEF3_NARY_UNION_VIA_PQ
  NUBPQ nary_union;
#else
  NUBQ nary_union;
#endif

  t3.start();

  //  int skip_shells = 0;
  int shells = N.number_of_volumes();
  Volume_const_iterator c = N.volumes_begin();
  ++c;
  for(;c!=N.volumes_end();++c) {
    std::cerr << "noch " << --shells << " shells" << std::endl;
//    if(shells == 805) CGAL_NEF_SETDTHREAD(223);
    //    if(skip_shells > 0) { --skip_shells; continue;}
    if(c->mark() == false) continue;

    Polyhedron P;
    convert_volume2polyhedron(N,c,P);
    Nef_polyhedron NP(P);

    Gausian_map G(NP,--NP.volumes_end());
    //    Gausian_map G(N, c);
    Gausian_map GcG;
    GcG.minkowski_sum(GC,G);

    Polyhedron tmp;
    gausian_map_to_polyhedron_3<Kernel, Polyhedron::HDS> Converter(GcG);
    tmp.delegate(Converter);
    CGAL_assertion(is_strongly_convex_3(tmp));
    Nef_polyhedron Ntmp(tmp);
    CGAL_assertion(Ntmp.is_valid());
    nary_union.add_polyhedron(Ntmp);
  }
  
  t3.stop();
  t4.start();
  Nef_polyhedron result = nary_union.get_union();
  t4.stop();
  t1.stop();

  //  std::cerr << result;

  std::cout << "Total runtime: " << t1.time() << std::endl;
  std::cout << "Decomposition: " << t2.time() << std::endl;
  std::cout << "Sum of convex Minkowski sums : " << t3.time() << std::endl;
  std::cout << "Union of subpolyhedra: " << t4.time() << std::endl;

  /*
  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();
*/

}
