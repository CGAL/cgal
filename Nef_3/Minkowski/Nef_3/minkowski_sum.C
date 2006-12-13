#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/Single_wall_creator2.h>
#include <CGAL/Nef_3/YVertical_wall_builder.h>
#include <CGAL/Nef_3/Reflex_edge_searcher.h>
#include <CGAL/Nef_3/Ray_hit_generator.h>
#include <CGAL/Nef_3/Ray_hit_generator2.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/Edge_sorter.h>
#include <CGAL/Nef_3/Edge_sorter2.h>
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

#include <CGAL/convexity_check_3.h>

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
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

typedef CGAL::Single_wall_creator<Nef_polyhedron>  Single_wall;
typedef CGAL::Single_wall_creator2<Nef_polyhedron> Single_wall2;
typedef CGAL::YVertical_wall_builder<Nef_polyhedron> YVertical_wall_builder;
typedef CGAL::Reflex_edge_searcher<Nef_polyhedron> Reflex_edge_searcher;
typedef Reflex_edge_searcher::Reflex_sedge_iterator Reflex_sedge_iterator;
typedef Reflex_edge_searcher::Container             Container;
typedef CGAL::Ray_hit_generator<Nef_polyhedron> Ray_hit;
typedef CGAL::Ray_hit_generator2<Nef_polyhedron> Ray_hit2;
typedef CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;
typedef CGAL::Edge_sorter<Nef_polyhedron, Container> Edge_sorter;
typedef CGAL::Edge_sorter2<Nef_polyhedron, Container> Edge_sorter2;


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

  External_structure_builder esb;

  Reflex_edge_searcher res(Sphere_point(1,0,0));
  N.delegate(res,false,false);

  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;

  Reflex_sedge_iterator rei;

  while(res.negative_rsedges_begin() != res.negative_rsedges_end()) {
  for(rei=res.negative_rsedges_begin(); rei!=res.negative_rsedges_end(); ++rei)
      CGAL_assertion(res.is_reflex_sedge(*rei));


  for(rei=res.negative_rsedges_begin(); rei!=res.negative_rsedges_end(); ++rei) {
      SHalfedge_handle se(*rei);
      Halfedge_handle split_edge;
    Ray_hit2 rh2a(Vector_3(-1,0,0),se->source()->source());
    N.delegate(rh2a);
    if(rh2a.split_edge(split_edge))
      res.handle_new_edge(split_edge);
    Ray_hit2 rh2b(Vector_3(-1,0,0),se->source()->twin()->source());
    N.delegate(rh2b);
    if(rh2b.split_edge(split_edge))
      res.handle_new_edge(split_edge);
  }
  
  Edge_sorter es(res.get_negative_rsedges());
  N.delegate(es);

  //  CGAL_NEF_SETDTHREAD(227*229*47);
//  CGAL_NEF_SETDTHREAD(233*229*227);
  for(rei=res.negative_rsedges_begin(); rei!=res.negative_rsedges_end(); ++rei) {
    Halfedge_handle e = (*rei)->source();
//    std::cerr << "handle reflex edge " << e->source()->point() << "->" 
//    	      << e->twin()->source()->point() << std::endl;
    CGAL_assertion(res.is_reflex_sedge(*rei));    
    if(e->point().hx() > 0)
      e = e->twin();
    Single_wall W(e,Vector_3(-1,0,0));
    N.delegate(W);
  }    

  N.delegate(esb);

  //  std::cerr << N;

//  Reflex_edge_searcher res2;
  N.delegate(res, false, false);

  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;
  }
  
  for(rei=res.positive_rsedges_begin(); rei!=res.positive_rsedges_end(); ++rei) {
      SHalfedge_handle se(*rei);
      Halfedge_handle split_edge;
    Ray_hit2 rh2a(Vector_3(1,0,0),se->source()->source());
    N.delegate(rh2a);
    if(rh2a.split_edge(split_edge))
      res.handle_new_edge(split_edge);
    Ray_hit2 rh2b(Vector_3(1,0,0),se->source()->twin()->source());
    N.delegate(rh2b);
    if(rh2b.split_edge(split_edge))
      res.handle_new_edge(split_edge);
  }

  Edge_sorter2 es2(res.get_positive_rsedges());
  N.delegate(es2);
  

//  QApplication b(argc, argv);
//  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* wb = 
//    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
//  b.setMainWidget(wb);
//  wb->show();
//  b.exec();


  //  CGAL_NEF_SETDTHREAD(229*233);

  rei=res.positive_rsedges_end();
  if(rei!=res.positive_rsedges_begin())
    do {
      --rei;
      Halfedge_handle e = (*rei)->source();
      //      std::cerr << "handle reflex edge " << e->source()->point() << "->" 
		//		<< e->twin()->source()->point() << std::endl;
      if(e->point().hx() < 0)
	e = e->twin();
      //      std::cerr << "handle reflex edge " << e->source()->point() << "->" 
      //		<< e->twin()->source()->point() << std::endl;
      Single_wall W(e,Vector_3(1,0,0));
      N.delegate(W);
    } while(rei!=res.positive_rsedges_begin());

  N.delegate(esb);

  N.delegate(res, false, false);
  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;


  YVertical_wall_builder Y;
  N.delegate(Y,false,false);

  YVertical_wall_builder::Vertical_redge_iterator vri=Y.pos_begin();
  for(; vri != Y.pos_end(); ++vri) {
    //    std::cerr << "pos: " << (*vri)->source()->point()
    //	      << "->" << (*vri)->twin()->source()->point() << std::endl;
    Single_wall2 W((*vri),Sphere_point(0,1,0));
    N.delegate(W);
  }
  for(vri = Y.neg_begin(); vri != Y.neg_end(); ++vri) {
    //    std::cerr << "neg: " << (*vri)->source()->point()
    //	      << "->" << (*vri)->twin()->source()->point() << std::endl;
    Single_wall2 W((*vri),Sphere_point(0,-1,0));
    N.delegate(W);
  }      

  N.delegate(esb);

  t2.stop();

  N.delegate(res, false, false);
  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;

  //  CGAL_NEF_SETDTHREAD(293);
  //  std::cerr << N;
  CGAL_assertion(N.is_valid());


//  QApplication b(argc, argv);
//  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* wb = 
//    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
//  b.setMainWidget(wb);
//  wb->show();
//  b.exec();

  typedef CGAL::Gausian_map<Kernel> Gausian_map;

  std::ifstream inc(argv[2]);
  Nef_polyhedron NC;
  inc >> NC;
  Gausian_map GC(NC, --NC.volumes_end());

#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
  std::list<Nef_polyhedron> queue;
#elif defined CGAL_NEF3_NARY_UNION_VIA_QUEUE
  std::list<Nef_polyhedron> queue;
#else
  typedef std::multimap<Nef_polyhedron::Size_type,Nef_polyhedron> PQ;
  typedef PQ::iterator      PQ_iterator;
  PQ pq;
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

#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
    queue.push_back(Ntmp);
#elif defined CGAL_NEF3_NARY_UNION_VIA_QUEUE
    queue.push_back(Ntmp);
#else
    pq.insert(make_pair(Ntmp.number_of_vertices(),Ntmp));
#endif
  }
  
  t3.stop();
  t4.start();
#ifdef CGAL_NEF3_NARY_UNION_VIA_SUMMUP
  Nef_polyhedron result(Nef_polyhedron::EMPTY);
  std::list<Nef_polyhedron>::iterator i1;
  for(i1=queue.begin(); i1!=queue.end(); ++i1) {
    std::cerr << queue.size() << " polyhedra in the queue (SUMMUP)" << std::endl;
    result += *i1;
  }
#elif defined CGAL_NEF3_NARY_UNION_VIA_QUEUE
  std::list<Nef_polyhedron>::iterator i1,i2;
  while(queue.size() > 1) {
    std::cerr << queue.size() << " polyhedra in the queue" << std::endl;

    i1 = i2 = queue.begin();
    ++i2;

    Nef_polyhedron Ntmp(*i1 + *i2);
    
    //    CGAL_assertion(Ntmp.is_valid());
    queue.pop_front();
    queue.pop_front();
    queue.push_back(Ntmp);
  }

  Nef_polyhedron result = *queue.begin();
#else
  PQ_iterator i1, i2; 
  while(pq.size() > 1) {
    i1 = i2 = pq.begin();
    ++i2;

    std::cerr << pq.size() << " polyhedra in the priority queue " << i1->first << "," << i2->first << std::endl;

    Nef_polyhedron N1(i1->second);
    Nef_polyhedron N2(i2->second);

    Nef_polyhedron Ntmp(N1 + N2);
 
    //    CGAL_assertion(Ntmp.is_valid());
    pq.erase(i1);
    pq.erase(i2);
    pq.insert(make_pair(Ntmp.number_of_vertices(),Ntmp));
  }

  Nef_polyhedron result = pq.begin()->second;

#endif

  t4.stop();
  t1.stop();

  //  std::cerr << result;

  std::cout << "Total runtime: " << t1.time() << std::endl;
  std::cout << "Decomposition: " << t2.time() << std::endl;
  std::cout << "Sum of convex Minkowski sums : " << t3.time() << std::endl;
  std::cout << "Union of subpolyhedra: " << t4.time() << std::endl;

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(result);
  a.setMainWidget(w);
  w->show();
  a.exec();

}
