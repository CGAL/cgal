#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/Single_wall_creator2.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <fstream>
#include <map>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
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
typedef Nef_polyhedron::Vector_3           Vector_3;
typedef Kernel::Plane_3            Plane_3;
typedef Kernel::Line_3             Line_3;
typedef Kernel::Point_3            Point_3;

typedef CGAL::Single_wall_creator<Nef_polyhedron> Single_wall;
typedef CGAL::Single_wall_creator2<Nef_polyhedron> Single_wall2;
typedef CGAL::Ray_hit_generator<Nef_polyhedron> Ray_hit;
typedef CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;

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

  Ray_hit rh(Vector_3(-1,0,0));
  N.delegate(rh);

  Halfedge_iterator e = D.halfedges_begin();
  for(;e != D.halfedges_end(); ++e) {
    if(e->is_twin() && e->twin() != Halfedge_handle()) {

      //      std::cerr << "edge: " << e->source()->point();
      //      std::cerr << "->" << e->twin()->source()->point() << std::endl;

      if(e->point().hx() > 0) {
	Single_wall W(e->twin(),Vector_3(-1,0,0));
	N.delegate(W);
      } else {
	Single_wall W(e,Vector_3(-1,0,0));
	N.delegate(W);
      }
    }
  }

  //  CGAL_NEF_SETDTHREAD(43);

  External_structure_builder esb;
  N.delegate(esb);

  Ray_hit rh2(Vector_3(1,0,0));
  N.delegate(rh2);

  for(e=D.halfedges_end();e!=D.halfedges_begin();--e) {  // Vorsicht mit begin() und end()
    if(e->is_twin() && e->twin() != Halfedge_handle()) {

      //      std::cerr << "edge2: " << e->source()->point() << "->" 
      //		<< e->twin()->source()->point() << std::endl;
      
      if(e->point().hx() < 0) {
	Single_wall W(e->twin(),Vector_3(1,0,0));
	N.delegate(W);
      } else {
	Single_wall W(e,Vector_3(1,0,0));
	N.delegate(W);
      }
    }
  }

  N.delegate(esb);
  
  //  std::cerr << N;

  SHalfedge_iterator se;
  for(se=D.shalfedges_begin();se!=D.shalfedges_end();++se) {
    Sphere_segment s(se->source()->point(), se->twin()->source()->point(), se->circle());
    if(se->incident_sface()->mark() == true && s.is_long() && se->circle().a() != 0) {

      /*
     std::cerr << "sedge at " << normalized(se->source()->source()->point()) 
		<< " in plane " << normalized(se->circle()) << std::endl;
      std::cerr << "sedge " << se->source()->point()
		<< "->" << se->twin()->source()->point() << std::endl;   
      */

      Plane_3 pl1(se->circle()), pl2(0,0,1,0);
      Line_3 l;
      CGAL::Object result = intersection(pl1,pl2);
      CGAL_assertion(assign(l,result));
    
      //      std::cerr << "intersection line " << l << std::endl;

      Vector_3 vec(l.to_vector());
      Sphere_point ip(CGAL::ORIGIN+vec);
      if(ip.hy() < 0) ip = ip.antipode();

      SHalfedge_handle sec = se;
      do {
	sec = sec->snext();
	vec = sec->source()->point() - CGAL::ORIGIN;
      } while(sec!=se && vec != Vector_3(1,0,0) && vec != Vector_3(-1,0,0));
      
      //      std::cerr << "senkrecht " << vec << std::endl;
      //      std::cerr << "sec " << sec->source()->point()
      //		<< "->" << sec->twin()->source()->point() << std::endl;

      if(s.has_on(ip) && s.source() != ip && s.target() != ip) {
	//	std::cerr << "intersection point " << ip << std::endl;

	Single_wall2 W(sec->source(), ip);
	N.delegate(W);
      }
     
      if(s.has_on(ip.antipode()) && 
	 s.source() != ip.antipode() && 
	 s.target() != ip.antipode()) {
	//	std::cerr << "antipode intersection point " << ip.antipode() << std::endl;
	Single_wall2 W(sec->source(), ip.antipode());
	N.delegate(W);
      }
    }
  }

  N.delegate(esb);

  t2.stop();

  //  std::cerr << N;
  CGAL_assertion(N.is_valid());

  /*
  QApplication b(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* wb = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
  b.setMainWidget(wb);
  wb->show();
  b.exec();
  */

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
    //    if(skip_shells > 0) { --skip_shells; continue;}
    if(c->mark() == false) continue;

    Gausian_map G(N, c);
    Gausian_map GcG;
    GcG.minkowski_sum(GC,G);

    Polyhedron tmp;
    gausian_map_to_polyhedron_3<Kernel, Polyhedron::HDS> Converter(GcG);
    tmp.delegate(Converter);
    // CGAL_assertion(is_strongly_convex_3(tmp));
    Nef_polyhedron Ntmp(tmp);
    //    CGAL_assertion(Ntmp.is_valid());

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
