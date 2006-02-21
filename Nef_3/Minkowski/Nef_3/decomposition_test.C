#include <algorithm>
#include <map>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/leda_integer.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/Single_wall_creator2.h>
#include <CGAL/Nef_3/Ray_hit_generator.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/Edge_sorter.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <fstream>

typedef leda_integer NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron;
typedef Nef_polyhedron::SNC_structure  SNC_structure;
typedef CGAL::SNC_decorator<SNC_structure>  SNC_decorator;
typedef SNC_structure::Halfedge_handle  Halfedge_handle;
typedef Nef_polyhedron::Halfedge_iterator  Halfedge_iterator;
typedef Nef_polyhedron::SHalfedge_iterator SHalfedge_iterator;
typedef Nef_polyhedron::SHalfedge_handle SHalfedge_handle;
typedef Nef_polyhedron::SFace_handle SFace_handle;
typedef Nef_polyhedron::Object_handle Object_handle;
typedef Nef_polyhedron::Sphere_segment Sphere_segment;
typedef Nef_polyhedron::Sphere_point   Sphere_point;

typedef Nef_polyhedron::Vector_3           Vector_3;
typedef Kernel::Plane_3            Plane_3;
typedef Kernel::Line_3             Line_3;

typedef CGAL::Single_wall_creator<Nef_polyhedron> Single_wall;
typedef CGAL::Single_wall_creator2<Nef_polyhedron> Single_wall2;
typedef CGAL::Ray_hit_generator<Nef_polyhedron> Ray_hit;
typedef CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;
typedef CGAL::Edge_sorter<Nef_polyhedron> Edge_sorter;
typedef Nef_polyhedron::SNC_structure  SNC_structure;
typedef SNC_structure::Sphere_map      Sphere_map;
typedef CGAL::SM_decorator<Sphere_map>          SM_decorator;  
typedef CGAL::SM_point_locator<SM_decorator>    SM_point_locator;

int main(int argc, char* argv[]) {
  
  CGAL_assertion(argc==2);
  std::ifstream in(argv[1]);
  Nef_polyhedron N;
  in >> N;

  SNC_decorator D(*const_cast<SNC_structure*>(N.sncp()));

  CGAL_NEF_SETDTHREAD(227*229*47);

#define CGAL_NEF3_DECOMPOSITION_WITH_VISUAL_WALLS
#ifdef CGAL_NEF3_DECOMPOSITION_WITH_VISUAL_WALLS

  Ray_hit rh(Vector_3(-1,0,0));
  N.delegate(rh);

  Halfedge_iterator e;
  for(e=D.halfedges_end();e != D.halfedges_begin(); --e) {
    if(e->is_twin() && e->twin() != Halfedge_handle()) {

      std::cerr << "edge: " << e->source()->point();
      std::cerr << "->" 
		<< e->twin()->source()->point() << std::endl;
      
      if(e->point().hx() > 0) {
	Single_wall W(e->twin(),Vector_3(-1,0,0));
	N.delegate(W);
      } else {
	Single_wall W(e,Vector_3(-1,0,0));
	N.delegate(W);
      }
    }
  }

  External_structure_builder esb;
  N.delegate(esb);

  Ray_hit rh2(Vector_3(1,0,0));
  N.delegate(rh2);

  //  Edge_sorter es;
  //  N.delegate(es);
  
  for(e=D.halfedges_end();e!=D.halfedges_begin();--e) {
    if(e->is_twin() && e->twin() != Halfedge_handle()) {

      std::cerr << "edge2: " << e->source()->point() << "->" 
		<< e->twin()->source()->point() << std::endl;
    
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
#else
  SHalfedge_iterator se;
  for(se=D.shalfedges_begin();se!=D.shalfedges_end();++se) {
    Sphere_segment s(se->source()->point(), se->twin()->source()->point(), se->circle());
    if(se->incident_sface()->mark() == true && s.is_long() && se->circle().a() != 0) {
      //      CGAL_assertion(se->circle().a() != 0);

      std::cerr << "sedge at " << normalized(se->source()->source()->point()) 
		<< " in plane " << normalized(se->circle()) << std::endl;
      std::cerr << "sedge " << se->source()->point()
		<< "->" << se->twin()->source()->point() << std::endl;      

      Plane_3 pl1(se->circle()), pl2(0,1,0,0);
      Line_3 l;
      CGAL::Object result = intersection(pl1,pl2);
      CGAL_assertion(assign(l,result));
    
      std::cerr << "intersection line " << l << std::endl;

      Vector_3 vec(l.to_vector());
      Sphere_point ip(CGAL::ORIGIN+vec);
      if(ip.hy() < 0) ip = ip.antipode();

      SHalfedge_handle sec = se;
      do {
	sec = sec->snext();
	vec = sec->source()->point() - CGAL::ORIGIN;
      } while(sec!=se && vec != Vector_3(1,0,0) && vec != Vector_3(-1,0,0));
      
      std::cerr << "senkrecht " << vec << std::endl;
      std::cerr << "sec " << sec->source()->point()
		<< "->" << sec->twin()->source()->point() << std::endl;

      if(s.has_on(ip) && s.source() != ip && s.target() != ip) {
	std::cerr << "intersection point " << ip << std::endl;
	Single_wall2 W(sec->source(), ip);
	N.delegate(W);
      }
     
      if(s.has_on(ip.antipode()) && 
	 s.source() != ip.antipode() && 
	 s.target() != ip.antipode()) {
	std::cerr << "antipode intersection point " << ip.antipode() << std::endl;
	Single_wall2 W(sec->source(), ip.antipode());
	N.delegate(W);
      }
    }
  }

  N.delegate(esb);
#endif


  //  std::cerr << N;
  /*
  SHalfedge_iterator se;
  for(se=D.shalfedges_begin();se!=D.shalfedges_end();++se) {
    Sphere_segment s(se->source()->point(), se->twin()->source()->point(), se->circle());
    if(se->incident_sface()->mark() == true && s.is_long() && se->circle().a() != 0) {
      //      CGAL_assertion(se->circle().a() != 0);

      std::cerr << "sedge at " << normalized(se->source()->source()->point()) 
		<< " in plane " << normalized(se->circle()) << std::endl;
      std::cerr << "sedge " << se->source()->point()
		<< "->" << se->twin()->source()->point() << std::endl;      

      Plane_3 pl1(se->circle()), pl2(0,0,1,0);
      Line_3 l;
      CGAL::Object result = intersection(pl1,pl2);
      CGAL_assertion(assign(l,result));
    
      std::cerr << "intersection line " << l << std::endl;

      Vector_3 vec(l.to_vector());
      Sphere_point ip(CGAL::ORIGIN+vec);
      if(ip.hy() < 0) ip = ip.antipode();

      SHalfedge_handle sec = se;
      do {
	sec = sec->snext();
	vec = sec->source()->point() - CGAL::ORIGIN;
      } while(sec!=se && vec != Vector_3(1,0,0) && vec != Vector_3(-1,0,0));
      
      std::cerr << "senkrecht " << vec << std::endl;
      std::cerr << "sec " << sec->source()->point()
		<< "->" << sec->twin()->source()->point() << std::endl;

      if(s.has_on(ip) && s.source() != ip && s.target() != ip) {
	std::cerr << "intersection point " << ip << std::endl;
	Single_wall2 W(sec->source(), ip);
	N.delegate(W);
      }
     
      if(s.has_on(ip.antipode()) && 
	 s.source() != ip.antipode() && 
	 s.target() != ip.antipode()) {
	std::cerr << "antipode intersection point " << ip.antipode() << std::endl;
	Single_wall2 W(sec->source(), ip.antipode());
	N.delegate(W);
      }
    }
  }

  //  CGAL_NEF_SETDTHREAD(43);
  N.delegate(esb);
  //  CGAL_NEF_SETDTHREAD(1);
  */

  std::cerr << N;

  QApplication a(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
  a.setMainWidget(w);
  w->show();
  a.exec();
}
