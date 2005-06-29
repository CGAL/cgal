#ifndef CGAL_NEF3_SINGLE_WALL_CREATOR_H
#define CGAL_NEF3_SINGLE_WALL_CREATOR_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_walls.h>
#include <CGAL/Nef_3/Ray_hit_generator.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Single_wall_creator : public Modifier_base<typename Nef_::SNC_and_PL> {
  
  typedef Nef_                                    Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL     SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure  SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>      Base;
  typedef CGAL::SNC_point_locator<Base>           SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>   SNC_intersection;
  typedef CGAL::SNC_constructor<SNC_structure>    SNC_constructor;

  typedef typename SNC_structure::Sphere_map      Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>          SM_decorator;  
  typedef CGAL::SM_point_locator<SM_decorator>    SM_point_locator; 
  typedef CGAL::SM_walls<Sphere_map>              SM_walls;
  typedef CGAL::Ray_hit_generator<Nef_polyhedron> Ray_hit;  

  typedef typename Base::Segment_3               Segment_3;
  typedef typename Base::Point_3                 Point_3;
  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Sphere_circle           Sphere_circle;
  typedef typename Base::Sphere_segment          Sphere_segment;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::SVertex_iterator        SVertex_iterator;

  typedef typename Base::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;

  Halfedge_handle ein;
  Vector_3 dir;
  
 public:
  Single_wall_creator(Halfedge_handle e, Vector_3 d)
    : ein(e), dir(d) {}

 private:
  int need_to_create_wall() const {

    int res = 0;

    if((ein->point() - CGAL::ORIGIN) == dir ||
       (ein->twin()->point() - CGAL::ORIGIN) == dir)
      return 4;

    std::cerr << "test 0 " << std::endl;
    
    SHalfedge_handle se;
    SM_point_locator PS(&*ein->source());
    Object_handle os = PS.locate(Sphere_point(dir));
    if(assign(se,os) && (se->source() == ein || se->twin()->source() == ein))
      return 4;
    
    std::cerr << "test 1 " << std::endl;

    SM_point_locator PT(&*ein->target());
    Object_handle ot = PT.locate(Sphere_point(dir));
    if(assign(se,ot) && (se->source() == ein->twin() || se->twin()->source() == ein->twin()))
      return 4;

    std::cerr << "test 2 " << std::endl;

    /*
    SFace_handle sfs, sft;
    if(assign(sfs,os) && sfs->mark() == false)
      res |= 1;
    if(assign(sft,os) && sft->mark() == false)
      res |= 2;

    std::cerr << "test 3 " << std::endl;
    */
  
    /*
    SVertex_handle sv;
    if(assign(sv,os)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send)
	if(sh->twin()->source() == ein) // tighter, i.e. sedge with circle(svs,ein)
	  res |= 1;
    }

    std::cerr << "test 4 " << std::endl;

    if(assign(sv,ot)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send)
	if(sh->twin()->source() == ein->twin()) // tighter, i.e. sedge with circle(svs,ein)
	  res |= 2;
    }
    */

    /*
    SVertex_handle sv;
    if(assign(sv,os)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send) {
	Sphere_circle c(sv->point(), ein->point());
	c = normalized(c);
	if(sh->circle() == c)
	  res |= 1;
      }
    }

    if(assign(sv,ot)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send) {
	Sphere_circle c(sv->point(), ein->twin()->point());
	c = normalized(c);
	if(sh->circle() == c)
	  res |= 2;
      }
    }
*/
    // TODO: some cases missing!!!
    
    std::cerr << "result " << res << std::endl;

    return res;
  }
    
 public:
  void operator()(SNC_and_PL& sncpl) {

    if(Sphere_point(CGAL::ORIGIN - dir) == ein->point() ||
       Sphere_point(CGAL::ORIGIN - dir) == ein->twin()->point())
      return;

    SNC_structure* sncp(sncpl.sncp);
    SNC_point_locator* pl(sncpl.pl);

    SNC_constructor C(*sncp,pl);

    int ntcw = need_to_create_wall();
    std::cerr << "ntcw" << ntcw << std::endl;
    CGAL_assertion(ntcw >=0 && ntcw <5);
    if(ntcw > 0) return;

    Vertex_handle origin[2];
    //    Vertex_handle opposite[2];
    origin[0] = ein->source();
    origin[1] = ein->target();
    SVertex_handle lateral_sv_tgt[2];
    lateral_sv_tgt[0] = ein;
    lateral_sv_tgt[1] = ein->twin();

   
    bool legal[2];
    Object_handle found_object[2];
    Sphere_point found_point[2];
    for(int i=0; i<2; ++i) {
      SM_walls SMW(&*origin[i]);
      SM_point_locator PL(&*origin[i]);
      Sphere_segment sphere_ray(lateral_sv_tgt[i]->point(), Sphere_point(dir));
      legal[i] = SMW.legal_direction(sphere_ray, found_object[i], found_point[i]);
    }

    SVertex_handle sv0, sv1;
    if(assign(sv0, found_object[0]) && assign(sv1, found_object[1])) {
      SHalfedge_around_svertex_circulator sh0(sv0->out_sedge()), send0(sh0);
      CGAL_For_all(sh0,send0)
	if(sh0->twin()->source() == ein) {
	  SHalfedge_around_svertex_circulator sh1(sv1->out_sedge()), send1(sh1);
	  CGAL_For_all(sh1,send1)
	    if(sh1->twin()->source() == ein->twin())
	      return;
	}
    }

    // TODO: check if cycle already exists

    if(!legal[0] || !legal[1]) 
      return;

    Sphere_circle c(ein->point(), Sphere_point(dir));
    c = normalized(c);

    SM_walls SMW(&*ein->source());
    Sphere_segment sphere_ray(ein->point(), ein->twin()->point(), c);
    SVertex_handle lateral_svertex = SMW.add_lateral_svertex(sphere_ray);
    SMW.add_sedge_between(ein, lateral_svertex, c);

    Ray_hit rh(sncp, pl, 3);
    Ray_3 r(lateral_svertex->source()->point(), lateral_svertex->point()-CGAL::ORIGIN);
    Vertex_handle v = rh.create_vertex_on_first_hit(r);
    
    while(v != ein->twin()->source()) {
      
      std::cerr << "current vertex " << v->point() << std::endl;

      std::cerr << "double coords" << CGAL::to_double(v->point().x())
		<< ", " << CGAL::to_double(v->point().y())
		<< ", " << CGAL::to_double(v->point().z()) << std::endl;

      SM_walls smw(&*v);
      SVertex_handle opp = smw.add_ray_svertex(lateral_svertex->point().antipode());
      opp->twin() = lateral_svertex;
      lateral_svertex->twin() = opp;
      pl->add_edge(lateral_svertex);

      std::cerr << "twins " << lateral_svertex->source()->point() << " + " << opp->source()->point() << std::endl;

      sphere_ray = Sphere_segment(lateral_svertex->point().antipode(), 
				  lateral_svertex->point(), c);
      lateral_svertex = smw.add_lateral_svertex(sphere_ray);
      smw.add_sedge_between(opp, lateral_svertex, c);

      // TODO: make use of existing edges along ray

      r = Ray_3(lateral_svertex->source()->point(), lateral_svertex->point()-CGAL::ORIGIN);
      v = rh.create_vertex_on_first_hit(r);
    }
    
    std::cerr << "last current vertex " << v->point() << std::endl;

    SM_walls smw(&*v);
    SVertex_handle opp = smw.add_ray_svertex(lateral_svertex->point().antipode());    
    opp->twin() = lateral_svertex;
    lateral_svertex->twin() = opp;
    pl->add_edge(lateral_svertex);

    smw.add_sedge_between(opp, ein->twin(), c);

    std::cerr << "final twins " << lateral_svertex->source()->point() << " + " << opp->source()->point() << std::endl;

    /*
    bool check_goal;
    if(!legal[0] || !legal[1]) 
      return;
    else {
      SM_walls SMW_src(&*origin[0]);
      lateral_sv_tgt[0] = SMW_src.add_svertex_into_object(found_point[0], found_object[0]);

      c = Sphere_circle(ein->point(), Sphere_point(dir));
      c = normalized(c);

      Sphere_circle test(ein->point(), lateral_sv_tgt[0]->point());
      CGAL_assertion(normalized(test) == c || ein->point().antipode() == lateral_sv_tgt[0]->point());

      SMW_src.add_sedge_between(ein, lateral_sv_tgt[0],c);
      check_goal = (found_point[1] == Sphere_point(dir)) ;
    }
    
    Point_3 goal;
    Object_handle goal_object;
    */

    /*
    for(int i=1; i<2; ++i) {
      if(!check_goal) continue;

      Vertex_handle v;
      Halfedge_handle e;
      Halffacet_handle f;
      Ray_3 r(origin[i]->point(),dir);
      std::cerr << "shoot ray " << r << std::endl;
      Object_handle o2 = pl->shoot(r);
      goal_object = o2;
      if(assign(f,o2)) {
	std::cerr << "Found facet " << f->plane() << std::endl;

	SNC_intersection I;
	I.does_intersect_internally(r, f, goal);

	std::cerr << "goal coords" << CGAL::to_double(goal.x())
		  << ", " << CGAL::to_double(goal.y())
		  << ", " << CGAL::to_double(goal.z()) << std::endl;

      } else if(assign(e,o2)) {
	std::cerr << "Found edge " << e->source()->point() 
		  << "->" << e->twin()->source()->point() << std::endl;

	SNC_intersection I;
	I.does_intersect_internally(r, Segment_3(e->source()->point(),
						 e->twin()->source()->point()),  
				    goal);
	std::cerr << "goal coords" << CGAL::to_double(goal.x())
		  << ", " << CGAL::to_double(goal.y())
		  << ", " << CGAL::to_double(goal.z()) << std::endl;

      } else if(assign(v,o2)) {
	
	std::cerr << "Found vertex " << v->point() << std::endl;
	std::cerr << "double coords" << CGAL::to_double(v->point().x())
		  << ", " << CGAL::to_double(v->point().y())
		  << ", " << CGAL::to_double(v->point().z()) << std::endl;

	opposite[1] = v;
	goal = v->point();
      } else {
	std::cerr << "Found nothing " << std::endl;
	opposite[i] = origin[i];
      }
    }
    */
    /*
    for(int i=0; i<2; ++i) {
      if(origin[i] != opposite[i] || origin[1-i] == opposite[1-i]) 
	continue;

      CGAL_assertion_msg(false, "special not changed, yet");

      std::cerr << "special " << std::endl;
      SM_walls SMW(&*opposite[i]);
      
      SVertex_handle src(i==0?ein:ein->twin());
      Sphere_circle c(Sphere_point(CGAL::ORIGIN - dir), src->point());
      Sphere_segment sphere_ray(src->point(), Sphere_point(CGAL::ORIGIN - dir), c);
      lateral_sv_tgt[i] = SMW.add_lateral_svertex(sphere_ray, true, Sphere_point(dir));
      SMW.add_sedge_between(src, lateral_sv_tgt[i]);    
    }
    */
    /*
    do {
      Ray_3 r(lateral_sv_tgt[0]->source()->point(),lateral_sv_tgt[0]->point()-CGAL::ORIGIN);

      CGAL_NEF_SETDTHREAD(503*509);
      Object_handle o2 = pl->shoot(r);
      CGAL_NEF_SETDTHREAD(1);

      std::cerr << "nachbehandlung: ray " << r << std::endl;
	std::cerr << "double coords" << CGAL::to_double(r.source().x())
		  << ", " << CGAL::to_double(r.source().y())
		  << ", " << CGAL::to_double(r.source().z())
	          << "->" << CGAL::to_double(r.to_vector().x()) 
	          << ", " << CGAL::to_double(r.to_vector().y()) 
	          << ", " << CGAL::to_double(r.to_vector().z()) << std::endl;

      Halffacet_handle f;
      Halfedge_handle e;
      Vertex_handle v;
      if(assign(f,o2)) {
	std::cerr << "Found facet " << f->plane() << std::endl;

	std::cerr << "double coords" << CGAL::to_double(f->plane().a())
		  << ", " << CGAL::to_double(f->plane().b())
		  << ", " << CGAL::to_double(f->plane().c())
	          << ", " << CGAL::to_double(f->plane().d()) << std::endl;

	Point_3 ip;
	SNC_intersection I;
	I.does_intersect_internally(r, f, ip);

	Vertex_handle v = C.create_from_facet(f,ip);
	pl->add_vertex(v);

	std::cerr << "double coords" << CGAL::to_double(v->point().x())
		  << ", " << CGAL::to_double(v->point().y())
		  << ", " << CGAL::to_double(v->point().z()) << std::endl;


	SM_walls SMW_tgt(&*v);
	
	SVertex_handle ray_sv_src = lateral_sv_tgt[0];
	SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN - lateral_sv_tgt[0]->point()));
	ray_sv_src->twin() = ray_sv_tgt;
	ray_sv_tgt->twin() = ray_sv_src;	
	
	SVertex_handle lateral_svertex = lateral_sv_tgt[0];
	Sphere_segment sphere_ray(Sphere_point(CGAL::ORIGIN - lateral_svertex->point()), 
				  lateral_svertex->point(),c);
	lateral_sv_tgt[0] = SMW_tgt.add_lateral_svertex(sphere_ray, true, Sphere_point(dir));

	SMW_tgt.add_sedge_between(ray_sv_tgt, lateral_sv_tgt[0],c);

      } else if(assign(e,o2)) {
	std::cerr << "Found edge " << e->source()->point() 
		  << "->" << e->twin()->source()->point() << std::endl;
	Point_3 ip;
	SNC_intersection I;
	I.does_intersect_internally(r, Segment_3(e->source()->point(),
						 e->twin()->source()->point()),  
				    ip);
	ip = normalized(ip);
	//	std::cerr << ip << "=?=" << lateral_sv_tgt[1]->source()->point() << std::endl;

	if(check_goal) {
	  Segment_3 seg(lateral_sv_tgt[0]->source()->point(), ip);
	  if(seg.has_on(goal)) {
	    
	    if(assign(f,goal_object)) {
	      std::cerr << "Found goal on facet " << std::endl;
	      
	      Vertex_handle v = opposite[1] = C.create_from_facet(f,goal);
	      pl->add_vertex(v);
	      
	    } else if(assign(e,o2)) {
	      std::cerr << "Found goal on edge " << std::endl;
	      
	      Vertex_handle v = opposite[1] = C.create_from_edge(e,goal);
	      pl->add_vertex(v);
	      
	      SVertex_iterator svi = v->svertices_begin();
	      SVertex_handle svf = svi;
	      SVertex_handle svb = ++svi;
	      
	      if(svf->point() == e->point()) {
		svb->twin() = e;
		svf->twin() = e->twin();
		e->twin()->twin() = svf;
		e->twin() = svb;
	      } else {
		svf->twin() = e;
		svb->twin() = e->twin();
		e->twin()->twin() = svb;
		e->twin() = svf;
	      }
	      
	    } else if(assign(v,o2)) {
	      
	      std::cerr << "Found goal on vertex " << v->point() << std::endl;
	      opposite[1] = v;
	    } else 
	      CGAL_assertion_msg(false, "wrong handle");
	    
	    SM_walls SMW_src(&*origin[1]);
	    SM_walls SMW_tgt(&*opposite[1]);
	    
	    SVertex_handle ray_sv_src = SMW_src.add_ray_svertex(Sphere_point(dir));
	    SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN -dir));
	    ray_sv_src->twin() = ray_sv_tgt;
	    ray_sv_tgt->twin() = ray_sv_src;
	    
	    lateral_sv_tgt[1] = SMW_tgt.add_ray_svertex(lateral_sv_tgt[0]->point().antipode()); 
	    
	    SMW_src.add_sedge_between(ray_sv_src, ein->twin(),c);
	    SMW_tgt.add_sedge_between(lateral_sv_tgt[1], ray_sv_tgt, c);
	    
	    lateral_sv_tgt[0]->twin() = lateral_sv_tgt[1];
	    lateral_sv_tgt[1]->twin() = lateral_sv_tgt[0];
	    return;
	  }
	}

	std::cerr << e->source()->point() << std::endl;
	Vertex_handle v = C.create_from_edge(e,ip);
	pl->add_vertex(v);
	std::cerr << "new vertex " << v->point() << std::endl;

	std::cerr << "double coords" << CGAL::to_double(v->point().x())
		  << ", " << CGAL::to_double(v->point().y())
		  << ", " << CGAL::to_double(v->point().z()) << std::endl;

	SVertex_iterator svi = v->svertices_begin();
	SVertex_handle svf = svi;
	SVertex_handle svb = ++svi;

	if(svf->point() == e->point()) {
	  svb->twin() = e;
	  svf->twin() = e->twin();
	  e->twin()->twin() = svf;
	  e->twin() = svb;
	} else {
	  svf->twin() = e;
	  svb->twin() = e->twin();
	  e->twin()->twin() = svb;
	  e->twin() = svf;
	}

	pl->add_edge(svf);
	pl->add_edge(svb);

	SM_walls SMW_tgt(&*v);
	
	SVertex_handle ray_sv_src = lateral_sv_tgt[0];
	SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN - lateral_sv_tgt[0]->point()));
	ray_sv_src->twin() = ray_sv_tgt;
	ray_sv_tgt->twin() = ray_sv_src;

	SVertex_handle lateral_svertex = lateral_sv_tgt[0];
	Sphere_segment sphere_ray(Sphere_point(CGAL::ORIGIN - lateral_svertex->point()), 
				  lateral_svertex->point(),c);
	lateral_sv_tgt[0] = SMW_tgt.add_lateral_svertex(sphere_ray, true, Sphere_point(dir));

	SMW_tgt.add_sedge_between(ray_sv_tgt, lateral_sv_tgt[0],c);
	
      } else if(assign(v,o2)) {
	std::cerr << "Found vertex " << v->point() << std::endl;
	
	std::cerr << "double coords" << CGAL::to_double(v->point().x())
		  << ", " << CGAL::to_double(v->point().y())
		  << ", " << CGAL::to_double(v->point().z()) << std::endl;	

	if(v == origin[1]) {
	  SM_walls SMW_tgt(&*v);
	  lateral_sv_tgt[1] = 
	    SMW_tgt.add_ray_svertex(CGAL::ORIGIN + (lateral_sv_tgt[0]->source()->point()-v->point()));

	  SMW_tgt.add_sedge_between(lateral_sv_tgt[1], ein->twin(), c);	  

	  lateral_sv_tgt[0]->twin() = lateral_sv_tgt[1];
	  lateral_sv_tgt[1]->twin() = lateral_sv_tgt[0];
	  return;
	}

	if(check_goal) {
	  Segment_3 seg(lateral_sv_tgt[0]->source()->point(),v->point());
	  if(seg.has_on(goal)) {
	    
	    if(assign(f,goal_object)) {
	      std::cerr << "Found goal on facet " << std::endl;
	      
	      Vertex_handle v = opposite[1] = C.create_from_facet(f,goal);
	      pl->add_vertex(v);
	      
	    } else if(assign(e,o2)) {
	      std::cerr << "Found goal on edge " << std::endl;
	      
	      Vertex_handle v = opposite[1] = C.create_from_edge(e,goal);
	      pl->add_vertex(v);
	      
	      SVertex_iterator svi = v->svertices_begin();
	      SVertex_handle svf = svi;
	      SVertex_handle svb = ++svi;
	      
	      if(svf->point() == e->point()) {
		svb->twin() = e;
		svf->twin() = e->twin();
		e->twin()->twin() = svf;
		e->twin() = svb;
	      } else {
		svf->twin() = e;
		svb->twin() = e->twin();
		e->twin()->twin() = svb;
		e->twin() = svf;
	      }
	      
	    } else if(assign(v,o2)) {
	      
	      std::cerr << "Found goal on vertex " << v->point() << std::endl;
	      opposite[1] = v;
	    } else 
	      CGAL_assertion_msg(false, "wrong handle");
	    
	    SM_walls SMW_src(&*origin[1]);
	    SM_walls SMW_tgt(&*opposite[1]);
	    
	    SVertex_handle ray_sv_src = SMW_src.add_ray_svertex(Sphere_point(dir));
	    SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN -dir));
	    ray_sv_src->twin() = ray_sv_tgt;
	    ray_sv_tgt->twin() = ray_sv_src;
	    
	    lateral_sv_tgt[1] = SMW_tgt.add_ray_svertex(lateral_sv_tgt[0]->point().antipode()); 
	    
	    SMW_src.add_sedge_between(ray_sv_src, ein->twin(),c);
	    SMW_tgt.add_sedge_between(lateral_sv_tgt[1], ray_sv_tgt, c);
	    
	    lateral_sv_tgt[0]->twin() = lateral_sv_tgt[1];
	    lateral_sv_tgt[1]->twin() = lateral_sv_tgt[0];
	    return;
	  }
	}

	SM_walls SMW_src(&*lateral_sv_tgt[0]->source());
	SM_walls SMW_tgt(&*v);
	
	SVertex_handle opp = SMW_tgt.add_ray_svertex(lateral_sv_tgt[0]->point().antipode());	
	opp->twin() = lateral_sv_tgt[0];
	lateral_sv_tgt[0]->twin() = opp;

	SVertex_handle lateral_svertex = lateral_sv_tgt[0];
	Sphere_segment sphere_ray(lateral_svertex->point().antipode(), 
				  lateral_svertex->point(),c);
	lateral_sv_tgt[0] = SMW_tgt.add_lateral_svertex(sphere_ray, true, Sphere_point(dir));
	SMW_tgt.add_sedge_between(opp, lateral_sv_tgt[0],c);	
      } else
	CGAL_assertion_msg(false, "wrong handle");
      
    } while(true);
    */    
    //    CGAL::SNC_io_parser<SNC_structure> O(std::cerr, *sncp, false);
    //    O.print();

  }
};

CGAL_END_NAMESPACE
#endif //CGAL_NEF3_SINGLE_WALL_CREATOR_H
