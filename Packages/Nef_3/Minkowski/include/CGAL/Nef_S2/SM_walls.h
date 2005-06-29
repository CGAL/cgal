#ifndef CGAL_SM_WALLS_H
#define CGAL_SM_WALLS_H

CGAL_BEGIN_NAMESPACE

template<typename SMap>
class SM_walls : SM_decorator<SMap> {

  typedef SMap                            Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>  Base;
  typedef Base                            SM_decorator;
  typedef CGAL::SM_point_locator<Base>    SM_point_locator;
  typedef typename Base::Sphere_point     Sphere_point;
  typedef typename Base::Sphere_circle    Sphere_circle;
  typedef typename Base::Sphere_segment   Sphere_segment;
  
  typedef typename Base::SVertex_handle   SVertex_handle;
  typedef typename Base::SFace_handle     SFace_handle;
  typedef typename Base::SHalfedge_handle SHalfedge_handle;
  typedef typename Base::SHalfloop_handle SHalfloop_handle;
  typedef typename Base::Object_handle    Object_handle;

  typedef typename Sphere_point::Vector_3       Vector_3;
  typedef typename Sphere_circle::Point_3        Point_3;
  typedef typename Sphere_circle::Plane_3        Plane_3;
  typedef typename Sphere_point::FT             FT;


  typedef typename SMap::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator ;
  
 public:
  SM_walls(Sphere_map* M) : Base(M) {}

  SHalfedge_handle find_cap(SVertex_handle sv, Sphere_point sp, Sphere_circle c) {

    std::cerr << "find_cap " << sv->source()->point() << ":" << sv->point() 
	      << " , sp : " << sp << std::endl;
    /*
    SHalfedge_handle se = sv->out_sedge();
    if(se != SHalfedge_handle())
      while(CGAL::spherical_orientation(cas(se)->twin()->source()->point(),
					sp,
					se->twin()->source()->point()) > 0)
	se = cas(se);    

    if(se != SHalfedge_handle()) {
      CGAL_assertion(Sphere_circle(sv->point(),sp) != 
		     Sphere_circle(sv->point(), se->twin()->source()->point()));
    }
    */
    
    SM_decorator SD(&*sv->source());
    if( SD.is_isolated(sv))
      return SHalfedge_handle();    

    //    Sphere_circle c(sv->point(), sp);
    std::cerr << "sv    " << sv->point() << std::endl;
    std::cerr << "c     " << c.orthogonal_vector() << std::endl;
    SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);

    Plane_3 pl(Point_3(0,0,0),Vector_3(sv->point()-CGAL::ORIGIN));
    Sphere_circle cc(pl);

    SHalfedge_around_svertex_circulator shnext(sh);      
    ++shnext;
    if(sh == shnext)
      return sh;

    CGAL_For_all(sh,send) {
      shnext =sh;
      ++shnext;

      std::cerr << "sh     " << sh->circle().orthogonal_vector() << std::endl;
      std::cerr << "shnext " << shnext->circle().orthogonal_vector() << std::endl;
      Sphere_segment seg(sh->circle().orthogonal_vector(), shnext->circle().orthogonal_vector(),cc);
      std::cerr << "seg " << seg << std::endl;
      if(seg.has_on(c.orthogonal_vector()))
	return sh;
    }
    CGAL_assertion_msg(false, "should not be executed");
    return SHalfedge_handle();
  }

  /*
  void insert_wall(Sphere_point sp, SVertex_handle sv) {
    SHalfedge_handle se = find_cap(sv, sp);
    SVertex_handle sv_new = new_svertex(sp);
    SHalfedge_handle se_new = new_shalfedge_pair(se, sv_new, AFTER);
    se_new->incident_sface() = se_new->twin()->incident_sface() = se->incident_sface();
    se_new->circle() = Sphere_circle(se_new->source()->point(), se_new->twin()->source()->point());
    se_new->twin()->circle() = se_new->circle().opposite();
  }
  */

  void insert_new_svertex_into_sedge(SVertex_handle sv, SHalfedge_handle se) {

    std::cerr << "insert new svertex into sedge " << sv->point() << " | " 
	      << se->source()->point() << "->" << se->twin()->source()->point() << std::endl;

    std::cerr << "double coords " << CGAL::to_double(sv->point().x())
	      << ", " << CGAL::to_double(sv->point().y())
	      << ", " << CGAL::to_double(sv->point().z()) << std::endl;

    CGAL_assertion(se->circle().has_on(sv->point()));

    SHalfedge_handle se_new = new_shalfedge_pair();
    SHalfedge_handle se_opp = se_new->twin();
    se_new->source() = sv;
    se_opp->source() = se->twin()->source();

    se_new->circle() = se->circle();
    CGAL_assertion(se_new->circle().has_on(se_new->source()->point()) &&
		   se_new->circle().has_on(se_opp->source()->point()));
    se_opp->circle() = se->twin()->circle();
    se->twin()->source() = sv;
    se_new->mark() = se_opp->mark() = se->mark();
    se_new->incident_sface() = se->incident_sface();
    se_opp->incident_sface() = se->twin()->incident_sface();

    se_new->snext() = se->snext();
    se->snext()->sprev() = se_new;
    se->snext() = se_new;
    se_new->sprev() = se;

    se_opp->sprev() = se->twin()->sprev();
    se->twin()->sprev()->snext() = se_opp;
    se->twin()->sprev() = se_opp;
    se_opp->snext() = se->twin();

    se_new->source()->out_sedge() = se_new;
    se_opp->source()->out_sedge() = se_opp;
  }

  void insert_new_svertex_into_sloop(SVertex_handle sv, SHalfloop_handle sl) {

    SHalfedge_handle se = new_shalfedge_pair(sv, sv);
    se->circle() = sl->circle();
    se->twin()->circle() = sl->twin()->circle();
    se->snext() = se->sprev() = se;
    se->twin()->snext() = se->twin()->sprev() = se->twin();
    se->incident_sface() = sl->incident_sface();
    se->twin()->incident_sface() = sl->twin()->incident_sface();
    se->mark() = se->twin()->mark() = sl->mark();
    store_sm_boundary_object(se,se->incident_sface());
    store_sm_boundary_object(se->twin(),se->twin()->incident_sface());

    unlink_as_loop(sl);
    unlink_as_loop(sl->twin());
    delete_loop_only();
  }

  bool legal_direction(Sphere_segment seg, Object_handle& o, Sphere_point& ip) {

    std::cerr << "legal_direction " << seg << std::endl;
    //    CGAL_NEF_SETDTHREAD(47);
    SM_point_locator P(sphere_map());
    o = P.ray_shoot(seg, ip, false, false);

    SVertex_handle sv;
    if(assign(sv, o)) {
      std::cerr << "  found svertex " << sv->point() << std::endl;
      return true;
    }

    SHalfedge_handle se;
    if(assign(se, o)) {
      std::cerr << "  found sedge " << ip << std::endl;
      return true;
    }

    SHalfloop_handle sl;
    if(assign(sl, o))
      CGAL_assertion_msg(false, "not implemented yet");

    SFace_handle sf;
    if(assign(sf, o))
      CGAL_assertion_msg(false, "wrong handle");

    ip = seg.target();
    o = P.locate(seg.target());

    if(assign(sf, o)) {
      std::cerr << "  found sface" << std::endl;
      if(sf->mark() == false)
	return false;
    }

    //    CGAL_NEF_SETDTHREAD(1);
    return true;
  }

  bool need_to_shoot(Sphere_point sp) {
    SM_point_locator pl(sphere_map());
    Object_handle o = pl.locate(sp);
    
    SVertex_handle sv;
    if(assign(sv, o))
      return false;
    
    SHalfedge_handle se;
    if(assign(se, o)) {
      sv = new_svertex(sp);
      insert_new_svertex_into_sedge(sv, se);
      return true;
    }
    
    SFace_handle sf;
    if(assign(sf, o)) {
      if(sf->mark() == false)
	return false;
      sv = new_svertex(sp);
      return true;
    }
    
    SHalfloop_handle sl;
    if(assign(sl, o))
      CGAL_assertion_msg(false, "not implemented yet");	
    
    CGAL_assertion_msg(false, "wrong handle");
    return false;
  }
  
  SVertex_handle add_ray_svertex(Sphere_point sp) {

    std::cerr << "add_ray_svertex " << sp << std::endl;

    SM_point_locator P(sphere_map());

    SM_decorator SD(sphere_map());
    SM_io_parser<SM_decorator>::dump(SD,std::cerr);

    //    CGAL_NEF_SETDTHREAD(47);
    Object_handle o = P.locate(sp);
    //    CGAL_NEF_SETDTHREAD(1);

    return add_svertex_into_object(sp,o);
  }
   
  SVertex_handle add_svertex_into_object(Sphere_point sp, Object_handle o) {

    std::cerr << "add_svertex_into_object " << sp << std::endl;

    SVertex_handle sv;
    SFace_handle sf;
    if(assign(sf, o)) {
      std::cerr << "  found sface" << std::endl;
      sv = new_svertex(sp);
      sv->incident_sface() = sf;
      return sv;
    }
    
    if(assign(sv, o)) {
      std::cerr << "  found svertex" << std::endl;
      return sv;
    }


    SHalfedge_handle se;
    if(assign(se, o)) {
      std::cerr << "  found sedge" <<std::endl;
      sv = new_svertex(sp);
      insert_new_svertex_into_sedge(sv, se);
      return sv;
    }
    
    SHalfloop_handle sl;
    if(assign(sl, o)) {
      std::cerr << " found sloop" << std::endl;
      sv = new_svertex(sp);
      insert_new_svertex_into_sloop(sv,sl);
      return sv;
    }

    CGAL_assertion_msg(false, "wrong handle");
    return SVertex_handle();
  }

  SVertex_handle add_lateral_svertex(Sphere_segment sphere_ray, 
				     bool compare_to_dir = false, 
				     const Sphere_point& dir = Sphere_point()) {

    std::cerr << "add_lateral_svertex " << sphere_ray << std::endl;

    CGAL_assertion(sphere_ray.source() != dir);

    Sphere_point sp1(sphere_ray.source());
    Sphere_point sp2(sphere_ray.target());
    std::cerr << "double coords " << CGAL::to_double(sp1.x())
	      << ", " << CGAL::to_double(sp1.y())
	      << ", " << CGAL::to_double(sp1.z())
	      << "->" << CGAL::to_double(sp2.x())
	      << ", " << CGAL::to_double(sp2.y())
	      << ", " << CGAL::to_double(sp2.z()) << std::endl;

    Sphere_point ip;
    SM_point_locator P(sphere_map());
    //    CGAL_NEF_SETDTHREAD(47);
    Object_handle o = P.ray_shoot(sphere_ray.source(), sphere_ray.sphere_circle(), ip);
    //    CGAL_NEF_SETDTHREAD(1);
    
    if(compare_to_dir) {
      Sphere_segment test_seg(sphere_ray.source(), ip, sphere_ray.sphere_circle());
      if(test_seg.has_on(dir))
	ip = dir;
    }

    SHalfedge_handle se;
    if(assign(se,o)) {
      std::cerr << "  split sedge" << std::endl;

      SVertex_handle sv = new_svertex(ip);
      insert_new_svertex_into_sedge(sv,se);
      
      return sv;    
    }
    
    SVertex_handle sv;
    if(assign(sv,o)) {
      std::cerr << "  found svertex " << sv->point() << std::endl;
      return sv;
    }

    SHalfloop_handle sl;
    if(assign(sl,o)) {
      std::cerr << "  found sloop " << std::endl;
      SVertex_handle sv = new_svertex(ip);
      
      insert_new_svertex_into_sloop(sv,sl);
      /*
      se = new_shalfedge_pair(sv, sv);
      se->circle() = sl->circle();
      se->twin()->circle() = sl->twin()->circle();
      se->snext() = se->sprev() = se;
      se->twin()->snext() = se->twin()->sprev() = se->twin();
      se->incident_sface() = sl->incident_sface();
      se->twin()->incident_sface() = sl->twin()->incident_sface();
      se->mark() = se->twin()->mark() = sl->mark();
      store_sm_boundary_object(se,se->incident_sface());
      store_sm_boundary_object(se->twin(),se->twin()->incident_sface());

      unlink_as_loop(sl);
      unlink_as_loop(sl->twin());
      delete_loop_only();
      */

      return sv;
    }
    
    CGAL_assertion_msg(false, "wrong handle");
    return SVertex_handle();
  }

  void add_sedge_between(SVertex_handle sv1, SVertex_handle sv2, 
			 Sphere_circle c = Sphere_circle()) { // = Sphere_circle(sv1->point(),sv2->point())) {
    std::cerr << "add sedges between " << sv1->point() 
	      << ", " << sv2->point() 
	      << " at " << sv1->source()->point() << std::endl;

    /*
    Sphere_circle test(Sphere_circle(sv1->point(),sv2->point()));
    CGAL_assertion(c == normalized(test) || 
		   sv1->point().antipode() == sv2->point());
    */

    if(c == Sphere_circle()) {
      c = Sphere_circle(sv1->point(), sv2->point());
      c = normalized(c);
    }

    SHalfedge_handle cap1 = find_cap(sv1,sv2->point(),c);
    if(cap1 != SHalfedge_handle()) CGAL_assertion(cap1->source()==sv1);
    SHalfedge_handle cap2 = find_cap(sv2,sv1->point(),c.opposite());
    if(cap2 != SHalfedge_handle()) CGAL_assertion(cap2->source()==sv2);

    SHalfedge_handle se_new;
    if(cap1 != SHalfedge_handle()) {
      if(cap2 != SHalfedge_handle())
	se_new = new_shalfedge_pair(cap1, cap2, AFTER, AFTER);
      else 
	se_new = new_shalfedge_pair(cap1, sv2, AFTER);
      se_new->incident_sface() = se_new->twin()->incident_sface() = cap1->incident_sface();
    } else {
      if(cap2 != SHalfedge_handle()) {
	se_new = new_shalfedge_pair(sv1, cap2, AFTER);
	se_new->incident_sface() = se_new->twin()->incident_sface() = cap2->incident_sface();
      } else {
	se_new = new_shalfedge_pair(sv1, sv2);
	se_new->incident_sface() = se_new->twin()->incident_sface() = sv1->incident_sface();
      }
    }
    
    se_new->mark() = se_new->twin()->mark() = true; // = se_new->incident_sface()->mark();

    std::cerr << sv1->point() << "->" << sv2->point() << "==" 
	      << se_new->source()->point() << "->" << se_new->twin()->source()->point() << std::endl;

    CGAL_assertion(sv1 == se_new->source() && 
		   sv2 == se_new->twin()->source());

    se_new->circle() = c;
    se_new->twin()->circle() = c.opposite();

    // does the new sedge split a facet ?
  }    
  /*
  SVertex_handle add_two(Sphere_point sp1, Sphere_point sp2) {
    std::cerr << "add_two " << sp1 << ", " << sp2 << std::endl;
    Sphere_point ip;
    SM_point_locator P(sphere_map());
    Object_handle o = P.ray_shoot(sp2,Sphere_circle(sp2,sp1),ip);

    SHalfedge_handle se;
    if(assign(se,o)) {
      //      CGAL_assertion_msg(false, "wrong handle");
      std::cerr << "split sedge" << std::endl;
  */
      /*
      SVertex_handle svt = se->source();
      typename SMap::SHalfedge_around_svertex_const_circulator sav(svt->out_sedge()), send(sav);      
      CGAL_For_all(sav,send) {
	std::cerr << sav->twin()->source()->point() << std::endl;
      }
      */
  /*
      SVertex_handle sv = new_svertex(ip);
      insert_new_svertex_into_sedge(sv,se);
  */  
      /*
      std::cerr << " " << std::endl;
      CGAL_For_all(sav,send) {
	std::cerr << sav->twin()->source()->point() << std::endl;
      }

      CGAL_For_all(sav,send) {
	typename SMap::SHalfedge_around_svertex_const_circulator snext(sav), snn(sav);
	++snext; ++snn; ++snn;
	typename SMap::SVertex_const_handle stwin1(sav->twin()->source());
	typename SMap::SVertex_const_handle stwin2(snext->twin()->source());
	typename SMap::SVertex_const_handle stwin3(snn->twin()->source());
	CGAL_assertion(spherical_orientation(stwin1->point(),
					     stwin2->point(),
					     stwin3->point()) >= 0);
      }
      */
      /*
      return sv;
    }

    SHalfloop_handle sl;
    CGAL_assertion_msg(!assign(sl,o), "not implemented yet");

    SVertex_handle sv;
    //    CGAL_assertion_msg(!assign(sv,o), "not implemented yet");
    if(assign(sv,o)) {
      insert_wall(sp2,sv);
      return sv;
    }

    CGAL_assertion_msg(false, "wrong handle");
    return SVertex_handle();
  }
  */
  /*
  SVertex_handle add_outgoing(Sphere_point sp, SVertex_handle sv) {
    SVertex_handle sv_new = new_svertex(sp);
    SHalfedge_handle se_cap = find_cap(sv,sp);
    SHalfedge_handle se_new = new_shalfedge_pair(se_cap, sv_new, AFTER);
    se_new->incident_sface() = se_new->twin()->incident_sface() = se_cap->incident_sface();
    se_new->circle() = Sphere_circle(se_new->source()->point(), se_new->twin()->source()->point());
    se_new->twin()->circle() = se_new->circle().opposite();
    return sv_new;
  }
  */
  /*
  void insert_new_svertex_into_sface(SFace_handle sf, const Sphere_point& sp) {
    SVertex_handle sv = new_svertex(sp);
    
    SFace_cycles_iterator sfc = sf.sface_cycles_begin();
    CGAL_assertion(sfc.is_shalfedge());
    SHalfedge_around_sface_circulator sec(SHalfedge_handle(sfc)), send(sec);
    CGAL_For_all(sec,send)
      new_shalfedge(sv,sec.source());
    SHalfedge_around_svertex_circulator sevc(sv.out_sedge()), eend(sevc);
    CGAL_For_all(sevc,eend) {
      SFace_handle sf_new = new_sface();
      sf_new.mark() = sf.mark();
      link_as_sface_cycle(sevc,sf_new);
    }
    delete_sface(sf);
  }
  */

  /*
  void create_opposite_vertex_on_facet(Halffacet_handle f, SVertex_handle sv_in, 
				       const Sphere_point& sp) {
    SVertex_handle sv = new_svertex(sp);
    SHalfedge_around_svertex_circulator sevc(sv_in.out_sedge()), eend(sevc);
    do { 
      Vector_3 or(f->plane()->orthogonal_vector());
      Sphere_point p(sevc.twin().source().point()),p_new;
      if(p.hx() != 0)
	p_new = Sphere_point(-(p.hy()*or.b()+p.hz()*or.c()),
			     or.a()*p.hy(), 
			     or.a()*p.hz(), 
			     or.a()*p.hw());
      else if(p.hy() !=0)
	p_new = Sphere_point(or.b()*p.hx(),
			     -(p.hx()*or.a()+p.hz()*or.c()),
			     or.b()*p.hz(), 
			     or.b()*p.hw());
      else if(p.hz() !=0)
	p_new = Sphere_point(or.c()*p.hx(),
			     or.c()*p.hy(),
			     -(p.hx()*or.a()+p.hy()*or.b()), 
			     or.c()*p.hw());
      SVertex_handle svh = new_svertex(p_new);
      svh->mark() = sevc->mark();
      --sevc;
    } while(sevc != eend);

    SVertex_iterator svi(svertices_begin());
    ++svi;
    SVertex_iterator sv_next(svi);
    ++sv_next;
    while(sv_next != svertices_end()) {
      new_shalfedge(*svi, *(++svi));
      new_shalfedge(*svi, sv); 
    }
    new_shalfedge(*svi,*(++svertices_begin()));
    new_shalfedge(*svi, sv); 
  }
  */

};

CGAL_END_NAMESPACE
#endif //CGAL_SM_WALLS_H
