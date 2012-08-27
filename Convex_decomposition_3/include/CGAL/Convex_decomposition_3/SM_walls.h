// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$ 
// $Id$
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_SM_WALLS_H
#define CGAL_CD3_SM_WALLS_H

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 227
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename SMap>
class SM_walls : SM_decorator<SMap> {

  typedef SMap                            Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>  Base;
  typedef Base                            SM_decorator;
  typedef CGAL::SM_point_locator<Base>    SM_point_locator;
  typedef typename Base::Sphere_point     Sphere_point;
  typedef typename Base::Sphere_circle    Sphere_circle;
  typedef typename Base::Sphere_direction Sphere_direction;
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


  typedef typename Base::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator ;
  typedef typename Base::SHalfedge_around_sface_circulator 
    SHalfedge_around_sface_circulator ;
  typedef typename Base::SFace_cycle_iterator
    SFace_cycle_iterator;
  
  using Base::new_svertex;
  using Base::link_as_face_cycle;
  using Base::unlink_as_face_cycle;
  using Base::link_as_isolated_vertex;
  using Base::unlink_as_isolated_vertex;
  using Base::new_shalfedge_pair;
  using Base::unlink_as_loop;
  using Base::is_isolated;
  using Base::is_sm_boundary_object;
  using Base::delete_face;

 public:
  SM_walls(Sphere_map* M) : Base(M) {
//   SM_decorator SD(sphere_map());
//   SM_io_parser<SM_decorator>::dump(SD,std::cerr);
  }

  SHalfedge_handle find_cap(SVertex_handle sv, Sphere_point sp, Sphere_circle c) {

    CGAL_NEF_TRACEN( "find_cap " << sv->source()->point() << ":" << sv->point() 
	      << " , sp : " << sp );
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
    CGAL_NEF_TRACEN( "sv    " << sv->point() );
    CGAL_NEF_TRACEN( "c     " << c.orthogonal_vector() );
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

      CGAL_NEF_TRACEN( "sh     " << sh->circle().orthogonal_vector() );
      CGAL_NEF_TRACEN( "shnext " << shnext->circle().orthogonal_vector() );
      Sphere_segment seg(sh->circle().orthogonal_vector(), shnext->circle().orthogonal_vector(),cc);
      CGAL_NEF_TRACEN( "seg " << seg );
      if(seg.has_on(c.orthogonal_vector()))
	return sh;
    }
    CGAL_error_msg( "should not be executed");
    return SHalfedge_handle();
  }

  void insert_new_svertex_into_sedge(SVertex_handle sv, SHalfedge_handle se) {
    
    CGAL_NEF_TRACEN( "insert new svertex into sedge " << sv->point() << " | " 
	      << se->source()->point() << "->" << se->twin()->source()->point() );

    CGAL_NEF_TRACEN( "double coords " << CGAL::to_double(sv->point().x())
	      << ", " << CGAL::to_double(sv->point().y())
	      << ", " << CGAL::to_double(sv->point().z()) );

    //    SM_decorator SD(sphere_map());
    //    SM_io_parser<SM_decorator>::dump(SD,std::cerr);

    CGAL_assertion(se->circle().has_on(sv->point()));

    SHalfedge_handle se_new = this->new_shalfedge_pair();
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

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    se_new->set_index(se->get_index());
    se_opp->set_index(se->twin()->get_index());
#endif

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

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    se->set_index(sl->get_index());
    se->twin()->set_index(sl->twin()->get_index());
#endif

    unlink_as_loop(sl);
    unlink_as_loop(sl->twin());

    link_as_face_cycle(se,se->incident_sface());
    link_as_face_cycle(se->twin(),se->twin()->incident_sface());

    this->delete_loop_only();
  }

  bool legal_direction(Sphere_segment seg, Object_handle& o, Sphere_point& ip) {

    CGAL_NEF_TRACEN( "legal_direction " << seg );
    SM_point_locator P(this->sphere_map());
    o = P.ray_shoot(seg, ip, false, false);

    SVertex_handle sv;
    if(assign(sv, o)) {
      CGAL_NEF_TRACEN( "  found svertex " << sv->point() );
      return true;
    }

    SHalfedge_handle se;
    if(assign(se, o)) {
      CGAL_NEF_TRACEN( "  found sedge " << ip );
      return true;
    }

    SHalfloop_handle sl;
    if(assign(sl, o)) {
      CGAL_NEF_TRACEN( "  found sloop " << ip );
      return true;
    }

    SFace_handle sf;
    if(assign(sf, o))
      CGAL_error_msg( "wrong handle");

    CGAL_NEF_TRACEN("did not find anything");

    ip = seg.target();
    o = P.locate(seg.target());
    /*
    if(assign(sv, o)) {
      CGAL_NEF_TRACEN( "  found svertex " << sv->point() );
      if(is_isolated(sv))
	return sv->incident_sface()->mark();
      else {
	bool collinear;
	Sphere_direction d(seg.sphere_circle().opposite());
	SHalfedge_handle e_res = P.out_wedge(sv, d, collinear);
	if(collinear) {
	  CGAL_NEF_TRACEN(" collinear ");
	  o = Object_handle(e_res);
	  return false;
	} else {
          if ( e_res->circle().has_on_negative_side(seg.source()) )
            e_res = e_res->sprev();
	  o = Object_handle(e_res->incident_sface());
	  CGAL_NEF_TRACEN("  sface " << e_res->incident_sface()->mark());
	  return e_res->incident_sface()->mark();
	}
      }
    }
    */
    if(assign(sf, o)) {
      CGAL_NEF_TRACEN( "  found sface" );
      return sf->mark();
    } 

    /*    
    if(assign(se,o))
      CGAL_error_msg("wrong handle");

    if(assign(sl,o))
      CGAL_error_msg("wrong handle");
    */
    //    CGAL_NEF_SETDTHREAD(1);
    return true;
  }

  bool need_to_shoot(Sphere_point sp, SVertex_handle& sv) {
    //    CGAL_NEF_SETDTHREAD(47);
    SM_point_locator pl(this->sphere_map());
    Object_handle o = pl.locate(sp);
    //        CGAL_NEF_SETDTHREAD(1);
	
    if(assign(sv, o))
      return false;
    
    SHalfedge_handle se;
    if(assign(se, o)) {
      sv = new_svertex(sp);
      sv->mark() = se->mark();
      insert_new_svertex_into_sedge(sv, se);
      return true;
    }
    
    SFace_handle sf;
    if(assign(sf, o)) {
      if(sf->mark() == false)
	return false;
      sv = new_svertex(sp);
      sv->mark() = sf->mark();
      link_as_isolated_vertex(sv, sf);
      return true;
    }
    
    SHalfloop_handle sl;
    if(assign(sl, o)) {
      sv = new_svertex(sp);
      sv->mark() = sl->mark();
      insert_new_svertex_into_sloop(sv, sl);
      return true;
    }
    
    CGAL_error_msg( "wrong handle");
    return false;
  }
  
  SVertex_handle add_ray_svertex(Sphere_point sp) {

    CGAL_NEF_TRACEN( "add_ray_svertex " << sp );

    SM_point_locator P(this->sphere_map());

    //    SM_decorator SD(this->sphere_map());
    //    SM_io_parser<SM_decorator>::dump(SD,std::cerr);

    Object_handle o = P.locate(sp);

    return add_svertex_into_object(sp,o);
  }
  
  SVertex_handle add_svertex_into_object(Sphere_point sp, Object_handle o) {

    CGAL_NEF_TRACEN( "add_svertex_into_object " << sp );

    //    SM_decorator SD(this->sphere_map());
    //    SM_io_parser<SM_decorator>::dump(SD,std::cerr);

    SVertex_handle sv;
    SFace_handle sf;
    if(assign(sf, o)) {
      CGAL_NEF_TRACEN( "  found sface with mark " << sf->mark());
      sv = new_svertex(sp);
      sv->mark() = sf->mark();
      sv->incident_sface() = sf;
      link_as_isolated_vertex(sv,sf);
      return sv;
    }
    
    if(assign(sv, o)) {
      CGAL_NEF_TRACEN( "  found svertex" );
      return sv;
    }

    SHalfedge_handle se;
    if(assign(se, o)) {
      CGAL_NEF_TRACEN( "  found sedge");
      sv = new_svertex(sp);
      sv->mark() = se->mark();
      insert_new_svertex_into_sedge(sv, se);
      return sv;
    }
    
    SHalfloop_handle sl;
    if(assign(sl, o)) {
      CGAL_NEF_TRACEN( " found sloop" );
      sv = new_svertex(sp);
      sv->mark() = sl->mark();
      insert_new_svertex_into_sloop(sv,sl);
      return sv;
    }

    CGAL_error_msg( "wrong handle");
    return SVertex_handle();
  }

  SVertex_handle add_lateral_svertex(Sphere_segment sphere_ray, 
				     bool compare_to_dir = false, 
				     const Sphere_point& dir = Sphere_point()) {

    CGAL_NEF_TRACEN( "add_lateral_svertex " << sphere_ray );

    //    CGAL_assertion(sphere_ray.source() != dir);

    Sphere_point sp1(sphere_ray.source());
    Sphere_point sp2(sphere_ray.target());
    CGAL_NEF_TRACEN( "double coords " << CGAL::to_double(sp1.x())
	      << ", " << CGAL::to_double(sp1.y())
	      << ", " << CGAL::to_double(sp1.z())
	      << "->" << CGAL::to_double(sp2.x())
	      << ", " << CGAL::to_double(sp2.y())
	      << ", " << CGAL::to_double(sp2.z()) );

    Sphere_point ip;
    SM_point_locator P(this->sphere_map());
    Object_handle o = P.ray_shoot(sphere_ray.source(), sphere_ray.sphere_circle(), ip);
    if(compare_to_dir && dir != sphere_ray.source() && dir != ip) {
      Sphere_segment test_seg(sphere_ray.source(), ip, sphere_ray.sphere_circle());
      if(test_seg.has_on(dir)) {
	SFace_handle sf;
	o = P.locate(dir);
	CGAL_assertion(assign(sf,o));
	SVertex_handle sv = new_svertex(Sphere_point(dir));
	sv->mark() = sf->mark();
	link_as_isolated_vertex(sv, sf);
	return sv;
      }
    }

    SHalfedge_handle se;
    if(assign(se,o)) {
      CGAL_NEF_TRACEN( "  split sedge" );

      SVertex_handle sv = new_svertex(ip);
      sv->mark() = se->mark();
      insert_new_svertex_into_sedge(sv,se);
      return sv;    
    }
    
    SVertex_handle sv;
    if(assign(sv,o)) {
      CGAL_NEF_TRACEN( "  found svertex " << sv->point() );
      return sv;
    }

    SHalfloop_handle sl;
    if(assign(sl,o)) {
      CGAL_NEF_TRACEN( "  found sloop " );
      SVertex_handle sv = new_svertex(ip);
      sv->mark() = sl->mark();
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
    
    CGAL_error_msg( "wrong handle");
    return SVertex_handle();
  }

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
  void add_sedge_between(SVertex_handle sv1, SVertex_handle sv2, 
			 int& index1, int& index2,
			 Sphere_circle c) {

#else
  void add_sedge_between(SVertex_handle sv1, SVertex_handle sv2, 
			 Sphere_circle c) { // = Sphere_circle(sv1->point(),sv2->point())) {
#endif
    CGAL_NEF_TRACEN( "add sedges between " << sv1->point() 
	      << ", " << sv2->point() 
	      << " at " << sv1->source()->point() );

    bool split_sface = true;

    if(is_isolated(sv1)) {
      split_sface = false;
      if(!is_sm_boundary_object(sv1)) {
	CGAL_NEF_TRACEN( "error " << sv1->point() << "at " << sv1->source()->point() );
      }
      unlink_as_isolated_vertex(sv1);
    }
    if(is_isolated(sv2)) {
      split_sface = false;
      if(!is_sm_boundary_object(sv2)) {
	CGAL_NEF_TRACEN( "error " << sv2->point() << "at " << sv2->source()->point() );
      }
      unlink_as_isolated_vertex(sv2);
    }

    SHalfedge_handle cap1 = find_cap(sv1,sv2->point(),c);
    if(cap1 != SHalfedge_handle()) CGAL_assertion(cap1->source()==sv1);
    SHalfedge_handle cap2 = find_cap(sv2,sv1->point(),c.opposite());
    if(cap2 != SHalfedge_handle()) CGAL_assertion(cap2->source()==sv2);

    if(split_sface && 
       cap1->incident_sface() == cap2->incident_sface()) {
      SHalfedge_around_sface_circulator sfc(cap1), send(sfc);
      CGAL_For_all(sfc,send)
	if(is_sm_boundary_object(sfc))
	  unlink_as_face_cycle(sfc);
    }

    /*
     bool same_sface;
     SHalfedge_handle entry;

     if(split_sface) {      
       same_sface = 
	 cap1->incident_sface() == cap2->incident_sface();

      std::cerr << "cap1 " << cap1->source()->point()
		 << "->" << cap1->twin()->source()->point() << std::endl;
       std::cerr << "cap2 " << cap2->source()->point()
		 << "->" << cap2->twin()->source()->point() << std::endl;

       if(same_sface) {
	 SHalfedge_around_sface_circulator sfc(cap1), send(sfc);
	 CGAL_For_all(sfc,send) {
	   std::cerr << "check " << sfc->source()->point()
		     << "->" << sfc->twin()->source()->point() << std::endl;
	   if(is_sm_boundary_object(sfc))
	     entry = sfc;
	   if(sfc == cap2)
	     same_sface = false;
	 }
       }
     }
    */
    SHalfedge_handle se_new;
    if(cap1 != SHalfedge_handle()) {
      if(cap2 != SHalfedge_handle())
	se_new = new_shalfedge_pair(cap1, cap2, this->AFTER, this->AFTER);
      else 
	se_new = new_shalfedge_pair(cap1, sv2, this->AFTER);
      se_new->incident_sface() = se_new->twin()->incident_sface() = cap1->incident_sface();
    } else {
      if(cap2 != SHalfedge_handle()) {
	se_new = new_shalfedge_pair(sv1, cap2, this->AFTER);
	se_new->incident_sface() = se_new->twin()->incident_sface() = cap2->incident_sface();
      } else {
	se_new = new_shalfedge_pair(sv1, sv2);
	se_new->incident_sface() = se_new->twin()->incident_sface() = sv1->incident_sface();
      }
    }
    
    se_new->mark() = se_new->twin()->mark() = se_new->incident_sface()->mark();
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    CGAL_assertion(index1==0 || index1!=index2);
    if(index1==0) {
      se_new->set_index();
      se_new->twin()->set_index();
      index1 = se_new->get_index();
      index2 = se_new->twin()->get_index();
    } else { 
      se_new->set_index(index1);
      se_new->twin()->set_index(index2);
    }
#endif 
    CGAL_NEF_TRACEN( sv1->point() << "->" << sv2->point() << "==" 
	      << se_new->source()->point() << "->" << se_new->twin()->source()->point() );

    CGAL_assertion(sv1 == se_new->source() && 
		   sv2 == se_new->twin()->source());

    se_new->circle() = c;
    se_new->twin()->circle() = c.opposite();

    //    SM_decorator SD(this->sphere_map());
    //    SM_io_parser<SM_decorator>::dump(SD,std::cerr);

    if(split_sface) {
      if(cap1->incident_sface() == cap2->incident_sface()) {
	SFace_handle sf_new = this->new_sface();
	SFace_handle sf_old = cap1->incident_sface();
      
	CGAL_NEF_TRACEN("sf_new->mark()=" << sf_old->mark());
	sf_new->mark() = sf_old->mark();
	CGAL_assertion(sf_old->mark());
	link_as_face_cycle(se_new, sf_new);
	link_as_face_cycle(se_new->twin(), sf_old);
      } else {
	/*
	SHalfedge_handle se = cap2;
	while(se->incident_sface() != cap1->incident_sface()) {
	  se->incident_sface() = cap1->incident_sface();
	  se=se->snext();
	}
	*/
	SFace_handle sf1 = cap1->incident_sface();
	delete_face(cap2->incident_sface());
	// TODO: some relinkings are redundant
	SHalfedge_around_sface_circulator hfc(cap1), hend(hfc);
	CGAL_For_all(hfc,hend) hfc->incident_sface() = sf1;	
      }
    }

    //    SM_decorator SD1(this->sphere_map());
    //    SM_io_parser<SM_decorator>::dump(SD1,std::cerr);

    // TODO: handle inner face cycles
  }    
};

} //namespace CGAL
#endif //CGAL_CD3_SM_WALLS_H
