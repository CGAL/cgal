// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_SINGLE_WALL_CREATOR_H
#define CGAL_CD3_SINGLE_WALL_CREATOR_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Convex_decomposition_3/SM_walls.h>
#include <CGAL/Convex_decomposition_3/Ray_hit_generator.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 229
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Single_wall_creator : public Modifier_base<typename Nef_::SNC_and_PL> {
  
  typedef Nef_                                    Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL     SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure  SNC_structure;
  typedef typename SNC_structure::Items           Items;
  typedef CGAL::SNC_decorator<SNC_structure>      Base;
  typedef CGAL::SNC_point_locator<Base>           SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>   SNC_intersection;
  typedef CGAL::SNC_constructor<Items, SNC_structure>
    SNC_constructor;

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
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::SVertex_iterator        SVertex_iterator;

  typedef typename Base::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;

  SVertex_handle ein;
  Vector_3 dir;
  SNC_structure* sncp;
  SNC_point_locator* pl;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
  int index1, index2;
#endif
  
 public:
  Single_wall_creator(SVertex_handle e, Vector_3 d)
    : ein(e), dir(d) 
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    , index1(0), index2(0)
#endif
{}

 public:
  bool need_to_create_wall() const {

    CGAL_assertion(Sphere_point(CGAL::ORIGIN - dir) != ein->point());
    CGAL_assertion(Sphere_point(CGAL::ORIGIN - dir) != ein->twin()->point());

    
    Vertex_handle origin[2];
    origin[0] = ein->source();
    origin[1] = ein->target();
    SVertex_handle lateral_sv_tgt[2];
    lateral_sv_tgt[0] = ein;
    lateral_sv_tgt[1] = ein->twin();

    // isolated ein not supported

    bool legal[2];
    Object_handle found_object[2];
    Sphere_point found_point[2];
    for(int i=0; i<2; ++i) {
      CGAL_NEF_TRACEN( "SM_walls " << origin[i]->point() );
      SM_walls SMW(&*origin[i]);
      SM_point_locator PL(&*origin[i]);
      Sphere_segment sphere_ray(lateral_sv_tgt[i]->point(), Sphere_point(dir));
      legal[i] = SMW.legal_direction(sphere_ray, found_object[i], found_point[i]);
    }

    SVertex_handle sv0, sv1;
    if(assign(sv0, found_object[0]) && assign(sv1, found_object[1])) {
      CGAL_NEF_TRACEN( "check " << sv0->point() << "+" << sv1->point() );
      SHalfedge_around_svertex_circulator sh0(sv0->out_sedge()), send0(sh0);
      CGAL_For_all(sh0,send0)
	if(sh0->twin()->source() == ein &&
	   Sphere_segment(sh0->source()->point(),
			  sh0->twin()->source()->point(),
			  sh0->circle()).is_short()) {
	  SHalfedge_around_svertex_circulator sh1(sv1->out_sedge()), send1(sh1);
	  CGAL_For_all(sh1,send1)
	    if(sh1->twin()->source() == ein->twin() &&
	       Sphere_segment(sh1->source()->point(),
			      sh1->twin()->source()->point(),
			      sh1->circle()).is_short()) {
	      CGAL_NEF_TRACEN( "did not process edge " );
	      CGAL_NEF_TRACEN( "check " << sh0->source()->point() << 
			       "->" << sh0->twin()->source()->point() );
	      CGAL_NEF_TRACEN( "check " << sh1->source()->point() << 
			       "->" << sh1->twin()->source()->point() );
	      return false;
	    }
	  CGAL_error_msg( "should not happen on one side only");
	}
    }

    // TODO: check if cycle already exists

    if(!legal[0] || !legal[1]) 
      return false;    

    return true;
  }
    
  SVertex_handle create_new_outer_cycle(SVertex_handle estart, Sphere_circle c) {

    CGAL_NEF_TRACEN( "create new outer cycle " << 
		     estart->source()->point() << " to " <<
		     estart->twin()->source()->point());
    CGAL_NEF_TRACEN( "double coords" << CGAL::to_double(estart->source()->point().x())
		     << ", " << CGAL::to_double(estart->source()->point().y())
		     << ", " << CGAL::to_double(estart->source()->point().z()) );
    CGAL_NEF_TRACEN( "double coords" << CGAL::to_double(estart->twin()->source()->point().x())
		     << ", " << CGAL::to_double(estart->twin()->source()->point().y())
		     << ", " << CGAL::to_double(estart->twin()->source()->point().z()) );

    SM_walls SMW(&*estart->source());
    Sphere_segment sphere_ray(estart->point(), estart->twin()->point(), c);
    SVertex_handle lateral_svertex = 
      SMW.add_lateral_svertex(sphere_ray);

    /*
    Sphere_segment test_seg(estart->point(), lateral_svertex->point(), c);
    if(lateral_svertex->point() != Sphere_point(1,0,0)) {
      //      std::cerr << lateral_svertex->point() << std::endl;
      CGAL_assertion(!test_seg.has_on(Sphere_point(1,0,0)));
    }
    if(lateral_svertex->point() != Sphere_point(-1,0,0)) {
      //      std::cerr << lateral_svertex->point() << std::endl;
      CGAL_assertion(!test_seg.has_on(Sphere_point(-1,0,0)));
    }
    */

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      SMW.add_sedge_between(estart, lateral_svertex, index1, index2, c);
#else
      SMW.add_sedge_between(estart, lateral_svertex, c);
#endif

    /* ???
    CGAL_assertion_code
      (Sphere_segment test2(estart->point(), 
			    lateral_svertex->point(), c));
    CGAL_assertion(test2.has_on(Sphere_point(dir)));
    CGAL_assertion(test2.has_on(Sphere_point(-dir)));
    */

    Ray_hit rh(sncp, pl, 3);
    Ray_3 r(lateral_svertex->source()->point(), 
	    lateral_svertex->point()-CGAL::ORIGIN);
    Vertex_handle v = rh.create_vertex_on_first_hit(r);

    while(v != estart->twin()->source()) {
      
      CGAL_NEF_TRACEN( "current vertex " << v->point() );

      CGAL_NEF_TRACEN( "double coords" << CGAL::to_double(v->point().x())
		       << ", " << CGAL::to_double(v->point().y())
		       << ", " << CGAL::to_double(v->point().z()) );

      CGAL_NEF_TRACEN( "SM_walls " << v->point() );
      SM_walls smw(&*v);
      SVertex_handle opp = smw.add_ray_svertex(lateral_svertex->point().antipode());
      SM_walls smwt(&*v);
      opp->twin() = lateral_svertex;
      lateral_svertex->twin() = opp;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      opp->set_index();
      lateral_svertex->set_index(opp->get_index());
#endif
      pl->add_edge(lateral_svertex);

      CGAL_NEF_TRACEN( "twins " << lateral_svertex->point() 
		       << " + " << opp->point() );

      sphere_ray = Sphere_segment(lateral_svertex->point().antipode(), 
				  lateral_svertex->point(), c);
      lateral_svertex = smw.add_lateral_svertex(sphere_ray);

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      smw.add_sedge_between(opp, lateral_svertex, index1, index2, c);
#else
      smw.add_sedge_between(opp, lateral_svertex, c);
#endif

      // TODO: make use of existing edges along ray
      r = Ray_3(lateral_svertex->source()->point(), lateral_svertex->point()-CGAL::ORIGIN);
      v = rh.create_vertex_on_first_hit(r);
    }
    
    CGAL_NEF_TRACEN( "last current vertex " << v->point() );

    CGAL_NEF_TRACEN( "SM_walls " << v->point() );
    SM_walls smw(&*v);
    SVertex_handle opp = 
      smw.add_ray_svertex(lateral_svertex->point().antipode());
    opp->twin() = lateral_svertex;
    lateral_svertex->twin() = opp;

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      opp->set_index();
      lateral_svertex->set_index(opp->get_index());
#endif
    pl->add_edge(lateral_svertex);

#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    smw.add_sedge_between(opp, estart->twin(), index1, index2, c);
#else
    smw.add_sedge_between(opp, estart->twin(), c);
#endif

    CGAL_NEF_TRACEN( "final twins " << lateral_svertex->source()->point() 
		     << " + " << opp->source()->point() );

    return lateral_svertex;
  }

  void insert_into_outer_cycle(SVertex_handle estart, Sphere_circle c) {

    CGAL_assertion(false);

    Ray_hit rh(sncp, pl, 3);
    SVertex_handle lateral_svertex = estart;
    Vertex_handle v = estart->twin()->source();

    do {
      CGAL_NEF_TRACEN( "current vertex " << v->point() );

      CGAL_NEF_TRACEN( "double coords" << CGAL::to_double(v->point().x())
		       << ", " << CGAL::to_double(v->point().y())
		       << ", " << CGAL::to_double(v->point().z()) );

      CGAL_NEF_TRACEN( "SM_walls " << v->point() );
      SM_walls smw(&*v);
      SVertex_handle opp = smw.add_ray_svertex(lateral_svertex->point().antipode());
      opp->twin() = lateral_svertex;
      lateral_svertex->twin() = opp;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      opp->set_index();
      lateral_svertex->set_index(opp->get_index());
#endif
      pl->add_edge(lateral_svertex);

      CGAL_NEF_TRACEN( "twins " << lateral_svertex->source()->point() 
		       << " + " << opp->source()->point() );

      Sphere_segment sphere_ray = Sphere_segment(lateral_svertex->point().antipode(), 
				  lateral_svertex->point(), c);
      lateral_svertex = smw.add_lateral_svertex(sphere_ray);
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    smw.add_sedge_between(opp, lateral_svertex, index1, index2, c);
#else
    smw.add_sedge_between(opp, lateral_svertex, c);
#endif

      // TODO: make use of existing edges along ray

      Ray_3 r = Ray_3(lateral_svertex->source()->point(), lateral_svertex->point()-CGAL::ORIGIN);
      //      CGAL_NEF_SETDTHREAD(503*509);
      v = rh.create_vertex_on_first_hit(r);
      //      CGAL_NEF_SETDTHREAD(1);
    } while(v != estart->source());
    
    CGAL_NEF_TRACEN( "last current vertex " << v->point() );

    CGAL_NEF_TRACEN( "SM_walls " << v->point() );
    SM_walls smw(&*v);
    SVertex_handle opp = smw.add_ray_svertex(lateral_svertex->point().antipode());    
    opp->twin() = lateral_svertex;
    lateral_svertex->twin() = opp;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      opp->set_index();
      lateral_svertex->set_index(opp->get_index());
#endif
    pl->add_edge(lateral_svertex);

    //    smw.add_sedge_between(opp, estart->twin(), c);

    //    CGAL_NEF_TRACEN( "final twins " << lateral_svertex->source()->point() 
    //		     << " + " << opp->source()->point() );
  }

 public:
  void operator()(SNC_and_PL& sncpl) {

    //    CGAL_NEF_SETDTHREAD(47*227*229*233);

    if(!need_to_create_wall()) return;

//    if(!Reflex_edge_searcher<Nef_polyhedron>::is_reflex_edge(ein))
//        return;

    sncp = sncpl.sncp;
    pl = sncpl.pl;

    //    SNC_io_parser<SNC_structure> O0(std::cerr,*sncp);
    //    O0.print();

    //    CGAL_NEF_SETDTHREAD(229*227);

    SVertex_handle target_svertex = ein->twin();
    Sphere_circle c(ein->point(), Sphere_point(dir));

    c = normalized(c);
    do {
      ein = target_svertex->twin(); // for subsequent runs of the loop
      SVertex_handle svopen = 
	create_new_outer_cycle(ein, c);

      if(ein->twin() != target_svertex) {
	// TODO: what indexes are needed here?
	SHalfedge_handle seopen = svopen->out_sedge();
	while(seopen->circle() == c || seopen->circle() == c.opposite())
	  seopen = seopen->sprev()->twin();
	CGAL_assertion(seopen->circle() != c &&
		       seopen->circle() != c.opposite());
	CGAL_NEF_TRACEN("found open sedge " << seopen->source()->point()
			<< "->" << seopen->twin()->source()->point()
			<< ":" << seopen->circle());
	// TODO: what if ein is isloated
	insert_into_outer_cycle(svopen, seopen->circle().opposite());
      }
	
    } while(ein->twin() != target_svertex);

  }
};

CGAL_END_NAMESPACE
#endif //CGAL_CD3_SINGLE_WALL_CREATOR_H
