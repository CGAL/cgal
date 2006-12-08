#ifndef CGAL_NEF3_SINGLE_WALL_CREATOR2_H
#define CGAL_NEF3_SINGLE_WALL_CREATOR2_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_walls.h>
#include <CGAL/Nef_S2/Normalizing.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 229
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Single_wall_creator2 : public Modifier_base<typename Nef_::SNC_and_PL> {
  
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>  SNC_intersection;
  typedef CGAL::SNC_constructor<SNC_structure>   SNC_constructor;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>         SM_decorator;  
  typedef CGAL::SM_point_locator<SM_decorator>   SM_point_locator; 
  typedef CGAL::SM_walls<Sphere_map>             SM_walls;
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
  Sphere_point spin;
  
 public:
  Single_wall_creator2(Halfedge_handle e, Sphere_point sp)
    : ein(e), spin(sp) {}

 private:
  bool need_to_create_wall() const {
    Sphere_circle c(ein->point(), spin);
    c=normalized(c);
    //    std::cerr << "need_to_create_wall " << ein->point() << ", " << spin << std::endl;
    //    std::cerr << "need_to_create_wall " << c << std::endl;
    
    SHalfedge_around_svertex_circulator svc(ein->out_sedge()), send(svc);
    CGAL_For_all(svc,send) {
      //      std::cerr << " " << svc->circle() << std::endl;
      if(normalized(svc->circle()) == c)
	return false;
    }
    return true;
  }
    
 public:
  void operator()(SNC_and_PL& sncpl) {

    if(!need_to_create_wall())
      return;

    SNC_structure* sncp(sncpl.sncp);
    SNC_point_locator* pl(sncpl.pl);
    SNC_constructor C(*sncp,pl);

    CGAL_NEF_TRACEN( "Single_wall_creator2: ein " << ein->source()->point() 
	      << "->" << ein->twin()->source()->point() );
    CGAL_NEF_TRACEN( "Single_wall_creator2: spin " << spin );

    SM_walls SMW_src(&*ein->source());
    SM_walls SMW_tgt(&*ein->twin()->source());
    Sphere_segment sphere_ray_src(ein->point(), spin);  
    Sphere_segment sphere_ray_tgt(ein->twin()->point(), spin);    
    SVertex_handle lateral_sv_tgt[2];
    lateral_sv_tgt[0] = SMW_src.add_lateral_svertex(sphere_ray_src);
    lateral_sv_tgt[1] = SMW_tgt.add_lateral_svertex(sphere_ray_tgt);
    
    CGAL_assertion(sphere_ray_src.sphere_circle() == sphere_ray_tgt.sphere_circle().opposite());

    SMW_src.add_sedge_between(ein, lateral_sv_tgt[0], sphere_ray_src.sphere_circle());
    SMW_tgt.add_sedge_between(ein->twin(), lateral_sv_tgt[1], sphere_ray_tgt.sphere_circle());

    Sphere_circle c(sphere_ray_src.sphere_circle());
    CGAL_assertion(c.a()==0 && c.b()==0);
    Ray_hit rh(sncp, pl);

    do {
      Ray_3 r(lateral_sv_tgt[0]->source()->point(),lateral_sv_tgt[0]->point()-CGAL::ORIGIN);

      CGAL_NEF_TRACEN( "nachbehandlung: ray " << r );
      CGAL_NEF_TRACEN( "double coords" << CGAL::to_double(r.source().x())
		<< ", " << CGAL::to_double(r.source().y())
		<< ", " << CGAL::to_double(r.source().z())
		<< "->" << CGAL::to_double(r.to_vector().x()) 
		<< ", " << CGAL::to_double(r.to_vector().y()) 
		<< ", " << CGAL::to_double(r.to_vector().z()) );

      Vertex_handle v = rh.create_vertex_on_first_hit(r);

      if(v == lateral_sv_tgt[1]->source()) {
	lateral_sv_tgt[0]->twin() = lateral_sv_tgt[1];
	lateral_sv_tgt[1]->twin() = lateral_sv_tgt[0];
	return;
      }

 	SM_walls SMW_tgt(&*v);	
	SVertex_handle opp = SMW_tgt.add_ray_svertex(lateral_sv_tgt[0]->point().antipode());
	opp->twin() = lateral_sv_tgt[0];
	lateral_sv_tgt[0]->twin() = opp;

	lateral_sv_tgt[0] = 
	  SMW_tgt.add_lateral_svertex(Sphere_segment(lateral_sv_tgt[0]->point().antipode(), 
						     lateral_sv_tgt[0]->point(),c));
	SMW_tgt.add_sedge_between(opp, lateral_sv_tgt[0], c);	

    } while(true);
  }
};

CGAL_END_NAMESPACE
#endif //CGAL_NEF3_SINGLE_WALL_CREATOR2_H
