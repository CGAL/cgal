#ifndef CGAL_NEF3_SINGLE_WALL_CREATOR_H
#define CGAL_NEF3_SINGLE_WALL_CREATOR_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_walls.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Single_wall_creator : public Modifier_base<typename Nef_::SNC_and_PL> {
  
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
  //typedef CGAL::Ray_shooter<SNC_structure>       Ray_shooter;
  

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
      return 3;

    std::cerr << "test 0 " << std::endl;
    
    SHalfedge_handle se;

    SM_point_locator PS(&*ein->source());
    Object_handle os = PS.locate(Sphere_point(dir));
    if(assign(se,os) && (se->source() == ein || se->twin()->source() == ein))
      return 3;
    
    std::cerr << "test 1 " << std::endl;

    SM_point_locator PT(&*ein->target());
    Object_handle ot = PT.locate(Sphere_point(dir));
    if(assign(se,ot) && (se->source() == ein->twin() || se->twin()->source() == ein->twin()))
      return 3;

    std::cerr << "test 2 " << std::endl;

    SFace_handle sfs, sft;
    if(assign(sfs,os)) {
      if(sfs->mark() == false)
	return 3;
      if(assign(sft,ot)) {
	if(sft->mark() == false)
	  return 3;
	if(sfs->volume() != sft->volume())
	  return 3;
      }
    } else if(assign(sft,ot) && sft->mark() == false)
      return 3;

    std::cerr << "test 3 " << std::endl;
    
    SVertex_handle sv;
    if(assign(sv,os)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send)
	if(sh->twin()->source() == ein) // tighther, i.e. sedge with circle(svs,ein)
	  res |= 1;
    }

    std::cerr << "test 4 " << std::endl;

    if(assign(sv,ot)) {
      SHalfedge_around_svertex_circulator sh(sv->out_sedge()), send(sh);
      CGAL_For_all(sh,send)
	if(sh->twin()->source() == ein->twin()) // tighther, i.e. sedge with circle(svs,ein)
	  res |= 2;
    }

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
    CGAL_assertion(ntcw >=0 && ntcw <4);
    if(ntcw == 3) return;

    Vertex_handle origin[2];
    Vertex_handle opposite[2];
    origin[0] = ein->source();
    origin[1] = ein->target();
    SVertex_handle lateral_sv_tgt[2];

    for(int i=0; i<2; ++i) {

      int tmp = (ntcw & (i+1));
      std::cerr << "ntcw " << ntcw << "&" << (i+1) << "=" << tmp << std::endl;

      /*
      if((ntcw & (i+1)) != 0)
	continue;
      */

      Vertex_handle v;
      Halfedge_handle e;
      Halffacet_handle f;
      Ray_3 r(origin[i]->point(),dir);
      std::cerr << "shoot ray " << r << std::endl;
      Object_handle o2 = pl->shoot(r);
      if(assign(f,o2)) {
	std::cerr << "Found facet " << std::endl;

	Point_3 ip;
	SNC_intersection I;
	I.does_intersect_internally(r, f, ip);

	Vertex_handle v = opposite[i] = C.create_from_facet(f,ip);
	
	SM_walls SMW_src(&*origin[i]);
	SM_walls SMW_tgt(&*opposite[i]);
	
	SVertex_handle ray_sv_src = SMW_src.add_ray_svertex(Sphere_point(dir));
	SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN -dir));
	ray_sv_src->twin() = ray_sv_tgt;
	ray_sv_tgt->twin() = ray_sv_src;	
	
	SVertex_handle lateral_svertex(i==0?ein:ein->twin());
	Sphere_segment sphere_ray(Sphere_point(CGAL::ORIGIN - dir), lateral_svertex->point());
	lateral_sv_tgt[i] = SMW_tgt.add_lateral_svertex(sphere_ray);

	SMW_src.add_sedge_between(ray_sv_src, lateral_svertex);
	SMW_tgt.add_sedge_between(ray_sv_tgt, lateral_sv_tgt[i]);

      } else if(assign(e,o2)) {
	std::cerr << "Found edge " << std::endl;

	Point_3 ip;
	SNC_intersection I;
	I.does_intersect_internally(r, Segment_3(e->source()->point(),
						 e->twin()->source()->point()),  
				    ip);

	Vertex_handle v = opposite[i] = C.create_from_edge(e,ip);
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
	  
	SM_walls SMW_src(&*origin[i]);
	SM_walls SMW_tgt(&*opposite[i]);
	
	SVertex_handle ray_sv_src = SMW_src.add_ray_svertex(Sphere_point(dir));
	SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN -dir));
	ray_sv_src->twin() = ray_sv_tgt;
	ray_sv_tgt->twin() = ray_sv_src;

	SVertex_handle lateral_svertex(i==0?ein:ein->twin());
	Sphere_segment sphere_ray(Sphere_point(CGAL::ORIGIN - dir), lateral_svertex->point());
	lateral_sv_tgt[i] = SMW_tgt.add_lateral_svertex(sphere_ray);

	SMW_src.add_sedge_between(ray_sv_src, lateral_svertex);
	SMW_tgt.add_sedge_between(ray_sv_tgt, lateral_sv_tgt[i]);

      } else if(assign(v,o2)) {
	
	std::cerr << "Found vertex " << v->point() << std::endl;
	opposite[i] = v;
	SM_walls SMW_src(&*origin[i]);
	SM_walls SMW_tgt(&*opposite[i]);
	
	SVertex_handle ray_sv_src = SMW_src.add_ray_svertex(Sphere_point(dir));
	SVertex_handle ray_sv_tgt = SMW_tgt.add_ray_svertex(Sphere_point(CGAL::ORIGIN -dir));
	ray_sv_src->twin() = ray_sv_tgt;
	ray_sv_tgt->twin() = ray_sv_src;
	
	SVertex_handle lateral_svertex(i==0?ein:ein->twin());
	Sphere_segment sphere_ray(Sphere_point(CGAL::ORIGIN - dir), lateral_svertex->point());
	lateral_sv_tgt[i] = SMW_tgt.add_lateral_svertex(sphere_ray);
	
	SMW_src.add_sedge_between(ray_sv_src, lateral_svertex);
	SMW_tgt.add_sedge_between(ray_sv_tgt, lateral_sv_tgt[i]);
	
      } else {
	std::cerr << "Found nothing " << std::endl;
	opposite[i] = origin[i];
      }
    }

    for(int i=0; i<2; ++i) {
      if(origin[i] != opposite[i] || origin[1-i] == opposite[1-i]) 
	continue;

      std::cerr << "special " << std::endl;
      SM_walls SMW(&*opposite[i]);
      
      SVertex_handle src(i==0?ein:ein->twin());
      Sphere_circle c(Sphere_point(CGAL::ORIGIN - dir), src->point());
      Sphere_segment sphere_ray(src->point(), Sphere_point(CGAL::ORIGIN - dir), c);
      lateral_sv_tgt[i] = SMW.add_lateral_svertex(sphere_ray);
      SMW.add_sedge_between(src, lateral_sv_tgt[i]);    
    }

    Sphere_circle c(Sphere_point(CGAL::ORIGIN - dir), lateral_sv_tgt[0]->point());
    do {
      Object_handle o2 = pl->shoot(Ray_3(lateral_sv_tgt[0]->source()->point(),
					 lateral_sv_tgt[0]->point()));
      
      Halffacet_handle f;
      Halfedge_handle e;
      Vertex_handle v;
      if(assign(f,o2))
	CGAL_assertion_msg(false, "wrong handle");
      else if(assign(e,o2)) {
	CGAL_assertion_msg(false, "not implemented yet");
      } else if(assign(v,o2)) {
	if(v == lateral_sv_tgt[1]->source()) {
	  lateral_sv_tgt[0]->twin() = lateral_sv_tgt[1];
	  lateral_sv_tgt[1]->twin() = lateral_sv_tgt[0];
	  return;
	}
	SM_walls SMW(&*lateral_sv_tgt[0]->source());
	SVertex_handle opp = SMW.add_ray_svertex(lateral_sv_tgt[0]->point().antipode());	
	opp->twin() = lateral_sv_tgt[0];
	lateral_sv_tgt[0]->twin() = opp;

	lateral_sv_tgt[0] = 
	  SMW.add_lateral_svertex(Sphere_segment(lateral_sv_tgt[0]->point().antipode(),
						 lateral_sv_tgt[0]->point(), c));
      }
    } while(true);
    
    //    CGAL::SNC_io_parser<SNC_structure> O(std::cerr, *sncp, false);
    //    O.print();

  }
};

CGAL_END_NAMESPACE
#endif //CGAL_NEF3_SINGLE_WALL_CREATOR_H
