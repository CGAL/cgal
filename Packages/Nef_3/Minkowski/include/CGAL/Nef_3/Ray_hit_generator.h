#ifndef CGAL_NEF3_RAY_HIT_GENERATOR_H
#define CGAL_NEF3_RAY_HIT_GENERATOR_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_S2/SM_walls.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Ray_hit_generator : public Modifier_base<typename Nef_::SNC_and_PL> {
  
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

  typedef typename Base::Segment_3               Segment_3;
  typedef typename Base::Point_3                 Point_3;
  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::SVertex_handle           SVertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::Vertex_iterator         Vertex_iterator;
  typedef typename Base::SVertex_iterator         SVertex_iterator;

  Vector_3 dir;
  SNC_structure* sncp;
  SNC_point_locator* pl;
  int mask;

 public:
  Ray_hit_generator(Vector_3 d = Vector_3()) : dir(d), mask(255) {}
  Ray_hit_generator(SNC_structure* sncin, SNC_point_locator* plin, int m) 
    : sncp(sncin), pl(plin), mask(m) {}
      
  Vertex_handle create_vertex_on_first_hit(const Ray_3& r) {

    std::cerr << "shoot ray in SNC " << r << std::endl;

    //    CGAL_NEF_SETDTHREAD(503*509);
    Object_handle o = pl->shoot(r, mask);
    //    CGAL_NEF_SETDTHREAD(1);

    Vertex_handle v;
    if(assign(v, o)) {
      std::cerr << "Found vertex " << v->point() << std::endl;
      return v;
    }

    Point_3 ip;
    SNC_intersection I;
    SNC_constructor C(*sncp,pl);

    Halfedge_handle e;
    if(assign(e, o)) {
      std::cerr << "Found edge " << e->source()->point() 
		<< "->" << e->twin()->source()->point() << std::endl;
      Segment_3 seg(e->source()->point(), e->twin()->source()->point());
      I.does_intersect_internally(r, seg, ip);
      ip = normalized(ip);
      v = C.create_from_edge(e,ip);
      pl->add_vertex(v);

      std::cerr << "new vertex " << ip << std::endl;

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

      return v;
    }

    Halffacet_handle f;
    if(assign(f, o)) {
      std::cerr << "Found facet " << std::endl;
      I.does_intersect_internally(r, f, ip);
      ip = normalized(ip);
      v = C.create_from_facet(f,ip);
      pl->add_vertex(v);

      std::cerr << "new vertex " << ip << std::endl;

      return v;
    }

    CGAL_assertion_msg(false, "ray should hit vertex, edge, or facet");
    return Vertex_handle();
  }

  void operator()(SNC_and_PL& sncpl) {
    
    sncp = sncpl.sncp;
    pl = sncpl.pl;

    Vertex_iterator vi;
    for(vi = sncp->vertices_begin(); vi != sncp->vertices_end(); ++vi) {
      SM_walls smw(&*vi);
      if(smw.need_to_shoot(Sphere_point(dir))) {
	Ray_3 r(vi->point(), dir);
	Vertex_handle v_new = create_vertex_on_first_hit(r);
	SM_walls smw(&*v_new);
	smw.add_ray_svertex(Sphere_point(-dir));
      }
    }
  }
};
  
CGAL_END_NAMESPACE
#endif //CGAL_NEF3_RAY_HIT_GENERATOR_H
