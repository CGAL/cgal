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
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;
  typedef typename Base::Object_handle           Object_handle;

  Halfedge_handle ein;
  Vector_3 dir;
  
 public:
  Single_wall_creator(Halfedge_handle e, Vector_3 d)
    : ein(e), dir(d) {}

  void operator()(SNC_and_PL& sncpl) {
    
    SNC_structure* sncp(sncpl.sncp);
    SNC_point_locator* pl(sncpl.pl);

    SNC_constructor C(*sncp,pl);

    Vertex_handle origin[2];
    Vertex_handle opposite[2];
    origin[0] = ein->source();
    origin[1] = ein->target();

    for(int i=0; i<2; ++i) {
      SM_point_locator P(&*origin[i]);
      Object_handle o = P.locate(Sphere_point(dir));
      SVertex_handle sv;
      SHalfedge_handle se;
      SHalfloop_handle sl;
      SFace_handle sf;
      if(assign(sv,o)) {
	opposite[i]=sv->twin()->source();
	std::cerr << " Found vertex directly !!!! " << std::endl;
      }
      else {
	Vertex_handle v;
	Halfedge_handle e;
	Halffacet_handle f;
	Ray_3 r(origin[i]->point(),dir);
	Object_handle o2 = pl->shoot(r);
	if(assign(f,o2))
	  std::cerr << "Found facet " << std::endl;
	else if(assign(e,o2)) {
	  std::cerr << "Found edge " << std::endl;
	  Point_3 ip;
	  SNC_intersection I;
	  I.does_intersect_internally(r, Segment_3(e->source()->point(),
						   e->twin()->source()->point()),  
				      ip);
	  opposite[i] = C.create_from_edge(e,ip);
	  SM_walls SMW(&*opposite[i]);
	  SMW.add_two(i==0?ein->point():ein->twin()->point(),
		      Sphere_point(CGAL::ORIGIN - dir));
	} else if(assign(v,o2)) {
	  std::cerr << "Found vertex " << std::endl;
	  opposite[i] = v;
	  SM_walls SMW(&*opposite[i]);
	  SMW.add_two(i==0?ein->point():ein->twin()->point(),
		      Sphere_point(CGAL::ORIGIN - dir));
	} else {
	  std::cerr << "Found nothing " << std::endl;
	  opposite[i] = origin[i];
	}
      }
    }

    SNC_point_locator* old_pl = pl;
    pl = pl->clone();
    sncpl.pl = pl;
    delete old_pl;
    C=SNC_constructor(*sncp,pl);
    C.clear_external_structure();
    C.build_external_structure();
  }
};

CGAL_END_NAMESPACE
#endif //CGAL_NEF3_SINGLE_WALL_CREATOR_H
