#include<CGAL/Nef_3/SNC_decorator.h>

template<typename SNC_>
class SNC_walls : SNC_decorator<SNC_>{
  
  typedef SNC_                                   SNC_strucutre;
  typedef SNC_decorator<SNC_structure>           Base;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typenmae Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;

  SNC_point_locator* pl;
  
  SNC_walls(SNC_structure& S, SNC_point_locator* spl) 
    : Base(S), pl(spl) {}

  /*
  Vector_3 dir;

  Unique_hash_map<Vertex_handle, Vertex_handle> opposite_up;
  Unique_hash_map<Vertex_handle, Vertex_handle> opposite_down;
  Unique_hash_map<Vertex_handle, Object_handle> object_up;
  Unique_hash_map<Vertex_handle, Object_handle> object_down;

  Vertex_handle create_opposite_vertex(Vertex_const_handle v, bool opp) {
    
    Vector_3 vec = opp ? dir.opposite : dir;
    Object_handle o;
    if(!opp)
      o = object_up[vi];
    else
      o = object_down[vi];

    SVertex_handle svh;
    if(assign(svh,o))
      return svh->twin()->source();

    int hit = 0;
    Object_handle o2 = pl()->shoot(vec);
    Vertex_handle vh;
    Halfedge_handle eh;
    Halffacet_handle fh;
    if(assign(fh,o2)) hit = 1;
    else if(assign(eh,o2)) hit = 2;
    else if(assign(vh,o2)) hit = 3;
    else return Vertex_handle();

    SHalfedge_handle seh;
    SFace_handle sfh;
    if(assign(sfh,o)) {
    SM_walls SW(vi);
      SW.insert_new_svertex_into_sface(sfh, Sphere_point(vec - CGAL::ORIGIN));
    } else if(assign(seh,o)) {
      SM_walls SW(vi);
      SW.insert_new_svertex_into_sedge(seh, Sphere_point(vec - CGAL::ORIGIN));
    } else CGAL_assertio_msg(false,"wrong object");

    switch(hit) {
    case 1  :
      Point_3 ip;
      SNC_intersection I;
      I.does_intersect_internally(Ray_3(v->point(),vec), fh,  ip);
      SM_walls SW(new_vertex(ip),fh->mark());
      SW.create_opposite_vertex_on_facet(fh, vi, Sphere_point(CGAL::ORIGIN-vec));
      break;
    case 2  :
      Point_3 ip;
      SNC_intersection I;
      I.does_intersect_internally(Ray_3(v->point(),vec), Segment(eh),  ip);     
      SM_walls SW(new_vertex(ip),fh->mark());
      SW.create_opposite_vertex_on_edge(eh, vi, Sphere_point(CGAL::ORIGIN-vec));
      break;
    case 3  :
      SM_walls SW(new_vertex(ip),fh->mark());
      SW.extend_vertex_by_inner_walls(vh, vi, Sphere_point(CGAL::ORIGIN-vec)); 
      break;
    default : CGAL_assertion_msg(false, "not implemented yet");
    }
    return Vertex_handle();
  }
  */

  void create_single_wall(Halfedge_handle ein, Vector_3 dir) {
    
    Vertex_handle origin[2];
    Vertex_handle opposite[2];
    origin[0] = ein->source();
    origin[1] = ein->target();

    for(int i=0; i<2; ++i) {
      SM_point_locator P(origin[i]);
      Object_handle o=P.shoot(dir);
      SVertex_handle sv;
      SHalfedge_handle se;
      SHalfloop_handle sl;
      SFace_handle sf;
      if(assign(sv,o)) {
	opposite[i]=sv->twin()->source();
	std::endl << " Found vertex directly !!!! " << std::endl;
      }
      else {
	Vertex_handle v;
	Halfedge_handle e;
	Halffacet_handle f;
	Object_handle o2 = pl()->shoot(dir);
	if(assign(f,o2))
	  std::cerr << "Found facet " << std::endl;
	else if(assign(e,o2)) 
	  std::cerr << "Found edge " << std::endl;
	else if(assign(v,o2))
	  std::cerr << "Found vertex " << std::endl;
	else
	  std::cerr << "Found nothing " << std::endl;
      }
    }
    
  }

  /*
  void create_walls() {

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, *this) {
      SM_point_locator P(vi);
      object_up[vi] = P.shoot(dir);
      object_up[vi] = P.shoot(dir.opposite());
    }

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, *this) {
      opposite_up[vi] = create_opposite_vertex(vi, dir);
      opposite_down[vi] = create_opposite_vertex(vi, dir.opposite());
    }

    Halfedge_iterator ei;
    CGAL_forall_vertices(ei, *this) {
      create_opposite_halfedge_up(opposite_up[ei->source()], 
				  opposite_up[ei->twin()->source()]);
      create_opposite_halfedge3_down(opposite_down[ei->source()], 
				     opposite_down[ei->twin()->source()]);
    }
  }
  */

};
