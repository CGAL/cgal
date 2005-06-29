#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/Gausian_map.h>
#include <CGAL/Modifier_base.h>

CGAL_BEGIN_NAMESPACE

template<typename Kernel_, typename Items_, typename Mark_>
  class gausian_map_to_nef_3 : public Modifier_base<SNC_structure<Kernel_,Items_,Mark_> > {
    
  typedef Kernel_                                 Kernel;
  typedef Items_                                  Items;
  typedef Mark_                                   Mark;
  typedef CGAL::SNC_structure<Kernel,Items,Mark>  SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>      SNC_decorator;
  typedef typename SNC_structure::Sphere_map      Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>          SM_decorator;
  typedef CGAL::Gausian_map<Kernel>               Gausian_map;

  typedef typename Gausian_map::SFace_const_iterator       SFace_const_iterator;
  typedef typename Gausian_map::SFace_const_handle         SFace_const_handle;
  typedef typename Gausian_map::SHalfedge_const_handle     SHalfedge_const_handle;
  typedef typename Gausian_map::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;

  typedef typename SNC_structure::Vertex_handle            Vertex_handle;
  typedef typename SNC_structure::SVertex_handle           SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle         SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle             SFace_handle;
  
  typedef typename SNC_structure::Sphere_circle            Sphere_circle;
  typedef typename SNC_structure::Sphere_point             Sphere_point;

  Gausian_map& G;

 public:
  gausian_map_to_nef_3(Gausian_map& Gin) : G(Gin) {}
    
  void operator()(SNC_structure& snc) {

    Unique_hash_map<SFace_const_handle, Vertex_handle> sface2vertex;
    
    SNC_decorator SD(snc);

    SFace_const_iterator sfi;
    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      sface2vertex[sfi] = snc.new_vertex(sfi->mark(),true);
    }

    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      Vertex_handle v = sface2vertex[sfi];
      SM_decorator SM(&*v);

      SHalfedge_const_handle sec = sfi->sface_cycles_begin();
      SHalfedge_around_sface_const_circulator sfc(sec), sfcend(sfc);
      
      SVertex_handle sv, sv_prev, sv_first;
      SHalfedge_handle se, se_prev, se_first;
      sv_first = 
	SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv_first->mark() = true;
      ++sfc;
      sv_prev = sv = 
	SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv->mark() = true;
      ++sfc;
      se_first = se_prev = SM.new_shalfedge_pair(sv_first, sv_prev);
      se_first->mark() = se_first->twin()->mark() = true;
      se_first->circle() = Sphere_circle(sv_first->point(), sv->point());
      se_first->circle() = normalized(se_first->circle());
      se_first->twin()->circle() = se_first->circle().opposite();

      CGAL_For_all(sfc,sfcend) {
	sv = SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
	sv->mark() = true;
	se = SM.new_shalfedge_pair(sv_prev, sv);
	se->mark() = se->twin()->mark() = true;
	se->circle() = Sphere_circle(sv_prev->point(), sv->point());
	se->circle() = normalized(se->circle());
	se->twin()->circle() = se->circle().opposite();
	se->sprev() = se_prev;
	se_prev->snext() = se;	
	sv_prev = sv;
	se_prev = se;
      }
     
      se = SM.new_shalfedge_pair(sv_prev, sv_first);
      se->mark() = se->twin()->mark() = true;
      se->circle() = Sphere_circle(sv_prev->point(), sv_first->point());
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();
      se->sprev() = se_prev;
      se_prev->snext() = se;
      se_first->sprev() = se;
      se->snext() = se_first;

      // orientation ?
      
      SFace_handle sf0 = SM.new_sface();
      SFace_handle sf1 = SM.new_sface();
      sf0->mark() = false;
      sf1->mark() = true;

      /*
      Sphere_point p1 = se->source()->point();
      Sphere_point p2 = se->snext()->source()->point();
      Sphere_point p3 = se->snext()->snext()->source()->point();
      if(spherical_orientation(p1,p2,p3) > 0)
	se = se->twin();
      */

      SM.link_as_face_cycle(se,sf0);
      SM.link_as_face_cycle(se->twin(),sf1);
    }
    
    //    SNC_io_parser<SNC_structure> O0(std::cerr,snc);
    //    O0.print();   
    //    CGAL_NEF_SETDTHREAD(43*31);
  }


};

CGAL_END_NAMESPACE
