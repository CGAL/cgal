#ifndef CGAL_NEF3_YVERTICAL_WALL_BUILDER_H
#define CGAL_NEF3_YVERTICAL_WALL_BUILDER_H

#include<CGAL/Nef_3/SNC_decorator.h>
#include<CGAL/Nef_3/Single_wall_creator2.h>

namespace CGAL {

template<typename Nef_>
class YVertical_wall_builder : public Modifier_base<typename Nef_::SNC_structure> {

  typedef Nef_                                            Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure          SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;
  typedef CGAL::Single_wall_creator2<Nef_polyhedron>   Single_wall2;

  typedef typename SNC_structure::Vertex_handle           Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle        Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle        SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle        SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle            SFace_handle;

  typedef typename SNC_structure::Volume_iterator         Volume_iterator;
  typedef typename SNC_structure::Halfedge_iterator       Halfedge_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;

  typedef typename SNC_structure::Vector_3                Vector_3;
  typedef typename SNC_structure::Point_3                 Point_3;

  typedef typename SNC_structure::Sphere_point            Sphere_point;
  typedef typename SNC_structure::Sphere_circle           Sphere_circle;
  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

  typedef typename std::list<Halfedge_handle>             Edge_list;
 public:
  typedef typename std::list<Halfedge_handle>::iterator   Vertical_redge_iterator;


 private:
  struct Reflex_edge_visitor {

    Edge_list& pos;
    Edge_list& neg;
    Sphere_point dir;

    Reflex_edge_visitor(Edge_list& p, Edge_list& n, Sphere_point dir_in)
        : pos(p), neg(n), dir(dir_in) {}

    void visit(Vertex_handle v) const {}
    void visit(Halfedge_handle e) const {}
    void visit(Halffacet_handle f) const {}
    void visit(SHalfloop_handle sl) const {}
    void visit(SFace_handle sf) const {}
    void visit(SHalfedge_handle se) const {
      int isrse = CGAL::is_reflex_sedge<SNC_structure>(se, dir);
      if((isrse&1)==1) pos.push_back(se->source());
      if((isrse&2)==2) neg.push_back(se->source());
      return;

      Halfedge_handle e = se->source();
      if(e->is_twin())
        return;
      SHalfedge_handle se2 = se->sprev()->twin();

      if(se2==se) {
        pos.push_back(e);
        neg.push_back(e);
      }

      Vector_3 vec1 = e->point() - CGAL::ORIGIN;
      Vector_3 vec2 = se->circle().orthogonal_vector();
      Sphere_point sp1 = CGAL::ORIGIN + cross_product(vec2,vec1);

      if(se2->circle().oriented_side(sp1) != ON_POSITIVE_SIDE)
        return;

      vec2 = se2->circle().orthogonal_vector();
      Sphere_point sp2 = CGAL::ORIGIN + cross_product(vec2,vec1);
      Sphere_segment s(sp1, sp2, Sphere_circle(sp2,sp1));
      //      std::cerr << "e:" << e->source()->point()
      //                << "->" << e->twin()->source()->point() << std::endl;
      //      std::cerr << "s:" << sp1 << "->" << sp2 << std::endl;
      CGAL_assertion(s.is_long());

      Sphere_point sp(dir);
      if(s.has_on(sp) && s.source() != sp && s.target() != sp)
        pos.push_back(e);

      sp = sp.antipode();
      if(s.has_on(sp) && s.source() != sp && s.target() != sp)
        neg.push_back(e);

    }
  };

  Sphere_point dir;
  Edge_list pos;
  Edge_list neg;

 public:
  YVertical_wall_builder(Sphere_point dir_in = Sphere_point(0,1,0))
    : dir (dir_in) {}

    /*
    void operator()(SNC_structure& snc) {
      Halfedge_iterator ei;
      CGAL_forall_halfedges(ei, snc) {
        if(ei->point() != Sphere_point(1,0,0)) continue;
        SHalfedge_around_svertex_circulator
          svc(ei->out_sedge()), send(svc);
        CGAL_For_all(svc, send) {
          if(!svc->incident_sface()->mark()) continue;
          if(!CGAL::is_reflex_sedge_in_any_direction<SNC_structure>(svc))
            continue;
          pos.push_back(svc);
        }
      }
    }
    */

  void operator()(SNC_structure& snc) {
    pos.clear();
    neg.clear();
    SNC_decorator D(snc);
    Reflex_edge_visitor rev(pos,neg,dir);
    Volume_iterator ci;
    for(ci=snc.volumes_begin(); ci!=snc.volumes_end(); ++ci)
      if(ci->mark())
        D.visit_shell_objects(SFace_handle(ci->shells_begin()),rev);
  }

  Vertical_redge_iterator pos_begin() { return pos.begin(); }
  Vertical_redge_iterator pos_end() { return pos.end(); }
  Vertical_redge_iterator neg_begin() { return neg.begin(); }
  Vertical_redge_iterator neg_end() { return neg.end(); }
};

} //namespace CGAL
#endif // CGAL_NEF3_YVERTICAL_WALL_BUILDER_H
