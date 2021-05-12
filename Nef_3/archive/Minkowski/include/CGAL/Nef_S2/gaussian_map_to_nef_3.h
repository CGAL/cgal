#ifndef CGAL_MINKOWSKI_GAUSSIAN_MAP_TO_NEF3_H
#define CGAL_MINKOWSKI_GAUSSIAN_MAP_TO_NEF3_H

#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/Gaussian_map.h>
#include <CGAL/Modifier_base.h>

namespace CGAL {

template<typename Nef3>
  class gaussian_map_to_nef_3 : public Modifier_base<typename Nef3::SNC_structure > {

  typedef typename Nef3::Kernel                   Kernel;
  typedef typename Nef3::SNC_structure            SNC_structure;
  typedef typename SNC_structure::Sphere_map      Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>          SM_decorator;
  typedef CGAL::Gaussian_map<Kernel>               Gaussian_map;

  typedef typename Gaussian_map::SFace_const_iterator       SFace_const_iterator;
  typedef typename Gaussian_map::SFace_const_handle         SFace_const_handle;
  typedef typename Gaussian_map::SHalfedge_const_iterator   SHalfedge_const_iterator;
  typedef typename Gaussian_map::SHalfedge_const_handle     SHalfedge_const_handle;
  typedef typename Gaussian_map::SVertex_const_iterator     SVertex_const_iterator;
  typedef typename Gaussian_map::SVertex_const_handle       SVertex_const_handle;
  typedef typename Gaussian_map::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;

  typedef typename SNC_structure::Vertex_handle            Vertex_handle;
  typedef typename SNC_structure::SVertex_handle           SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle         SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle             SFace_handle;

  typedef typename SNC_structure::Sphere_circle            Sphere_circle;
  typedef typename SNC_structure::Sphere_point             Sphere_point;

  const Gaussian_map& G;

 public:
  gaussian_map_to_nef_3(const Gaussian_map& Gin) : G(Gin) {}

  void operator()(SNC_structure& snc) {

    snc.clear();

#ifdef CGAL_NEF_INDEXED_ITEMS
    CGAL::Unique_hash_map<SHalfedge_const_handle, int> SE2i;
    SHalfedge_const_iterator sei;
    CGAL_forall_sedges(sei, G) {
      SE2i[sei] = Index_generator::get_unique_index();
      SE2i[sei->twin()] = SE2i[sei];
    }

    CGAL::Unique_hash_map
      <SVertex_const_handle, std::pair<int, int> > SV2i;
    SVertex_const_iterator svi;
    CGAL_forall_svertices(svi, G)
      SV2i[svi] = std::pair<int, int>
      (Index_generator::get_unique_index(),
       Index_generator::get_unique_index());
#endif

    CGAL::Unique_hash_map<SFace_const_handle, Vertex_handle> sface2vertex;
    SFace_const_iterator sfi;
    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      sface2vertex[sfi] = snc.new_vertex(sfi->mark().point(),
                                         sfi->mark().boolean());
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
      sv_first->mark() = sfc->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
      sv_first->set_index(SE2i[sfc]);
#endif
      ++sfc;
      sv_prev = sv =
        SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv->mark() = sfc->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
      sv->set_index(SE2i[sfc]);
#endif
      se_first = se_prev = SM.new_shalfedge_pair(sv_first, sv_prev);
      se_first->mark() = se_first->twin()->mark() = sfc->source()->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
      se_first->set_index(SV2i[sfc->source()].first);
      se_first->twin()->set_index(SV2i[sfc->source()].second);
#endif
      se_first->circle() = Sphere_circle(sv_first->point(), sv->point());
      se_first->circle() = normalized(se_first->circle());
      se_first->twin()->circle() = se_first->circle().opposite();

      ++sfc;
      CGAL_For_all(sfc,sfcend) {
        sv = SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
        sv->mark() = sfc->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
        sv->set_index(SE2i[sfc]);
#endif
        se = SM.new_shalfedge_pair(sv_prev, sv);
        se->mark() = se->twin()->mark() = sfc->source()->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
        se->set_index(SV2i[sfc->source()].first);
        se->twin()->set_index(SV2i[sfc->source()].second);
#endif
        se->circle() = Sphere_circle(sv_prev->point(), sv->point());
        se->circle() = normalized(se->circle());
        se->twin()->circle() = se->circle().opposite();
        se->sprev() = se_prev;
        se_prev->snext() = se;
        sv_prev = sv;
        se_prev = se;
      }

      se = SM.new_shalfedge_pair(sv_prev, sv_first);
      se->mark() = se->twin()->mark() = sfc->source()->mark().boolean();
#ifdef CGAL_NEF_INDEXED_ITEMS
      se->set_index(SV2i[sfc->source()].first);
      se->twin()->set_index(SV2i[sfc->source()].second);
#endif
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
      SM.link_as_face_cycle(se,sf0);
      SM.link_as_face_cycle(se->twin(),sf1);
    }

    //    SNC_io_parser<SNC_structure> O0(std::cerr,snc);
    //    O0.print();
    //    CGAL_NEF_SETDTHREAD(43*31);
  }


};

} //namespace CGAL
#endif // CGAL_MINKOWSKI_GAUSSIAN_MAP_TO_NEF3_H
