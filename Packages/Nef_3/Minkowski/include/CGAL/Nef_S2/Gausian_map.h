#ifndef CGAL_NEF_GAUSIAN_MAP
#define CGAL_NEF_GAUSIAN_MAP

#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>

CGAL_BEGIN_NAMESPACE

template <class K>
class Gausian_map : public CGAL::SM_decorator<CGAL::Sphere_map<CGAL::Sphere_geometry<K>,
                                                               CGAL::SM_items, 
                                                               typename K::Point_3> > {

  typedef CGAL::Sphere_geometry<K>                        Kernel;
  typedef typename Kernel::Point_3                        Mark;
  typedef CGAL::Sphere_map<Kernel,CGAL::SM_items,Mark>    Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                  SM_decorator;
  typedef SM_decorator                                    Base;
  typedef CGAL::SM_overlayer<SM_decorator>                SM_overlayer;
  typedef typename Kernel::Sphere_circle                  Sphere_circle;
 public:
  typedef typename Sphere_map::SVertex_handle             SVertex_handle;
  typedef typename Sphere_map::SHalfedge_handle           SHalfedge_handle;
  typedef typename Sphere_map::SFace_handle               SFace_handle;
  typedef typename Sphere_map::SFace_iterator             SFace_iterator;
  typedef typename Sphere_map::SFace_const_iterator       SFace_const_iterator;
  typedef typename Sphere_map::SHalfedge_iterator         SHalfedge_iterator;
  typedef typename Sphere_map::SVertex_iterator           SVertex_iterator;
  typedef typename Sphere_map::SVertex_const_iterator     SVertex_const_iterator;
  typedef typename Sphere_map::SHalfedge_around_svertex_const_circulator
    SHalfedge_around_svertex_const_circulator;
  typedef typename Sphere_map::SHalfedge_around_sface_circulator
    SHalfedge_around_sface_circulator;

  struct VECTOR_ADDITION { 
    Mark operator()(const Mark& b1, const Mark& b2) const {
      return b1+(b2-CGAL::ORIGIN); 
    } 
  };

 public:
  Gausian_map() : Base(new Sphere_map) {}

  template<typename NK> 
    Gausian_map(const CGAL::Nef_polyhedron_3<NK>& N3) : Base(new Sphere_map) {
    
    typedef CGAL::Nef_polyhedron_3<NK> Nef_polyhedron_3;
    typedef typename Nef_polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_const_iterator
      Halffacet_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;
    typedef typename Nef_polyhedron_3::Volume_const_handle
      Volume_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle
      Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle
      SHalfedge_const_handle;
    
    Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;

    Volume_const_handle c(N3.volumes_begin());

    Halffacet_const_iterator f;
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      SVertex_handle sv = new_svertex(f->plane().orthogonal_vector());
      sv->mark() = Mark(0,0,0);
      Facet2SVertex[f] = sv;
    }
    
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hc(se), hend(hc);
      do {
	Halfedge_const_handle e = hc->source();
	SHalfedge_handle set = Edge2SEdge[e->twin()];
	if(set == SHalfedge_handle())
	  Edge2SEdge[e] = new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	else {
	  link_as_target_and_append(Facet2SVertex[f], set);
	  set->mark() = set->twin()->mark() = Mark(0,0,0);
	  set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	  set->twin()->circle() = set->circle().opposite();
	  Edge2SEdge[e] = set->twin();
	}
	++hc;
      } while(hc != hend);
    }
    
    Vertex_const_iterator v;
    CGAL_forall_vertices(v,N3) {
      typename Nef_polyhedron_3::Nef_polyhedron_S2 SD(N3.get_sphere_map(v));
      Halfedge_const_handle e(SD.svertices_begin());
      SHalfedge_handle se = Edge2SEdge[e];
      SFace_handle sf = new_sface();
      sf->mark() = v->point();
      link_as_face_cycle(se,sf);
    }
  }
  
  template<typename PK> 
    Gausian_map(const CGAL::Polyhedron_3<PK>& P) : Base(new Sphere_map) {
    
    typedef CGAL::Polyhedron_3<PK> Polyhedron_3;
    typedef typename Polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Polyhedron_3::Facet_const_iterator
      Facet_const_iterator;
    typedef typename Polyhedron_3::Halfedge_around_facet_const_circulator
      Halfedge_around_facet_const_circulator;
    typedef typename Polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Polyhedron_3::Facet_const_handle
      Facet_const_handle;
    
    Unique_hash_map<Facet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;

    Facet_const_iterator f;
    for(f = P.facets_begin(); f != P.facets_end(); ++f) {
      SVertex_handle sv = new_svertex(f->plane().orthogonal_vector());
      sv->mark() = Mark(0,0,0);
      Facet2SVertex[f] = sv;
    }

    for(f = P.facets_begin(); f != P.facets_end(); ++f) {
      Halfedge_around_facet_const_circulator hc(f->facet_begin()),hend(hc);
      do {
	Halfedge_const_handle e = hc;
	SHalfedge_handle set = Edge2SEdge[e->opposite()];
	if(set == SHalfedge_handle())
	  Edge2SEdge[e] = new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	else {
	  link_as_target_and_append(Facet2SVertex[f], set,1);
	  set->mark() = set->twin()->mark() = Mark(0,0,0);
	  set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	  set->twin()->circle() = set->circle().opposite();
	  Edge2SEdge[e] = set->twin();
	}
	++hc;
      } while(hc != hend);
    }
   
    Vertex_const_iterator v;
    for(v = P.vertices_begin(); v != P.vertices_end(); ++v) {
      Halfedge_const_handle e(v->halfedge());
      SHalfedge_handle se = Edge2SEdge[e];
      SFace_handle sf = new_sface();
      sf->mark() = v->point();
      link_as_face_cycle(se,sf);
    }
  }

    void simplify() {
      CGAL_NEF_TRACEN("simplify");

      typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
      CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem(NULL);
      CGAL::Union_find< SFace_handle> UF;
  
      SFace_iterator f;
      CGAL_forall_sfaces(f,*this) {
	Pitem[f] = UF.make_set(f);
	clear_face_cycle_entries(f);
      }
      
      SHalfedge_iterator e;
      for(e = this->shalfedges_begin(); e != this->shalfedges_end(); ++e) { 
	if (e->is_twin() ) continue;
	CGAL_NEF_TRACEN("can simplify ? " << PH(e));
	CGAL_NEF_TRACEN(mark(e) << " " << mark(face(e)) << " " << mark(face(twin(e))));
	if (mark(face(e)) == mark(face(twin(e)))) {
	  CGAL_NEF_TRACEN("deleting "<<PH(e));
	  if ( !UF.same_set(Pitem[face(e)],
			    Pitem[face(twin(e))]) ) {
	    
	    UF.unify_sets( Pitem[face(e)],
			   Pitem[face(twin(e))] );
	    CGAL_NEF_TRACEN("unioning disjoint faces");
	  }
	  
	  CGAL_NEF_TRACEN("is_closed_at_source " << is_closed_at_source(e) << 
			  " " << is_closed_at_source(twin(e)));
       	  delete_edge_pair(e);
	}
      }
      
      CGAL::Unique_hash_map<SHalfedge_handle,bool> linked(false);
      for (e = this->shalfedges_begin(); e != this->shalfedges_end(); ++e) {
	if ( linked[e] ) continue;
	SHalfedge_around_sface_circulator hfc(e),hend(hfc);
	SFace_handle f = *(UF.find( Pitem[face(e)]));
	CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
	store_sm_boundary_object(e,f);
      }
      
      SVertex_iterator v;
      for(v = this->svertices_begin(); v != this->svertices_end(); ++v) {
	if ( is_isolated(v) ) {
	  delete_vertex_only(v);
	  continue;
	}
	if ( has_outdeg_two(v)) {
	  merge_edge_pairs_at_target(previous(first_out_edge(v))); 
	} 
      }
      
      for (f = this->sfaces_begin(); f != this->sfaces_end(); ++f) { 
	Union_find_handle pit = Pitem[f];
	if ( UF.find(pit) != pit ) {
	  CGAL_NEF_TRACEN("delete face " << &*f);
	  delete_face_only(f);
	}
      }
    }    
    
    void minkowski_sum(const Gausian_map& G1, const Gausian_map& G2) {
      SM_overlayer O(sphere_map());
      O.subdivide(G1.sphere_map(), G2.sphere_map());
      VECTOR_ADDITION va;
      O.select(va);
      simplify();
    }

    void dump() {
      SM_io_parser<Base>::dump(*this,std::cerr);
    }

    //    ~Gausian_map() { delete (Base*) this; }
};

template<typename Kernel>
std::ostream& operator<<(std::ostream& out, const CGAL::Gausian_map<Kernel>& G) {
  out << "OFF" << std::endl;
  out << G.number_of_sfaces() << " " << G.number_of_svertices() << " 0" << std::endl;
  
  typedef typename CGAL::Gausian_map<Kernel>::SFace_const_iterator SFace_const_iterator;
  CGAL::Unique_hash_map<SFace_const_iterator, int> SFace2int;
  
  int i=0;
  SFace_const_iterator sf;
  CGAL_forall_sfaces(sf, G) {
    SFace2int[sf] = i++;
    out << CGAL::to_double(sf->mark().x()) << " " 
	<< CGAL::to_double(sf->mark().y()) << " " 
	<< CGAL::to_double(sf->mark().z()) << std::endl;
  }

  typename CGAL::Gausian_map<Kernel>::SVertex_const_iterator sv;
  CGAL_forall_svertices(sv,G) {
    typename CGAL::Gausian_map<Kernel>::SHalfedge_around_svertex_const_circulator 
      svc(G.first_out_edge(sv)),
      svc1(svc),
      svend(svc);
    out << std::distance(++svc1,svend)+1;
    CGAL_For_all(svc,svend)
      out << " " << SFace2int[svc->incident_sface()];
    out << std::endl;
  }

  return out;
}

CGAL_END_NAMESPACE
#endif // CGAL_NEF_GAUSIAN_MAP
