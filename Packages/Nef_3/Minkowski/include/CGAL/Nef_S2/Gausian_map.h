#ifndef CGAL_NEF_GAUSIAN_MAP
#define CGAL_NEF_GAUSIAN_MAP

#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_polyhedron_S2.h>

#include <CGAL/IO/Qt_widget_Nef_S2.h>
#include <qapplication.h>

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


  template<typename Nef_polyhedron_3>
  class SVertex_creator {     

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator   
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

    typedef typename Nef_polyhedron_3::Point_3                  Point_3;

    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;
    typedef CGAL::Unique_hash_map<Vertex_const_handle, bool> Vertex2bool_hash;
    typedef CGAL::Unique_hash_map<Halfedge_const_handle, bool> Edge2bool_hash;

    SM_decorator SM;
    Facet2SVertex_hash& Facet2SVertex;
    Vertex2bool_hash& omit_vertex;
    Edge2bool_hash& omit_edge;

  public:
    SVertex_creator(Sphere_map* smap, Facet2SVertex_hash& F2SV, Vertex2bool_hash& V2b, Edge2bool_hash& E2b)
      : SM(smap), Facet2SVertex(F2SV), omit_vertex(V2b), omit_edge(E2b) {}

  private:
    bool svertex_exists(const Point_3& p, SVertex_handle& sv) {
      SVertex_iterator svi;
      CGAL_forall_svertices(svi, SM) {
	if((Point_3) svi->point() == p) {
	  sv = svi;
	  return true;
	}
      }
      return false;
    }

  public:
      void visit(Vertex_const_handle v) {}
      void visit(Halfedge_const_handle e) {}
      void visit(SHalfedge_const_handle se) {}
      void visit(SHalfloop_const_handle sl) {}
      void visit(SFace_const_handle sf) {}
  
      void visit(Halffacet_const_handle f) {
	Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
	SHalfedge_const_handle se(fc);
	SHalfedge_around_facet_const_circulator hc(se), hend(hc);

	CGAL_For_all(hc,hend) {
	  std::cerr << "circles " << hc->circle() << " + " << hc->sprev()->circle() << std::endl;
	  if(hc->sprev()->circle() == hc->circle()) {
	    omit_vertex[hc->source()->source()] = omit_vertex[hc->source()->twin()->source()] = true;
	    omit_edge[hc->source()] = omit_edge[hc->source()->twin()] = true;
	  }
	  if(hc->source()->point() == hc->twin()->source()->point().antipode())
	    omit_edge[hc->source()] = omit_edge[hc->source()->twin()] = true;
	}

	std::cerr << "SVertex_creator " << f->plane() << std::endl;
	SVertex_handle sv;
	if(!svertex_exists(CGAL::ORIGIN + f->twin()->plane().orthogonal_vector(), sv)) {
	  sv = SM.new_svertex(f->twin()->plane().orthogonal_vector());
	  sv->mark() = Mark(0,0,0);
	}
	Facet2SVertex[f] = sv;
      }
  };

  template<typename Nef_polyhedron_3>
    class SEdge_creator {

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef typename SM_decorator::SHalfedge_around_svertex_circulator
      SHalfedge_around_svertex_circulator;

    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator   
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;
    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;
    typedef CGAL::Unique_hash_map<Halfedge_const_handle, bool> Edge2bool_hash;

    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;
    Facet2SVertex_hash& Facet2SVertex;
    Edge2bool_hash& omit_edge;

  public:
    SEdge_creator(Sphere_map* smap, Edge2SEdge_hash& E2SE, Facet2SVertex_hash& F2SV,Edge2bool_hash& E2b) : 
      SM(smap), Edge2SEdge(E2SE), Facet2SVertex(F2SV), omit_edge(E2b) {}
      
  private:
    bool sedge_exists(SVertex_handle s, SVertex_handle t, SHalfedge_handle& se) {
      SHalfedge_around_svertex_circulator sh(s->out_sedge()), send(sh);
      CGAL_For_all(sh, send)
	if(sh->twin()->source() == t) {
	  se = sh;
	  return true;
	}
      return false;
    }
    
  public:
    void visit(Vertex_const_handle v) {}
    void visit(Halfedge_const_handle e) {}
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
    void visit(SFace_const_handle sf) {}
    
    void visit(Halffacet_const_handle f) {
      std::cerr << "SEdge_creator " << f->plane() << std::endl;
      Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hc(se), hend(hc);
      do{--hc;}
      while(!omit_edge[hc->source()] && hc != hend);
      hend = hc;
      do {
	if(!omit_edge[hc->source()]) {
	  std::cerr << "edge " << hc->source()->source()->point() 
		    << ":" << hc->source()->point() << std::endl;

	  Halfedge_const_handle e = hc->source();
	  SHalfedge_handle set = Edge2SEdge[e->twin()];
	  if(set == SHalfedge_handle())
	    Edge2SEdge[e] = SM.new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	  else {
	    SM.link_as_target_and_append(Facet2SVertex[f], set);
	    set->mark() = set->twin()->mark() = Mark(0,0,0);
	    set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	    set->twin()->circle() = set->circle().opposite();
	    Edge2SEdge[e] = set->twin();
	  }
	  
	  CGAL::SM_io_parser<SM_decorator> O(std::cerr,SM); 
	  O.print();
	}
	++hc;
      } while(hc != hend);
    }
  };
  

  template<typename Nef_polyhedron_3>
    class SFace_creator {     
  
    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;
    typedef CGAL::Unique_hash_map<Vertex_const_handle, bool> Vertex2bool_hash;

    const Nef_polyhedron_3& N3;
    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;
    Vertex2bool_hash& omit_vertex;

  public:
    SFace_creator(const Nef_polyhedron_3& N, Sphere_map* smap, Edge2SEdge_hash& E2SE, Vertex2bool_hash& V2b) : 
      N3(N), SM(smap), Edge2SEdge(E2SE), omit_vertex(V2b) {}

      void visit(Halfedge_const_handle e) {}
      void visit(SHalfedge_const_handle se) {}
      void visit(SHalfloop_const_handle sl) {}
      void visit(SFace_const_handle sf) {}
      void visit(Halffacet_const_handle f) {}

      void visit(Vertex_const_handle v) {
	std::cerr << "SFace_creator " << v->point() << std::endl;

	if(omit_vertex[v]) {
	  std::cerr << "omit " << v->point() << std::endl;
	  return;
	}

        typename Nef_polyhedron_3::Nef_polyhedron_S2 SD(N3.get_sphere_map(v));
	typename Nef_polyhedron_3::Halfedge_const_iterator ei(SD.svertices_begin());
	SHalfedge_handle se = Edge2SEdge[ei];
	while(se == SHalfedge_handle()) {
	  ++ei;
	  se = Edge2SEdge[ei];
	}

	CGAL_assertion(ei != SD.svertices_end());

	SFace_handle sf = SM.new_sface();
	sf->mark() = v->point();
	SM.link_as_face_cycle(se,sf);
      }
  };


  struct VECTOR_ADDITION { 
    Mark operator()(const Mark& b1, const Mark& b2) const {
      return b1+(b2-CGAL::ORIGIN); 
    } 
  };

 public:
  Gausian_map() : Base(new Sphere_map) {}

  template<typename NK> 
    Gausian_map(const CGAL::Nef_polyhedron_3<NK>& N3,
		typename CGAL::Nef_polyhedron_3<NK>::Volume_const_iterator c) : Base(new Sphere_map) {

    typedef CGAL::Nef_polyhedron_3<NK> Nef_polyhedron_3;
    typedef typename Nef_polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_const_iterator
      Halffacet_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;
    typedef typename Nef_polyhedron_3::Vertex_const_handle
      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Volume_const_handle
      Volume_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle
      Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle
      SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle
      SFace_const_handle;

    Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;
    Unique_hash_map<Vertex_const_handle, bool> Vertex2bool(false);
    Unique_hash_map<Halfedge_const_handle, bool> Edge2bool(false);

    SFace_const_handle sf = c->shells_begin();

    SVertex_creator<Nef_polyhedron_3> create_svertices(sphere_map(), Facet2SVertex, Vertex2bool, Edge2bool);
    SEdge_creator<Nef_polyhedron_3>   create_sedges(sphere_map(), Edge2SEdge, Facet2SVertex, Edge2bool);
    SFace_creator<Nef_polyhedron_3>   create_sfaces(N3, sphere_map(), Edge2SEdge, Vertex2bool);

    N3.visit_shell_objects(sf, create_svertices);
    N3.visit_shell_objects(sf, create_sedges);
    N3.visit_shell_objects(sf, create_sfaces);

    CGAL::SM_io_parser<SM_decorator> O(std::cerr,*this); 
    O.print();
  }

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

    Volume_const_handle c(--N3.volumes_end());

    Halffacet_const_iterator f;
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      SVertex_handle sv = new_svertex(f->twin()->plane().orthogonal_vector());
      sv->mark() = Mark(0,0,0);
      Facet2SVertex[f] = sv;
    }
    
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
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

      CGAL::SM_io_parser<SM_decorator> O(std::cerr,*this); 
      O.print();

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

    void visualize() {
      int argc=1;
      char* argv[argc];
      argv[0] = "Gaussian Map Viewer";

      typedef typename CGAL::Nef_polyhedron_S2<K,CGAL::SM_items,Mark> Nef_polyhedron_S2;
      typedef typename Nef_polyhedron_S2::Const_decorator SM_const_decorator;
      SM_const_decorator SMCD(sphere_map());

      QApplication a(argc, argv);
      CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>* w = 
	new CGAL::Qt_widget_Nef_S2<Nef_polyhedron_S2>(SMCD);
      a.setMainWidget(w);
      w->show();
      a.exec();     
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
