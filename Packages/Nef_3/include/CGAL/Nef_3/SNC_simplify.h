#ifndef CGAL_SNC_SIMPLIFY_H
#define CGAL_SNC_SIMPLIFY_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#ifdef CGAL_NEF_DEBUG
#include <CGAL/Nef_3/SNC_io_parser.h>
#endif

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 41
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

template<typename SNC_structure>
class SNC_simplify : public SNC_decorator<SNC_structure> {
  
  typedef CGAL::SNC_simplify<SNC_structure>             Self;
  typedef CGAL::SNC_decorator<SNC_structure>            SNC_decorator;
  typedef SNC_decorator                                 Base;
  typedef typename SNC_structure::Sphere_map            Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                SM_decorator;
#ifdef CGAL_NEF_DEBUG
  typedef CGAL::SNC_io_parser<SNC_structure>            SNC_io_parser;
#endif
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator; 
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;  
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator; 

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
                                  SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator 
                                  SHalfedge_around_sface_circulator;

  typedef typename SNC_structure::Halffacet_cycle_iterator 
                                  Halffacet_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_iterator 
                                  SFace_cycle_iterator;

  typedef typename Base::Sphere_kernel Sphere_kernel;
  typedef typename Base::Sphere_point Sphere_point;
  typedef typename Base::Sphere_segment Sphere_segment;
  typedef typename Base::Sphere_circle Sphere_circle;
  typedef typename Base::Sphere_direction Sphere_direction;

#ifdef CGAL_NEF_DEBUG
  SNC_io_parser *IO;
#endif

 public:
  SNC_simplify(SNC_structure& sncs) : Base(sncs), IO(NULL) {}

  char PSE(SHalfedge_handle h)
    /* prints a sphere segment */ {
    SNC_decorator D;
    SM_decorator SD;

    TRACE(IO->index(h)<<" @ "<<IO->index(D.vertex(h))<<
	  " "<<D.point(D.vertex(h))<<
	  ", prev " <<IO->index(D.previous(h))<<
	  ", next " <<IO->index(D.next(h))<<
	  ", sprev "<<IO->index(SD.previous(h))<<
	  ", snext "<<IO->index(SD.next(h))<<
	  ", twin " <<IO->index(D.twin(h)));
    return ' ';
  }

  char PFC(SHalfedge_handle e)
    /* prints a facet cycle */ {
    TRACEN("--> Facet cycle begin");
    SHalfedge_around_facet_circulator c(e), cend(c);
    CGAL_For_all(c, cend) {
      TRACEN(PSE(c));
    }
    TRACE("--> Facet cycle end"); 
    return ' ';
  }

  char PFB(Halffacet_handle f)
    /* prints facet boundary entry points */ {
    Halffacet_cycle_iterator fc;
    CGAL_forall_facet_cycles_of(fc, f) {
      if(fc.is_shalfedge()) {
	SHalfedge_handle(e);
	TRACE(' '<<IO->index(e)); }
      else if(fc.is_shalfloop()) {
	SHalfloop_handle l(fc);
	TRACE(' '<<IO->index(l)<<"(sl)");
      }
    }
    return '.';
  }

  char PSFB(SFace_handle f)
    /* prints sphere face boundary entry points */ {
    SFace_cycle_iterator it;
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_shalfedge() ) TRACE(IO->index(SHalfedge_handle(it))<<' ');
    TRACE(", ");
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_svertex() ) TRACE(IO->index(SVertex_handle(it))<<' ');
    TRACE(", ");
    CGAL_forall_sface_cycles_of(it,f)
      if ( it.is_shalfloop() ) TRACE(IO->index(SHalfloop_handle(it)));
    return '.';
  }

  typedef typename Union_find< Volume_handle>::handle UFH_volume;
  typedef typename Union_find< Halffacet_handle>::handle UFH_facet;
  typedef typename Union_find< SFace_handle>::handle UFH_sface;

  void remove_f_including_all_edge_uses_in_its_boundary_cycles
    ( Halffacet_handle f,
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf )
    /* removes f and its boundary cycles, and merges up the sphere facets
       incident to them. */ {
    SNC_decorator D;
    Halffacet_cycle_iterator fc;
    CGAL_forall_facet_cycles_of(fc, f) {
      if(fc.is_shalfedge() ) {
	SHalfedge_handle e(fc);
	SHalfedge_around_facet_circulator u(e), eend(e);
	CGAL_For_all(u, eend) {
	  SFace_handle fu = D.sface(u), ftu = D.sface(D.twin(u));
	  TRACEN("sfUNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	  merge_sets( fu, ftu, hash, uf);
	  SM_decorator SD(&*D.vertex(u));
	  TRACEN("removing "<<IO->index(u)<<" & "<<IO->index(SD.twin(u)));
	  Halfedge_handle src(SD.source(u)), tgt(SD.target(u));
	  if ( SD.is_closed_at_source(u) ) 
	    SD.set_face( src, fu);
	  if ( SD.is_closed_at_source( SD.twin(u)) ) 
 	    SD.set_face( tgt, fu);
	  /* TO VERIFY: does is_closed_at_source(u) imply is_isolated(src)?
	     if it is true, the svertex face update is not necesary. */
	  SD.delete_edge_pair(u); 
	  if( SD.is_isolated(src))
	    //	    SD.delete_vertex_only(src);
	    SD.set_face(src,fu);
	  SM_decorator SD2(&*D.vertex(tgt));
	  if( SD2.is_isolated(tgt))
	    // SD.delete_vertex_only(tgt);
	    SD2.set_face(tgt,fu);
	  /* TO VERIFY: can both svertices be isolated at the same time? */
      }
      }
      else if(fc.is_shalfloop()) {
	SHalfloop_handle l(fc);
	// this code is currenlty not used, but it is potentially need 
	// in the future, e.g for complex marks or a relative interior 
	// function 
	SFace_handle fu = D.sface(l), ftu = D.sface(D.twin(l));
	TRACEN("UNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	merge_sets( fu, ftu, hash, uf);
	SM_decorator SD(&*D.vertex(l));
	TRACEN("removing "<<IO->index(l)<<" & "<<IO->index(SD.twin(l)));
	SD.delete_loop_only();
      }
    }
    TRACEN("removing "<<IO->index(f)<<" & "<<IO->index(D.twin(f)));
    this->sncp()->delete_halffacet_pair(f);
    return;
  }

  bool is_part_of_volume(Vertex_handle v) 
    /* determines if a vertex v is part of a volume, cheking if its local
       graph is trivial (only one sface with no boundary). */  {
    SM_decorator SD(&*v);
    CGAL_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
    if( is_empty_range( SD.svertices_begin(), SD.svertices_end()) &&
	is_empty_range( SD.shalfedges_begin(), SD.shalfedges_end()) &&
	!SD.has_shalfloop())
      return true;
    return false;
  }

  bool is_part_of_facet(Vertex_handle v) 
    /* determines if a vertex v is part of a the relative interior of a 
       facet, checking if its local graph consists just of a sloop and
       two incident sfaces. */ {
    SM_decorator SD(&*v);
    CGAL_assertion( !is_empty_range( SD.svertices_begin(),
					  SD.svertices_end()) ||
			 is_empty_range( SD.shalfedges_begin(),
					 SD.shalfedges_end()));
    return( SD.has_shalfloop() &&
	    is_empty_range( SD.svertices_begin(), SD.svertices_end()));
  }

  bool is_part_of_edge(Vertex_handle v) {
    /* determines if a vertex v is part of a edge, checking at its local 
       graph for exactly two antipodal vertices  */

    SM_decorator SD(&*v);
    if(SD.has_shalfloop())
      return false;
    if(SD.svertices_begin() == SD.svertices_end())
      return false;
    if(++(SD.svertices_begin()) == SD.svertices_end())
      return false;
    
    TRACE(SNC_decorator(*this->sncp()).point(v)<<" is in edge interior? ");
    SVertex_iterator sv(SD.svertices_begin());
    SVertex_handle p1(sv++), p2(sv++);
    if( sv != SD.svertices_end())
      return false;
    
    TRACE("has two svertices ");
    Sphere_point sp1(SD.point(p1)), sp2(SD.point(p2));
    return (sp1 == sp2.antipode());
  }

  bool vertex_simplification(bool snc_computed = true) {
    bool simplified = false;

#ifdef CGAL_NEF_DEBUG
    delete IO;
    IO = new SNC_io_parser(std::cerr, *this->sncp());
#endif

    SNC_decorator D(*this->sncp());

    Vertex_iterator v = (*this->sncp()).vertices_begin();
    while( v != (*this->sncp()).vertices_end()) {
      SM_decorator SD(&*v);
      Vertex_iterator v_next(v);
      v_next++;
      if( is_part_of_volume(v)) {
	TRACEN("mark("<<IO->index(v)<<")="<<D.mark(v)<<", "<<
	       "mark("<<IO->index(D.volume(SD.sfaces_begin()))<<")="<<
	       SD.mark(SD.sfaces_begin()));
	if(D.mark(v) == SD.mark(SD.sfaces_begin())) {
	  TRACEN("removing isolated vertex "<<IO->index(v));
	  this->sncp()->delete_vertex(v);
	  simplified = true;
	}
      }
      else if( is_part_of_facet(v)) {
	if( D.mark(v) == SD.mark(SD.shalfloop())) {
	  TRACEN("removing "<<IO->index(v)<<
		 " on facet ");
	  this->sncp()->delete_vertex(v);
	  simplified = true;
	}
      }
      else if( is_part_of_edge(v)) {
	SVertex_iterator sv(SD.svertices_begin());
	Halfedge_handle e1(sv++), e2(sv++);
	CGAL_assertion( sv == SD.svertices_end());
	if( D.mark(e1) == D.mark(v) && D.mark(v) == D.mark(e2)) {
	  TRACEN("merging "<<IO->index(e1)<<" & "<<IO->index(e2)<<
		 " in "<<IO->index(v));
	  if(snc_computed)
	    merge_halfedge_pairs( e1, e2);
	  else
	    this->sncp()->delete_vertex(v);
	  simplified = true;
	}
      }
      v = v_next;
    }
    return simplified;
  }
  
  bool simplify() {

    bool update_facets  =  false;
    bool update_sfaces  =  false;
    bool update_volumes =  false;

    TRACEN(">>> simplifying");
    SNC_decorator D(*this->sncp());

#ifdef CGAL_NEF_DEBUG
    delete IO;
    IO = new SNC_io_parser(std::cerr, *this->sncp());
#endif
    
    Unique_hash_map< Volume_handle, UFH_volume> hash_volume;
    Unique_hash_map< Halffacet_handle, UFH_facet> hash_facet;
    Unique_hash_map< SFace_handle, UFH_sface> hash_sface;
    Union_find< Volume_handle> uf_volume;
    Union_find< Halffacet_handle> uf_facet;
    Union_find< SFace_handle> uf_sface;

    /* We discard  the information about boundary entry points, first
       on volumes, facets on sfacets.  Since during the volumes simplification
       is required the remotion of facet cycles, the information about those
       cycles is keep until the this simplification step is performed. */

    this->sncp()->clear_boundary();

    Volume_iterator c;
    CGAL_forall_volumes( c, *this->sncp()) {
      hash_volume[c] = uf_volume.make_set(c);
      this->sncp()->reset_object_list(c->shell_entry_objects());
    }
    SFace_iterator sf;
    CGAL_forall_sfaces( sf, *this->sncp()) {
      hash_sface[sf] = uf_sface.make_set(sf);
      this->sncp()->reset_sm_object_list(sf->boundary_entry_objects());
    }

    /* 
     * Volumes simplification 
     */

    Halffacet_handle f(D.halffacets_begin());
    while( f != D.halffacets_end() && f->is_twin())
      f++;
    while( f != D.halffacets_end()) {
      CGAL_assertion( !f->is_twin());
      Halffacet_iterator f_next(f);
      do
	f_next++;
      while( f_next != D.halffacets_end() && f_next->is_twin());
      CGAL_assertion( f != D.twin(f));
      Volume_handle c1 = D.volume(f), c2 = D.volume(D.twin(f));
      TRACEN(" mark("<<IO->index(c1)<<")="<<D.mark(c1)<<
      	     " mark("<<IO->index(f) <<")="<<D.mark(f) <<
	     " mark("<<IO->index(c2)<<")="<<D.mark(c2)<<
	     " is_twin(f)="<<f->is_twin());
      if( D.mark(c1) == D.mark(f) && D.mark(f) == D.mark(c2)
	  && D.is_standard(f)) {
	merge_sets( c1, c2, hash_volume, uf_volume);
	remove_f_including_all_edge_uses_in_its_boundary_cycles
	  (f, hash_sface, uf_sface);
	update_sfaces = update_volumes = true;
	TRACEN("UNION of "<<IO->index(c1)<<" & "<<IO->index(c2));
      }
      f = f_next;
    }
    
    CGAL_forall_halffacets( f, *this->sncp()) {
      hash_facet[f] = uf_facet.make_set(f);
      this->sncp()->reset_object_list(f->boundary_entry_objects());
    }
    
    /* 
     * Edges simplification
     */

    Halfedge_iterator e(D.halfedges_begin());
    while( e != D.halfedges_end() && e->is_twin())
      e++;
    while( e != (*this->sncp()).halfedges_end()) {
      CGAL_assertion( !e->is_twin());
      Halfedge_iterator e_next(e);
      do 
	e_next++;
      while( e_next != D.halfedges_end() && e_next->is_twin());
      SM_decorator SD(&*D.source(e));
      if( SD.is_isolated(e)) {
	if(D.mark(e) == D.mark(D.volume(D.sface(e)))) {
	  TRACEN("removing pair "<<IO->index(e)<<' '<<IO->index(D.twin(e)));
	  this->sncp()->delete_halfedge_pair(e);
	  update_facets = true;
	}
      } 
      else { 
	if( D.has_outdeg_two(e)) {
	  SHalfedge_handle e1(SD.first_out_edge(e)); 
	  SHalfedge_handle e2(SD.cyclic_adj_succ(e1));
	  if( SD.circle(e1)==SD.circle(SD.twin(e2)) &&
	      D.mark(e1)==D.mark(e) && D.mark(e)==D.mark(e2)) {
	    Halffacet_handle f1(D.facet(e1)); 
	    Halffacet_handle f2(D.facet(e2));
	    TRACEN("UNION of "<<IO->index(f1)<<" & "<<IO->index(D.twin(f2))<<
		   " ("<<IO->index(D.twin(f1))<<" & "<<IO->index(f2)<<")");
	    merge_sets( f1, D.twin(f2), hash_facet, uf_facet);
	    merge_sets( D.twin(f1), f2, hash_facet, uf_facet);
	    TRACEN("BEFORE "<<PFC(e1)<<std::endl<<PFC(D.twin(e2)));
	    TRACEN("removing "<<IO->index(e));
	    remove_edge_and_merge_facet_cycles(e);
	    update_facets = true;
	    // TRACEN("AFTER "<<PFC(e1)); // e1 not valid after sloop creation
	  }
	}
      }
      e = e_next;
    }

    update_facets = vertex_simplification() || update_facets;
    purge_no_find_objects(hash_volume, hash_facet, hash_sface, uf_volume, 
			  uf_facet, uf_sface);
    create_boundary_links_forall_sfaces( hash_sface, uf_sface);
    create_boundary_links_forall_facets( hash_facet, uf_facet);
    create_boundary_links_forall_volumes( hash_volume, uf_volume);
    
    TRACEN(">>> simplifying done ");

    return update_sfaces || update_facets || update_volumes;
  }
   
  void remove_edge_and_merge_facet_cycles( Halfedge_handle e) {
     SNC_decorator D(*this->sncp());
     CGAL_assertion( D.has_outdeg_two(e));
     Halfedge_handle et = D.twin(e);
     CGAL_assertion( D.has_outdeg_two(et));
     SM_decorator SD1(&*D.vertex(e));
     SM_decorator SD2(&*D.vertex(et));
     TRACEN("source " << IO->index(D.vertex(e)));
     TRACEN("target " << IO->index(D.vertex(et)));
     SHalfedge_handle e1 = SD1.first_out_edge(e);
     SHalfedge_handle e2 = SD2.next(D.previous(e1));
     merge_sedges_at_target_and_remove_svertex( D.twin(e1), e);
     merge_sedges_at_target_and_remove_svertex( D.twin(e2), et);
   }

   void merge_sedges_at_target_and_remove_svertex( SHalfedge_handle s1,
						   SVertex_handle v) {
     SNC_decorator D(*this->sncp());
     SM_decorator SD(&*D.vertex(v));
     CGAL_assertion( SD.target(s1) == v);
     SHalfedge_handle s2(SD.next(s1));
     CGAL_assertion( SD.source(s2) == v);
     TRACEN("s1 = " << IO->index(s1));
     TRACEN("s2 = " << IO->index(s2));
     if( s1 == s2) {
       TRACEN(IO->index(s1)<<'('<<IO->index(D.twin(s2))<<") to sloop");
       SD.convert_edge_to_loop(s1);
       CGAL_assertion(SD.shalfloop() != SHalfloop_handle());
       D.add_sloop_to_facet( SD.shalfloop(), D.facet(s1));
       TRACEN(IO->index(s2)<<" removed");
     }
     else {
       CGAL_assertion( D.has_outdeg_two(v));
       D.link_as_prev_next_pair( s1, D.next(s2));
       TRACEN(IO->index(s1)<<" "<<IO->index(D.next(s2))<<" linked.");
       D.link_as_prev_next_pair( D.twin(D.next(s2)), D.twin(s1));
       TRACEN(IO->index(D.twin(D.next(s2)))<<" "<<
	      IO->index(D.twin(s1))<<" linked.");
       SD.merge_edge_pairs_at_target( s1); // s2 is removed
       TRACEN(IO->index(s2)<<" removed");
     }
   }

   void merge_halfedge_pairs( SVertex_handle p, SVertex_handle q) {
     SNC_decorator D(*this->sncp());
     CGAL_assertion( D.vertex(p) == D.vertex(q));
     Vertex_handle v(D.vertex(p)); 
     CGAL_assertion( is_part_of_edge(v));
     SM_decorator SD(&*v);
     SHalfedge_around_svertex_circulator s(SD.first_out_edge(p)), se(s);
     CGAL_For_all( s, se) {
       D.link_as_prev_next_pair( D.previous(s), D.next(s));
       D.link_as_prev_next_pair( D.previous(SD.twin(s)), D.next(SD.twin(s)));
     }
     D.make_twins( D.twin(p), D.twin(q));
     SD.delete_vertex(p);
     SD.delete_vertex(q);
     this->sncp()->delete_vertex(v);
   }

   void purge_no_find_objects( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash_volume,
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash_facet,
      Unique_hash_map< SFace_handle, UFH_sface>& hash_sface,
      Union_find< Volume_handle>& uf_volume,
      Union_find< Halffacet_handle>& uf_facet,
      Union_find< SFace_handle>& uf_sface ) {
     SNC_decorator D(*this->sncp());
     SFace_iterator sf;
     std::list<SFace_handle> sflist;
     CGAL_forall_sfaces( sf, *this->sncp()) {
       if( uf_sface.find(hash_sface[sf]) != hash_sface[sf]) {
	 TRACEN("no find object "<<IO->index(sf));
	 sflist.push_back(sf);
       }
     }

     typename std::list<SFace_handle>::const_iterator sfli;
     for(sfli = sflist.begin(); sfli != sflist.end(); sfli++){     
       SM_decorator SD(&*D.vertex(*sfli));
       SD.delete_face_only(*sfli);
     }

     Halffacet_iterator f;
     std::list<Halffacet_handle> flist;
     CGAL_forall_facets( f, *this->sncp()) {
       TRACEN("facet "<<IO->index(f));
       if( uf_facet.find(hash_facet[f]) != hash_facet[f]) {
	 TRACEN("no find object "<<IO->index(f));
	 flist.push_back(f);
       }
     }
     
     typename std::list<Halffacet_handle>::const_iterator fli;
     for(fli = flist.begin(); fli != flist.end(); fli++)
       this->sncp()->delete_halffacet_pair(*fli);

     Volume_iterator c;
     std::list<Volume_handle> clist;
     CGAL_forall_volumes( c, *this->sncp()) {
       if( uf_volume.find(hash_volume[c]) != hash_volume[c]) {
	 TRACEN("no find object "<<IO->index(c));
	 clist.push_back(c);
       }
     }

     typename std::list<Volume_handle>::const_iterator cli;
     for(cli = clist.begin(); cli != clist.end(); cli++){     
       this->sncp()->delete_volume(*cli);
     }
   }

  void create_boundary_links_forall_sfaces(
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf ) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this->sncp());
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e, *this->sncp()) {
      if( linked[e])
	continue;
      SM_decorator SD(&*D.vertex(e));
      SFace_handle sf = *(uf.find(hash[D.sface(e)]));
      CGAL_assertion( sf != SFace_handle());
      SHalfedge_around_sface_circulator c(e), cend(c);
      CGAL_For_all( c, cend) {
	SD.set_face(c, sf);
	linked[c] = true;
      }
      SD.store_sm_boundary_object( e, sf);
    }

    SVertex_handle sv;
    CGAL_forall_svertices(sv, *this->sncp()) {
      SM_decorator SD(&*D.vertex(sv));
      if( SD.is_isolated(sv)) {
	SFace_handle sf = *(uf.find(hash[SD.face(sv)])); 
	CGAL_assertion( sf != SFace_handle());
	SD.set_face( sv, sf);
	SD.store_sm_boundary_object( sv, sf);
      }
    }

    SHalfloop_handle sl;
    CGAL_forall_shalfloops(sl, *this->sncp()) {
      SM_decorator SD(&*D.vertex(sl));
      SD.store_sm_boundary_object( sl, SD.face(sl));
    }
  }

  void create_boundary_links_forall_facets(
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash,
      Union_find< Halffacet_handle>& uf) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this->sncp());
    SHalfedge_iterator u;
    CGAL_forall_shalfedges(u, *this->sncp()) {
      if( linked[u])
	continue;
      /* set find(f) as incident facet of every edge use on the cycle of u */
      SHalfedge_handle u_min = u;
      Halffacet_handle f = *(uf.find(hash[D.facet(u)]));
      SHalfedge_around_facet_circulator c(u), cend(c);
      CGAL_For_all( c, cend) {
	D.set_facet( c, f);
	if( lexicographically_xyz_smaller(D.point(D.vertex(c)), D.point(D.vertex(u_min))))
	  u_min = c;
	linked[c] = true;
      }
      /* store the edge use at the lexicographicaly minimum facet vertex, as
	 a cycle entry of f.  The outermost cycle is stored at first
	 on the facet's cycles list. */
      if( is_empty_range( f->boundary_entry_objects().begin(),
			  f->boundary_entry_objects().end())) {
	D.store_boundary_object( u_min, f);
	TRACEN("new outer cycle min. vertex: "<<D.point(D.vertex(u_min)));
      }
      else {
	SHalfedge_handle f_sedge;
	CGAL_assertion( CGAL::assign( f_sedge, 
				     f->boundary_entry_objects().front()));
	CGAL::assign( f_sedge, f->boundary_entry_objects().front());
	if( lexicographically_xyz_smaller(D.point(D.vertex(u_min)), D.point(D.vertex(f_sedge))))
	  D.store_as_first_boundary_object( u_min, f);
	else
	  D.store_boundary_object( u_min, f);
      }
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops( l, *this->sncp()) {
      Halffacet_handle f = *(uf.find(hash[D.facet(l)]));
      D.set_facet( l, f);
      D.store_boundary_object( l, f);
    }
  }

  void create_boundary_links_forall_volumes( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash,
      Union_find< Volume_handle>& uf) {
    typedef typename SNC_decorator::template Shell_volume_setter<SNC_decorator> Volume_setter;
    //   typedef Unique_hash_map< SFace_handle, bool> SFace_map;
    //  SFace_map linked(false);

#ifdef CGAL_NEF_DEBUG
    delete IO;
    IO = new SNC_io_parser(std::cerr, *this->sncp());
#endif

    SNC_decorator D(*this->sncp());
    Volume_setter setter(D);

    SFace_iterator sf;
    Volume_handle c;
    CGAL_forall_sfaces(sf, *this->sncp()) {
      TRACEN("SFace " << IO->index(sf));
      if( setter.is_linked(sf)) continue;
      c = *(uf.find(hash[D.volume(sf)]));
      TRACEN("Volume " << IO->index(c));
      setter.set_volume(c);
      D.visit_shell_objects( sf, setter );      
      D.store_boundary_object( sf, c);
    }
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_SNC_STRUCTURE_H
