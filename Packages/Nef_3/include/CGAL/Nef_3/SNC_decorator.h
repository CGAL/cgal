// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_decorator.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan ert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_decorator.h         decorator for global SNC structure 
// ============================================================================
#ifndef CGAL_SNC_DECORATOR_H
#define CGAL_SNC_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_SM_decorator.h>

#undef _DEBUG
#define _DEBUG 13
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure_>
class SNC_decorator 
{ typedef SNC_structure_ SNC_structure;
  typedef SNC_decorator<SNC_structure_> Self;
  typedef CGAL::SNC_SM_decorator<SNC_structure_> SM_decorator;
  SNC_structure* sncp_;
public:
#define USING(t) typedef typename SNC_structure_::t t
  USING(Vertex_iterator);
  USING(Vertex_handle);
  USING(Vertex_const_handle);
  USING(Halfedge_iterator);
  USING(Halfedge_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_iterator); 
  USING(Halffacet_handle);
  USING(Halffacet_const_handle);
  USING(Volume_iterator);
  USING(Volume_handle);
  USING(Volume_const_handle);
  USING(SVertex_iterator);  
  USING(SVertex_handle);
  USING(SVertex_const_handle);
  USING(SHalfedge);
  USING(SHalfedge_iterator);
  USING(SHalfedge_handle);
  USING(SHalfedge_const_handle);
  USING(SHalfloop_iterator);
  USING(SHalfloop_handle);
  USING(SHalfloop_const_handle);
  USING(SFace_iterator);  
  USING(SFace_handle);
  USING(SFace_const_handle);
  USING(SHalfedge_const_iterator); 
  USING(Object_handle);
  USING(SObject_handle);
  USING(SFace_cycle_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Shell_entry_iterator);
  USING(Point_3);
  USING(Plane_3);
  USING(Segment_3);
  USING(Line_3);
  USING(Vector_3);
  USING(Direction_3);
  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Mark);
  USING(Size_type);
  USING(SHalfedge_around_facet_const_circulator);
  USING(SHalfedge_around_facet_circulator);
#undef USING
  typedef void* GenPtr;

  SNC_decorator() : sncp_() {}
  SNC_decorator(SNC_structure& W) : sncp_(&W) {}
  SNC_structure* sncp() const { return sncp_; }

  Vertex_handle vertex( Halfedge_handle e) const
  { return e->center_vertex_; }
  Vertex_const_handle vertex( Halfedge_const_handle e) const
  { return e->center_vertex_; }
  Halfedge_handle twin( Halfedge_handle e) const
  { return e->twin_; }
  Vertex_handle source( Halfedge_handle e) const
  { return e->center_vertex_; }
  Vertex_handle target( Halfedge_handle e) const
  { return source(twin(e)); }
  SFace_handle sface( Halfedge_handle e) const
  { return e->incident_sface_; }
  SFace_const_handle sface( Halfedge_const_handle e) const
  { return e->incident_sface_; }
  /* SVertex queries*/

  Vertex_handle vertex(SHalfedge_handle e) const
  { return vertex(e->source_); }
  SHalfedge_handle twin(SHalfedge_handle e) const
  { return e->twin_; }
  Vertex_handle source(SHalfedge_handle e) const
  { return e->source_->center_vertex_; }
  Vertex_handle source(SHalfedge e) const
  { return e.source_->center_vertex_; }
  Vertex_handle target(SHalfedge_handle e) const
  { return e->twin_->source_->twin_->center_vertex_; }
  SHalfedge_handle previous(SHalfedge_handle e) const
  { return e->prev_; }
  SHalfedge_handle next(SHalfedge_handle e) const
  { return e->next_; }
  Halffacet_handle facet(SHalfedge_handle e) const
  { return e->incident_facet_; }
  SFace_handle sface(SHalfedge_handle e) const
  { return e->incident_sface_; }
  Halfedge_handle ssource(SHalfedge_handle e) const
  { return e->source_; }
  Halfedge_handle starget(SHalfedge_handle e) const
  { return e->twin_->source_; }
  /* SHalfedge queries */

  const Mark& mark(SHalfedge_handle e) const
  /*{\Mop returns the mark associated to |e| as
  a sphere object. This is temporary information!!!}*/
  { return e->tmp_mark(); }

  const Sphere_circle& tmp_circle(SHalfedge_handle e) const
  { return e->tmp_circle(); }
  const Sphere_circle& tmp_circle(SHalfloop_handle l) const
  { return l->tmp_circle(); }

  std::string debug(SHalfedge_handle e) const
  { std::ostrstream os; set_pretty_mode(os);
    os << "sedge-use " << point(source(e)) << point(target(e))<<'\0';
    std::string res(os.str()); os.freeze(0); return res; 
  }

  SHalfloop_handle twin( SHalfloop_handle l) const
  { return l->twin_; }
  Halffacet_handle facet( SHalfloop_handle l) const
  { return l->incident_facet_; }
  Vertex_handle vertex( SHalfloop_handle l) const
  { return l->incident_sface_->center_vertex_; }
  SFace_handle sface( SHalfloop_handle l) const
  { return l->incident_sface_; }
  SFace_const_handle sface( SHalfloop_const_handle l) const
  { return l->incident_sface_; }
  /* SHalfloop queries */

  Vertex_handle vertex(SFace_handle f) const
  { return f->center_vertex_; }
  Volume_handle volume(SFace_handle f) const
  { return f->incident_volume_; }
  /* SHalffacet queries */

  Halffacet_handle twin(Halffacet_handle f) const
  { return f->twin_; }
  Volume_handle volume(Halffacet_handle f) const
    { return f->volume_; }
  Volume_const_handle volume(Halffacet_const_handle f) const
    { return f->volume_; }
  /* Halffacet queries */

  SFace_handle adjacent_sface(Halffacet_handle f) const {
    Halffacet_cycle_iterator fc(f->facet_cycles_begin());
    CGAL_assertion( fc != f->facet_cycles_end());
    SHalfedge_handle se;
    if ( assign(se, fc) ) { 
      CGAL_assertion( facet(se) == f);
      CGAL_assertion( sface(se) != SFace_handle());
      CGAL_assertion( volume(sface(twin(se))) == volume(f));
      return sface(twin(se));
    } 
    else 
      CGAL_assertion_msg( 0, "Facet outer cycle entry point"
			     "is not an SHalfedge? ");
    return SFace_handle(); // never reached
  }

  // attributes::
  Point_3& point(Vertex_handle v) const
  { return v->point(); }
  const Point_3& point(Vertex_const_handle v) const
  { return v->point(); }

  Sphere_point tmp_point(Halfedge_handle e) const
  { return e->tmp_point(); }
  Sphere_point calc_point(Halfedge_handle e) const
  { CGAL_nef3_assertion(twin(e)!=Halfedge_handle());
    Point_3 ps(point(source(e)));
    Point_3 pt(point(target(e)));
    return Sphere_point(pt-ps);
  }

  Segment_3 segment(Halfedge_handle e) const
  { return Segment_3(point(source(e)),
		     point(target(e))); }

  Plane_3& plane(Halffacet_handle f) const
  { return f->plane(); }

  Vector_3 orthogonal_vector(Halffacet_handle f) const
  { return f->plane().orthogonal_vector(); }

  Mark& mark(Vertex_handle v) const
  { return v->mark(); }
  Mark& mark(Halfedge_handle e) const
  { return e->mark(); }
  Mark& mark(Halffacet_handle f) const
  { return f->mark(); }
  Mark& mark(Volume_handle c) const
  { return c->mark(); }

  GenPtr& info(Vertex_handle v) const
  { return v->info(); }

  Mark& tmp_mark(Halfedge_handle e) const
  { return e->tmp_mark(); }

  template <typename H>
  bool is_boundary_object(H h) const
  { return sncp()->is_boundary_object(h); }

  template <typename H>
  void store_boundary_object(H h, Halffacet_handle f) const
  { f->boundary_entry_objects_.push_back(Object_handle(h));
    sncp()->store_boundary_item(h, --(f->facet_cycles_end()));
  }

  template <typename H>
  void undo_boundary_object(H h, Halffacet_handle f) const
  { CGAL_assertion(sncp()->is_boundary_object(h));
    Halffacet_cycle_iterator it = sncp()->boundary_item(h);
    sncp()->undef_boundary_item(h);
    f->boundary_entry_objects_.erase(it);
  }

  void link_as_facet_cycle(SHalfedge_handle e, Halffacet_handle f) const
  { SHalfedge_around_facet_circulator hfc(e), hend(hfc);
    CGAL_For_all(hfc,hend) hfc->incident_facet_ = f;
    store_boundary_object(e,f);
  } 

  void link_as_interior_loop(SHalfloop_handle l, Halffacet_handle f) const
  { l->incident_facet_ = f;
    store_boundary_object(l,f);
  } 

  template <typename H>
  void make_twins(H h1, H h2) const
  { h1->twin_ = h2; h2->twin_ = h1; }

  void link_as_prev_next_pair(SHalfedge_handle e1, SHalfedge_handle e2) const
  { e1->next_ = e2; e2->prev_ = e1; } 

  void store_boundary_object(SFace_handle h, Volume_handle c) const
  { c->shell_entry_objects_.push_back(Object_handle(h));
    sncp()->store_boundary_item(h, --(c->shells_end()));
  }

  struct Shell_volume_setter {
    const SNC_decorator D;
    Volume_handle c;
    Shell_volume_setter(const SNC_decorator& Di, Volume_handle ci)
      : D(Di), c(ci) {}
    void visit(SFace_handle h) { D.set_volume(h, c); }
    void visit(Vertex_handle h) { /* empty */ }
    void visit(Halfedge_handle h) { /* empty */ }
    void visit(Halffacet_handle h ) { D.set_volume(h, c); }
    void set_volume(Volume_handle ci) { c = ci; }
  };

  void link_as_outer_shell( SFace_handle f, Volume_handle c ) const {
    CGAL_assertion(c->shell_entry_objects_.size() == 0);
    Shell_volume_setter Setter(*this, c);
    visit_shell_objects( f, Setter );
    TRACEN("Volume "<<&*c<<", outer shell "<<&*f);
    store_boundary_object( f, c );
  }

  void link_as_inner_shell( SFace_handle f, Volume_handle c ) const {
    // CGAL_assertion(c->shell_entry_objects_.size() > 0);
    Shell_volume_setter Setter(*this, c);
    visit_shell_objects( f, Setter );
    TRACEN("Volume "<<&*c<<", inner shell "<<&*f);
    store_boundary_object( f, c );
  }

  void set_volume(Halffacet_handle h, Volume_handle c) const
  { h->volume_ = c; CGAL_assertion(h->volume_ == c); }
  void set_volume(SFace_handle h, Volume_handle c) const 
  { h->incident_volume_ = c; CGAL_assertion(h->incident_volume_ == c); }

  template <typename Visitor>
  void visit_shell_objects(SFace_handle f, Visitor& V) const;

  Vertex_iterator   vertices_begin()   { return sncp()->vertices_begin(); }
  Vertex_iterator   vertices_end()     { return sncp()->vertices_end(); }
  Halfedge_iterator halfedges_begin()  { return sncp()->halfedges_begin(); }
  Halfedge_iterator halfedges_end()    { return sncp()->halfedges_end(); }
  Halffacet_iterator    halffacets_begin()     { return sncp()->halffacets_begin(); }
  Halffacet_iterator    halffacets_end()       { return sncp()->halffacets_end(); }
  Volume_iterator   volumes_begin()    { return sncp()->volumes_begin(); }
  Volume_iterator   volumes_end()      { return sncp()->volumes_end(); }

  Shell_entry_iterator  shells_begin(Volume_handle c) {
	return sncp()->shells_begin(c);
  }
  Shell_entry_iterator  shells_end(Volume_handle c) {
	return sncp()->shells_end(c);
  }

  Size_type number_of_vertices() const  
  { return sncp()->number_of_vertices(); }
  Size_type number_of_halfedges() const 
  { return sncp()->number_of_halfedges(); }
  Size_type number_of_edges() const     
  { return sncp()->number_of_edges(); }
  Size_type number_of_halffacets() const    
  { return sncp()->number_of_halffacets();}
  Size_type number_of_facets() const    
  { return sncp()->number_of_facets();}
  Size_type number_of_volumes() const   
  { return sncp()->number_of_volumes();}

};




/* visiting shell objects:

Objects are marked as done, when placed in the output list.  We have
to maintain a stack of sface candidates (the spherical rubber sectors
that provide connectivity at the local graphs of vertices) and facet
candiates (the plane pieces in three space also providing
connectivity). Note that we have to take care about the orientation of
sobjects and facets. We have to take care that (1) the search along
the shell extends along the whole shell structure (2) does not visit
any object twice, and (3) all 3-space objects have to be reported and
this also just once.

The facets and sfaces are marked |done| when they are put into their
corresponding queues thus each such object is visited exactly once
when taken out of the queue. 

When an sface |sf| is taken out of the queue |SFaceCandiates| its
boundary structure is examined and all 2-skeleton objects (vertices
and edges of 3-space) that are incident to the volume represented by
|sf| are reported. Facets are reported when they are taken out of
|FacetCandiates|.

*/

template <typename EW>
template <typename Visitor>
void SNC_decorator<EW>::
visit_shell_objects(SFace_handle f, Visitor& V) const
{ 
  typedef SM_decorator::SHalfedge_around_sface_circulator 
    SHalfedge_around_sface_circulator;
  std::list<SFace_handle> SFaceCandidates;
  std::list<Halffacet_handle> FacetCandidates;
  CGAL::Generic_handle_map<bool> Done(false);
  SFaceCandidates.push_back(f);  Done[f] = true;
  while ( true ) {
    if ( SFaceCandidates.empty() && FacetCandidates.empty() ) break;
    if ( !FacetCandidates.empty() ) {
      Halffacet_handle f = *FacetCandidates.begin();
      FacetCandidates.pop_front();
      V.visit(f); // report facet
      Halffacet_cycle_iterator fc;
      CGAL_nef3_forall_facet_cycles_of(fc,f) {
        SHalfedge_handle e; SHalfloop_handle l;
        if ( assign(e,fc) ) { 
          SHalfedge_around_facet_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) { e = twin(ec);
            if ( Done[sface(e)] ) continue;
            SFaceCandidates.push_back(sface(e));
            Done[sface(e)] = true;
          }
        } else if ( assign(l,fc) ) { l = twin(l);
          if ( Done[sface(l)] ) continue;
          SFaceCandidates.push_back(sface(l));
          Done[sface(l)] = true;
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
      }
    }
    if ( !SFaceCandidates.empty() ) {
      SFace_handle sf = *SFaceCandidates.begin();
      SFaceCandidates.pop_front();
      V.visit(sf);
      if ( !Done[vertex(sf)] )
        V.visit(vertex(sf)); // report vertex
      SFace_cycle_iterator fc;
      CGAL_nef3_forall_sface_cycles_of(fc,sf) {
        SVertex_handle v; SHalfedge_handle e; SHalfloop_handle l;
        if ( assign(e,fc) ) {
          SHalfedge_around_sface_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) { 
            v = starget(ec);
            if ( !Done[v] ) {
              V.visit(v); // report edge
              Done[v] = Done[twin(v)] = true;
            }
            Halffacet_handle f = facet(twin(ec));
            if ( Done[f] ) continue;
            FacetCandidates.push_back(f); Done[f] = true;
          }
        } else if ( assign(v,fc) ) {
          if ( Done[v] ) continue; 
          V.visit(v); // report edge
          Done[v] = Done[twin(v)] = true;
          // note that v is isolated, thus twin(v) is isolated too
          SFace_handle fo = sface(twin(v));
          if ( Done[fo] ) continue;
          SFaceCandidates.push_back(fo); Done[fo] = true;
        } else if ( assign(l,fc) ) {
          Halffacet_handle f = facet(twin(l));
          if ( Done[f] ) continue;
          FacetCandidates.push_back(f);  Done[f] = true;
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
      }
    }
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_DECORATOR_H

