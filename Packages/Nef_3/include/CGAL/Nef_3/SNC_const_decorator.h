// Copyright (c) 1997-2002  Max-Planck-Institute Saarbrucken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_CONST_DECORATOR_H
#define CGAL_SNC_CONST_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/Normalizing.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_io_parser.h>
#undef _DEBUG
#define _DEBUG 191
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure_>
class SNC_const_decorator { 
  typedef SNC_structure_                            Base;
  typedef SNC_structure_                            SNC_structure;
  typedef SNC_const_decorator<SNC_structure>        Self;
  typedef SNC_SM_const_decorator<SNC_structure>     SM_const_decorator;
  const SNC_structure* sncp_;

public:
#define USING(t) typedef typename SNC_structure_::t t
  USING(Vertex_const_iterator);
  USING(Vertex_const_handle);
  USING(Halfedge_const_iterator);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_iterator); 
  USING(Halffacet_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_iterator);
  USING(Volume_handle);
  USING(Volume_const_handle);
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
  USING(SHalfedge_around_facet_const_circulator);
  USING(SHalfedge_around_facet_circulator);
  USING(SFace_cycle_const_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Halffacet_cycle_const_iterator);
  USING(Shell_entry_iterator);
  USING(Shell_entry_const_iterator);
  USING(Kernel);
  USING(FT);
  USING(RT);
  USING(Point_3);
  USING(Plane_3);
  USING(Segment_3);
  USING(Line_3);
  USING(Vector_3);
  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Mark);
  USING(Size_type);
  USING(Infi_box);
#undef USING

#define DECUSING(t) typedef typename SM_const_decorator::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
#undef DECUSING

 public:
  typedef void* GenPtr;

  SNC_const_decorator() : sncp_(0) {}
  SNC_const_decorator(const SNC_structure& W) : sncp_(&W) {}
  const SNC_structure* sncp() const { return sncp_; }

  Vertex_const_handle vertex( Halfedge_const_handle e) const
  { return e->center_vertex_; }
  Halfedge_const_handle twin( Halfedge_const_handle e) const
  { return e->twin_; }
  Vertex_const_handle source( Halfedge_const_handle e) const
  { return e->center_vertex_; }
  Vertex_const_handle target( Halfedge_const_handle e) const
  { return source(twin(e)); }  
  SFace_const_handle sface( Halfedge_const_handle e) const
  { return e->incident_sface_; }
  /* SVertex queries*/

  Vertex_const_handle vertex(SHalfedge_const_handle e) const
  { return vertex(e->source_); }
  SHalfedge_const_handle twin(SHalfedge_const_handle e) const
  { return e->twin_; }
  Vertex_const_handle source(SHalfedge_const_handle e) const
  { return e->source_->center_vertex_; }
  Vertex_const_handle source(SHalfedge e) const
  { return e.source_->center_vertex_; }
  Vertex_const_handle target(SHalfedge_const_handle e) const
  { return e->twin_->source_->twin_->center_vertex_; }
  SHalfedge_const_handle previous(SHalfedge_const_handle e) const
  { return e->prev_; }
  SHalfedge_const_handle next(SHalfedge_const_handle e) const
  { return e->next_; }
  Halffacet_const_handle facet(SHalfedge_const_handle e) const
  { return e->incident_facet_; }
  SFace_const_handle sface(SHalfedge_const_handle e) const
  { return e->incident_sface_; }
  Halfedge_const_handle ssource(SHalfedge_const_handle e) const
  { return e->source_; }
  Halfedge_const_handle starget(SHalfedge_const_handle e) const
  { return e->twin_->source_; }
  /* SHalfedge queries */

  SHalfloop_const_handle twin( SHalfloop_const_handle l) const
  { return l->twin_; }
  Halffacet_const_handle facet( SHalfloop_const_handle l) const
  { return l->incident_facet_; }
  Vertex_const_handle vertex( SHalfloop_const_handle l) const
  { return l->incident_sface_->center_vertex_; }
  SFace_const_handle sface( SHalfloop_const_handle l) const
  { return l->incident_sface_; }
  /* SHalfloop queries */

  Vertex_const_handle vertex(SFace_const_handle f) const
  { return f->center_vertex_; }
  Volume_const_handle volume(SFace_const_handle f) const
  { return f->incident_volume_; }
  /* SHalffacet queries */

  Halffacet_const_handle twin(Halffacet_const_handle f) const
  { return f->twin_; }
  Volume_const_handle volume(Halffacet_const_handle f) const
    { return f->volume_; }
  /* Halffacet queries */

  SFace_const_handle adjacent_sface(Halffacet_const_handle f) const {
    Halffacet_cycle_const_iterator fc(f->facet_cycles_begin());
    CGAL_nef3_assertion( fc != f->facet_cycles_end());
    SHalfedge_const_handle se;
    if ( assign(se, fc) ) { 
      CGAL_nef3_assertion( facet(se) == f);
      CGAL_nef3_assertion( sface(se) != SFace_const_handle());
      CGAL_nef3_assertion( volume(sface(twin(se))) == volume(f));
      return sface(twin(se));
    } 
    else 
      CGAL_nef3_assertion_msg( 0, "Facet outer cycle entry point"
			     "is not an SHalfedge? ");
    return SFace_const_handle(); // never reached
  }

  // attributes::
  const Point_3& point(Vertex_const_handle v) const
  { return v->point(); }

  Sphere_point tmp_point(Halfedge_const_handle e) const
  { return e->tmp_point(); }
  Sphere_point calc_point(Halfedge_const_handle e) const
  { CGAL_nef3_assertion(twin(e)!=Halfedge_const_handle());
    Point_3 ps(point(source(e)));
    Point_3 pt(point(target(e)));
    return Sphere_point(pt-ps);
  }

  Segment_3 segment(Halfedge_const_handle e) const
  { return Segment_3(point(source(e)),
		     point(target(e))); }

  const Plane_3 plane(Halffacet_const_handle f) const
  { return f->plane(); }

  Vector_3 orthogonal_vector(Halffacet_const_handle f) const
  { return f->plane().orthogonal_vector(); }

  Mark mark(Vertex_const_handle v) const
  { return v->mark(); }
  Mark mark(Halfedge_const_handle e) const
  { return e->mark(); }
  Mark mark(Halffacet_const_handle f) const
  { return f->mark(); }
  Mark mark(Volume_const_handle c) const
  { return c->mark(); }

  template <typename H>
  bool is_boundary_object(H h) const
  { return sncp()->is_boundary_object(h); }

  bool is_infbox_vertex( Vertex_const_handle v) const {
    return !Infi_box::is_standard(v->point());
  }

  /* returns true if |f| is part of the infinimaximal box.*/
  bool is_infbox_facet(Halffacet_const_handle f) const {
    return !Infi_box::is_standard(f->plane());
  }

  template <typename Visitor>
  void visit_shell_objects(SFace_const_handle f, Visitor& V) const;

  Vertex_const_iterator   vertices_begin()   { return sncp()->vertices_begin(); }
  Vertex_const_iterator   vertices_end()     { return sncp()->vertices_end(); }
  Halfedge_const_iterator halfedges_begin()  { return sncp()->halfedges_begin(); }
  Halfedge_const_iterator halfedges_end()    { return sncp()->halfedges_end(); }
  Halffacet_const_iterator halffacets_begin(){ return sncp()->halffacets_begin(); }
  Halffacet_const_iterator halffacets_end()  { return sncp()->halffacets_end(); }
  Volume_const_iterator   volumes_begin()    { return sncp()->volumes_begin(); }
  Volume_const_iterator   volumes_end()      { return sncp()->volumes_end(); }

  Shell_entry_const_iterator shells_begin(Volume_const_handle c) {
    return c->shells_begin();
  }
  Shell_entry_const_iterator shells_end(Volume_const_handle c) {
    return c->shells_end();
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

template <typename EW>
template <typename Visitor>
void SNC_const_decorator<EW>::
visit_shell_objects(SFace_const_handle f, Visitor& V) const
{ 
  typedef typename SM_const_decorator::SHalfedge_around_sface_const_circulator 
    SHalfedge_around_sface_const_circulator;
  std::list<SFace_const_handle> SFaceCandidates;
  std::list<Halffacet_const_handle> FacetCandidates;
  CGAL::Unique_hash_map<SFace_const_handle,bool> DoneSF(false);
  CGAL::Unique_hash_map<Vertex_const_handle,bool> DoneV(false);
  CGAL::Unique_hash_map<SVertex_const_handle,bool> DoneSV(false);
  CGAL::Unique_hash_map<Halffacet_const_handle,bool> DoneF(false);
  SFaceCandidates.push_back(f);  DoneSF[f] = true;
  while ( true ) {
    if ( SFaceCandidates.empty() && FacetCandidates.empty() ) break;
    if ( !FacetCandidates.empty() ) {
      Halffacet_const_handle f = *FacetCandidates.begin();
      FacetCandidates.pop_front();
      V.visit(f); // report facet
      Halffacet_cycle_const_iterator fc;
      CGAL_nef3_forall_facet_cycles_of(fc,f) {
        SHalfedge_handle e; SHalfloop_handle l;
        if ( assign(e,fc) ) {
	  SHalfedge_const_handle she;
          SHalfedge_around_facet_const_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) { she = twin(ec);
            if ( DoneSF[sface(she)] ) continue;
            SFaceCandidates.push_back(sface(she));
            DoneSF[sface(she)] = true;
          }
        } else if ( assign(l,fc) ) { 
	  SHalfloop_const_handle ll = twin(l);
          if ( DoneSF[sface(ll)] ) continue;
          SFaceCandidates.push_back(sface(ll));
          DoneSF[sface(ll)] = true;
        } else CGAL_nef3_assertion_msg(0,"Damn wrong handle.");
      }
    }
    if ( !SFaceCandidates.empty() ) {
      SFace_const_handle sf = *SFaceCandidates.begin();
      SFaceCandidates.pop_front();
      V.visit(sf);
      if ( !DoneV[vertex(sf)] )
        V.visit(vertex(sf)); // report vertex
      DoneV[vertex(sf)] = true;
      //      SVertex_const_handle sv;
      SM_const_decorator SD(vertex(sf));
      /*      
      CGAL_nef3_forall_svertices(sv,SD){
	if(SD.is_isolated(sv) && !DoneSV[sv])
	  V.visit(sv);
      }
      */
      SFace_cycle_const_iterator fc;
      CGAL_nef3_forall_sface_cycles_of(fc,sf) {
        SVertex_handle v; SHalfedge_handle e; SHalfloop_handle l;
        if ( assign(e,fc) ) {
	  SHalfedge_around_sface_const_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) { 
            SVertex_const_handle vv = starget(ec);
            if ( !SD.is_isolated(vv) && !DoneSV[vv] ) {
              V.visit(vv); // report edge
              DoneSV[vv] = DoneSV[twin(vv)] = true;
            }
            Halffacet_const_handle f = facet(twin(ec));
            if ( DoneF[f] ) continue;
            FacetCandidates.push_back(f); DoneF[f] = true;
          }
        } else if ( assign(v,fc) ) {
          if ( DoneSV[v] ) continue; 
          V.visit(v); // report edge
	  V.visit(twin(v));
          DoneSV[v] = DoneSV[twin(v)] = true;
	  CGAL_assertion(SD.is_isolated(v));
	  SFaceCandidates.push_back(sface(twin(v)));
	  DoneSF[sface(twin(v))]=true;
          // note that v is isolated, thus twin(v) is isolated too
	  //	  SM_const_decorator SD;
	  //	  SFace_const_handle fo;
	  //	  fo = sface(twin(v));
	  /*
	  if(SD.is_isolated(v)) 
	    fo = source(v)->sfaces_begin();
	  else
	    fo = sface(twin(v));
	  */
        } else if ( assign(l,fc) ) {
          Halffacet_const_handle f = facet(twin(l));
          if ( DoneF[f] ) continue;
          FacetCandidates.push_back(f);  DoneF[f] = true;
        } else CGAL_nef3_assertion_msg(0,"Damn wrong handle.");
      }
    }
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONST_DECORATOR_H

