// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
#include <list>

#undef _DEBUG
#define _DEBUG 191
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template<typename S> class SM_const_decorator;

template <typename SNC_structure_>
class SNC_const_decorator { 
  typedef SNC_structure_                            Base;
  typedef SNC_structure_                            SNC_structure;
  typedef SNC_const_decorator<SNC_structure>        Self;
  typedef typename SNC_structure::Sphere_map        Sphere_map;
  typedef CGAL::SM_const_decorator<Sphere_map>      SM_const_decorator;
  const SNC_structure* sncp_;

  typedef typename SNC_structure::SHalfedge  SHalfedge;

public:
  typedef typename SNC_structure::Object_handle   Object_handle;
  typedef typename SNC_structure::Object_const_iterator Object_const_iterator;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::Volume_const_handle Volume_const_handle;
  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle;
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::Vertex_const_iterator Vertex_const_iterator;
  typedef typename SNC_structure::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename SNC_structure::Halffacet_const_iterator Halffacet_const_iterator; 
  typedef typename SNC_structure::Volume_const_iterator Volume_const_iterator;
  typedef typename SNC_structure::SVertex_const_iterator SVertex_const_iterator;
  typedef typename SNC_structure::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename SNC_structure::SHalfloop_const_iterator SHalfloop_const_iterator;
  typedef typename SNC_structure::SFace_const_iterator SFace_const_iterator;

  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Shell_entry_const_iterator Shell_entry_const_iterator;

  typedef typename SNC_structure::Kernel Kernel;
  typedef typename SNC_structure::FT FT;
  typedef typename SNC_structure::RT RT;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Line_3 Line_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Aff_transformation_3 Aff_transformation_3;

  typedef typename SNC_structure::Sphere_kernel Sphere_kernel;
  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;
  typedef typename SNC_structure::Sphere_direction Sphere_direction;

  typedef typename SNC_structure::Size_type Size_type;
  typedef typename SNC_structure::Mark Mark;
  typedef typename SNC_structure::Infi_box Infi_box;

  typedef typename SM_const_decorator::SHalfedge_around_svertex_const_circulator 
                                       SHalfedge_around_svertex_const_circulator;
  typedef typename SM_const_decorator::SHalfedge_around_sface_const_circulator 
                                       SHalfedge_around_sface_const_circulator;

 public:
  typedef void* GenPtr;

  SNC_const_decorator() : sncp_(0) {}
  SNC_const_decorator(const SNC_structure& W) : sncp_(&W) {}

 protected:  
  void set_snc(const SNC_structure& W) {
    sncp_ = &W;
  }

 public:
  const SNC_structure* sncp() const { return sncp_; }
  
  Vertex_const_handle vertex( Halfedge_const_handle e) const
  { return e->center_vertex(); }
  Halfedge_const_handle twin( Halfedge_const_handle e) const
  { return e->twin(); }
  Vertex_const_handle source( Halfedge_const_handle e) const
  { return e->center_vertex(); }
  Vertex_const_handle target( Halfedge_const_handle e) const
  { return source(twin(e)); }  
  SFace_const_handle sface( Halfedge_const_handle e) const
  { return e->incident_sface(); }
  /* SVertex queries*/

  Vertex_const_handle vertex(SHalfedge_const_handle e) const
  { return vertex(e->source()); }
  SHalfedge_const_handle twin(SHalfedge_const_handle e) const
  { return e->twin(); }
  Vertex_const_handle source(SHalfedge_const_handle e) const
  { return e->source()->center_vertex(); }
  Vertex_const_handle source(SHalfedge e) const
  { return e.source()->center_vertex(); }
  Vertex_const_handle target(SHalfedge_const_handle e) const
  { return e->twin()->source()->twin()->center_vertex(); }
  SHalfedge_const_handle previous(SHalfedge_const_handle e) const
  { return e->prev(); }
  SHalfedge_const_handle next(SHalfedge_const_handle e) const
  { return e->next(); }
  Halffacet_const_handle facet(SHalfedge_const_handle e) const
  { return e->incident_facet(); }
  SFace_const_handle sface(SHalfedge_const_handle e) const
  { return e->incident_sface(); }
  Halfedge_const_handle ssource(SHalfedge_const_handle e) const
  { return e->source(); }
  Halfedge_const_handle starget(SHalfedge_const_handle e) const
  { return e->twin()->source(); }
  /* SHalfedge queries */

  SHalfloop_const_handle twin( SHalfloop_const_handle l) const
  { return l->twin(); }
  Halffacet_const_handle facet( SHalfloop_const_handle l) const
  { return l->incident_facet(); }
  Vertex_const_handle vertex( SHalfloop_const_handle l) const
  { return l->incident_sface()->center_vertex(); }
  SFace_const_handle sface( SHalfloop_const_handle l) const
  { return l->incident_sface(); }
  /* SHalfloop queries */

  Vertex_const_handle vertex(SFace_const_handle f) const
  { return f->center_vertex(); }
  Volume_const_handle volume(SFace_const_handle f) const
  { return f->incident_volume(); }
  /* SFace queries */

  Halffacet_const_handle twin(Halffacet_const_handle f) const
  { return f->twin(); }
  Volume_const_handle volume(Halffacet_const_handle f) const
    { return f->volume(); }
  /* Halffacet queries */

  SFace_const_handle adjacent_sface(Halffacet_const_handle f) const {
    Halffacet_cycle_const_iterator fc(f->facet_cycles_begin());
    CGAL_assertion( fc != f->facet_cycles_end());
    SHalfedge_const_handle se;
    if ( assign(se, fc) ) { 
      CGAL_assertion( facet(se) == f);
      CGAL_assertion( sface(se) != SFace_const_handle());
      CGAL_assertion( volume(sface(twin(se))) == volume(f));
      return sface(twin(se));
    } 
    else 
      CGAL_assertion_msg( 0, "Facet outer cycle entry point"
			     "is not an SHalfedge? ");
    return SFace_const_handle(); // never reached
  }

  const Point_3& point(Vertex_const_handle v) const
  { return v->point(); }

  Vector_3 vector(Halfedge_const_handle e) const
  { return Vector_3(e->point_); }

  Segment_3 segment(Halfedge_const_handle e) const
  { return Segment_3(point(source(e)),
		     point(target(e))); }

  const Plane_3 plane(Halffacet_const_handle f) const
  { return f->plane(); }

  Mark mark(Vertex_const_handle v) const
  { return v->mark(); }
  Mark mark(Halfedge_const_handle e) const
  { return e->mark(); }
  Mark mark(Halffacet_const_handle f) const
  { return f->mark(); }
  Mark mark(Volume_const_handle c) const
  { return c->mark(); }

  template <typename Visitor>
  void visit_shell_objects(SFace_const_handle f, Visitor& V) const;

  Vertex_const_iterator   vertices_begin() const { 
    return sncp()->vertices_begin(); }
  Vertex_const_iterator   vertices_end()   const { 
    return sncp()->vertices_end(); }
  Halfedge_const_iterator halfedges_begin()const {
    return sncp()->halfedges_begin(); }
  Halfedge_const_iterator halfedges_end()  const { 
    return sncp()->halfedges_end(); }
  Halffacet_const_iterator halffacets_begin() const { 
    return sncp()->halffacets_begin(); }
  Halffacet_const_iterator halffacets_end() const { 
    return sncp()->halffacets_end(); }
  Volume_const_iterator   volumes_begin() const   { 
    return sncp()->volumes_begin(); }
  Volume_const_iterator   volumes_end()   const   { 
    return sncp()->volumes_end(); }

  Shell_entry_const_iterator shells_begin(Volume_const_handle c) const {
    return c->shells_begin();
  }
  Shell_entry_const_iterator shells_end(Volume_const_handle c) const {
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

  static bool is_standard(const Vertex_const_handle v) {
    return Infi_box::is_standard(v->point());
  }
  static bool is_standard(const Halffacet_const_handle f) {
    return Infi_box::is_standard(f->plane());
  }
  static bool is_standard_kernel() { return Infi_box::standard_kernel(); }
  static bool is_extended_kernel() { return Infi_box::extended_kernel(); }
  static void set_size_of_infimaximal_box(const typename Infi_box::NT& size) { 
    Infi_box::set_size_of_infimaximal_box(size); 
  }
};

template <typename EW>
template <typename Visitor>
void SNC_const_decorator<EW>::
visit_shell_objects(SFace_const_handle f, Visitor& V) const
{ 
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
      CGAL_forall_facet_cycles_of(fc,f) {
        SHalfedge_handle e;
	SHalfloop_handle l;
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
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
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
      SM_const_decorator SD(&*vertex(sf));
      /*      
      CGAL_forall_svertices(sv,SD){
	if(SD.is_isolated(sv) && !DoneSV[sv])
	  V.visit(sv);
      }
      */
      SFace_cycle_const_iterator fc;
      CGAL_forall_sface_cycles_of(fc,sf) {
        SVertex_handle v;
	SHalfedge_handle e;
	SHalfloop_handle l;
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
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
      }
    }
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONST_DECORATOR_H

