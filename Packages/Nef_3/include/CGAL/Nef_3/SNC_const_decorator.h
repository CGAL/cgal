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
#include <CGAL/Nef_3/SNC_decorator_traits.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
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
  typedef SNC_decorator_const_traits<SNC_structure>  Decorator_traits;

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
  SNC_const_decorator(const SNC_const_decorator& S) {
    sncp_ = S.sncp_;
  }

 protected:  
  void set_snc(const SNC_structure& W) {
    sncp_ = &W;
  }

 public:
  const SNC_structure* sncp() const { return sncp_; }
  
  static Vertex_const_handle vertex( Halfedge_const_handle e)
  { return e->center_vertex(); }
  static Halfedge_const_handle twin( Halfedge_const_handle e)
  { return e->twin(); }
  static Vertex_const_handle source( Halfedge_const_handle e)
  { return e->center_vertex(); }
  static Vertex_const_handle target( Halfedge_const_handle e)
  { return source(twin(e)); }  
  static SFace_const_handle sface( Halfedge_const_handle e)
  { return e->incident_sface(); }
  /* SVertex queries*/

  static Vertex_const_handle vertex(SHalfedge_const_handle e)
  { return vertex(e->source()); }
  static SHalfedge_const_handle twin(SHalfedge_const_handle e)
  { return e->twin(); }
  static Vertex_const_handle source(SHalfedge_const_handle e)
  { return e->source()->center_vertex(); }
  static Vertex_const_handle source(SHalfedge e)
  { return e.source()->center_vertex(); }
  static Vertex_const_handle target(SHalfedge_const_handle e)
  { return e->twin()->source()->twin()->center_vertex(); }
  static SHalfedge_const_handle previous(SHalfedge_const_handle e)
  { return e->prev(); }
  static SHalfedge_const_handle next(SHalfedge_const_handle e)
  { return e->next(); }
  static Halffacet_const_handle facet(SHalfedge_const_handle e)
  { return e->incident_facet(); }
  static SFace_const_handle sface(SHalfedge_const_handle e)
  { return e->incident_sface(); }
  static Halfedge_const_handle ssource(SHalfedge_const_handle e)
  { return e->source(); }
  static Halfedge_const_handle starget(SHalfedge_const_handle e)
  { return e->twin()->source(); }
  /* SHalfedge queries */

  static SHalfloop_const_handle twin( SHalfloop_const_handle l)
  { return l->twin(); }
  static Halffacet_const_handle facet( SHalfloop_const_handle l)
  { return l->incident_facet(); }
  static Vertex_const_handle vertex( SHalfloop_const_handle l)
  { return l->incident_sface()->center_vertex(); }
  static SFace_const_handle sface( SHalfloop_const_handle l)
  { return l->incident_sface(); }
  /* SHalfloop queries */

  static Vertex_const_handle vertex(SFace_const_handle f)
  { return f->center_vertex(); }
  static Volume_const_handle volume(SFace_const_handle f)
  { return f->incident_volume(); }
  /* SFace queries */

  static Halffacet_const_handle twin(Halffacet_const_handle f)
  { return f->twin(); }
  static Volume_const_handle volume(Halffacet_const_handle f)
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

  static const Point_3& point(Vertex_const_handle v)
  { return v->point(); }

  static Vector_3 vector(Halfedge_const_handle e)
  { return Vector_3(e->point()-CGAL::ORIGIN); }

  static Segment_3 segment(Halfedge_const_handle e)
  { return Segment_3(point(source(e)),
		     point(target(e))); }

  static const Plane_3 plane(Halffacet_const_handle f)
  { return f->plane(); }

  static Mark mark(Vertex_const_handle v)
  { return v->mark(); }
  static Mark mark(Halfedge_const_handle e)
  { return e->mark(); }
  static Mark mark(Halffacet_const_handle f)
  { return f->mark(); }
  static Mark mark(Volume_const_handle c)
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

  bool is_bounded() const {
    if(is_standard_kernel())
      return true;
    int i = 0;
    Halffacet_const_handle hf;
    CGAL_forall_facets(hf, *sncp()) {
      if(!Infi_box::is_standard(plane(hf)))
	++i;
    }
    CGAL_assertion(i>=6);
    return (i == 6);
  }

  static bool is_bounded(Halffacet_const_handle f) {
    Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
    CGAL_assertion(fc.is_shalfedge());
    SHalfedge_const_handle sh(fc);
    SHalfedge_around_facet_const_circulator fcc(sh), fend(fcc);
    CGAL_For_all(fcc,fend)
      if(!is_standard(fcc->source()->source()))
	return false;
    return true;
  }

  static bool is_standard(Vertex_const_handle v) {
    return Infi_box::is_standard(v->point());
  }
  static bool is_standard(Halffacet_const_handle f) {
    return Infi_box::is_standard(f->plane());
  }
  static bool is_standard_kernel() { return Infi_box::standard_kernel(); }
  static bool is_extended_kernel() { return Infi_box::extended_kernel(); }
  static void set_size_of_infimaximal_box(const typename Infi_box::NT& size) { 
    Infi_box::set_size_of_infimaximal_box(size); 
  }

  typedef SM_point_locator<SM_const_decorator>         SM_point_locator;

  Halffacet_const_handle get_visible_facet( const Vertex_const_handle v, 
				      const Ray_3& ray) const
    /*{\Mop when one shoot a ray |ray| in order to find the facet below to
      an object, and vertex |v| is hit, we need to choose one of the facets
      in the adjacency list of |v| such that it could be 'seen' from the
      piercing point of the |ray| on the sphere map on |v|.  We make it just
      locating the sphere facet |sf| pierced by |ray| and taking the adjacent 
      facet to one of the sphere segments on the boundary of |sf|.
      \precondition |ray| target is on |v| and the intersection between
      |ray| and the 2-skeleton incident to v is empty. }*/ {

    Halffacet_const_handle f_visible;
    CGAL_assertion( ray.source() != point(v));
    CGAL_assertion( ray.has_on(point(v)));
    Sphere_point sp(ray.source() - point(v));
    TRACEN( "Locating "<<sp <<" in "<<point(v));
    CGAL_assertion(Infi_box::degree(sp.hx()) < 2 && 
		   Infi_box::degree(sp.hy()) < 2 && 
		   Infi_box::degree(sp.hz()) < 2 && 
		   Infi_box::degree(sp.hw()) == 0);
    sp = Infi_box::simplify(sp);
    TRACEN( "Locating "<<sp <<" in "<<point(v));
    SM_point_locator L(&*v);
    Object_handle o = L.locate(sp);

    SFace_const_handle sf;
    //    CGAL_assertion(assign(sf,o));
    //    assign(sf,o);
    if(!assign( sf, o)) {
      CGAL_assertion_msg( 0, "it is not possible to decide which one is a visible facet (if any)");
      return Halffacet_const_handle();
    }
    /*    SM_decorator SD;
    if(sncp()->halfedges_begin() == sncp()->halfedges_end() || 
       SD.is_isolated(sncp()->halfedges_begin())) 
      return Halffacet_const_handle();  
    */

    SFace_cycle_const_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) {
	TRACEN( "no adjacent facet found.");
	f_visible =  Halffacet_const_handle();
    }
    else {
      SHalfedge_const_handle se; 
      SHalfloop_const_handle sl;
      SVertex_const_handle sv;
      if (fc.is_shalfedge()) {
	se = SHalfedge_const_handle(fc);
	TRACEN( "adjacent facet found (SEdges cycle).");
	TRACEN("se"<<PH(se));
	SM_const_decorator SD;
	TRACEN(plane(facet(se))<<"/"<<plane(facet(SD.next(se)))<<"/"<<plane(facet(SD.next(SD.next(se)))));
	f_visible = facet(twin(se));
	TRACEN("f_visible"<<plane(f_visible));
      }
      else if (fc.is_shalfloop()) {
	sl = SHalfloop_const_handle(sl);
	SM_const_decorator SD;
	TRACEN( "adjacent facet found (SHalfloop cycle)."<<SD.circle(sl) 
		<< " with facet "<<plane(facet(sl)));
	f_visible = facet(twin(sl));
	TRACEN("f_visible"<<plane(f_visible));
      }
      else if(fc.is_svertex()) {
	sv = SVertex_const_handle(fc);
#ifdef _DEBUG
	// TODO: is there any warranty that the outter facet cycle enty point is always at first
	// in the cycles list?
	++fc; while( fc != fce)  { CGAL_assertion( assign( sv, fc)); ++fc; }
#endif
	TRACEN( "no adjacent facets were found (but incident edge(s)).");
	f_visible = Halffacet_const_handle();
      }
      else
	CGAL_assertion_msg(0,"Damn wrong handle");
    }
    return f_visible;
  }

  Halffacet_const_handle get_visible_facet( const Halfedge_const_handle e,
                                      const Ray_3& ray) const {
   //{\Mop when one shoot a ray |ray| in order to find the facet below to
   //  an object, and an edge |e| is hit, we need to choose one of the two 
   //  facets in the adjacency list of |e| that could be 'seen'  from the
   //  piercing point of the |ray| on the local (virtual) view  of |e|
   //  \precondition |ray| target belongs to |e|. } 

    SM_const_decorator SD(&*source(e));
    if( SD.is_isolated(e))
      return Halffacet_const_handle();
    
    // We search for the plane in the adjacency list of e, which is closest
    // to the ray. The cross product of the direction of e and the orthogonal
    // vector of a plane gives us a vector vec0/vec1 on the plane of the facet 
    // and orthogonal to e, pointing inside of the facet.

    Vector_3 ev(segment(e).to_vector()), rv(ray.to_vector());
    SHalfedge_around_svertex_const_circulator sh(SD.first_out_edge(e));
    Halffacet_const_handle res = facet(sh); 
    Vector_3 vec0(cross_product(ev,plane(res).orthogonal_vector()));
    /* // probably incorrect assertion
    CGAL_assertion_code
      (Sphere_segment _ess( SD.point(SD.source(sh)), 
			    SD.point(SD.source(next(sh))),
			    SD.circle(sh)));
    CGAL_assertion( _ess.has_on(vec0));
    */
    SHalfedge_around_svertex_const_circulator send(sh);
    TRACEN("initial face candidate "<< plane(res)<<" with vector  "<<vec0);

    // We compare the vectors vec0/vec1 of the facets. The one that is nearest 
    // to pointing in the opposite direction of the ray, is chosen. The 
    // respective facet is the nearest to the ray.

    sh++;
    CGAL_For_all(sh,send) {
      Vector_3 vec1(cross_product(ev,plane(facet(sh)).orthogonal_vector()));
      TRACEN("test face candidate "<< plane(facet(sh))<<" with vector  "<<vec1);
      FT sk0(rv*vec0),  sk1(rv*vec1);
      if(sk0<FT(0) && sk1>FT(0))
        continue;
      if(sk0>FT(0) && sk1<FT(0)) {
        res = facet(sh); 
	vec0 = vec1;
	continue;
      }

      // We have to comapare the two skalar products sk0 and sk1. Therefore 
      // we have to normalize the input vectors vec0 and vec1, which means
      // that we have to divide them by their lengths len0 and len1. 
      // To cicumvent irrational numbers, we sqaure the whole inequality.

      FT len0 = vec0.x()*vec0.x()+vec0.y()*vec0.y()+vec0.z()*vec0.z();
      FT len1 = vec1.x()*vec1.x()+vec1.y()*vec1.y()+vec1.z()*vec1.z();
      FT diff = len0*sk1*sk1 - len1*sk0*sk0;

      if((sk0>FT(0) && diff<FT(0)) || (sk0<FT(0) && diff>FT(0))) {
        res = facet(sh);
	vec0 = vec1;
      }
    }
 
    // We have to check which of the two halffacet is visible from
    // the ray. 

    if(rv*plane(res).orthogonal_vector() > FT(0))
      res = twin(res);

    TRACEN("return "<<plane(res));
    return res; // never reached
  }

  Halffacet_const_handle get_visible_facet( const Halffacet_const_handle f,
				      const Ray_3& ray) const 
    /*{\Mop when one shoot a ray |ray| in order to find the facet below to
      an object, and a facet |f| is hit, we need to choose the right facet
      from the halffacet pair |f| that  could be 'seen'  from the
      piercing point of the |ray| on the local (virtual) view  of |f|.
      \precondition |ray| target belongs to |f| and the intersection between
      |ray| and is not coplanar with |f|. }*/ {

    Halffacet_const_handle f_visible = f;
    // CGAL_assertion( !plane(f_visible).has_on(ray.source()));
    if( plane(f_visible).has_on_negative_side(ray.source()))
      f_visible = twin(f);
    CGAL_assertion( plane(f_visible).has_on_positive_side(ray.source()));
    return f_visible;
  }

  Halffacet_const_handle get_visible_facet( const Vertex_const_handle v, 
				      const Segment_3& ray) const 
    /*{\Mop when one shoots a ray |ray| in order to find the facet below to
      an object, and vertex |v| is hit, we need to choose one of the facets
      in the adjacency list of |v| such that it could be 'seen' from the
      piercing point of the |ray| on the sphere map on |v|.  We make it just
      locating the sphere facet |sf| pierced by |ray| and taking the adjacent 
      facet to one of the sphere segments on the boundary of |sf|.
      \precondition |ray| target is on |v| and the intersection between
      |ray| and the 2-skeleton incident to v is empty. }*/ {

    Halffacet_const_handle f_visible;
    CGAL_assertion( ray.source() != point(v));
    CGAL_assertion( ray.has_on(point(v)));
    Sphere_point sp(ray.source() - point(v));
    TRACEN( "Locating "<<sp <<" in "<<point(v));
    CGAL_assertion(Infi_box::degree(sp.hx()) < 2 && 
		   Infi_box::degree(sp.hy()) < 2 && 
		   Infi_box::degree(sp.hz()) < 2 && 
		   Infi_box::degree(sp.hw()) == 0);
    sp = Infi_box::simplify(sp);
    TRACEN( "Locating "<<sp <<" in "<<point(v));
    SM_point_locator L(v);
    Object_handle o = L.locate(sp);

    SFace_const_handle sf;
    CGAL_assertion(assign(sf,o));
    assign( sf, o);

    SFace_cycle_const_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) {
	TRACEN( "no adjacent facets were found.");
	f_visible =  Halffacet_const_handle();
    }
    else {
      SHalfedge_const_handle se; 
      SHalfloop_const_handle sl;
      SVertex_const_handle sv;
      if ( assign( se, fc)) {
	TRACEN( "adjacent facet found (SEdges cycle).");
	TRACEN("se"<<PH(se));
	f_visible = facet(twin(se));
	TRACEN("f_visible"<<&f_visible);
      }
      else if ( assign( sl, fc)) {
	TRACEN( "adjacent facet found (SHalfloop cycle).");
	f_visible = facet(twin(sl));
      }
      else if( assign( sv, fc)) {
	TRACEN( "no adjacent facets were found (but incident edge(s)).");
	f_visible = Halffacet_const_handle();
      }
      else
	CGAL_assertion_msg(0,"Damn wrong handle");
    }
    return f_visible;
  }

  Halffacet_const_handle get_visible_facet( const Halfedge_const_handle e,
                                      const Segment_3& ray) const {
    //{\Mop when one shoot a ray |ray| in order to find the facet below to
    //  an object, and an edge |e| is hit, we need to choose one of the two 
   //   facets in the adjacency list of |e| that could be 'seen'  from the
   //   piercing point of the |ray| on the local (virtual) view  of |e|
   //   \precondition |ray| target belongs to |e|. } 

    CGAL_assertion(false);

    SM_const_decorator SD;
    if( SD.is_isolated(e))
      return Halffacet_const_handle();

    Halffacet_const_handle res = facet(sh);    

    Vector_3 ed(segment(e).to_vector());
    Vector_3 ev(segment(e).to_vector()), rv(ray.to_vector());
    SHalfedge_around_svertex_const_circulator sh(SD.first_out_edge(e)), send(sh);
    Vector_3 vec0(cross_product(ev,plane(res).orthogonal_vector()));
    TRACEN("initial face candidate "<< plane(res));
   
    sh++;
    CGAL_For_all(sh,send) {
    Vector_3 vec1(cross_product(ev,plane(sh).orthogonal_vector()));
      RT sk0(rv*vec0),  sk1(rv*vec1);
      if(sk0<0 && sk1>0)
        continue;
      if(sk0>0 && sk1<0) {
        res = facet(sh);
	continue;
      }

      RT len0 = vec0.x()*vec0.x()+vec0.y()*vec0.y()+vec0.z()*vec0.z();
      RT len1 = vec1.x()*vec1.x()+vec1.y()*vec1.y()+vec1.z()*vec1.z();
      RT sq0 = sk0 * sk0;
      RT sq1 = sk1 * sk1;     
      RT diff = len0*sq1 - len1*sq0;

      if((sk0 > 0 && diff<0) || (sk0 < 0 && diff>0))
        res = facet(sh);
    }

    return Halffacet_const_handle(); // never reached
  }

  Halffacet_const_handle get_visible_facet( const Halffacet_const_handle f,
				      const Segment_3& ray) const 
    //{\Mop when one shoot a ray |ray| in order to find the facet below to
    //  an object, and a facet |f| is hit, we need to choose the right facet
    //  from the halffacet pair |f| that  could be 'seen'  from the
    //  piercing point of the |ray| on the local (virtual) view  of |f|.
    //  \precondition |ray| target belongs to |f| and the intersection between
    //  |ray| and is not coplanar with |f|. } 
    {
      Halffacet_const_handle f_visible = f;
      CGAL_assertion( !plane(f_visible).has_on(ray.source()));
      if( plane(f_visible).has_on_negative_side(ray.source()))
	f_visible = twin(f);
      CGAL_assertion( plane(f_visible).has_on_positive_side(ray.source()));
      return f_visible;
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
	    V.visit(SHalfedge_const_handle(ec));
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
	  V.visit(l);
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

