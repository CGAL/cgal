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
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_io_parser.h>
#include <CGAL/Nef_3/SNC_ray_shoter.h>
#ifdef  SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#undef _DEBUG
#define _DEBUG 19
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename SNC_structure_>
class SNC_decorator { 
  typedef SNC_structure_ SNC_structure;
  typedef SNC_decorator<SNC_structure>        Self;
  typedef SNC_constructor<SNC_structure>      SNC_constructor;
  typedef SNC_ray_shoter<SNC_structure>       SNC_ray_shoter;
  typedef SNC_SM_decorator<SNC_structure>     SM_decorator;
  typedef SNC_SM_overlayer<SNC_structure>     SM_overlayer;
  typedef SNC_SM_point_locator<SNC_structure> SM_point_locator;
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
  USING(SHalfedge_around_facet_const_circulator);
  USING(SHalfedge_around_facet_circulator);
  USING(SFace_cycle_iterator);
  USING(SFace_cycle_const_iterator);
  USING(Halffacet_cycle_iterator);
  USING(Halffacet_cycle_const_iterator);
  USING(Shell_entry_iterator);
  USING(Shell_entry_const_iterator);
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
#undef USING

#define DECUSING(t) typedef typename SM_decorator::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
  DECUSING(SHalfedge_around_svertex_circulator);
#undef DECUSING

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
  { std::stringstream os; set_pretty_mode(os);
    os << "sedge-use " << point(source(e)) << point(target(e))<<'\0';
    return os.str();
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
    CGAL_nef3_assertion( fc != f->facet_cycles_end());
    SHalfedge_handle se;
    if ( assign(se, fc) ) { 
      CGAL_nef3_assertion( facet(se) == f);
      CGAL_nef3_assertion( sface(se) != SFace_handle());
      CGAL_nef3_assertion( volume(sface(twin(se))) == volume(f));
      return sface(twin(se));
    } 
    else 
      CGAL_nef3_assertion_msg( 0, "Facet outer cycle entry point"
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
  void store_as_first_boundary_object(H h, Halffacet_handle f) const
  { f->boundary_entry_objects_.push_front(Object_handle(h));
    sncp()->store_boundary_item(h, --(f->facet_cycles_end()));
  }

  template <typename H>
  void undo_boundary_object(H h, Halffacet_handle f) const
  { CGAL_nef3_assertion(sncp()->is_boundary_object(h));
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

  template <typename H>
  void undo_boundary_object(H h, Volume_handle c) const
  { CGAL_nef3_assertion(sncp()->is_boundary_object(h));
    Shell_entry_iterator it = sncp()->boundary_item(h);
    sncp()->undef_boundary_item(h);
    c->shell_entry_objects_.erase(it);
  }

  template <typename H>
    void store_boundary_object(H h, Volume_handle c) const
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
    CGAL_nef3_assertion(c->shell_entry_objects_.size() == 0);
    Shell_volume_setter Setter(*this, c);
    visit_shell_objects( f, Setter );
    TRACEN("Volume "<<&*c<<", outer shell "<<&*f);
    store_boundary_object( f, c );
  }

  void link_as_inner_shell( SFace_handle f, Volume_handle c ) const {
    // CGAL_nef3_assertion(c->shell_entry_objects_.size() > 0);
    Shell_volume_setter Setter(*this, c);
    visit_shell_objects( f, Setter );
    TRACEN("Volume "<<&*c<<", inner shell "<<&*f);
    store_boundary_object( f, c );
  }

  struct Shell_mark_setter {
    const SNC_decorator D;
    Mark m;
    Shell_mark_setter(const SNC_decorator& Di, Mark mi)
      : D(Di), m(mi) {}
    void visit(SFace_handle h)        { /* empty */ }
    void visit(Vertex_handle h)       { D.mark(h) = m; }
    void visit(Halfedge_handle h)     { D.mark(h) = m; }
    void visit(Halffacet_handle h)    { D.mark(h) = m; }
    void set_volume(Volume_handle ci) { /* empty */ }
  };

  void clear_outer_box_marks() {
    SObject_handle o = shells_begin(volumes_begin());
    SFace_handle sf;
    CGAL_nef3_assertion( assign( sf, o));
    CGAL_assertion( sf != sncp()->sfaces_end());
    assign( sf, o);
    Shell_mark_setter Setter( *this, false);
    visit_shell_objects( sf, Setter );
  }

  template <class H> void set_facet(H h, Halffacet_handle f) const 
    { h->incident_facet_ = f; }
  void set_volume(Halffacet_handle h, Volume_handle c) const
    { h->volume_ = c; }
  void set_volume(SFace_handle h, Volume_handle c) const 
    { h->incident_volume_ = c; }

  void add_sloop_to_facet(SHalfloop_handle l, Halffacet_handle f) const {
    SM_decorator SD(vertex(l));
    Sphere_circle facet_plane(plane(f));
    if( facet_plane == SD.circle(l)) {
      l->incident_facet_ = f;
      SD.twin(l)->incident_facet_ = twin(f);
    } else {
      CGAL_nef3_assertion( facet_plane.opposite() == SD.circle(l));
      l->incident_facet_ = twin(f);
      SD.twin(l)->incident_facet_ = f;
    }
  }

  /* returns true if |f| is part of the infinimaximal box.*/
  bool is_infbox_facet(Halffacet_handle f) const {
    return (volume(f) == sncp()->volumes_begin() ||
	    volume(twin(f)) == sncp()->volumes_begin());
  }

  /* returns true when |v| has outdegree two.*/
  bool has_outdeg_two(SVertex_handle v) const {
    SM_decorator SD;
    if( SD.is_isolated(v))
      return false;
    SHalfedge_handle e1 = SD.first_out_edge(v);
    SHalfedge_handle e2 = SD.cyclic_adj_succ(e1);
    return( e1!=e2 && SD.cyclic_adj_succ(e2)==e1);
  }

  Halffacet_handle get_visible_facet( const Vertex_handle v, 
				      const Segment_3& ray) const 
    /*{\Mop when one shot a ray |ray| in order to find the facet below to
      an object, and vertex |v| is hit, we need to choose one of the facets
      in the adjacency list of |v| such that it could be 'seen' from the
      piercing point of the |ray| on the sphere map on |v|.  We make it just
      locating the sphere facet |sf| pierced by |ray| and taking the adjacent 
      facet to one of the sphere segments on the boundary of |sf|.
      \precondition |ray| target is on |v| and the intersection between
      |ray| and the 2-skeleton incident to v is empty. }*/ {
    Halffacet_handle f_visible;
    CGAL_nef3_assertion( ray.target() == point(v));
    Sphere_point sp(ray.source() - point(v));
    TRACEN( "Locating " << sp <<" in " << point(v));
    SM_point_locator L(v);
    SObject_handle o = L.locate(sp);
    SFace_const_handle sf;
    CGAL_nef3_assertion( assign( sf, o));
    assign( sf, o);
    SFace_cycle_const_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) {
	TRACEN( "no adjacent facet found.");
	f_visible =  Halffacet_handle();
    }
    else {
      SHalfedge_handle se; 
      SHalfloop_handle sl;
      if ( assign( se, fc)) {
	TRACEN( "adjacent facet found (SEdges cycle).");
	f_visible = facet(twin(se));
      }
      else if ( assign( sl, fc)) {
	TRACEN( "adjacent facet found (SHalfloop cycle).");
	f_visible = facet(twin(sl));
      }
      else 
	CGAL_nef3_assertion_msg(0, "Damn, wrong handle.");
    }
    return f_visible;
  }

  Halffacet_handle get_visible_facet( const Halfedge_handle e,
				      const Segment_3& ray) const
    /*{\Mop when one shot a ray |ray| in order to find the facet below to
      an object, and an edge |e| is hit, we need to choose one of the two 
      facets in the adjacency list of |e| that could be 'seen'  from the
      piercing point of the |ray| on the local (virtual) view  of |e|
      \precondition |ray| target belongs to |e|. }*/ {
    CGAL_nef3_assertion( segment(e).has_on( ray.target()));
    SM_decorator SD;
    if( SD.is_isolated(e))
      return Halffacet_handle();
    Direction_3 ed(segment(e).direction());
    Vector_3 ev(ed), rv(ed);
    SHalfedge_around_svertex_circulator sh(SD.first_out_edge(e)), sg(sh);
    Vector_3 h(plane(facet(twin(sh))).orthogonal_vector());
    TRACEN("initial face candidate "<<&*facet(twin(sh)));
    sg++;
    while ( true ) {
      Vector_3 g(plane(facet(twin(sg))).orthogonal_vector());
      if( CGAL_NTS is_positive( cross_product(g, ev) * h)) {
	if( CGAL_NTS is_negative( rv * g))
	  return facet(twin(sh));
	else {
	  sh = sg;
	  h = g;
	  TRACEN("new candidate "<<&*facet(twin(sh)));
	}
      }
      else
	return facet(twin(sh));
    }
    return Halffacet_handle(); // never reached
  }
  
  Halffacet_handle get_visible_facet( const Halffacet_handle f,
				      const Segment_3& ray) const 
    /*{\Mop when one shot a ray |ray| in order to find the facet below to
      an object, and a facet |f| is hit, we need to choose the right facet
      from the halffacet pair |f| that  could be 'seen'  from the
      piercing point of the |ray| on the local (virtual) view  of |f|.
      \precondition |ray| target belongs to |f| and the intersection between
      |ray| and is not coplanar with |f|. }*/ {
    Halffacet_handle f_visible = f;
    CGAL_nef3_assertion( !plane(f_visible).has_on(ray.source()));
    if( plane(f_visible).has_on_negative_side(ray.source()))
      f_visible = twin(f);
    CGAL_nef3_assertion( plane(f_visible).has_on_positive_side(ray.source()));
    return f_visible;
  }

  template <typename Selection>
    Vertex_handle binop_local_views( Vertex_handle v0, Vertex_handle v1,
				     const Selection& BOP, SNC_structure& rsnc)
    /*{\opOverlays two spheres maps.}*/ {
    typedef SNC_SM_io_parser<SNC_structure> SNC_SM_io_parser;
    SNC_SM_io_parser IO0( std::cerr, v0);
    SNC_SM_io_parser IO1( std::cerr, v1);
    TRACEN(" sphere maps before local binary operation");
    TRACEN(v0->debug());
    TRACEN(v1->debug());
    IO0.print();
    IO1.print();
    CGAL_assertion( point(v0) == point(v1));
    Vertex_handle v01 = rsnc.new_vertex( point(v0), BOP( mark(v0),mark(v1)));
    TRACEN("  binop result on vertex "<<&*v01<<" on "<<&*(v01->sncp()));
    SM_overlayer O(v01);
    O.subdivide( v0, v1);
    O.select( BOP);
    O.simplify();
    O.check_integrity_and_topological_planarity();

    TRACEN(" result sphere map:");
    SNC_SM_io_parser IO01( std::cerr, v01);
    TRACEN(v01->debug());
    IO01.print();
    TRACEN(" sphere maps after local binary operation");
    IO0.print();
    IO1.print();
#ifdef SM_VISUALIZOR
    typedef SNC_SM_visualizor<SNC_structure> SMV;
    CGAL::OGL::add_sphere();
    SMV V0(v0, CGAL::OGL::spheres_.back());
    V0.draw_map();
    SMV V1(v1, CGAL::OGL::spheres_.back());
    V1.draw_map();
    SMV V01(v01, CGAL::OGL::spheres_.back());
    V01.draw_map();
    CGAL::OGL::start_viewer();
    TRACEN("any key to continue...");
    char c;
    std::cin >> c;
#endif
    return v01;
  }

  Vertex_handle create_local_view_on( const Point_3& p, Halfedge_handle e) {
    SNC_constructor C(*sncp());
    return C.create_from_edge( e, p);
  }

  Vertex_handle create_local_view_on( const Point_3& p, Halffacet_handle f) {
    SNC_constructor C(*sncp());
    return C.create_from_facet( f, p);
  }

  Vertex_handle create_local_view_on( const Point_3& p, Volume_handle c) {
    Vertex_handle v = sncp()->new_vertex( p, mark(c));
    SM_decorator SD(v);
    SFace_handle f = SD.new_face();
    SD.mark(f) = mark(c);
    SM_point_locator PL(v);
    PL.init_marks_of_halfspheres(); // necessary to init default marks
    return v;
  }

  Vertex_handle qualify_with_respect( const Point_3 p,
				      SNC_structure& P1i,
				      SNC_structure& result)
    /*{\op }*/ {
    SNC_ray_shoter rs(P1i);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Volume_handle c;
    SNC_decorator D(result);
    Object_handle o = rs.locate(p);
    if( assign( v, o)) {
      TRACEN("<-> vertex local view on "<<point(v));
      return v;
    }
    else if( assign( e, o)) {
      TRACEN("<-> edge local view on "<<p);
      return D.create_local_view_on( p, e);
    }
    else if( assign( f, o)) {
      TRACEN("<-> facet local view on "<<p);
      return D.create_local_view_on( p, f);
    }
    else if( assign( c, o)) {
      TRACEN("<-> volume local view on "<<p);
      return D.create_local_view_on( p, c);
    }
    else CGAL_nef3_assertion_msg(0, "Where is the point then?");
    return Vertex_handle(); // never reached
  }

  template <typename Selection>
    void binary_operation( SNC_structure& snc1i, 
			   const Selection& BOP,
			   SNC_structure& result)
    /*{\opPerforms a binary operation defined on |BOP| between two
      SNC structures.  The input structures are not modified and the
      result of the operation is stored in |result|.
      \precondition: the structure |result| is empty.}*/ {
    typedef Unique_hash_map<Vertex_handle, bool> Hash_map;
    CGAL_nef3_assertion( result.is_empty());
    Hash_map Ignore(false);
    Vertex_iterator v0, v1;
    
    TRACEN("=> for all v0 in snc0, qualify v0 with respect snc1");

    TRACEN("vertices on snc0:");
    CGAL_nef3_forall_vertices( v0, *sncp()) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc1:");
    CGAL_nef3_forall_vertices( v0, snc1i) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc01:");
    CGAL_nef3_forall_vertices( v0, result) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    
    TRACEN("number of vertices on snc0 sn1 snc01: "<<
	   sncp()->number_of_vertices()<<' '<<
	   snc1i.number_of_vertices()<<' '<<
	   result.number_of_vertices());

    CGAL_nef3_forall_vertices( v0, *sncp()) {
      CGAL_nef3_assertion(!Ignore[v0]);
      v1 = qualify_with_respect( point(v0), snc1i, result);
      TRACEN("=> overlay of vertices v0 "<<&*v0<<" v1 "<<&*v1);
      binop_local_views( v0, v1, BOP, result);
      if( v1->sncp() == &result) /* if v1 is a copy */
	result.delete_vertex(v1);
      else
	Ignore[v1] = true;

      TRACEN("vertices on snc0 sn1 snc01: "<<
	     sncp()->number_of_vertices()<<' '<<
	     snc1i.number_of_vertices()<<' '<<
	     result.number_of_vertices());
    }

    TRACEN("=> for all v1 in snc1, qualify v1 with respect snc0");

    TRACEN("vertices on snc0:");
    CGAL_nef3_forall_vertices( v0, *sncp()) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc1:");
    CGAL_nef3_forall_vertices( v0, snc1i) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc01:");
    CGAL_nef3_forall_vertices( v0, result) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);

    CGAL_nef3_forall_vertices( v1, snc1i) {
      if( Ignore[v1]) continue;
      v0 = qualify_with_respect( point(v1), *sncp(), result);
      TRACEN("=> overlay of vertices v1 "<<&*v1<<" v0 "<<&*v0);
      binop_local_views( v0, v1, BOP, result);
      CGAL_nef3_assertion( v0->sncp() == &result);
      result.delete_vertex(v0);

      TRACEN("vertices on snc0 sn1 snc01: "<<
	     sncp()->number_of_vertices()<<' '<<
	     snc1i.number_of_vertices()<<' '<<
	     result.number_of_vertices());
    }

    TRACEN("=> edge facet intersection");

    TRACEN("vertices on snc0:");
    CGAL_nef3_forall_vertices( v0, *sncp()) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc1:");
    CGAL_nef3_forall_vertices( v0, snc1i) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);
    TRACEN("vertices on snc01:");
    CGAL_nef3_forall_vertices( v0, result) TRACEN(point(v0)<<&*(v0->sncp_));
    TRACEN("end vertices"<<std::endl);

    SNC_ray_shoter rs(*sncp());

    Halfedge_iterator e0, e1;
    Halffacet_iterator f0, f1;
    CGAL_nef3_forall_edges( e0, *sncp()) { 
      CGAL_nef3_forall_facets( f1, snc1i) { 
	Point_3 ip;
	if( rs.does_intersect_internally( segment(e0), f1, ip )) {
	  TRACEN(" edge0 face1 intersection...");
	  v0 = qualify_with_respect( ip, *sncp(), result);
	  v1 = qualify_with_respect( ip, snc1i, result);
	  binop_local_views( v0, v1, BOP, result);
	  result.delete_vertex(v0);
	  result.delete_vertex(v1);
	}
      }
    }
    CGAL_nef3_forall_edges( e1, snc1i) { 
      CGAL_nef3_forall_facets( f0, *sncp()) { 
	Point_3 ip;
	if( rs.does_intersect_internally( segment(e1), f0, ip )) {
	  TRACEN(" edge1 face0 intersection...");
	  v1 = qualify_with_respect( ip, snc1i, result);
	  v0 = qualify_with_respect( ip, *sncp(), result);
	  binop_local_views( v0, v1, BOP, result);
	  result.delete_vertex(v0);
	  result.delete_vertex(v1);
	}
      }
    }
    TRACEN("=> edge edge intersection");
    CGAL_nef3_forall_edges( e0, *sncp()) { 
      CGAL_nef3_forall_edges( e1, snc1i) { 
	Point_3 ip;
	if( rs.does_intersect_internally( segment(e0), segment(e1), ip )) {
	  TRACEN(" edge0 edge1 intersection...");
	  Vertex_handle v0, v1;
	  v0 = qualify_with_respect( ip, *sncp(), result);
	  v1 = qualify_with_respect( ip, snc1i, result);
	  binop_local_views( v0, v1, BOP, result);
	  result.delete_vertex(v0);
	  result.delete_vertex(v1);
	}
      }
    }
    TRACEN("=> resultant vertices: ");
    CGAL_nef3_forall_vertices( v0, result) {
      TRACEN(&*v0<<" "<<point(v0)<<&*(v0->sncp_));
    }
    TRACEN("=> pre-construction result");
    SNC_io_parser<SNC_structure> O(std::cout, result);
    O.print();

    // remove vertices whose local view is not that of a vertex
    Vertex_iterator vi, vin;
    for( vi = result.vertices_begin(); vi != result.vertices_end(); vi = vin) {
      vin = vi;
      vin++;
      SM_decorator SD(vi);
      if( (result.is_part_of_volume(vi) && 
	   mark(vi) == SD.mark(SD.sfaces_begin())) 
	  ||
	  (result.is_part_of_facet(vi) && 
	   mark(vi) == SD.mark(SD.shalfloop())) 
	  ||
	  (result.is_part_of_edge(vi) &&
	   mark(vi) == SD.mark(SD.svertices_begin())))
	result.delete_vertex(vi);
    }

    // synthesis of spatial structure
    SNC_constructor C(result);
    C.pair_up_halfedges();
    C.link_shalfedges_to_facet_cycles();
    C.categorize_facet_cycles_and_create_facets();
    C.create_volumes();
    result.simplify();

    TRACEN("=> construction completed, result: ");
    SNC_io_parser<SNC_structure> Op(std::cout, result);
    Op.print();

    TRACEN("=> end binary operation. ");
  }

  template <typename Visitor>
  void visit_shell_objects(SFace_handle f, Visitor& V) const;

  Vertex_iterator   vertices_begin()   { return sncp()->vertices_begin(); }
  Vertex_iterator   vertices_end()     { return sncp()->vertices_end(); }
  Halfedge_iterator halfedges_begin()  { return sncp()->halfedges_begin(); }
  Halfedge_iterator halfedges_end()    { return sncp()->halfedges_end(); }
  Halffacet_iterator halffacets_begin(){ return sncp()->halffacets_begin(); }
  Halffacet_iterator halffacets_end()  { return sncp()->halffacets_end(); }
  Volume_iterator   volumes_begin()    { return sncp()->volumes_begin(); }
  Volume_iterator   volumes_end()      { return sncp()->volumes_end(); }

  Shell_entry_iterator shells_begin(Volume_handle c) {
    return c->shells_begin();
  }
  Shell_entry_iterator shells_end(Volume_handle c) {
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
  typedef typename SM_decorator::SHalfedge_around_sface_circulator 
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
        } else CGAL_nef3_assertion_msg(0,"Damn wrong handle.");
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
        } else CGAL_nef3_assertion_msg(0,"Damn wrong handle.");
      }
    }
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_DECORATOR_H

