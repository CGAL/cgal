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
// WARRANTY OF DESISGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_SNC_DECORATOR_H
#define CGAL_SNC_DECORATOR_H

//#define CGAL_NEF3_SM_VISUALIZOR
//#define CGAL_NEF3_DUMP_SPHERE_MAPS
//#define CGAL_NEF3_DUMP_SNC_OPERATORS

#include <CGAL/basic.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_3/SNC_simplify.h>
#include <CGAL/Nef_3/binop_intersection_tests.h>
#include <CGAL/Nef_3/SNC_decorator_traits.h>

#include <CGAL/IO/Verbose_ostream.h>
#ifdef  CGAL_NEF3_SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // CGAL_NEF3_SM_VISUALIZOR
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 19
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

#if defined (CGAL_NEF3_TIMER_OVERLAY) || defined(CGAL_NEF3_TIMER_INTERSECTION)
CGAL::Timer timer_overlay;
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
int number_of_point_location_queries;
int number_of_ray_shooting_queries;
CGAL::Timer timer_point_location;
CGAL::Timer timer_ray_shooting;
#endif 

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
int number_of_edge_facet_overlays=0;
int number_of_clones=0;
int number_of_sphere_sweeps=0;
CGAL::Timer timer_sphere_sweeps;
#endif

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
int number_of_plane_sweeps=0;
CGAL::Timer timer_plane_sweeps;
#endif

#ifdef CGAL_NEF3_DUMP_STATISTICS
int number_of_intersections;
int number_of_intersection_candidates;
#endif

template <typename Map>
class SNC_decorator : public SNC_const_decorator<Map> { 
 public:
  typedef Map SNC_structure;
  typedef typename Map::Sphere_map                     Sphere_map;
  typedef CGAL::SNC_decorator<SNC_structure>           Self;
  typedef CGAL::SNC_const_decorator<SNC_structure>     Base;
  typedef Base                                         SNC_const_decorator;
  typedef CGAL::SNC_constructor<SNC_structure>         SNC_constructor;
  typedef CGAL::SM_decorator<Sphere_map>               SM_decorator;
  typedef CGAL::SM_const_decorator<Sphere_map>         SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<SM_decorator>         SM_overlayer;
  typedef CGAL::SM_point_locator<SM_decorator>         SM_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>        SNC_intersection;
  typedef CGAL::SNC_point_locator<SNC_decorator>       SNC_point_locator;
  typedef CGAL::SNC_simplify<SNC_structure>            SNC_simplify;
  SNC_structure* sncp_;

  typedef SNC_decorator_traits<SNC_structure>  Decorator_traits;

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

  typedef typename SNC_structure::SHalfedge_around_facet_circulator SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
 
  typedef typename SNC_structure::Vertex Vertex;
  typedef typename SNC_structure::Halfedge Halfedge;
  typedef typename SNC_structure::Halffacet Halffacet;
  typedef typename SNC_structure::Volume Volume;
  typedef typename SNC_structure::SVertex SVertex;
  typedef typename SNC_structure::SHalfedge SHalfedge;
  typedef typename SNC_structure::SHalfloop SHalfloop;
  typedef typename SNC_structure::SFace SFace;

  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::Object_iterator Object_iterator;

  typedef typename Base::Vertex_const_handle Vertex_const_handle;
  typedef typename Base::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Base::Halffacet_const_handle Halffacet_const_handle;
  typedef typename Base::Volume_const_handle Volume_const_handle;
  typedef typename Base::SVertex_const_handle SVertex_const_handle;
  typedef typename Base::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Base::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Base::SFace_const_handle SFace_const_handle;

  typedef typename Base::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Base::Halfedge_const_iterator Halfedge_const_iterator;
  typedef typename Base::Halffacet_const_iterator Halffacet_const_iterator;
  typedef typename Base::Volume_const_iterator Volume_const_iterator;
  typedef typename Base::SVertex_const_iterator SVertex_const_iterator;
  typedef typename Base::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Base::SHalfloop_const_iterator SHalfloop_const_iterator;
  typedef typename Base::SFace_const_iterator SFace_const_iterator;

  typedef typename Base::SHalfedge_around_facet_const_circulator SHalfedge_around_facet_const_circulator;
  typedef typename Base::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename Base::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
  typedef typename Base::Shell_entry_const_iterator Shell_entry_const_iterator;

  typedef typename Base::Kernel Kernel;
  typedef typename Base::FT FT;
  typedef typename Base::RT RT;

  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Segment_3 Segment_3;
  typedef typename Base::Ray_3 Ray_3;
  typedef typename Base::Line_3 Line_3;
  typedef typename Base::Plane_3 Plane_3;
  typedef typename Base::Vector_3 Vector_3;

  typedef typename Base::Sphere_kernel Sphere_kernel;
  typedef typename Base::Sphere_point Sphere_point;
  typedef typename Base::Sphere_segment Sphere_segment;
  typedef typename Base::Sphere_circle Sphere_circle;
  typedef typename Base::Sphere_direction Sphere_direction;

  typedef typename Base::Size_type Size_type;
  typedef typename Base::Mark Mark;
  typedef typename Base::Infi_box Infi_box;
  typedef typename Base::Aff_transformation_3 Aff_transformation_3;

  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
                                 SHalfedge_around_svertex_circulator;
  typedef typename SM_decorator::SHalfedge_around_sface_circulator 
                                 SHalfedge_around_sface_circulator;

  struct points_lt {
    bool operator()(Point_3& p1, Point_3& p2) const {
      return CGAL::lexicographically_xyz_smaller(p1,p2);
    }
  };

  enum {NO_SNC, WITH_SNC};

 public:
  typedef void* GenPtr;

  SNC_decorator() : Base(), sncp_() {}
  SNC_decorator(SNC_structure& W) 
    : Base(W), sncp_(&W) {}
  SNC_decorator(const Self& S) {
    sncp_ = S.sncp_;
  }

  SNC_structure* sncp() const { 
    CGAL_assertion( sncp_ != NULL);
    return sncp_; 
  }

  void set_snc(SNC_structure& W) {
    sncp_ = &W;
  }

  std::string debug(SHalfedge_handle e) const
  { std::stringstream os; set_pretty_mode(os);
    os << "sedge-use " << e->source()->source()->point() 
       << e->twin()->source()->twin()->source()->point() <<'\0';
    return os.str();
  }

  SFace_handle adjacent_sface(Halffacet_handle f) const {
    Halffacet_cycle_iterator fc(f->facet_cycles_begin());
    CGAL_assertion( fc != f->facet_cycles_end());
    if ( fc.is_shalfedge() ) {
      SHalfedge_handle se(fc);
      CGAL_assertion( se->facet() == f);
      CGAL_assertion( se->incident_sface() != SFace_handle());
      CGAL_assertion( se->twin()->incident_sface()->volume() == f->incident_volume());
      return se->twin()->incident_sface();
    } 
    else 
      CGAL_assertion_msg( 0, "Facet outer cycle entry point"
			     "is not an SHalfedge? ");
    return SFace_handle(); // never reached
  }

  static Vector_3 to_vector(Halfedge_handle e) {  // rename to to_vector
    return Vector_3(e->vector()-CGAL::ORIGIN);
  }

  static Segment_3 segment(Halfedge_handle e)
    { return Segment_3(e->source()->point(),
		       e->twin()->source()->point()); }

  GenPtr& info(Vertex_handle v) const
  { return v->info(); }

  template <typename H>
  bool is_boundary_object(H h) const
  { return sncp()->is_boundary_object(h); }

  template <typename H>
  void store_boundary_object(H h, Halffacet_handle f) const
  { f->boundary_entry_objects().push_back(Object_handle(h));
    sncp()->store_boundary_item(h, --(f->facet_cycles_end()));
  }

  template <typename H>
  void store_as_first_boundary_object(H h, Halffacet_handle f) const
  { f->boundary_entry_objects().push_front(Object_handle(h));
    sncp()->store_boundary_item(h, --(f->facet_cycles_end()));
  }

  template <typename H>
  void undo_boundary_object(H h, Halffacet_handle f) const
  { CGAL_assertion(sncp()->is_boundary_object(h));
    Halffacet_cycle_iterator it = sncp()->boundary_item(h);
    sncp()->undef_boundary_item(h);
    f->boundary_entry_objects().erase(it);
  }

  void link_as_facet_cycle(SHalfedge_handle e, Halffacet_handle f) const
  { SHalfedge_around_facet_circulator hfc(e), hend(hfc);
    CGAL_For_all(hfc,hend) hfc->facet() = f;
    store_boundary_object(e,f);
  } 

  void link_as_interior_loop(SHalfloop_handle l, Halffacet_handle f) const
  { l->facet() = f;
    store_boundary_object(l,f);
  } 

  template <typename H>
  void make_twins(H h1, H h2) const
  { h1->twin() = h2; h2->twin() = h1; }

  void link_as_prev_next_pair(SHalfedge_handle e1, SHalfedge_handle e2) const
  { e1->next() = e2; e2->prev() = e1; } 

  template <typename H>
  void undo_boundary_object(H h, Volume_handle c) const
  { CGAL_assertion(sncp()->is_boundary_object(h));
    Shell_entry_iterator it = sncp()->boundary_item(h);
    sncp()->undef_boundary_item(h);
    c->shell_entry_objects().erase(it);
  }

  template <typename H>
    void store_boundary_object(H h, Volume_handle c, bool at_front = false) const
  { 
    if(at_front) {
      c->shell_entry_objects().push_front(Object_handle(h));
      sncp()->store_boundary_item(h, c->shells_begin());
    }
    else {
      c->shell_entry_objects().push_back(Object_handle(h));
      sncp()->store_boundary_item(h, --(c->shells_end()));
    }
  }

  template<typename SNCD_>
  struct Shell_volume_setter {
    const SNCD_ D;
    Volume_handle c;
    typedef Unique_hash_map< SFace_handle, bool> SFace_map;
    SFace_map linked;
    Shell_volume_setter(const SNCD_& Di)
      : D(Di), linked(false) {}
    void visit(SFace_handle h) { 
      CGAL_NEF_TRACEN(h->center_vertex()->point()); 
      D.set_volume(h, c); 
      linked[h] = true;
    }
    void visit(Vertex_handle h) { /* empty */ }
    void visit(Halfedge_handle h) { /* empty */ }
    void visit(Halffacet_handle h ) { D.set_volume(h, c); }
    void visit(SHalfedge_handle se) {}
    void visit(SHalfloop_handle sl) {}
    void set_volume(Volume_handle ci) { c = ci; }
    bool is_linked(SFace_handle h) {return linked[h];}
  };

  void link_as_outer_shell( SFace_handle f, Volume_handle c ) const {
    //    CGAL_assertion(c->shell_entry_objects().size() == 0);
    Shell_volume_setter<SNC_decorator> Setter(*this);
    Setter.set_volume(c);
    visit_shell_objects( f, Setter );
    CGAL_NEF_TRACEN("Volume "<<&*c<<", outer shell "<<&*f);
    store_boundary_object( f, c );
  }

  void link_as_inner_shell( SFace_handle f, Volume_handle c ) const {
    // CGAL_assertion(c->shell_entry_objects().size() > 0);
    Shell_volume_setter<SNC_decorator> Setter(*this);
    Setter.set_volume(c);
    visit_shell_objects( f, Setter );
    CGAL_NEF_TRACEN("Volume "<<&*c<<", inner shell "<<&*f);
    store_boundary_object( f, c);
  }

  template <class H> void set_facet(H h, Halffacet_handle f) const 
    { h->facet() = f; }
  void set_volume(Halffacet_handle h, Volume_handle c) const
    { h->incident_volume() = c; }
  void set_volume(SFace_handle h, Volume_handle c) const 
    { h->volume() = c; }

  void add_sloop_to_facet(SHalfloop_handle l, Halffacet_handle f) const {
    SM_decorator SD(&*l->incident_sface()->center_vertex());
    Sphere_circle facet_plane(f->plane());
    if( facet_plane == l->circle()) {
      l->facet() = f;
      l->twin()->facet() = f->twin();
    } else {
      CGAL_assertion( facet_plane.opposite() == l->circle());
      l->facet() = f->twin();
      l->twin()->facet() = f;
    }
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
				      const Ray_3& ray) const
    /*{\Mop when one shoot a ray |ray| in order to find the facet below to
      an object, and vertex |v| is hit, we need to choose one of the facets
      in the adjacency list of |v| such that it could be 'seen' from the
      piercing point of the |ray| on the sphere map on |v|.  We make it just
      locating the sphere facet |sf| pierced by |ray| and taking the adjacent 
      facet to one of the sphere segments on the boundary of |sf|.
      \precondition |ray| target is on |v| and the intersection between
      |ray| and the 2-skeleton incident to v is empty. }*/ {

    Halffacet_handle f_visible;
    CGAL_assertion( ray.source() != v->point());
    CGAL_assertion( ray.has_on(v->point()));
    Sphere_point sp(ray.source() - v->point());
    CGAL_NEF_TRACEN( "Locating "<<sp <<" in "<< v->point());
    CGAL_assertion(Infi_box::degree(sp.hx()) < 2 && 
		   Infi_box::degree(sp.hy()) < 2 && 
		   Infi_box::degree(sp.hz()) < 2 && 
		   Infi_box::degree(sp.hw()) == 0);
    sp = Infi_box::simplify(sp);
    CGAL_NEF_TRACEN( "Locating "<<sp <<" in "<< v->point());
    SM_point_locator L(&*v);
    Object_handle o = L.locate(sp);

    SFace_handle sf;
    if(!CGAL::assign(sf,o)) {
      SHalfedge_handle se;
      if(CGAL::assign(se,o))
	CGAL_NEF_TRACEN("on sedge " << PH(se) << " on facet " << se->facet()->plane());
      SVertex_handle sv;
      if(CGAL::assign(sv,o))
	CGAL_NEF_TRACEN("on svertex " << sv->point());
      CGAL_assertion_msg( 0, "it is not possible to decide which one is a visible facet (if any)");
      return Halffacet_handle();
    }

    SFace_cycle_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) {
	CGAL_NEF_TRACEN( "no adjacent facet found.");
	f_visible =  Halffacet_handle();
    }
    else {
      if ( fc.is_shalfedge()) {
	SHalfedge_handle se(fc);
	CGAL_NEF_TRACEN( "adjacent facet found (SEdges cycle).");
	CGAL_NEF_TRACEN("se"<<PH(se));
	CGAL_NEF_TRACEN(se->facet()->plane() <<"/"<<
	       se->snext()->facet()->plane()  <<"/"<< 
	       se->snext()->snext()->facet()->plane());
	f_visible = se->twin()->facet();
	CGAL_NEF_TRACEN("f_visible"<<f_visible->plane());
      }
      else if ( fc.is_shalfloop()) {
	SHalfloop_handle sl(fc);
	CGAL_NEF_TRACEN( "adjacent facet found (SHalfloop cycle)."<< sl->circle() 
			 << " with facet "<<sl->facet()->plane());
	f_visible = sl->twin()->facet();
	CGAL_NEF_TRACEN("f_visible"<<f_visible->plane());
      }
      else if(fc.is_svertex()) {
#ifdef CGAL_NEF_DEBUG
	// TODO: is there any warranty that the outter facet cycle enty point is always at first
	// in the cycles list?
	++fc; while( fc != fce)  { CGAL_assertion( fc.is_svertex()); ++fc; }
#endif
	CGAL_NEF_TRACEN( "no adjacent facets were found (but incident edge(s)).");
	f_visible = Halffacet_handle();
      }
      else
	CGAL_assertion_msg(0,"Damn wrong handle");
    }
    return f_visible;
  }

  Halffacet_handle get_visible_facet( const Halfedge_handle e,
                                      const Ray_3& ray) const {
   //{\Mop when one shoot a ray |ray| in order to find the facet below to
   //  an object, and an edge |e| is hit, we need to choose one of the two 
   //  facets in the adjacency list of |e| that could be 'seen'  from the
   //  piercing point of the |ray| on the local (virtual) view  of |e|
   //  \precondition |ray| target belongs to |e|. } 

    SM_decorator SD(&*e->source());
    if( SD.is_isolated(e))
      return Halffacet_handle();
    
    // We search for the plane in the adjacency list of e, which is closest
    // to the ray. The cross product of the direction of e and the orthogonal
    // vector of a plane gives us a vector vec0/vec1 on the plane of the facet 
    // and orthogonal to e, pointing inside of the facet.

    Vector_3 ev(segment(e).to_vector()), rv(ray.to_vector());
    SHalfedge_around_svertex_circulator sh(SD.first_out_edge(e));
    Halffacet_handle res = sh->facet(); 
    Vector_3 vec0(cross_product(ev,res->plane().orthogonal_vector()));
    /* // probably incorrect assertion
    CGAL_assertion_code
      (Sphere_segment _ess( sh->source()->source()->point(),
			    sh->next()->source()->source()->point(),
			    sh->circle());
    CGAL_assertion( _ess.has_on(vec0));
    */
    SHalfedge_around_svertex_circulator send(sh);
    CGAL_NEF_TRACEN("initial face candidate "<< res->plane()<<" with vector  "<<vec0);

    // We compare the vectors vec0/vec1 of the facets. The one that is nearest 
    // to pointing in the opposite direction of the ray, is chosen. The 
    // respective facet is the nearest to the ray.

    sh++;
    CGAL_For_all(sh,send) {
      Vector_3 vec1(cross_product(ev,sh->facet()->plane().orthogonal_vector()));
      CGAL_NEF_TRACEN("test face candidate "<< sh->facet()->plane()<<" with vector  "<<vec1);
      FT sk0(rv*vec0),  sk1(rv*vec1);
      if(sk0<=FT(0) && sk1>=FT(0))
	continue;
      if(sk0>=FT(0) && sk1<=FT(0)) {
        res = sh->facet(); 
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

      // if sk0<0 (and therefore sk1<0) both vectors point in a good direction.
      // Therefore we take the one pointing more in the good direction.
      // if sk0>0 (and therefore sk1>0) both vectors point in a bad direction.
      // Therefore we take the one pointing less in the bad direction.

      if((sk0>FT(0) && diff<FT(0)) || (sk0<FT(0) && diff>FT(0))) {
        res = sh->facet();
	vec0 = vec1;
      }
    }
 
    // We have to check which of the two halffacet is visible from
    // the ray. 

    if(rv*res->plane().orthogonal_vector() > FT(0))
      res = res->twin();

    CGAL_NEF_TRACEN("return "<<res->plane());
    return res; // never reached
  }

  Halffacet_handle get_visible_facet( const Halffacet_handle f,
				      const Ray_3& ray) const 
    /*{\Mop when one shoot a ray |ray| in order to find the facet below to
      an object, and a facet |f| is hit, we need to choose the right facet
      from the halffacet pair |f| that  could be 'seen'  from the
      piercing point of the |ray| on the local (virtual) view  of |f|.
      \precondition |ray| target belongs to |f| and the intersection between
      |ray| and is not coplanar with |f|. }*/ {
    
    CGAL_NEF_TRACEN("get visible facet " << ray << ", " << f->plane() 
		    << " has on source " << f->plane().has_on(ray.source()));
    Halffacet_handle f_visible = f;
    // CGAL_assertion( !f_visible->plane().has_on(ray.source()));
    if( f_visible->plane().has_on_negative_side(ray.source()))
      f_visible = f->twin();
    CGAL_assertion( f_visible->plane().has_on_positive_side(ray.source()));
    return f_visible;
  }

  Halffacet_handle get_visible_facet( const Vertex_handle v, 
				      const Segment_3& ray) const 
    /*{\Mop when one shoots a ray |ray| in order to find the facet below to
      an object, and vertex |v| is hit, we need to choose one of the facets
      in the adjacency list of |v| such that it could be 'seen' from the
      piercing point of the |ray| on the sphere map on |v|.  We make it just
      locating the sphere facet |sf| pierced by |ray| and taking the adjacent 
      facet to one of the sphere segments on the boundary of |sf|.
      \precondition |ray| target is on |v| and the intersection between
      |ray| and the 2-skeleton incident to v is empty. }*/ {

    Halffacet_handle f_visible;
    CGAL_assertion( ray.source() != v->point());
    CGAL_assertion( ray.has_on(v->point()));
    Sphere_point sp(ray.source() - v->point());
    CGAL_NEF_TRACEN( "Locating "<<sp <<" in "<<v->point());
    CGAL_assertion(Infi_box::degree(sp.hx()) < 2 && 
		   Infi_box::degree(sp.hy()) < 2 && 
		   Infi_box::degree(sp.hz()) < 2 && 
		   Infi_box::degree(sp.hw()) == 0);
    sp = Infi_box::simplify(sp);
    CGAL_NEF_TRACEN( "Locating "<<sp <<" in "<< v->point());
    SM_point_locator L(v);
    Object_handle o = L.locate(sp);

    SFace_handle sf;
    CGAL_assertion(CGAL::assign(sf,o));
    CGAL::assign(sf,o);

    SFace_cycle_iterator fc = sf->sface_cycles_begin(),
      fce = sf->sface_cycles_end();
    if( is_empty_range( fc, fce)) {
	CGAL_NEF_TRACEN( "no adjacent facets were found.");
	f_visible =  Halffacet_handle();
    }
    else {
      if (fc.is_shalfege()) {
	SHalfedge_handle se(fc);
	CGAL_NEF_TRACEN( "adjacent facet found (SEdges cycle).");
	CGAL_NEF_TRACEN("se"<<PH(se));
	f_visible = se->twin()->facet();
	CGAL_NEF_TRACEN("f_visible"<<&f_visible);
      }
      else if (fc.is_shalfloop()) {
	SHalfloop_handle sl(fc);
	CGAL_NEF_TRACEN( "adjacent facet found (SHalfloop cycle).");
	f_visible = sl->twin()->facet();
      }
      else if(fc.is_svertex()) {
	CGAL_NEF_TRACEN( "no adjacent facets were found (but incident edge(s)).");
	f_visible = Halffacet_handle();
      }
      else
	CGAL_assertion_msg(0,"Damn wrong handle");
    }
    return f_visible;
  }

  /*
  Halffacet_handle get_visible_facet( const Halfedge_handle e,
                                      const Segment_3& ray) const {
    //{\Mop when one shoot a ray |ray| in order to find the facet below to
    //  an object, and an edge |e| is hit, we need to choose one of the two 
   //   facets in the adjacency list of |e| that could be 'seen'  from the
   //   piercing point of the |ray| on the local (virtual) view  of |e|
   //   \precondition |ray| target belongs to |e|. } 

    CGAL_assertion(false);

    SM_decorator SD;
    if( SD.is_isolated(e))
      return Halffacet_handle();

    Halffacet_handle res = sh->facet();    

    Vector_3 ed(segment(e).to_vector());
    Vector_3 ev(segment(e).to_vector()), rv(ray.to_vector());
    SHalfedge_around_svertex_circulator sh(SD.first_out_edge(e)), send(sh);
    Vector_3 vec0(cross_product(ev,res->plane().orthogonal_vector()));
    CGAL_NEF_TRACEN("initial face candidate "<< res->plane());
   
    sh++;
    CGAL_For_all(sh,send) {
    Vector_3 vec1(cross_product(ev,sh->plane().orthogonal_vector()));
      RT sk0(rv*vec0),  sk1(rv*vec1);
      if(sk0<0 && sk1>0)
        continue;
      if(sk0>0 && sk1<0) {
        res = sh->facet();
	continue;
      }

      RT len0 = vec0.x()*vec0.x()+vec0.y()*vec0.y()+vec0.z()*vec0.z();
      RT len1 = vec1.x()*vec1.x()+vec1.y()*vec1.y()+vec1.z()*vec1.z();
      RT sq0 = sk0 * sk0;
      RT sq1 = sk1 * sk1;     
      RT diff = len0*sq1 - len1*sq0;

      if((sk0 > 0 && diff<0) || (sk0 < 0 && diff>0))
        res = sh->facet();
    }

    return Halffacet_handle(); // never reached
  }
  */

  Halffacet_handle get_visible_facet( const Halffacet_handle f,
				      const Segment_3& ray) const 
    //{\Mop when one shoot a ray |ray| in order to find the facet below to
    //  an object, and a facet |f| is hit, we need to choose the right facet
    //  from the halffacet pair |f| that  could be 'seen'  from the
    //  piercing point of the |ray| on the local (virtual) view  of |f|.
    //  \precondition |ray| target belongs to |f| and the intersection between
    //  |ray| and is not coplanar with |f|. } 
    {
      Halffacet_handle f_visible = f;
      CGAL_assertion( !f_visible()->plane().has_on(ray.source()));
      if( f_visible()->plane().has_on_negative_side(ray.source()))
	f_visible = f->twin();
      CGAL_assertion( f_visible()->plane().has_on_positive_side(ray.source()));
      return f_visible;
    }

  template <typename Selection>
    Vertex_handle binop_local_views( Vertex_const_handle v0, Vertex_const_handle v1,
				     const Selection& BOP, SNC_structure& rsnc)
    /*{\opOverlays two spheres maps.}*/ {
  
#ifdef CGAL_NEF3_DUMP_SPHERE_MAPS
    typedef SM_io_parser<SM_decorator> SM_io_parser;
    SM_io_parser IO0( std::cerr, v0);
    SM_io_parser IO1( std::cerr, v1);
    CGAL_NEF_TRACEN(" sphere maps before local binary operation");
    CGAL_NEF_TRACEN(v0->debug());
    CGAL_NEF_TRACEN(v1->debug());
    IO0.debug();    IO1.debug();
    IO0.print();    IO1.print();
#endif // CGAL_NEF3_DUMP_SPHERE_MAPS
 
   CGAL_assertion( v0->point() == v1->point());
    Vertex_handle v01 = rsnc.new_vertex( v0->point(), BOP( v0->mark(),v1->mark()));
    //    cerr<<"BOP Vertex "<<v0->mark()<<" "<<v1->mark()<<std::endl;
    CGAL_NEF_TRACEN("  binop result on vertex "<<&*v01<<" on "<<&*(v01->sncp()));
    SM_overlayer O(&*v01);
    O.subdivide( &*v0, &*v1);
    O.select( BOP);
    O.simplify();

#ifdef CGAL_NEF3_DUMP_SPHERE_MAPS
    CGAL_NEF_TRACEN(" result sphere map:");
    SM_io_parser IO01( std::cerr, v01);
    CGAL_NEF_TRACEN(v01->debug());
    IO01.print();
    CGAL_NEF_TRACEN(" sphere maps after local binary operation");
    IO0.debug();
    IO1.debug();
#endif //CGAL_NEF3_DUMP_SPHERE_MAPS

#ifdef CGAL_NEF3_SM_VISUALIZOR
    typedef SNC_SM_visualizor<SNC_structure> SMV;
    CGAL::OGL::add_sphere();
    SMV V0(v0, CGAL::OGL::spheres_.back());
    V0.draw_map();
    SMV V1(v1, CGAL::OGL::spheres_.back());
    V1.draw_map();
    SMV V01(v01, CGAL::OGL::spheres_.back());
    V01.draw_map();
    CGAL::OGL::start_viewer();
    CGAL_NEF_TRACEN("any key to continue...");
    char c;
    std::cin >> c;
#endif // CGAL_NEF3_VISUALIZOR
    return v01;
  }

  Vertex_handle create_local_view_on( const Point_3& p, Halfedge_const_handle e) {
    SNC_constructor C(*sncp());
    return C.create_from_edge( e, p);
  }
  
  Vertex_handle create_local_view_on( const Point_3& p, Halffacet_const_handle f) {
    SNC_constructor C(*sncp());
    return C.create_from_facet( f, p);
  }

  Vertex_handle create_local_view_on( const Point_3& p, Volume_const_handle c) {
    Vertex_handle v = sncp()->new_vertex( p, c->mark());
    SM_decorator SD(&*v);
    SFace_handle f = SD.new_sface();
    f->mark() = c->mark();
    CGAL_NEF_TRACEN("volume "<<&*c<<" marked as "<<c->mark()); 
    //    SM_point_locator PL(v);
    //    PL.init_marks_of_halfspheres(); // necessary to init default marks
    return v;
  }

  Vertex_handle 
    qualify_with_respect( const Point_3 p, 
			  const SNC_point_locator* pl1,
			  SNC_structure& result) {
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Volume_handle c;
    Self D(result);

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    number_of_point_location_queries++;
    timer_point_location.start();
#endif
    Object_handle o = pl1->locate(p);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    timer_point_location.stop();
#endif
    if( CGAL::assign( v, o)) {
      CGAL_NEF_TRACEN("<-> vertex local view on "<< v->point());
      return v;
    }
    else if( CGAL::assign( e, o)) {
      CGAL_NEF_TRACEN("<-> edge local view of "<<p<<" on "<<&*e);
      return D.create_local_view_on( p, e);
    }
    else if( CGAL::assign( f, o)) {
      CGAL_NEF_TRACEN("<-> facet local view of "<<p<<" on "<<&*f);
      return D.create_local_view_on( p, f);
    }
    else if( CGAL::assign( c, o)) {
      CGAL_NEF_TRACEN("<-> volume local view of "<<p<<" on "<<&*c);
      return D.create_local_view_on( p, c);
    }
    else CGAL_assertion_msg(0, "Where is the point then?");
    return Vertex_handle(); // never reached
  }

  /*
  Vertex_handle 
    qualify_with_respect( const Point_3& p, 
	                  Object_handle o,
			  const SNC_structure& snc) {
    Halfedge_const_handle e;
    Halffacet_handle f;
    if( CGAL::assign( e, o))
      return create_local_view_on( p, e);
    else if( CGAL::assign( f, o))
      return create_local_view_on( p, f);
    else 
      CGAL_assertion_msg( 0, "wrong handle");
    return Vertex_handle(); // never reached
  }
      */

  template <typename SNC_decorator, typename Selection>
  class Intersection_call_back : 
    public SNC_point_locator::Intersection_call_back 
  {
    typedef typename SNC_decorator::Decorator_traits Decorator_traits;
    typedef typename Decorator_traits::Halfedge_handle Halfedge_handle;
    typedef typename Decorator_traits::Halffacet_handle Halffacet_handle;
    
  public:
    Intersection_call_back( SNC_structure& s0, SNC_structure& s1,
			    const Selection& _bop, SNC_structure& r, 
			    bool invert_order = false) :
      snc0(s0), snc1(s1), bop(_bop), result(r),
      inverse_order(invert_order) {}

      void operator()(Halfedge_handle e0, Object_handle o1, const Point_3& ip)
      const {

#ifdef CGAL_NEF3_DUMP_STATISTICS
      ++number_of_intersections;
#endif

      Halfedge_handle e;
      Halffacet_handle f;
      
      Point_3 p(normalized(ip));

      CGAL_NEF_TRACEN("Intersection_call_back: intersection reported on " << p << " (normalized: " << normalized(p) << " )");
#ifdef CGAL_NEF_DEBUG
      CGAL_NEF_TRACEN("edge 0 has source " << e0->source()->point() << " and direction " << e0->vector());
      if( CGAL::assign( e, o1)) {
	CGAL_NEF_TRACEN("edge 1 has source " << e->source()->point() << " and direction " << e->vector());
      }
      else if( CGAL::assign( f, o1)) {
	CGAL_NEF_TRACEN("face 1 has plane equation " << f->plane());
      }
      else 
      	CGAL_assertion_msg( 0, "wrong handle");
#endif      

#if defined (CGAL_NEF3_TIMER_OVERLAY) || (CGAL_NEF3_TIMER_INTERSECTION)
      timer_overlay.start();
#endif

      if( CGAL::assign( e, o1)) {
	Self D(result);
	Vertex_handle v0, v1;
	v0 = D.create_local_view_on( p, e0);
	v1 = D.create_local_view_on( p, e);
	if( inverse_order)
	  std::swap( v0, v1);
	D.binop_local_views( v0, v1, bop, result);
	result.delete_vertex(v0);
	result.delete_vertex(v1);
      }
      else if( CGAL::assign( f, o1)) {
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
	Self D(result);
	Vertex_handle v0, v1;
	v0 = D.create_local_view_on( p, e0);
	v1 = D.create_local_view_on( p, f);
	if( inverse_order)
	  std::swap( v0, v1);
	D.binop_local_views( v0, v1, bop, result);
	result.delete_vertex(v0);
	result.delete_vertex(v1);	
#else
	SNC_constructor C(result);
	Sphere_map* M0 = C.create_edge_facet_overlay(e0, f, p, bop, inverse_order);
	SM_overlayer O(M0);
	O.simplify();
#endif
      }
      else 
	CGAL_assertion_msg( 0, "wrong handle");

#if defined (CGAL_NEF3_TIMER_OVERLAY) || (CGAL_NEF3_TIMER_INTERSECTION)
      timer_overlay.stop();
#endif

    }
  private:
    const SNC_structure& snc0;
    const SNC_structure& snc1;
    const Selection& bop;
    SNC_structure& result;
    bool inverse_order;
  };

  template <typename Selection>
    void binary_operation( SNC_point_locator* pl0,
			   const SNC_structure& snc1,
			   const SNC_point_locator* pl1,
			   const SNC_structure& snc2,
			   const SNC_point_locator* pl2,
			   const Selection& BOP)
    /*{\opPerforms a binary operation defined on |BOP| between two
      SNC structures.  The input structures are not modified and the
      result of the operation is stored in |result|.
      \precondition: the structure |result| is empty.}*/ {
    CGAL_assertion( sncp()->is_empty());
    CGAL_assertion( pl1 != NULL && pl2 != NULL);

    //    Progress_indicator_clog v_qualifying
    //      (sncp()->number_of_vertices()+snc1i.number_of_vertices(),
    //       "binary_operator: qualifying vertices...");

#ifdef CGAL_NEF3_TIMER_BINARY_OPERATION
    CGAL::Timer timer_binary_operation;
    timer_binary_operation.start();
#endif 

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    timer_sphere_sweeps.reset();
    number_of_sphere_sweeps=0;
    number_of_edge_facet_overlays=0;
    number_of_clones=0;
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION    
    timer_point_location.reset();
    number_of_point_location_queries=0;
#endif

#ifdef CGAL_NEF3_TIMER_OVERLAY
      timer_overlay.reset();
#endif

#ifdef CGAL_NEF3_DUMP_STATISTICS
      number_of_intersections=0;
      number_of_intersection_candidates=0;
#endif

    Unique_hash_map<Vertex_const_handle, bool> ignore(false);
    Vertex_const_iterator v0;

    //    CGAL_NEF_SETDTHREAD(19*43*131*509);
    CGAL_NEF_TRACEN("=> binary operation");
#ifdef CGAL_NEF3_DUMP_SNC_OPERATORS
    CGAL_NEF_TRACEN("=> first operand:");
    SNC_io_parser<SNC_structure> O0(std::cout, snc1);
    O0.print();
    CGAL_NEF_TRACEN("=> second operand:");
    SNC_io_parser<SNC_structure> O1(std::cout, snc2);
    O1.print();
#endif // CGAL_NEF3_DUMP_SNC_OPERATORS

    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());

    CGAL_NEF_TRACEN("=> for all v0 in snc1, qualify v0 with respect snc2");
    CGAL_forall_vertices( v0, snc1) {
      //      v_qualifying++;
      CGAL_assertion(!ignore[v0]);
      Point_3 p0(v0->point());
      Vertex_handle v;
      Halfedge_handle e;
      Halffacet_handle f;
      Volume_handle c;
      CGAL_NEF_TRACEN("Locating point " << p0);

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      ++number_of_point_location_queries;
      timer_point_location.start();
#endif 
      Object_handle o = pl2->locate(p0);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      timer_point_location.stop();
#endif 
    
#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.start();
#endif
      if( CGAL::assign( v, o)) {
	CGAL_NEF_TRACEN("p0 found on vertex");
	binop_local_views( v0, v, BOP, *sncp());
	ignore[v] = true;
      }
      else if( CGAL::assign( e, o)) {
	CGAL_NEF_TRACEN("p0 found on edge");
	Vertex_handle v1 = create_local_view_on( p0, e);
	binop_local_views( v0, v1, BOP, *sncp());
	sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( f, o)) {
	CGAL_NEF_TRACEN("p0 found on facet");
	Vertex_handle v1 = create_local_view_on( p0, f);
	binop_local_views( v0, v1, BOP, *sncp());
	sncp()->delete_vertex(v1);
      }
      else if( CGAL::assign( c, o)) {
	CGAL_NEF_TRACEN("p0 found on volume with mark " << c->mark());
#ifdef CGAL_NEF3_OVERLAY_IF_NEEDED_OFF
        if(true) {
#else
	if( BOP( true, c->mark()) != BOP( false, c->mark())) {
#endif
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
	  Vertex_handle v1 = create_local_view_on( p0, c);
	  binop_local_views( v0, v1, BOP, *sncp());
	  sncp()->delete_vertex(v1);
#else
	  SNC_constructor C(*sncp());
	  Vertex_handle v1 = C.clone_SM(v0);
	  SM_decorator SM(&*v1);
	  SM.change_marks(BOP, c->mark());
	  SM_overlayer O(&*v1);
	  O.simplify();
#endif
	} else {
	  CGAL_NEF_TRACEN("vertex in volume deleted " << std::endl << 
		 "  vertex: " <<  v0->point() << std::endl << 
		 "  mark of volume: " << c->mark()); 
	}
      }
      else CGAL_assertion_msg( 0, "wrong handle");

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.stop();
#endif
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());

    CGAL_NEF_TRACEN("=> for all v1 in snc1, qualify v1 with respect snc0");
    CGAL_forall_vertices( v0, snc2) {
      //      v_qualifying++;
      if(ignore[v0]) continue;
      Point_3 p1(v0->point());
      Halfedge_handle e;
      Halffacet_handle f;
      Volume_handle c;
      CGAL_NEF_TRACEN("Locating point " << p1);

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      number_of_point_location_queries++;
      timer_point_location.start();
#endif 
      Object_handle o = pl1->locate(p1);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
      timer_point_location.stop();
#endif 
      
      CGAL_assertion_code(Vertex_handle v);
      CGAL_assertion( !CGAL::assign( v, o));

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.start();
#endif
      if( CGAL::assign( e, o)) {
	CGAL_NEF_TRACEN("p1 found on edge");
	Vertex_handle v1 = create_local_view_on( p1, e);
	binop_local_views( v1, v0, BOP, *sncp());
	sncp()->delete_vertex(v1);
      } 
      else if( CGAL::assign( f, o)) {
	CGAL_NEF_TRACEN("p1 found on facet");
	Vertex_handle v1 = create_local_view_on( p1, f);
	binop_local_views( v1, v0, BOP, *sncp());
	sncp()->delete_vertex(v1);
      } 
      else if( CGAL::assign( c, o)) {
	CGAL_NEF_TRACEN("p1 found on volume with mark " << c->mark());
#ifdef CGAL_NEF3_OVERLAY_IF_NEEDED_OFF
        if(true) {
#else
	if( BOP( c->mark(), true) != BOP( c->mark(), false)) {
#endif
#ifdef CGAL_NEF3_OVERLAY_BY_HAND_OFF
	Vertex_handle v1 = create_local_view_on( p1, c);
	binop_local_views( v1, v0, BOP, *sncp());
	sncp()->delete_vertex(v1);	  
#else
	  SNC_constructor C(*sncp());
	  Vertex_handle v1 = C.clone_SM(v0);
	  SM_decorator SM(&*v1);
	  SM.change_marks(c->mark(), BOP);
	  SM_overlayer O(&*v1);
	  O.simplify();
#endif
	} else {
	  CGAL_NEF_TRACEN("vertex in volume deleted " << std::endl << 
		 "  vertex: " <<  v0->point() << std::endl << 
		 "  mark of volume: " << c->mark()); 
	}
      }
      else CGAL_assertion_msg( 0, "wrong handle");

#if defined(CGAL_NEF3_TIMER_OVERLAY)
      timer_overlay.stop();
#endif
    }

    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());

    // Each time the intersect method of the point locator for the 
    // SNC structure finds an intersection between the segment defined
    // by an edge on the other SNC structure, the call back method is
    // called with the intersecting objects and the intersection point.
    // The responsability of the call back functor is to construct the
    // local view on the intersection point on both SNC structures,
    // overlay them and add the resulting sphere map to the result.

    // CGAL_NEF_SETDTHREAD(19*509*43*131);

    //    Progress_indicator_clog ee_intersections
    //      (sncp()->number_of_edges(),
    //       "binary_operator: finding edge-edge intersections...");

    Intersection_call_back<SNC_decorator, Selection> call_back0
      ( const_cast<SNC_structure&>(snc1), const_cast<SNC_structure&>(snc2), BOP, *sncp());
    Intersection_call_back<SNC_decorator, Selection> call_back1 
      ( const_cast<SNC_structure&>(snc2), const_cast<SNC_structure&>(snc2), BOP, *sncp(), true);

#ifdef CGAL_NEF3_TIMER_INTERSECTION
    double split_intersection = timer_overlay.time();
    CGAL::Timer timer_intersection;
    timer_intersection.start();
#endif

    // choose between intersection algorithms
#ifdef CGAL_NEF3_INTERSECTION_BY_KDTREE
    Halfedge_iterator e0, e1;
    /*
    CGAL_NEF_TRACEN("=> finding edge-edge intersections...");
    CGAL_forall_edges( e0, snc1) {
      //      ee_intersections++;
      pl2->intersect_with_edges( e0, call_back0);
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());

    //    Progress_indicator_clog ef_intersections
    //      (sncp()->number_of_edges()+snc1i.number_of_edges(),
    //       "binary_operator: finding edge-facet intersections...");

    CGAL_NEF_TRACEN("number of vertices (so far...) = "<< sncp()->number_of_vertices());
    CGAL_NEF_TRACEN("=> finding edge0-facet1 intersections...");
    CGAL_forall_edges( e0, snc1) {
      //      ef_intersections++;
      pl2->intersect_with_facets( e0, call_back0);
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());
    */

    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc1))
      pl2->intersect_with_edges_and_facets(e0,call_back0);


    CGAL_NEF_TRACEN("=> finding edge1-facet0 intersections...");
    CGAL_forall_edges( e1,const_cast<SNC_structure&>(snc2)) {
      //      ef_intersections++;
      pl1->intersect_with_facets( e1, call_back1);
    }
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());
#elif defined CGAL_NEF3_INTERSECTION_NAIVE
    
    CGAL::SNC_point_locator_naive<SNC_decorator> pln1;
    CGAL::SNC_point_locator_naive<SNC_decorator> pln2;
    pln1.initialize(const_cast<SNC_structure*>(&snc1));
    pln2.initialize(const_cast<SNC_structure*>(&snc2));
    Halfedge_iterator e0;   
    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc1))
      pln2.intersect_with_edges_and_facets(e0,call_back0);
    CGAL_forall_edges(e0,const_cast<SNC_structure&>(snc2))
      pln1.intersect_with_facets( e0, call_back1);

#else
    CGAL_NEF_TRACEN("intersection by fast box intersection");
        binop_intersection_test_segment_tree<SNC_decorator> binop_box_intersection;
        binop_box_intersection(call_back0, call_back1, 
			       const_cast<SNC_structure&>(snc1), 
			       const_cast<SNC_structure&>(snc2));
#endif

#ifdef CGAL_NEF3_TIMER_INTERSECTION
    timer_intersection.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_intersection: "
		<< timer_intersection.time()-timer_overlay.time()+split_intersection 
		<< std::endl;
#endif

    CGAL_NEF_TRACEN("=> resultant vertices (before simplification): ");
    CGAL_assertion_code(CGAL_forall_vertices( v0, *sncp()) 
			  CGAL_NEF_TRACEN(&*v0<<" "<<v0->point()));

#ifdef CGAL_NEF3_TIMER_SIMPLIFICATION
    CGAL::Timer timer_simplification;
    timer_simplification.start();
#endif

    SNC_simplify simp(*sncp());
    simp.vertex_simplification(NO_SNC);

#ifdef CGAL_NEF3_TIMER_SIMPLIFICATION
    timer_simplification.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_simplification: " 
		<< timer_simplification.time() << std::endl;
#endif

    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "<< sncp()->number_of_vertices());

    CGAL_NEF_TRACEN("=> resultant vertices (after simplification): ");
    CGAL_assertion_code(CGAL_forall_vertices( v0, *sncp()) 
			  CGAL_NEF_TRACEN(&*v0<<" "<<v0->point()));
#ifdef CGAL_NEF3_DUMP_SNC_OPERATORS
    CGAL_NEF_TRACEN("=> pre-construction result");
    SNC_io_parser<SNC_structure> O(std::cerr, *sncp());
    O.print();
#endif
    
#ifndef CGAL_NEF3_NO_EXTERNAL
    //    CGAL_assertion( !sncp()->is_empty());
    SNC_constructor C(*sncp(), pl0);
    C.build_external_structure();
#endif

#ifdef CGAL_NEF3_DUMP_SNC_OPERATORS
    CGAL_NEF_TRACEN("=> construction completed, result: ");
    SNC_io_parser<SNC_structure> Op(std::cout, *sncp());
    Op.print();
#endif

#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_plane_sweeps: " 
		<< number_of_plane_sweeps << std::endl;
      std::cout << "Runtime_plane_sweeps: "
		<< timer_plane_sweeps.time() << std::endl;
    }
#endif

    CGAL_assertion(!simp.simplify());
    CGAL_NEF_TRACEN("=> end binary operation. ");

#ifdef CGAL_NEF3_DUMP_STATISTICS
    if(cgal_nef3_timer_on) {
      std::cout << "Vertices_in_object_A: "
		<< snc1.number_of_vertices() << std::endl;
      std::cout << "Vertices_in_object_B: "
		<< snc2.number_of_vertices() << std::endl;
      std::cout << "Number_of_intersections: "
		<< number_of_intersections << std::endl;
      std::cout << "Number_of_intersection_candidates: "
		<< number_of_intersection_candidates << std::endl;    
      std::cout << "Vertices_in_Result: "
		<< sncp()->number_of_vertices() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_edge_facet_overlays: " 
		<< number_of_edge_facet_overlays << std::endl;
      std::cout << "Number_of_clones: " 
		<< number_of_clones << std::endl;
      std::cout << "Number_of_sphere_sweeps: " 
		<< number_of_sphere_sweeps << std::endl;
      std::cout << "Runtime_sphere_sweeps: "
		<< timer_sphere_sweeps.time() << std::endl;
    }
#endif

#if defined (CGAL_NEF3_TIMER_SPHERE_SWEEPS) && defined (CGAL_NEF3_TIMER_PLANE_SWEEPS)
    if(cgal_nef3_timer_on) {
      std::cout << "Runtime_all_sweeps: "
		<< timer_sphere_sweeps.time()+timer_plane_sweeps.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_OVERLAY
    if(cgal_nef3_timer_on) {
      std::cout << "Runtime_overlay: " 
		<< timer_overlay.time() << std::endl;
    }
#endif

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    if(cgal_nef3_timer_on) {
      std::cout << "Number_of_ray_shooting_queries: "
		<< number_of_ray_shooting_queries << std::endl;
      std::cout << "Runtime_ray_shooting: "
		<< timer_ray_shooting.time() << std::endl;
      std::cout << "Number_of_point_location_queries: "
		<< number_of_point_location_queries << std::endl;
      std::cout << "Runtime_point_location: "
		<< timer_point_location.time() << std::endl;
      std::cout << "Number_of_kd_tree_queries: "
		<< number_of_point_location_queries + 
	           number_of_ray_shooting_queries << std::endl;
      std::cout << "Runtime_kd_tree_queries: "
		<< timer_point_location.time() + 
                   timer_ray_shooting.time() << std::endl;
    }
#endif 

#ifdef CGAL_NEF3_TIMER_BINARY_OPERATION
    if(cgal_nef3_timer_on) {
      timer_binary_operation.stop();
      std::cout << "Runtime_binary_operation: " 
		<< timer_binary_operation.time() << std::endl;
    }
#endif 
    }
  

  bool is_valid( bool verb = false, int level = 0) {
    
    Verbose_ostream verr(verb);
    verr << "begin CGAL::SNC_decorator<...>::is_valid( verb=true, "
      "level = " << level << "):" << std::endl;
    
    int max = number_of_vertices() 
      + number_of_halfedges() 
      + number_of_halffacets() 
      + number_of_volumes()
      + 2 * sncp()->number_of_shalfedges()
      + sncp()->number_of_shalfloops()
      + sncp()->number_of_sfaces();
    int count = 0;

    bool valid = true;
    Vertex_iterator vi;
    std::list<Point_3> Points(false);   // durch hashmap ersetzen    
    CGAL_forall_vertices(vi,*this) {
      if(!valid) break;

      valid = valid && (vi->sncp()==sncp());
      valid = valid && vi->is_valid(verb, level);

      SM_decorator SD(&*vi);
      Unique_hash_map<SVertex_handle, bool> SVvisited(false);
      Unique_hash_map<SHalfedge_handle, bool> SEvisited(false);
      Unique_hash_map<SFace_handle, bool> SFvisited(false);

      valid = valid && SD.is_valid(SVvisited,SEvisited,SFvisited, verb, level);

      Points.push_back(vi->point());

      //      Mark m1 = SD.mark_of_halfsphere(-1);  
      //      Mark m2 = SD.mark_of_halfsphere(+1);
      //      SM_point_locator PL(vi);
      //      PL.init_marks_of_halfspheres();
      //      valid = valid && (SD.mark_of_halfsphere(-1) == m1);
      //      valid = valid && (SD.mark_of_halfsphere(+1) == m2);
     
      //      CGAL_NEF_TRACEN(m1 << "," << m2 << " should have been " << 
      //	     SD.mark_of_halfsphere(-1) << "," << SD.mark_of_halfsphere(+1));

      valid = valid && (++count <= max);
    }

   verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    Points.sort(points_lt());
    typename std::list<Point_3>::const_iterator li1, li2;
    li2 = li1 = Points.begin();
    if(!Points.empty()) {
      li2++;
      while(valid && li2 != Points.end()) {
	valid = valid && (*li1++ != *li2++);
      }
    }

    verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;
    
    Halfedge_iterator he;
    CGAL_forall_halfedges(he,*this) {

      valid = valid && he->is_valid(verb, level);
      valid = valid && (he->twin()!=he);
      valid = valid && (he->twin()->twin()==he);

      if(he->is_twin()) {
	SM_decorator S1(&*he->source());
	SM_decorator S2(&*he->twin()->source());
	SHalfedge_handle se1(S1.first_out_edge(he));
	SHalfedge_handle se2(S2.first_out_edge(he->twin()));
	if(se1 != NULL && se2 != NULL) {
	  SHalfedge_handle start1(se1);
	  SHalfedge_handle start2(se2->twin()->snext());
	  while(se1->facet() != se2->twin()->facet() && se2 != start2)
	    se2 = se2->sprev()->twin();

	  start2 = se2;
	  do {
	    se1 = se1->twin()->snext();
	    se2 = se2->sprev()->twin();
	    valid = valid && (se1->facet() == se2->facet()->twin());
	  } while(se1 != start1);
	  valid = valid && (se2 == start2);
	}

	//	Line_3 supporting_line(he->source()->point(), he->point());
	//	valid = valid && supporting_line.has_on(he->twin()->source()->point());
      }
      valid = valid && (++count <= max);
    }

    verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    Unique_hash_map<SHalfedge_handle, bool> SEinUniqueFC(false);
    Halffacet_iterator hfi;
    CGAL_forall_halffacets(hfi,*this) {
      valid = valid && hfi->is_valid(verb, level);
      valid = valid && (hfi->twin()!=hfi);
      valid = valid && (hfi->twin()->twin()==hfi);
      valid = valid && (hfi->plane().opposite() == hfi->twin()->plane()); 
     
      Halffacet_cycle_iterator hfci;
      CGAL_forall_facet_cycles_of(hfci,hfi) { 
	if(hfci.is_shalfedge()) {
	  SHalfedge_handle sheh(hfci);
	  valid = valid && (sheh != SHalfedge_handle());
	  SHalfedge_around_facet_circulator shec1(sheh), shec2(shec1);
       	  CGAL_For_all(shec1, shec2) {
	    CGAL_assertion(!SEinUniqueFC[shec1]);
	    SEinUniqueFC[shec1] = true;
	    Plane_3 p_ref(shec1->source()->source()->point(),shec1->circle().opposite().orthogonal_vector());
	    valid = valid && Infi_box::check_point_on_plane(shec1->source()->source()->point(), hfi->plane());
	    valid = valid && (normalized(hfi->plane()) == normalized(p_ref));
	    valid = valid && (sheh->incident_sface()->volume() == 
			      hfi->twin()->incident_volume());
	  }
	}
	else if(hfci.is_shalfloop())
	  valid = valid && (SHalfloop_handle(hfci) != SHalfloop_handle()); 
	else
	  valid = false;
	valid = valid && (++count <= max);	
      }
      
      valid = valid && (++count <= max);
    }

    verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    Volume_iterator voli;
    CGAL_forall_volumes(voli,*this) {
      
      if(number_of_vertices() > 0)
	valid = valid && (voli->is_valid(verb, level));

      Shell_entry_iterator si;
      CGAL_forall_shells_of(si,voli) {
	valid = valid && (si != Shell_entry_iterator() &&
			  SFace_handle(si) != SFace_handle() && 
			  SFace_handle(si) != NULL);
	valid = valid && (++count <= max);
      }
      valid = valid && (++count <= max);
    }

   verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    SHalfedge_iterator she;
    CGAL_forall_shalfedges(she,*sncp()) {
      valid = valid && (she->next() != she);
      valid = valid && (she->prev() != she);
      valid = valid && (she->next()->prev() == she);
      valid = valid && (she->prev()->next() == she);
      valid = valid && (she->next()->twin()->next() == she->twin());
      valid = valid && (she->prev()->twin()->prev() == she->twin());
      valid = valid && (she->facet()->twin() == she->twin()->facet());
      valid = valid && (she->facet()->mark() == she->mark());
      valid = valid && (she->twin()->facet()->incident_volume() == 
			she->incident_sface()->volume());
      valid = valid && (++count <= max);
    }

   verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    SHalfloop_iterator shl;
    CGAL_forall_shalfloops(shl,*sncp()){
      SM_decorator SD;
      valid = valid && (shl->facet()->mark() == shl->mark());
      //      valid = valid && (sh1->twin()->facet()->plane() == sh1->circle()); 
      //      valid = valid && (sh1->twin()->facet()->volume() == sh1->sface()->volume());
    }

    verr << "CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    SFace_iterator sf;
    CGAL_forall_sfaces(sf,*sncp()) {
      SM_decorator SD;
      valid = valid && (sf->volume()->mark() == sf->mark());
      SFace_cycle_iterator sfc;
      for(sfc=sf->sface_cycles_begin();sfc!=sf->sface_cycles_end();++sfc)
	if(sfc.is_shalfedge())
	  valid = valid && (sf==SHalfedge_handle(sfc)->incident_sface());
    }

    verr << "end of CGAL::SNC_decorator<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;

    return valid;
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
  SVertex_iterator   svertices_begin() { return sncp()->svertices_begin(); }
  SVertex_iterator   svertices_end()   { return sncp()->svertices_end(); }
  SHalfedge_iterator shalfedges_begin(){ return sncp()->shalfedges_begin(); }
  SHalfedge_iterator shalfedges_end()  { return sncp()->shalfedges_end(); }
  SHalfloop_iterator shalfloops_begin(){ return sncp()->shalfloops_begin(); }
  SHalfloop_iterator shalfloops_end()  { return sncp()->shalfloops_end(); }
  SFace_iterator   sfaces_begin() { return sncp()->sfaces_begin(); }
  SFace_iterator   sfaces_end()   { return sncp()->sfaces_end(); }

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
      CGAL_forall_facet_cycles_of(fc,f) {
        if (fc.is_shalfedge() ) {
	  SHalfedge_handle e(fc);
          SHalfedge_around_facet_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) { e = ec->twin();
            if ( Done[e->incident_sface()] ) continue;
            SFaceCandidates.push_back(e->incident_sface());
            Done[e->incident_sface()] = true;
          }
        } else if (fc.is_shalfloop()) {
	  SHalfloop_handle l(fc);
	  l = l->twin();
          if ( Done[l->incident_sface()] ) continue;
          SFaceCandidates.push_back(l->incident_sface());
          Done[l->incident_sface()] = true;
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
      }
    }
    if ( !SFaceCandidates.empty() ) {
      SFace_handle sf = *SFaceCandidates.begin();
      SFaceCandidates.pop_front();
      V.visit(sf); // report sface
      if ( !Done[sf->center_vertex()] )
        V.visit(sf->center_vertex()); // report vertex
      Done[sf->center_vertex()] = true;
      //      SVertex_handle sv;
      SM_decorator SD(&*sf->center_vertex());
      /*      
      CGAL_forall_svertices(sv,SD){
	if(SD.is_isolated(sv) && !Done[sv])
	  V.visit(sv);
      }
      */
      SFace_cycle_iterator fc;
      CGAL_forall_sface_cycles_of(fc,sf) {
        if ( fc.is_shalfedge() ) {
	  SHalfedge_handle e(fc);
	  SHalfedge_around_sface_circulator ec(e),ee(e);
          CGAL_For_all(ec,ee) {
	    V.visit(SHalfedge_handle(ec));
            SVertex_handle v = ec->twin()->source();
            if ( !SD.is_isolated(v) && !Done[v] ) {
              V.visit(v); // report edge
              Done[v] = Done[v->twin()] = true;
            }
            Halffacet_handle f = ec->twin()->facet();
            if ( Done[f] ) continue;
            FacetCandidates.push_back(f); Done[f] = true;
          }
        } else if ( fc.is_svertex() ) {
	  SVertex_handle v(fc);
          if ( Done[v] ) continue; 
          V.visit(v); // report edge
	  V.visit(v->twin());
          Done[v] = Done[v->twin()] = true;
	  CGAL_assertion(SD.is_isolated(v));
	  SFaceCandidates.push_back(v->twin()->incident_sface());
	  Done[v->twin()->incident_sface()]=true;
          // note that v is isolated, thus v->twin() is isolated too
	  //	  SM_decorator SD;
	  //	  SFace_handle fo;
	  //	  fo = v->twin()->incident_sface();
	  /*
	  if(SD.is_isolated(v)) 
	    fo = v->source()->sfaces_begin();
	  else
	    fo = v->twin()->incident_sface();
	  */
	  //	  if ( Done[fo] ) continue;
	  //	  SFaceCandidates.push_back(fo); Done[fo] = true;
        } else if (fc.is_shalfloop()) {
	  SHalfloop_handle l(fc);
	  V.visit(l);
          Halffacet_handle f = l->twin()->facet();
          if ( Done[f] ) continue;
          FacetCandidates.push_back(f);  Done[f] = true;
        } else CGAL_assertion_msg(0,"Damn wrong handle.");
      }
    }
  }
}



CGAL_END_NAMESPACE
#endif //CGAL_SNC_DECORATOR_H

