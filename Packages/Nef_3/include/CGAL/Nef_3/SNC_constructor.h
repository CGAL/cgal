// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related docuhmentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_constructor.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_constructor.h       construction of basic SNCs and global construction
// ============================================================================
#ifndef CGAL_SNC_CONSTRUCTOR_H
#define CGAL_SNC_CONSTRUCTOR_H

#include <CGAL/basic.h>
#include <CGAL/functional.h> 
#include <CGAL/function_objects.h> 
#include <CGAL/Circulator_project.h>
#include <CGAL/Nef_3/Normalizing.h>
#include <CGAL/Nef_3/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <CGAL/Nef_3/SNC_ray_shooter.h>
#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#include <map>
#include <list>
#undef _DEBUG
#define _DEBUG 43
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge,Decorator> Self;
  Point p; int i; Edge e;
  Decorator& D;
  Halfedge_key(Point pi, int ii, Edge ei, Decorator& Di ) : 
    p(pi), i(ii), e(ei), D(Di) {}
  Halfedge_key(const Self& k) : p(k.p), i(k.i), e(k.e), D(k.D) {}
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge,Decorator> Key;
  typedef typename Point::R R;
  typedef typename R::Vector_3 Vector;
  typedef typename R::Direction_3 Direction;
  bool operator()( const Key& k1, const Key& k2) const { 
    if ( k1.p == k2.p ) 
      return (k1.i < k2.i);
    /* previous code: 
       else return CGAL::lexicographically_xyz_smaller(k1.p,k2.p); */
    Direction l(Vector(CGAL::ORIGIN, k1.D.tmp_point(k1.e)));
    if( k1.i < 0) l = -l;
    return (Direction( k2.p - k1.p) == l); 
  }
};

template <typename Point, typename Edge, class Decorator>
std::ostream& operator<<(std::ostream& os, 
                         const Halfedge_key<Point,Edge,Decorator>& k )
{ os << k.p << " " << k.i; return os; }

template <typename R>
int sign_of(const CGAL::Plane_3<R>& h)
{ if ( h.a() != 0 ) return CGAL_NTS sign(h.a());
  if ( h.b() != 0 ) return CGAL_NTS sign(h.b());
  return CGAL_NTS sign(h.c());
}

struct Plane_lt {
  template <typename R>
  bool operator()(const CGAL::Plane_3<R>& h1,
                  const CGAL::Plane_3<R>& h2) const
  { 
    typedef typename R::RT     RT;
    RT diff = h1.a()-h2.a();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.b()-h2.b();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.c()-h2.c();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.d()-h2.d(); return CGAL_NTS sign(diff) < 0;
  }
};

// ----------------------------------------------------------------------------
// SNC_constructor 
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_constructor}{SNC}{overlay functionality}{O}}*/

template <typename SNC_structure_>
class SNC_constructor : public SNC_decorator<SNC_structure_>
{ 
public:
  typedef SNC_structure_ SNC_structure;
  typedef typename SNC_structure_::Sphere_kernel        Sphere_kernel;
  typedef typename SNC_structure_::Kernel               Kernel;
  typedef CGAL::SNC_constructor<SNC_structure>          Self;
  typedef CGAL::SNC_decorator<SNC_structure>            Base;
  typedef CGAL::SNC_decorator<SNC_structure>            SNC_decorator;
  typedef CGAL::SNC_ray_shooter<SNC_structure>          SNC_ray_shooter;
  typedef CGAL::SNC_FM_decorator<SNC_structure>         FM_decorator;
  typedef CGAL::SNC_SM_decorator<SNC_structure>         SM_decorator;
  typedef CGAL::SNC_SM_overlayer<SNC_structure>         SM_overlayer;
  typedef CGAL::SNC_SM_point_locator<SNC_structure>     SM_point_locator;
  typedef CGAL::SNC_SM_const_decorator<SNC_structure>   SM_const_decorator;

  #define USING(t) typedef typename SNC_structure::t t
  USING(Vertex);
  USING(Halfedge);
  USING(Halffacet);
  USING(Volume);
  
  USING(Vertex_iterator);
  USING(Halfedge_iterator);
  USING(Halffacet_iterator);
  USING(Volume_iterator);

  USING(Vertex_handle);
  USING(Halfedge_handle);
  USING(Halffacet_handle);
  USING(Volume_handle);

  USING(Vertex_const_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_handle);

  USING(SVertex_iterator);
  USING(SHalfedge_iterator);
  USING(SFace_iterator);
  USING(SHalfloop_iterator);

  USING(SVertex);
  USING(SHalfedge);
  USING(SFace);
  USING(SHalfloop);

  USING(SVertex_handle);
  USING(SHalfedge_handle);
  USING(SFace_handle);
  USING(SHalfloop_handle);

  USING(SVertex_const_handle); 
  USING(SHalfedge_const_handle); 
  USING(SHalfloop_const_handle); 
  USING(SFace_const_handle); 

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
  USING(Vector_3);
  USING(Direction_3);
  USING(Segment_3);
  USING(Line_3);
  USING(Plane_3);

  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Sphere_direction);

  USING(Mark);
  #undef USING

  #define DECUSING(t) typedef typename SM_decorator::t t
  DECUSING(SHalfedge_around_svertex_const_circulator);
  DECUSING(SHalfedge_around_svertex_circulator);
  #undef DECUSING

  typedef void* GenPtr;

  typedef CGAL::Unique_hash_map<SFace_const_handle,unsigned int>  
                                                         Shell_number_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> Shell_closed_hash;

  struct Shell_explorer {
    const SNC_decorator& D;
    Shell_number_hash&  Shell;
    Shell_closed_hash& Closed;
    SFace_visited_hash& Done;
    Vertex_handle v_min;
    int n;

    Shell_explorer(const SNC_decorator& Di, Shell_number_hash& Si, 
                   Shell_closed_hash& Sc, SFace_visited_hash& Vi) 
      : D(Di), Shell(Si), Closed(Sc), Done(Vi), n(0) {}

    void visit(SFace_handle h) { 
      TRACEN("visit sf "<<D.point(D.vertex(h)));
      Shell[h]=n;
      Done[h]=true;
    }
    void visit(Vertex_handle h) { 
      TRACEN("visit v  "<<D.point(h));
      if ( CGAL::lexicographically_xyz_smaller(
           D.point(h),D.point(v_min)) ) 
	v_min = h; 
    }
    void visit(Halfedge_handle h) { 
      TRACEN("visit he "<<D.point(D.source(h)));
      SM_decorator SD(D.vertex(h));
      SFace_handle sf = D.source(h)->sfaces_begin();
      if( Closed[sf] ) {
	if( SD.is_isolated(h) || SD.first_out_edge(h) == SD.last_out_edge(h) ) {
	  Closed[sf] = false;
	}
      }
    }
    void visit(Halffacet_handle h) { /* do nothing */ }

    Vertex_handle& minimal_vertex() { return v_min; }

    void increment_shell_number() { 
      TRACEN("leaving shell "<<n);
      ++n; 
    }
  };

  SNC_constructor(SNC_structure& W) : Base(W) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/


  Vertex_handle create_box_corner(int x, int y, int z,
                                  bool space=true, bool boundary=true) const; 

  Vertex_handle create_extended_box_corner(int x, int y, int z,
                                  bool space=true, bool boundary=true) const;

  /*{\Mop produces the sphere map representing thp,e box corner in
          direction $(x,y,z)$.}*/

  Vertex_handle create_from_facet(Halffacet_handle f,
				  const Point_3& p) const; 
  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $f$. \precond $p$ is part of $f$.}*/

  Vertex_handle create_from_edge(Halfedge_handle e,
				 const Point_3& p) const; 
  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $e$. \precond $p$ is part of $e$.}*/

  void pair_up_halfedges() const;
  /*{\Mop pairs all halfedge stubs to create the edges in 3-space.}*/

  void link_shalfedges_to_facet_cycles() const;
  /*{\Mop creates all non-trivial facet cycles from sedges. 
  \precond |pair_up_halfedges()| was called before.}*/

  void categorize_facet_cycles_and_create_facets() const;
  /*{\Mop collects all facet cycles incident to a facet and creates
  the facets. \precond |link_shalfedges_to_facet_cycles()| was called
  before.}*/

  void create_volumes() const;
  /*{\Mop collects all shells incident to a volume and creates the
  volumes.  \precond |categorize_facet_cycles_and_creating_facets()| was
  called before.}*/

  Volume_handle determine_volume( SFace_handle sf, 
                const std::vector< Vertex_handle>& MinimalVertex, 
                const Shell_number_hash&  Shell ) const
    /*{\Mop determines the volume |C| that a shell |S| pointed by |sf| 
      belongs to.  \precondition |S| separates the volume |C| from an enclosed
      volume.}*/ {
    TRACEN("determine volume");
    Vertex_handle v_min = MinimalVertex[Shell[sf]];
    Halffacet_handle f_below = get_facet_below(v_min);
    if ( f_below == Halffacet_handle())
      return Base(*this).volumes_begin();
    Volume_handle c = volume(f_below);
    if( c != Volume_handle()) {
      TRACE( "Volume " << &*c << " hit ");
      TRACEN("(Shell #" << Shell[adjacent_sface(f_below)] << ")");
      return c;
    }
    SFace_handle sf_below = adjacent_sface(f_below);
    TRACE( "Shell not assigned to a volume hit ");
    TRACEN( "(Inner shell #" << Shell[sf_below] << ")");
    c = determine_volume( sf_below, MinimalVertex, Shell);
    link_as_inner_shell( sf_below, c);
    return c;
  }

  Halffacet_handle get_facet_below( Vertex_handle vi) const 
     /*{\Mop determines the facet below a vertex |vi| via ray shooting. }*/ {
    Halffacet_handle f_below;
    Point_3 p = point(vi);
    if(!Kernel::is_standard(p))
      return Halffacet_handle();

    Segment_3 ray( p, Kernel::epoint( 0, p.hx()[0], 0,p.hy()[0], 1,0, p.hw()[0]));

    SNC_ray_shooter rs(*sncp());
    Object_handle o = rs.shoot(ray);
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    TRACEN("get_facet_below");
    if( assign(v, o)) {
      TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if( f_below == Halffacet_handle())
	f_below = get_facet_below(v);
    }
    else if( assign(e, o)) {
      TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if( f_below == Halffacet_handle())
	f_below = get_facet_below(vertex(e));
    }
    else if( assign(f, o)) {
      TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_nef3_assertion( f_below != Halffacet_handle());
    }
    else { TRACEN("no facet below found..."); }
    return f_below;
  }
  
  void create_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, Plane_3 h, bool boundary) const;
  void create_corner_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, int max, bool boundary) const;
  void create_degenerate_corner_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, int min, int max, bool boundary) const;

}; // SNC_constructor<SNC>

template <typename SNC_>
void
SNC_constructor<SNC_>::
create_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, Plane_3 h, bool boundary) const { 

  int max = 0;
  if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hy()))
    max = 1;
  if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hz()))
    max = 2;

  TRACEN("create frame point " << p << std::endl << sp1 << std::endl << sp2);
  Vertex_handle v=sncp()->new_vertex(p, boundary);
  SM_decorator SD(v); 

  TRACEN("create spoints");
  Sphere_point SP[4];
  switch(max) {
  case 0: SP[0] = Sphere_point(1,0,0); break;
  case 1: SP[0] = Sphere_point(0,1,0); break;
  case 2: SP[0] = Sphere_point(0,0,1); break;
  default: CGAL_nef3_assertion_msg(0,"wrong value");
  }
  
  SP[1]=sp1;
  SP[2]=sp2;
  
  if (spherical_orientation(SP[0],SP[1],SP[2]) < 0) {
    SP[3] = SP[0];
    SP[0] = -Vector_3(SP[3]);
  }
  else
    SP[3] = -Vector_3(SP[0]);


  typedef typename Kernel::RT RT;
  RT delta = h.a()*SP[0].hx()+h.b()*SP[0].hy()+h.c()*SP[0].hz();
  CGAL_nef3_assertion(delta !=0);
  bool fmark0 = (delta <  0);

  TRACEN("create svertices");
  SVertex_handle sv[4];
  for(int i=0; i<4; ++i) {
    sv[i] = SD.new_vertex(SP[i]);
    mark(sv[i]) = boundary;
  }
  
  TRACEN("create sedges");
  SHalfedge_handle she[5];
  for(int si=0; si<3;++si) {
    she[si]=SD.new_edge_pair(sv[si], sv[(si+1)%3]);
    SD.circle(she[si])= Sphere_circle(Plane_3(SP[si],SP[(si+1)%3],Point_3(0,0,0)));
    SD.circle(SD.twin(she[si])) =  SD.circle(she[si]).opposite();
    TRACEN("circles " << SD.circle(she[si]) << "    " << SD.circle(she[si]).opposite() << ": " << SP[si] << "          " << SP[(si+1)%3]);	  
    SD.mark(she[si]) = boundary;
  }

  she[3] = SD.new_edge_pair(sv[3],she[1],-1);
  she[4] = SD.new_edge_pair(sv[3],SD.twin(she[1]),1);
  
  SD.circle(SD.twin(she[3]))= Sphere_circle(Plane_3(SP[1],SP[3],Point_3(0,0,0)));
  SD.circle(she[3]) =  SD.circle(SD.twin(she[3])).opposite();	  
  TRACEN("circles2 " << SD.circle(she[3]) << "    " << SD.circle(she[3]).opposite());	  
  SD.mark(she[3]) = boundary;
  // SD.mark(SD.twin(she[3])) = !fmark0;
  
  SD.circle(she[4])= Sphere_circle(Plane_3(SP[3],SP[2],Point_3(0,0,0)));
  SD.circle(SD.twin(she[4])) =  SD.circle(she[4]).opposite();
  TRACEN("circles3 " << SD.circle(she[4]) << "    " << SD.circle(she[4]).opposite());	  	  
  SD.mark(she[4]) = boundary;
  // SD.mark(SD.twin(she[4])) = !fmark0;
  
  SFace_handle sf[3];
  for(int i=0; i<3; ++i)
    sf[i] = SD.new_face();
  SD.mark(sf[0])= fmark0;
  SD.mark(sf[1])=!fmark0;
  SD.mark(sf[2])= 0;

  SD.link_as_face_cycle(she[0],sf[0]);
  SD.link_as_face_cycle(SD.twin(she[1]),sf[1]);
  SD.link_as_face_cycle(SD.twin(she[0]),sf[2]);

  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
}

template <typename SNC_>
void
SNC_constructor<SNC_>::
create_corner_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, int max, bool boundary) const { 

  TRACEN("create corner frame point " << p << std::endl << sp1 << std::endl << sp2);
  Vertex_handle v=sncp()->new_vertex(p, boundary);
  SM_decorator SD(v); 

  TRACEN("create sfaces");
  SFace_handle sf[3];
  for(int i=0; i<3; ++i)
    sf[i] = SD.new_face();
  SD.mark(sf[0])=SD.mark(sf[1])=1;
  SD.mark(sf[2])=0;

  typename Kernel::RT::NT vp[3];
  vp[0] = -v->point().hx()[1];
  vp[1] = -v->point().hy()[1];
  vp[2] = -v->point().hz()[1];
  
  TRACEN("create spoints");
  Sphere_point SP[5];
  switch(max) {
  case 0: SP[0] = Sphere_point(0,vp[1],0); SP[3]= Sphere_point(0,0,vp[2]); SP[4] = Sphere_point(vp[0],0,0); break;
  case 1: SP[0] = Sphere_point(vp[0],0,0); SP[3]= Sphere_point(0,0,vp[2]); SP[4] = Sphere_point(0,vp[1],0); break;
  case 2: SP[0] = Sphere_point(vp[0],0,0); SP[3]= Sphere_point(0,vp[1],0); SP[4] = Sphere_point(0,0,vp[2]); break;
  default: CGAL_nef3_assertion_msg(0,"wrong value");
  }
  
  if (spherical_orientation(SP[0],Sphere_point(sp1),Sphere_point(sp2)) > 0) {
    SP[1] = sp1;
    SP[2] = sp2;
  }
  else {
    SP[1] = sp2;
    SP[2] = sp1;
  }

  if (spherical_orientation(SP[3],SP[0],SP[1]) < 0) {
    Sphere_point sx = SP[3];
    SP[3] = SP[0];
    SP[0] = sx;
  }

  TRACEN("create svertices");
  SVertex_handle sv[5];
  for(int i=0; i<5; ++i) {
    TRACEN("svertex " << SP[i]);
    sv[i] = SD.new_vertex(SP[i]);
    mark(sv[i]) = boundary;
  }
  
  TRACEN("create sedges");
  SHalfedge_handle she[6];
  for(int si=0; si<4;++si) {
    she[si]=SD.new_edge_pair(sv[si], sv[(si+1)%4]);
    SD.circle(she[si])= Sphere_circle(Plane_3(SP[si],SP[(si+1)%4],Point_3(0,0,0)));
    SD.circle(SD.twin(she[si])) =  SD.circle(she[si]).opposite();
    //    TRACEN("circles " << SD.circle(she[si]) << "    " << SD.circle(she[si]).opposite());
    SD.mark(she[si]) = boundary;
  }

  TRACEN("1");
  she[4] = SD.new_edge_pair(sv[4],she[1],1);
  she[5] = SD.new_edge_pair(sv[4],SD.twin(she[1]),-1);
  
  TRACEN("2");
  SD.circle(SD.twin(she[4]))= Sphere_circle(Plane_3(SP[1],SP[4],Point_3(0,0,0)));
  SD.circle(she[4]) =  SD.circle(SD.twin(she[4])).opposite();
  //  TRACEN("circles2 " << SD.circle(she[4]) << "    " << SD.circle(she[4]).opposite());	  
  SD.mark(she[4]) = boundary;
  
  TRACEN("3");
  SD.circle(she[5])= Sphere_circle(Plane_3(SP[4],SP[2],Point_3(0,0,0)));
  SD.circle(SD.twin(she[5])) =  SD.circle(she[5]).opposite();
  //  TRACEN("circles2 " << SD.circle(she[5]) << "    " << SD.circle(she[5]).opposite());		  
  SD.mark(she[5]) = 1;
  
  TRACEN("link faces");
  SD.link_as_face_cycle(she[0],sf[0]);
  SD.link_as_face_cycle(SD.twin(she[0]),sf[1]);
  SD.link_as_face_cycle(she[1],sf[2]);

  TRACEN("fpwc sedges");
  SHalfedge_handle se;
  CGAL_nef3_forall_shalfedges(se,SM_decorator(v))
    TRACEN(PH(se));

  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
}

template <typename SNC_>
void
SNC_constructor<SNC_>::
create_degenerate_corner_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, int min, int max, bool boundary) const { 

  TRACEN("create degenerate corner frame point " << p << std::endl << sp1 << std::endl << sp2);
  Vertex_handle v=sncp()->new_vertex(p, boundary);
  SM_decorator SD(v); 

  TRACEN("create sfaces");
  SFace_handle sf[3];
  for(int i=0; i<3; ++i)
    sf[i] = SD.new_face();
  SD.mark(sf[0])=SD.mark(sf[1])=1;
  SD.mark(sf[2])=0;

  typename Kernel::RT::NT vp[3];
  vp[0] = -v->point().hx()[1];
  vp[1] = -v->point().hy()[1];
  vp[2] = -v->point().hz()[1];
  
  TRACEN("create spoints");
  Sphere_point SP[4];

  switch(max) { 
  case 0: SP[0] = Sphere_point(vp[0],0,0); break;
  case 1: SP[0] = Sphere_point(0,vp[1],0); break;
  case 2: SP[0] = Sphere_point(0,0,vp[2]); break;
  default: CGAL_nef3_assertion_msg(0,"wrong value \"max\"");
  }

  switch(min+max) {
  case 1: SP[3] = Sphere_point(0,0,vp[2]); break;
  case 2: SP[3] = Sphere_point(0,vp[1],0); break;
  case 3: SP[3] = Sphere_point(vp[0],0,0); break;
  default: CGAL_nef3_assertion_msg(0,"wrong value \"min+max\"");
  }
  
  if (spherical_orientation(SP[0],Sphere_point(sp1),Sphere_point(sp2)) > 0) {
    SP[1] = sp1;
    SP[2] = sp2;
  }
  else {
    SP[1] = sp2;
    SP[2] = sp1;
  }

  TRACEN("create svertices");
  SVertex_handle sv[4];
  for(int i=0; i<4; ++i) {
    sv[i] = SD.new_vertex(SP[i]);
    mark(sv[i]) = boundary;
  }
  
  TRACEN("create sedges");
  SHalfedge_handle she[5];
  for(int si=0; si<3;++si) {
    she[si]=SD.new_edge_pair(sv[si], sv[(si+1)%3]);
    SD.circle(she[si])= Sphere_circle(Plane_3(SP[si],SP[(si+1)%3],Point_3(0,0,0)));
    SD.circle(SD.twin(she[si])) =  SD.circle(she[si]).opposite();
    TRACEN("circles " << SD.circle(she[si]) << "    " << SD.circle(she[si]).opposite() << ": " << SP[si] << "          " << SP[(si+1)%3]);	  
    SD.mark(she[si]) = boundary;
  }

  she[3] = SD.new_edge_pair(sv[3],she[1],-1);
  she[4] = SD.new_edge_pair(sv[3],SD.twin(she[1]),1);
  
  SD.circle(SD.twin(she[3]))= Sphere_circle(Plane_3(SP[1],SP[3],Point_3(0,0,0)));
  SD.circle(she[3]) =  SD.circle(SD.twin(she[3])).opposite();	  
  TRACEN("circles2 " << SD.circle(she[3]) << "    " << SD.circle(she[3]).opposite());	  
  SD.mark(she[3]) = boundary;
  
  SD.circle(she[4])= Sphere_circle(Plane_3(SP[3],SP[2],Point_3(0,0,0)));
  SD.circle(SD.twin(she[4])) =  SD.circle(she[4]).opposite();
  TRACEN("circles3 " << SD.circle(she[4]) << "    " << SD.circle(she[4]).opposite());	  	  
  SD.mark(she[4]) = boundary;

  SD.link_as_face_cycle(SD.twin(she[0]),sf[0]);
  SD.link_as_face_cycle(SD.twin(she[1]),sf[1]);
  SD.link_as_face_cycle(she[0],sf[2]);

  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
}

// ----------------------------------------------------------------------------
// create_box_corner()
// creates the local graph at the corner of a cube in direction (x,y,z)
// 'space' specifies if the bounded volume is selected.
// 'boundary' specifies if the boundary of the box is selected 

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_box_corner(int x, int y, int z, bool space, bool boundary) const { 
  CGAL_nef3_assertion(CGAL_NTS abs(x) == CGAL_NTS abs(y) &&
		      CGAL_NTS abs(y) == CGAL_NTS abs(z));
  TRACEN("  constructing box corner on "<<Point_3(x,y,z)<<"...");
  Vertex_handle v = sncp()->new_vertex( Point_3(x, y, z), boundary);
  SM_decorator SD(v);
  Sphere_point sp[] = { Sphere_point(-x, 0, 0), 
			Sphere_point(0, -y, 0), 
			Sphere_point(0, 0, -z) };
  /* create box vertices */
  SVertex_handle sv[3];
  for(int vi=0; vi<3; ++vi) {
    sv[vi] = SD.new_vertex(sp[vi]);
    mark(sv[vi]) = boundary;
  }
  /* create facet's edge uses */
  Sphere_segment ss[3];
  SHalfedge_handle she[3];
  for(int si=0; si<3; ++si) {
    she[si] = SD.new_edge_pair(sv[si], sv[(si+1)%3]);
    ss[si] = Sphere_segment(sp[si],sp[(si+1)%3]);
    SD.circle(she[si]) = ss[si].sphere_circle();
    SD.circle(SD.twin(she[si])) = ss[si].opposite().sphere_circle();
    SD.mark(she[si]) = boundary;
  }
  /* create facets */
  SFace_handle fi = SD.new_face();
  SFace_handle fe = SD.new_face();
  SD.link_as_face_cycle(she[0], fi);
  SD.link_as_face_cycle(SD.twin(she[0]), fe);
  /* set face mark */
  Sphere_point p1 = SD.point(SD.source(she[0]));
  Sphere_point p2 = SD.point(SD.target(she[0]));
  Sphere_point p3 = SD.point(SD.target(SD.next(she[0])));
  if ( spherical_orientation(p1,p2,p3) > 0 ) {
    SD.mark(fi) = space;
    SD.mark(fe) = 0;
  }
  else {
    SD.mark(fi) = 0;
    SD.mark(fe) = space;
  }

  // SD.mark_of_halfsphere(-1) = (x<0 && y>0 && z>0);
  // SD.mark_of_halfsphere(+1) = (x>0 && y>0 && z<0);
  /* TODO: to check if the commented code above could be wrong */
  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
  return v;
}

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_extended_box_corner(int x, int y, int z, bool space, bool boundary) const { 

  CGAL_nef3_assertion(CGAL_NTS abs(x) == CGAL_NTS abs(y) &&
		      CGAL_NTS abs(y) == CGAL_NTS abs(z));

  TRACEN("  constructing box corner on "<<Point_3(x,y,z)<<"...");
  Point_3 p = SNC_::Kernel::epoint(x,0,y,0,z,0,1);
  Vertex_handle v = sncp()->new_vertex(p , boundary);
  SM_decorator SD(v);
  Sphere_point sp[] = { Sphere_point(-x, 0, 0), 
			Sphere_point(0, -y, 0), 
			Sphere_point(0, 0, -z) };
  
  /* create box vertices */
  SVertex_handle sv[3];
  for(int vi=0; vi<3; ++vi) {
    sv[vi] = SD.new_vertex(sp[vi]);
    mark(sv[vi]) = boundary;
  }
  /* create facet's edge uses */
  //  Sphere_segment ss[3];
  SHalfedge_handle she[3];
  for(int si=0; si<3; ++si)
    she[si] = SD.new_edge_pair(sv[si], sv[(si+1)%3]);
  
  for(int i=0; i<3;++i) {
    SD.circle(she[i]) = Sphere_circle(Plane_3(sp[i],sp[(i+1)%3],Point_3(0,0,0)));
    SD.circle(SD.twin(she[i])) =  SD.circle(she[i]).opposite();
    SD.mark(she[i]) = boundary;
  }

  /* create facets */
  SFace_handle fi = SD.new_face();
  SFace_handle fe = SD.new_face();
  SD.link_as_face_cycle(she[0], fi);
  SD.link_as_face_cycle(SD.twin(she[0]), fe);

  Sphere_point p1 = SD.point(SD.source(she[0]));
  Sphere_point p2 = SD.point(SD.target(she[0]));
  Sphere_point p3 = SD.point(SD.target(SD.next(she[0])));
  if ( spherical_orientation(p1,p2,p3) > 0 ) {
    SD.mark(fi) = space;
    SD.mark(fe) = 0;
  }
  else {
    SD.mark(fi) = 0;
    SD.mark(fe) = space;
  }

  // SD.mark_of_halfsphere(-1) = (x<0 && y>0 && z>0);
  // SD.mark_of_halfsphere(+1) = (x>0 && y>0 && z<0);
  /* TODO: to check if the commented code above could be wrong */
  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
  return v;
}


// ----------------------------------------------------------------------------
// create_from_facet() 
// Creates the local graph of a facet f at point p.
// Precondition is that p ist part of f.

#define CGAL_NEF3_BUGGY_CODE

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_from_facet(Halffacet_handle f, const Point_3& p) const
{ 
  /* TODO: CGAL_nef3_assertion(FM_decorator(f).contains(p));*/
  Vertex_handle v = sncp()->new_vertex( p, mark(f));
  point(v) = p;
  Sphere_circle c(plane(f)); // circle through origin parallel to h
  SM_decorator D(v);
  SHalfloop_handle l = D.new_loop_pair();
  SFace_handle f1 = D.new_face(), f2 = D.new_face();
  D.link_as_loop(l,f1);
  D.link_as_loop(twin(l),f2);

  D.circle(l) = c; 
  D.circle(twin(l)) = c.opposite();
  D.mark(f1) = mark(volume(f));
  D.mark(f2) = mark(volume(twin(f)));
  D.mark(l) = mark(f);
#ifdef CGAL_NEF3_BUGGY_CODE
  Sphere_point q(0,-1,0);
  CGAL::Oriented_side os = c.oriented_side(q);
  switch ( os ) {
    case ON_POSITIVE_SIDE: 
      D.mark_of_halfsphere(-1) = D.mark_of_halfsphere(+1) = true;
      break;
    case ON_NEGATIVE_SIDE:
      D.mark_of_halfsphere(-1) = D.mark_of_halfsphere(+1) = false;
      break;
    case ON_ORIENTED_BOUNDARY:
      if ( c.a()<=0 && c.c()>=0 ) // normal(c) dx<=0&&dz>=0
        D.mark_of_halfsphere(+1) = true;
      if ( c.a()>=0 && c.c()<=0 ) // normal(c) dx<=0&&dz>=0
        D.mark_of_halfsphere(-1) = true;
  }
  /* TODO: to find why the code chuck above is wrong */
#endif // CGAL_NEF3_BUGGY_CODE
  SM_point_locator L(v);
  L.init_marks_of_halfspheres();
  return v;
}


// ----------------------------------------------------------------------------
// create_from_edge()
// Creates the local graph of an edge e at point p.
// Precondition is that p ist part of segment(e).

template <typename SNC_>
typename SNC_::Vertex_handle 
SNC_constructor<SNC_>::
create_from_edge(Halfedge_handle e,
		 const Point_3& p) const
{ CGAL_nef3_assertion(segment(e).has_on(p));
  Vertex_handle v = sncp()->new_vertex( p, mark(e));
  SM_decorator D(v);
  SM_const_decorator E(source(e));
  Sphere_point ps = calc_point(e);
  SVertex_handle v1 = D.new_vertex(ps);
  SVertex_handle v2 = D.new_vertex(ps.antipode());
  D.mark(v1) = D.mark(v2) = mark(e);
  bool first = true;

  //  SETDTHREAD(19*43*131);
  SM_decorator EE(vertex(e));
  SHalfedge_handle eee;
  TRACEN("---------------------" << point(vertex(e)));
  CGAL_nef3_forall_shalfedges(eee,EE)
    TRACEN("|" << EE.circle(eee) <<"|" << EE.mark(eee) << " " << EE.mark(EE.face(eee)));
  TRACEN(" ");

  SHalfedge_around_svertex_const_circulator ec1(E.out_edges(e)), ee(ec1);
  SHalfedge_handle e1,e2;
  CGAL_For_all(ec1,ee) {
    if (first) e1 = D.new_edge_pair(v1,v2);
    else       e1 = D.new_edge_pair(e1, e2, SM_decorator::AFTER, 
				            SM_decorator::BEFORE);
    e2 = D.twin(e1); 
    first = false;
  }

  CGAL_nef3_forall_shalfedges(eee,D)
    TRACEN("|" << D.circle(eee));
  TRACEN(" ");


  ec1 = E.out_edges(e);
  SHalfedge_around_svertex_circulator ec2(D.out_edges(v1));
  CGAL_For_all(ec1,ee) {
    TRACEN("|" << E.circle(ec1) <<"|" << E.mark(ec1) << " " << E.mark(E.face(ec1)));
    D.mark(ec2) = E.mark(ec1);
    D.circle(ec2) = E.circle(ec1);
    D.circle(D.twin(ec2)) = E.circle(E.twin(ec1));
    SFace_handle f = D.new_face();
    D.link_as_face_cycle(ec2,f);
    D.mark(f) = E.mark(E.face(ec1));
    ++ec2;
  }
  TRACEN(" ");
  SM_point_locator L(v);
  L.init_marks_of_halfspheres();

  CGAL_nef3_forall_shalfedges(eee,D)
    TRACEN("|" << D.circle(eee) <<"|" << D.mark(eee) << " " << D.mark(D.face(eee)));
  TRACEN("---------------------");

  return v;
}

// ----------------------------------------------------------------------------
// pair_up_halfedges()
// Starting from all local graphs of all vertices of a nef polyhedron
// we pair up all halfedges to halfedge pairs. 

template <typename SNC_>
void SNC_constructor<SNC_>::
pair_up_halfedges() const
{ 
  //  SETDTHREAD(43*61);
  TRACEN(">>>>>pair_up_halfedges");
  typedef Halfedge_key< Point_3, Halfedge_handle, SNC_decorator>
    Halfedge_key;
  typedef Halfedge_key_lt< Point_3, Halfedge_handle, SNC_decorator> 
    Halfedge_key_lt;
  typedef std::list<Halfedge_key>  Halfedge_list;

  typedef CGAL::Pluecker_line_3<Kernel> Pluecker_line_3;
  typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
  typedef std::map< Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
    Pluecker_line_map;

  SNC_decorator D(*this);
  Pluecker_line_map M;
  Pluecker_line_map M2;

  Halfedge_iterator e;
  CGAL_nef3_forall_halfedges(e,*sncp()) {
    Point_3 p = point(vertex(e));
    Pluecker_line_3 l(p, p + Vector_3(tmp_point(e)));
    int inverted;
    l = categorize(l,inverted);
    if(D.is_edge_on_infibox(e))
      M2[l].push_back(Halfedge_key(p,inverted,e,D));
    else
      M[l].push_back(Halfedge_key(p,inverted,e,D));
    TRACEN(" ("<<p<<") ("<<p+ Vector_3(tmp_point(e))<<") "<< " (" << tmp_point(e) << ") " << &*e<<" |"<<l << " " << inverted);
  }

  typename Pluecker_line_map::iterator it;

  CGAL_nef3_forall_iterators(it,M2) {
    it->second.sort(Halfedge_key_lt());
    TRACEN("search opposite  "<<it->first<< "\n    "); 
    typename Halfedge_list::iterator itl;
    CGAL_nef3_forall_iterators(itl,it->second) {
      Halfedge_handle e1 = itl->e;
      ++itl; CGAL_nef3_assertion(itl != it->second.end());
      Halfedge_handle e2 = itl->e;
      TRACEN(point(vertex(e1)) << "-" << point(vertex(e2)));
      CGAL_nef3_assertion(tmp_point(e1)==tmp_point(e2).antipode());
      make_twins(e1,e2);
      CGAL_nef3_assertion(mark(e1)==mark(e2));

      // discard temporary sphere_point ?
    }
  }

  CGAL_nef3_forall_iterators(it,M) {
    it->second.sort(Halfedge_key_lt());
    TRACEN("search opposite  "<<it->first<< "\n    "); 
    typename Halfedge_list::iterator itl;
    CGAL_nef3_forall_iterators(itl,it->second) {
      Halfedge_handle e1 = itl->e;
      ++itl; CGAL_nef3_assertion(itl != it->second.end());
      Halfedge_handle e2 = itl->e;
      // TRACEN("e1="<<tmp_point(e1)<<"@"<<point(vertex(e1))<<
      //     " & e2="<<tmp_point(e2)<<"@"<<point(vertex(e2)));
      CGAL_nef3_assertion(tmp_point(e1)==tmp_point(e2).antipode());
      make_twins(e1,e2);
      CGAL_nef3_assertion(mark(e1)==mark(e2));

      // discard temporary sphere_point ?
    }
  }
}

// ----------------------------------------------------------------------------
// link_shalfedges_to_facet_cycles()
// links all edge-uses to facets cycles within the corresponding planes

template <typename SNC_>
void SNC_constructor<SNC_>::
link_shalfedges_to_facet_cycles() const
{
  TRACEN(">>>>>link_shalfedges_to_facet_cycles");
  Halfedge_iterator e;
  CGAL_nef3_forall_edges(e,*sncp()) {
    Halfedge_iterator et = twin(e);
    SM_decorator D(vertex(e)), Dt(vertex(et));
    if ( D.is_isolated(e) ) continue;
    SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
    SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);

    CGAL_For_all(cet,cete)
      if ( Dt.circle(cet) == D.circle(ce).opposite() && 
           twin(Dt.source(cet)) == D.source(ce) ) 
          break;

    /*    DEBUG 

    if( Dt.circle(cet) != D.circle(ce).opposite() )
      TRACEN("assertion failed!");
      
      TRACEN("vertices " << point(vertex(e)) << "    " << point(vertex(et)) << std::endl);
      
      SHalfedge_around_svertex_circulator sc(D.first_out_edge(e));
      SHalfedge_around_svertex_circulator sct(Dt.first_out_edge(et));
      CGAL_For_all(sc,cee)
	TRACEN("sseg@E addr="<<&*sc<<
	       " src="<<D.point(D.source(sc))<<
	       " tgt="<<D.point(D.target(sc))<<endl<<
	       " circle=" << D.circle(sc));
      CGAL_For_all(sct,cete)
	TRACEN("sseg@ET addr="<<&*sct<<
	       " src="<<Dt.point(Dt.source(sct))<<
	       " tgt="<<Dt.point(Dt.target(sct))<<endl<<
	       " circle=" << Dt.circle(sct));
#ifdef SM_VISUALIZOR
      typedef CGAL::SNC_SM_visualizor<SNC_structure> SMV;
      CGAL::OGL::add_sphere();
      SMV V(vertex(e), CGAL::OGL::spheres_.back());
      V.draw_map();
      SMV Vt(vertex(et), CGAL::OGL::spheres_.back());
      Vt.draw_map();
      CGAL::OGL::start_viewer();
      char c;
      cin >> c;
#endif
DEBUG*/

    CGAL_nef3_assertion( Dt.circle(cet) == D.circle(ce).opposite() ); 
    CGAL_nef3_assertion( twin(Dt.source(cet)) == D.source(ce)); 
    SNC_io_parser<SNC_structure> Op(std::cerr, *sncp());
    CGAL_For_all(ce,cee) { 
      CGAL_nef3_assertion(ce->tmp_mark()==cet->tmp_mark());
      link_as_prev_next_pair(Dt.twin(cet),ce);
      link_as_prev_next_pair(D.twin(ce),cet);
      --cet; // ce moves ccw, cet moves cw
    }
  }
}

// ----------------------------------------------------------------------------
// categorize_facet_cycles_and_create_facets()
// sweeping all edge-uses we categorize facet cycle incidence, create
// the facet objects and assign the facet cycles.

template <typename SNC_>
void SNC_constructor<SNC_>::
categorize_facet_cycles_and_create_facets() const
{ TRACEN(">>>>>categorize_facet_cycles_and_create_facets");

  typedef std::list<SObject_handle> SObject_list;
  typedef std::map<Plane_3, SObject_list, Plane_lt> 
    Map_planes;

  Map_planes M;
  SHalfedge_iterator e;
  CGAL_nef3_forall_shalfedges(e,*sncp()) {
    Sphere_circle c(tmp_circle(e));
    Plane_3 h = c.plane_through(point(vertex(e)));
    SM_decorator SD(vertex(e));
    TRACEN(point(target(SD.source(e))) << " - " << point(vertex(e)) << " - " << point(target(SD.target(e))) << " has plane " << h << " has circle " << tmp_circle(e) << " has signum " << sign_of(h));
    if ( sign_of(h)<0 ) continue;
    //   CGAL_nef3_assertion( h == normalized(h));
    M[normalized(h)].push_back(SObject_handle(twin(e)));
    TRACEN(" normalized as " << normalized(h)); 
  }
  SHalfloop_iterator l;
  CGAL_nef3_forall_shalfloops(l,*sncp()) {
    Sphere_circle c(tmp_circle(l));
    Plane_3 h = c.plane_through(point(vertex(l))); 
    if ( sign_of(h)<0 ) continue;
    // CGAL_nef3_assertion( h == normalized(h));
    M[normalized(h)].push_back(SObject_handle(twin(l)));
  }

  typename Map_planes::iterator it;
  CGAL_nef3_forall_iterators(it,M) { 
    TRACEN("  plane "<<it->first<<"             "<<(it->first).point());
    FM_decorator D(*sncp());
    D.create_facet_objects(it->first,it->second.begin(),it->second.end());
  }
}

// ----------------------------------------------------------------------------
// create_volumes()
// categorizes all shells and creates volume objects.

template <typename SNC_>
void SNC_constructor<SNC_>::
create_volumes() const
{ 
  TRACEN(">>>>>create_volumes");
  Shell_number_hash  Shell(-1);
  Shell_closed_hash Closed(true);
  SFace_visited_hash Done(false);
  Shell_explorer V(*this,Shell,Closed,Done);
  std::vector<Vertex_handle> MinimalVertex;
  std::vector<SFace_handle> EntrySFace;
  SFace_iterator f;
  /* First, we classify all the Shere Faces per Shell.  For each Shell we
     determine its minimum lexicographyly vertex and we check wheter the
     Shell encloses a region (closed surface) or not. */
  CGAL_nef3_forall_sfaces(f,*sncp()) { 
    if ( Done[f] ) 
      continue;
    V.minimal_vertex() = vertex(f);
    visit_shell_objects(f,V);

    MinimalVertex.push_back(V.minimal_vertex());
    EntrySFace.push_back(f);
    V.increment_shell_number();
  }
  /* then, we determine the Shells which correspond to Volumes via a ray
     shootting in the direction (-1,0,0) over the Sphere_map of the minimal 
     vertex.  The Shell corresponds to a Volume if the object hit belongs 
     to another Shell. */
  sncp()->new_volume(); // outermost volume (nirvana)
  for( unsigned int i = 0; i < MinimalVertex.size(); ++i) {
    Vertex_handle v = MinimalVertex[i];
    TRACEN( "Shell #" << i << " minimal vertex: " << point(v));
    SM_point_locator D(v);
    SObject_handle o = D.locate(Sphere_point(-1,0,0));
    SFace_const_handle sfc;
    if( !assign(sfc, o) || Shell[sfc] != i) { /*UNTESTED CASE: !assign(sfc,o)*/
      SFace_handle f = EntrySFace[i];
      CGAL_nef3_assertion( Shell[EntrySFace[i]] == i );
      if( Closed[f] ) {
	TRACEN("Shell #" << i << " is closed");
	SM_decorator SD(v);
	Volume_handle c = sncp()->new_volume();
	mark(c) = SD.mark(f);
	link_as_outer_shell(f, c );
	TRACE( "Shell #" << i <<" linked as outer shell");
	TRACEN( "(sface" << (assign(sfc,o)?"":" not") << " hit case)");
      }
    }
  }
  /* finaly, we go through all the Shells which do not correspond to a Volume 
     and we assign them to its enclosing Volume determined via a facet below
     check. */
  CGAL_nef3_forall_sfaces(f,*sncp()) {
    if ( volume(f) != Volume_handle() ) 
      continue;
    TRACEN( "Inner shell #" << Shell[f] << " volume?");
    Volume_handle c = determine_volume( f, MinimalVertex, Shell );
    link_as_inner_shell( f, c );
  }
}


CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONSTRUCTOR_H

