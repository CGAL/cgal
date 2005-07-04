// Copyright (c) 1997-2002,2005 Max-Planck-Institute Saarbruecken (Germany).
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
//                 Ralf Osbild     <osbild@mpi-sb.mpg.de>
#ifndef CGAL_NEF_POLYHEDRON_3_H
#define CGAL_NEF_POLYHEDRON_3_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Nef_3/SNC_items.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
//#include <CGAL/Nef_3/SNC_binop.h>
//#include <CGAL/Nef_3/SNC_walker.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_3/SNC_SM_explorer.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Modifier_base.h>

#ifdef CGAL_NEF3_CGAL_NEF3_SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // CGAL_NEF3_SM_VISUALIZOR

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/polyhedron_3_to_nef_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/assertions.h>

#include <list> // || (circulator_size(c) != 2 && !result));

// RO: includes for "vertex cycle to Nef" constructor
#include <CGAL/Nef_3/vertex_cycle_to_nef_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/normal_vector_newell_3.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 11
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename K, typename I, typename M> class Nef_polyhedron_3;
template <typename K, typename I, typename M> class Nef_polyhedron_3_rep;

template <typename K, typename I, typename M>
std::ostream& operator<<(std::ostream& os, Nef_polyhedron_3<K,I,M>& NP);

template <typename K, typename I, typename M>
std::istream& operator>>(std::istream& os, Nef_polyhedron_3<K,I,M>& NP);


template <typename K, typename I, typename M>
class Nef_polyhedron_3_rep 
{ 
  typedef Nef_polyhedron_3_rep<K,I,M>                  Self;
  friend class Nef_polyhedron_3<K,I,M>;
 public:
  typedef CGAL::SNC_structure<K,I,M>                   SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>           SNC_decorator;
  typedef CGAL::SNC_const_decorator<SNC_structure>     SNC_const_decorator;
  //  typedef CGAL::SNC_binop<SNC_structure>               SNC_binop;
  typedef CGAL::SNC_constructor<SNC_structure>         SNC_constructor;
  //typedef CGAL::SNC_walker<SNC_structure>              SNC_walker;
  typedef CGAL::SNC_point_locator<SNC_decorator> SNC_point_locator;
  typedef CGAL::SNC_simplify<SNC_structure>            SNC_simplify;
#ifdef CGAL_NEF3_POINT_LOCATOR_NAIVE
  typedef CGAL::SNC_point_locator_naive<SNC_decorator> SNC_point_locator_default;
#else
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_decorator> SNC_point_locator_default;
#endif

  typedef typename SNC_structure::Sphere_map       Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>           SM_decorator;
  typedef CGAL::SM_const_decorator<Sphere_map>     SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<SM_decorator>     SM_overlayer;
  typedef CGAL::SM_point_locator<SNC_structure>    SM_point_locator;
  typedef CGAL::SM_io_parser<SM_decorator>         SM_io_parser;

#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef CGAL::SNC_SM_visualizor<SNC_structure>       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR

 private:
  SNC_structure snc_;
  SNC_point_locator* pl_;
  
 public:
  Nef_polyhedron_3_rep() : snc_(), pl_() {}
  ~Nef_polyhedron_3_rep() { 
    CGAL_NEF_TRACEN( "Nef_polyhedron_3_rep: destroying SNC structure "<<&snc_<<
	    ", point locator "<<pl_);
    snc_.clear(); 
    delete pl_; 
  }
};

/*{\Manpage {Nef_polyhedron_3} {T} {Nef Polyhedra in Space}{N}}*/

/*{\Mdefinition
An instance of data type |\Mname| is a subset of 3-space which is the
result of forming complements and intersections starting from a set |H| of
halfspaces. |\Mtype| is closed under all binary set opertions |intersection|,
|union|, |difference|, |complement| and under the topological operations
|boundary|, |closure|, and |interior|.}*/

template <typename Kernel_, typename Items_ = SNC_items, typename Mark_ = bool>
class Nef_polyhedron_3 : public CGAL::Handle_for< Nef_polyhedron_3_rep<Kernel_, Items_, Mark_> >, 
			 public SNC_const_decorator<SNC_structure<Kernel_,Items_,Mark_> >
{ 
 public:
  /*{\Mtypes 7}*/  
  
  typedef Kernel_                                     Kernel;
  typedef Kernel_                                     Traits;
  typedef Items_                                      Items;
  typedef Mark_                                       Mark;
  typedef Nef_polyhedron_3<Kernel, Items, Mark>       Self;
  typedef Nef_polyhedron_3<Kernel, Items, Mark>       Nef_polyhedron;
  typedef Handle_for< Nef_polyhedron_3_rep<Kernel, Items, Mark> >   Base;
  typedef typename Kernel::Point_3                    Point_3;
  typedef typename Kernel::Plane_3                    Plane_3;
  typedef typename Kernel::Vector_3                   Vector_3;
  typedef typename Kernel::Segment_3                  Segment_3;
  typedef typename Kernel::Aff_transformation_3       Aff_transformation_3;
 
  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  typedef enum { EMPTY=0, COMPLETE=1 } Content;
  /*{\Menum construction selection}*/

  typedef enum { DEFAULT, NAIVE, WALKING, SPATIAL_SUBDIVISION  } Location_mode;
  /*{\Menum selection flag for the point location mode.}*/

protected: 
  struct AND { Mark operator()(const Mark& b1, const Mark& b2, bool inverted=false) const { return b1&&b2; } };
  struct OR { Mark operator()(const Mark& b1, const Mark& b2, bool inverted=false) const { return b1||b2; } };
  struct DIFF { Mark operator()(const Mark& b1, const Mark& b2, bool inverted=false) const { 
    if(inverted) return !b1&&b2; return b1&&!b2; } };
  struct XOR { Mark operator()(const Mark& b1, const Mark& b2, bool inverted=false) const 
    { return (b1&&!b2)||(!b1&&b2); } };

 public:
  typedef Nef_polyhedron_3_rep<Kernel,Items, Mark>    Nef_rep;
  typedef typename Nef_rep::SNC_structure       SNC_structure;
 protected:
  typedef typename Nef_rep::SNC_decorator       SNC_decorator;
  typedef typename Nef_rep::SNC_const_decorator SNC_const_decorator;
  typedef typename Nef_rep::SNC_constructor     SNC_constructor;
 //  typedef typename Nef_rep::SNC_binop           SNC_binop;
  //typedef typename Nef_rep::SNC_walker          SNC_walker;
  typedef typename Nef_rep::SNC_point_locator   SNC_point_locator;
  typedef typename Nef_rep::SNC_point_locator_default 
    SNC_point_locator_default;

  typedef typename Nef_rep::SM_decorator        SM_decorator;
  typedef typename Nef_rep::SM_const_decorator  SM_const_decorator;
  typedef typename Nef_rep::SM_overlayer        SM_overlayer;
  typedef typename Nef_rep::SM_point_locator    SM_point_locator;
  typedef typename Nef_rep::SM_io_parser        SM_io_parser;
  typedef typename Nef_rep::SNC_simplify        SNC_simplify;
#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef typename Nef_rep::SM_visualizor       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR

 typedef typename Nef_rep::Sphere_map                Sphere_map;
 public:
 typedef Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>  Nef_polyhedron_S2;
 protected:

  SNC_structure& snc() { return this->ptr()->snc_; } 
  const SNC_structure& snc() const { return this->ptr()->snc_; } 

  SNC_point_locator*& pl() { return this->ptr()->pl_; }
  const SNC_point_locator* pl() const { return this->ptr()->pl_; }

  friend std::ostream& operator<< <>
      (std::ostream& os, Nef_polyhedron_3<Kernel,Items, Mark>& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Nef_polyhedron_3<Kernel,Items, Mark>& NP);

  typedef typename SNC_decorator::Vertex_handle    Vertex_handle;
  typedef typename SNC_decorator::Halfedge_handle  Halfedge_handle;
  typedef typename SNC_decorator::Halffacet_handle
                                                   Halffacet_handle;
  typedef typename SNC_decorator::Volume_handle    Volume_handle;

 public:
  typedef typename SNC_structure::Sphere_point                 Sphere_point;
  typedef typename SNC_structure::Sphere_segment               Sphere_segment;
  typedef typename SNC_structure::Sphere_circle                Sphere_circle;
  typedef typename SNC_structure::Vertex_base                  Vertex;
  typedef typename SNC_structure::Halfedge_base                Halfedge;
  typedef typename SNC_structure::Halffacet_base               Halffacet;
  typedef typename SNC_structure::Volume_base                  Volume;
  typedef typename SNC_structure::Vertex_const_handle          Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle        Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle       Halffacet_const_handle;
  typedef typename SNC_structure::Volume_const_handle          Volume_const_handle;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_svertex_const_circulator 
                                  SHalfedge_around_svertex_const_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
                                  SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator 
                                  SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_const_circulator 
                                  SHalfedge_around_sface_const_circulator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator     
                                  Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Infi_box                     Infi_box;
  typedef typename SNC_structure::Size_type Size_type;
  typedef Size_type                         size_type;

  typedef typename Kernel::RT                       RT;

 public:
  typedef typename SM_decorator::SVertex_handle    SVertex_handle;
  typedef typename SM_decorator::SHalfedge_handle  SHalfedge_handle;
  typedef typename SM_decorator::SFace_handle      SFace_handle;
  typedef typename SM_decorator::SVertex_const_handle
                                                   SVertex_const_handle;
  typedef typename SM_decorator::SHalfedge_const_handle
                                                   SHalfedge_const_handle;
  typedef typename SM_decorator::SHalfloop_const_handle
                                                   SHalfloop_const_handle;
  typedef typename SM_decorator::SFace_const_handle
                                                   SFace_const_handle;
  typedef typename SNC_decorator::Vertex_iterator  Vertex_iterator;
  typedef typename SNC_decorator::Halfedge_iterator
                                                   Halfedge_iterator;
  typedef typename SNC_decorator::Halffacet_iterator
                                                   Halffacet_iterator;
  typedef typename SNC_structure::Shell_entry_iterator
                                                   Shell_entry_iterator;
  typedef typename SNC_decorator::Volume_iterator  Volume_iterator;
  typedef typename SNC_structure::Vertex_const_iterator
                                                    Vertex_const_iterator;
  typedef typename SNC_structure::Halfedge_const_iterator 
                                                   Halfedge_const_iterator;
  typedef typename SNC_structure::Halffacet_const_iterator     
                                                   Halffacet_const_iterator;
  typedef typename SNC_structure::Volume_const_iterator     
                                                   Volume_const_iterator;
  typedef typename SNC_structure::Shell_entry_const_iterator
                                                   Shell_entry_const_iterator;
  typedef typename SM_decorator::SVertex_iterator  SVertex_iterator;
  typedef typename SM_decorator::SHalfedge_iterator
                                                   SHalfedge_iterator;
  typedef typename SM_decorator::SFace_iterator    SFace_iterator;
  typedef typename SM_decorator::SVertex_const_iterator
                                                   SVertex_const_iterator;
  typedef typename SM_decorator::SHalfedge_const_iterator 
                                                   SHalfedge_const_iterator;
  typedef typename SM_decorator::SHalfloop_const_iterator 
                                                   SHalfloop_const_iterator;
  typedef typename SM_decorator::SFace_const_iterator     
                                                   SFace_const_iterator;
  typedef typename SNC_decorator::SFace_cycle_const_iterator     
                                                   SFace_cycle_const_iterator;

 protected: 
  void initialize_infibox_vertices(Content space) {
    SNC_constructor C(snc()); 
    Infi_box::initialize_infibox_vertices(C, space == COMPLETE);
  }

  void check_h_for_intersection_of_12_cube_edges_and_add_vertices
  (const Plane_3& p);
  void create_intersection_vertex_of_h_and_e();
  void init_cube_vertices_depending_on_h(const Plane_3& p);
  void add_h_to_local_view_of_v();
  
 public:
  void build_external_structure() {
    SNC_constructor C(snc(), pl());
    C.build_external_structure();
  }

 public:
  /*{\Mcreation 3}*/

  Nef_polyhedron_3( Content space = EMPTY,
		    SNC_point_locator* _pl = new SNC_point_locator_default);
		   
  /*{\Mcreate creates an instance |\Mvar| of type |\Mname|
  and initializes it to the empty set if |space == EMPTY|
  and to the whole space if |space == COMPLETE|.}*/

  Nef_polyhedron_3(const Plane_3& p, 
		   Boundary b = INCLUDED,
		   SNC_point_locator* _pl = new SNC_point_locator_default);
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the
  halfspace on the positive side of |p| including |p| if |b==INCLUDED|,
  excluding |p| if |b==EXCLUDED|.}*/

  Nef_polyhedron_3(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) : Base(N1) {
    set_snc(snc());
  } 

  Nef_polyhedron_3& operator=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) { 
    Base::operator=(N1);
    set_snc(snc());
    return (*this); 
  }

  ~Nef_polyhedron_3() { 
    CGAL_NEF_TRACEN("~Nef_polyhedron_3: destructor called for snc "<<&snc()<<
	   ", pl "<<pl());
  }

   // RO: "vertex cycle to Nef" constructor (main part)
   // II input iterator; KN kernel of normal (may differ from Nef kernel)
   template <class II, class KN>
   Nef_polyhedron_3 (II v_first, II v_last,
                    const CGAL::Vector_3<KN> &normal, bool verb = false)
   {  CGAL_NEF_TRACEN("construction from vertex cycle (main part)");

      // project and triangulate vertices,
      // convert result to Nef_polyhedron
      CGAL_precondition (!CGAL::is_empty_range (v_first, v_last));
      bool is_nef = vertex_cycle_to_nef_3<Nef_polyhedron> (snc(),
         v_first, v_last, normal, verb);
      if (is_nef)
      {
	 // TO DO:
	 // Wie kann der eigene point_locator pl() eingebunden werden?
	 // Wie kann der Konstruktor umgangen werden?
         typedef CGAL::SNC_point_locator_by_spatial_subdivision
            <CGAL::SNC_decorator<SNC_structure> >    Point_locator;

         Point_locator Pl;
         CGAL::SNC_constructor<SNC_structure> con (snc(), &Pl);
         con.build_external_structure();
         *this = Nef_polyhedron(snc(), &Pl);
      }
      else
      {  *this = Nef_polyhedron();
      }
      set_snc (snc());
      CGAL_expensive_postcondition (is_valid());
   }

   // RO: "vertex cycle to Nef" constructor (normal computation)
   template <class II>
   Nef_polyhedron_3 (II v_first, II v_last, bool verb = false)
   {  CGAL_NEF_TRACEN("construction from vertex cycle (normal computation)");

      // compute normal vector
      CGAL_precondition (!CGAL::is_empty_range (v_first, v_last));
      CGAL::Vector_3<typename II::value_type::R> normal;
      normal_vector_newell_3 (v_first, v_last, normal);

      // call "main" constructor
      *this = Nef_polyhedron_3 (v_first, v_last, normal, verb);
      set_snc (snc());
   }

  typedef Polyhedron_3< Kernel> Polyhedron;
  
 template <class T1, class T2,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template <class T31, class T32, class T33>
#endif
           class T3, class T4 >
 Nef_polyhedron_3( CGAL::Polyhedron_3<T1,T2,T3,T4>& P,
		   SNC_point_locator* _pl = new SNC_point_locator_default) {
    CGAL_NEF_TRACEN("construction from Polyhedron_3");
    SNC_structure rsnc;
    *this = Nef_polyhedron_3(rsnc, _pl, false);
    initialize_infibox_vertices(EMPTY);
    polyhedron_3_to_nef_3
      <CGAL::Polyhedron_3<T1,T2,T3,T4>, SNC_structure>( P, snc());
    build_external_structure();
    simplify();
    set_snc(snc());
    CGAL_assertion(orientation() == 1);
  }
  
 protected:  
  template <class HDS>
  class Build_polyhedron : public CGAL::Modifier_base<HDS> {
    
    class Visitor {

      const Object_index<Vertex_const_iterator>& VI;
      Polyhedron_incremental_builder_3<HDS>& B;
      SNC_const_decorator& D;
      
    public:
      Visitor(Polyhedron_incremental_builder_3<HDS>& BB,
	      SNC_const_decorator& sd,
	      Object_index<Vertex_const_iterator>& vi) : VI(vi), B(BB), D(sd){}

      void visit(Halffacet_const_handle opposite_facet) {

	CGAL_NEF_TRACEN("Build_polyhedron: visit facet " << D.plane(opposite_facet));
 
	CGAL_assertion(Infi_box::is_standard(D.plane(opposite_facet)));
	
	SHalfedge_const_handle se;
	Halffacet_cycle_const_iterator fc;
     	
	Halffacet_const_handle f = opposite_facet->twin();

	B.begin_facet();
	fc = f->facet_cycles_begin();
	se = SHalfedge_const_handle(fc);
	CGAL_assertion(se!=0);
	SHalfedge_around_facet_const_circulator hc_start(se);
	SHalfedge_around_facet_const_circulator hc_end(hc_start);
	CGAL_For_all(hc_start,hc_end) {
	  CGAL_NEF_TRACEN("   add vertex " << hc_start->source()->center_vertex()->point());
	  B.add_vertex_to_facet(VI[hc_start->source()->center_vertex()]);
	}
	B.end_facet();
      }

      void visit(SFace_const_handle s) {}
      void visit(Halfedge_const_handle e) {}
      void visit(Vertex_const_handle v) {}
      void visit(SHalfedge_const_handle se) {}
      void visit(SHalfloop_const_handle sl) {}
    };

  public:

    SNC_const_decorator& scd;
    Object_index<Vertex_const_iterator> VI;

    Build_polyhedron(SNC_const_decorator& s) : 
      scd(s), VI(s.vertices_begin(),s.vertices_end(),'V') {}
    
    void operator()( HDS& hds) {

      Polyhedron_incremental_builder_3<HDS> B(hds, true);

      int skip_volumes;
      if(Infi_box::extended_kernel()) {
	B.begin_surface(scd.number_of_vertices()-8, 
			scd.number_of_facets()-6,
			scd.number_of_edges()-12);
	skip_volumes = 2;
      }
      else {
	B.begin_surface(scd.number_of_vertices(), 
			scd.number_of_facets(),
			scd.number_of_edges());
	skip_volumes = 1;
      }
      
      int vertex_index = 0;
      Vertex_const_iterator v;
      CGAL_forall_vertices(v,scd) {
	if(Infi_box::is_standard(v->point())) {
	  VI[v]=vertex_index++;
	  B.add_vertex(v->point());
	}
      }     
      
      Visitor V(B,scd,VI);
      Volume_const_handle c;
      CGAL_forall_volumes(c,scd)
	if(skip_volumes-- <= 0)
	  scd.visit_shell_objects(SFace_const_handle(c->shells_begin()),V);
      B.end_surface();
    }

  };

 public:
 void delegate( Modifier_base<SNC_structure>& modifier) {
   // calls the `operator()' of the `modifier'. Precondition: The
   // `modifier' returns a consistent representation.
   modifier(snc());
   // build_external_structure(); // TO DO: conflict with CGAL::Mark_bounded_volumes
   simplify();
   CGAL_expensive_postcondition( is_valid());
 }

 struct SNC_and_PL {
   SNC_structure* sncp;
   SNC_point_locator* pl;

   SNC_and_PL(SNC_structure* s, SNC_point_locator* p) : sncp(s), pl(p) {}
 };

 void delegate( Modifier_base<SNC_and_PL>& modifier) {
   // calls the `operator()' of the `modifier'. Precondition: The
   // `modifier' returns a consistent representation.
   SNC_and_PL sncpl(&snc(),pl());
   modifier(sncpl);
   pl() = sncpl.pl;
   CGAL_expensive_postcondition( is_valid());
 }
 
 protected:
   int orientation() {
    // Function does not work correctly with every Nef_polyhedron. 
    // It works correctly if v_max is 2-manifold
    CGAL_NEF_TRACEN("orientation");
    
    typedef typename SM_decorator::SHalfedge_around_sface_circulator
      SHalfedge_around_sface_circulator;

    if ( this->number_of_vertices() == 0)
      return 0;
    
    // Find top-most vertex v_max (top-most also works fine for terrains)
    SNC_decorator D(snc());
    bool first = true;
    Vertex_handle v_max;
    Vertex_handle vh;
    CGAL_forall_vertices(vh, D) {
      CGAL_NEF_TRACEN("  is_standard " << vh->point() << " " << Infi_box::is_standard(D.point(vh)));
      if (Infi_box::is_standard(D.point(vh)) && 
	  (first || D.point(vh).z() > D.point(v_max).z())) {
	first = false;
	v_max = vh;
      }
    }
    CGAL_NEF_TRACEN("  vmax " << v_max->point());
    CGAL_assertion_msg(!first, "off file is empty");
    SM_decorator SD(&*v_max);

    // Add up z-coord. of all normalized normal vectors of the incident facets
    
    double z=1.0;
    first = true;
    SFace_handle sf;
    CGAL_forall_sfaces(sf, SD) {
      if(Volume_const_handle(sf->volume()) != Infi_box::getNirvana(snc())) continue;
      CGAL_assertion_msg(first, "function shall not be used on this polyhedron");
      SHalfedge_around_sface_circulator estart(sf->sface_cycles_begin()), 
	                                eend(estart);
      CGAL_For_all(estart,eend) {
	Sphere_circle c(SD.circle(estart));
	CGAL_NEF_TRACEN("  circle " << c); 
	double delta_z = CGAL::to_double(c.orthogonal_vector().hz());
	CGAL_NEF_TRACEN("  delta_z " << delta_z);
	Point_3 target = ORIGIN + c.orthogonal_vector();
	Segment_3 s(Point_3(0,0,0), target);
	delta_z /= CGAL::sqrt(CGAL::to_double(s.squared_length()));
	z += delta_z;	
	CGAL_NEF_TRACEN("  z " << z);
      }
      first = false;
    }
    return sign(z);
  }

 public:

  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  void convert_to_Polyhedron(Polyhedron& P) {
    CGAL_precondition(is_simple());
    Build_polyhedron<HalfedgeDS> bp(*this);    
    P.delegate(bp);
  }

 //  void dump(bool sorted = false, std::ostream& os = std::cout) 
 //  { SNC_io_parser::dump( snc(), os, sorted); }

  bool is_valid( bool verb = false, int level = 0) {
    // checks the combinatorial consistency.
    Verbose_ostream verr(verb);
    verr << "begin CGAL::Nef_polyhedron_3<...>::is_valid( verb=true, "
      "level = " << level << "):" << std::endl;

    SNC_decorator D(snc());
    bool valid = D.is_valid(verb, level);
    verr << "end of CGAL::Nef_polyhedron_3<...>::is_valid(): structure is "
	 << ( valid ? "valid." : "NOT VALID.") << std::endl;
    return valid;
  }

  bool is_simple() {

    Halfedge_iterator e;
    CGAL_forall_edges(e,snc())
      if(!is_edge_2manifold(e))
	return false;

    Vertex_iterator v;
    CGAL_forall_vertices(v,snc())
      if(!is_vertex_2manifold(v))
	return false;

    Halffacet_iterator f;
    CGAL_forall_halffacets(f,snc())
      if(!is_facet_simple(f))
	return false;

    return true;
  }
 
 bool is_convex() const {
   
   Vertex_const_iterator v;
   CGAL_forall_vertices(v, *this) {

     SM_const_decorator SD(&*v);
     if(std::distance(SD.sfaces_begin(),SD.sfaces_end())!=2)
       return false;

     if(!Infi_box::is_standard(v->point())) continue;

     SFace_const_iterator sf;
     CGAL_forall_sfaces(sf,SD) {
       if(sf->volume() == Infi_box::getNirvana(snc())) continue;
       if(std::distance(sf->sface_cycles_begin(),sf->sface_cycles_end())!=1)
	 return false;
       SFace_cycle_const_iterator sfi(sf->sface_cycles_begin());
       if(!sfi.is_shalfedge())
	 return false;
       SHalfedge_const_handle se(sf->sface_cycles_begin());
       SHalfedge_around_sface_const_circulator sec(se),send(sec);
       CGAL_For_all(sec,send)
	 if(spherical_orientation(sec->source()->point(),
				  sec->snext()->source()->point(),
				  sec->snext()->snext()->source()->point())<0) {
	   std::cerr << "vertex at " << v->point() << " is not convex" << std::endl;
	   return false;
	 }
     }
   }
   return true;
 }

 private:  
  bool is_edge_2manifold(const Halfedge_handle& e) {

    SM_decorator SD;
    SHalfedge_around_svertex_circulator c(SD.first_out_edge(e)), c2(c);

    if(c == 0) {
      CGAL_assertion(circulator_size(c) !=2);
      return false;
    }

    if(++c == c2){
      CGAL_assertion(circulator_size(c) !=2);
      return false;
    }

    if(++c != c2) { 
      CGAL_assertion(circulator_size(c) !=2);
      return false;
    }
    
    CGAL_assertion(circulator_size(c) == 2);
    return true;
  }
 
  bool is_vertex_2manifold(const Vertex_handle& v) {
     
    SFace_iterator sfi(v->sfaces_begin());
    if (++sfi != v->sfaces_last())
      return false;

    return true;
  }

  bool is_facet_simple(const Halffacet_const_handle& f) {
    
    bool found_first = false;
    Halffacet_cycle_const_iterator it; 
    CGAL_forall_facet_cycles_of(it,f)
      if (found_first || !it.is_shalfedge())
	return false;
      else
	found_first = true;
   
    return true;
  }

 public:
   
  void clear(Content space = EMPTY)
    { *this = Nef_polyhedron_3(space, pl()->clone()); }
  /*{\Mop makes |\Mvar| the empty set if |space == EMPTY| and the
  full space if |space == COMPLETE|.}*/

 bool is_empty() const {
   /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
   if(Infi_box::extended_kernel()) 
     return this->number_of_vertices() == 8 &&
            this->number_of_edges() == 12 &&
            this->number_of_facets() == 6 &&
            this->number_of_volumes() == 2 &&
            mark(++this->volumes_begin()) == false;

   else 
     return this->number_of_vertices() == 0 &&
            this->number_of_edges() == 0 &&
            this->number_of_facets() == 0 &&
            this->number_of_volumes() == 1 &&
            mark(this->volumes_begin()) == false;
  }

 bool is_space() const {
  /*{\Mop returns true if |\Mvar| is the whole space, false otherwise.}*/
   if(Infi_box::extended_kernel()) 
     return this->number_of_vertices() == 8 &&
            this->number_of_edges() == 12 &&
            this->number_of_facets() == 6 &&
            this->number_of_volumes() == 2 &&
            mark(++this->volumes_begin()) == true;

   else 
     return this->number_of_vertices() == 0 &&
            this->number_of_edges() == 0 &&
            this->number_of_facets() == 0 &&
            this->number_of_volumes() == 1 &&
            mark(this->volumes_begin()) == true;
  }

  /*{\Xtext \headerline{Destructive Operations}}*/

 protected:
  void clone_rep() { *this = Nef_polyhedron_3<Kernel,Items, Mark>(snc(), pl()); }
  void empty_rep() { 
    SNC_structure rsnc;
    *this = Nef_polyhedron_3<Kernel,Items, Mark>(rsnc, new SNC_point_locator_default);
  }

 public:
  Nef_polyhedron_3( const SNC_structure& W, 
		    SNC_point_locator* _pl = new SNC_point_locator_default,
		    bool clone_pl = true,
		    bool clone_snc = true);
  /*{\Xcreate makes |\Mvar| a new object.  If |cloneit==true| then the
  underlying structure of |W| is copied into |\Mvar|.}*/
  // TODO: granados: define behavior when clone=false

  /*{\Moperations 4 3 }*/

  void simplify() {
    SNC_simplify simp(snc());
    bool simplified = simp.simplify();
    CGAL_NEF_TRACEN( "simplify(): structure simplified? "<<simplified);
    
    if( simplified) {
#ifdef CGAL_NEF3_UPDATE_K3TREE_AFTER_SIMPLIFICATION
      /*debug*/ snc().print_statistics();
      Unique_hash_map<Vertex_handle, bool> 
	V(false, snc().number_of_vertices());
      Unique_hash_map<Halfedge_handle, bool> 
	E(false, snc().number_of_halfedges());
      Unique_hash_map<Halffacet_handle, bool> 
	F(false, snc().number_of_halffacets());
      Vertex_iterator v;
      Halfedge_iterator e;
      Halffacet_iterator f;
      CGAL_forall_vertices( v, snc()) V[Vertex_handle(v)] = true;
      CGAL_forall_halfedges( e, snc()) E[Halfedge_handle(e)] = true;
      CGAL_forall_halffacets( f, snc()) F[Halffacet_handle(f)] = true;
      bool updated = pl()->update( V, E, F);
      CGAL_NEF_TRACEN("simplify(): point locator structure updated? " << updated);
#else
      SNC_point_locator* old_pl = pl();
      pl() = pl()->clone();
      pl()->initialize(&snc());
      delete old_pl;
#endif
    }
  }

 public:
  Nef_polyhedron_S2 get_sphere_map(Vertex_const_handle v) const {
    return Nef_polyhedron_S2(*v);
  }

  void extract_complement();
  /*{\Xop converts |\Mvar| to its complement. }*/
  void extract_interior();
  /*{\Xop converts |\Mvar| to its interior. }*/
  void extract_boundary();
  /*{\Xop converts |\Mvar| to its boundary. }*/

  void extract_closure()
  /*{\Xop converts |\Mvar| to its closure. }*/
  { CGAL_NEF_TRACEN("extract closure");
    if( this->is_shared()) clone_rep();
    extract_complement();
    extract_interior();
    extract_complement();
  }

  void extract_regularization()
  /*{\Xop converts |\Mvar| to its regularization. }*/
  { CGAL_NEF_TRACEN("extract regularization");
    if( this->is_shared()) clone_rep();
    extract_interior();
    extract_closure();
  }

  /*{\Mtext \headerline{Constructive Operations}}*/

  Nef_polyhedron_3<Kernel,Items, Mark> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane. }*/
  { Nef_polyhedron_3<Kernel,Items, Mark> res = *this;
    res.extract_complement();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> interior() const
  /*{\Mop    returns the interior of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items, Mark> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> closure() const
  /*{\Mop returns the closure of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items, Mark> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> boundary() const
  /*{\Mop returns the boundary of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items, Mark> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> regularization() const
  /*{\Mop    returns the regularized polyhedron (closure of 
    the interior).}*/
  { Nef_polyhedron_3<Kernel,Items, Mark> res = *this;
    res.extract_regularization();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark>
  intersection(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
    /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/ {
    CGAL_NEF_TRACEN(" intersection between nef3 "<<&*this<<" and "<<&N1);
    AND _and;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, AND> tests_impl;
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    SNC_decorator D( res.snc());
    D.binary_operation(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _and);
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> 
  join(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/ { 
    CGAL_NEF_TRACEN(" join between nef3 "<<&*this<<" and "<<&N1);
    OR _or;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, OR> tests_impl;
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    SNC_decorator D( res.snc());
    D.binary_operation(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _or);
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark> 
  difference(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/ { 
    CGAL_NEF_TRACEN(" difference between nef3 "<<&*this<<" and "<<&N1);
    DIFF _diff;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, DIFF> tests_impl;
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    SNC_decorator D( res.snc());
    D.binary_operation(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _diff);
    return res;
  }    

  Nef_polyhedron_3<Kernel,Items, Mark> 
  symmetric_difference(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/ {
    CGAL_NEF_TRACEN(" symmetic difference between nef3 "<<&*this<<" and "<<&N1);
    XOR _xor;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, XOR> tests_impl;
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    SNC_decorator D( res.snc());
    D.binary_operation(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _xor);
    return res;
  }


  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement}. There are also the corresponding
  modification operations |*=,+=,-=,^=|.}*/

  Nef_polyhedron_3<Kernel,Items, Mark>  operator*(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const 
  { return intersection(N1); }

  Nef_polyhedron_3<Kernel,Items, Mark>  operator+(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return join(N1); }

  Nef_polyhedron_3<Kernel,Items, Mark>  operator-(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return difference(N1); }

  Nef_polyhedron_3<Kernel,Items, Mark>  operator^(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return symmetric_difference(N1); }

  Nef_polyhedron_3<Kernel,Items, Mark>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_3<Kernel,Items, Mark>& operator*=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1)
  { *this = intersection(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items, Mark>& operator+=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1)
  { *this = join(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items, Mark>& operator-=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1)
  { *this = difference(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items, Mark>& operator^=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1)
  { *this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality.}*/

  bool operator==(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { CGAL_NEF_TRACEN(" equality comparision between nef3 "<<&*this<<" and "<<&N1);
    return symmetric_difference(N1).is_empty(); }

  bool operator!=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { CGAL_NEF_TRACEN(" inequality comparision between nef3 "<<&*this<<" and "<<&N1);
    return !operator==(N1); }  

  bool operator<(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return !N1.difference(*this).is_empty() && difference(N1).is_empty(); } 

  bool operator>(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return difference(*this).is_empty() && !difference(N1).is_empty(); } 

  bool operator<=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return difference(N1).is_empty(); } 

  bool operator>=(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) const
  { return N1.difference(*this).is_empty(); } 
 /*
  void transform( const Aff_transformation_3& aff) {
    
    // precondition: the polyhedron as a bounded boundary
    // (needs to be explicitly tested at some time)

    if( is_shared())
      clone_rep();
    // only linear transform for the origin-centered sphere maps
    Aff_transformation_3 linear( aff.hm(0,0), aff.hm(0,1), aff.hm(0,2),
				 aff.hm(1,0), aff.hm(1,1), aff.hm(1,2),
				 aff.hm(2,0), aff.hm(2,1), aff.hm(2,2),
				 aff.hm(3,3));
    
    SNC_decorator deco( snc());
    
    list<Vertex_handle> vertex_list;
    Vertex_iterator vi;
    CGAL_forall_vertices( vi, snc()) {
      
      CGAL_NEF_TRACEN("transform vertex ");
      if (!is_standard(vi)) {
	vertex_list.push_bach(vi);
      } else {
	vi->point() = vi->point().transform( aff);
	SM_decorator sdeco(&*vi);
	sdeco.transform( linear);
      }
    }

    CGAL_forall_halffacets(fi, snc()) {
      if (!deco.is_standard( fi)) {
	Halffacet_cycle_iterator fc = fi.facet_cycles_begin();
	CGAL_assertion(fc.is_shalfedge());
	SHalfedge_around_facet_circulator fcc(SHalfedge(fc)), fend(fcc);
	CGAL_For_all(fcc,fend) {
	  deco.add_interim_points(point(source(fcc)),point(target(fcc)));
	}
      }
    }

    typename list<Vertex_handle>::iterator li;
    for(li = vertex_list.begin(); li != vertex_list.end(); li++){
      deco.transform_sphere(*li);
    }
        
    if(!is_finite()) {
      build_external_structure();
    } else {
      Halffacet_iterator fi;
      CGAL_forall_halffacets(fi,snc()) {
	fi->plane() = fi->plane().transform( aff);
      }    
    }

    if(aff.homogeneous(0,1) != 0 ||
       aff.homogeneous(0,2) != 0 ||
       aff.homogeneous(1,0) != 0 ||
       aff.homogeneous(1,2) != 0 ||
       aff.homogeneous(2,0) != 0 ||
       aff.homogeneous(2,1) != 0 ||
       aff.homogeneous(0,0) != aff.homogeneous(1,1) ||
       aff.homogeneous(0,0) != aff.homogeneous(2,2)) {

      SNC_point_locator* old_pl = pl();
      pl() = pl()->clone();
      pl()->initialize(&snc());
      delete old_pl;   
    }
    else
      pl()->transform(aff);

  }
*/

 bool is_90degree_rotation(const Aff_transformation_3& aff) const {
   if(aff.hm(0,3) != 0) return false;
   if(aff.hm(1,3) != 0) return false;
   if(aff.hm(2,3) != 0) return false;
   if(CGAL_NTS abs(aff.hm(0,0)) + 
      CGAL_NTS abs(aff.hm(0,1)) + 
      CGAL_NTS abs(aff.hm(0,2)) != aff.hm(3,3)) return false;
   if(CGAL_NTS abs(aff.hm(1,0)) + 
      CGAL_NTS abs(aff.hm(1,1)) + 
      CGAL_NTS abs(aff.hm(1,2)) != aff.hm(3,3)) return false;
   if(CGAL_NTS abs(aff.hm(2,0)) + 
      CGAL_NTS abs(aff.hm(2,1)) + 
      CGAL_NTS abs(aff.hm(2,2)) != aff.hm(3,3)) return false;
   if(CGAL_NTS abs(aff.hm(0,0)) + 
      CGAL_NTS abs(aff.hm(1,0)) + 
      CGAL_NTS abs(aff.hm(2,0)) != aff.hm(3,3)) return false;
   if(CGAL_NTS abs(aff.hm(0,1)) + 
      CGAL_NTS abs(aff.hm(1,1)) + 
      CGAL_NTS abs(aff.hm(2,1)) != aff.hm(3,3)) return false;
   if(CGAL_NTS abs(aff.hm(0,2)) + 
      CGAL_NTS abs(aff.hm(1,2)) + 
      CGAL_NTS abs(aff.hm(2,2)) != aff.hm(3,3)) return false;
   return true;
 }

 bool is_scaling(const Aff_transformation_3& aff) const {
   if(aff.hm(0,3) != 0) return false;
   if(aff.hm(1,3) != 0) return false;
   if(aff.hm(2,3) != 0) return false;   
   if(aff.hm(0,1) != 0) return false;
   if(aff.hm(0,2) != 0) return false;
   if(aff.hm(1,0) != 0) return false;   
   if(aff.hm(1,2) != 0) return false;
   if(aff.hm(2,0) != 0) return false;
   if(aff.hm(2,1) != 0) return false;
   if(aff.hm(0,0) != aff.hm(1,1)) return false;
   if(aff.hm(0,0) != aff.hm(2,2)) return false;
   return true;
 }

  void transform( const Aff_transformation_3& aff) {
    
    CGAL_precondition(aff.is_even());

    // precondition: the polyhedron as a bounded boundary
    // (needs to be explicitly tested at some time)

    if( this->is_shared())
      clone_rep();
    // only linear transform for the origin-centered sphere maps
    Aff_transformation_3 linear( aff.hm(0,0), aff.hm(0,1), aff.hm(0,2),
				 aff.hm(1,0), aff.hm(1,1), aff.hm(1,2),
				 aff.hm(2,0), aff.hm(2,1), aff.hm(2,2),
				 aff.hm(3,3));
    
    SNC_constructor cstr(snc());
    
    std::list<Vertex_handle> vertex_list;
    std::list<Vertex_handle> corner_list;
    std::list<Vertex_handle> delete_list;
    typename std::list<Vertex_handle>::iterator li;
    typename std::list<Vertex_handle>::iterator li2;

    bool ninety = is_90degree_rotation(aff);
    bool scale = is_scaling(aff);

    Vertex_iterator vi;
    CGAL_forall_vertices( vi, snc()) {
      
      CGAL_NEF_TRACEN("transform vertex ");
      if(scale) {
	if(is_standard(vi))
	  vi->point() = vi->point().transform( aff);
	else if(!Infi_box::is_infibox_corner(vi->point())) {
	  vi->point() = normalized(Infi_box::normalize_transformed_vertex(vi->point().transform(aff)));
	}
      } else if (!is_standard(vi) && !ninety) {
	if(Infi_box::is_infibox_corner(vi->point()))
	  corner_list.push_back(vi);
	vertex_list.push_back(vi);
      } else {
	vi->point() = vi->point().transform( aff);
	SM_decorator sdeco(&*vi);
	sdeco.transform( linear);
      }
    }

    if(!this->is_bounded() && !ninety && !scale) {
      Halffacet_iterator fi;
      CGAL_forall_facets(fi, snc()) {
	if(!is_standard(fi) || is_bounded(fi)) continue;
	Plane_3 pt = fi->plane();
	pt = pt.transform(aff);
	std::list<Point_3> points(Infi_box::find_points_of_box_with_plane(cstr,pt));
	std::list<Vertex_handle> newVertices;
	newVertices = Infi_box::create_vertices_on_infibox(cstr,
							   pt, points, fi->mark(), 
							   fi->twin()->incident_volume()->mark(), 
							   fi->incident_volume()->mark());

	for(li = newVertices.begin(); li != newVertices.end(); ++li) {
	  if(Infi_box::is_infibox_corner(point(*li))) {
	    li2 = corner_list.begin();
	    while(li2 != corner_list.end() && point(*li2) != point(*li)) ++li2;
	    CGAL_assertion(li2 != corner_list.end());
	    delete_list.push_back(*li2);
	    *li2 = *li;
	  }
	}
      }
      
      for(li = vertex_list.begin(); li != vertex_list.end();++li) {
	SM_decorator SD(&**li);
	if(Infi_box::is_complex_facet_infibox_intersection(**li)) {
	  Halffacet_handle hf[2];
	  int i=0;
	  SHalfedge_iterator sei;
	  CGAL_forall_sedges(sei,SD) {
	    if(!Infi_box::is_sedge_on_infibox(sei)) {
	      hf[i] = sei->facet();
	      if(hf[i]->is_twin()) hf[i] = hf[i]->twin();
	      ++i;
	    }
	    if(i>1)
	      break;
	  }
	}
      }

      cstr.clear_external_structure();
      for(li = vertex_list.begin(); li != vertex_list.end();++li){
	if(Infi_box::is_complex_facet_infibox_intersection(**li)) {
	  Vertex_handle v2;
	  Vertex_handle v1 = cstr.create_for_infibox_overlay(*li);
	  v1->point() = normalized(Infi_box::normalize_transformed_vertex(point(*li).transform(aff)));
	  SM_decorator sdeco(&*v1);
	  sdeco.transform(linear);	    
	  switch(Infi_box::type_of_infibox_point(v1->point())) {
	  case 1: 
	    v2 = cstr.create_from_point_on_infibox_facet(v1->point()); 
	    break;
	  case 2: 
	    v2 = cstr.create_from_point_on_infibox_edge(v1->point());
	    break;
	  case 3: 
	    v2 = cstr.create_from_point_on_infibox_vertex(v1->point());
	    li2 = corner_list.begin();
	    while(li2 != corner_list.end() && point(*li2) != v2->point()) ++li2;
	    if(li2 != corner_list.end())
	      delete_list.push_back(*li2);
	    break;
	  default: CGAL_assertion_msg(false, "wrong value");
	  }
	  Vertex_handle v = snc().new_vertex(v1->point(), mark(*li));
	  SM_overlayer O(&*v);
	  O.subdivide(&*v1,&*v2);
	  AND _and;
	  O.select(_and);
	  O.simplify();
	  snc().delete_vertex(v1);
	  snc().delete_vertex(v2);
	}
	
	if(Infi_box::is_infibox_corner(point(*li))) {
	  SM_decorator SD(&**li);
	  if(SD.number_of_svertices() < 4)
	    continue;
	  li2 = corner_list.begin();
	  while(li2 != corner_list.end() && point(*li2) != point(*li)) ++li2;
	  CGAL_assertion(li2 != corner_list.end());
	  if(*li == *li2) {
	    delete_list.push_back(*li2);
	    *li2 = cstr.create_from_point_on_infibox_vertex(point(*li));
	  }
	} else 
	  snc().delete_vertex(*li);	  
      }

      for(li = delete_list.begin(); li != delete_list.end(); ++li)
	snc().delete_vertex(*li);

      while(cstr.erase_redundant_vertices());
      cstr.correct_infibox_sedge_marks();

      build_external_structure();
      cstr.correct_infibox_sface_marks();

      SNC_point_locator* old_pl = pl();
      pl() = pl()->clone();
      pl()->initialize(&snc());
      delete old_pl;   

    } else {
      Halffacet_iterator fi;
      CGAL_forall_halffacets(fi,snc()) {
	if(is_standard(fi) || ninety)
	  fi->plane() = fi->plane().transform( aff);
      }    

      if(aff.homogeneous(0,1) != 0 ||
	 aff.homogeneous(0,2) != 0 ||
	 aff.homogeneous(1,0) != 0 ||
	 aff.homogeneous(1,2) != 0 ||
	 aff.homogeneous(2,0) != 0 ||
	 aff.homogeneous(2,1) != 0 ||
	 aff.homogeneous(0,0) != aff.homogeneous(1,1) ||
	 aff.homogeneous(0,0) != aff.homogeneous(2,2) ||
	 !this->is_bounded()) {
	   SNC_point_locator* old_pl = pl();
	   pl() = pl()->clone();
	   pl()->initialize(&snc());
	   delete old_pl;   
	 }
      else pl()->transform(aff); 
    }
  }
  
  /*{\Mtext \headerline{Exploration}
  As Nef polyhedra are the result of forming complements 
  and intersections starting from a set |H| of halfspaces which are
  defined by oriented planes in three space. The corresponding 
  structure is represented by an extended wuerzburg structure 
  $W = (V,E,F,C)$. For topological queries within |W| the following 
  types and operations allow exploration access to this structure.}*/

  /*{\Mtypes 3}*/
    
    typedef CGAL::SNC_SM_explorer<SM_const_decorator>  SM_explorer;

    SM_explorer SMexplorer(Vertex_const_handle v) const { 
      SM_const_decorator SMCD(&*v);
      return SM_explorer(SMCD); 
    }

  typedef typename SNC_structure::Object_list Object_list;
  typedef typename SNC_structure::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying
  plane map. The kind of object |(vertex, halfedge, face)| can 
  be determined and the object can be assigned to a corresponding
  handle by the three functions:\\
  |bool assign(Vertex_const_handle& h, Object_handle)|\\
  |bool assign(Edge_const_handle& h, Object_handle)|\\
  |bool assign(Facet_const_handle& h, Object_handle)|\\
  |bool assign(Volume_const_handle& h, Object_handle)|\\
  where each function returns |true| iff the assignment to
  |h| was done.}*/

  /*{\Moperations 3 1 }*/

  bool contains(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the set
  represented by |\Mvar|.}*/
    // { SNC_point_locator PL(snc()); return PL.mark(h);} 
    { CGAL_assertion_msg( 0, "not implemented."); return false;}

  bool contained_in_boundary(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the $2$-skeleton
  of |\Mvar|.}*/
  { Vertex_const_handle v;
    Halfedge_const_handle e;
    Halffacet_const_handle f;
    return  ( assign(v,h) || assign(e,h) || assign(f,h) );
  }

  Object_handle locate(const Point_3& p) const
  /*{\Mop  returns a generic handle |h| to an object (vertex, edge, facet,
  volume) of the underlying SNC which contains the point |p| in its relative 
  interior. The point |p| is contained in the set represented by |\Mvar| if 
  |\Mvar.contains(h)| is true.}*/ {
    CGAL_NEF_TRACEN( "locating point...");
    CGAL_assertion( pl() != NULL);
    Object_handle o = pl()->locate(p);
    
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    Volume_handle c;
    if(assign(v,o)) return Vertex_const_handle(v);
    if(assign(e,o)) return Halfedge_const_handle(e);
    if(assign(f,o)) return Halffacet_const_handle(f);
    if(assign(c,o)) return Volume_const_handle(c);
    return Object_handle();
  }

  /*{\Mimplementation Nef polyhedra are implemented on top of an
  extended Wuerzburg structure data structure (EWS) and use linear
  space in the number of vertices, edges and facets.  Operations like
  empty take constant time. The operations clear, complement, interior,
  closure, boundary, regularization, input and output take linear
  time. All binary set operations and comparison operations take time
  $O(N^2)$ where $N$ is the size of the output plus the size of the
  input.

  The point location operations run in linear query time without any 
  preprocessing.}*/

  /*{\Mexample Nef polyhedra are parameterized by a so called extended
  geometric kernel. There's currently only one such kernel based on a
  homogeneous representation of extended points called
  |Extended_homogeneous<NT>|.  The kernel is parameterized by a
  multiprecision integer type. The member types of |Nef_polyhedron_3<
  Extended_homogeneous_3<NT> >| map to corresponding types of the CGAL
  geometry kernel (e.g. |Nef_polyhedron::Plane_3| equals
  |CGAL::Homogeneous<leda_integer>::Plane_3| in the example below).
  \begin{Mverb}
  #include <CGAL/basic.h>
  #include <CGAL/leda_integer.h>
  #include <CGAL/Extended_homogeneous.h>
  #include <CGAL/Nef_polyhedron_3.h>

  using namespace CGAL;
  typedef  Extended_homogeneous<leda_integer>   Extended_kernel;
  typedef  Nef_polyhedron_3<Extended_kernel>    Nef_polyhedron;
  typedef  Nef_polyhedron::Plane_3              Plane_3;

  int main()
  {
    Nef_polyhedron N1(Plane_3(1,0,0,0));
    Nef_polyhedron N2(Plane_3(0,1,0,0), Nef_polyhedron::EXCLUDED);
    Nef_polyhedron N3 = N1 * N2; // line (*)
    return 0;
  }
  \end{Mverb}
  After line (*) |N3| is the intersection of |N1| and |N2|.}*/

  std::size_t bytes() {
    // bytes used for the Nef_polyhedron_3.
    return sizeof(Self) + (snc().bytes() - sizeof(SNC_structure));
  }

  std::size_t bytes_reduced() {
    // bytes used for the Nef_polyhedron_3.
    std::cout << sizeof(Self) + (snc().bytes_reduced2() - sizeof(SNC_structure)) << std::endl;
    return sizeof(Self) + (snc().bytes_reduced() - sizeof(SNC_structure));
  }

}; // end of Nef_polyhedron_3

template <typename Kernel, typename Items, typename Mark>
Nef_polyhedron_3<Kernel,Items, Mark>::
Nef_polyhedron_3( Content space, SNC_point_locator* _pl) {
  CGAL_NEF_TRACEN("construction from empty or space.");
  SNC_structure rsnc;
  empty_rep();
  set_snc(snc());
  pl() = _pl;
  if(Infi_box::extended_kernel()) {
    initialize_infibox_vertices(space);
    build_external_structure();
  } else {
    build_external_structure();
    SNC_decorator D(snc());
    D.mark(D.volumes_begin()) = (space == COMPLETE) ? 1 : 0;
  }
}

template <typename Kernel, typename Items, typename Mark>
Nef_polyhedron_3<Kernel,Items, Mark>::
Nef_polyhedron_3(const Plane_3& h, Boundary b, SNC_point_locator* _pl) {
  CGAL_NEF_TRACEN("construction from plane "<<h);
  empty_rep();
  set_snc(snc());
  SNC_constructor C(snc(), pl());
  Infi_box::create_vertices_of_box_with_plane(C,h,(b==INCLUDED));
  pl() = _pl;
  build_external_structure();
}
 
template <typename Kernel, typename Items, typename Mark>
Nef_polyhedron_3<Kernel,Items, Mark>::
Nef_polyhedron_3( const SNC_structure& W, SNC_point_locator* _pl, 
		  bool clone_pl,
		  bool clone_snc) {
  CGAL_assertion( clone_snc == true || clone_pl == false);
  // TODO: granados: define behavior when clone=false
  //  CGAL_NEF_TRACEN("construction from an existing SNC structure (clone="<<clone<<")"); 

  this->copy_on_write();
  if(clone_snc) {
    snc() = W;
    set_snc(snc());
  }
  if(clone_pl) {
    pl() = _pl->clone();
    pl()->initialize(&snc());
  } 
  else
    pl() = _pl;
}

template <typename Kernel, typename Items, typename Mark>
void
Nef_polyhedron_3<Kernel,Items, Mark>::
extract_complement() {
  CGAL_NEF_TRACEN("extract complement");
  if( this->is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_forall_vertices(v,D){
    D.mark(v) = !D.mark(v); 
    SM_decorator SM(&*v);
    SM.extract_complement();
  }
  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) D.mark(f) = !D.mark(f); 
 
  Volume_iterator c;
  CGAL_forall_volumes(c,D) 
    //    if(!(Infi_box::extended_kernel && c==D.volumes_begin()))
      D.mark(c) = !D.mark(c);
}

template <typename Kernel, typename Items, typename Mark>
void
Nef_polyhedron_3<Kernel,Items, Mark>::
extract_interior() {
  CGAL_NEF_TRACEN("extract interior");
  if (this->is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_forall_vertices(v,D){
    D.mark(v) = false;
    SM_decorator SM(&*v);
    SM.extract_interior();
  }
  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) D.mark(f) = false;

  simplify();
}

template <typename Kernel, typename Items, typename Mark>
void
Nef_polyhedron_3<Kernel,Items, Mark>::
extract_boundary() {
  CGAL_NEF_TRACEN("extract boundary");
  if (this->is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_forall_vertices(v,D) {
    D.mark(v) = true;
    SM_decorator SM(&*v);
    SM.extract_boundary();
  }
  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) D.mark(f) = true;
  Volume_iterator c;
  CGAL_forall_volumes(c,D) D.mark(c) = false;
  simplify();
}

CGAL_END_NAMESPACE

#endif //CGAL_NEF_POLYHEDRON_3_H
