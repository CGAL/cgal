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
//#include <CGAL/Nef_3/SNC_walker.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_3/SNC_SM_explorer.h>
#include <CGAL/Nef_polyhedron_S2.h>

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

#undef _DEBUG
#define _DEBUG 11
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename K, typename I> class Nef_polyhedron_3;
template <typename K, typename I> class Nef_polyhedron_3_rep;

template <typename K, typename I>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_3<K,I>&); 
template <typename K, typename I>
std::istream& operator>>(std::istream&, Nef_polyhedron_3<K,I>&); 

template <typename K, typename I>
class Nef_polyhedron_3_rep 
{ 
  typedef Nef_polyhedron_3_rep<K,I>                     Self;
  friend class Nef_polyhedron_3<K,I>;
 public:
  typedef CGAL::SNC_structure<K,I>                     SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>           SNC_decorator;
  typedef CGAL::SNC_const_decorator<SNC_structure>     SNC_const_decorator;
  typedef CGAL::SNC_constructor<SNC_structure>         SNC_constructor;
  //typedef CGAL::SNC_walker<SNC_structure>              SNC_walker;
  typedef CGAL::SNC_point_locator<SNC_structure>       SNC_point_locator;
  typedef CGAL::SNC_simplify<SNC_structure>            SNC_simplify;
#ifdef CGAL_NEF3_POINT_LOCATOR_NAIVE
  typedef CGAL::SNC_point_locator_naive<SNC_structure> SNC_point_locator_default;
#else
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_structure> SNC_point_locator_default;
#endif

  typedef typename SNC_structure::Sphere_map       Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>           SM_decorator;
  typedef CGAL::SM_const_decorator<Sphere_map>     SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<SM_decorator>     SM_overlayer;
  typedef CGAL::SM_point_locator<SNC_structure>    SM_point_locator;
  typedef CGAL::SM_io_parser<SNC_structure>        SM_io_parser;

#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef CGAL::SNC_SM_visualizor<SNC_structure>       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR

 private:
  SNC_structure snc_;
  SNC_point_locator* pl_;
  
 public:
  Nef_polyhedron_3_rep() : snc_(), pl_() {}
  ~Nef_polyhedron_3_rep() { 
    TRACEN( "Nef_polyhedron_3_rep: destroying SNC structure "<<&snc_<<
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

template <typename Kernel_, typename Items_ = SNC_items>
class Nef_polyhedron_3 : public CGAL::Handle_for< Nef_polyhedron_3_rep<Kernel_, Items_> >, 
			 public SNC_const_decorator<SNC_structure<Kernel_,Items_> >
{ 
 public:
  /*{\Mtypes 7}*/  
  
  typedef Kernel_                                     Kernel;
  typedef Items_                                      Items;
  typedef Nef_polyhedron_3<Kernel, Items>             Self;
  typedef Handle_for< Nef_polyhedron_3_rep<Kernel, Items> >   Base;
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
  struct AND { bool operator()(bool b1, bool b2, bool inverted=false) const { return b1&&b2; } };
  struct OR { bool operator()(bool b1, bool b2, bool inverted=false) const { return b1||b2; } };
  struct DIFF { bool operator()(bool b1, bool b2, bool inverted=false) const { 
    if(inverted) return !b1&&b2; return b1&&!b2; } };
  struct XOR { bool operator()(bool b1, bool b2, bool inverted=false) const 
    { return (b1&&!b2)||(!b1&&b2); } };

 public:
  typedef Nef_polyhedron_3_rep<Kernel,Items>    Nef_rep;
  typedef typename Nef_rep::SNC_structure       SNC_structure;
 protected:
  typedef typename Nef_rep::SNC_decorator       SNC_decorator;
  typedef typename Nef_rep::SNC_constructor     SNC_constructor;
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
 typedef Nef_polyhedron_S2<Kernel,Items,Sphere_map>  Nef_polyhedron_S2;
 protected:
 typedef bool   Mark;

  SNC_structure& snc() { return ptr()->snc_; } 
  const SNC_structure& snc() const { return ptr()->snc_; } 

  SNC_point_locator*& pl() { return ptr()->pl_; }
  const SNC_point_locator*& pl() const { return ptr()->pl_; }

  friend std::ostream& operator<< <>
      (std::ostream& os, Nef_polyhedron_3<Kernel,Items>& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Nef_polyhedron_3<Kernel,Items>& NP);

  typedef typename SNC_decorator::Vertex_handle    Vertex_handle;
  typedef typename SNC_decorator::Halfedge_handle  Halfedge_handle;
  typedef typename SNC_decorator::Halffacet_handle
                                                   Halffacet_handle;
  typedef typename SNC_decorator::Volume_handle    Volume_handle;

  typedef typename SNC_structure::Sphere_point                 Sphere_point;
  typedef typename SNC_structure::Sphere_segment               Sphere_segment;
  typedef typename SNC_structure::Sphere_circle                Sphere_circle;
  typedef typename SNC_structure::Vertex_const_handle          Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle        Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle       Halffacet_const_handle;
  typedef typename SNC_structure::Volume_const_handle          Volume_const_handle;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
                                  SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator 
                                  SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::Halffacet_cycle_iterator     Halffacet_cycle_iterator;
  typedef typename SNC_structure::Infi_box                     Infi_box;

  typedef typename Kernel::RT                       RT;

 public:
  typedef typename SM_decorator::SVertex_handle    SVertex_handle;
  typedef typename SM_decorator::SHalfedge_handle  SHalfedge_handle;
  typedef typename SM_decorator::SFace_handle      SFace_handle;
  typedef typename SM_decorator::SVertex_const_handle
                                                   SVertex_const_handle;
  typedef typename SM_decorator::SHalfedge_const_handle
                                                   SHalfedge_const_handle;
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
  typedef typename SM_decorator::SFace_const_iterator     
                                                   SFace_const_iterator;

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

  Nef_polyhedron_3(const Nef_polyhedron_3<Kernel,Items>& N1) : Base(N1) {
    set_snc(snc());
  } 

  Nef_polyhedron_3& operator=(const Nef_polyhedron_3<Kernel,Items>& N1) { 
    Base::operator=(N1);
    set_snc(snc());
    return (*this); 
  }

  ~Nef_polyhedron_3() { 
    TRACEN("~Nef_polyhedron_3: destructor called for snc "<<&snc()<<
	   ", pl "<<pl());
  }
  
  typedef Polyhedron_3< Kernel> Polyhedron;
  Nef_polyhedron_3( Polyhedron& P, 
		    SNC_point_locator* _pl = new SNC_point_locator_default) {
    TRACEN("construction from Polyhedron_3");
    SNC_structure rsnc;
    *this = Nef_polyhedron_3(rsnc, _pl, false);
    initialize_infibox_vertices(EMPTY);
    polyhedron_3_to_nef_3
      < Polyhedron, SNC_structure>( P, snc());
    build_external_structure();
    simplify();
    set_snc(snc());
    //    CGAL_assertion(orientation() == 1);
  }
  
 protected:  
  template <class HDS>
  class Build_polyhedron : public CGAL::Modifier_base<HDS> {
    
    class Visitor {

      const Object_index<Vertex_iterator>& VI;
      Polyhedron_incremental_builder_3<HDS>& B;
      SNC_decorator& D;
      
    public:
      Visitor(Polyhedron_incremental_builder_3<HDS>& BB,
	      SNC_decorator& sd,
	      Object_index<Vertex_iterator>& vi) : VI(vi), B(BB), D(sd){}

      void visit(Halffacet_handle opposite_facet) {

	TRACEN("Build_polyhedron: visit facet " << D.plane(opposite_facet));
 
	CGAL_assertion(Infi_box::is_standard(D.plane(opposite_facet)));
	
	SHalfedge_handle se;
	Halffacet_cycle_iterator fc;
     	
	Halffacet_handle f = D.twin(opposite_facet);

	B.begin_facet();
	fc = f->facet_cycles_begin();
	assign(se,fc);
	CGAL_assertion(se!=0);
	SHalfedge_around_facet_circulator hc_start(se);
	SHalfedge_around_facet_circulator hc_end(hc_start);
	CGAL_For_all(hc_start,hc_end) {
	  TRACEN("   add vertex " << D.vertex(hc_start)->point());
	  B.add_vertex_to_facet(VI[D.vertex(hc_start)]);
	}
	B.end_facet();
      }

      void visit(SFace_handle s) {}
      void visit(Halfedge_handle e) {}
      void visit(Vertex_handle v) {}
    };

  public:

    SNC_structure& snc;
    Object_index<Vertex_iterator> VI;

    Build_polyhedron(SNC_structure& s) : 
      snc(s), VI(s.vertices_begin(),s.vertices_end(),'V') {}
    
    void operator()( HDS& hds) {

      SNC_decorator D(snc);

      Polyhedron_incremental_builder_3<HDS> B(hds, true);

      int skip_volumes;
      if(Infi_box::extended_kernel()) {
	B.begin_surface(D.number_of_vertices()-8, 
			D.number_of_facets()-6,
			D.number_of_edges()-12);
	skip_volumes = 2;
      }
      else {
	B.begin_surface(D.number_of_vertices(), 
			D.number_of_facets(),
			D.number_of_edges());
	skip_volumes = 1;
      }
      
      int vertex_index = 0;
      Vertex_iterator v;
      CGAL_forall_vertices(v,snc) {
	if(Infi_box::is_standard(v->point())) {
	  VI[v]=vertex_index++;
	  B.add_vertex(v->point());
	}
      }     
      
      Visitor V(B,D,VI);
      Volume_handle c;
      CGAL_forall_volumes(c,snc)
	if(skip_volumes-- <= 0)
	  D.visit_shell_objects(SFace_handle(c->shells_begin()),V);
      B.end_surface();
    }

  };

  int orientation() {
    // Function does not work correctly with every Nef_polyhedron. 
    // It works correctly if v_max is 2-manifold

    typedef typename SM_decorator::SHalfedge_around_sface_circulator
      SHalfedge_around_sface_circulator;

    if ( number_of_vertices() == 0)
      return 0;
    
    // Find top-most vertex v_max (top-most also works fine for terrains)
    SNC_decorator D(snc());
    Vertex_handle v_max = D.vertices_begin();
    Vertex_handle vh;
    CGAL_forall_vertices(vh, D) {
      if (Infi_box::is_standard(D.point(vh)) && 
	  D.point(vh).z() > D.point(v_max).z())
	v_max = vh;
    }
    SM_decorator SD(&*v_max);

    // Add up z-coord. of all normalized normal vectors of the incident facets
    
    double z=1.0;
    bool first = true;
    SFace_handle sf;
    CGAL_forall_sfaces(sf, SD) {
      if(D.volume(sf) != Infi_box::getNirvana(snc())) continue;
      if(!first) return 0; // orientation can't be decided via v_max
      SHalfedge_around_sface_circulator estart(sf->sface_cycles_begin()), 
	                                eend(estart);
      CGAL_For_all(estart,eend) {
	Sphere_circle c(SD.circle(estart));
	double delta_z = CGAL::to_double(c.orthogonal_vector().hz());
	Point_3 target = ORIGIN + c.orthogonal_vector();
	Segment_3 s(Point_3(0,0,0), target);
	delta_z /= sqrt(to_double(s.squared_length()));
	z += delta_z;	
      }
      first = false;
    }
    
    return sign(z);
  }

 public:

  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  void convert_to_Polyhedron(Polyhedron& P) {
    CGAL_precondition(is_simple());
    Build_polyhedron<HalfedgeDS> bp(snc());    
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
      
    if (++(v->sfaces_begin()) != v->sfaces_last())
      return false;

    return true;
  }

  bool is_facet_simple(const Halffacet_handle& f) {
    
    bool found_first = false;
    Halffacet_cycle_iterator it; 
    CGAL_forall_facet_cycles_of(it,f)
      if (found_first || !it.is_shalfedge())
	return false;
      else
	found_first = true;
   
    return true;
  }

 public:
   
  void clear(Content space = EMPTY)
    { *this = Nef_polyhedron_3(space, pl->clone()); }
  /*{\Mop makes |\Mvar| the empty set if |space == EMPTY| and the
  full space if |space == COMPLETE|.}*/

  bool is_empty() //const
    /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/ {
    return snc().has_bbox_only();
  }

  bool is_space() //const
  /*{\Mop returns true if |\Mvar| is the whole space, false otherwise.}*/
  //{ SM_const_decorator D(pm());
  //  typename PM_const_decorator::Volume_const_iterator v = D.volumes_begin();
  { //SNC_structure snc_non_const = snc(); // const_decorator not implemented
    SNC_decorator D(snc());
    Volume_iterator v = D.volumes_begin(); // it should be const_iterator
    return (D.number_of_vertices()==8 &&
            D.number_of_edges()==12 &&
            D.number_of_facets()==6 && 
	    D.number_of_volumes()==2 && 
            D.mark(++v) == true);
  }

  /*{\Xtext \headerline{Destructive Operations}}*/

 protected:
  void clone_rep() { *this = Nef_polyhedron_3<Kernel,Items>(snc(), pl()); }
  void empty_rep() { 
    SNC_structure rsnc;
    *this = Nef_polyhedron_3<Kernel,Items>(rsnc, new SNC_point_locator_default);
  }

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
    TRACEN( "simplify(): structure simplified? "<<simplified);
    
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
      TRACEN("simplify(): point locator structure updated? " << updated);
#else
      SNC_point_locator* old_pl = pl();
      pl() = pl()->clone();
      pl()->initialize(&snc());
      delete old_pl;
#endif
    }
  }

 public:
  Nef_polyhedron_S2 get_sphere_map(Vertex_const_handle v) {
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
  { TRACEN("extract closure");
    if( is_shared()) clone_rep();
    extract_complement();
    extract_interior();
    extract_complement();
  }

  void extract_regularization()
  /*{\Xop converts |\Mvar| to its regularization. }*/
  { TRACEN("extract regularization");
    if( is_shared()) clone_rep();
    extract_interior();
    extract_closure();
  }

  /*{\Mtext \headerline{Constructive Operations}}*/

  Nef_polyhedron_3<Kernel,Items> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane. }*/
  { Nef_polyhedron_3<Kernel,Items> res = *this;
    res.extract_complement();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> interior() const
  /*{\Mop    returns the interior of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> closure() const
  /*{\Mop returns the closure of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> boundary() const
  /*{\Mop returns the boundary of |\Mvar|. }*/
  { Nef_polyhedron_3<Kernel,Items> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> regularization() const
  /*{\Mop    returns the regularized polyhedron (closure of 
             the interior).}*/
  { Nef_polyhedron_3<Kernel,Items> res = *this;
    res.extract_regularization();
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> intersection(Nef_polyhedron_3<Kernel,Items>& N1 )
    /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/ {
    TRACEN(" intersection between nef3 "<<&*this<<" and "<<&N1);
    AND _and;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, AND> tests_impl;
    SNC_structure rsnc;
    SNC_decorator D( snc(), pl());
    Nef_polyhedron_3<Kernel,Items> res(rsnc, new SNC_point_locator_default, false);
    D.binary_operation( N1.snc(), N1.pl(), _and, res.snc(), res.pl());
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> join(Nef_polyhedron_3<Kernel,Items>& N1) 
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/ { 
    TRACEN(" join between nef3 "<<&*this<<" and "<<&N1);
    OR _or;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, OR> tests_impl;
    SNC_structure rsnc;
    SNC_decorator D( snc(), pl());
    Nef_polyhedron_3<Kernel,Items> res(rsnc, new SNC_point_locator_default, false);
    D.binary_operation( N1.snc(), N1.pl(), _or, res.snc(), res.pl());
    //delete rpl; // TODO: analize how the improve the Nef_3 constructor so this instruction is not needed
    return res;
  }

  Nef_polyhedron_3<Kernel,Items> difference(Nef_polyhedron_3<Kernel,Items>& N1)
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/ { 
    TRACEN(" difference between nef3 "<<&*this<<" and "<<&N1);
    DIFF _diff;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, DIFF> tests_impl;
    SNC_structure rsnc;
    SNC_decorator D( snc(), pl());
    Nef_polyhedron_3<Kernel,Items> res(rsnc, new SNC_point_locator_default, false);
    D.binary_operation( N1.snc(), N1.pl(), _diff, res.snc(), res.pl());
    return res;
  }    

  Nef_polyhedron_3<Kernel,Items> symmetric_difference(Nef_polyhedron_3<Kernel,Items>& N1)
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/ {
    TRACEN(" symmetic difference between nef3 "<<&*this<<" and "<<&N1);
    XOR _xor;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, XOR> tests_impl;
    SNC_structure rsnc;
    SNC_decorator D( snc(), pl());
    Nef_polyhedron_3<Kernel,Items> res(rsnc, new SNC_point_locator_default, false);
    D.binary_operation( N1.snc(), N1.pl(), _xor, res.snc(), res.pl());
    return res;
  }


  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement}. There are also the corresponding
  modification operations |*=,+=,-=,^=|.}*/

  Nef_polyhedron_3<Kernel,Items>  operator*(const Nef_polyhedron_3<Kernel,Items>& N1) const 
  { return intersection(N1); }

  Nef_polyhedron_3<Kernel,Items>  operator+(Nef_polyhedron_3<Kernel,Items>& N1) 
  { return join(N1); }

  Nef_polyhedron_3<Kernel,Items>  operator-(Nef_polyhedron_3<Kernel,Items>& N1) 
  { return difference(N1); }

  Nef_polyhedron_3<Kernel,Items>  operator^(Nef_polyhedron_3<Kernel,Items>& N1) 
  { return symmetric_difference(N1); }

  Nef_polyhedron_3<Kernel,Items>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_3<Kernel,Items>& operator*=(Nef_polyhedron_3<Kernel,Items>& N1)
  { *this = intersection(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items>& operator+=(Nef_polyhedron_3<Kernel,Items>& N1)
  { *this = join(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items>& operator-=(Nef_polyhedron_3<Kernel,Items>& N1)
  { *this = difference(N1); return *this; }

  Nef_polyhedron_3<Kernel,Items>& operator^=(Nef_polyhedron_3<Kernel,Items>& N1)
  { *this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality.}*/

  bool operator==(Nef_polyhedron_3<Kernel,Items>& N1) 
  { TRACEN(" equality comparision between nef3 "<<&*this<<" and "<<&N1);
    return symmetric_difference(N1).is_empty(); }

  bool operator!=(Nef_polyhedron_3<Kernel,Items>& N1)
  { TRACEN(" inequality comparision between nef3 "<<&*this<<" and "<<&N1);
    return !operator==(N1); }  

  bool operator<( Nef_polyhedron_3<Kernel,Items>& N1) 
  { return !N1.difference(*this).is_empty() && difference(N1).is_empty(); } 

  bool operator>( Nef_polyhedron_3<Kernel,Items>& N1)    
  { return difference(*this).is_empty() && !difference(N1).is_empty(); } 

  bool operator<=( Nef_polyhedron_3<Kernel,Items>& N1) 
  { return difference(N1).is_empty(); } 

  bool operator>=( Nef_polyhedron_3<Kernel,Items>& N1) 
  { return N1.difference(*this).is_empty(); } 

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
    
    Vertex_iterator vi;
    CGAL_forall_vertices( vi, snc()) {
      
      TRACEN("transform vertex ");
      if (is_standard(vi)) {
	vi->point() = vi->point().transform( aff);
	SM_decorator sdeco(&*vi);
	sdeco.transform( linear);
      }
    }

    Halffacet_iterator fi;
    CGAL_forall_halffacets(fi,snc()) {
      if ( deco.is_standard( fi))
        fi->plane() = fi->plane().transform( aff);
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
  
  /*{\Mtext \headerline{Exploration}
  As Nef polyhedra are the result of forming complements 
  and intersections starting from a set |H| of halfspaces which are
  defined by oriented planes in three space. The corresponding 
  structure is represented by an extended wuerzburg structure 
  $W = (V,E,F,C)$. For topological queries within |W| the following 
  types and operations allow exploration access to this structure.}*/

  /*{\Mtypes 3}*/
    
    typedef typename Nef_rep::SNC_const_decorator      SNC_const_decorator;
    typedef typename Nef_rep::SM_const_decorator       SM_const_decorator;
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
    { CGAL_assertion_msg( 0, "not implemented."); }

  bool contained_in_boundary(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the $2$-skeleton
  of |\Mvar|.}*/
  { Vertex_const_handle v;
    Halfedge_const_handle e;
    Halffacet_const_handle f;
    return  ( assign(v,h) || assign(e,h) || assign(f,h) );
  }

  Object_handle locate(const Point_3& p)
  /*{\Mop  returns a generic handle |h| to an object (vertex, edge, facet,
  volume) of the underlying SNC which contains the point |p| in its relative 
  interior. The point |p| is contained in the set represented by |\Mvar| if 
  |\Mvar.contains(h)| is true.}*/ {
    TRACEN( "locating point...");
    CGAL_assertion( pl() != NULL);
    return pl()->locate(p);
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
  |Extended_homogeneous_3<NT>|.  The kernel is parameterized by a
  multiprecision integer type. The member types of |Nef_polyhedron_3<
  Extended_homogeneous_3<NT> >| map to corresponding types of the CGAL
  geometry kernel (e.g. |Nef_polyhedron::Plane_3| equals
  |CGAL::Homogeneous<leda_integer>::Plane_3| in the example below).
  \begin{Mverb}
  #include <CGAL/basic.h>
  #include <CGAL/leda_integer.h>
  #include <CGAL/Extended_homogeneous_3.h>
  #include <CGAL/Nef_polyhedron_3.h>

  using namespace CGAL;
  typedef  Extended_homogeneous_3<leda_integer> Extended_kernel;
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
    cout << sizeof(Self) + (snc().bytes_reduced2() - sizeof(SNC_structure)) << std::endl;
    return sizeof(Self) + (snc().bytes_reduced() - sizeof(SNC_structure));
  }

}; // end of Nef_polyhedron_3

template <typename Kernel, typename Items>
Nef_polyhedron_3<Kernel,Items>::
Nef_polyhedron_3( Content space, SNC_point_locator* _pl) {
  TRACEN("construction from empty or space.");
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

template <typename Kernel, typename Items>
Nef_polyhedron_3<Kernel,Items>::
Nef_polyhedron_3(const Plane_3& h, Boundary b, SNC_point_locator* _pl) {
  TRACEN("construction from plane "<<h);
  empty_rep();
  set_snc(snc());
  SNC_constructor C(snc(), pl());
  Infi_box::create_vertices_of_box_with_plane(C,h,(b==INCLUDED));
  pl() = _pl;
  build_external_structure();
}
 
template <typename Kernel, typename Items>
Nef_polyhedron_3<Kernel,Items>::
Nef_polyhedron_3( const SNC_structure& W, SNC_point_locator* _pl, 
		  bool clone_pl,
		  bool clone_snc) {
  CGAL_assertion( clone_snc == true || clone_pl == false);
  // TODO: granados: define behavior when clone=false
  //  TRACEN("construction from an existing SNC structure (clone="<<clone<<")"); 
  copy_on_write();
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

template <typename Kernel, typename Items>
void
Nef_polyhedron_3<Kernel,Items>::
extract_complement() {
  TRACEN("extract complement");
  if( is_shared()) clone_rep();
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

template <typename Kernel, typename Items>
void
Nef_polyhedron_3<Kernel,Items>::
extract_interior() {
  TRACEN("extract interior");
  if (is_shared()) clone_rep();
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

template <typename Kernel, typename Items>
void
Nef_polyhedron_3<Kernel,Items>::
extract_boundary() {
  TRACEN("extract boundary");
  if (is_shared()) clone_rep();
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
