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
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_3/SNC_ray_shooter.h>
//#include <CGAL/Nef_3/SNC_walker.h>
#ifdef CGAL_NEF3_VISUALIZOR
#include <CGAL/Nef_3/SNC_visualizor_OGL.h>
#endif // CGAL_NEF3_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_decorator.h>
#include <CGAL/Nef_3/SNC_SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_io_parser.h>
#include <CGAL/Nef_3/SNC_explorer.h>
#include <CGAL/Nef_3/SNC_SM_explorer.h>

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

template <typename T> class Nef_polyhedron_3;
template <typename T> class Nef_polyhedron_3_rep;

template <typename T>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_3<T>&); 
template <typename T>
std::istream& operator>>(std::istream&, Nef_polyhedron_3<T>&); 

template <typename T>
class Nef_polyhedron_3_rep : public Ref_counted
{ 
  typedef Nef_polyhedron_3_rep<T>                     Self;
  friend class Nef_polyhedron_3<T>;
  typedef CGAL::SNC_structure<T>                       SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>           SNC_decorator;
  typedef CGAL::SNC_const_decorator<SNC_structure>     SNC_const_decorator;
  typedef CGAL::SNC_constructor<SNC_structure>         SNC_constructor;
  typedef CGAL::SNC_ray_shooter<SNC_structure> SNC_ray_shooter; // DEPRECATED
  //typedef CGAL::SNC_walker<SNC_structure>              SNC_walker;
  typedef CGAL::SNC_io_parser<SNC_structure>           SNC_io_parser;
  typedef CGAL::SNC_point_locator<SNC_structure>       SNC_point_locator;
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_structure>
    //typedef CGAL::SNC_point_locator_naive<SNC_structure>
    SNC_point_locator_default;

#ifdef CGAL_NEF3_VISUALIZOR
  typedef CGAL::SNC_visualizor_OGL<SNC_structure>      SNC_visualizor;
#endif // CGAL_NEF3_VISUALIZOR
  typedef CGAL::SNC_SM_decorator<SNC_structure>        SM_decorator;
  typedef CGAL::SNC_SM_const_decorator<SNC_structure>  SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<SNC_structure>        SM_overlayer;
  typedef CGAL::SNC_SM_point_locator<SNC_structure>    SM_point_locator;
  typedef CGAL::SNC_SM_io_parser<SNC_structure>        SM_io_parser;

#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef CGAL::SNC_SM_visualizor<SNC_structure>       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR

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

template <typename T>
class Nef_polyhedron_3 : public CGAL::Handle_for< Nef_polyhedron_3_rep<T> >
{ 
public:
  typedef T Extended_kernel; 
  static  T EK; // static extended kernel

  /*{\Mtypes 7}*/  
  typedef Nef_polyhedron_3<T>                   Self;
  typedef Handle_for< Nef_polyhedron_3_rep<T> > Base;
  typedef typename T::Kernel                    Kernel;
  typedef typename T::Point_3                   Point_3;
  typedef typename T::Plane_3                   Plane_3;
  typedef typename T::Vector_3                  Vector_3;
  typedef typename Kernel::Aff_transformation_3      Aff_transformation_3;

  typedef bool Mark;
  /*{\Xtypemember marking set membership or exclusion.}*/
 
  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  typedef enum { EMPTY=0, COMPLETE=1 } Content;
  /*{\Menum construction selection}*/

  typedef enum { DEFAULT, NAIVE, WALKING, SPATIAL_SUBDIVISION  } Location_mode;
  /*{\Menum selection flag for the point location mode.}*/

protected:
  struct AND { bool operator()(bool b1, bool b2) const { return b1&&b2; } };
  struct OR { bool operator()(bool b1, bool b2) const { return b1||b2; } };
  struct DIFF { bool operator()(bool b1, bool b2) const { return b1&&!b2; } };
  struct XOR { bool operator()(bool b1, bool b2) const 
    { return (b1&&!b2)||(!b1&&b2); } };

  typedef Nef_polyhedron_3_rep<T>               Nef_rep;
  typedef typename Nef_rep::SNC_structure       SNC_structure;
  typedef typename Nef_rep::SNC_decorator       SNC_decorator;
  typedef typename Nef_rep::SNC_constructor     SNC_constructor;
  typedef typename Nef_rep::SNC_ray_shooter SNC_ray_shooter; // DEPRECATED
  //typedef typename Nef_rep::SNC_walker          SNC_walker;
  typedef typename Nef_rep::SNC_io_parser       SNC_io_parser;
  typedef typename Nef_rep::SNC_point_locator   SNC_point_locator;
  typedef typename Nef_rep::SNC_point_locator_default 
    SNC_point_locator_default;

#ifdef CGAL_NEF3_VISUALIZOR
  typedef typename Nef_rep::SNC_visualizor      SNC_visualizor;
#endif // CGAL_NEF3_VISUALIZOR
  typedef typename Nef_rep::SM_decorator        SM_decorator;
  typedef typename Nef_rep::SM_const_decorator  SM_const_decorator;
  typedef typename Nef_rep::SM_overlayer        SM_overlayer;
  typedef typename Nef_rep::SM_point_locator    SM_point_locator;
  typedef typename Nef_rep::SM_io_parser        SM_io_parser;
#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef typename Nef_rep::SM_visualizor       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR

  SNC_structure& snc() { return ptr()->snc_; } 
  const SNC_structure& snc() const { return ptr()->snc_; } 

  SNC_point_locator*& pl() { return ptr()->pl_; }
  const SNC_point_locator*& pl() const { return ptr()->pl_; }

  friend std::ostream& operator<< <>
      (std::ostream& os, Nef_polyhedron_3<T>& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Nef_polyhedron_3<T>& NP);

  typedef typename SNC_decorator::Vertex_handle    Vertex_handle;
  typedef typename SNC_decorator::Halfedge_handle  Halfedge_handle;
  typedef typename SNC_decorator::Halffacet_handle
                                                   Halffacet_handle;
  typedef typename SNC_decorator::Volume_handle    Volume_handle;

#define USING(t) typedef typename SNC_structure::t t
  USING(Sphere_point);
  USING(Sphere_segment);
  USING(Sphere_circle);
  USING(Vertex_const_handle);
  USING(Halfedge_const_handle);
  USING(Halffacet_const_handle);
  USING(Volume_const_handle);
  USING(SHalfedge_around_svertex_circulator);
  USING(SHalfedge_around_facet_circulator);
  USING(SHalfedge_around_facet_const_circulator);
  USING(Halffacet_cycle_iterator);
  USING(Infi_box);
# undef USING

  typedef typename Kernel::RT                       RT;

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

  void clear_box_marks() 
    /* unset all frame marks */ {
    SNC_decorator D(snc());
    Volume_iterator c = D.volumes_begin(); 
    D.mark(c) = false;
    D.clear_outer_box_marks();
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

  Nef_polyhedron_3(const Nef_polyhedron_3<T>& N1) : Base(N1) {} 

  Nef_polyhedron_3& operator=(const Nef_polyhedron_3<T>& N1)
  { Base::operator=(N1); return (*this); }

  ~Nef_polyhedron_3() { 
    TRACEN("~Nef_polyhedron_3: destructor called for snc "<<&snc()<<
	   ", pl "<<pl());
  }
  
  typedef Polyhedron_3< Kernel> Polyhedron;
  Nef_polyhedron_3( Polyhedron& P, 
		    SNC_point_locator* _pl = new SNC_point_locator_default) {
    TRACEN("construction from Polyhedron_3");
    initialize_infibox_vertices(EMPTY);
    polyhedron_3_to_nef_3
      < Polyhedron, SNC_structure, SM_point_locator>( P, snc());
    pl() = _pl;
    build_external_structure();
    simplify();
  }
  
  Nef_polyhedron_3( const char* filename, 
		    SNC_point_locator* _pl = new SNC_point_locator_default) {
    bool OK = true;
    std::ifstream in(filename);
    OK = OK && in;
    OK = OK && snc().load(in);
    if(!OK) 
      std::cerr << "Failure while loading data" << std::endl;
    pl() = _pl;
    pl()->initialize(&snc());
  }

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
      
      B.begin_surface(D.number_of_vertices()-8, 
		      D.number_of_facets()-6,
		      D.number_of_edges()-12);
      
      int vertex_index = -8;
      Vertex_iterator v;
      CGAL_nef3_forall_vertices(v,snc) {
	if(vertex_index >= 0) {
	  VI[v]=vertex_index;
	  B.add_vertex(v->point());
	}
	vertex_index++;
      }     
        
      int outer_volume = 0;
      //Shell_entry_iterator it;
      Visitor V(B,D,VI);
      Volume_handle c;
      CGAL_nef3_forall_volumes(c,snc) {
	if(outer_volume++ > 1)
	  D.visit_shell_objects(SFace_handle(c->shells_begin()),V);
      }
   
      B.end_surface();
    }

  };
  
  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  void convert_to_Polyhedron(Polyhedron& P) {
       
    CGAL_precondition(is_simple());
    Build_polyhedron<HalfedgeDS> bp(snc());    
    P.delegate(bp);
  }

  void dump(std::ostream& os = std::cerr) { SNC_io_parser::dump( snc(), os); }

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
    CGAL_nef3_forall_edges(e,snc())
      if(!is_edge_2manifold(e))
	return false;

    Vertex_iterator v;
    CGAL_nef3_forall_vertices(v,snc())
      if(!is_vertex_2manifold(v))
	return false;

    Halffacet_iterator f;
    CGAL_nef3_forall_halffacets(f,snc())
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
    CGAL_nef3_forall_facet_cycles_of(it,f)
      if (found_first || !it.is_shalfedge())
	return false;
      else
	found_first = true;
   
    return true;
  }

 public:

  void visualize() { 
#ifdef CGAL_NEF3_VISUALIZOR
    SNC_visualizor sncv( snc());
    sncv.draw();
    //OGL::polyhedra_.back().debug();
    OGL::start_viewer();
#endif // CGAL_NEF3_VISUALIZOR
  }

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
  void clone_rep() { *this = Nef_polyhedron_3<T>(snc(), pl()); }

  Nef_polyhedron_3( const SNC_structure& W, 
		    SNC_point_locator* _pl = new SNC_point_locator_default,
		    bool clone = true);
  /*{\Xcreate makes |\Mvar| a new object.  If |cloneit==true| then the
  underlying structure of |W| is copied into |\Mvar|.}*/
  // TODO: granados: define behavior when clone=false

  /*{\Moperations 4 3 }*/

  void simplify() {
    bool simplified = snc().simplify();
    TRACEN( "simpliy(): structure simplified? "<<simplified);
    
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
      CGAL_nef3_forall_vertices( v, snc()) V[Vertex_handle(v)] = true;
      CGAL_nef3_forall_halfedges( e, snc()) E[Halfedge_handle(e)] = true;
      CGAL_nef3_forall_halffacets( f, snc()) F[Halffacet_handle(f)] = true;
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

  template< class Selection >
  Nef_polyhedron_3<T> bin_op( Nef_polyhedron_3<T>& N1,
                              const Selection& BOP,
                              binop_intersection_tests<SNC_decorator,
                                                       Selection>& tests_impl )
  {
    SNC_structure rsnc;
    SNC_decorator D(snc());
    D.binary_operation( N1.snc(), BOP, rsnc, tests_impl );
    return Nef_polyhedron_3<T> (rsnc);
  }
  
 public:
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

  Nef_polyhedron_3<T> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_complement();
    return res;
  }

  Nef_polyhedron_3<T> interior() const
  /*{\Mop    returns the interior of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_3<T> closure() const
  /*{\Mop returns the closure of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_3<T> boundary() const
  /*{\Mop returns the boundary of |\Mvar|. }*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_3<T> regularization() const
  /*{\Mop    returns the regularized polyhedron (closure of 
             the interior).}*/
  { Nef_polyhedron_3<T> res = *this;
    res.extract_regularization();
    return res;
  }

  Nef_polyhedron_3<T>
  intersection(Nef_polyhedron_3<T>& N1,
               CGAL::binop_intersection_tests<SNC_decorator, AND>& tests_impl )
  {
    TRACEN(" intersection between nef3 "<<&*this<<" and "<<&N1);
    return bin_op( N1, AND(), tests_impl );
  }

  Nef_polyhedron_3<T>
  join(Nef_polyhedron_3<T>& N1,
       CGAL::binop_intersection_tests<SNC_decorator, OR>& tests_impl )
  {
    TRACEN(" join between nef3 "<<&*this<<" and "<<&N1);
    return bin_op( N1, OR(), tests_impl );
  }

  Nef_polyhedron_3<T>
  difference(Nef_polyhedron_3<T>& N1,
             CGAL::binop_intersection_tests<SNC_decorator, DIFF>& tests_impl )
  {
    TRACEN(" difference between nef3 "<<&*this<<" and "<<&N1);
    return bin_op( N1, DIFF(), tests_impl );
  }

  Nef_polyhedron_3<T>
  symmetric_difference(Nef_polyhedron_3<T>& N1,
                       CGAL::binop_intersection_tests<SNC_decorator,
                                                      XOR>& tests_impl )
  {
    TRACEN(" symmetric difference between nef3 "<<&*this<<" and "<<&N1);
    return bin_op( N1, XOR(), tests_impl );
  }

  Nef_polyhedron_3<T> intersection(Nef_polyhedron_3<T>& N1 )
    /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/ {
    TRACEN(" intersection between nef3 "<<&*this<<" and "<<&N1);
    AND _and;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, AND> tests_impl;
    SNC_structure rsnc;
    //<<<<<<< Nef_polyhedron_3.h
    SNC_point_locator *rpl = pl()->clone();
    SNC_decorator D( snc(), pl());
    D.binary_operation( N1.snc(), N1.pl(), _and, rsnc, rpl);
    Nef_polyhedron_3<T> res(rsnc, rpl);
    //=======
    //SNC_decorator D(snc());
    //D.binary_operation( N1.snc(), _and, rsnc, tests_impl );
    //Nef_polyhedron_3<T> res(rsnc);
    //>>>>>>> 1.43
    //    res.clear_box_marks();
    return res;
  }

  Nef_polyhedron_3<T> join(Nef_polyhedron_3<T>& N1) 
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/ { 
    TRACEN(" join between nef3 "<<&*this<<" and "<<&N1);
    OR _or;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, OR> tests_impl;
    SNC_structure rsnc;
    //<<<<<<< Nef_polyhedron_3.h
    SNC_point_locator *rpl = pl()->clone();
    SNC_decorator D( snc(), pl());
    D.binary_operation( N1.snc(), N1.pl(), _or, rsnc, rpl);
    Nef_polyhedron_3<T> res(rsnc, rpl);
    //delete rpl; // TODO: analize how the improve the Nef_3 constructor so this instruction is not needed
    //=======
    //SNC_decorator D(snc());
    //D.binary_operation( N1.snc(), _or, rsnc, tests_impl);
    //Nef_polyhedron_3<T> res(rsnc);
    //>>>>>>> 1.43
    //    res.clear_box_marks();
    return res;
  }

  Nef_polyhedron_3<T> difference(Nef_polyhedron_3<T>& N1)
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/ { 
    TRACEN(" difference between nef3 "<<&*this<<" and "<<&N1);
    DIFF _diff;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, DIFF> tests_impl;
    SNC_structure rsnc;
    //<<<<<<< Nef_polyhedron_3.h
    SNC_point_locator *rpl = pl()->clone();
    SNC_decorator D( snc(), pl());
    D.binary_operation( N1.snc(), N1.pl(), _diff, rsnc, rpl);
    Nef_polyhedron_3<T> res(rsnc, rpl);
    //=======
    //SNC_decorator D(snc());
    //D.binary_operation( N1.snc(), _diff, rsnc, tests_impl);
    //Nef_polyhedron_3<T> res(rsnc);
    //>>>>>>> 1.43
    //    res.clear_box_marks();
    return res;
  }    

  Nef_polyhedron_3<T> symmetric_difference(Nef_polyhedron_3<T>& N1)
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/ {
    TRACEN(" symmetic difference between nef3 "<<&*this<<" and "<<&N1);
    XOR _xor;
    //CGAL::binop_intersection_tests_allpairs<SNC_decorator, XOR> tests_impl;
    SNC_structure rsnc;
    //<<<<<<< Nef_polyhedron_3.h
    SNC_point_locator *rpl = pl()->clone();
    SNC_decorator D( snc(), pl());
    D.binary_operation( N1.snc(), N1.pl(), _xor, rsnc, rpl);
    Nef_polyhedron_3<T> res(rsnc, rpl);
    //=======
    //SNC_decorator D(snc());
    //D.binary_operation( N1.snc(), _xor, rsnc, tests_impl);
    //Nef_polyhedron_3<T> res(rsnc);
    //>>>>>>> 1.43
    //    res.clear_box_marks();
    return res;
  }


  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement}. There are also the corresponding
  modification operations |*=,+=,-=,^=|.}*/

  Nef_polyhedron_3<T>  operator*(Nef_polyhedron_3<T>& N1) 
  { return intersection(N1); }

  Nef_polyhedron_3<T>  operator+(Nef_polyhedron_3<T>& N1) 
  { return join(N1); }

  Nef_polyhedron_3<T>  operator-(Nef_polyhedron_3<T>& N1) 
  { return difference(N1); }

  Nef_polyhedron_3<T>  operator^(Nef_polyhedron_3<T>& N1) 
  { return symmetric_difference(N1); }

  Nef_polyhedron_3<T>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_3<T>& operator*=(Nef_polyhedron_3<T>& N1)
  { this = intersection(N1); return *this; }

  Nef_polyhedron_3<T>& operator+=(Nef_polyhedron_3<T>& N1)
  { this = join(N1); return *this; }

  Nef_polyhedron_3<T>& operator-=(Nef_polyhedron_3<T>& N1)
  { this = difference(N1); return *this; }

  Nef_polyhedron_3<T>& operator^=(Nef_polyhedron_3<T>& N1)
  { this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality.}*/

  bool operator==(Nef_polyhedron_3<T>& N1) 
  { TRACEN(" equality comparision between nef3 "<<&*this<<" and "<<&N1);
    return symmetric_difference(N1).is_empty(); }

  bool operator!=(Nef_polyhedron_3<T>& N1)
  { TRACEN(" inequality comparision between nef3 "<<&*this<<" and "<<&N1);
    return !operator==(N1); }  

  bool operator<( Nef_polyhedron_3<T>& N1) 
  { return !N1.difference(*this).is_empty() && difference(N1).is_empty(); } 

  bool operator>( Nef_polyhedron_3<T>& N1)    
  { return difference(*this).is_empty() && !difference(N1).is_empty(); } 

  bool operator<=( Nef_polyhedron_3<T>& N1) 
  { return difference(N1).is_empty(); } 

  bool operator>=( Nef_polyhedron_3<T>& N1) 
  { return N1.difference(*this).is_empty(); } 

  void transform( const Aff_transformation_3& aff) {
    
    //      SETDTHREAD(11*23);
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
    CGAL_nef3_forall_vertices( vi, snc()) {
      
      TRACEN("transform vertex ");
      if ( ! deco.is_infbox_vertex(vi)) {
	vi->point() = vi->point().transform( aff);
	SM_decorator sdeco(vi);
	sdeco.transform( linear);
      }
    }
    deco.compute_all_marks_of_halfspheres();
    
    Halffacet_iterator fi;
    CGAL_nef3_forall_halffacets(fi,snc()) {
      if ( ! deco.is_infbox_facet( fi))
                fi->plane() = fi->plane().transform( aff);
    }
    
    SNC_point_locator* old_pl = pl();
    pl() = pl()->clone();
    pl()->initialize(&snc());
    delete old_pl;
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
    typedef CGAL::SNC_explorer<SNC_const_decorator>    SNC_explorer;
    typedef typename Nef_rep::SM_const_decorator       SM_const_decorator;
    typedef CGAL::SNC_SM_explorer<SM_const_decorator>  SM_explorer;

    SNC_explorer SNCexplorer() const { 
      SNC_const_decorator SCD(snc());
      return SNC_explorer(SCD); 
    }

    SM_explorer SMexplorer(Vertex_const_handle v) const { 
      SM_const_decorator SMCD(v);
      return SM_explorer(SMCD); 
    }

  //typedef CGAL::SNC_explorer<SNC_decorator,T> SNC_explorer; //not implemented
  /*{\Mtypemember a decorator to examine the underlying plane map. 
  See the manual page of |EW_explorer|.}*/

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
    { CGAL_nef3_assertion_msg( 0, "not implemented."); }

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

  //EW_explorer explorer() const // not implemented, replaced by decorator?
  /*{\Mop returns a decorator object which allows read-only access of
  the underlying three-dimensional subdivision structure. 
  See the manual page |EW_explorer| for its usage.}*/
  //{ return EW_explorer(snc,EK); }

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

template <typename T>
T Nef_polyhedron_3<T>::EK = T();

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3( Content space, SNC_point_locator* _pl) {
  TRACEN("construction from empty or space.");
  pl() = _pl;
  if(Infi_box::extended_Kernel()) {
    SNC_constructor C(snc());
    Infi_box::initialize_infibox_vertices(C,space == COMPLETE);
    //  initialize_simple_cube_vertices(space);
    build_external_structure();
  }
}

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3(const Plane_3& h, Boundary b, SNC_point_locator* _pl) {
  TRACEN("construction from plane "<<h);
  SNC_constructor C(snc(), pl());
  Infi_box::create_vertices_of_box_with_plane(C,h,(b==INCLUDED));
  pl() = _pl;
  build_external_structure();
}

template <typename T>
Nef_polyhedron_3<T>::
Nef_polyhedron_3( const SNC_structure& W, SNC_point_locator* _pl, bool clone) {
  CGAL_nef3_assertion( clone == true);
  // TODO: granados: define behavior when clone=false
  TRACEN("construction from an existing SNC structure (clone="<<clone<<")");
  if (clone) { 
    snc() = W;
    pl() = _pl->clone();
    pl()->initialize(&snc());
  }
}

template <typename T>
void Nef_polyhedron_3<T>::extract_complement() {
  TRACEN("extract complement");
  if( is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_nef3_forall_vertices(v,D){
    D.mark(v) = !D.mark(v); 
    SM_decorator SM(v);
    SM.extract_complement();
    SM.mark_of_halfsphere(-1) = !SM.mark_of_halfsphere(-1);    
    SM.mark_of_halfsphere(+1) = !SM.mark_of_halfsphere(+1);    
  }
  Halffacet_iterator f;
  CGAL_nef3_forall_facets(f,D) D.mark(f) = !D.mark(f); 
 
  Volume_iterator c;
  CGAL_nef3_forall_volumes(c,D) 
    //    if(!(Infi_box::extended_Kernel && c==D.volumes_begin()))
      D.mark(c) = !D.mark(c);
}


template <typename T>
void Nef_polyhedron_3<T>::extract_interior() {
  TRACEN("extract interior");
  if (is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_nef3_forall_vertices(v,D){
    D.mark(v) = false;
    SM_decorator SM(v);
    SM.extract_interior();
  }
  Halffacet_iterator f;
  CGAL_nef3_forall_facets(f,D) D.mark(f) = false;

  simplify();
}

template <typename T>
void Nef_polyhedron_3<T>::extract_boundary() {
  TRACEN("extract boundary");
  if (is_shared()) clone_rep();
  SNC_decorator D(snc());
  Vertex_iterator v;
  CGAL_nef3_forall_vertices(v,D) {
    D.mark(v) = true;
    SM_decorator SM(v);
    SM.extract_boundary();
    SM.mark_of_halfsphere(-1) = SM.mark_of_halfsphere(+1) = false;
  }
  Halffacet_iterator f;
  CGAL_nef3_forall_facets(f,D) D.mark(f) = true;
  Volume_iterator c;
  CGAL_nef3_forall_volumes(c,D) D.mark(c) = false;
  simplify();
}

template <typename T>
std::ostream& operator<<
 (std::ostream& os, Nef_polyhedron_3<T>& NP)
{
  typedef typename Nef_polyhedron_3<T>::SNC_structure SNC_structure;
  CGAL::SNC_io_parser<SNC_structure> O(os, NP.snc());
  O.print();
  return os;
}

template <typename T>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_3<T>& NP)
{
  typedef typename Nef_polyhedron_3<T>::SNC_decorator SNC_decorator;
  CGAL::SNC_io_parser<SNC_decorator> I(is, NP.snc());
  I.read();
  I.check_integrity();
  return is;
}

CGAL_END_NAMESPACE

#endif //CGAL_NEF_POLYHEDRON_3_H


