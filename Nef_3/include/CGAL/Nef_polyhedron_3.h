// Copyright (c) 1997-2002,2005 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
//                 Ralf Osbild     <osbild@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NEF_POLYHEDRON_3_H
#define CGAL_NEF_POLYHEDRON_3_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Nef_3/Default_items.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_external_structure.h>
#include <CGAL/Nef_3/Combine_with_halfspace.h>
#ifdef CGAL_NEF_VISUAL_HULL
#include <CGAL/Nef_3/Binary_operation_vs.h>
#else
#include <CGAL/Nef_3/Binary_operation.h>
#endif
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_SM_explorer.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Nef_3/Mark_bounded_volumes.h>

#ifdef CGAL_NEF3_POINT_LOCATOR_NAIVE
#include <CGAL/Nef_3/SNC_ray_shooter.h>
#endif

#ifdef CGAL_NEF3_CGAL_NEF3_SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // CGAL_NEF3_SM_VISUALIZOR

#ifdef CGAL_NEF3_OLD_VISUALIZATION 
#include <CGAL/Nef_3/Visualizor_OpenGL_3.h>
#endif // CGAL_NEF3_OLD_VISUALIZATION 

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/polyhedron_3_to_nef_3.h>
#include <CGAL/Nef_3/shell_to_nef_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/assertions.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <list>

// RO: includes for "vertex cycle to Nef" constructor
#include <CGAL/Nef_3/vertex_cycle_to_nef_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/normal_vector_newell_3.h>

#ifdef CGAL_NEF_VISUAL_HULL
#include <CGAL/Nef_3/Modifying_binary_operation_vs.h>
#endif

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 11
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

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
  typedef CGAL::SNC_structure<K,I,M>                      SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;
  typedef CGAL::SNC_const_decorator<SNC_structure>        SNC_const_decorator;
  typedef CGAL::Binary_operation<SNC_structure>           Binary_operation;
  typedef CGAL::SNC_constructor<I, SNC_structure>         SNC_constructor;
  typedef CGAL::SNC_external_structure<I, SNC_structure>  SNC_external_structure;
  typedef CGAL::SNC_point_locator<SNC_decorator> SNC_point_locator;
  typedef CGAL::SNC_simplify<I, SNC_structure>            SNC_simplify;
#ifdef CGAL_NEF3_POINT_LOCATOR_NAIVE
  typedef CGAL::SNC_point_locator_naive<SNC_decorator> SNC_point_locator_default;
#else
  typedef CGAL::SNC_point_locator_by_spatial_subdivision<SNC_decorator> SNC_point_locator_default;
#endif

  typedef typename SNC_structure::Sphere_map       Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>           SM_decorator;
  typedef CGAL::SM_const_decorator<Sphere_map>     SM_const_decorator;
  typedef CGAL::SNC_SM_overlayer<I, SM_decorator>  SM_overlayer;
  typedef CGAL::SM_point_locator<SNC_structure>    SM_point_locator;

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

template <typename Kernel_, typename Items_ = typename CGAL::Default_items<Kernel_>::Items, typename Mark_ = bool>
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

#ifndef _MSC_VER
  // VC++ has a problem to digest the following typedef,
  // and does not need the using statements -- AF
  // The left and right part of these typedefs have the same name. It is
  // very important to qualify the left part with the CGAL:: namespace, no
  // to confuse g++. -- Laurent Rineau, 2010/09/13
  typedef CGAL::SNC_structure<Kernel,Items,Mark> SNC_structure;
  typedef CGAL::SNC_const_decorator<SNC_structure> SNC_const_decorator;
  using SNC_const_decorator::set_snc;
  using SNC_const_decorator::is_standard;
  using SNC_const_decorator::is_bounded;
#endif

  struct Polylines_tag {};

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  typedef enum { EMPTY=0, COMPLETE=1 } Content;
  /*{\Menum construction selection}*/

  typedef enum { DEFAULT, NAIVE, WALKING, SPATIAL_SUBDIVISION  } Location_mode;
  /*{\Menum selection flag for the point location mode.}*/

protected: 
  struct AND { Mark operator()(const Mark& b1, const Mark& b2, bool /* inverted */ =false) const { return b1&&b2; } };
  struct OR { Mark operator()(const Mark& b1, const Mark& b2, bool /* inverted */ =false) const { return b1||b2; } };
  struct DIFF { Mark operator()(const Mark& b1, const Mark& b2, bool inverted=false) const { 
    if(inverted) return !b1&&b2; return b1&&!b2; } };
  struct XOR { Mark operator()(const Mark& b1, const Mark& b2, bool /* inverted */ =false) const 
    { return (b1&&!b2)||(!b1&&b2); } };

 public:
  typedef Nef_polyhedron_3_rep<Kernel,Items, Mark>    Nef_rep;

  typedef typename Nef_rep::SM_decorator        SM_decorator;
  typedef typename Nef_rep::SM_const_decorator  SM_const_decorator;
 protected:
  typedef typename Nef_rep::SNC_decorator       SNC_decorator;
  typedef typename Nef_rep::SNC_constructor     SNC_constructor;
  typedef typename Nef_rep::SNC_external_structure SNC_external_structure;
  typedef typename Nef_rep::Binary_operation    Binary_operation;
  typedef typename Nef_rep::SNC_point_locator   SNC_point_locator;
  typedef typename Nef_rep::SNC_point_locator_default
    SNC_point_locator_default;
  typedef CGAL::Combine_with_halfspace<SNC_structure, SNC_point_locator> 
          Combine_with_halfspace;
public:
 enum Intersection_mode { 
	 CLOSED_HALFSPACE = Combine_with_halfspace::CLOSED_HALFSPACE, 
     OPEN_HALFSPACE = Combine_with_halfspace::OPEN_HALFSPACE, 
     PLANE_ONLY = Combine_with_halfspace::PLANE_ONLY};

protected: 
  typedef typename Nef_rep::SM_overlayer        SM_overlayer;
  typedef typename Nef_rep::SM_point_locator    SM_point_locator;
  typedef typename Nef_rep::SNC_simplify        SNC_simplify;
#ifdef CGAL_NEF3_SM_VISUALIZOR
  typedef typename Nef_rep::SM_visualizor       SM_visualizor;
#endif // CGAL_NEF3_SM_VISUALIZOR
#ifdef CGAL_NEF3_OLD_VISUALIZATION 
  typedef CGAL::Nef_Visualizor_OpenGL_3<Nef_polyhedron_3> Visualizor;
#endif // CGAL_NEF3_OLD_VISUALIZATION 

 typedef typename Nef_rep::Sphere_map                Sphere_map;
 public:
 typedef CGAL::Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map> Nef_polyhedron_S2;
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
  typedef typename SNC_structure::Halffacet_cycle_iterator     
                                  Halffacet_cycle_iterator;
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
  typedef typename SM_decorator::SHalfloop_iterator
                                                   SHalfloop_iterator;
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

  typedef typename SNC_decorator::Association  Association;


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
    SNC_external_structure es(snc(), pl());
    es.build_external_structure();
  }

 public:
  /*{\Mcreation 3}*/

  Nef_polyhedron_3( Content space = EMPTY);
		   
  /*{\Mcreate creates an instance |\Mvar| of type |\Mname|
  and initializes it to the empty set if |space == EMPTY|
  and to the whole space if |space == COMPLETE|.}*/

  explicit Nef_polyhedron_3(const Plane_3& p, Boundary b = INCLUDED);
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the
  halfspace on the negative side of |p| including |p| if |b==INCLUDED|,
  excluding |p| if |b==EXCLUDED|.}*/

  Nef_polyhedron_3(const Nef_polyhedron_3<Kernel,Items, Mark>& N1) 
 : Base(N1) , SNC_const_decorator() {
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
         SNC_external_structure es(snc(), &Pl);
         es.build_external_structure();
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
  
 template<typename Items, typename SNC_structure>
 class Sphere_map_creator {
   typedef typename SNC_structure::SM_decorator     SM_decorator;
   typedef typename SNC_structure::Vertex_handle    Vertex_handle;
   typedef typename SNC_structure::SVertex_handle   SVertex_handle;
   typedef typename SNC_structure::SFace_handle     SFace_handle;
   typedef typename SNC_structure::Sphere_point     Sphere_point;
   
   public:
   Sphere_map_creator() {}
   
   template<typename point_iterator>
   void create_end_sphere_map(SNC_structure& snc,
			      point_iterator cur,
			      point_iterator prev) {
     Vertex_handle v(snc.new_vertex(*cur, true));
     SM_decorator SM(&*v);
     SVertex_handle sv(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
				      true));
     SFace_handle sf(v->new_sface());
     SM.link_as_isolated_vertex(sv,sf);
   }
   
   template<typename point_iterator>
   void create_sphere_map(SNC_structure& snc,
			  point_iterator cur,
			  point_iterator prev,
			  point_iterator next) {
     Vertex_handle v(snc.new_vertex(*cur, true));
     SM_decorator SM(&*v);
     SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
				       true));
     SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
				       true));      
     SFace_handle sf(v->new_sface());
     SM.link_as_isolated_vertex(sv1,sf);
     SM.link_as_isolated_vertex(sv2,sf);
   }
 };
 
 template<typename SNC_structure>
 class Sphere_map_creator<CGAL::SNC_indexed_items, SNC_structure> {
   typedef typename SNC_structure::SM_decorator     SM_decorator;
   typedef typename SNC_structure::Vertex_handle    Vertex_handle;
   typedef typename SNC_structure::SVertex_handle   SVertex_handle;
   typedef typename SNC_structure::SFace_handle     SFace_handle;
   typedef typename SNC_structure::Sphere_point     Sphere_point;
   
   bool first;
   int index;
 public:
   Sphere_map_creator() : first(true), index(0) {}
     
     template<typename point_iterator>
       void create_end_sphere_map(SNC_structure& snc,
				  point_iterator cur,
				  point_iterator prev) {
       Vertex_handle v(snc.new_vertex(*cur, true));
       SM_decorator SM(&*v);
       SVertex_handle sv(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
					true));
       SFace_handle sf(v->new_sface());
       SM.link_as_isolated_vertex(sv,sf);
       if(first) {
	 sv->set_index();
	 index = sv->get_index();
	 first = false;
       } else
	 sv->set_index(index);
     }
     
     template<typename point_iterator>
       void create_sphere_map(SNC_structure& snc,
			      point_iterator cur,
			      point_iterator prev,
			      point_iterator next) {
       Vertex_handle v(snc.new_vertex(*cur, true));
       SM_decorator SM(&*v);
       SVertex_handle sv1(v->new_svertex(Sphere_point(ORIGIN+(*prev-*cur)),
					 true));
       SVertex_handle sv2(v->new_svertex(Sphere_point(ORIGIN+(*next-*cur)),
					 true));      
       SFace_handle sf(v->new_sface());
       SM.link_as_isolated_vertex(sv1,sf);
       SM.link_as_isolated_vertex(sv2,sf);
       sv1->set_index(index);
       sv2->set_index();
       index = sv2->get_index();
     }
 };
 
 template <typename InputIterator>
 Nef_polyhedron_3(InputIterator begin, InputIterator end, Polylines_tag) {
   typedef typename std::iterator_traits<InputIterator>::value_type
     point_iterator_pair;
   typedef typename point_iterator_pair::first_type
     point_iterator;

   empty_rep();
   set_snc(snc());
   initialize_infibox_vertices(EMPTY);

   point_iterator pbegin, pend, pnext, pprev;
   Sphere_map_creator<Items, SNC_structure> smc;
   for(;begin != end; ++begin) {
     pend = begin->second;
     pprev = pnext = pbegin = begin->first;
     ++pnext;
     CGAL_assertion(pnext != pend);
     smc.create_end_sphere_map(snc(),pbegin,pnext);
     for(++pbegin,++pnext; pnext!=pend; ++pbegin,++pprev,++pnext)
       smc.create_sphere_map(snc(),pbegin,pprev,pnext);
     smc.create_end_sphere_map(snc(),pbegin,pprev);
   }
   build_external_structure();
   simplify();
 }

 template <class T1, class T2,
           template <class T31, class T32, class T33>
           class T3, class T4 >
 Nef_polyhedron_3( CGAL::Polyhedron_3<T1,T2,T3,T4>& P) {
    CGAL_NEF_TRACEN("construction from Polyhedron_3");
    SNC_structure rsnc;
    *this = Nef_polyhedron_3(rsnc, new SNC_point_locator_default, false);
    initialize_infibox_vertices(EMPTY);
    polyhedron_3_to_nef_3
      <CGAL::Polyhedron_3<T1,T2,T3,T4>, SNC_structure>( P, snc());
    build_external_structure();
    simplify();
    CGAL::Mark_bounded_volumes<Nef_polyhedron_3> mbv(true);
    delegate(mbv);
    set_snc(snc());
  }
  
 Nef_polyhedron_3(const Nef_polyhedron& N, 
		  SFace_const_iterator sf) 
 {
   SNC_structure rsnc;
   *this = Nef_polyhedron_3(rsnc, new SNC_point_locator_default, false);
   initialize_infibox_vertices(EMPTY);
   shell_to_nef_3(N, sf, snc());
   build_external_structure();
   simplify();
   CGAL::Mark_bounded_volumes<Nef_polyhedron_3> mbv(true);
   delegate(mbv);
   set_snc(snc());
 }


 protected: 

  template<typename Kernel>
  class Triangulation_handler2 {

    typedef typename CGAL::Triangulation_vertex_base_2<Kernel>               Vb;
    typedef typename CGAL::Constrained_triangulation_face_base_2<Kernel>     Fb;
    typedef typename CGAL::Triangulation_data_structure_2<Vb,Fb>             TDS;
    typedef typename CGAL::Constrained_triangulation_2<Kernel,TDS>           CT;

    typedef typename CT::Face_handle           Face_handle;
    typedef typename CT::Vertex_handle         CTVertex_handle;
    typedef typename CT::Finite_faces_iterator Finite_face_iterator;
    typedef typename CT::Edge                  Edge;

    CT ct;
    CGAL::Unique_hash_map<Face_handle, bool> visited;
    CGAL::Unique_hash_map<CTVertex_handle, Vertex_const_handle> ctv2v;
    Finite_face_iterator fi;
    Plane_3 supporting_plane;

  public:
    Triangulation_handler2(Halffacet_const_handle f) : 
      visited(false), supporting_plane(f->plane()) {

      Halffacet_cycle_const_iterator fci;
      for(fci=f->facet_cycles_begin(); fci!=f->facet_cycles_end(); ++fci) {
	if(fci.is_shalfedge()) {
          SHalfedge_around_facet_const_circulator sfc(fci), send(sfc);
	  CGAL_For_all(sfc,send) {
            CGAL_NEF_TRACEN("  insert point" << sfc->source()->source()->point());
	    CTVertex_handle ctv = ct.insert(sfc->source()->source()->point());
	    ctv2v[ctv] = sfc->source()->source();
          }
        }
      }

      for(fci=f->facet_cycles_begin(); fci!=f->facet_cycles_end(); ++fci) {
	if(fci.is_shalfedge()) {
          SHalfedge_around_facet_const_circulator sfc(fci), send(sfc);
	  CGAL_For_all(sfc,send) {
            CGAL_NEF_TRACEN("  insert constraint" << sfc->source()->source()->point()
	                     << "->" << sfc->source()->twin()->source()->point());
	    ct.insert_constraint(sfc->source()->source()->point(),
	                         sfc->source()->twin()->source()->point());
          }
        }
      }
      CGAL_assertion(ct.is_valid());

      CGAL_NEF_TRACEN("number of finite triangles " << ct.number_of_faces());

      typename CT::Face_handle infinite = ct.infinite_face();
      typename CT::Vertex_handle ctv = infinite->vertex(1);
      if(ct.is_infinite(ctv)) ctv = infinite->vertex(2);
      CGAL_assertion(!ct.is_infinite(ctv));

      typename CT::Face_handle opposite;
      typename CT::Face_circulator vc(ctv,infinite);
      do { opposite = vc++;
      } while(!ct.is_constrained(typename CT::Edge(vc,vc->index(opposite))));
      typename CT::Face_handle first = vc;

      CGAL_assertion(!ct.is_infinite(first));
      traverse_triangulation(first, first->index(opposite));

      fi = ct.finite_faces_begin();
    }

    void traverse_triangulation(Face_handle f, int parent) {
      visited[f] = true;
      if(!ct.is_constrained(Edge(f,ct.cw(parent))) && !visited[f->neighbor(ct.cw(parent))]) {
	Face_handle child(f->neighbor(ct.cw(parent)));
	traverse_triangulation(child, child->index(f));
      } 
      if(!ct.is_constrained(Edge(f,ct.ccw(parent))) && !visited[f->neighbor(ct.ccw(parent))]) {
	Face_handle child(f->neighbor(ct.ccw(parent)));
	traverse_triangulation(child, child->index(f));
      } 
    } 
 
    template<typename Triangle_3>
    bool get_next_triangle(Triangle_3& tr) {
      while(fi != ct.finite_faces_end() && visited[fi] == false) ++fi;
      if(fi == ct.finite_faces_end()) return false;
      tr = Triangle_3(fi->vertex(0)->point(), fi->vertex(1)->point(), fi->vertex(2)->point());
      ++fi;
      return true;
    }

    bool same_orientation(Plane_3 p1) const {
      if(p1.a() != 0)
	return CGAL::sign(p1.a()) == CGAL::sign(supporting_plane.a());
      if(p1.b() != 0)
	return CGAL::sign(p1.b()) == CGAL::sign(supporting_plane.b());
      return CGAL::sign(p1.c()) == CGAL::sign(supporting_plane.c());
    }

    template<typename PIB, typename Index>
    void handle_triangles(PIB& pib, Index& VI) {
      while(fi != ct.finite_faces_end() && visited[fi] == false) ++fi;
      while(fi != ct.finite_faces_end()) {
	Plane_3 plane(fi->vertex(0)->point(),
		      fi->vertex(1)->point(),
		      fi->vertex(2)->point());
	pib.begin_facet();
	if(same_orientation(plane)) {
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(0)]]);
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(1)]]);
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(2)]]);
	} else {
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(0)]]);
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(2)]]);
	  pib.add_vertex_to_facet(VI[ctv2v[fi->vertex(1)]]);
	}
	pib.end_facet();
	do {
	  ++fi;
	} while(fi != ct.finite_faces_end() && visited[fi] == false);
      }
    }
  };
 
  template <class HDS>
  class Build_polyhedron : public CGAL::Modifier_base<HDS> {
    
    class Visitor {

      typedef typename CGAL::Projection_traits_xy_3<Kernel>       XY;
      typedef typename CGAL::Projection_traits_yz_3<Kernel>       YZ;
      typedef typename CGAL::Projection_traits_xz_3<Kernel>       XZ;

      const Object_index<Vertex_const_iterator>& VI;
      Polyhedron_incremental_builder_3<HDS>& B;
      const SNC_const_decorator& D;
      
    public:
      Visitor(Polyhedron_incremental_builder_3<HDS>& BB,
	      const SNC_const_decorator& sd,
	      Object_index<Vertex_const_iterator>& vi) : VI(vi), B(BB), D(sd){}

      void visit(Halffacet_const_handle opposite_facet) {

	CGAL_NEF_TRACEN("Build_polyhedron: visit facet " << opposite_facet->plane());
 
	CGAL_assertion(Infi_box::is_standard(opposite_facet->plane()));
	
	SHalfedge_const_handle se;
	Halffacet_cycle_const_iterator fc;
     	
	Halffacet_const_handle f = opposite_facet->twin();

	SHalfedge_around_facet_const_circulator 
	  sfc1(f->facet_cycles_begin()), sfc2(sfc1);
	
	if(++f->facet_cycles_begin() != f->facet_cycles_end() ||
	   ++(++(++sfc1)) != sfc2) {
	  Vector_3 orth = f->plane().orthogonal_vector();
	  int c = CGAL::abs(orth[0]) > CGAL::abs(orth[1]) ? 0 : 1;
	  c = CGAL::abs(orth[2]) > CGAL::abs(orth[c]) ? 2 : c;
	  
	  if(c == 0) {
	    Triangulation_handler2<YZ> th(f);
	    th.handle_triangles(B, VI);
	  } else if(c == 1) {
	    Triangulation_handler2<XZ> th(f);
	    th.handle_triangles(B, VI);
	  } else if(c == 2) {
	    Triangulation_handler2<XY> th(f);
	    th.handle_triangles(B, VI);
	  } else
	    CGAL_error_msg( "wrong value");

	} else {

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
      }

      void visit(SFace_const_handle) {}
      void visit(Halfedge_const_handle) {}
      void visit(Vertex_const_handle) {}
      void visit(SHalfedge_const_handle) {}
      void visit(SHalfloop_const_handle) {}
    };

  public:

    const SNC_const_decorator& scd;
    Object_index<Vertex_const_iterator> VI;

    Build_polyhedron(const SNC_const_decorator& s) : 
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
			2*scd.number_of_vertices()-4,
			3*scd.number_of_vertices()-6);
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
      if ( B.error() ) B.rollback();
    }

  };

  template <class HDS>
  class Build_polyhedron2 : public CGAL::Modifier_base<HDS> {
    
    class Find_holes {

      Unique_hash_map<Vertex_const_handle, bool>& omit_vertex;
      int nov, nof;

    public:
      Find_holes(Unique_hash_map<Vertex_const_handle, bool>& omit_vertex_) 
	: omit_vertex(omit_vertex_), nov(0), nof(0) {}

      void visit(Halffacet_const_handle f) {
	++nof;
	Halffacet_cycle_const_iterator fc = f->facet_cycles_begin();
	for(++fc; fc != f->facet_cycles_end(); ++fc) {
	  if(fc.is_shalfedge()) {
	    --nof;
	    SHalfedge_around_facet_const_circulator 
	      sfc(fc), send(sfc);
	    CGAL_For_all(sfc, send) {
	      omit_vertex[sfc->source()->source()] = true;
	      --nov;
	    }
	  } else if(fc.is_shalfloop()) {
	    SHalfloop_const_handle sl(fc);
	    omit_vertex[sl->incident_sface()->center_vertex()];
	    --nov;
	  } else
	    CGAL_error_msg( "wrong handle type");
	}
      }

      void visit(Vertex_const_handle) { ++nov; }
      void visit(SFace_const_handle) {}
      void visit(Halfedge_const_handle) {}
      void visit(SHalfedge_const_handle) {}
      void visit(SHalfloop_const_handle) {}

      int number_of_vertices() const {
	return nov;
      }

      int number_of_facets() const {
	return nof;
      }
    };

    class Add_vertices {
      
      Polyhedron_incremental_builder_3<HDS>& B;
      Unique_hash_map<Vertex_const_handle, bool>& omit_vertex;
      Object_index<Vertex_const_iterator>& VI;      
      int vertex_index;

    public:
      Add_vertices(Polyhedron_incremental_builder_3<HDS>& B_,
		   Unique_hash_map<Vertex_const_handle, bool>& omit_vertex_,
		   Object_index<Vertex_const_iterator>& VI_) 
	: B(B_), omit_vertex(omit_vertex_), VI(VI_), vertex_index(0) {}
	
      void visit(Vertex_const_handle v) {
	if(omit_vertex[v]) return;
	VI[v]=vertex_index++;
	B.add_vertex(v->point());
      }

      void visit(Halffacet_const_handle) {}
      void visit(SFace_const_handle) {}
      void visit(Halfedge_const_handle) {}
      void visit(SHalfedge_const_handle) {}
      void visit(SHalfloop_const_handle) {}

    };

    class Visitor {

      const Object_index<Vertex_const_iterator>& VI;
      Polyhedron_incremental_builder_3<HDS>& B;
      const Unique_hash_map<Vertex_const_handle, bool>& omit_vertex;
      SNC_const_decorator& D;
      
    public:
      Visitor(Polyhedron_incremental_builder_3<HDS>& BB,
	      const Unique_hash_map<Vertex_const_handle, bool>& omit_vertex_,
	      SNC_const_decorator& sd,
	      Object_index<Vertex_const_iterator>& vi) 
	: VI(vi), B(BB), omit_vertex(omit_vertex_), D(sd){}

      void visit(Halffacet_const_handle opposite_facet) {

	CGAL_NEF_TRACEN("Build_polyhedron: visit facet " << opposite_facet->plane());
 
	CGAL_assertion(Infi_box::is_standard(opposite_facet->plane()));
	
	SHalfedge_const_handle se;
	Halffacet_cycle_const_iterator fc;
     	
	Halffacet_const_handle f = opposite_facet->twin();

	fc = f->facet_cycles_begin();
	se = SHalfedge_const_handle(fc);
	CGAL_assertion(se!=0);
	if(omit_vertex[se->source()->source()]) return;
	B.begin_facet();
	SHalfedge_around_facet_const_circulator hc_start(se);
	SHalfedge_around_facet_const_circulator hc_end(hc_start);
	CGAL_For_all(hc_start,hc_end) {
	  CGAL_NEF_TRACEN("   add vertex " << hc_start->source()->center_vertex()->point());
	  B.add_vertex_to_facet(VI[hc_start->source()->center_vertex()]);
	}
	B.end_facet();
      }

      void visit(SFace_const_handle) {}
      void visit(Halfedge_const_handle) {}
      void visit(Vertex_const_handle) {}
      void visit(SHalfedge_const_handle) {}
      void visit(SHalfloop_const_handle) {}
    };

  public:

    SFace_const_handle sf;
    SNC_const_decorator& scd;
    Object_index<Vertex_const_iterator> VI;
    Unique_hash_map<Vertex_const_handle, bool> omit_vertex;

    Build_polyhedron2(SFace_const_handle sf_, SNC_const_decorator& s) : 
      sf(sf_), scd(s), VI(s.vertices_begin(),s.vertices_end(),'V'), 
      omit_vertex(false) {}
    
      void operator()(HDS& hds) {

      Polyhedron_incremental_builder_3<HDS> B(hds, true);
      
      Find_holes F(omit_vertex);
      scd.visit_shell_objects(sf, F);

      B.begin_surface(F.number_of_vertices(), 
		      F.number_of_facets(),
		      F.number_of_vertices()+F.number_of_facets()-2);

      Add_vertices A(B,omit_vertex, VI);
      scd.visit_shell_objects(sf, A);

      Visitor V(B,omit_vertex, scd,VI);
      scd.visit_shell_objects(sf, V);
      B.end_surface();
      if ( B.error() ) B.rollback();
    }

  };


 public:
 void delegate( Modifier_base<SNC_structure>& modifier, 
		bool compute_external = false, 
		bool do_simplify = true) {

   // calls the `operator()' of the `modifier'. Precondition: The
   // `modifier' returns a consistent representation.
   if( this->is_shared()) clone_rep();
   modifier(snc());
   if(compute_external) {
     SNC_external_structure es(snc());
     es.clear_external_structure();
     
     build_external_structure();
   }
   if(do_simplify)
     simplify();
   CGAL_expensive_postcondition( is_valid());
 }

 struct SNC_and_PL {
   SNC_structure* sncp;
   SNC_point_locator* pl;

   SNC_and_PL(SNC_structure* s, SNC_point_locator* p) : sncp(s), pl(p) {}
 };

 void delegate( Modifier_base<SNC_and_PL>& modifier, 
		bool compute_external = false,
		bool do_simplify = false) {
   // calls the `operator()' of the `modifier'. Precondition: The
   // `modifier' returns a consistent representation.
   if( this->is_shared()) clone_rep();
   SNC_and_PL sncpl(&snc(),pl());
   modifier(sncpl);
   pl() = sncpl.pl;
   if(compute_external) {
     SNC_external_structure es(snc());
     es.clear_external_structure();
     build_external_structure();
   }
   if(do_simplify)
     simplify();
   CGAL_expensive_postcondition( is_valid());
 }
 
 public:

 template<typename Polyhedron>
 void convert_to_Polyhedron(Polyhedron& P) const {
   convert_to_polyhedron(P);
 }

 template<typename Polyhedron>
 void convert_to_polyhedron(Polyhedron& P) const {
   typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
   CGAL_precondition(is_simple());
   P.clear();
   Build_polyhedron<HalfedgeDS> bp(*this);    
   P.delegate(bp);
 }

 template<typename Polyhedron>
 void convert_inner_shell_to_polyhedron(SFace_const_iterator sf, Polyhedron& P) {
   typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
   P.clear();
   Build_polyhedron2<HalfedgeDS> bp(sf, *this);
   P.delegate(bp);
 }

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

  bool is_simple() const {

    Halfedge_const_iterator e;
    CGAL_forall_edges(e,snc())
      if(!is_edge_2manifold(e))
	return false;

    CGAL_NEF_TRACEN("there is no edge with non-manifold situation");

    Vertex_const_iterator v;
    CGAL_forall_vertices(v,snc())
      if(!is_vertex_2manifold(v))
	return false;

    CGAL_NEF_TRACEN("there is no vertex with non-manifold situation");
/*
    Halffacet_iterator f;
    CGAL_forall_halffacets(f,snc())
      if(!is_facet_simple(f))
	return false;

    CGAL_NEF_TRACEN("there are no holes");
*/
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
  bool is_edge_2manifold(const Halfedge_const_handle& e) const {

    SM_decorator SD;
    SHalfedge_around_svertex_const_circulator c(SD.first_out_edge(e)), c2(c);

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
 
  bool is_vertex_2manifold(const Vertex_const_handle& v) const {
     
    SFace_const_iterator sfi(v->sfaces_begin());
    if (++sfi != v->sfaces_last())
      return false;

    return true;
  }

  bool is_facet_simple(const Halffacet_const_handle& f) const {
    
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
#ifdef CGAL_NEF3_OLD_VISUALIZATION   
  void visualize() { 
    Visualizor sncv(*this);
    sncv.draw();
    //OGL::polyhedra_.back().debug();
    OLDOGL::start_viewer();
  }
#endif // CGAL_NEF3_OLD_VISUALIZATION   
   
  void clear(Content space = EMPTY)
    { *this = Nef_polyhedron_3(space); }
  /*{\Mop makes |\Mvar| the empty set if |space == EMPTY| and the
  full space if |space == COMPLETE|.}*/

 bool is_empty() const {
   /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
   if(Infi_box::extended_kernel()) 
     return this->number_of_vertices() == 8 &&
            this->number_of_edges() == 12 &&
            this->number_of_facets() == 6 &&
            this->number_of_volumes() == 2 &&
            (++this->volumes_begin())->mark() == false;

   else 
     return this->number_of_vertices() == 0 &&
            this->number_of_edges() == 0 &&
            this->number_of_facets() == 0 &&
            this->number_of_volumes() == 1 &&
            (this->volumes_begin())->mark() == false;
  }

 bool is_space() const {
  /*{\Mop returns true if |\Mvar| is the whole space, false otherwise.}*/
   if(Infi_box::extended_kernel()) 
     return this->number_of_vertices() == 8 &&
            this->number_of_edges() == 12 &&
            this->number_of_facets() == 6 &&
            this->number_of_volumes() == 2 &&
            (++this->volumes_begin())->mark() == true;

   else 
     return this->number_of_vertices() == 0 &&
            this->number_of_edges() == 0 &&
            this->number_of_facets() == 0 &&
            this->number_of_volumes() == 1 &&
            (this->volumes_begin())->mark() == true;
  }

  /*{\Xtext \headerline{Destructive Operations}}*/

 protected:
  void clone_rep() { *this = Nef_polyhedron_3<Kernel,Items, Mark>(snc(), pl()); }
  void empty_rep() { 
    SNC_structure rsnc;
    delete pl();
    *this = Nef_polyhedron_3<Kernel,Items, Mark>(rsnc, new SNC_point_locator_default,false);
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
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    Binary_operation bo( res.snc());
    bo(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _and);
    return res;
  }

  Nef_polyhedron_3<Kernel,Items, Mark>
   intersection(const Plane_3& plane, 
		Intersection_mode im) const {
    AND _and;
    SNC_structure rsnc;
    Nef_polyhedron_3<Kernel,Items, Mark> res(rsnc, new SNC_point_locator_default, false);
    Combine_with_halfspace cwh(res.snc(), res.pl());
    cwh.combine_with_halfspace(snc(), plane, _and, 
							   static_cast<typename Combine_with_halfspace::Intersection_mode>(im));
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
    Binary_operation bo(res.snc());
    bo(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _or);
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
    Binary_operation bo(res.snc());
    bo(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _diff);
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
    Binary_operation bo(res.snc());
    bo(res.pl(), snc(), pl(), N1.snc(), N1.pl(), _xor);
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
  { 
#ifdef CGAL_NEF_VISUAL_HULL
    CGAL_NEF_SETDTHREAD(19*43*71);
    std::cerr << "visual hull code " << std::endl;
    std::cerr << *this << std::endl;
    std::cerr << const_cast<Nef_polyhedron&>(N1) << std::endl;
    AND _and;
    typename CGAL::Modifying_binary_operation<SNC_structure> 
      mbo(this->snc());
    mbo(const_cast<SNC_structure&>(N1.snc()), N1.pl(), pl(), _and);
    return *this;
#else
    *this = intersection(N1); return *this; 
#endif
  }

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

 void change_orientation(bool full = false) {

   Halffacet_handle ftemp;
   Volume_handle vtemp;
   SVertex_handle svtemp;
   SHalfedge_handle setemp;
   SFace_handle sftemp;

   SVertex_iterator sv;
   CGAL_forall_svertices(sv, snc()) {
     sv->out_sedge() = sv->out_sedge()->twin();
   }

   SHalfedge_iterator se;
   CGAL_forall_shalfedges(se, snc()) {
     if(se->is_twin()) {
       svtemp = se->source();
       se->source() = se->twin()->source();
       se->twin()->source() = svtemp;

       if(full) {
	 ftemp = se->facet();
	 se->facet() = se->twin()->facet();
	 se->twin()->facet() = ftemp;
       }
       //       sftemp = se->incident_sface();
       //       se->incident_sface() = se->twin()->incident_sface();
       //       se->twin()->incident_sface() = sftemp;
     }
       
     setemp = se->sprev();
     se->sprev() = se->snext();
     se->snext() = setemp;

     se->circle() = se->circle().opposite();

     if(full) {
       setemp = se->prev();
       se->prev() = se->next();
       se->next() = setemp;
     }
   }

   if(full) {
     Halffacet_iterator f;
     CGAL_forall_facets(f, snc()) {
       vtemp  = f->incident_volume();
       f->incident_volume() = f->twin()->incident_volume();
       f->twin()->incident_volume() = vtemp;
       Halffacet_cycle_iterator fc(f->facet_cycles_begin()),
	 fct(f->twin()->facet_cycles_begin());
       while(fc!=f->facet_cycles_end()) {
	 CGAL_assertion(fct!=f->twin()->facet_cycles_end());
	 if(fc.is_shalfedge()) {
	   CGAL_assertion(fct.is_shalfedge());
	   setemp = fc;
	   *fc = *fct;
	   *fct = make_object(setemp);
	 }
	 ++fc;
	 ++fct;
       }
     }
   
     CGAL_forall_halffacets(f, snc()) {
       Halffacet_cycle_iterator fc(f->facet_cycles_begin());
       for(;fc!=f->facet_cycles_end();++fc) {
	 if(fc.is_shalfedge()) {
	   setemp = fc;
	   SHalfedge_around_facet_circulator hfc(setemp),hend(hfc);
	   ++hfc;
	   CGAL_For_all(hfc,hend) {
	     if ( CGAL::lexicographically_xyz_smaller(hfc->source()->source()->point(),
						      setemp->source()->source()->point()))
	       setemp = hfc;
	   }
	   *fc = make_object(setemp);
	 }
       }
     }
   }
 }

  void transform( const Aff_transformation_3& aff) {
    
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
	  if(Infi_box::is_infibox_corner((*li)->point())) {
	    li2 = corner_list.begin();
	    while(li2 != corner_list.end() && (*li2)->point() != (*li)->point()) ++li2;
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

      Association A;
      SNC_external_structure es(snc());
      es.clear_external_structure();
      for(li = vertex_list.begin(); li != vertex_list.end();++li){
	if(Infi_box::is_complex_facet_infibox_intersection(**li)) {
	  Vertex_handle v2;
	  Vertex_handle v1 = cstr.create_for_infibox_overlay(*li);
	  v1->point() = normalized(Infi_box::normalize_transformed_vertex((*li)->point().transform(aff)));
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
	    while(li2 != corner_list.end() && (*li2)->point() != v2->point()) ++li2;
	    if(li2 != corner_list.end())
	      delete_list.push_back(*li2);
	    break;
	  default: CGAL_error_msg( "wrong value");
	  }
	  Vertex_handle v = snc().new_vertex(v1->point(), (*li)->mark());
	  SM_overlayer O(&*v);
	  O.subdivide(&*v1, &*v2, A);
	  AND _and;
	  O.select(_and);
	  O.simplify(A);
	  snc().delete_vertex(v1);
	  snc().delete_vertex(v2);
	}
	
	if(Infi_box::is_infibox_corner((*li)->point())) {
	  SM_decorator SD(&**li);
	  if(SD.number_of_svertices() < 4)
	    continue;
	  li2 = corner_list.begin();
	  while(li2 != corner_list.end() && (*li2)->point() != (*li)->point()) ++li2;
	  CGAL_assertion(li2 != corner_list.end());
	  if(*li == *li2) {
	    delete_list.push_back(*li2);
	    *li2 = cstr.create_from_point_on_infibox_vertex((*li)->point());
	  }
	} else 
	  snc().delete_vertex(*li);	  
      }

      for(li = delete_list.begin(); li != delete_list.end(); ++li)
	snc().delete_vertex(*li);

      if(!aff.is_even())
	change_orientation();

      while(cstr.erase_redundant_vertices()) ;
      cstr.correct_infibox_sedge_marks();

      build_external_structure();
      cstr.correct_infibox_sface_marks();

      // are the upcoming lines necessary? 
      SNC_point_locator* old_pl = pl();
      pl() = pl()->clone();
      pl()->initialize(&snc());
      delete old_pl;   

    } else {
      Halffacet_iterator fi;
      CGAL_forall_halffacets(fi,snc()) {
	if(is_standard(fi) || ninety) {
	  fi->plane() = fi->plane().transform( aff);
#ifdef CGAL_NEF3_FACET_WITH_BOX 
	  typedef typename Halffacet::Box Box;
	  bool first = true;
	  Halffacet_cycle_iterator cycle_it = fi->facet_cycles_begin();
	  if( cycle_it.is_shalfedge() ) {
	    SHalfedge_iterator edge_it(cycle_it);
	    SHalfedge_around_facet_circulator
	      start( edge_it ), end( edge_it );
	    CGAL_For_all( start, end ) {
	      const Point_3& p = start->source()->source()->point();
	      typename Kernel::FT q[3];
	      q[0] = p.x();
	      q[1] = p.y();
	      q[2] = p.z();
	      if(first) {
		fi->b = Box(q,q);
		first = false;
	      } else
		fi->b.extend(q);
	    }
	  } else
	    CGAL_error_msg( "is facet first cycle a SHalfloop?"); 
#endif
	}
      }    

      if(!aff.is_even())
	change_orientation(true);

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

    SNC_constructor C(snc());
    C.assign_indices(); 
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

  bool contains(Object_handle /*h*/) const
  /*{\Mop  returns true iff the object |h| is contained in the set
  represented by |\Mvar|.}*/
    // { SNC_point_locator PL(snc()); return PL.mark(h);} 
    { CGAL_error_msg( "not implemented."); return false;}

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
    if(assign(v,o)) return make_object(Vertex_const_handle(v));
    if(assign(e,o)) return make_object(Halfedge_const_handle(e));
    if(assign(f,o)) return make_object(Halffacet_const_handle(f));
    if(assign(c,o)) return make_object(Volume_const_handle(c));

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
Nef_polyhedron_3( Content space) {
  CGAL_NEF_TRACEN("construction from empty or space.");
  empty_rep();
  set_snc(snc());
  if(Infi_box::extended_kernel()) {
    initialize_infibox_vertices(space);
    build_external_structure();
  } else {
    build_external_structure();
    snc().volumes_begin()->mark() = (space == COMPLETE) ? 1 : 0;
  }
}

template <typename Kernel, typename Items, typename Mark>
Nef_polyhedron_3<Kernel,Items, Mark>::
Nef_polyhedron_3(const Plane_3& h, Boundary b) {
  CGAL_NEF_TRACEN("construction from plane "<<h);
  empty_rep();
  set_snc(snc());
  SNC_constructor C(snc());
  Infi_box::create_vertices_of_box_with_plane(C,h,(b==INCLUDED));
  build_external_structure();
  /*
  if(Infi_box::extended_kernel()) {
    SNC_structure snc1, snc2;
    SNC_point_locator* pl1 = new SNC_point_locator_default; 
    SNC_point_locator* pl2 = new SNC_point_locator_default; 

    SNC_constructor c1(snc1); 
    Infi_box::initialize_infibox_vertices(c1, true);
    SNC_external_structure es1(snc1, pl1);
    es1.build_external_structure();

    SNC_constructor c2(snc2);
    c2.create_vertices_for_halfspace(h, b);
    SNC_external_structure es2(snc2, pl2);
    es2.pair_up_halfedges();
    es2.link_shalfedges_to_facet_cycles();
    c2.create_facets_and_volumes_of_halfspace(h);
    pl2->initialize(&snc2);

    AND _and;
    Binary_operation bo(snc());
    bo(pl(), snc1, pl1, snc2, pl2, _and);
    
    delete pl1;
    delete pl2;
  } else
    CGAL_error_msg
      ("Constructor is only available with extended kernels");
  */
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
    v->mark() = !v->mark(); 
    SM_decorator SM(&*v);
    SM.extract_complement();
  }

  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) f->mark() = !f->mark(); 
 
  Volume_iterator c;
  CGAL_forall_volumes(c,D) 
    //    if(!(Infi_box::extended_kernel && c==D.volumes_begin()))
      c->mark() = !c->mark();
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
    v->mark() = false;
    SM_decorator SM(&*v);
    SM.extract_interior();
  }
  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) f->mark() = false;

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
    v->mark() = true;
    SM_decorator SM(&*v);
    SM.extract_boundary();
  }
  Halffacet_iterator f;
  CGAL_forall_halffacets(f,D) f->mark() = true;
  Volume_iterator c;
  CGAL_forall_volumes(c,D) c->mark() = false;
  simplify();
}

} //namespace CGAL

#endif //CGAL_NEF_POLYHEDRON_3_H
