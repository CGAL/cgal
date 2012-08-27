// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYHEDRON_2_H
#define CGAL_NEF_POLYHEDRON_2_H

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4800) // complaint about performance in std::map where we can't do anything
#endif                          


#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Nef_2/HDS_items.h>
#include <CGAL/HalfedgeDS_default.h>

#include <CGAL/Is_extended_kernel.h>
#include <CGAL/Nef_2/PM_explorer.h>
#include <CGAL/Nef_2/PM_decorator.h>
#include <CGAL/Nef_2/PM_io_parser.h>
#include <CGAL/Nef_2/PM_overlayer.h>
//#include <CGAL/Nef_2/PM_transformer.h>
#include <CGAL/Nef_2/PM_point_locator.h>
#include <CGAL/Nef_2/Bounding_box_2.h>
#include <vector>
#include <list>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 11
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename T, typename I, typename M> class Nef_polyhedron_2;
template <typename T, typename I, typename M> class Nef_polyhedron_2_rep;

template <typename T, typename I, typename M>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_2<T,I,M>&); 

template <typename T, typename I, typename M>
std::istream& operator>>(std::istream&, Nef_polyhedron_2<T,I,M>&);

template <typename T, typename Items, typename Mark_>
class Nef_polyhedron_2_rep 
{ 
  typedef Nef_polyhedron_2_rep<T,Items,Mark_> Self;
  friend class Nef_polyhedron_2<T,Items,Mark_>;

  struct HDS_traits {
    typedef typename T::Point_2 Point;
    typedef Mark_                Mark;
  };

public: // gcc-3.3 otherwise claims that Decorator in Polyhedron_2 is private
  typedef CGAL_HALFEDGEDS_DEFAULT<HDS_traits,Items> Plane_map;
  typedef CGAL::PM_const_decorator<Plane_map>           Const_decorator;
  typedef CGAL::PM_decorator<Plane_map>                 Decorator;
  typedef CGAL::PM_naive_point_locator<Decorator,T>     Slocator;
  typedef CGAL::PM_point_locator<Decorator,T>           Locator;
  typedef CGAL::PM_overlayer<Decorator,T>               Overlayer;

private:

  Plane_map pm_; 
  Locator* pl_;
  
  void init_locator() 
  { 
    if ( !pl_ ) 
      pl_ = new Locator(pm_); 
  }

  void clear_locator() 
  { 
    if ( pl_ ) { 
      delete pl_; 
      pl_=0; 
    } 
  }

public:
  Nef_polyhedron_2_rep() 
    : pm_(), pl_(0) 
  {}

  Nef_polyhedron_2_rep(const Self& ) 
    : pm_(), pl_(0) 
  {}

  ~Nef_polyhedron_2_rep() 
  { 
    pm_.clear(); 
    clear_locator(); 
  }
};

/*{\Moptions print_title=yes }*/ 
/*{\Manpage {Nef_polyhedron_2}{T}{Nef Polyhedra in the Plane}{N}}*/

/*{\Mdefinition
An instance of data type |\Mname| is a subset of the plane that is
the result of forming complements and intersections starting from a
finite set |H| of half-spaces. |\Mtype| is closed under all binary set
operations |intersection|, |union|, |difference|, |complement| and
under the topological operations |boundary|, |closure|, and
|interior|.

The template parameter |T| is specified via an extended kernel
concept. |T| must be a model of the concept |ExtendedKernelTraits_2|.
}*/

template <typename T, typename Items_=HDS_items, typename Mark_=bool>
class Nef_polyhedron_2 
  : public Handle_for< Nef_polyhedron_2_rep<T,Items_,Mark_> >
{ 
public:
typedef T Extended_kernel;
static  T EK; // static extended kernel

  /*{\Mtypes 7}*/
  typedef Mark_ Mark;
  /*{\Xtypemember marking set membership or exclusion.}*/
  typedef Items_ Items;
  typedef Nef_polyhedron_2<T,Items,Mark> Self;
  typedef Handle_for< Nef_polyhedron_2_rep<T,Items,Mark> > Base;
  typedef typename T::Point_2   Extended_point;
  typedef typename T::Segment_2 Extended_segment;

  typedef typename T::Standard_line_2 Line;
  /*{\Mtypemember the oriented lines modeling half-planes}*/
  typedef typename T::Standard_point_2 Point;
  /*{\Mtypemember the affine points of the plane.}*/
  typedef typename T::Standard_direction_2 Direction;
  /*{\Mtypemember directions in our plane.}*/
  typedef typename T::Standard_aff_transformation_2  Aff_transformation;
  /*{\Mtypemember affine transformations of the plane.}*/


  // types for choosing the right constructor
  struct Polylines {};
  struct Polygons {};

  enum Operation { JOIN=0 };

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  enum Content { EMPTY=0, COMPLETE=1 };
  /*{\Menum construction selection}*/

  static const Polylines POLYLINES;
  static const Polygons POLYGONS;
protected:
  struct AND { bool operator()(bool b1, bool b2)  const { return b1&&b2; }  };
  struct OR { bool operator()(bool b1, bool b2)   const { return b1||b2; }  };
  struct DIFF { bool operator()(bool b1, bool b2) const { return b1&&!b2; } };
  struct XOR { bool operator()(bool b1, bool b2)  const 
                                           { return (b1&&!b2)||(!b1&&b2); } };

  typedef Nef_polyhedron_2_rep<T,Items,Mark>           Nef_rep;
  typedef typename Nef_rep::Plane_map       Plane_map;
  typedef typename Nef_rep::Decorator       Decorator;
  typedef typename Nef_rep::Const_decorator Const_decorator;
  typedef typename Nef_rep::Overlayer       Overlayer;
  //typedef typename Nef_rep::T               Transformer;
  typedef typename Nef_rep::Slocator        Slocator;
  typedef typename Nef_rep::Locator         Locator;

  using Base::ptr;
  using Base::is_shared;

  Plane_map& pm() { return ptr()->pm_; } 
  const Plane_map& pm() const { return ptr()->pm_; } 

  friend std::ostream& operator<< <>
      (std::ostream& os, const Nef_polyhedron_2<T,Items,Mark>& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Nef_polyhedron_2<T,Items,Mark>& NP);

  typedef typename Decorator::Vertex_handle         Vertex_handle;
  typedef typename Decorator::Halfedge_handle       Halfedge_handle;
  typedef typename Decorator::Face_handle           Face_handle;
  typedef typename Decorator::Vertex_const_handle   Vertex_const_handle;
  typedef typename Decorator::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Decorator::Face_const_handle     Face_const_handle;

  typedef typename Decorator::Vertex_iterator       Vertex_iterator;
  typedef typename Decorator::Halfedge_iterator     Halfedge_iterator;
  typedef typename Decorator::Face_iterator         Face_iterator;
  typedef typename Const_decorator::Vertex_const_iterator   
                                                    Vertex_const_iterator;
  typedef typename Const_decorator::Halfedge_const_iterator 
                                                    Halfedge_const_iterator;
  typedef typename Const_decorator::Face_const_iterator     
                                                    Face_const_iterator;

  typedef Bounding_box_2<typename Is_extended_kernel<Extended_kernel>::value_type, 
                         Extended_kernel> Box_2;

  struct Except_frame_box_edges {
    Decorator D_; 
    Face_handle f_;

    Except_frame_box_edges(Plane_map& P) 
      : D_(P), f_(D_.faces_begin()) 
    {}
    
    bool 
    operator()(Halfedge_handle e, const Tag_true& ) const
    { 
      return D_.face(e)==f_ || D_.face(D_.twin(e))==f_; 
    }

    bool 
    operator()(Halfedge_handle /*e*/, const Tag_false& ) const
    { 
      return false;
    }
    bool
    operator()(Halfedge_handle e) const
    {
      return this->operator()(e, typename Is_extended_kernel<Extended_kernel>::value_type());
    }

  };

  friend struct Except_frame_box_edges;

  typedef std::list<Extended_segment>      ES_list;
  typedef typename ES_list::const_iterator ES_iterator;

  void fill_with_frame_segs(ES_list& L, const Tag_true& ) const
  /*{\Xop fills the list with the four segments which span our frame,
     the convex hull of SW,SE,NW,NE.}*/
  { L.push_back(Extended_segment(EK.SW(),EK.NW()));
    L.push_back(Extended_segment(EK.SW(),EK.SE()));
    L.push_back(Extended_segment(EK.NW(),EK.NE()));
    L.push_back(Extended_segment(EK.SE(),EK.NE()));
  }

  void fill_with_frame_segs(ES_list& , const Tag_false& ) const
  {}

  void fill_with_frame_segs(ES_list& L) const
  { 

    fill_with_frame_segs(L, typename Is_extended_kernel<Extended_kernel>::value_type());
  }

  struct Link_to_iterator {
    const Decorator& D;
    Halfedge_handle _e;
    Vertex_handle   _v;
    ES_iterator     _it;
    Mark            _m;
    Link_to_iterator(const Decorator& d, ES_iterator it, Mark m) : 
      D(d), _e(), _v(), _it(it), _m(m) {}

    void supporting_segment(Halfedge_handle e, ES_iterator it) 
    { if ( it == _it ) _e = e; D.mark(e) = _m; }
    void trivial_segment(Vertex_handle v, ES_iterator it) 
    { if ( it == _it ) _v = v; D.mark(v) = _m; }
    void starting_segment(Vertex_handle v, ES_iterator) 
    { D.mark(v) = _m; }
    void passing_segment(Vertex_handle v, ES_iterator) 
    { D.mark(v) = _m; }
    void ending_segment(Vertex_handle v, ES_iterator) 
    { D.mark(v) = _m; }

  };

  template<typename IT>
  struct From_intersecting_polygons {

    Unique_hash_map<Halfedge_handle,IT>& halfedge2iterator;

    From_intersecting_polygons(Unique_hash_map<Halfedge_handle,IT>& e2i) 
      : halfedge2iterator(e2i) {}

    void supporting_segment(Halfedge_handle e, IT it) 
    { 
      halfedge2iterator[e->opposite()] = 
	halfedge2iterator[e] = it; e->mark() = true;}      

    void trivial_segment(Vertex_handle v, IT) 
    { v->mark() = true; }
    void starting_segment(Vertex_handle v, IT) 
    { v->mark() = true; }
    void passing_segment(Vertex_handle v, IT) 
    { v->mark() = true; }
    void ending_segment(Vertex_handle v, IT) 
    { v->mark() = true; }
  };  

  friend struct Link_to_iterator;

  void clear_outer_face_cycle_marks(const Tag_true&) 
  { // unset all frame marks
    Decorator D(pm());
    Face_iterator f = D.faces_begin(); 
    D.mark(f) = false;
    Halfedge_handle e = D.holes_begin(f);
    D.set_marks_in_face_cycle(e, false);
  }

  void clear_outer_face_cycle_marks(const Tag_false&)
  {}

  void clear_outer_face_cycle_marks()
  {
    clear_outer_face_cycle_marks(typename Is_extended_kernel<Extended_kernel>::value_type());
  }

public:
  /*{\Mcreation 3}*/
  Nef_polyhedron_2(Content plane = EMPTY) : Base(Nef_rep())
  /*{\Mcreate creates an instance |\Mvar| of type |\Mname|
  and initializes it to the empty set if |plane == EMPTY|
  and to the whole plane if |plane == COMPLETE|.}*/
  {
    ES_list L;
    fill_with_frame_segs(L);
    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), false);
    D.create(L.begin(),L.end(),I);
    D.mark(--D.faces_end()) = bool(plane);
  }


  Nef_polyhedron_2(const Line& l, Boundary line = INCLUDED) : Base(Nef_rep())
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the half-plane
  left of |l| including |l| if |line==INCLUDED|, excluding |l| if 
  |line==EXCLUDED|.}*/  
  {   CGAL_NEF_TRACEN("Nconstruction from line "<<l);
    ES_list L;
    fill_with_frame_segs(L);
    if(check_tag(typename Is_extended_kernel<Extended_kernel>::value_type())) {
      Extended_point ep1 = EK.construct_opposite_point(l);
      Extended_point ep2 = EK.construct_point(l);
      L.push_back(EK.construct_segment(ep1,ep2));
    }
    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), false);
    D.create(L.begin(),L.end(),I);
    if(check_tag(typename Is_extended_kernel<Extended_kernel>::value_type())) {
      CGAL_assertion( I._e != Halfedge_handle() );
      Halfedge_handle el = I._e;
      if ( D.point(D.target(el)) != EK.target(L.back()) )
	el = D.twin(el);
      D.mark(D.face(el)) = true;
      D.mark(el) = bool(line);
    } else {
      D.mark(--D.faces_end()) = bool(EMPTY);
      std::cerr << "Constructor not available with standard kernel. "
                   " Returned empty polygon!" << std::endl;
    }
  }


  template <class Forward_iterator>
  Nef_polyhedron_2(Forward_iterator it, Forward_iterator end,
    Boundary b = INCLUDED) : Base(Nef_rep())
  /*{\Mcreate creates a Nef polyhedron |\Mvar| from the simple polygon
  |P| spanned by the list of points in the iterator range |[it,end)| and
  including its boundary if |b = INCLUDED| and excluding the boundary
  otherwise. |Forward_iterator| has to be an iterator with value type
  |Point|. This construction expects that |P| is simple. The degenerate
  cases where |P| contains no point, one point or spans just one segment
  (two points) are correctly handled. In all degenerate cases there's
  only one unbounded face adjacent to the degenerate polygon. If |b ==
  INCLUDED| then |\Mvar| is just the boundary. If |b == EXCLUDED| then
  |\Mvar| is the whole plane without the boundary.}*/
  { 
    ES_list L;
    fill_with_frame_segs(L);
    bool empty = false;  
    if (it != end) 
      {
        Extended_point ef, ep = ef = EK.construct_point(*it);
        Forward_iterator itl=it; ++itl;
        if (itl == end) // case only one point
          L.push_back(EK.construct_segment(ep,ep));
        else { // at least one segment
          while( itl != end ) {
            Extended_point en = EK.construct_point(*itl);
            L.push_back(EK.construct_segment(ep,en));
            ep = en; ++itl;
          }
          L.push_back(EK.construct_segment(ep,ef));
        }
      }

    else empty = true;
    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), true);
    D.create(L.begin(),L.end(),I);
    if ( empty ) {
      D.mark(--D.faces_end()) = !bool(b); return; }
    CGAL_assertion( I._e != Halfedge_handle() || I._v != Vertex_handle() );

    if ( EK.is_degenerate(L.back()) ) {
      // its a point
      CGAL_assertion(I._v != Vertex_handle());
      D.mark(D.face(I._v)) = !bool(b); D.mark(I._v) = b;
    } else {
      // at least one segment
      Halfedge_handle el = I._e;
      if ( D.point(D.target(el)) != EK.target(L.back()) )
	el = D.twin(el);  
      D.set_marks_in_face_cycle(el,bool(b));
      unsigned int n = 
        check_tag(typename Is_extended_kernel<Extended_kernel>::value_type()) ? 2 : 1;
      if ( D.number_of_faces() > n ) D.mark(D.face(el)) = true;
      else                           D.mark(D.face(el)) = !bool(b);
    }

    clear_outer_face_cycle_marks();
  }

  // The constructor which takes an iterator range of polygons
  template <class Forward_iterator>
  Nef_polyhedron_2(Forward_iterator pit, Forward_iterator pend,	      
		   Polygons, Operation op = JOIN) : Base(Nef_rep()) { 

    CGAL_assertion(op==JOIN);

    typedef typename std::iterator_traits<Forward_iterator>::value_type
      iterator_pair;
    typedef typename iterator_pair::first_type point_iterator;
    point_iterator it, itl, end;

    ES_list L;

    fill_with_frame_segs(L);
    for(;pit != pend; ++pit) {
      it = pit->first;
      end = pit->second;
      if (it != end) {
        Extended_point ef, ep = ef = EK.construct_point(*it);
        itl=it; ++itl;
        if (itl == end) // case only one point
          L.push_back(EK.construct_segment(ep,ep));
        else { // at least one segment
          while( itl != end ) {
            Extended_point en = EK.construct_point(*itl);
            L.push_back(EK.construct_segment(ep,en));
            ep = en; ++itl;
          }
          L.push_back(EK.construct_segment(ep,ef));
        }
      }
    }

    Overlayer D(pm());
    Unique_hash_map<Halfedge_handle,ES_iterator> e2i;
    From_intersecting_polygons<ES_iterator> fip(e2i);
    D.create(L.begin(),L.end(),fip);

    Face_handle outer_face;
    if(check_tag(typename Is_extended_kernel<Extended_kernel>::value_type()))
      outer_face = ++D.faces_begin();
    else
      outer_face = D.faces_begin();
    Halfedge_handle e;
    for(e=D.halfedges_begin(); e!=D.halfedges_end(); ++e) {
      if(&*e < &*(D.twin(e)) && EK.is_standard(D.source(e)->point())) {
	ES_iterator eit = e2i[e];
	if(lexicographically_xy_smaller(EK.standard_point(eit->source()),
					EK.standard_point(eit->target()))) {
	  if(lexicographically_xy_smaller(EK.standard_point(D.source(D.twin(e))->point()),
					  EK.standard_point(D.source(e)->point())))
	    e = D.twin(e);
	} else
	  if(lexicographically_xy_smaller(EK.standard_point(D.source(e)->point()),
					  EK.standard_point(D.source(D.twin(e))->point())))
	    e = D.twin(e);
	if(D.face(e) != outer_face)
	  D.mark(D.face(e)) = true;
      }
    }
    
    D.simplify(Except_frame_box_edges(pm()));
    clear_outer_face_cycle_marks();
  }


  // The constructor which takes an iterator range of polylines
  template <class Forward_iterator>
  Nef_polyhedron_2(Forward_iterator pit, Forward_iterator pend,	      
		   Polylines) : Base(Nef_rep()) { 

    typedef typename std::iterator_traits<Forward_iterator>::value_type 
      iterator_pair;
    typedef typename iterator_pair::first_type point_iterator;
    point_iterator it, itl, end;

    ES_list L;

    fill_with_frame_segs(L);
    for(;pit != pend; ++pit) {
      it = pit->first;
      end = pit->second;
      if (it != end) {
        Extended_point ep  = EK.construct_point(*it);
        itl=it; ++itl;
        if (itl == end) // case only one point
          L.push_back(EK.construct_segment(ep,ep));
        else { // at least one segment
          while( itl != end ) {
            Extended_point en = EK.construct_point(*itl);
            L.push_back(EK.construct_segment(ep,en));
            ep = en;
	    ++itl;
          }
        }
      }
    }

    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), true);
    D.create(L.begin(),L.end(),I, Overlayer::POLYLINE);
    
    clear_outer_face_cycle_marks();
  }

  Nef_polyhedron_2(const Nef_polyhedron_2<T,Items,Mark>& N1) : Base(N1) {}
  Nef_polyhedron_2& operator=(const Nef_polyhedron_2<T,Items,Mark>& N1)
  { Base::operator=(N1); return (*this); }
  ~Nef_polyhedron_2() {}



  template <class Forward_iterator>
  Nef_polyhedron_2(Forward_iterator first, Forward_iterator beyond, 
    double p) : Base(Nef_rep())
  /*{\Xcreate creates a random Nef polyhedron from the arrangement of
  the set of lines |S = set[first,beyond)|. The cells of the arrangement
  are selected uniformly at random with probability $p$. \precond $0 < p
  < 1$.}*/
  { CGAL_assertion(0<p && p<1);
    ES_list L; fill_with_frame_segs(L);
    while ( first != beyond ) {
      Extended_point ep1 = EK.construct_opposite_point(*first);
      Extended_point ep2 = EK.construct_point(*first);
      L.push_back(EK.construct_segment(ep1,ep2)); ++first;
    }
    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), false);
    D.create(L.begin(),L.end(),I);

    boost::rand48 rng;
    boost::uniform_real<> dist(0,1);
    boost::variate_generator<boost::rand48&, boost::uniform_real<> > get_double(rng,dist);

    Vertex_iterator v; Halfedge_iterator e; Face_iterator f;
    for (v = D.vertices_begin(); v != D.vertices_end(); ++v)
      D.mark(v) = ( get_double() < p ? true : false );
    for (e = D.halfedges_begin(); e != D.halfedges_end(); ++(++e))
      D.mark(e) = ( get_double() < p ? true : false );
    for (f = D.faces_begin(); f != D.faces_end(); ++f)
      D.mark(f) = ( get_double() < p ? true : false );
    D.simplify(Except_frame_box_edges(pm()));
    clear_outer_face_cycle_marks(); 
  }



  protected:
  Nef_polyhedron_2(const Plane_map& H, bool clone=true) : Base(Nef_rep()) 
  /*{\Xcreate makes |\Mvar| a new object.  If |clone==true| then the
  underlying structure of |H| is copied into |\Mvar|.}*/
  { if (clone) {
      Decorator D(pm()); // a decorator working on the rep plane map
      D.clone(H);        // cloning H into pm()
    }
  }
  void clone_rep() { *this = Nef_polyhedron_2<T,Items,Mark>(pm()); }

  /*{\Moperations 4 3 }*/
  public:

  void clear(Content plane = EMPTY)
  { *this = Nef_polyhedron_2(plane); }
  /*{\Mop makes |\Mvar| the empty set if |plane == EMPTY| and the
  full plane if |plane == COMPLETE|.}*/

  bool is_empty() const
  /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
  { Const_decorator D(pm());
    Face_const_iterator f = D.faces_begin();
    if(check_tag(typename Is_extended_kernel<Extended_kernel>::
		 value_type()))
      return (D.number_of_vertices()==4 &&
	      D.number_of_edges()==4 &&
	      D.number_of_faces()==2 &&
	      D.mark(++f) == false);
    else
      return (D.number_of_vertices()==0 &&
	      D.number_of_edges()==0 &&
	      D.number_of_faces()==1 &&
		D.mark(f) == false);    
  }

  bool is_plane() const
  /*{\Mop returns true if |\Mvar| is the whole plane, false otherwise.}*/
  { Const_decorator D(pm());
    Face_const_iterator f = D.faces_begin();
    if(check_tag(typename Is_extended_kernel<Extended_kernel>::
         value_type()))
      return (D.number_of_vertices()==4 &&
              D.number_of_edges()==4 &&
              D.number_of_faces()==2 &&
              D.mark(++f) == true);
    else
      return (D.number_of_vertices()==0 &&
          D.number_of_edges()==0 &&
          D.number_of_faces()==1 &&
          D.mark(f) == true);
  }

  void extract_complement()
  { CGAL_NEF_TRACEN("extract complement");
  if ( this->is_shared() ) {
	  clone_rep();
  }
    Overlayer D(pm());
    Vertex_iterator v, vend = D.vertices_end();
    for(v = D.vertices_begin(); v != vend; ++v)      D.mark(v) = !D.mark(v);
    Halfedge_iterator e, eend = D.halfedges_end();
    for(e = D.halfedges_begin(); e != eend; ++(++e)) D.mark(e) = !D.mark(e);
    Face_iterator f, fend = D.faces_end();
    for(f = D.faces_begin(); f != fend; ++f)         D.mark(f) = !D.mark(f);
    clear_outer_face_cycle_marks();
  }

  void extract_interior()
  { CGAL_NEF_TRACEN("extract interior");
    if ( this->is_shared() ) clone_rep();
    Overlayer D(pm());
    Vertex_iterator v, vend = D.vertices_end();
    for(v = D.vertices_begin(); v != vend; ++v)      D.mark(v) = false;
    Halfedge_iterator e, eend = D.halfedges_end();
    for(e = D.halfedges_begin(); e != eend; ++(++e)) D.mark(e) = false;
    D.simplify(Except_frame_box_edges(pm()));
  }


  void extract_boundary()
  { CGAL_NEF_TRACEN("extract boundary");
    if ( this->is_shared() ) clone_rep();
    Overlayer D(pm());
    Vertex_iterator v, vend = D.vertices_end();
    for(v = D.vertices_begin(); v != vend; ++v)      D.mark(v) = true;
    Halfedge_iterator e, eend = D.halfedges_end();
    for(e = D.halfedges_begin(); e != eend; ++(++e)) D.mark(e) = true;
    Face_iterator f, fend = D.faces_end();
    for(f = D.faces_begin(); f != fend; ++f)         D.mark(f) = false;
    clear_outer_face_cycle_marks();
    D.simplify(Except_frame_box_edges(pm()));
  }

  void extract_closure()
  /*{\Xop converts |\Mvar| to its closure. }*/
  { CGAL_NEF_TRACEN("extract closure");
    extract_complement();
    extract_interior();
    extract_complement();
  }

  void extract_regularization()
  /*{\Xop converts |\Mvar| to its regularization. }*/
  { CGAL_NEF_TRACEN("extract regularization");
    extract_interior();
    extract_closure();
  }

  /*{\Mtext \headerline{Constructive Operations}}*/

  Nef_polyhedron_2<T,Items,Mark> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane.}*/
  { Nef_polyhedron_2<T,Items,Mark> res = *this;
    res.extract_complement();
    return res;
  }


  Nef_polyhedron_2<T,Items,Mark> interior() const
  /*{\Mop returns the interior of |\Mvar|.}*/
  { Nef_polyhedron_2<T,Items,Mark> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_2<T,Items,Mark> closure() const
  /*{\Mop returns the closure of |\Mvar|.}*/
  { Nef_polyhedron_2<T,Items,Mark> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_2<T,Items,Mark> boundary() const
  /*{\Mop returns the boundary of |\Mvar|.}*/
  { Nef_polyhedron_2<T,Items,Mark> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_2<T,Items,Mark> regularization() const
  /*{\Mop returns the regularized polyhedron (closure of interior).}*/
  { Nef_polyhedron_2<T,Items,Mark> res = *this;
    res.extract_regularization();
    return res;
  }


  Nef_polyhedron_2<T,Items,Mark> intersection(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/
  { Nef_polyhedron_2<T,Items,Mark> res(pm(),false); // empty, no frame
    Overlayer D(res.pm());
    D.subdivide(pm(),N1.pm());
    AND _and; D.select(_and);
    res.clear_outer_face_cycle_marks();
    D.simplify(Except_frame_box_edges(res.pm()));
    return res;
  }


  Nef_polyhedron_2<T,Items,Mark> join(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/
  { Nef_polyhedron_2<T,Items,Mark> res(pm(),false); // empty, no frame
    Overlayer D(res.pm());
    D.subdivide(pm(),N1.pm());
    OR _or; D.select(_or);
    res.clear_outer_face_cycle_marks();
    D.simplify(Except_frame_box_edges(res.pm()));
    return res;
  }

  Nef_polyhedron_2<T,Items,Mark> difference(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/
  { Nef_polyhedron_2<T,Items,Mark> res(pm(),false); // empty, no frame
    Overlayer D(res.pm());
    D.subdivide(pm(),N1.pm());
    DIFF _diff; D.select(_diff);
    res.clear_outer_face_cycle_marks();
    D.simplify(Except_frame_box_edges(res.pm()));
    return res;
  }    

  Nef_polyhedron_2<T,Items,Mark> symmetric_difference(
    const Nef_polyhedron_2<T,Items,Mark>& N1) const
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/
  { Nef_polyhedron_2<T,Items,Mark> res(pm(),false); // empty, no frame
    Overlayer D(res.pm());
    D.subdivide(pm(),N1.pm());
    XOR _xor; D.select(_xor);
    res.clear_outer_face_cycle_marks();
    D.simplify(Except_frame_box_edges(res.pm()));
    return res;
  }

  #if 0
  Nef_polyhedron_2<T,Items,Mark> transform(const Aff_transformation& t) const
  /*{\Mop returns $t(|\Mvar|)$.}*/
  { Nef_polyhedron_2<T,Items,Mark> res(pm()); // cloned
    Transformer PMT(res.pm());
    PMT.transform(t);
    return res;
  }
  #endif


  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement} respectively. There are also the
  corresponding modification operations |*=,+=,-=,^=|.}*/

  Nef_polyhedron_2<T,Items,Mark>  operator*(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return intersection(N1); }

  Nef_polyhedron_2<T,Items,Mark>  operator+(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return join(N1); }

  Nef_polyhedron_2<T,Items,Mark>  operator-(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return difference(N1); }

  Nef_polyhedron_2<T,Items,Mark>  operator^(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return symmetric_difference(N1); }

  Nef_polyhedron_2<T,Items,Mark>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_2<T,Items,Mark>& operator*=(const Nef_polyhedron_2<T,Items,Mark>& N1)
  { *this = intersection(N1); return *this; }

  Nef_polyhedron_2<T,Items,Mark>& operator+=(const Nef_polyhedron_2<T,Items,Mark>& N1)
  { *this = join(N1); return *this; }

  Nef_polyhedron_2<T,Items,Mark>& operator-=(const Nef_polyhedron_2<T,Items,Mark>& N1)
  { *this = difference(N1); return *this; }

  Nef_polyhedron_2<T,Items,Mark>& operator^=(const Nef_polyhedron_2<T,Items,Mark>& N1)
  { *this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality, respectively.}*/

  bool operator==(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return symmetric_difference(N1).is_empty(); }

  bool operator!=(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return !operator==(N1); }  

  bool operator<=(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return difference(N1).is_empty(); } 

  bool operator<(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return difference(N1).is_empty() && !N1.difference(*this).is_empty(); } 

  bool operator>=(const Nef_polyhedron_2<T,Items,Mark>& N1) const
  { return N1.difference(*this).is_empty(); } 

  bool operator>(const Nef_polyhedron_2<T,Items,Mark>& N1) const   
  { return N1.difference(*this).is_empty() && !difference(N1).is_empty(); } 


  /*{\Mtext \headerline{Exploration - Point location - Ray shooting}
  As Nef polyhedra are the result of forming complements 
  and intersections starting from a set |H| of half-spaces that are
  defined by oriented lines in the plane, they can be represented by
  an attributed plane map $M = (V,E,F)$. For topological queries
  within |M| the following types and operations allow exploration
  access to this structure.}*/

  /*{\Mtypes 3}*/
  typedef Const_decorator Topological_explorer;

  typedef CGAL::PM_explorer<Const_decorator,T> Explorer;
  /*{\Mtypemember a decorator to examine the underlying plane map. 
  See the manual page of |Explorer|.}*/

  typedef typename Locator::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying
  plane map. The kind of object |(vertex, halfedge, face)| can 
  be determined and the object can be assigned to a corresponding
  handle by the three functions:\\
  |bool assign(Vertex_const_handle& h, Object_handle)|\\
  |bool assign(Halfedge_const_handle& h, Object_handle)|\\
  |bool assign(Face_const_handle& h, Object_handle)|\\
  where each function returns |true| iff the assignment to
  |h| was done.}*/

  enum Location_mode { DEFAULT, NAIVE, LMWT };
  /*{\Menum selection flag for the point location mode.}*/


  /*{\Moperations 3 1 }*/

  void init_locator() const 
  { const_cast<Self*>(this)->ptr()->init_locator(); }
  const Locator& locator() const 
  { CGAL_assertion(ptr()->pl_); return *(ptr()->pl_); }


  bool contains(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the set
  represented by |\Mvar|.}*/
  { Slocator PL(pm()); return PL.mark(h); }

  bool contained_in_boundary(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the $1$-skeleton
  of |\Mvar|.}*/
  { Vertex_const_handle v;
    Halfedge_const_handle e;
    return  ( CGAL::assign(v,h) || CGAL::assign(e,h) );
  }

  Object_handle locate(const Point& p, Location_mode m = DEFAULT) const
  /*{\Mop  returns a generic handle |h| to an object (face, halfedge, vertex) 
  of the underlying plane map that contains the point |p| in its relative 
  interior. The point |p| is contained in the set represented by |\Mvar| if 
  |\Mvar.contains(h)| is true. The location mode flag |m| allows one to choose
  between different point location strategies.}*/
  { 
    if (m == DEFAULT || m == LMWT) {
      init_locator();
      Extended_point ep = EK.construct_point(p);
      return locator().locate(ep);
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      Extended_segment s(EK.construct_point(p),
			 PL.point(PL.vertices_begin()));
      return PL.locate(s); 
    }
    CGAL_error_msg("location mode not implemented.");
    return Object_handle();
  }

  struct INSET {
    const Const_decorator& D;
    INSET(const Const_decorator& Di) : D(Di) {}
    bool operator()(Vertex_const_handle v) const { return D.mark(v); }
    bool operator()(Halfedge_const_handle e) const { return D.mark(e); }
    bool operator()(Face_const_handle f) const { return D.mark(f); }
  };

  friend struct INSET;

  Object_handle ray_shoot(const Point& p, const Direction& d, 
                          Location_mode m = DEFAULT) const
  /*{\Mop returns a handle |h| with |\Mvar.contains(h)| that can be
  converted to a |Vertex_/Halfedge_/Face_const_handle| as described
  above. The object returned is intersected by the ray starting in |p|
  with direction |d| and has minimal distance to |p|.  The operation
  returns the null handle |NULL| if the ray shoot along |d| does not hit
  any object |h| of |\Mvar| with |\Mvar.contains(h)|. The location mode
  flag |m| allows one to choose between different point location
  strategies.}*/
  { 
    Extended_point ep, eq;
    if(!check_tag(typename Is_extended_kernel<Extended_kernel>::value_type())) {
      Const_decorator D(pm());
      Box_2 b(D.vertices_begin(), D.vertices_end());
      ep = EK.construct_point(p);
      eq = b.intersection_ray_bbox(p,d);
    } else {
      ep = EK.construct_point(p);
      eq = EK.construct_point(p,d);      
    }
      
    if (m == DEFAULT || m == LMWT) {
      init_locator();
      return locator().ray_shoot(EK.construct_segment(ep,eq),
                                 INSET(locator())); 
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      return PL.ray_shoot(EK.construct_segment(ep,eq),INSET(PL));
    }
    CGAL_error_msg("location mode not implemented.");
    return Object_handle();
  }

  struct INSKEL {
    bool operator()(Vertex_const_handle) const { return true; }
    bool operator()(Halfedge_const_handle) const { return true; }
    bool operator()(Face_const_handle) const { return false; }
  };

  Object_handle ray_shoot_to_boundary(const Point& p, const Direction& d, 
                Location_mode m = DEFAULT) const
  /*{\Mop returns a handle |h| that can be converted to a
  |Vertex_/Halfedge_const_handle| as described above. The object
  returned is part of the $1$-skeleton of |\Mvar|, intersected by the
  ray starting in |p| with direction |d| and has minimal distance to
  |p|.  The operation returns the null handle |NULL| if the ray shoot
  along |d| does not hit any $1$-skeleton object |h| of |\Mvar|. The
  location mode flag |m| allows one to choose between different point
  location strategies.}*/
  { 
    Extended_point ep, eq;
    if(!check_tag(typename Is_extended_kernel<Extended_kernel>::value_type())) {
      Const_decorator D(pm());
      Box_2 b(D.vertices_begin(), D.vertices_end());
      ep = EK.construct_point(p);
      eq = b.intersection_ray_bbox(p,d);
    } else {
      ep = EK.construct_point(p);
      eq = EK.construct_point(p,d);      
    }
      
    if (m == DEFAULT || m == LMWT) {
      init_locator();
      return locator().ray_shoot(EK.construct_segment(ep,eq),
                                 INSKEL()); 
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      return PL.ray_shoot(EK.construct_segment(ep,eq),INSKEL());
    }
    CGAL_error_msg("location mode not implemented.");
    return Object_handle();
  }


  Explorer explorer() const { return Explorer(pm(),EK); }
  /*{\Mop returns a decorator object which allows read-only access of
  the underlying plane map. See the manual page |Explorer| for its 
  usage.}*/


  /*{\Mimplementation Nef polyhedra are implemented on top of a halfedge
  data structure and use linear space in the number of vertices, edges
  and facets.  Operations like |empty| take constant time. The
  operations |clear|, |complement|, |interior|, |closure|, |boundary|,
  |regularization|, input and output take linear time. All binary set
  operations and comparison operations take time $O(n \log n)$ where $n$
  is the size of the output plus the size of the input.

  The point location and ray shooting operations are implemented in
  two flavors. The |NAIVE| operations run in linear query time without
  any preprocessing, the |DEFAULT| operations (equals |LMWT|) run in
  sub-linear query time, but preprocessing is triggered with the first
  operation. Preprocessing takes time $O(N^2)$, the sub-linear point
  location time is either logarithmic when LEDA's persistent
  dictionaries are present or if not then the point location time is
  worst-case linear, but experiments show often sublinear runtimes.  Ray
  shooting equals point location plus a walk in the constrained
  triangulation overlayed on the plane map representation. The cost of
  the walk is proportional to the number of triangles passed in
  direction |d| until an obstacle is met. In a minimum weight
  triangulation of the obstacles (the plane map representing the
  polyhedron) the theory provides a $O(\sqrt{n})$ bound for the number
  of steps. Our locally minimum weight triangulation approximates the
  minimum weight triangulation only heuristically (the calculation of
  the minimum weight triangulation is conjectured to be NP hard). Thus
  we have no runtime guarantee but a strong experimental motivation for
  its approximation.}*/

  /*{\Mexample Nef polyhedra are parameterized by a so-called extended
  geometric kernel. There are three kernels, one based on a homogeneous
  representation of extended points called |Extended_homogeneous<RT>|
  where |RT| is a ring type providing additionally a |gcd| operation and
  one based on a cartesian representation of extended points called
  |Extended_cartesian<NT>| where |NT| is a field type, and finally
  |Filtered_extended_homogeneous<RT>| (an optimized version of the
  first).

  The member types of |Nef_polyhedron_2< Extended_homogeneous<NT> >|
  map to corresponding types of the CGAL geometry kernel
  (e.g. |Nef_polyhedron::Line| equals
  |CGAL::Homogeneous<leda_integer>::Line_2| in the example below).
  \begin{Mverb}
  #include <CGAL/basic.h>
  #include <CGAL/leda_integer.h>
  #include <CGAL/Extended_homogeneous.h>
  #include <CGAL/Nef_polyhedron_2.h>

  using namespace CGAL;
  typedef  Extended_homogeneous<leda_integer> Extended_kernel;
  typedef  Nef_polyhedron_2<Extended_kernel>  Nef_polyhedron;
  typedef  Nef_polyhedron::Line               Line;

  int main()
  {
    Nef_polyhedron N1(Line(1,0,0));
    Nef_polyhedron N2(Line(0,1,0), Nef_polyhedron::EXCLUDED);
    Nef_polyhedron N3 = N1 * N2; // line (*)
    return 0;
  }
  \end{Mverb}
  After line (*) |N3| is the intersection of |N1| and |N2|.}*/


}; // end of Nef_polyhedron_2

template <typename T, typename Items, typename Mark>
T Nef_polyhedron_2<T,Items,Mark>::EK;


template <typename T, typename Items, typename Mark>
const typename Nef_polyhedron_2<T,Items,Mark>::Polygons Nef_polyhedron_2<T,Items,Mark>::POLYGONS = typename Nef_polyhedron_2<T,Items,Mark>::Polygons();

template <typename T, typename Items, typename Mark>
const typename Nef_polyhedron_2<T,Items,Mark>::Polylines Nef_polyhedron_2<T,Items,Mark>::POLYLINES = typename Nef_polyhedron_2<T,Items,Mark>::Polylines();

template <typename T, typename Items, typename Mark>
std::ostream& operator<<
 (std::ostream& os, const Nef_polyhedron_2<T,Items,Mark>& NP)
{
  os << "Nef_polyhedron_2<" << NP.EK.output_identifier() << ">\n";
  typedef typename Nef_polyhedron_2<T,Items,Mark>::Decorator Decorator;
  CGAL::PM_io_parser<Decorator> O(os, NP.pm()); O.print();
  return os;
}

template <typename T, typename Items, typename Mark>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_2<T,Items,Mark>& NP)
{
  typedef typename Nef_polyhedron_2<T,Items,Mark>::Decorator Decorator;
  CGAL::PM_io_parser<Decorator> I(is, NP.pm()); 
  if (I.check_sep("Nef_polyhedron_2<") &&
      I.check_sep(NP.EK.output_identifier()) &&
      I.check_sep(">")) I.read();
  else {
    std::cerr << "Nef_polyhedron_2 input corrupted." << std::endl;
    NP = Nef_polyhedron_2<T,Items,Mark>();
  }
  typename Nef_polyhedron_2<T,Items,Mark>::Topological_explorer D(NP.explorer());
  D.check_integrity_and_topological_planarity();
  return is;
}



} //namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif //CGAL_NEF_POLYHEDRON_2_H
