// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_polyhedron_2.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Nef_polyhedron_2.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Nef polyhedra in the plane
// ============================================================================

#ifndef CGAL_NEF_POLYHEDRON_2_H
#define CGAL_NEF_POLYHEDRON_2_H

#if defined(_MSC_VER) || defined(__BORLANDC__)
#define CGAL_SIMPLE_HDS
#endif

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#ifndef CGAL_SIMPLE_HDS
#include <CGAL/Nef_2/HDS_items.h>
#include <CGAL/HalfedgeDS_default.h>
#else
#include <CGAL/Nef_2/HalfedgeDS_default_MSC.h>
#endif
#include <CGAL/Nef_2/PM_explorer.h>
#include <CGAL/Nef_2/PM_decorator.h>
#include <CGAL/Nef_2/PM_io_parser.h>
#include <CGAL/Nef_2/PM_overlayer.h>
//#include <CGAL/Nef_2/PM_transformer.h>
#include <CGAL/Nef_2/PM_point_locator.h>
#include <vector>
#include <list>

#undef _DEBUG
#define _DEBUG 11
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE

template <typename T> class Nef_polyhedron_2;
template <typename T> class Nef_polyhedron_2_rep;

template <typename T>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_2<T>&); 
template <typename T>
std::istream& operator>>(std::istream&, Nef_polyhedron_2<T>&);
template <typename T>
class Nef_polyhedron_2_rep : public Ref_counted
{ typedef Nef_polyhedron_2_rep<T> Self;
  friend class Nef_polyhedron_2<T>;
#ifndef CGAL_SIMPLE_HDS
  struct HDS_traits {
    typedef typename T::Point_2 Point;
    typedef bool                Mark;
  };
  typedef CGAL_HALFEDGEDS_DEFAULT<HDS_traits,HDS_items> Plane_map;
  typedef CGAL::PM_const_decorator<Plane_map>           Const_decorator;
  typedef CGAL::PM_decorator<Plane_map>                 Decorator;
  typedef CGAL::PM_naive_point_locator<Decorator,T>     Slocator;
  typedef CGAL::PM_point_locator<Decorator,T>           Locator;
  typedef CGAL::PM_overlayer<Decorator,T>               Overlayer;

#else
  struct HDS_traits {
    typedef typename T::Point_2 Point;
    typedef bool                Mark;
  };
  typedef CGAL::HalfedgeDS_default_MSC<HDS_traits>  Plane_map;
  typedef CGAL::PM_const_decorator<Plane_map>       Const_decorator;
  typedef CGAL::PM_decorator<Plane_map>             Decorator;
  typedef CGAL::PM_naive_point_locator<Decorator,T> Slocator;
  typedef CGAL::PM_point_locator<Decorator,T>       Locator;
  typedef CGAL::PM_overlayer<Decorator,T>           Overlayer;

#endif
  //typedef CGAL::PM_transformer<Decorator,T> Transformer;
  Plane_map pm_; Locator* pl_;
  
  void init_locator() 
  { if ( !pl_ ) pl_ = new Locator(pm_); }
  void clear_locator() 
  { if ( pl_ ) delete pl_; pl_=0; }
public:
  Nef_polyhedron_2_rep() : Ref_counted(), pm_(), pl_(0) {}
  Nef_polyhedron_2_rep(const Self& R) : Ref_counted(), pm_(), pl_(0) {}
  ~Nef_polyhedron_2_rep() { pm_.clear(); clear_locator(); }
};

/*{\Moptions print_title=yes }*/ 
/*{\Manpage {Nef_polyhedron_2}{T}{Nef Polyhedra in the Plane}{N}}*/

/*{\Mdefinition
An instance of data type |\Mname| is a subset of the plane that is
the result of forming complements and intersections starting from a
finite set |H| of halfspaces. |\Mtype| is closed under all binary set
operations |intersection|, |union|, |difference|, |complement| and
under the topological operations |boundary|, |closure|, and
|interior|.

The template parameter |T| is specified via an extended kernel
concept. |T| must be a model of the concept |ExtendedKernelTraits_2|.
}*/

template <typename T>
class Nef_polyhedron_2 : public Handle_for< Nef_polyhedron_2_rep<T> >
{ 
public:
typedef T Extended_kernel;
static  T EK; // static extended kernel

  /*{\Mtypes 7}*/
  typedef Nef_polyhedron_2<T> Self;
  typedef Handle_for< Nef_polyhedron_2_rep<T> > Base;
  typedef typename T::Point_2   Extended_point;
  typedef typename T::Segment_2 Extended_segment;

  typedef typename T::Standard_line_2 Line;
  /*{\Mtypemember the oriented lines modeling halfplanes}*/
  typedef typename T::Standard_point_2 Point;
  /*{\Mtypemember the affine points of the plane.}*/
  typedef typename T::Standard_direction_2 Direction;
  /*{\Mtypemember directions in our plane.}*/
  typedef typename T::Standard_aff_transformation_2  Aff_transformation;
  /*{\Mtypemember affine transformations of the plane.}*/

  typedef bool Mark;
  /*{\Xtypemember marking set membership or exclusion.}*/

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  enum Content { EMPTY=0, COMPLETE=1 };
  /*{\Menum construction selection}*/

protected:
  struct AND { bool operator()(bool b1, bool b2)  const { return b1&&b2; }  };
  struct OR { bool operator()(bool b1, bool b2)   const { return b1||b2; }  };
  struct DIFF { bool operator()(bool b1, bool b2) const { return b1&&!b2; } };
  struct XOR { bool operator()(bool b1, bool b2)  const 
                                           { return (b1&&!b2)||(!b1&&b2); } };

  typedef Nef_polyhedron_2_rep<T>           Nef_rep;
  typedef typename Nef_rep::Plane_map       Plane_map;
  typedef typename Nef_rep::Decorator       Decorator;
  typedef typename Nef_rep::Const_decorator Const_decorator;
  typedef typename Nef_rep::Overlayer       Overlayer;
  //typedef typename Nef_rep::T               Transformer;
  typedef typename Nef_rep::Slocator        Slocator;
  typedef typename Nef_rep::Locator         Locator;

  Plane_map& pm() { return ptr->pm_; } 
  const Plane_map& pm() const { return ptr->pm_; } 

  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
      (std::ostream& os, const Nef_polyhedron_2<T>& NP);
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
      (std::istream& is, Nef_polyhedron_2<T>& NP);

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

  typedef std::list<Extended_segment>      ES_list;
  typedef typename ES_list::const_iterator ES_iterator;

  void fill_with_frame_segs(ES_list& L) const
  /*{\Xop fills the list with the four segments which span our frame,
     the convex hull of SW,SE,NW,NE.}*/
  { L.push_back(Extended_segment(EK.SW(),EK.NW()));
    L.push_back(Extended_segment(EK.SW(),EK.SE()));
    L.push_back(Extended_segment(EK.NW(),EK.NE()));
    L.push_back(Extended_segment(EK.SE(),EK.NE()));
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

  void clear_outer_face_cycle_marks() 
  { // unset all frame marks
    Decorator D(pm());
    Face_iterator f = D.faces_begin(); 
    D.mark(f) = false;
    Halfedge_handle e = D.holes_begin(f);
    D.set_marks_in_face_cycle(e, false);
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
    D.mark(++D.faces_begin()) = bool(plane);
  }


  Nef_polyhedron_2(const Line& l, Boundary line = INCLUDED) : Base(Nef_rep())
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the halfplane
  left of |l| including |l| if |line==INCLUDED|, excluding |l| if 
  |line==EXCLUDED|.}*/  
  {   TRACEN("Nconstruction from line "<<l);
    ES_list L;
    fill_with_frame_segs(L);
    Extended_point ep1 = EK.construct_opposite_point(l);
    Extended_point ep2 = EK.construct_point(l);
    L.push_back(EK.construct_segment(ep1,ep2));
    Overlayer D(pm());
    Link_to_iterator I(D, --L.end(), false);
    D.create(L.begin(),L.end(),I);
    CGAL_assertion( I._e != Halfedge_handle() );
    Halfedge_handle el = I._e;
    if ( D.point(D.target(el)) != EK.target(L.back()) )
      el = D.twin(el);
    D.mark(D.face(el)) = true;
    D.mark(el) = bool(line);
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
      D.mark(++D.faces_begin()) = !bool(b); return; }
    CGAL_assertion( I._e != Halfedge_handle() || I._v != Vertex_handle() );
    if ( EK.is_degenerate(L.back()) ) {
      CGAL_assertion(I._v != Vertex_handle());
      D.mark(D.face(I._v)) = !bool(b); D.mark(I._v) = b;
    } else {
      Halfedge_handle el = I._e;
      if ( D.point(D.target(el)) != EK.target(L.back()) )
        el = D.twin(el);  
      D.set_marks_in_face_cycle(el,bool(b));
      if ( D.number_of_faces() > 2 ) D.mark(D.face(el)) = true;
      else                           D.mark(D.face(el)) = !bool(b);
    }
    clear_outer_face_cycle_marks();


  }

  Nef_polyhedron_2(const Nef_polyhedron_2<T>& N1) : Base(N1) {}
  Nef_polyhedron_2& operator=(const Nef_polyhedron_2<T>& N1)
  { Base::operator=(N1); return (*this); }
  ~Nef_polyhedron_2() {}

  protected:
  Nef_polyhedron_2(const Plane_map& H, bool clone=true) : Base(Nef_rep()) 
  /*{\Xcreate makes |\Mvar| a new object.  If |clone==true| then the
  underlying structure of |H| is copied into |\Mvar|.}*/
  { if (clone) {
      Decorator D(pm()); // a decorator working on the rep plane map
      D.clone(H);        // cloning H into pm()
    }
  }
  void clone_rep() { *this = Nef_polyhedron_2<T>(pm()); }

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
    return (D.number_of_vertices()==4 &&
            D.number_of_edges()==4 &&
            D.number_of_faces()==2 &&
            D.mark(++f) == false);
  }

  bool is_plane() const
  /*{\Mop returns true if |\Mvar| is the whole plane, false otherwise.}*/
  { Const_decorator D(pm());
    Face_const_iterator f = D.faces_begin();
    return (D.number_of_vertices()==4 &&
            D.number_of_edges()==4 &&
            D.number_of_faces()==2 &&
            D.mark(++f) == true);
  }

  void extract_complement()
  { TRACEN("extract complement");
    if ( ptr->is_shared() ) clone_rep();
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
  { TRACEN("extract interior");
    if ( ptr->is_shared() ) clone_rep();
    Overlayer D(pm());
    Vertex_iterator v, vend = D.vertices_end();
    for(v = D.vertices_begin(); v != vend; ++v)      D.mark(v) = false;
    Halfedge_iterator e, eend = D.halfedges_end();
    for(e = D.halfedges_begin(); e != eend; ++(++e)) D.mark(e) = false;
    D.simplify();
  }


  void extract_boundary()
  { TRACEN("extract boundary");
    if ( ptr->is_shared() ) clone_rep();
    Overlayer D(pm());
    Vertex_iterator v, vend = D.vertices_end();
    for(v = D.vertices_begin(); v != vend; ++v)      D.mark(v) = true;
    Halfedge_iterator e, eend = D.halfedges_end();
    for(e = D.halfedges_begin(); e != eend; ++(++e)) D.mark(e) = true;
    Face_iterator f, fend = D.faces_end();
    for(f = D.faces_begin(); f != fend; ++f)         D.mark(f) = false;
    clear_outer_face_cycle_marks();
    D.simplify();
  }

  void extract_closure()
  /*{\Xop converts |\Mvar| to its closure. }*/
  { TRACEN("extract closure");
    extract_complement();
    extract_interior();
    extract_complement();
  }

  void extract_regularization()
  /*{\Xop converts |\Mvar| to its regularization. }*/
  { TRACEN("extract regularization");
    extract_interior();
    extract_closure();
  }

  /*{\Mtext \headerline{Constructive Operations}}*/

  Nef_polyhedron_2<T> complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane.}*/
  { Nef_polyhedron_2<T> res = *this;
    res.extract_complement();
    return res;
  }


  Nef_polyhedron_2<T> interior() const
  /*{\Mop returns the interior of |\Mvar|.}*/
  { Nef_polyhedron_2<T> res = *this;
    res.extract_interior();
    return res;
  }

  Nef_polyhedron_2<T> closure() const
  /*{\Mop returns the closure of |\Mvar|.}*/
  { Nef_polyhedron_2<T> res = *this;
    res.extract_closure();
    return res;
  }

  Nef_polyhedron_2<T> boundary() const
  /*{\Mop returns the boundary of |\Mvar|.}*/
  { Nef_polyhedron_2<T> res = *this;
    res.extract_boundary();
    return res;
  }

  Nef_polyhedron_2<T> regularization() const
  /*{\Mop returns the regularized polyhedron (closure of interior).}*/
  { Nef_polyhedron_2<T> res = *this;
    res.extract_regularization();
    return res;
  }


  Nef_polyhedron_2<T> intersection(const Nef_polyhedron_2<T>& N1) const
  /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/
  { Nef_polyhedron_2<T> res(pm(),false); // empty, no frame
    Overlayer PMO(res.pm());
    PMO.subdivide(pm(),N1.pm());
    AND _and; PMO.select(_and);
    res.clear_outer_face_cycle_marks();
    PMO.simplify();
    return res;
  }


  Nef_polyhedron_2<T> join(const Nef_polyhedron_2<T>& N1) const
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/
  { Nef_polyhedron_2<T> res(pm(),false); // empty, no frame
    Overlayer PMO(res.pm());
    PMO.subdivide(pm(),N1.pm());
    OR _or; PMO.select(_or);
    res.clear_outer_face_cycle_marks();
    PMO.simplify();
    return res;
  }

  Nef_polyhedron_2<T> difference(const Nef_polyhedron_2<T>& N1) const
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/
  { Nef_polyhedron_2<T> res(pm(),false); // empty, no frame
    Overlayer PMO(res.pm());
    PMO.subdivide(pm(),N1.pm());
    DIFF _diff; PMO.select(_diff);
    res.clear_outer_face_cycle_marks();
    PMO.simplify();
    return res;
  }    

  Nef_polyhedron_2<T> symmetric_difference(
    const Nef_polyhedron_2<T>& N1) const
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$ 
          |T - \Mvar|. }*/
  { Nef_polyhedron_2<T> res(pm(),false); // empty, no frame
    Overlayer PMO(res.pm());
    PMO.subdivide(pm(),N1.pm());
    XOR _xor; PMO.select(_xor);
    res.clear_outer_face_cycle_marks();
    PMO.simplify();
    return res;
  }

  #if 0
  Nef_polyhedron_2<T> transform(const Aff_transformation& t) const
  /*{\Mop returns $t(|\Mvar|)$.}*/
  { Nef_polyhedron_2<T> res(pm()); // cloned
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

  Nef_polyhedron_2<T>  operator*(const Nef_polyhedron_2<T>& N1) const
  { return intersection(N1); }

  Nef_polyhedron_2<T>  operator+(const Nef_polyhedron_2<T>& N1) const
  { return join(N1); }

  Nef_polyhedron_2<T>  operator-(const Nef_polyhedron_2<T>& N1) const
  { return difference(N1); }

  Nef_polyhedron_2<T>  operator^(const Nef_polyhedron_2<T>& N1) const
  { return symmetric_difference(N1); }

  Nef_polyhedron_2<T>  operator!() const
  { return complement(); }
   
  Nef_polyhedron_2<T>& operator*=(const Nef_polyhedron_2<T>& N1)
  { this = intersection(N1); return *this; }

  Nef_polyhedron_2<T>& operator+=(const Nef_polyhedron_2<T>& N1)
  { this = join(N1); return *this; }

  Nef_polyhedron_2<T>& operator-=(const Nef_polyhedron_2<T>& N1)
  { this = difference(N1); return *this; }

  Nef_polyhedron_2<T>& operator^=(const Nef_polyhedron_2<T>& N1)
  { this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality, respectively.}*/

  bool operator==(const Nef_polyhedron_2<T>& N1) const
  { return symmetric_difference(N1).is_empty(); }

  bool operator!=(const Nef_polyhedron_2<T>& N1) const
  { return !operator==(N1); }  

  bool operator<=(const Nef_polyhedron_2<T>& N1) const
  { return difference(N1).is_empty(); } 

  bool operator<(const Nef_polyhedron_2<T>& N1) const
  { return difference(N1).is_empty() && !N1.difference(*this).is_empty(); } 

  bool operator>=(const Nef_polyhedron_2<T>& N1) const
  { return N1.difference(*this).is_empty(); } 

  bool operator>(const Nef_polyhedron_2<T>& N1) const   
  { return N1.difference(*this).is_empty() && !difference(N1).is_empty(); } 


  /*{\Mtext \headerline{Exploration - Point location - Ray shooting}
  As Nef polyhedra are the result of forming complements 
  and intersections starting from a set |H| of halfspaces that are
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

  void init_locator() const { ptr->init_locator(); }
  const Locator& locator() const 
  { assert(ptr->pl_); return *(ptr->pl_); }


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
      ptr->init_locator();
      Extended_point ep = EK.construct_point(p);
      return locator().locate(ep);
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      Extended_segment s(EK.construct_point(p),
                         PL.point(PL.vertices_begin()));
      return PL.locate(s); 
    }
    CGAL_assertion_msg(0,"location mode not implemented.");
    return Object_handle();
  }

  struct INSET {
    const Const_decorator& D;
    INSET(const Const_decorator& Di) : D(Di) {}
    bool operator()(Vertex_const_handle v) const { return D.mark(v); }
    bool operator()(Halfedge_const_handle e) const { return D.mark(e); }
    bool operator()(Face_const_handle f) const { return D.mark(f); }
  };

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
    if (m == DEFAULT || m == LMWT) {
      ptr->init_locator();
      Extended_point ep = EK.construct_point(p), 
                     eq = EK.construct_point(p,d);
      return locator().ray_shoot(EK.construct_segment(ep,eq),
                                 INSET(locator())); 
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      Extended_point ep = EK.construct_point(p), 
                     eq = EK.construct_point(p,d);
      return PL.ray_shoot(EK.construct_segment(ep,eq),INSET(PL));
    }
    CGAL_assertion_msg(0,"location mode not implemented.");
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
    if (m == DEFAULT || m == LMWT) {
      ptr->init_locator();
      Extended_point ep = EK.construct_point(p), 
                     eq = EK.construct_point(p,d);
      return locator().ray_shoot(EK.construct_segment(ep,eq),INSKEL()); 
    } else if (m == NAIVE) {
      Slocator PL(pm(),EK);
      Extended_point ep = EK.construct_point(p), 
                     eq = EK.construct_point(p,d);
      return PL.ray_shoot(EK.construct_segment(ep,eq),INSKEL()); 
    }
    CGAL_assertion_msg(0,"location mode not implemented.");
    return Object_handle();
  }


  Explorer explorer() const { return Explorer(pm(),EK); }
  /*{\Mop returns a decorator object which allows read-only access of
  the underlying plane map. See the manual page |Explorer| for its 
  usage.}*/


  /*{\Mtext\headerline{Input and Output}
  A Nef polyhedron |\Mvar| can be visualized in a |Window_stream W|. The 
  output operator is defined in the file 
  |CGAL/IO/Nef_\-poly\-hedron_2_\-Win\-dow_\-stream.h|.
  }*/

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

template <typename T>
T Nef_polyhedron_2<T>::EK;


template <typename T>
std::ostream& operator<<
 (std::ostream& os, const Nef_polyhedron_2<T>& NP)
{
  os << "Nef_polyhedron_2<" << NP.EK.output_identifier() << ">\n";
  typedef typename Nef_polyhedron_2<T>::Decorator Decorator;
  CGAL::PM_io_parser<Decorator> O(os, NP.pm()); O.print();
  return os;
}

template <typename T>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_2<T>& NP)
{
  typedef typename Nef_polyhedron_2<T>::Decorator Decorator;
  CGAL::PM_io_parser<Decorator> I(is, NP.pm()); 
  if (I.check_sep("Nef_polyhedron_2<") &&
      I.check_sep(NP.EK.output_identifier()) &&
      I.check_sep(">")) I.read();
  else {
    std::cerr << "Nef_polyhedron_2 input corrupted." << std::endl;
    NP = Nef_polyhedron_2<T>();
  }
  typename Nef_polyhedron_2<T>::Topological_explorer D(NP.explorer());
  D.check_integrity_and_topological_planarity();
  return is;
}


CGAL_END_NAMESPACE

#undef CGAL_SIMPLE_HDS
#endif //CGAL_NEF_POLYHEDRON_2_H


