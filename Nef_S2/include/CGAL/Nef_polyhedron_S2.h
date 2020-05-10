// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYHEDRON_S2_H
#define CGAL_NEF_POLYHEDRON_S2_H

#include <CGAL/license/Nef_S2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_S2/SM_overlayer.h>
#include <CGAL/Modifier_base.h>

#include <vector>
#include <list>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 53
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename K, typename I, typename Mk, typename M> class Nef_polyhedron_S2;
template <typename K, typename I, typename Mk, typename M> class Nef_polyhedron_S2_rep;
template <typename K, typename I, typename Mk> class Nef_polyhedron_3;
class SNC_items;

template <typename K, typename I, typename Mk, typename M>
std::ostream& operator<<(std::ostream&, const Nef_polyhedron_S2<K,I,Mk,M>&);
template <typename K, typename I, typename Mk, typename M>
std::istream& operator>>(std::istream&, Nef_polyhedron_S2<K,I,Mk,M>&);


template <typename K, typename I, typename Mk, typename M>
class Nef_polyhedron_S2_rep {

  typedef Nef_polyhedron_S2_rep<K,I,Mk,M>        Self;
  friend class Nef_polyhedron_S2<K,I,Mk,M>;

 public:
  typedef CGAL::Sphere_geometry<K>                     Sphere_kernel;
  typedef Mk                                           Mark;
  typedef M                                            Sphere_map;
  typedef CGAL::SM_const_decorator<Sphere_map>         Const_decorator;
  typedef CGAL::SM_decorator<Sphere_map>               Decorator;
  typedef CGAL::SM_overlayer<Decorator>                Overlayer;
  typedef CGAL::SM_point_locator<Const_decorator>      Locator;

 private:
  Sphere_map sm_;

public:
  Nef_polyhedron_S2_rep() : sm_() {}
  Nef_polyhedron_S2_rep(const Self&) : sm_() {}
  ~Nef_polyhedron_S2_rep() { sm_.clear(); }
};

/*{\Moptions print_title=yes }*/
/*{\Manpage {Nef_polyhedron_S2}{K}
{Nef Polyhedra in the sphere surface}{N}}*/

/*{\Mdefinition An instance of data type |\Mname| is a subset of $S_2$
that is the result of forming complements and intersections starting
from a finite set |H| of half-spaces. |\Mtype| is closed under all
binary set operations |intersection|, |union|, |difference|,
|complement| and under the topological operations |boundary|,
|closure|, and |interior|.

The template parameter |Kernel| is specified via a kernel concept.
|Kernel| must be a model of the concept |NefSphereKernelTraits_2|.
}*/

template <typename Kernel_, typename Items_ = SM_items, typename Mark_ = bool,
          typename Map_ = Sphere_map<Sphere_geometry<Kernel_>,Items_, Mark_> >
class Nef_polyhedron_S2 : public Handle_for< Nef_polyhedron_S2_rep<Kernel_,Items_,Mark_,Map_> >,
                          public Nef_polyhedron_S2_rep<Kernel_,Items_,Mark_,Map_>::Const_decorator {

  using Nef_polyhedron_S2_rep<Kernel_,Items_,Mark_,Map_>::Const_decorator::set_sm;

public:
  /*{\Mtypes 7}*/
  typedef Items_                                              Items;
  typedef Kernel_                                             Kernel;
  typedef Map_                                                Sphere_map;
  typedef Mark_                                               Mark;
  typedef Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>     Self;
  typedef Nef_polyhedron_S2_rep<Kernel,Items,Mark,Sphere_map> Rep;
  typedef Handle_for< Nef_polyhedron_S2_rep<Kernel,Items,Mark,Sphere_map> >  Base;
  typedef typename Rep::Sphere_kernel                         Sphere_kernel;
//  typedef typename Rep::Sphere_map                            Sphere_map;

  typedef typename Sphere_kernel::Sphere_point   Sphere_point;
  /*{\Mtypemember points in the sphere surface.}*/

  typedef typename Sphere_kernel::Sphere_segment Sphere_segment;
  /*{\Mtypemember segments in the sphere surface.}*/

  typedef typename Sphere_kernel::Sphere_circle  Sphere_circle;
  /*{\Mtypemember oriented great circles modeling spherical half-spaces}*/

  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;


//  typedef typename Rep::Mark Mark;
  /*{\Xtypemember marking set membership or exclusion.}*/

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  /*{\Menum construction selection.}*/

  enum Content { EMPTY=0, COMPLETE=1 };
  /*{\Menum construction selection}*/

  const Sphere_map& sphere_map() const { return this->ptr()->sm_; }
protected:
  Sphere_map& sphere_map() { return this->ptr()->sm_; }

  struct AND { bool operator()(const Mark& b1, const Mark& b2)  const { return b1&&b2; }  };
  struct OR { bool operator()(const Mark& b1, const Mark& b2)   const { return b1||b2; }  };
  struct DIFF { bool operator()(const Mark& b1, const Mark& b2) const { return b1&&!b2; } };
  struct XOR { bool operator()(const Mark& b1, const Mark& b2)  const
               { return (b1&&!b2)||(!b1&&b2); } };

  typedef Nef_polyhedron_S2_rep<Kernel,Items,Mark,Sphere_map>  Nef_rep;
  typedef typename Nef_rep::Decorator                     Decorator;
public:
  typedef typename Nef_rep::Const_decorator               Const_decorator;
protected:
  typedef typename Nef_rep::Overlayer                     Overlayer;
  typedef typename Nef_rep::Locator                       Locator;

  friend std::ostream& operator<< <>
      (std::ostream& os, const Self& NP);
  friend std::istream& operator>> <>
      (std::istream& is, Self& NP);

public:
  typedef typename Decorator::SVertex_handle         SVertex_handle;
  typedef typename Decorator::SHalfedge_handle       SHalfedge_handle;
  typedef typename Decorator::SHalfloop_handle       SHalfloop_handle;
  typedef typename Decorator::SFace_handle           SFace_handle;

  typedef typename Sphere_map::SVertex_base          SVertex;
  typedef typename Sphere_map::SHalfedge_base        SHalfedge;
  typedef typename Sphere_map::SHalfloop             SHalfloop;
  typedef typename Sphere_map::SFace_base            SFace;

  typedef typename Decorator::SVertex_const_handle   SVertex_const_handle;
  typedef typename Decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Decorator::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Decorator::SFace_const_handle     SFace_const_handle;

  typedef typename Decorator::SVertex_iterator       SVertex_iterator;
  typedef typename Decorator::SHalfedge_iterator     SHalfedge_iterator;
  typedef typename Decorator::SHalfloop_iterator     SHalfloop_iterator;
  typedef typename Decorator::SFace_iterator         SFace_iterator;

  typedef typename Const_decorator::SVertex_const_iterator
                                                    SVertex_const_iterator;
  typedef typename Const_decorator::SHalfedge_const_iterator
                                                    SHalfedge_const_iterator;
  typedef typename Const_decorator::SHalfloop_const_iterator
                                                    SHalfloop_const_iterator;
  typedef typename Const_decorator::SFace_const_iterator
                                                    SFace_const_iterator;
  typedef typename Const_decorator::Size_type Size_type;
  typedef Size_type size_type;

  typedef std::list<Sphere_segment>  SS_list;
  typedef typename SS_list::const_iterator SS_iterator;

  friend class Nef_polyhedron_3<Kernel, SNC_items, Mark>;

public:
  /*{\Mcreation 3}*/

  Nef_polyhedron_S2(Content sphere = EMPTY) : Base(Nef_rep())
  /*{\Mcreate creates an instance |\Mvar| of type |\Mname|
  and initializes it to the empty set if |sphere == EMPTY|
  and to the whole sphere if |sphere == COMPLETE|.}*/
  {
    set_sm(&sphere_map());
    Decorator D(&sphere_map());
    SFace_handle sf=D.new_sface();
    sf->mark() = bool(sphere);
  }


  Nef_polyhedron_S2(const Sphere_circle& c,
                    Boundary circle = INCLUDED) : Base(Nef_rep()) {
  /*{\Mcreate creates a Nef polyhedron |\Mvar| containing the half-sphere
  left of |c| including |c| if |circle==INCLUDED|, excluding |c| if
  |circle==EXCLUDED|.}*/

    set_sm(&sphere_map());
    CGAL_NEF_TRACEN("Nef_polyhedron_S2(): construction from circle "<<c);
    Decorator D(&sphere_map());
    Overlayer O(&sphere_map());
    O.create(c);
    SHalfloop_handle h = D.shalfloop();
    if ( h->circle() != c ) h = h->twin();
    h->incident_sface()->mark() = true;
    h->mark() = h->twin()->mark() = bool(circle);
  }


  template <class Forward_iterator>
  Nef_polyhedron_S2(Forward_iterator first, Forward_iterator beyond,
    Boundary b = INCLUDED) : Base(Nef_rep())
  /*{\Mcreate creates a Nef polyhedron |\Mvar| from the set of sphere
    segments in the iterator range |[first,beyond)|. If the set of sphere
    segments is a simple polygon that separates the sphere surface
    into two regions, then the polygonal region that is left of the
    segment |*first| is selected. The polygonal region includes its
    boundary if |b = INCLUDED| and excludes the boundary
    otherwise. |Forward_iterator| has to be an iterator with value
    type |Sphere_segment|.}*/
  { CGAL_NEF_TRACEN("Nef_polyhedron_S2(): creation from segment range");
    CGAL_assertion(first!=beyond);
    set_sm(&sphere_map());
    Overlayer D(&sphere_map());
    Sphere_segment s = *first;
    D.create_from_segments(first,beyond);
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,D) {
      Sphere_circle c(e->circle());
      if ( c == s.sphere_circle() ) break;
    }
    if ( e != SHalfedge_iterator() ) {
      if ( e->circle() != s.sphere_circle() ) e = e->twin();
      CGAL_assertion( e->circle() == s.sphere_circle() );
      D.set_marks_in_face_cycle(e,bool(b));
      if ( D.number_of_sfaces() > 2 ) e->incident_sface()->mark() = true;
      else                            e->incident_sface()->mark() = !bool(b);
      return;
    }
    D.simplify();
  }

  Nef_polyhedron_S2(const Self& N1) : Base(N1), Const_decorator() {
    set_sm(&sphere_map());
  }
  Nef_polyhedron_S2& operator=(const Self& N1)
  { Base::operator=(N1); set_sm(&sphere_map()); return (*this); }
  ~Nef_polyhedron_S2() {}

  template <class Forward_iterator>
  Nef_polyhedron_S2(Forward_iterator first, Forward_iterator beyond,
    double p) : Base(Nef_rep())
  /*{\Xcreate creates a random Nef polyhedron from the arrangement of
  the set of circles |S = set[first,beyond)|. The cells of the arrangement
  are selected uniformly at random with probability $p$. \precond $0 < p
  < 1$.}*/
  { CGAL_assertion(0<=p && p<=1);
    CGAL_assertion(first!=beyond);
    set_sm(&sphere_map());
    Overlayer D(&sphere_map());
    D.create_from_circles(first, beyond); D.simplify();

    boost::rand48 rng;
    boost::uniform_real<> dist(0,1);
    boost::variate_generator<boost::rand48&, boost::uniform_real<> > get_double(rng,dist);

    SVertex_iterator v; SHalfedge_iterator e; SFace_iterator f;
    CGAL_forall_svertices(v,D)
      v->mark() = ( get_double() < p ? true : false );
    CGAL_forall_shalfedges(e,D)
      e->mark() = ( get_double() < p ? true : false );
    CGAL_forall_sfaces(f,D)
      f->mark() = ( get_double() < p ? true : false );
    D.simplify();
  }

 void delegate( Modifier_base<Sphere_map>& modifier) {
   // calls the `operator()' of the `modifier'. Precondition: The
   // `modifier' returns a consistent representation.
   modifier(sphere_map());
   //   CGAL_expensive_postcondition( is_valid());
 }

//protected:
  Nef_polyhedron_S2(const Sphere_map& H, bool clone=true) : Base(Nef_rep())
  /*{\Xcreate makes |\Mvar| a new object.  If |clone==true| then the
  underlying structure of |H| is copied into |\Mvar|.}*/
{
    if(clone)
      this->ptr()->sm_ = H;
    set_sm(&sphere_map());
  }

  void clone_rep() { *this = Self(sphere_map()); }

  /*{\Moperations 4 3 }*/
  public:

  void clear(Content plane = EMPTY)
  { *this = Nef_polyhedron_S2(plane); }
  /*{\Mop makes |\Mvar| the empty set if |plane == EMPTY| and the
  full plane if |plane == COMPLETE|.}*/

  bool is_empty() const
  /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
  { Const_decorator D(&sphere_map());
    CGAL_NEF_TRACEN("is_empty()"<<*this);
    SFace_const_iterator f = D.sfaces_begin();
    return (D.number_of_svertices()==0 &&
            D.number_of_sedges()==0 &&
            D.number_of_sloops()==0 &&
            D.number_of_sfaces()==1 &&
            f->mark() == false);
  }

  bool is_plane() const
  /*{\Mop returns true if |\Mvar| is the whole plane, false otherwise.}*/
  { Const_decorator D(&sphere_map());
    SFace_const_iterator f = D.sfaces_begin();
    return (D.number_of_svertices()==0 &&
            D.number_of_sedges()==0 &&
            D.number_of_sloops()==0 &&
            D.number_of_sfaces()==1 &&
            f->mark() == true);
  }

  bool is_sphere() const
  {
    return is_plane();
  }

  void extract_complement()
  { CGAL_NEF_TRACEN("extract complement");
    if ( this->is_shared() ) clone_rep();
    Overlayer D(&sphere_map());
    SVertex_iterator v;
    SHalfedge_iterator e;
    SFace_iterator f;
    CGAL_forall_svertices(v,D) v->mark() = !v->mark();
    CGAL_forall_sedges(e,D) e->mark() = !e->mark();
    CGAL_forall_sfaces(f,D) f->mark() = !f->mark();

    if ( D.has_shalfloop() )
      D.shalfloop()->mark() =
        D.shalfloop()->twin()->mark() =
        !D.shalfloop()->mark();
  }

  void extract_interior()
  { CGAL_NEF_TRACEN("extract interior");
    if ( this->is_shared() ) clone_rep();
    Overlayer D(&sphere_map());
    SVertex_iterator v;
    SHalfedge_iterator e;
    CGAL_forall_svertices(v,D) v->mark() = false;
    CGAL_forall_sedges(e,D) e->mark() = false;
    if ( D.has_shalfloop() ) D.shalfloop()->mark() = false;
    D.simplify();
  }


  void extract_boundary()
  { CGAL_NEF_TRACEN("extract boundary");
    if ( this->is_shared() ) clone_rep();
    Overlayer D(&sphere_map());
    SVertex_iterator v;
    SHalfedge_iterator e;
    SFace_iterator f;
    CGAL_forall_svertices(v,D) v->mark() = true;
    CGAL_forall_sedges(e,D)    e->mark() = true;
    CGAL_forall_sfaces(f,D)    f->mark() = false;
    if ( D.has_shalfloop() )       D.shalfloop()->mark() = D.shalfoop()->twin()->mark() = true;
    D.simplify();
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

  Self complement() const
  /*{\Mop returns the complement of |\Mvar| in the plane.}*/
  { Self res = *this;
    res.extract_complement();
    return res;
  }


  Self interior() const
  /*{\Mop returns the interior of |\Mvar|.}*/
  { Self res = *this;
    res.extract_interior();
    return res;
  }

  Self closure() const
  /*{\Mop returns the closure of |\Mvar|.}*/
  { Self res = *this;
    res.extract_closure();
    return res;
  }

  Self boundary() const
  /*{\Mop returns the boundary of |\Mvar|.}*/
  { Self res = *this;
    res.extract_boundary();
    return res;
  }

  Self regularization() const
  /*{\Mop returns the regularized polyhedron (closure of interior).}*/
  { Self res = *this;
    res.extract_regularization();
    return res;
  }


  Self intersection(const Self& N1) const
  /*{\Mop returns |\Mvar| $\cap$ |N1|. }*/
  { Self res(sphere_map(),false); // empty
    Overlayer D(&res.sphere_map());
    D.subdivide(&sphere_map(),&N1.sphere_map());
    AND _and; D.select(_and); D.simplify();
    return res;
  }


  Self join(const Self& N1) const
  /*{\Mop returns |\Mvar| $\cup$ |N1|. }*/
  { Self res(sphere_map(),false); // empty
    Overlayer D(&res.sphere_map());
    D.subdivide(&sphere_map(),&N1.sphere_map());
    OR _or; D.select(_or); D.simplify();
    return res;
  }

  Self difference(const Self& N1) const
  /*{\Mop returns |\Mvar| $-$ |N1|. }*/
  { Self res(sphere_map(),false); // empty
    Overlayer D(&res.sphere_map());
    D.subdivide(&sphere_map(),&N1.sphere_map());
    DIFF _diff; D.select(_diff); D.simplify();
    return res;
  }

  Self symmetric_difference(
    const Self& N1) const
  /*{\Mop returns the symmectric difference |\Mvar - T| $\cup$
          |T - \Mvar|. }*/
  { Self res(sphere_map(),false); // empty
    Overlayer D(&res.sphere_map());
    D.subdivide(&sphere_map(),&N1.sphere_map());
    XOR _xor; D.select(_xor); D.simplify();
    return res;
  }

  /*{\Mtext Additionally there are operators |*,+,-,^,!| which
  implement the binary operations \emph{intersection}, \emph{union},
  \emph{difference}, \emph{symmetric difference}, and the unary
  operation \emph{complement} respectively. There are also the
  corresponding modification operations |*=,+=,-=,^=|.}*/

  Self  operator*(const Self& N1) const
  { return intersection(N1); }

  Self  operator+(const Self& N1) const
  { return join(N1); }

  Self  operator-(const Self& N1) const
  { return difference(N1); }

  Self  operator^(const Self& N1) const
  { return symmetric_difference(N1); }

  Self  operator!() const
  { return complement(); }

  Self& operator*=(const Self& N1)
  { *this = intersection(N1); return *this; }

  Self& operator+=(const Self& N1)
  { *this = join(N1); return *this; }

  Self& operator-=(const Self& N1)
  { *this = difference(N1); return *this; }

  Self& operator^=(const Self& N1)
  { *this = symmetric_difference(N1); return *this; }

  /*{\Mtext There are also comparison operations like |<,<=,>,>=,==,!=|
  which implement the relations subset, subset or equal, superset, superset
  or equal, equality, inequality, respectively.}*/

  bool operator==(const Self& N1) const
  { return symmetric_difference(N1).is_empty(); }

  bool operator!=(const Self& N1) const
  { return !operator==(N1); }

  bool operator<=(const Self& N1) const
  { return difference(N1).is_empty(); }

  bool operator<(const Self& N1) const
  { return difference(N1).is_empty() && !N1.difference(*this).is_empty(); }

  bool operator>=(const Self& N1) const
  { return N1.difference(*this).is_empty(); }

  bool operator>(const Self& N1) const
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

  typedef Const_decorator Explorer;
  /*{\Mtypemember a decorator to examine the underlying plane map.
  See the manual page of |Explorer|.}*/

  typedef typename Locator::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of the underlying
  plane map. The kind of object |(vertex, halfedge, face)| can
  be determined and the object can be assigned to a corresponding
  handle by the three functions:\\
  |bool assign(SVertex_const_handle& h, Object_handle)|\\
  |bool assign(SHalfedge_const_handle& h, Object_handle)|\\
  |bool assign(SFace_const_handle& h, Object_handle)|\\
  where each function returns |true| iff the assignment to
  |h| was done.}*/


  /*{\Moperations 3 1 }*/

  bool contains(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the set
  represented by |\Mvar|.}*/
  { Locator PL(&sphere_map()); return PL.mark(h); }

  bool contained_in_boundary(Object_handle h) const
  /*{\Mop  returns true iff the object |h| is contained in the $1$-skeleton
  of |\Mvar|.}*/
  { SVertex_const_handle v;
    SHalfedge_const_handle e;
    return  ( CGAL::assign(v,h) || CGAL::assign(e,h) );
  }


  Object_handle locate(const Sphere_point& p) const
  /*{\Mop  returns a generic handle |h| to an object (face, halfedge, vertex)
  of the underlying plane map that contains the point |p| in its relative
  interior. The point |p| is contained in the set represented by |\Mvar| if
  |\Mvar.contains(h)| is true. The location mode flag |m| allows one to choose
  between different point location strategies.}*/
  {
    Locator PL(&sphere_map());
    return PL.locate(p);
  }

  struct INSET {
    const Const_decorator& D;
    INSET(const Const_decorator& Di) : D(Di) {}
    bool operator()(SVertex_const_handle v) const { return v->mark(); }
    bool operator()(SHalfedge_const_handle e) const { return e->mark(); }
    bool operator()(SHalfloop_const_handle l) const { return l->mark(); }
    bool operator()(SFace_const_handle f) const { return f->mark(); }
  };

  Object_handle ray_shoot(const Sphere_point& p,
                          const Sphere_direction& d) const
  /*{\Mop returns a handle |h| with |\Mvar.contains(h)| that can be
  converted to a |SVertex_/SHalfedge_/SFace_const_handle| as described
  above. The object returned is intersected by the ray starting in |p|
  with direction |d| and has minimal distance to |p|.  The operation
  returns the null handle |nullptr| if the ray shoot along |d| does not hit
  any object |h| of |\Mvar| with |\Mvar.contains(h)|.}*/
  {
    Locator PL(&sphere_map());
    return PL.ray_shoot(p,d,INSET(PL));
  }

  struct INSKEL {
    bool operator()(SVertex_const_handle) const { return true; }
    bool operator()(SHalfedge_const_handle) const { return true; }
    bool operator()(SHalfloop_const_handle) const { return true; }
    bool operator()(SFace_const_handle) const { return false; }
  };

  Object_handle ray_shoot_to_boundary(const Sphere_point& p,
                                      const Sphere_direction& d) const
  /*{\Mop returns a handle |h| that can be converted to a
  |SVertex_/SHalfedge_const_handle| as described above. The object
  returned is part of the $1$-skeleton of |\Mvar|, intersected by the
  ray starting in |p| with direction |d| and has minimal distance to
  |p|.  The operation returns the null handle |nullptr| if the ray shoot
  along |d| does not hit any $1$-skeleton object |h| of |\Mvar|. The
  location mode flag |m| allows one to choose between different point
  location strategies.}*/
  {
    Locator PL(&sphere_map());
    return PL.ray_shoot(p,d,INSKEL());
  }


  //  Explorer explorer() const
  /*{\Mop returns a decorator object which allows read-only access of
  the underlying plane map. See the manual page |Explorer| for its
  usage.}*/
  //  { return Explorer(const_cast<Sphere_map*>(&sphere_map())); }

  /*{\Mtext\headerline{Input and Output}
  A Nef polyhedron |\Mvar| can be visualized in an open GL window. The
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
  the naive way. The operations run in linear query time without
  any preprocessing.}*/

  /*{\Mexample Nef polyhedra are parameterized by a standard CGAL
  kernel.

  \begin{Mverb}
  #include <CGAL/Homogeneous.h>
  #include <CGAL/leda_integer.h>
  #include <CGAL/Nef_polyhedron_S2.h>

  using namespace CGAL;
  typedef  Homogeneous<leda_integer>   Kernel;
  typedef  SM_items<Kernel, bool>      SM_items;
  typedef  Nef_polyhedron_S2<SM_items> Nef_polyhedron;
  typedef  Nef_polyhedron::Sphere_circle Sphere_circle;

  int main()
  {
    Nef_polyhedron N1(Sphere_circle(1,0,0));
    Nef_polyhedron N2(Sphere_circle(0,1,0), Nef_polyhedron::EXCLUDED);
    Nef_polyhedron N3 = N1 * N2; // line (*)
    return 0;
  }
  \end{Mverb}
  After line (*) |N3| is the intersection of |N1| and |N2|.}*/


}; // end of Nef_polyhedron_S2

template <typename Kernel,typename Items,typename Mark, typename Sphere_map>
std::ostream& operator<<
 (std::ostream& os, const Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>& NP)
{
  os << "Nef_polyhedron_S2\n";
  typedef typename Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>::Explorer Decorator;
  CGAL::SM_io_parser<Decorator> O(os, Decorator(&NP.sphere_map()));
  O.print();
  return os;
}

template <typename Kernel,typename Items, typename Mark, typename Sphere_map>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>& NP)
{
  typedef typename Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>::Decorator Decorator;
  CGAL::SM_io_parser<Decorator> I(is, Decorator(&NP.sphere_map()));
  //  if ( I.check_sep("Nef_polyhedron_S2") )
  I.read();
  /*
  else {
    std::cerr << "Nef_polyhedron_S2 input corrupted." << std::endl;
    NP = Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>();
  }
  */
  /*
  typename Nef_polyhedron_S2<Kernel,Items,Mark,Sphere_map>::Topological_explorer D(NP.explorer());
  D.check_integrity_and_topological_planarity();
  */
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_NEF_POLYHEDRON_S2_H
