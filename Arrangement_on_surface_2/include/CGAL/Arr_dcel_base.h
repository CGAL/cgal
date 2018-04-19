// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 (based on old version by: Iddo Hanniel and Oren Nechushtan)

#ifndef CGAL_ARR_DCEL_BASE_H
#define CGAL_ARR_DCEL_BASE_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The definition of the base DCEL class for planar arrangements and its
 * peripheral records.
 */

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <list>
#include <map>
#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/In_place_list.h>
#include <CGAL/function_objects.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Arrangement_2/Arrangement_2_iterators.h>
#include <CGAL/assertions.h>


namespace CGAL {

inline void* _clean_pointer(const void* p)
{
  CGAL_static_assertion(sizeof(void*) == sizeof(size_t));
  const size_t  mask = ~1;
  const size_t  val = (reinterpret_cast<size_t>(p) & mask);

  return (reinterpret_cast<void*>(val));
}

inline void* _set_lsb(const void* p)
{
  const size_t  mask = 1;
  const size_t  val = (reinterpret_cast<size_t>(p) | mask);
  return (reinterpret_cast<void*>( val));
}

inline bool _is_lsb_set(const void* p)
{
  const size_t  mask = 1;
  const size_t  val = reinterpret_cast<size_t>(p);
  return ((val & mask) != 0);
}

/*! \class
 * Base vertex class.
 */
template <typename Point_> class Arr_vertex_base {
public:
  typedef Point_       Point;

  /*! \struct
   * An auxiliary structure for rebinding the vertex with a new point class.
   */
  template<typename PNT> struct rebind { typedef Arr_vertex_base<PNT> other; };

protected:
  void* p_inc;  // An incident halfedge pointing at the vertex,
                // or the isolated vertex information (in case it is
                // isolated). The LSB of the pointer indicates whether
                // the vertex is isolated.
  Point* p_pt;  // The point associated with the vertex.
  char pss[2];  // The x and y parameter spaces (condensed in two bytes).

public:
  /*! Default constructor. */
  Arr_vertex_base() :
    p_inc(NULL),
    p_pt(NULL)
  { pss[0] = pss[1] = static_cast<char>(CGAL::ARR_INTERIOR); }

  /*! Destructor. */
  virtual ~Arr_vertex_base() {}

  /*! Check if the point pointer is NULL. */
  bool has_null_point() const { return (p_pt == NULL); }

  /*! Get the point (const version). */
  const Point& point() const
  {
    CGAL_assertion(p_pt != NULL);
    return (*p_pt);
  }

  /*! Get the point (non-const version). */
  Point& point()
  {
    CGAL_assertion(p_pt != NULL);
    return (*p_pt);
  }

  /*! Set the point (may be a NULL point). */
  void set_point(Point* p) { p_pt = p; }

  /*! Get the boundary type in x. */
  Arr_parameter_space parameter_space_in_x() const
  { return (Arr_parameter_space(pss[0])); }

  /*! Get the boundary type in y. */
  Arr_parameter_space parameter_space_in_y() const
  { return (Arr_parameter_space(pss[1])); }

  /*! Set the boundary conditions of the vertex. */
  void set_boundary(Arr_parameter_space ps_x, Arr_parameter_space ps_y)
  {
    pss[0] = static_cast<char>(ps_x);
    pss[1] = static_cast<char>(ps_y);
    return;
  }

  /*! Assign from another vertex. */
  virtual void assign(const Arr_vertex_base<Point>& v)
  {
    p_pt = v.p_pt;
    pss[0] = v.pss[0];
    pss[1] = v.pss[1];
  }
};

/*! \class
 * Base halfedge class.
 */
template <typename X_monotone_curve_> class Arr_halfedge_base {
public:
  typedef X_monotone_curve_  X_monotone_curve;

  /*! \struct
   * An auxiliary structure for rebinding the halfedge with a new curve class.
   */
  template<typename XCV>
  struct rebind { typedef Arr_halfedge_base<XCV> other; };

protected:
  void* p_opp;  // The opposite halfedge.
  void* p_prev; // The previous halfedge in the component boundary.
  void* p_next; // The next halfedge in the component boundary.

  void* p_v;    // The incident vertex (the target of the halfedge).
                // The LSB of this pointer is used to store the
                // direction of the halfedge.
  void* p_comp; // The component this halfedge belongs to: the incident
                // face for outer CCBs and the inner CCB information for
                // inner CCBs. The LSB of the pointer indicates whether
                // the halfedge lies on the boundary of an inner CCB.

  X_monotone_curve* p_cv; // The associated x-monotone curve.

public:
  /*! Default constructor */
  Arr_halfedge_base() :
    p_opp(NULL),
    p_prev(NULL),
    p_next(NULL),
    p_v(NULL),
    p_comp(NULL),
    p_cv(NULL)
  {}

  /*! Destructor. */
  virtual ~Arr_halfedge_base() {}

  /*! Check if the curve pointer is NULL. */
  bool has_null_curve() const { return (p_cv == NULL); }

  /*! Get the x-monotone curve (const version). */
  const X_monotone_curve& curve() const
  {
    CGAL_precondition(p_cv != NULL);
    return (*p_cv);
  }

  /*! Get the x-monotone curve (non-const version). */
  X_monotone_curve& curve()
  {
    CGAL_precondition(p_cv != NULL);
    return (*p_cv);
  }

  /*! Set the x-monotone curve. */
  void set_curve(X_monotone_curve* c)
  {
    p_cv = c;

    // Set the curve for the opposite halfedge as well.
    Arr_halfedge_base<X_monotone_curve>* opp =
      reinterpret_cast<Arr_halfedge_base<X_monotone_curve>* >(p_opp);

    opp->p_cv = c;
  }

  /*! Assign from another halfedge. */
  virtual void assign(const Arr_halfedge_base<X_monotone_curve>& he)
  { p_cv = he.p_cv; }
};

/*!
 * Base face class.
 */
class Arr_face_base
{
public:
  typedef std::list<void*>                      Outer_ccbs_container;
  typedef Outer_ccbs_container::iterator        Outer_ccb_iterator;
  typedef Outer_ccbs_container::const_iterator  Outer_ccb_const_iterator;

  typedef std::list<void*>                      Inner_ccbs_container;
  typedef Inner_ccbs_container::iterator        Inner_ccb_iterator;
  typedef Inner_ccbs_container::const_iterator  Inner_ccb_const_iterator;

  typedef std::list<void*>                      Isolated_vertices_container;
  typedef Isolated_vertices_container::iterator Isolated_vertex_iterator;
  typedef Isolated_vertices_container::const_iterator
                                                Isolated_vertex_const_iterator;

protected:
  enum {
    IS_UNBOUNDED = 1,
    IS_FICTITIOUS = 2
  };

  int                          flags;      // Face flags.
  Outer_ccbs_container         outer_ccbs; // The outer CCBs of the faces.
  Inner_ccbs_container         inner_ccbs; // The inner CCBs of the face.
  Isolated_vertices_container  iso_verts;  // The isolated vertices inside
                                           // the face.
public:
  /*! Default constructor. */
  Arr_face_base() : flags(0) {}

  /*! Destructor. */
  virtual ~Arr_face_base() {}

  /*! Check if the face is unbounded. */
  bool is_unbounded() const { return ((flags & IS_UNBOUNDED) != 0); }

  /*! Set the face as bounded or unbounded. */
  void set_unbounded(bool unbounded)
  { flags = (unbounded) ? (flags | IS_UNBOUNDED) : (flags & ~IS_UNBOUNDED); }

  /*! Check if the face is fictitious. */
  bool is_fictitious() const { return ((flags & IS_FICTITIOUS) != 0); }

  /*! Set the face as fictitious or valid. */
  void set_fictitious(bool fictitious)
  { flags = (fictitious) ? (flags | IS_FICTITIOUS) : (flags & ~IS_FICTITIOUS); }

  /*! Assign from another face. */
  virtual void assign(const Arr_face_base& f) { flags = f.flags; }
};

// Forward declarations:
template <class V, class H, class F> class Arr_vertex;
template <class V, class H, class F> class Arr_halfedge;
template <class V, class H, class F> class Arr_face;
template <class V, class H, class F> class Arr_outer_ccb;
template <class V, class H, class F> class Arr_inner_ccb;
template <class V, class H, class F> class Arr_isolated_vertex;

/*! \class
 * The default arrangement DCEL vertex class.
 */
template <class V, class H, class F>
class Arr_vertex : public V, public In_place_list_base<Arr_vertex<V,H,F> >
{
public:

  typedef V                           Base;
  typedef Arr_vertex<V,H,F>           Vertex;
  typedef Arr_halfedge<V,H,F>         Halfedge;
  typedef Arr_isolated_vertex<V,H,F>  Isolated_vertex;

  /*! Default constructor. */
  Arr_vertex() {}

  /*! Check if the vertex is isolated. */
  bool is_isolated() const
  {
    // Note that we use the LSB of the p_inc pointer as a Boolean flag.
    return (_is_lsb_set(this->p_inc));
  }

  /*! Get an incident halfedge (const version). */
  const Halfedge* halfedge() const
  {
    CGAL_precondition(! is_isolated());
    return (reinterpret_cast<const Halfedge*>(this->p_inc));
  }

  /*! Get an incident halfedge (non-const version). */
  Halfedge* halfedge()
  {
    CGAL_precondition(! is_isolated());
    return (reinterpret_cast<Halfedge*>(this->p_inc));
  }

  /*! Set an incident halfedge (for non-isolated vertices). */
  void set_halfedge(Halfedge* he)
  {
    // Set the halfedge pointer and reset the LSB.
    this->p_inc = he;
  }

  /*! Get the isolated vertex information (const version). */
  const Isolated_vertex* isolated_vertex() const
  {
    CGAL_precondition(is_isolated());
    return (reinterpret_cast<const Isolated_vertex*>(_clean_pointer
                                                     (this->p_inc)));
  }

  /*! Get the isolated vertex information (non-const version). */
  Isolated_vertex* isolated_vertex()
  {
    CGAL_precondition(is_isolated());
    return (reinterpret_cast<Isolated_vertex*>(_clean_pointer(this->p_inc)));
  }

  /*! Set the isolated vertex information. */
  void set_isolated_vertex(Isolated_vertex* iv)
  {
    // Set the isolated vertex-information pointer and set its LSB.
    this->p_inc = _set_lsb(iv);
  }
};

/*! \class
 * The default arrangement DCEL halfedge class.
 */
template <class V, class H, class F>
class Arr_halfedge : public H,
                     public In_place_list_base<Arr_halfedge<V,H,F> >
{
public:
  typedef H                     Base;
  typedef Arr_vertex<V,H,F>     Vertex;
  typedef Arr_halfedge<V,H,F>   Halfedge;
  typedef Arr_face<V,H,F>       Face;
  typedef Arr_outer_ccb<V,H,F>  Outer_ccb;
  typedef Arr_inner_ccb<V,H,F>  Inner_ccb;

  /*! Default constructor. */
  Arr_halfedge() {}

  /*! Get the opposite halfedge (const version). */
  const Halfedge* opposite () const
  { return (reinterpret_cast<const Halfedge*>(this->p_opp)); }

  /*! Get the opposite halfedge (non-const version). */
  Halfedge* opposite() { return (reinterpret_cast<Halfedge*>(this->p_opp)); }

  /*! Sets the opposite halfedge. */
  void set_opposite(Halfedge* he) { this->p_opp = he; }

  /*! Get the direction of the halfedge. */
  Arr_halfedge_direction direction() const
  {
    // Note that we use the LSB of the p_v pointer as a Boolean flag.
    if (_is_lsb_set(this->p_v)) return (ARR_LEFT_TO_RIGHT);
    else return (ARR_RIGHT_TO_LEFT);
  }

  /*! Set the direction of the edge (and of its opposite halfedge). */
  void set_direction(Arr_halfedge_direction dir)
  {
    Halfedge* opp = reinterpret_cast<Halfedge*>(this->p_opp);

    if (dir == ARR_LEFT_TO_RIGHT) {
      this->p_v = _set_lsb(this->p_v);
      opp->p_v = _clean_pointer(opp->p_v);
    }
    else {
      this->p_v = _clean_pointer(this->p_v);
      opp->p_v = _set_lsb(opp->p_v);
    }
  }

  /*! Get the previous halfedge along the chain (const version). */
  const Halfedge* prev() const
  { return (reinterpret_cast<const Halfedge*>(this->p_prev)); }

  /*! Get the previous halfedge along the chain (const version). */
  Halfedge* prev() { return (reinterpret_cast<Halfedge*>(this->p_prev)); }

  /*! Set the previous halfedge along the chain. */
  void set_prev(Halfedge* he)
  {
    this->p_prev = he;
    he->p_next = this;
  }

  /*! Get the next halfedge along the chain (const version). */
  const Halfedge* next() const
  { return (reinterpret_cast<const Halfedge*>(this->p_next)); }

  /*! Get the next halfedge along the chain (const version). */
  Halfedge* next() { return (reinterpret_cast<Halfedge*>(this->p_next)); }

  /*! Set the next halfedge along the chain. */
  void set_next(Halfedge* he)
  {
    this->p_next = he;
    he->p_prev = this;
  }

  /*! Get the target vertex (const version). */
  const Vertex* vertex() const
  { return (reinterpret_cast<const Vertex*>(_clean_pointer(this->p_v))); }

  /*! Get the target vertex (non-const version). */
  Vertex* vertex()
  { return (reinterpret_cast<Vertex*>(_clean_pointer(this->p_v))); }

  /*! Set the target vertex. */
  void set_vertex(Vertex* v)
  {
    // Set the vertex pointer, preserving the content of the LSB.
    if (_is_lsb_set(this->p_v)) this->p_v = _set_lsb(v);
    else this->p_v = v;
  }

  /*! Check whether the halfedge lies on the boundary of an outer CCB. */
  bool is_on_outer_ccb() const { return (!_is_lsb_set(this->p_comp)); }

  /*! Get an incident outer CCB (const version).
   * \pre The edge does not lie on an inner CCB.
   */
  const Outer_ccb* outer_ccb() const
  {
    CGAL_precondition(! is_on_inner_ccb());
    return (reinterpret_cast<const Outer_ccb*>(this->p_comp));
  }

  /*! Get an incident outer CCB (non-const version).
   * \pre The edge does not lie on an inner CCB.
   */
  Outer_ccb* outer_ccb()
  {
    CGAL_precondition(! is_on_inner_ccb());
    return (reinterpret_cast<Outer_ccb*>(this->p_comp));
  }

  /*! Set the incident outer CCB. */
  void set_outer_ccb(Outer_ccb *oc)
  {
    // Set the component pointer and reset its LSB.
    this->p_comp = oc;
  }

  /*! Check whether the halfedge lies on the boundary of an inner CCB. */
  bool is_on_inner_ccb() const { return (_is_lsb_set(this->p_comp)); }

  /*! Get an incident inner CCB (const version).
   * \pre The edge lies on an inner CCB.
   */
  const Inner_ccb* inner_ccb() const
  {
    CGAL_precondition(is_on_inner_ccb());
    return (reinterpret_cast<const Inner_ccb*>(_clean_pointer(this->p_comp)));
  }

  /*! Get an incident inner CCB (non-const version).
   * \pre The edge lies on an inner CCB.
   */
  Inner_ccb* inner_ccb()
  {
    CGAL_precondition(is_on_inner_ccb());
    return (reinterpret_cast<Inner_ccb*>(_clean_pointer(this->p_comp)));
  }

  /*! Set the incident inner CCB. */
  void set_inner_ccb(Inner_ccb *ic)
  {
    // Set the component pointer and set its LSB.
    this->p_comp = _set_lsb(ic);
  }
};

/*! \class
 * The default arrangement DCEL face class.
 */
template <class V, class H, class F>
class Arr_face : public F,
                 public In_place_list_base<Arr_face<V,H,F> >
{
public:
  typedef F                            Base;
  typedef Arr_vertex<V,H,F>            Vertex;
  typedef Arr_halfedge<V,H,F>          Halfedge;
  typedef Arr_face<V,H,F>              Face;
  typedef Arr_outer_ccb<V,H,F>         Outer_ccb;
  typedef Arr_inner_ccb<V,H,F>         Inner_ccb;
  typedef Arr_isolated_vertex<V,H,F>   Isolated_vertex;

  typedef Inner_ccb                    Hole;

private:
  typedef Cast_function_object<void*, Halfedge*> _Ccb_to_halfedge_cast;
  // typedef Cast_function_object<const void*, const Halfedge*>
  //                                             _Const_ccb_to_halfedge_cast;
  typedef _Ccb_to_halfedge_cast                  _Const_ccb_to_halfedge_cast;

public:

  /*! Default constructor. */
  Arr_face()
  {}

  // Definition of the outer CCB iterators:
  typedef Iterator_project<typename F::Outer_ccb_iterator,
                           _Ccb_to_halfedge_cast>   Outer_ccb_iterator;

  typedef Iterator_project<typename F::Outer_ccb_const_iterator,
                           _Const_ccb_to_halfedge_cast>
                                                    Outer_ccb_const_iterator;

  /*! Get the number of outer CCBs the face has. */
  size_t number_of_outer_ccbs() const { return (this->outer_ccbs.size()); }

  /*! Get an iterator for the first outer CCB of the face. */
  Outer_ccb_iterator outer_ccbs_begin() { return (this->outer_ccbs.begin()); }

  /*! Get a past-the-end iterator for the outer CCBs inside the face. */
  Outer_ccb_iterator outer_ccbs_end() { return (this->outer_ccbs.end()); }

  /*! Get an const iterator for the first outer CCB inside the face. */
  Outer_ccb_const_iterator outer_ccbs_begin() const
  { return (this->outer_ccbs.begin()); }

  /*! Get a const past-the-end iterator for the outer CCBs inside the face. */
  Outer_ccb_const_iterator outer_ccbs_end() const
  { return (this->outer_ccbs.end()); }

  /*! Add an outer CCB to the face. */
  void add_outer_ccb(Outer_ccb *oc, Halfedge *h)
  { oc->set_iterator(this->outer_ccbs.insert(this->outer_ccbs.end(), h)); }

  /*! Erase an outer CCB of the face. */
  void erase_outer_ccb(Outer_ccb *oc)
  { this->outer_ccbs.erase(oc->iterator().current_iterator()); }

  // Definition of the inner CCB iterators:
  typedef Iterator_project<typename F::Inner_ccb_iterator,
                           _Ccb_to_halfedge_cast>   Inner_ccb_iterator;

  typedef Iterator_project<typename F::Inner_ccb_const_iterator,
                           _Const_ccb_to_halfedge_cast>
                                                    Inner_ccb_const_iterator;

  typedef Inner_ccb_iterator                        Hole_iterator;
  typedef Inner_ccb_const_iterator                  Hole_const_iterator;

  /*! Get the number of inner CCBs the face has. */
  size_t number_of_inner_ccbs() const { return (this->inner_ccbs.size()); }

  /*! Get an iterator for the first inner CCB of the face. */
  Inner_ccb_iterator inner_ccbs_begin() { return (this->inner_ccbs.begin()); }

  /*! Get a past-the-end iterator for the inner CCBs inside the face. */
  Inner_ccb_iterator inner_ccbs_end() { return (this->inner_ccbs.end()); }

  /*! Get an const iterator for the first inner CCB inside the face. */
  Inner_ccb_const_iterator inner_ccbs_begin() const
  { return (this->inner_ccbs.begin()); }

  /*! Get a const past-the-end iterator for the inner CCBs inside the face. */
  Inner_ccb_const_iterator inner_ccbs_end() const
  { return (this->inner_ccbs.end()); }

  /*! Add an inner CCB to the face. */
  void add_inner_ccb(Inner_ccb* ic, Halfedge* h)
  { ic->set_iterator(this->inner_ccbs.insert(this->inner_ccbs.end(), h)); }

  /*! Erase an inner CCB of the face. */
  void erase_inner_ccb(Inner_ccb* ic)
  { this->inner_ccbs.erase(ic->iterator().current_iterator()); }

  /*! Move all inner CCBs (holes) from the face to another. */
  Inner_ccb_iterator splice_inner_ccbs(Arr_face& other)
  {
    const bool was_empty = this->inner_ccbs.empty();
    typename Base::Inner_ccbs_container::iterator previous =
      this->inner_ccbs.end();
    if (!was_empty) --previous;
    this->inner_ccbs.splice(this->inner_ccbs.end(), other.inner_ccbs);
    if (was_empty) previous = this->inner_ccbs.begin();
    else ++previous;
    for (typename Base::Inner_ccbs_container::iterator it = previous;
         it != this->inner_ccbs.end(); ++it)
    {
      Inner_ccb* ccb = static_cast<Halfedge*>(*it)->inner_ccb();
      ccb->set_iterator(it);
      ccb->set_face(this);
    }
    return previous;
  }

  // Backward compatibility:
  size_t number_of_holes() const { return number_of_inner_ccbs(); }
  Hole_iterator holes_begin() { return inner_ccbs_begin(); }
  Hole_iterator holes_end() { return inner_ccbs_end(); }
  Hole_const_iterator holes_begin() const { return inner_ccbs_begin(); }
  Hole_const_iterator holes_end() const { return inner_ccbs_end(); }

  // Definition of the isloated vertices iterators:
  typedef I_Dereference_iterator<
    typename F::Isolated_vertex_iterator,
    Vertex,
    typename F::Isolated_vertex_iterator::difference_type,
    typename F::Isolated_vertex_iterator::iterator_category>
                                              Isolated_vertex_iterator;

  typedef I_Dereference_const_iterator<
    typename F::Isolated_vertex_const_iterator,
    typename F::Isolated_vertex_iterator,
    Vertex,
    typename F::Isolated_vertex_iterator::difference_type,
    typename F::Isolated_vertex_iterator::iterator_category>
                                              Isolated_vertex_const_iterator;

  /*! Get the number of isloated vertices inside the face. */
  size_t number_of_isolated_vertices() const
  { return (this->iso_verts.size()); }

  /*! Get an iterator for the first isloated vertex inside the face. */
  Isolated_vertex_iterator isolated_vertices_begin()
  { return (this->iso_verts.begin()); }

  /*! Get a past-the-end iterator for the isloated vertices inside the face. */
  Isolated_vertex_iterator isolated_vertices_end()
  { return (this->iso_verts.end()); }

  /*! Get an const iterator for the first isloated vertex inside the face. */
  Isolated_vertex_const_iterator isolated_vertices_begin() const
  { return (this->iso_verts.begin()); }

  /*! Get a const past-the-end iterator for the isloated vertices inside the
   * face. */
  Isolated_vertex_const_iterator isolated_vertices_end() const
  { return (this->iso_verts.end()); }

  /*! Add an isloated vertex inside the face. */
  void add_isolated_vertex(Isolated_vertex *iv, Vertex* v)
  { iv->set_iterator(this->iso_verts.insert(this->iso_verts.end(), v)); }

  /*! Erase an isloated vertex from the face. */
  void erase_isolated_vertex(Isolated_vertex *iv)
  { this->iso_verts.erase(iv->iterator().current_iterator()); }

  /*! Move all isolated vertices from the face to another. */
  Isolated_vertex_iterator splice_isolated_vertices(Arr_face& other)
  {
    const bool was_empty = this->iso_verts.empty();
    typename Base::Isolated_vertices_container::iterator previous =
      this->iso_verts.end();
    if (!was_empty) --previous;
    this->iso_verts.splice(this->iso_verts.end(), other.iso_verts);
    if (was_empty) previous = this->iso_verts.begin();
    else ++previous;
    for (typename Base::Isolated_vertices_container::iterator it = previous;
         it != this->iso_verts.end(); ++it)
    {
      Isolated_vertex* iv = static_cast<Vertex*>(*it)->isolated_vertex();
      iv->set_iterator(it);
      iv->set_face(this);
    }
    return previous;
  }
};

/*! \class
 * Representation of an outer CCB.
 */
template <class V, class H, class F>
class Arr_outer_ccb : public In_place_list_base<Arr_outer_ccb<V,H,F> > {
public:
  typedef Arr_outer_ccb<V,H,F>               Self;
  typedef Arr_halfedge<V,H,F>                Halfedge;
  typedef Arr_face<V,H,F>                    Face;
  typedef typename Face::Outer_ccb_iterator  Outer_ccb_iterator;

private:
  Face* p_f;                  // The face the contains the CCB in its interior.
  Outer_ccb_iterator  iter;   // The outer CCB identifier.
  bool iter_is_not_singular;

public:
  /*! Default constructor. */
  Arr_outer_ccb() : p_f(NULL), iter_is_not_singular(false) {}

  /*! Copy constructor. */
  Arr_outer_ccb(const Arr_outer_ccb& other) :
    p_f(other.p_f), iter_is_not_singular(other.iter_is_not_singular)
  { if (other.iter_is_not_singular) iter = other.iter; }

  /*! Get a halfedge along the component (const version). */
  const Halfedge* halfedge() const { return (*iter); }

  /*! Get a halfedge along the component (non-const version). */
  Halfedge* halfedge() { return (*iter); }

  /*! Set a representative halfedge for the component. */
  void set_halfedge(Halfedge* he) { *iter = he; }

  /*! Get the incident face (const version). */
  const Face* face() const { return (p_f); }

  /*! Get the incident face (non-const version). */
  Face* face() { return (p_f); }

  /*! Set the incident face. */
  void set_face(Face* f) { p_f = f; }

  /*! Get the iterator (const version). */
  Outer_ccb_iterator iterator() const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Get the iterator (non-const version). */
  Outer_ccb_iterator iterator()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Set the outer CCB iterator. */
  void set_iterator(Outer_ccb_iterator it)
  {
    iter = it;
    iter_is_not_singular = true;
  }
};

/*! \class
 * Representation of an inner CCB.
 */
template <class V, class H, class F>
class Arr_inner_ccb : public In_place_list_base<Arr_inner_ccb<V,H,F> >
{
public:
  typedef Arr_inner_ccb<V,H,F>               Self;
  typedef Arr_halfedge<V,H,F>                Halfedge;
  typedef Arr_face<V,H,F>                    Face;
  typedef typename Face::Inner_ccb_iterator  Inner_ccb_iterator;

private:
  Face* p_f;                  // The face the contains the CCB in its interior.
  Inner_ccb_iterator  iter;   // The inner CCB identifier.
  bool iter_is_not_singular;

public:
  /*! Default constructor. */
  Arr_inner_ccb() : p_f(NULL), iter_is_not_singular(false) {}

  /*! Copy constructor. */
  Arr_inner_ccb(const Arr_inner_ccb& other) :
    p_f(other.p_f), iter_is_not_singular(other.iter_is_not_singular)
  { if (other.iter_is_not_singular) iter = other.iter; }

  /*! Get a halfedge along the component (const version). */
  const Halfedge* halfedge() const { return (*iter); }

  /*! Get a halfedge along the component (non-const version). */
  Halfedge* halfedge() { return (*iter); }

  /*! Set a representative halfedge for the component. */
  void set_halfedge(Halfedge *he) { *iter = he; }

  /*! Get the incident face (const version). */
  const Face* face() const { return (p_f); }

  /*! Get the incident face (non-const version). */
  Face* face() { return (p_f); }

  /*! Set the incident face. */
  void set_face(Face* f) { p_f = f; }

  /*! Get the iterator (const version). */
  Inner_ccb_iterator iterator() const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Get the iterator (non-const version). */
  Inner_ccb_iterator iterator()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iter);
  }

  /*! Set the inner CCB iterator. */
  void set_iterator(Inner_ccb_iterator it)
  {
    iter = it;
    iter_is_not_singular = true;
  }
};

/*! \class
 * Representation of an isolated vertex.
 */
template <class V, class H, class F>
class Arr_isolated_vertex :
public In_place_list_base<Arr_isolated_vertex<V,H,F> > {
public:
  typedef Arr_isolated_vertex<V,H,F>                Self;
  typedef Arr_face<V,H,F>                           Face;
  typedef typename Face::Isolated_vertex_iterator   Isolated_vertex_iterator;

private:
  Face* p_f;                             // The containing face.
  Isolated_vertex_iterator   iv_it;     // The isolated vertex identifier.
  bool iter_is_not_singular;

public:
  /*! Default constructor. */
  Arr_isolated_vertex() : p_f(NULL), iter_is_not_singular(false) {}

  /*! Copy constructor. */
  Arr_isolated_vertex(const Arr_isolated_vertex& other) :
    p_f(other.p_f), iter_is_not_singular(other.iter_is_not_singular)
  { if (other.iter_is_not_singular) iv_it = other.iv_it; }

  /*! Get the containing face (const version). */
  const Face* face() const { return (p_f); }

  /*! Get the containing face (non-const version). */
  Face* face() { return (p_f); }

  /*! Set the incident face, the one that contains the isolated vertex. */
  void set_face(Face* f) { p_f = f; }

  /*! Get the isolated vertex iterator (const version). */
  Isolated_vertex_iterator iterator() const
  {
    CGAL_assertion(iter_is_not_singular);
    return (iv_it);
  }

  /*! Get the isolated vertex iterator (non-const version). */
  Isolated_vertex_iterator iterator()
  {
    CGAL_assertion(iter_is_not_singular);
    return (iv_it);
  }

  /*! Set the isolated vertex iterator. */
  void set_iterator(Isolated_vertex_iterator iv)
  {
    iv_it = iv;
    iter_is_not_singular = true;
  }
};

/*! \class
 * The arrangement DCEL class.
 */
template <class V, class H, class F,
          class Allocator = CGAL_ALLOCATOR(int) >
class Arr_dcel_base {
public:
  // Define the vertex, halfedge and face types.
  typedef Arr_dcel_base<V,H,F>        Self;
  typedef Arr_vertex<V,H,F>           Vertex;
  typedef Arr_halfedge<V,H,F>         Halfedge;
  typedef Arr_face<V,H,F>             Face;
  typedef Arr_outer_ccb<V,H,F>        Outer_ccb;
  typedef Arr_inner_ccb<V,H,F>        Inner_ccb;
  typedef Arr_isolated_vertex<V,H,F>  Isolated_vertex;

  typedef Inner_ccb                   Hole;

protected:
  // The vetices, halfedges and faces are stored in three in-place lists.
  typedef In_place_list<Vertex, false>           Vertex_list;
  typedef In_place_list<Halfedge, false>         Halfedge_list;
  typedef In_place_list<Face, false>             Face_list;
  typedef In_place_list<Outer_ccb, false>        Outer_ccb_list;
  typedef In_place_list<Inner_ccb, false>        Inner_ccb_list;
  typedef In_place_list<Isolated_vertex, false>  Iso_vert_list;

#ifdef CGAL_CXX11
    typedef std::allocator_traits<Allocator> Allocator_traits;
    typedef typename Allocator_traits::template rebind_alloc<Vertex>          Vertex_allocator;
    typedef typename Allocator_traits::template rebind_alloc<Halfedge>        Halfedge_allocator;
    typedef typename Allocator_traits::template rebind_alloc<Face>            Face_allocator;
    typedef typename Allocator_traits::template rebind_alloc<Outer_ccb>       Outer_ccb_allocator;
    typedef typename Allocator_traits::template rebind_alloc<Inner_ccb>       Inner_ccb_allocator;
    typedef typename Allocator_traits::template rebind_alloc<Isolated_vertex> Iso_vert_allocator;
#else // not CGAL_CXX11
  // Vertex allocator.
  typedef typename Allocator::template rebind<Vertex>    Vertex_alloc_rebind;
  typedef typename Vertex_alloc_rebind::other            Vertex_allocator;

  // Halfedge allocator.
  typedef typename Allocator::template rebind<Halfedge>  Halfedge_alloc_rebind;
  typedef typename Halfedge_alloc_rebind::other          Halfedge_allocator;

  // Face allocator.
  typedef typename Allocator::template rebind<Face>      Face_alloc_rebind;
  typedef typename Face_alloc_rebind::other              Face_allocator;

  // Outer CCB allocator.
  typedef typename Allocator::template rebind<Outer_ccb> Out_ccb_alloc_rebind;
  typedef typename Out_ccb_alloc_rebind::other           Outer_ccb_allocator;

  // Inner CCB allocator.
  typedef typename Allocator::template rebind<Inner_ccb> In_ccb_alloc_rebind;
  typedef typename In_ccb_alloc_rebind::other            Inner_ccb_allocator;

  // Isolated vertex allocator.
  typedef typename Allocator::template rebind<Isolated_vertex>
                                                         Iso_vert_alloc_rebind;
  typedef typename Iso_vert_alloc_rebind::other          Iso_vert_allocator;
#endif // not CGAL_CXX11

public:
  typedef typename Halfedge_list::size_type              Size;
  typedef typename Halfedge_list::size_type              size_type;
  typedef typename Halfedge_list::difference_type        difference_type;
  typedef typename Halfedge_list::difference_type        Difference;
  typedef std::bidirectional_iterator_tag                iterator_category;

protected:

  Vertex_list         vertices;             // The vertices container.
  Halfedge_list       halfedges;            // The halfedges container.
  Face_list           faces;                // The faces container.
  Outer_ccb_list      out_ccbs;             // The outer CCBs.
  Inner_ccb_list      in_ccbs;              // The inner CCBs.
  Iso_vert_list       iso_verts;            // The isolated vertices.

  Vertex_allocator    vertex_alloc;         // An allocator for vertices.
  Halfedge_allocator  halfedge_alloc;       // An allocator for halfedges.
  Face_allocator      face_alloc;           // An allocator for faces.
  Outer_ccb_allocator out_ccb_alloc;        // An allocator for outer CCBs.
  Inner_ccb_allocator in_ccb_alloc;         // An allocator for inner CCBs.
  Iso_vert_allocator  iso_vert_alloc;       // Allocator for isolated vertices.

public:
  // Definitions of iterators.
  typedef typename Vertex_list::iterator              Vertex_iterator;
  typedef typename Halfedge_list::iterator            Halfedge_iterator;
  typedef typename Face_list::iterator                Face_iterator;
  typedef CGAL::N_step_adaptor_derived<Halfedge_iterator, 2>
                                                      Edge_iterator;

  // Definitions of const iterators.
  typedef typename Vertex_list::const_iterator        Vertex_const_iterator;
  typedef typename Halfedge_list::const_iterator      Halfedge_const_iterator;
  typedef typename Face_list::const_iterator          Face_const_iterator;
  typedef CGAL::N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                      Edge_const_iterator;

private:
  // Copy constructor - not supported.
  Arr_dcel_base(const Self&);

  // Assignment operator - not supported.
  Self& operator=(const Self&);

public:
  /// \name Construction and destruction.
  //@{
  /*! Default constructor. */
  Arr_dcel_base() {}

  /*! Destructor. */
  ~Arr_dcel_base() { delete_all(); }
  //@}

  /// \name The DCEL size.
  //@{
  /*! Get the number of DCEL vertices. */
  Size size_of_vertices() const { return (vertices.size()); }

  /*! Get the number of DCEL halfedges (twice the number of edges). */
  Size size_of_halfedges() const { return (halfedges.size()); }

  /*! Get the number of DCEL faces. */
  Size size_of_faces() const { return (faces.size()); }

  /*! Get the number of outer CCBs. */
  Size size_of_outer_ccbs() const { return (out_ccbs.size()); }

  /*! Get the number of inner CCBs. */
  Size size_of_inner_ccbs() const { return (in_ccbs.size()); }

  /*! Get the number of isolated vertices. */
  Size size_of_isolated_vertices() const { return (iso_verts.size()); }
  //@}

  /// \name Obtaining iterators.
  //@{
  Vertex_iterator   vertices_begin()  { return vertices.begin(); }
  Vertex_iterator   vertices_end()    { return vertices.end(); }
  Iterator_range<Prevent_deref<Vertex_iterator> >
  vertex_handles()
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end());
  }
  Halfedge_iterator halfedges_begin() { return halfedges.begin();}
  Halfedge_iterator halfedges_end()   { return halfedges.end(); }
  Iterator_range<Prevent_deref<Halfedge_iterator> >
  halfedge_handles()
  {
    return make_prevent_deref_range(halfedges_begin(), halfedges_end());
  }
  Face_iterator     faces_begin()     { return faces.begin(); }
  Face_iterator     faces_end()       { return faces.end(); }
  Iterator_range<Prevent_deref<Face_iterator> >
  face_handles()
  {
    return make_prevent_deref_range(faces_begin(), faces_end());
  }
  Edge_iterator     edges_begin()     { return halfedges.begin(); }
  Edge_iterator     edges_end()       { return halfedges.end(); }
  Iterator_range<Prevent_deref<Edge_iterator> >
  edge_handles()
  {
    return make_prevent_deref_range(edges_begin(), edges_end());
  }
  //@}

  /// \name Obtaining constant iterators.
  //@{
  Vertex_const_iterator   vertices_begin() const { return vertices.begin(); }
  Vertex_const_iterator   vertices_end() const { return vertices.end(); }
  Iterator_range<Prevent_deref<Vertex_const_iterator> >
  vertex_handles() const
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end());
  }
  Halfedge_const_iterator halfedges_begin() const { return halfedges.begin(); }
  Halfedge_const_iterator halfedges_end() const { return halfedges.end(); }
  Iterator_range<Prevent_deref<Halfedge_const_iterator> >
  halfedge_handles() const
  {
    return make_prevent_deref_range(halfedges_begin(), halfedges_end());
  }
  Face_const_iterator     faces_begin() const { return faces.begin(); }
  Face_const_iterator     faces_end() const { return faces.end(); }
  Iterator_range<Prevent_deref<Face_const_iterator> >
  face_handles() const
  {
    return make_prevent_deref_range(faces_begin(), faces_end());
  }
  Edge_const_iterator     edges_begin() const { return halfedges.begin(); }
  Edge_const_iterator     edges_end() const { return halfedges.end(); }
  Iterator_range<Prevent_deref<Edge_const_iterator> >
  edge_handles() const
  {
    return make_prevent_deref_range(edges_begin(), edges_end());
  }
  //@}

  // \name Creation of new DCEL features.
  //@{
  /*! Create a new vertex. */
  Vertex* new_vertex()
  {
    Vertex* v = vertex_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Vertex_allocator>::construct(vertex_alloc,v);
#else
    vertex_alloc.construct(v, Vertex());
#endif
    vertices.push_back(*v);
    return v;
  }

  /*! Create a new pair of opposite halfedges. */
  Halfedge* new_edge()
  {
    // Create two new halfedges.
    Halfedge* h1 = _new_halfedge();
    Halfedge* h2 = _new_halfedge();

    // Pair them together.
    h1->set_opposite(h2);
    h2->set_opposite(h1);

    return (h1);
  }

  /*! Create a new face. */
  Face* new_face()
  {
    Face* f = face_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Face_allocator>::construct(face_alloc, f);
#else
    face_alloc.construct(f, Face());
#endif
    faces.push_back (*f);
    return(f);
  }

  /*! Create a new outer CCB. */
  Outer_ccb* new_outer_ccb()
  {
    Outer_ccb* oc = out_ccb_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Outer_ccb_allocator>::construct(out_ccb_alloc, oc);
#else
    out_ccb_alloc.construct(oc, Outer_ccb());
#endif
    out_ccbs.push_back(*oc);
    return (oc);
  }

  /*! Create a new inner CCB. */
  Inner_ccb* new_inner_ccb()
  {
    Inner_ccb* ic = in_ccb_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Inner_ccb_allocator>::construct(in_ccb_alloc, ic);
#else
    in_ccb_alloc.construct(ic, Inner_ccb());
#endif
    in_ccbs.push_back(*ic);
    return (ic);
  }

  /*! Create a new isolated vertex. */
  Isolated_vertex* new_isolated_vertex()
  {
    Isolated_vertex* iv = iso_vert_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Iso_vert_allocator>::construct(iso_vert_alloc, iv);
#else
    iso_vert_alloc.construct(iv, Isolated_vertex());
#endif
    iso_verts.push_back(*iv);
    return (iv);
  }
  //@}

  /// \name Deletion of DCEL features.
  //@{
  /*! Delete an existing vertex. */
  void delete_vertex(Vertex* v)
  {
    vertices.erase(v);
#ifdef CGAL_CXX11
    std::allocator_traits<Vertex_allocator>::destroy(vertex_alloc, v);
#else
    vertex_alloc.destroy(v);
#endif
    vertex_alloc.deallocate(v,1);
  }

  /*! Delete an existing pair of opposite halfedges. */
  void delete_edge(Halfedge *h)
  {
    Halfedge* h_opp = h->opposite();
    _delete_halfedge(h);
    _delete_halfedge(h_opp);
  }

  /*! Delete an existing face. */
  void delete_face(Face* f)
  {
    faces.erase(f);
#ifdef CGAL_CXX11
    std::allocator_traits<Face_allocator>::destroy(face_alloc, f);
#else
    face_alloc.destroy(f);
#endif
    face_alloc.deallocate(f, 1);
  }

  /*! Delete an existing outer CCB. */
  void delete_outer_ccb(Outer_ccb* oc)
  {
    out_ccbs.erase(oc);
#ifdef CGAL_CXX11
    std::allocator_traits<Outer_ccb_allocator>::destroy(out_ccb_alloc, oc);
#else
    out_ccb_alloc.destroy(oc);
#endif
    out_ccb_alloc.deallocate(oc, 1);
  }

  /*! Delete an existing inner CCB. */
  void delete_inner_ccb(Inner_ccb* ic)
  {
    in_ccbs.erase(ic);
#ifdef CGAL_CXX11
    std::allocator_traits<Inner_ccb_allocator>::destroy(in_ccb_alloc, ic);
#else
    in_ccb_alloc.destroy(ic);
#endif
    in_ccb_alloc.deallocate(ic, 1);
  }

  /*! Delete an existing isolated vertex. */
  void delete_isolated_vertex(Isolated_vertex* iv)
  {
    iso_verts.erase(iv);
#ifdef CGAL_CXX11
    std::allocator_traits<Iso_vert_allocator>::destroy(iso_vert_alloc, iv);
#else
    iso_vert_alloc.destroy(iv);
#endif
    iso_vert_alloc.deallocate(iv, 1);
  }

  /*! Delete all DCEL features. */
  void delete_all()
  {
    // Free all vertices.
    Vertex_iterator vit = vertices.begin(), v_curr;
    while (vit != vertices.end()) {
      v_curr = vit;
      ++vit;
      delete_vertex(&(*v_curr));
    }

    // Free all halfedges.
    Halfedge_iterator  hit = halfedges.begin(), h_curr;
    while (hit != halfedges.end()) {
      h_curr = hit;
      ++hit;
      _delete_halfedge(&(*h_curr));
    }

    // Free all faces.
    Face_iterator fit = faces.begin(), f_curr;
    while (fit != faces.end()) {
      f_curr = fit;
      ++fit;
      delete_face(&(*f_curr));
    }

    // Free all outer CCBs.
    typename Outer_ccb_list::iterator ocit = out_ccbs.begin(), oc_curr;
    while (ocit != out_ccbs.end()) {
      oc_curr = ocit;
      ++ocit;
      delete_outer_ccb(&(*oc_curr));
    }

    // Free all inner CCBs.
    typename Inner_ccb_list::iterator icit = in_ccbs.begin(), ic_curr;
    while (icit != in_ccbs.end()) {
      ic_curr = icit;
      ++icit;
      delete_inner_ccb(&(*ic_curr));
    }

    // Free all isolated vertices.
    typename Iso_vert_list::iterator ivit = iso_verts.begin(), iv_curr;
    while (ivit != iso_verts.end()) {
      iv_curr = ivit;
      ++ivit;
      delete_isolated_vertex(&(*iv_curr));
    }
  }
  //@}

  /*! Assign our DCEL the contents of another DCEL.
   */
  void assign(const Self& dcel)
  {
    // Clear the current contents of the DCEL.
    delete_all();

    // Create duplicated of the DCEL features and map the features of the
    // given DCEL to their corresponding duplicates.
    typedef std::map<const Vertex*, Vertex*>                    Vertex_map;
    typedef std::map<const Halfedge*, Halfedge*>                Halfedge_map;
    typedef std::map<const Face*, Face*>                        Face_map;
    typedef std::map<const Outer_ccb*, Outer_ccb*>              Outer_ccb_map;
    typedef std::map<const Inner_ccb*, Inner_ccb*>              Inner_ccb_map;
    typedef std::map<const Isolated_vertex*, Isolated_vertex*>  Iso_vert_map;

    Vertex_map v_map;
    Vertex_const_iterator vit;
    Vertex* dup_v;

    for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit) {
      dup_v = new_vertex();
      dup_v->assign(*vit);
      v_map.insert(typename Vertex_map::value_type(&(*vit), dup_v));
    }

    Halfedge_map he_map;
    Halfedge_const_iterator hit;
    Halfedge* dup_h;

    for (hit = dcel.halfedges_begin(); hit != dcel.halfedges_end(); ++hit) {
      dup_h = _new_halfedge();
      dup_h->assign(*hit);
      he_map.insert(typename Halfedge_map::value_type(&(*hit), dup_h));
    }

    Face_map f_map;
    Face_const_iterator fit;
    Face* dup_f;

    for (fit = dcel.faces_begin(); fit != dcel.faces_end(); ++fit) {
      dup_f = new_face();
      dup_f->assign(*fit);
      f_map.insert(typename Face_map::value_type(&(*fit), dup_f));
    }

    Outer_ccb_map oc_map;
    typename Outer_ccb_list::const_iterator ocit;
    Outer_ccb* dup_oc;

    for (ocit = dcel.out_ccbs.begin(); ocit != dcel.out_ccbs.end(); ++ocit) {
      dup_oc = new_outer_ccb();
      oc_map.insert(typename Outer_ccb_map::value_type(&(*ocit), dup_oc));
    }

    Inner_ccb_map ic_map;
    typename Inner_ccb_list::const_iterator icit;
    Inner_ccb* dup_ic;

    for (icit = dcel.in_ccbs.begin(); icit != dcel.in_ccbs.end(); ++icit) {
      dup_ic = new_inner_ccb();
      ic_map.insert(typename Inner_ccb_map::value_type(&(*icit), dup_ic));
    }

    Iso_vert_map iv_map;
    typename Iso_vert_list::const_iterator ivit;
    Isolated_vertex* dup_iv;

    for (ivit = dcel.iso_verts.begin(); ivit != dcel.iso_verts.end(); ++ivit) {
      dup_iv = new_isolated_vertex();
      iv_map.insert(typename Iso_vert_map::value_type(&(*ivit), dup_iv));
    }

    // Update the vertex records.
    const Vertex* v;
    const Halfedge* h;
    const Face* f;
    const Outer_ccb* oc;
    const Inner_ccb* ic;
    const Isolated_vertex* iv;

    for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit) {
      v = &(*vit);
      dup_v = (v_map.find(v))->second;

      if (v->is_isolated()) {
        // Isolated vertex - set its information.
        iv = v->isolated_vertex();
        dup_iv = (iv_map.find(iv))->second;

        dup_v->set_isolated_vertex(dup_iv);
      }
      else {
        // Regular vertex - set its incident halfedge.
        h = v->halfedge();
        dup_h = (he_map.find(h))->second;

        dup_v->set_halfedge(dup_h);
      }
    }

    // Update the halfedge records.
    const Halfedge* opp;
    const Halfedge* prev;
    const Halfedge* next;
    Halfedge* dup_opp;
    Halfedge* dup_prev;
    Halfedge* dup_next;

    for (hit = dcel.halfedges_begin(); hit != dcel.halfedges_end(); ++hit) {
      h = &(*hit);
      v = h->vertex();
      opp = h->opposite();
      prev = h->prev();
      next = h->next();

      dup_h = (he_map.find(h))->second;
      dup_v = (v_map.find(v))->second;
      dup_opp = (he_map.find(opp))->second;
      dup_prev = (he_map.find(prev))->second;
      dup_next = (he_map.find(next))->second;

      dup_h->set_vertex(dup_v);
      dup_h->set_opposite(dup_opp);
      dup_h->set_prev(dup_prev);
      dup_h->set_next(dup_next);
      dup_h->set_direction(h->direction());

      if (h->is_on_inner_ccb()) {
        // The halfedge lies on an inner CCB - set its inner CCB record.
        ic = h->inner_ccb();
        dup_ic = (ic_map.find(ic))->second;
        dup_h->set_inner_ccb(dup_ic);
      }
      else {
        // The halfedge lies on an outer CCB - set its outer CCB record.
        oc = h->outer_ccb();
        dup_oc = (oc_map.find(oc))->second;
        dup_h->set_outer_ccb(dup_oc);
      }
    }

    // Update the face records, along with the CCB records and isolated vertex
    // records.
    typename Face::Outer_ccb_const_iterator out_ccb_it;
    typename Face::Inner_ccb_const_iterator in_ccb_it;
    typename Face::Isolated_vertex_const_iterator iso_vert_it;
    const Halfedge* hccb;
    const Vertex* iso_vert;
    Halfedge* dup_hccb;
    Vertex* dup_iso_vert;

    for (fit = dcel.faces_begin(); fit != dcel.faces_end(); ++fit) {
      f = &(*fit);
      dup_f = (f_map.find(f))->second;
      dup_f->set_unbounded(f->is_unbounded());
      dup_f->set_fictitious(f->is_fictitious());

      // Assign the outer CCBs of the face.
      for (out_ccb_it = f->outer_ccbs_begin();
           out_ccb_it != f->outer_ccbs_end(); ++out_ccb_it)
      {
        hccb = *out_ccb_it;

        dup_hccb = (he_map.find(hccb))->second;
        dup_oc = dup_hccb->outer_ccb();

        dup_oc->set_face(dup_f);
        dup_f->add_outer_ccb(dup_oc, dup_hccb);
      }

      // Assign the inner CCBs of the face.
      for (in_ccb_it = f->inner_ccbs_begin();
           in_ccb_it != f->inner_ccbs_end(); ++in_ccb_it)
      {
        hccb = *in_ccb_it;

        dup_hccb = (he_map.find(hccb))->second;
        dup_ic = dup_hccb->inner_ccb();

        dup_ic->set_face(dup_f);
        dup_f->add_inner_ccb(dup_ic, dup_hccb);
      }

      // Assign the isolated vertices.
      for (iso_vert_it = f->isolated_vertices_begin();
           iso_vert_it != f->isolated_vertices_end(); ++iso_vert_it)
      {
        iso_vert = &(*iso_vert_it);

        dup_iso_vert = (v_map.find(iso_vert))->second;
        dup_iv = dup_iso_vert->isolated_vertex();

        dup_iv->set_face(dup_f);
        dup_f->add_isolated_vertex(dup_iv, dup_iso_vert);
      }
    }
  }

protected:
  /*! Create a new halfedge. */
  Halfedge* _new_halfedge()
  {
    Halfedge* h = halfedge_alloc.allocate(1);
#ifdef CGAL_CXX11
    std::allocator_traits<Halfedge_allocator>::construct(halfedge_alloc, h);
#else
    halfedge_alloc.construct(h, Halfedge());
#endif
    halfedges.push_back(*h);
    return (h);
  }

  /*! Delete an existing halfedge. */
  void _delete_halfedge(Halfedge* h)
  {
    halfedges.erase(h);
#ifdef CGAL_CXX11
    std::allocator_traits<Halfedge_allocator>::destroy(halfedge_alloc,h);
#else
    halfedge_alloc.destroy(h);
#endif
    halfedge_alloc.deallocate(h, 1);
  }
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
