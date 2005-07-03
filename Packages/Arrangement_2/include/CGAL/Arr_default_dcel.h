// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 (based on old version by: Iddo Hanniel and Oren Nechushtan)
#ifndef CGAL_ARR_DEFAULT_DCEL_H
#define CGAL_ARR_DEFAULT_DCEL_H

/*! \file
 * The definition of the Arr_default_dcel<Traits> class.
 */

#include <CGAL/basic.h>
#include <list>
#include <map>
#include <CGAL/N_step_adaptor.h>
#include <CGAL/In_place_list.h>
#include <CGAL/HalfedgeDS_iterator.h>
#include <CGAL/Arrangement_2/Arrangement_2_iterators.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * Base vertex class. 
 */
template <class Point_> class Arr_vertex_base 
{
public:

  typedef Point_       Point;

protected:

  void       *p_he;   // An incident halfedge pointing at the vertex.
  Point      *p_pt;   // The point associated with the vertex.

public:

  /*! Default constructor. */
  Arr_vertex_base() :
    p_he (NULL),
    p_pt (NULL)
  {}
  
  /*! Destructor. */
  virtual ~Arr_vertex_base() {}

  /*! Get the point (const version). */
  const Point& point() const 
  { 
    return (*p_pt);
  }

  /*! Get the point (non-const version). */
  Point& point() 
  { 
    return (*p_pt);
  }

  /*! Set the point. */
  void set_point (Point *p) 
  {
    p_pt = p;
  }

  /*! Assign from another vertex. */
  virtual void assign (const Arr_vertex_base<Point>& v)
  {
    p_pt = v.p_pt;
  }
};

/*! \class
 * Base halfedge class.
 */
template <class X_monotone_curve_> class Arr_halfedge_base 
{
public:

  typedef X_monotone_curve_  X_monotone_curve;

protected:

  void       *p_opp;   // The opposite halfedge.
  void       *p_prev;  // The previous halfedge in the component boundary.
  void       *p_next;  // The next halfedge in the component boundary.

  void       *p_v;     // The incident vertex (the target of the halfedge).
  void       *p_f;     // The incident face (to the left of the halfedge).
  
  X_monotone_curve *p_cv; // The associated x-monotone curve.

public:

  /*! Default constructor */
  Arr_halfedge_base() :
    p_opp (NULL),
    p_prev (NULL),
    p_next (NULL),
    p_v (NULL),
    p_f (NULL),
    p_cv (NULL)
  {}

  /*! Destructor. */
  virtual ~Arr_halfedge_base()
  {}

  /*! Get the x-monotone curve (const version). */
  const X_monotone_curve& curve() const 
  {
    return (*p_cv);
  }

  /*! Get the x-monotone curve (non-const version). */
  X_monotone_curve& curve () 
  {
    return (*p_cv);
  }


  /*! Set the x-monotone curve. */
  void set_curve (X_monotone_curve* c)
  { 
    p_cv = c;

    // Set the curve for the opposite halfedge as well.
    Arr_halfedge_base<X_monotone_curve>* opp =
      reinterpret_cast<Arr_halfedge_base<X_monotone_curve>* > (p_opp);

    opp->p_cv = c;
  }

  /*! Assign from another halfedge. */
  virtual void assign (const Arr_halfedge_base<X_monotone_curve>& he)
  {
    p_cv = he.p_cv;
  }
};

/*!
 * Base face class.
 */
class Arr_face_base
{
public:

  typedef std::list<void*>                  Holes_container;
  typedef Holes_container::iterator         Holes_iterator;
  typedef Holes_container::const_iterator   Holes_const_iterator;

  typedef std::list<void*>                  Isolated_vertices_container;
  typedef Isolated_vertices_container::iterator
                                            Isolated_vertices_iterator;
  typedef Isolated_vertices_container::const_iterator
                                            Isolated_vertices_const_iterator;

protected:

  void           *p_he;        // An incident halfedge along the face boundary.
  Holes_container              holes;      // The holes inside the face.
  Isolated_vertices_container  iso_verts;  // The isolated vertices inside
                                           // the face.

public:

  /*! Default constructor. */
  Arr_face_base() :
    p_he (NULL),
    holes()
  {}

  /*! Destructor. */
  virtual ~Arr_face_base()
  {}

  /*! Assign from another face (does nothing). */
  virtual void assign (const Arr_face_base& )
  {}
};

// Forward declarations:
template <class V, class H, class F> class Arr_vertex;
template <class V, class H, class F> class Arr_halfedge;
template <class V, class H, class F> class Arr_face;

/*! \class
 * The default arrangement DCEL vertex class.
 */
template <class V, class H, class F>
class Arr_vertex : public V,
                   public In_place_list_base<Arr_vertex<V,H,F> >
{
public:

  typedef V                     Base;
  typedef Arr_vertex<V,H,F>     Vertex;
  typedef Arr_halfedge<V,H,F>   Halfedge;
  typedef Arr_face<V,H,F>       Face;

  /*! Default constructor. */
  Arr_vertex() 
  {}

  /*! Get an incident halfedge (const version). */
  const Halfedge* halfedge () const
  {
    return (reinterpret_cast<const Halfedge*>(p_he));
  }

  /*! Get an incident halfedge (non-const version). */
  Halfedge* halfedge ()
  {
    return (reinterpret_cast<Halfedge*>(p_he));
  }

  /*! Set an incident halfedge. */
  void set_halfedge (Halfedge* he)
  { 
    p_he = he;
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

  /*! Default constructor. */
  Arr_halfedge()
  {}

  /*! Get the opposite halfedge (const version). */
  const Halfedge* opposite () const
  { 
    return (reinterpret_cast<const Halfedge*>(p_opp));
  }
  
  /*! Get the opposite halfedge (non-const version). */
  Halfedge* opposite ()
  { 
    return (reinterpret_cast<Halfedge*>(p_opp));
  }

  /*! Sets the opposite halfedge. */
  void set_opposite (Halfedge* he) 
  { 
    p_opp = he;
  }

  /*! Get the previous halfedge along the chain (const version). */
  const Halfedge* prev () const
  {
    return (reinterpret_cast<const Halfedge*>(p_prev));
  }

  /*! Get the previous halfedge along the chain (const version). */
  Halfedge* prev ()
  {
    return (reinterpret_cast<Halfedge*>(p_prev));
  }

  /*! Set the previous halfedge along the chain. */
  void set_prev (Halfedge* he)
  {
    p_prev = he;
    he->p_next = this;
  }

  /*! Get the next halfedge along the chain (const version). */
  const Halfedge* next () const
  {
    return (reinterpret_cast<const Halfedge*>(p_next));
  }

  /*! Get the next halfedge along the chain (const version). */
  Halfedge* next ()
  {
    return (reinterpret_cast<Halfedge*>(p_next));
  }

  /*! Set the next halfedge along the chain. */
  void set_next (Halfedge* he)
  {
    p_next = he;
    he->p_prev = this;
  }

  /*! Get the target vertex (const version). */
  const Vertex* vertex () const 
  { 
    return (reinterpret_cast<const Vertex*>(p_v));
  }

  /*! Get the target vertex (non-const version). */
  Vertex* vertex ()
  { 
    return (reinterpret_cast<Vertex*>(p_v));
  }

  /*! Set the target vertex. */
  void set_vertex (Vertex* v)
  {
    p_v = v;
  }

  /*! Get the incident face (const version). */
  const Face* face () const 
  {     
    return (reinterpret_cast<const Face*>(p_f));
  }

  /*! Get the incident face (non-const version). */
  Face* face () 
  {
    return (reinterpret_cast<Face*>(p_f));
  }

  /*! Set the incident face. */
  void set_face (Face* f)
  { 
    p_f = f;
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

  typedef F                     Base;
  typedef Arr_vertex<V,H,F>     Vertex;
  typedef Arr_halfedge<V,H,F>   Halfedge;
  typedef Arr_face<V,H,F>       Face;

  /*! Default constructor. */
  Arr_face()
  {}

  /*! Get an incident halfedge (const version). */
  const Halfedge * halfedge() const
  {
    return (reinterpret_cast<const Halfedge*>(p_he));
  }

  /*! Get an incident halfedge (non-const version). */
  Halfedge * halfedge()
  {
    return (reinterpret_cast<Halfedge*>(p_he));
  }

  /*! Set an incident halfedge. */
  void set_halfedge (Halfedge* he)
  {
    p_he = he;
  }

  // Define the hole iterators:
  typedef I_HalfedgeDS_iterator<
    typename F::Holes_iterator,
    Halfedge*,
    typename F::Holes_iterator::difference_type,
    typename F::Holes_iterator::iterator_category>       Holes_iterator;

  typedef I_HalfedgeDS_const_iterator<
    typename F::Holes_const_iterator,
    typename F::Holes_iterator,
    const Halfedge*,
    typename F::Holes_const_iterator::difference_type,
    typename F::Holes_const_iterator::iterator_category> Holes_const_iterator;

  /*! Add a hole inside the face. */
  void add_hole (Halfedge* h)
  {
    holes.push_back (h);
  }

  /*! Erase a hole from the face. */
  void erase_hole (Holes_iterator hit)
  {
    holes.erase (hit.current_iterator());
  }

  /*! Get an iterator for the first hole inside the face. */
  Holes_iterator holes_begin()
  {
    return holes.begin();
  }

  /*! Get a past-the-end iterator for the holes inside the face. */
  Holes_iterator holes_end()
  {
    return holes.end();
  }

  /*! Get an const iterator for the first hole inside the face. */
  Holes_const_iterator holes_begin() const
  {
    return holes.begin();
  }

  /*! Get a const past-the-end iterator for the holes inside the face. */
  Holes_const_iterator holes_end() const
  {
    return holes.end();
  }

  // Define the isloated vertices iterators:
  typedef I_Dereference_iterator<
    typename F::Isolated_vertices_iterator,
    Vertex,
    typename F::Isolated_vertices_iterator::difference_type,
    typename F::Isolated_vertices_iterator::iterator_category>
                                              Isolated_vertices_iterator;
  
  typedef I_Dereference_const_iterator<
    typename F::Isolated_vertices_const_iterator,
    typename F::Isolated_vertices_iterator,
    Vertex,
    typename F::Isolated_vertices_iterator::difference_type,
    typename F::Isolated_vertices_iterator::iterator_category>
                                              Isolated_vertices_const_iterator;

  /*! Add an isloated vertex inside the face. */
  void add_isolated_vertex (Vertex* v)
  {
    iso_verts.push_back (v);
  }

  /*! Erase an isloated vertex from the face. */
  void erase_isolated_vertex (Isolated_vertices_iterator ivit)
  {
    iso_verts.erase (ivit.current_iterator());
  }

  /*! Get an iterator for the first isloated vertex inside the face. */
  Isolated_vertices_iterator isolated_vertices_begin()
  {
    return iso_verts.begin();
  }

  /*! Get a past-the-end iterator for the isloated vertices inside the face. */
  Isolated_vertices_iterator isolated_vertices_end()
  {
    return iso_verts.end();
  }

  /*! Get an const iterator for the first isloated vertex inside the face. */
  Isolated_vertices_const_iterator isolated_vertices_begin() const
  {
    return iso_verts.begin();
  }

  /*! Get a const past-the-end iterator for the isloated vertices inside the
   * face. */
  Isolated_vertices_const_iterator isolated_vertices_end() const
  {
    return iso_verts.end();
  }
};

/*! \class
 * The arrangement DCEL class.
 */
template <class V, class H, class F,
          class Allocator = CGAL_ALLOCATOR(int) >
class Arr_dcel
{
public:

  // Define the vertex, halfedge and face types.
  typedef Arr_dcel<V,H,F>       Self;
  typedef Arr_vertex<V,H,F>     Vertex;
  typedef Arr_halfedge<V,H,F>   Halfedge;
  typedef Arr_face<V,H,F>       Face;
  
protected:

  // The vetices, halfedges and faces are stored in three in-place lists.
  typedef In_place_list<Vertex, false>   Vertex_list;
  typedef In_place_list<Halfedge, false> Halfedge_list;
  typedef In_place_list<Face, false>     Face_list;

  // Vertex allocator.
  typedef typename Allocator::template rebind<Vertex>    Vertex_alloc_rebind;
  typedef typename Vertex_alloc_rebind::other            Vertex_allocator;

  // Halfedge allocator.
  typedef typename Allocator::template rebind<Halfedge>  Halfedge_alloc_rebind;
  typedef typename Halfedge_alloc_rebind::other          Halfedge_allocator;

  // Face allocator.
  typedef typename Allocator::template rebind<Face>      Face_alloc_rebind;
  typedef typename Face_alloc_rebind::other              Face_allocator;

public:

  typedef typename Halfedge_list::size_type             Size;
  typedef typename Halfedge_list::size_type             size_type;
  typedef typename Halfedge_list::difference_type       difference_type;
  typedef typename Halfedge_list::difference_type       Difference;
  typedef std::bidirectional_iterator_tag               iterator_category;

protected:

  Vertex_list         vertices;             // The vertices container.
  Halfedge_list       halfedges;            // The halfedges container.
  Face_list           faces;                // The faces container.

  Vertex_allocator    vertex_alloc;         // An allocator for vertices.
  Halfedge_allocator  halfedge_alloc;       // An allocator for halfedges.
  Face_allocator      face_alloc;           // An allocator for faces.

public:

  // Definitions of iterators.
  typedef typename Vertex_list::iterator              Vertex_iterator;
  typedef typename Halfedge_list::iterator            Halfedge_iterator;
  typedef typename Face_list::iterator                Face_iterator;
  typedef CGAL::N_step_adaptor<Halfedge_iterator, 2>  Edge_iterator;
  
  // Definitions of const iterators.
  typedef typename Vertex_list::const_iterator        Vertex_const_iterator;
  typedef typename Halfedge_list::const_iterator      Halfedge_const_iterator;
  typedef typename Face_list::const_iterator          Face_const_iterator;
  typedef CGAL::N_step_adaptor<Halfedge_const_iterator,
                               2>                     Edge_const_iterator;
private:

  // Copy Constructor - not supported.
  Arr_dcel (const Self& ) ;

  // Assignment operator - not supported.
  Self& operator= (const Self& );

public:

  /// \name Construction and destruction.
  //@{
  /*! Default constructor. */
  Arr_dcel ()
  {}
  
  /*! Destructor. */
  ~Arr_dcel ()
  {
    delete_all();
  }
  //@}

  /// \name The DCEL size.
  //@{
  /*! Get the number of DCEL vertices. */
  Size size_of_vertices () const
  { 
    return (vertices.size());
  }

  /*! Get the number of DCEL halfedges (twice the number of edges). */
  Size size_of_halfedges () const
  {
    return (halfedges.size());
  }

  /*! Get the number of DCEL faces. */
  Size size_of_faces() const
  {
    return (faces.size());
  }
  //@}

  /// \name Obtaining iterators.
  //@{
  Vertex_iterator   vertices_begin()  { return vertices.begin(); }
  Vertex_iterator   vertices_end()    { return vertices.end(); }
  Halfedge_iterator halfedges_begin() { return halfedges.begin();}
  Halfedge_iterator halfedges_end()   { return halfedges.end(); }
  Face_iterator     faces_begin()     { return faces.begin(); }
  Face_iterator     faces_end()       { return faces.end(); }
  Edge_iterator     edges_begin()     { return halfedges.begin(); }  
  Edge_iterator     edges_end()       { return halfedges.end(); }  
  //@}

  /// \name Obtaining constant iterators.
  //@{
  Vertex_const_iterator   vertices_begin() const { return vertices.begin(); }
  Vertex_const_iterator   vertices_end() const { return vertices.end(); }
  Halfedge_const_iterator halfedges_begin() const { return halfedges.begin(); }
  Halfedge_const_iterator halfedges_end() const { return halfedges.end(); }
  Face_const_iterator     faces_begin() const { return faces.begin(); }
  Face_const_iterator     faces_end() const { return faces.end(); }
  Edge_const_iterator     edges_begin() const { return halfedges.begin(); }  
  Edge_const_iterator     edges_end() const { return halfedges.end(); }  
  //@}

  // \name Creation of new DCEL features.
  //@{
  /*! Create a new vertex. */
  Vertex* new_vertex()
  {
    Vertex     *v = vertex_alloc.allocate (1);

    vertex_alloc.construct (v, Vertex());
    vertices.push_back (*v);
    return v;
  }
  
  /*! Create a new pair of opposite halfedges. */
  Halfedge* new_edge() 
  {
    // Create two new halfedges.
    Halfedge   *h1 = _new_halfedge ();
    Halfedge   *h2 = _new_halfedge ();

    // Pair them together.
    h1->set_opposite (h2);
    h2->set_opposite (h1);

    return (h1);
  }

  /*! Create a new face. */
  Face* new_face()
  {
    Face       *f = face_alloc.allocate (1);
    
    face_alloc.construct (f, Face());
    faces.push_back (*f);
    return (f);
  }
  //@}

  /// \name Deletion of DCEL features.
  //@{
  /*! Delete an existing vertex. */
  void delete_vertex (Vertex * v)
  {
    vertices.erase (v);
    vertex_alloc.destroy (v);
    vertex_alloc.deallocate (v,1);
  }
  
  /*! Delete an existing pair of opposite halfedges. */
  void delete_edge (Halfedge * h) 
  {
    Halfedge   *h_opp = h->opposite();

    _delete_halfedge (h);
    _delete_halfedge (h_opp);
  }

  /*! Delete an existing face. */
  void delete_face(Face * f)
  {
    faces.erase (f);
    face_alloc.destroy (f);
    face_alloc.deallocate (f, 1);
  }
  
  /*! Delete all DCEL features. */
  void delete_all() 
  {
    // Free all vertices.
    Vertex_iterator    vit = vertices.begin(), v_curr;

    while (vit != vertices.end())
    {
      v_curr = vit;
      ++vit;
      delete_vertex (&(*v_curr));
    }

    // Free all halfedges.
    Halfedge_iterator  hit = halfedges.begin(), h_curr;

    while (hit != halfedges.end())
    {
      h_curr = hit;
      ++hit;
      _delete_halfedge (&(*h_curr));
    }

    // Free all faces.
    Face_iterator      fit = faces.begin(), f_curr;

    while (fit != faces.end())
    {
      f_curr = fit;
      ++fit;
      delete_face (&(*f_curr));
    }
  }
  //@}

  /*!
   * Assign our DCEL the contents of another DCEL.
   * \param dcel The DCEL to be copied.
   * \param uf A pointer to the unbounded face in this DCEL.
   * \return A pointer to the unbounded face in the copied DCEL.
   */
  Face* assign (const Self& dcel, const Face *uf)
  {
    // Clear the current contents of the DCEL.
    delete_all();

    // Create duplicated of the DCEL features and map the features of the
    // given DCEL to their corresponding duplicates.
    typedef std::map<const Vertex*, Vertex*>     Vertex_map;
    typedef std::map<const Halfedge*, Halfedge*> Halfedge_map;
    typedef std::map<const Face*, Face*>         Face_map;

    Vertex_map                v_map;
    Vertex_const_iterator     vit;
    Vertex                   *dup_v;

    for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
    {
      dup_v = new_vertex();
      dup_v->assign (*vit);
      v_map.insert (Vertex_map::value_type (&(*vit), dup_v));
    }

    Halfedge_map              he_map;
    Halfedge_const_iterator   hit;
    Halfedge                 *dup_h;

    for (hit = dcel.halfedges_begin(); hit != dcel.halfedges_end(); ++hit)
    {
      dup_h = _new_halfedge();
      dup_h->assign (*hit);
      he_map.insert (Halfedge_map::value_type(&(*hit), dup_h));
    }

    Face_map                  f_map;
    Face_const_iterator       fit;
    Face                     *dup_f;

    for (fit = dcel.faces_begin(); fit != dcel.faces_end(); ++fit)
    {
      dup_f = new_face();
      dup_f->assign (*fit);
      f_map.insert (Face_map::value_type(&(*fit), dup_f));
    }

    // Update the vertex records.
    const Vertex             *v;
    const Halfedge           *h;
    const Face               *f;
    
    for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
    {
      v = &(*vit);
      h = v->halfedge();

      dup_v = (v_map.find (v))->second;
      dup_h = (he_map.find (h))->second;

      dup_v->set_halfedge (dup_h);
    }

    // Update the halfedge records.
    const Halfedge           *opp, *prev, *next;
    Halfedge                 *dup_opp, *dup_prev, *dup_next;

    for (hit = dcel.halfedges_begin(); hit != dcel.halfedges_end(); ++hit)
    {
      h = &(*hit);
      v = h->vertex();
      f = h->face();
      opp = h->opposite();
      prev = h->prev();
      next = h->next();

      dup_h = (he_map.find (h))->second;
      dup_v = (v_map.find (v))->second;
      dup_f = (f_map.find (f))->second;
      dup_opp = (he_map.find (opp))->second;
      dup_prev = (he_map.find (prev))->second;
      dup_next = (he_map.find (next))->second;

      dup_h->set_vertex (dup_v);
      dup_h->set_face (dup_f);
      dup_h->set_opposite (dup_opp);
      dup_h->set_prev (dup_prev);
      dup_h->set_next (dup_next);
    }

    // Update the face records.
    typename Face::Holes_const_iterator              holes_it;
    typename Face::Isolated_vertices_const_iterator  iso_verts_it;
    const Halfedge                      *hole;
    const Vertex                        *iso_vert;
    Halfedge                            *dup_hole;
    Vertex                              *dup_iso_vert;

    for (fit = dcel.faces_begin(); fit != dcel.faces_end(); ++fit)
    {
      f = &(*fit);
      h = f->halfedge();

      // Set the pointer to the outer boundary edge (may be NULL in case that
      // the current face f is the unbounded face).
      dup_f = (f_map.find (f))->second;
      if (h != NULL)
        dup_h = (he_map.find (h))->second;
      else
        dup_h = NULL;

      dup_f->set_halfedge (dup_h);

      // Assign the holes.
      for (holes_it = f->holes_begin();
           holes_it != f->holes_end(); ++holes_it)
      {
        hole = *holes_it;

        dup_hole = (he_map.find (hole))->second;
        dup_f->add_hole (dup_hole);
      }

      // Assign the isolated vertices.
      for (iso_verts_it = f->isolated_vertices_begin();
           iso_verts_it != f->isolated_vertices_end(); ++iso_verts_it)
      {
        iso_vert = &(*iso_verts_it);

        dup_iso_vert = (v_map.find (iso_vert))->second;
        dup_f->add_isolated_vertex (dup_iso_vert);
      }
    }

    // Return the unbounded face in the copied DCEL.
    return ((f_map.find (uf))->second);
  }

protected:

  /*! Create a new halfedge. */
  Halfedge * _new_halfedge ()
  {
    Halfedge   *h = halfedge_alloc.allocate (1);

    halfedge_alloc.construct (h, Halfedge());
    halfedges.push_back (*h);
    return (h);
  }

  /*! Delete an existing halfedge. */
  void _delete_halfedge (Halfedge* h)
  {
    halfedges.erase (h);
    halfedge_alloc.destroy (h);
    halfedge_alloc.deallocate (h, 1);
  }
};

/*! \class
 * The default arrangement DCEL class.
 * The Traits parameters corresponds to a geometric traits class, which 
 * defines the Point_2 and X_monotone_curve_2 types.
 */
template <class Traits_>
class Arr_default_dcel :
  public Arr_dcel<Arr_vertex_base<typename Traits_::Point_2>,
                  Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
                  Arr_face_base>
{
public:

  /*! Default constructor. */
  Arr_default_dcel()
  {}
};

CGAL_END_NAMESPACE

#endif 
