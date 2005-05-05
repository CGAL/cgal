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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_2_ITERATORS_H
#define CGAL_ARRANGEMENT_2_ITERATORS_H

/*! \file
 * Definitions of the arrangement iterator and circulators.
 */

CGAL_BEGIN_NAMESPACE

// Forward declerations.
template <class Traits, class Dcel> class _Vertex_handle;
template <class Traits, class Dcel> class _Vertex_const_handle;
template <class Traits, class Dcel> class _Vertex_iterator;
template <class Traits, class Dcel> class _Vertex_const_iterator;

template <class Traits, class Dcel> class _Halfedge_handle;
template <class Traits, class Dcel> class _Halfedge_const_handle;
template <class Traits, class Dcel> class _Halfedge_iterator;
template <class Traits, class Dcel> class _Halfedge_const_iterator;
template <class Traits, class Dcel> class _Edge_iterator;
template <class Traits, class Dcel> class _Edge_const_iterator;

template <class Traits, class Dcel> class _Face_handle;
template <class Traits, class Dcel> class _Face_const_handle;
template <class Traits, class Dcel> class _Face_iterator;
template <class Traits, class Dcel> class _Face_const_iterator;

template <class Traits, class Dcel> class _Halfedge_around_vertex_circulator;
template <class Traits, class Dcel> class 
                                    _Halfedge_around_vertex_const_circulator;
template <class Traits, class Dcel> class _Ccb_halfedge_circulator;
template <class Traits, class Dcel> class _Ccb_halfedge_const_circulator;
template <class Traits, class Dcel> class _Holes_iterator;
template <class Traits, class Dcel> class _Holes_const_iterator;

template <class Traits, class Dcel>  class Arrangement_2;

// --------------------------------------------------------------------------
// Circulators
// --------------------------------------------------------------------------

/*! \class
 * Circulator for the edges around a given vertex.
 */
template <class Traits_, class Dcel_> 
class _Halfedge_around_vertex_circulator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Halfedge_around_vertex_const_circulator<Traits_, Dcel_>;
  friend class _Vertex_handle<Traits_, Dcel_>;

protected:

  typedef Dcel_                                             Dcel;
  typedef _Halfedge_around_vertex_circulator<Traits_, Dcel> Self;
  typedef typename Dcel::Halfedge                           Halfedge;

  Halfedge    *p_he;           // The current halfedge.

  /*! Constructor from a halfedge pointer. */
  _Halfedge_around_vertex_circulator (Halfedge *e) :
    p_he (e)
  {}

public:

  /*! Default constructor. */
  _Halfedge_around_vertex_circulator () :
    p_he (NULL)
  {}

  /*! Equality operator. */
  bool operator== (const Self& circ) const
  {
    return (p_he == circ.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& circ) const
  {
    return (p_he != circ.p_he);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    p_he = p_he->next()->opposite();
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    p_he = p_he->next()->opposite();
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    p_he = p_he->opposite()->previous();
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    p_he = p_he->opposite()->previous();
    return (temp);
  }

  /*! Get the halfedge pointer. */
  Halfedge* halfedge () const
  {
    return (p_he);
  }
};

/*! \class
 * Const circulator for the edges around a given vertex.
 */
template <class Traits_, class Dcel_> 
class _Halfedge_around_vertex_const_circulator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Vertex_const_handle<Traits_, Dcel_>;

protected:

  typedef Dcel_                                                   Dcel;
  typedef _Halfedge_around_vertex_const_circulator<Traits_, Dcel> Self;
  typedef typename Dcel::Halfedge                                 Halfedge;

  const Halfedge    *p_he;           // The current halfedge.

  /*! Constructor from a halfedge pointer. */
  _Halfedge_around_vertex_const_circulator (const Halfedge *e) :
    p_he (e)
  {}

public:

  /*! Default constructor. */
  _Halfedge_around_vertex_const_circulator () :
    p_he (NULL)
  {}

  /*! Constructor from a non-const circulator. */
  _Halfedge_around_vertex_const_circulator
      (const _Halfedge_around_vertex_circulator<Traits_, Dcel>& circ) :
    p_he (circ.p_he)
  {}

  /*! Equality operator. */
  bool operator== (const Self& circ) const
  {
    return (p_he == circ.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& circ) const
  {
    return (p_he != circ.p_he);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    p_he = p_he->next()->opposite();
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    p_he = p_he->next()->opposite();
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    p_he = p_he->opposite()->previous();
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    p_he = p_he->opposite()->previous();
    return (temp);
  }

  /*! Get the halfedge pointer. */
  const Halfedge* halfedge () const
  {
    return (p_he);
  }
};

/*! \class
 * Circulator for the edges around the boundary of a connected component.
 */
template <class Traits_, class Dcel_> 
class _Ccb_halfedge_circulator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Ccb_halfedge_const_circulator<Traits_, Dcel_>;
  friend class _Face_handle<Traits_, Dcel_>;
  friend class _Holes_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                    Dcel;
  typedef _Ccb_halfedge_circulator<Traits_, Dcel>  Self;
  typedef typename Dcel::Halfedge                  Halfedge;

  Halfedge    *p_he;           // The current halfedge.

  /*! Constructor from a halfedge pointer. */
  _Ccb_halfedge_circulator (Halfedge *e) :
    p_he (e)
  {}

public:

  /*! Default constructor. */
  _Ccb_halfedge_circulator () :
    p_he (NULL)
  {}

  /*! Equality operator. */
  bool operator== (const Self& circ) const
  {
    return (p_he == circ.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& circ) const
  {
    return (p_he != circ.p_he);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    p_he = p_he->next();
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    p_he = p_he->next();
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    p_he = p_he->previous();
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    p_he = p_he->previous();
    return (temp);
  }

  /*! Get the halfedge pointer. */
  Halfedge* halfedge () const
  {
    return (p_he);
  }
};

/*! \class
 * Const circulator for the edges around the boundary of a connected component.
 */
template <class Traits_, class Dcel_> 
class _Ccb_halfedge_const_circulator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Face_const_handle<Traits_, Dcel_>;
  friend class _Holes_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                          Dcel;
  typedef _Ccb_halfedge_const_circulator<Traits_, Dcel>  Self;
  typedef typename Dcel::Halfedge                        Halfedge;

  const Halfedge    *p_he;           // The current halfedge.

  /*! Constructor from a halfedge pointer. */
  _Ccb_halfedge_const_circulator (const Halfedge *e) :
    p_he (e)
  {}

public:

  /*! Default constructor. */
  _Ccb_halfedge_const_circulator () :
    p_he (NULL)
  {}

  /*! Constructor from a non-const circulator. */
  _Ccb_halfedge_const_circulator
      (const _Ccb_halfedge_circulator<Traits_, Dcel>& circ) :
    p_he (circ.p_he)
  {}

  /*! Equality operator. */
  bool operator== (const Self& circ) const
  {
    return (p_he == circ.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& circ) const
  {
    return (p_he != circ.p_he);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    p_he = p_he->next();
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    p_he = p_he->next();
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    p_he = p_he->previous();
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    p_he = p_he->previous();
    return (temp);
  }

  /*! Get the halfedge pointer. */
  const Halfedge* halfedge () const
  {
    return (p_he);
  }
};

/*! \class
 * Iterator for the holes inside a face.
 */
template <class Traits_, class Dcel_> 
class _Holes_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Holes_const_iterator<Traits_, Dcel_>;   
  friend class _Face_handle<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Holes_iterator<Traits_, Dcel>         Self;
  typedef typename Dcel::Face::Holes_iterator    Holes_iter;

  Holes_iter     holes_it;           // The stored iterator.

  /*! Constructor from a DCEL holes iterator. */
  _Holes_iterator (const Holes_iter& it) :
    holes_it (it)
  {}

public:

public:

  typedef unsigned int                            Size;
  typedef unsigned int                            size_type;
  typedef typename Holes_iter::difference_type    difference_type;
  typedef typename Holes_iter::difference_type    Difference;
  typedef std::bidirectional_iterator_tag         iterator_category;

  /*! Default constructor. */
  _Holes_iterator ()
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (holes_it == iter.holes_it);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (holes_it != iter.holes_it);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++holes_it;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++holes_it;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --holes_it;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --holes_it;
    return (temp);
  }

  /*! Get a circulator for the boundary of the current hole. */
  _Ccb_halfedge_circulator<Traits_, Dcel> operator* () const
  {
    _Ccb_halfedge_circulator<Traits_, Dcel>     circ (*holes_it);
    return (circ);
  }
};

/*! \class
 * Const iterator for the holes inside a face.
 */
template <class Traits_, class Dcel_>
class _Holes_const_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Face_const_handle<Traits_, Dcel_>;

protected:

  typedef Dcel_                                      Dcel;
  typedef _Holes_const_iterator<Traits_, Dcel>       Self;
  typedef typename Dcel::Face::Holes_const_iterator  Holes_const_iter;

  Holes_const_iter     holes_it;           // The stored iterator.

  /*! Constructor from a DCEL holes iterator. */
  _Holes_const_iterator (const Holes_const_iter& it) :
    holes_it (it)
  {}

public:

  typedef unsigned int                               Size;
  typedef unsigned int                               size_type;
  typedef typename Holes_const_iter::difference_type difference_type;
  typedef typename Holes_const_iter::difference_type Difference;
  typedef std::bidirectional_iterator_tag            iterator_category;

  /*! Default constructor. */
  _Holes_const_iterator ()
  {}

  /*! Constructor from a non-const iterator. */
  _Holes_const_iterator (const _Holes_iterator<Traits_, Dcel>& iter) :
    holes_it (iter.holes_it)
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (holes_it == iter.holes_it);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (holes_it != iter.holes_it);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++holes_it;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++holes_it;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --holes_it;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --holes_it;
    return (temp);
  }

  /*! Get a const circulator for the boundary of the current hole. */
  _Ccb_halfedge_const_circulator<Traits_, Dcel> operator* () const
  {
    _Ccb_halfedge_const_circulator<Traits_, Dcel>     circ (*holes_it);
    return (circ);  
  }
};

// --------------------------------------------------------------------------
// Handles
// --------------------------------------------------------------------------

/*! \class
 * The vertex handle class.
 */
template <class Traits_, class Dcel_>
class _Vertex_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Vertex_const_handle<Traits_, Dcel_>;
  friend class _Halfedge_handle<Traits_, Dcel_>;
  friend class _Vertex_iterator<Traits_, Dcel_>;

protected:

  typedef Traits_                         Traits;
  typedef Dcel_                           Dcel;
  typedef _Vertex_handle<Traits, Dcel>    Self;
  typedef typename Dcel::Size             Size;
  typedef typename Dcel::Vertex           Vertex;
  typedef typename Dcel::Halfedge         Halfedge;

  Vertex        *p_v;     // The vertex associated with the handle.

  /*! Constructor from a vertex pointer. */
  _Vertex_handle (Vertex *v) :
    p_v (v)
  {}

public:

  /*! Default constructor. */
  _Vertex_handle () :
    p_v (NULL)
  {}

  /*! Equality operator. */
  bool operator== (const Self& vh) const
  {
    return (p_v == vh.p_v);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& vh) const
  {
    return (p_v != vh.p_v);
  }

  /*! Get the point associated with the vertex. */
  const typename Traits::Point_2& point () const
  {
    return (p_v->point());
  }

  /*!
   * Set the point associated with the vertex.
   * \param p The point.
   * \pre p must be geometrically equal to the current point the
   *      vertex stores.
   */
  void set_point (const typename Traits::Point_2& p)
  {
    CGAL_precondition_code(
      Traits     traits;
    );
    CGAL_precondition (traits.points_equal (p_v->point(), p));

    // Set the point (not that point() returns a Point_2 reference).
    p_v->point() = p;
  }

  /*! Get the vertex degree (the number of incident halfedges). */
  Size degree () const
  {
    // Go over all incident halfedges and count them.
    Halfedge   *first = p_v->halfedge();
    Halfedge   *curr = first;
    Size        deg = 0;

    if (first == NULL)
      return (deg);

    do
    {
      deg++;
      curr = curr->next()->opposite();
    } while (curr != first);

    return (deg);
  }

  /*! Get a circulator for the incident halfedges. */
  _Halfedge_around_vertex_circulator<Traits, Dcel> incident_halfedges () const
  {
    _Halfedge_around_vertex_circulator<Traits, Dcel>   circ (p_v->halfedge());
    return (circ);
  } 
};

/*! \class
 * The vertex const handle class.
 */
template <class Traits_, class Dcel_>
class _Vertex_const_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Halfedge_const_handle<Traits_, Dcel_>;
  friend class _Vertex_const_iterator<Traits_, Dcel_>;

protected:

  typedef Traits_                             Traits;
  typedef Dcel_                               Dcel;
  typedef _Vertex_const_handle<Traits, Dcel>  Self;
  typedef typename Dcel::Size                 Size;
  typedef typename Dcel::Vertex               Vertex;
  typedef typename Dcel::Halfedge             Halfedge;

  const Vertex        *p_v;     // The vertex associated with the handle.

  /*! Constructor from a vertex pointer. */
  _Vertex_const_handle (const Vertex *v) :
    p_v (v)
  {}

public:

  /*! Default constructor. */
  _Vertex_const_handle () :
    p_v (NULL)
  {}

  /* Constructor from a non-const handle. */
  _Vertex_const_handle (const _Vertex_handle<Traits, Dcel>& vh) :
    p_v (vh.p_v)
  {}

  /*! Equality operator. */
  bool operator== (const Self& vh) const
  {
    return (p_v == vh.p_v);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& vh) const
  {
    return (p_v != vh.p_v);
  }

  /*! Get the point associated with the vertex. */
  const typename Traits::Point_2& point () const
  {
    return (p_v->point());
  }

  /*! Get the vertex degree (the number of incident halfedges). */
  Size degree () const
  {
    // Go over all incident halfedges and count them.
    const Halfedge     *first = p_v->halfedge();
    const Halfedge     *curr = first;
    Size                deg = 0;

    if (first == NULL)
      return (deg);

    do
    {
      deg++;
      curr = curr->next()->opposite();
    } while (curr != first);

    return (deg);
  }

  /*! Get a circulator for the incident halfedges. */
  _Halfedge_around_vertex_const_circulator<Traits, Dcel>
    incident_halfedges () const
  {
    _Halfedge_around_vertex_const_circulator<Traits, Dcel>
      circ (p_v->halfedge());
    
    return (circ);
  } 
};

/*! \class
 * The face handle class.
 */
template <class Traits_, class Dcel_>
class _Face_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Face_const_handle<Traits_, Dcel_>;
  friend class _Halfedge_handle<Traits_, Dcel_>;
  friend class _Face_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                           Dcel;
  typedef _Face_handle<Traits_, Dcel>     Self;
  typedef typename Dcel::Face             Face;

  Face          *p_f;     // The face associated with the handle.

  /*! Constructor from a face pointer. */
  _Face_handle (Face *f) :
    p_f (f)
  {}

public:

  /*! Default constructor. */
  _Face_handle () :
    p_f (NULL)
  {}

  /*! Equality operator. */
  bool operator== (const Self& fh) const
  {
    return (p_f == fh.p_f);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& fh) const
  {
    return (p_f != fh.p_f);
  }

  /*! Check whether the face is unbounded. */
  bool is_unbounded () const
  {
    // Check whether the outer-boundary edge exists or not.
    return (p_f->halfedge() == NULL);
  }

  /*! Get a circulator for the outer boundary of the face. */
  _Ccb_halfedge_circulator<Traits_, Dcel> outer_ccb () const
  {
    _Ccb_halfedge_circulator<Traits_, Dcel>   circ (p_f->halfedge());
    return (circ);
  }

  /*! Get an iterator for the first hole inside the face. */
  _Holes_iterator<Traits_, Dcel> holes_begin () const
  {
    _Holes_iterator<Traits_, Dcel>            iter (p_f->holes_begin());
    return (iter);
  }

  /*! Get a past-the-end iterator for holes inside the face. */
  _Holes_iterator<Traits_, Dcel> holes_end () const
  {
    _Holes_iterator<Traits_, Dcel>            iter (p_f->holes_end());
    return (iter);
  }
};

/*! \class
 * The face const handle class.
 */
template <class Traits_, class Dcel_>
class _Face_const_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Halfedge_const_handle<Traits_, Dcel_>;
  friend class _Face_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                              Dcel;
  typedef _Face_const_handle<Traits_, Dcel>  Self;
  typedef typename Dcel::Face                Face;

  const Face          *p_f;     // The face associated with the handle.

  /*! Constructor from a face pointer. */
  _Face_const_handle (const Face *f) :
    p_f (f)
  {}

public:

  /*! Default constructor. */
  _Face_const_handle () :
    p_f (NULL)
  {}

  /*! Constructor from a non-const handle. */
  _Face_const_handle (const _Face_handle<Traits_, Dcel>& fh) :
    p_f (fh.p_f)
  {}

  /*! Equality operator. */
  bool operator== (const Self& fh) const
  {
    return (p_f == fh.p_f);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& fh) const
  {
    return (p_f != fh.p_f);
  }

  /*! Check whether the face is unbounded. */
  bool is_unbounded () const
  {
    // Check whether the outer-boundary edge exists or not.
    return (p_f->halfedge() == NULL);
  }

  /*! Get a circulator for the outer boundary of the face. */
  _Ccb_halfedge_const_circulator<Traits_, Dcel> outer_ccb () const
  {
    _Ccb_halfedge_const_circulator<Traits_, Dcel> circ (p_f->halfedge());
    return (circ);
  }

  /*! Get an iterator for the first hole inside the face. */
  _Holes_const_iterator<Traits_, Dcel> holes_begin () const
  {
    _Holes_const_iterator<Traits_, Dcel>          iter (p_f->holes_begin());
    return (iter);
  }

  /*! Get a past-the-end iterator for holes inside the face. */
  _Holes_const_iterator<Traits_, Dcel> holes_end () const
  {
    _Holes_const_iterator<Traits_, Dcel>          iter (p_f->holes_end());
    return (iter);
  }
};


/*! \class
 * The halfedge handle class.
 */
template <class Traits_, class Dcel_>
class _Halfedge_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Halfedge_const_handle<Traits_, Dcel_>;

protected:

  typedef Traits_                         Traits;
  typedef Dcel_                           Dcel;
  typedef _Halfedge_handle<Traits, Dcel>  Self;
  typedef typename Dcel::Halfedge         Halfedge;

  Halfedge          *p_he;     // The halfedge associated with the handle.

public:

  /*! Default constructor. */
  _Halfedge_handle () :
    p_he (NULL)
  {}

  /*! Constructor from a halfedge pointer. */
  _Halfedge_handle (Halfedge *e) :
    p_he (e)
  {}

  /*! Equality operator. */
  bool operator== (const Self& hh) const
  {
    return (p_he == hh.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& hh) const
  {
    return (p_he != hh.p_he);
  }

  /*! Get the curve associated with the halfedge. */
  const typename Traits::X_monotone_curve_2& curve () const
  {
    return (p_he->curve());
  }

  /*!
   * Set the curve associated with the halfedge.
   * \param cv The x-monotone curve.
   * \pre cv must be geometrically equal to the current curve the
   *      halfedge stores.
   */
  void set_curve (const typename Traits::X_monotone_curve_2& cv)
  {
    CGAL_precondition_code(
      Traits     traits;
    );
    CGAL_precondition (traits.curves_equal (p_he->curve(), cv));

    // Set the curve (not that curve() returns an X_monotone_curve_2
    // reference). Also notice it is not necessary to set the curve
    // of the twin halfedge, as it stores a pointer to the same
    // X_monotone_curve_2 object we have just modified.
    p_v->curve() = cv;
  }

  /*! Get the source vertex. */
  _Vertex_handle<Traits, Dcel> source () const
  {
    // Return the incident vertex of the opposite halfedge.
    _Vertex_handle<Traits, Dcel>      vh (p_he->opposite()->vertex());
    return (vh);
  }

  /*! Get the target vertex. */
  _Vertex_handle<Traits, Dcel> target () const
  {
    // Return the incident vertex of this halfedge.
    _Vertex_handle<Traits, Dcel>      vh (p_he->vertex());
    return (vh);
  }

  /*! Get the incident face. */
  _Face_handle<Traits, Dcel> face () const
  {
    _Face_handle<Traits, Dcel>        fh (p_he->face());
    return (fh);
  }

  /*! Get the twin halfedge. */
  Self twin () const
  {
    Self                      hh (p_he->opposite());
    return (hh);
  }

  /*! Get the previous halfedge along the component boundary. */
  Self previous () const
  {
    Self                      hh (p_he->previous());
    return (hh);
  }

  /*! Get the next halfedge along the component boundary. */
  Self next () const
  {
    Self                      hh (p_he->next());
    return (hh);
  }
};

// Dereference operators:
template <class Traits, class Dcel>
_Halfedge_handle<Traits, Dcel> operator*
    (const _Halfedge_around_vertex_circulator<Traits, Dcel>& circ)
{
  _Halfedge_handle<Traits, Dcel>   hh (circ.halfedge());
  return (hh);
}

template <class Traits, class Dcel>
_Halfedge_handle<Traits, Dcel> operator* 
    (const _Ccb_halfedge_circulator<Traits, Dcel>& circ)
{
  _Halfedge_handle<Traits, Dcel>   hh (circ.halfedge());
  return (hh);
}

/*! \class
 * The halfedge const handle class.
 */
template <class Traits_, class Dcel_>
class _Halfedge_const_handle
{
  friend class Arrangement_2<Traits_, Dcel_>;

protected:

  typedef Traits_                               Traits;
  typedef Dcel_                                 Dcel;
  typedef _Halfedge_const_handle<Traits, Dcel>  Self;
  typedef typename Dcel::Halfedge               Halfedge;

  const Halfedge     *p_he;     // The halfedge associated with the handle.

public:

  /*! Default constructor. */
  _Halfedge_const_handle () :
    p_he (NULL)
  {}

  /*! Constructor from a non-const handle. */
  _Halfedge_const_handle (const _Halfedge_handle<Traits, Dcel>& hh) :
    p_he (hh.p_he)
  {}

  /*! Constructor from a halfedge pointer. */
  _Halfedge_const_handle (const Halfedge *e) :
    p_he (e)
  {}

  /*! Equality operator. */
  bool operator== (const Self& hh) const
  {
    return (p_he == hh.p_he);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& hh) const
  {
    return (p_he != hh.p_he);
  }

  /*! Get the curve associated with the halfedge. */
  const typename Traits::X_monotone_curve_2& curve () const
  {
    return (p_he->curve());
  }

  /*! Get the source vertex. */
  _Vertex_const_handle<Traits, Dcel> source () const
  {
    // Return the incident vertex of the opposite halfedge.
    _Vertex_const_handle<Traits, Dcel>      vh (p_he->opposite()->vertex());
    return (vh);
  }

  /*! Get the target vertex. */
  _Vertex_const_handle<Traits, Dcel> target () const
  {
    // Return the incident vertex of this halfedge.
    _Vertex_const_handle<Traits, Dcel>      vh (p_he->vertex());
    return (vh);
  }

  /*! Get the incident face. */
  _Face_const_handle<Traits, Dcel> face () const
  {
    _Face_const_handle<Traits, Dcel>        fh (p_he->face());
    return (fh);
  }

  /*! Get the twin halfedge. */
  Self twin () const
  {
    Self                      hh (p_he->opposite());
    return (hh);
  }

  /*! Get the previous halfedge along the component boundary. */
  Self previous () const
  {
    Self                      hh (p_he->previous());
    return (hh);
  }

  /*! Get the next halfedge along the component boundary. */
  Self next () const
  {
    Self                      hh (p_he->next());
    return (hh);
  }
};

// Dereference operators:
template <class Traits, class Dcel>
_Halfedge_const_handle<Traits, Dcel> operator* 
    (const _Halfedge_around_vertex_const_circulator<Traits, Dcel>& circ)
{
  _Halfedge_const_handle<Traits, Dcel>   hh (circ.halfedge());
  return (hh);
}

template <class Traits, class Dcel>
_Halfedge_const_handle<Traits, Dcel> operator*
    (const _Ccb_halfedge_const_circulator<Traits, Dcel>& circ)
{
  _Halfedge_const_handle<Traits, Dcel>   hh (circ.halfedge());
  return (hh);
}

// --------------------------------------------------------------------------
// Iterators
// --------------------------------------------------------------------------

/*! \class
 * Iterator for the DCEL vertices.
 */
template <class Traits_, class Dcel_> 
class _Vertex_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Vertex_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Vertex_iterator<Traits_, Dcel>        Self;
  typedef typename Dcel::Vertex                  Vertex;
  typedef typename Dcel::Vertex_iterator         Vertex_iter;

  Vertex_iter     vit;           // The stored iterator.

  /*! Constructor from a DCEL vertex iterator. */
  _Vertex_iterator (const Vertex_iter& it) :
    vit (it)
  {}

public:

  typedef typename Vertex_iter::size_type         Size;
  typedef typename Vertex_iter::size_type         size_type;
  typedef typename Vertex_iter::difference_type   difference_type;
  typedef typename Vertex_iter::difference_type   Difference;
  typedef std::bidirectional_iterator_tag         iterator_category;

  /*! Default constructor. */
  _Vertex_iterator ()
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (vit == iter.vit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (vit != iter.vit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++vit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++vit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --vit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --vit;
    return (temp);
  }

  /*! Get a handle to the current vertex. */
  _Vertex_handle<Traits_, Dcel> operator* () const
  {
    _Vertex_handle<Traits_, Dcel>     vh (&(*vit));
    return (vh);
  }

  /*! Get a pointer to the current vertex. */
  Vertex* operator-> () const
  {
    return (&(*vit));
  }
};

/*! \class
 * Const iterator for the DCEL vertices.
 */
template <class Traits_, class Dcel_>
class _Vertex_const_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Vertex_const_iterator<Traits_, Dcel>  Self;
  typedef typename Dcel::Vertex                  Vertex;
  typedef typename Dcel::Vertex_const_iterator   Vertex_const_iter;

  Vertex_const_iter     vit;           // The stored iterator.

  /*! Constructor from a DCEL vertex iterator. */
  _Vertex_const_iterator (const Vertex_const_iter& it) :
    vit (it)
  {}

public:

  typedef typename Vertex_const_iter::size_type       Size;
  typedef typename Vertex_const_iter::size_type       size_type;
  typedef typename Vertex_const_iter::difference_type difference_type;
  typedef typename Vertex_const_iter::difference_type Difference;
  typedef std::bidirectional_iterator_tag             iterator_category;

  /*! Default constructor. */
  _Vertex_const_iterator ()
  {}

  /*! Constructor from a non-const iterator. */
  _Vertex_const_iterator (const _Vertex_iterator<Traits_, Dcel>& iter) :
    vit (iter.vit)
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (vit == iter.vit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (vit != iter.vit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++vit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++vit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --vit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --vit;
    return (temp);
  }

  /*! Get a handle to the current vertex. */
  _Vertex_const_handle<Traits_, Dcel> operator* () const
  {
    _Vertex_const_handle<Traits_, Dcel>     vh (&(*vit));
    return (vh);
  }

  /*! Get a pointer to the current vertex. */
  const Vertex* operator-> () const
  {
    return (&(*vit));
  }
};

/*! \class
 * Iterator for the DCEL halfedges.
 */
template <class Traits_, class Dcel_> 
class _Halfedge_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Halfedge_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Halfedge_iterator<Traits_, Dcel>      Self;
  typedef typename Dcel::Halfedge                Halfedge;
  typedef typename Dcel::Halfedge_iterator       Halfedge_iter;

  Halfedge_iter     hit;           // The stored iterator.

  /*! Constructor from a DCEL halfedge iterator. */
  _Halfedge_iterator (const Halfedge_iter& it) :
    hit (it)
  {}

public:

  typedef typename Halfedge_iter::size_type       Size;
  typedef typename Halfedge_iter::size_type       size_type;
  typedef typename Halfedge_iter::difference_type difference_type;
  typedef typename Halfedge_iter::difference_type Difference;
  typedef std::bidirectional_iterator_tag         iterator_category;

  /*! Default constructor. */
  _Halfedge_iterator ()
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (hit == iter.hit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (hit != iter.hit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++hit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++hit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --hit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --hit;
    return (temp);
  }

  /*! Get a handle to the current halfedge. */
  _Halfedge_handle<Traits_, Dcel> operator* () const
  {
    _Vertex_handle<Traits_, Dcel>     hh (&(*hit));
    return (hh);
  }

  /*! Get a pointer to the current halfedge. */
  Halfedge* operator-> () const
  {
    return (&(*hit));
  }
};

/*! \class
 * Const iterator for the DCEL halfedges.
 */
template <class Traits_, class Dcel_>
class _Halfedge_const_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

protected:

  typedef Dcel_                                    Dcel;
  typedef _Halfedge_const_iterator<Traits_, Dcel>  Self;
  typedef typename Dcel::Halfedge                  Halfedge;
  typedef typename Dcel::Halfedge_const_iterator   Halfedge_const_iter;

  Halfedge_const_iter     hit;           // The stored iterator.

  /*! Constructor from a DCEL halfedge iterator. */
  _Halfedge_const_iterator (const Halfedge_const_iter& it) :
    hit (it)
  {}

public:

  typedef typename Halfedge_const_iter::size_type       Size;
  typedef typename Halfedge_const_iter::size_type       size_type;
  typedef typename Halfedge_const_iter::difference_type difference_type;
  typedef typename Halfedge_const_iter::difference_type Difference;
  typedef std::bidirectional_iterator_tag               iterator_category;

  /*! Default constructor. */
  _Halfedge_const_iterator ()
  {}

  /*! Constructor from a non-const iterator. */
  _Halfedge_const_iterator (const _Halfedge_iterator<Traits_, Dcel>& iter) :
    hit (iter.hit)
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (hit == iter.hit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (hit != iter.hit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++hit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++hit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --hit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --hit;
    return (temp);
  }

  /*! Get a handle to the current halfedge. */
  _Halfedge_const_handle<Traits_, Dcel> operator* () const
  {
    _Halfedge_const_handle<Traits_, Dcel>     hh (&(*hit));
    return (hh);  
  }

  /*! Get a pointer to the current halfedge. */
  const Halfedge* operator-> () const
  {
    return (&(*hit));
  }
};

/*! \class
 * Iterator for the DCEL edges.
 */
template <class Traits_, class Dcel_>
class _Edge_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Edge_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Edge_iterator<Traits_, Dcel>          Self;
  typedef typename Dcel::Halfedge                Halfedge;
  typedef typename Dcel::Edge_iterator           Edge_iter;

  Edge_iter     eit;           // The stored iterator.

  /*! Constructor from a DCEL edge iterator. */
  _Edge_iterator (const Edge_iter& it) :
    eit (it)
  {}

public:

  typedef typename Edge_iter::size_type           Size;
  typedef typename Edge_iter::size_type           size_type;
  typedef typename Edge_iter::difference_type     difference_type;
  typedef typename Edge_iter::difference_type     Difference;
  typedef std::bidirectional_iterator_tag         iterator_category;

  /*! Default constructor. */
  _Edge_iterator ()
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (eit == iter.eit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (eit != iter.eit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++eit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++eit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --eit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --eit;
    return (temp);
  }

  /*! Get a handle to the current halfedge. */
  _Halfedge_handle<Traits_, Dcel> operator* () const
  {
    _Halfedge_handle<Traits_, Dcel>     hh (&(*eit));
    return (hh);
  }

  /*! Get a pointer to the current halfedge. */
  Halfedge* operator-> () const
  {
    return (&(*eit));
  }
};

/*! \class
 * Const iterator for the DCEL edges.
 */
template <class Traits_, class Dcel_>
class _Edge_const_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Edge_const_iterator<Traits_, Dcel>    Self;
  typedef typename Dcel::Halfedge                Halfedge;
  typedef typename Dcel::Edge_const_iterator     Edge_const_iter;

  Edge_const_iter     eit;           // The stored iterator.

  /*! Constructor from a DCEL edge iterator. */
  _Edge_const_iterator (const Edge_const_iter& it) :
    eit (it)
  {}

public:

  typedef typename Edge_const_iter::size_type         Size;
  typedef typename Edge_const_iter::size_type         size_type;
  typedef typename Edge_const_iter::difference_type   difference_type;
  typedef typename Edge_const_iter::difference_type   Difference;
  typedef std::bidirectional_iterator_tag             iterator_category;

  /*! Default constructor. */
  _Edge_const_iterator ()
  {}

  /*! Constructor from a non-const iterator. */
  _Edge_const_iterator (const _Edge_iterator<Traits_, Dcel>& iter) :
    eit (iter.eit)
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (eit == iter.eit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (eit != iter.eit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++eit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++eit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --eit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --eit;
    return (temp);
  }

  /*! Get a handle to the current halfedge. */
  _Halfedge_const_handle<Traits_, Dcel> operator* () const
  {
    _Halfedge_const_handle<Traits_, Dcel>     hh (&(*eit));
    return (hh);
  }

  /*! Get a pointer to the current halfedge. */
  const Halfedge* operator-> () const
  {
    return (&(*eit));
  }
};

/*! \class
 * Iterator for the DCEL faces.
 */
template <class Traits_, class Dcel_>
class _Face_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

  friend class _Face_const_iterator<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Face_iterator<Traits_, Dcel>          Self;
  typedef typename Dcel::Face                    Face;
  typedef typename Dcel::Face_iterator           Face_iter;

  Face_iter     fit;           // The stored iterator.

  /*! Constructor from a DCEL face iterator. */
  _Face_iterator (const Face_iter& it) :
    fit (it)
  {}

public:

  typedef typename Face_iter::size_type           Size;
  typedef typename Face_iter::size_type           size_type;
  typedef typename Face_iter::difference_type     difference_type;
  typedef typename Face_iter::difference_type     Difference;
  typedef std::bidirectional_iterator_tag         iterator_category;

  /*! Default constructor. */
  _Face_iterator ()
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (fit == iter.fit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (fit != iter.fit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++fit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++fit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --fit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --fit;
    return (temp);
  }

  /*! Get a handle to the current face. */
  _Face_handle<Traits_, Dcel> operator* () const
  {
    _Face_handle<Traits_, Dcel>     fh (&(*fit));
    return (fh);  
  }

  /*! Get a pointer to the current face. */
  Face* operator-> () const
  {
    return (&(*fit));
  }
};

/*! \class
 * Const iterator for the DCEL faces.
 */
template <class Traits_, class Dcel_>
class _Face_const_iterator
{
  friend class Arrangement_2<Traits_, Dcel_>;

protected:

  typedef Dcel_                                  Dcel;
  typedef _Face_const_iterator<Traits_, Dcel>    Self;
  typedef typename Dcel::Face                    Face;
  typedef typename Dcel::Face_const_iterator     Face_const_iter;

  Face_const_iter     fit;           // The stored iterator.

  /*! Constructor from a DCEL face iterator. */
  _Face_const_iterator (const Face_const_iter& it) :
    fit (it)
  {}

public:

  typedef typename Face_const_iter::size_type         Size;
  typedef typename Face_const_iter::size_type         size_type;
  typedef typename Face_const_iter::difference_type   difference_type;
  typedef typename Face_const_iter::difference_type   Difference;
  typedef std::bidirectional_iterator_tag             iterator_category;

  /*! Default constructor. */
  _Face_const_iterator ()
  {}

  /*! Constructor from a non-const iterator. */
  _Face_const_iterator (const _Face_iterator<Traits_, Dcel>& iter) :
    fit (iter.fit)
  {}

  /*! Equality operator. */
  bool operator== (const Self& iter) const
  {
    return (fit == iter.fit);
  }

  /*! Inequality operator. */
  bool operator!= (const Self& iter) const
  {
    return (fit != iter.fit);
  }

  /*! Increment operator (prefix). */
  Self& operator++ ()
  {
    ++fit;
    return (*this);
  }

  /*! Increment operator (postfix). */
  Self operator++ (int )
  {
    Self    temp = *this;
    ++fit;
    return (temp);
  }

  /*! Decrement operator (prefix). */
  Self& operator-- ()
  {
    --fit;
    return (*this);
  }

  /*! Decrement operator (postfix). */
  Self operator-- (int )
  {
    Self    temp = *this;
    --fit;
    return (temp);
  }

  /*! Get a handle to the current vertex. */
  _Face_const_handle<Traits_, Dcel> operator* () const
  {
    _Face_const_handle<Traits_, Dcel>     fh (&(*fit));
    return (fh);
  }

  /*! Get a pointer to the current face. */
  const Face* operator-> () const
  {
    return (&(*fit));
  }
};

CGAL_END_NAMESPACE

#endif
