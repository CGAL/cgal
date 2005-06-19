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
//                 (based on old version by: Iddo Hanniel,
//                                           Eyal Flato,
//                                           Oren Nechushtan,
//                                           Ester Ezra,
//                                           Shai Hirsch,
//                                           and Eugene Lipovetsky)
#ifndef CGAL_ARRANGEMENT_2_H
#define CGAL_ARRANGEMENT_2_H

/*! \file
 * The header file for the Arrangement_2<Traits,Dcel> class.
 */

#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arrangement_2_iterators.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/In_place_list.h>
#include <map>

CGAL_BEGIN_NAMESPACE

/*! \class
 * The arrangement class, representing planar subdivisions induced by
 * a set of arbitrary planar curves. 
 * The Traits parameter corresponds to a traits class that defines the
 * Point_2 and X_monotone_curve_2 types and implements the geometric
 * predicates and constructions for the family of curves it defines.
 * The Dcel parameter should be a model of the ArrDcel concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <class Traits_, 
          class Dcel_ = Arr_default_dcel<Traits_> > 
class Arrangement_2
{
public:

  typedef Traits_                               Traits_2;
  typedef Dcel_                                 Dcel;
  typedef Arrangement_2<Traits_2,Dcel>          Self;

  typedef typename Traits_2::Point_2            Point_2;
  typedef typename Traits_2::X_monotone_curve_2 X_monotone_curve_2;

  typedef typename Dcel::Size                   Size;

protected:

  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;

  typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;

  typedef typename Dcel::Vertex                 Vertex;
  typedef typename Dcel::Halfedge               Halfedge;
  typedef typename Dcel::Face                   Face;

public:

  // Handles, iterators and circulators.
  typedef _Vertex_handle<Traits_2, Dcel>              Vertex_handle;
  typedef _Vertex_const_handle<Traits_2, Dcel>        Vertex_const_handle;
  typedef _Vertex_iterator<Traits_2, Dcel>            Vertex_iterator;
  typedef _Vertex_const_iterator<Traits_2, Dcel>      Vertex_const_iterator;

  typedef _Halfedge_handle<Traits_2, Dcel>            Halfedge_handle;
  typedef _Halfedge_const_handle<Traits_2, Dcel>      Halfedge_const_handle;
  typedef _Halfedge_iterator<Traits_2, Dcel>          Halfedge_iterator;
  typedef _Halfedge_const_iterator<Traits_2, Dcel>    Halfedge_const_iterator;
  typedef _Edge_iterator<Traits_2, Dcel>              Edge_iterator;
  typedef _Edge_const_iterator<Traits_2, Dcel>        Edge_const_iterator;

  typedef _Face_handle<Traits_2, Dcel>                Face_handle;
  typedef _Face_const_handle<Traits_2, Dcel>          Face_const_handle;
  typedef _Face_iterator<Traits_2, Dcel>              Face_iterator;
  typedef _Face_const_iterator<Traits_2, Dcel>        Face_const_iterator;

  typedef _Halfedge_around_vertex_circulator<Traits_2, Dcel>
                                      Halfedge_around_vertex_circulator;
  typedef _Halfedge_around_vertex_const_circulator<Traits_2, Dcel>
                                      Halfedge_around_vertex_const_circulator;
  typedef _Ccb_halfedge_circulator<Traits_2, Dcel>        
                                      Ccb_halfedge_circulator;
  typedef _Ccb_halfedge_const_circulator<Traits_2, Dcel>  
                                      Ccb_halfedge_const_circulator;
  typedef _Holes_iterator<Traits_2, Dcel>             Holes_iterator;
  typedef _Holes_const_iterator<Traits_2, Dcel>       Holes_const_iterator;
  typedef _Isolated_vertices_iterator<Traits_2, Dcel>
                                      Isolated_vertices_iterator;
  typedef _Isolated_vertices_const_iterator<Traits_2, Dcel>
                                      Isolated_vertices_const_iterator;

protected:

  /*!
   * \class Representation of a point object stored in the points' container.
   */
  typedef typename Traits_2::Point_2                  Base_point_2;

  class Stored_point_2 : public Base_point_2,
                         public In_place_list_base<Stored_point_2>
  {
  public:

    /*! Default constructor. */
    Stored_point_2 ()
    {}

    /*! Constructor from a point. */
    Stored_point_2 (const Base_point_2& p) :
      Base_point_2 (p)
    {}
  };

  /*!
   * \class Representation of an x-monotone curve object stored in the 
   *        curves' container.
   */  
  class Stored_curve_2 : public X_monotone_curve_2,
                         public In_place_list_base<Stored_curve_2>
  {
  public:

    /*! Default constructor. */
    Stored_curve_2 ()
    {}

    /*! Constructor from an x-monotone curve. */
    Stored_curve_2 (const X_monotone_curve_2& cv) :
      X_monotone_curve_2 (cv)
    {}
  };

  typedef CGAL_ALLOCATOR(Stored_point_2)          Points_alloc;
  typedef CGAL_ALLOCATOR(Stored_curve_2)          Curves_alloc;

  typedef In_place_list<Stored_point_2, false>    Points_container;
  typedef In_place_list<Stored_curve_2, false>    X_curves_container;

  typedef Arr_observer<Self>                      Observer;
  typedef std::list<Observer*>                    Observers_container;
  typedef typename Observers_container::iterator  Observers_iterator;
  typedef typename Observers_container::reverse_iterator  
                                                  Observers_rev_iterator;

  // Data members:
  Dcel                dcel;         // The DCEL representing the arrangement.
  Face               *un_face;      // The unbounded face of the DCEL.
  Points_container    points;       // Container for the points that
                                    // correspond to the vertices.
  Points_alloc        points_alloc; // Allocator for the points.
  X_curves_container  curves;       // Container for the x-monotone curves
                                    // that correspond to the edges.
  Curves_alloc        curves_alloc; // Allocator for the curves.
  Observers_container observers;    // Storing pointers to existing observers.
  Traits_wrapper_2   *traits;       // The traits wrapper.
  bool                own_traits;   // Inidicate whether we should evetually
                                    // free the traits object.

public:

  /// \name Constructors.
  //@{

  /*! Default constructor */
  Arrangement_2 ();

  /*! Copy constructor. */
  Arrangement_2 (const Self& arr);

  /*! Constructor given a traits object. */
  Arrangement_2 (const Traits_2 *tr);
  //@}

  /// \name Assignment functions.
  //@{

  /*! Assignment operator. */
  Self& operator= (const Self& arr);

  /*! Assign an arrangement. */
  void assign (const Self& arr);
  //@}

  /// \name Destruction functions..
  //@{

  /*! Destructor. */
  virtual ~Arrangement_2 ();

  /*! Clear the arrangement. */
  void clear();
  //@}

  /*! Access the traits object. */
  const Traits_2* get_traits () const
  {
    return (traits);
  }

  /// \name Access the arrangement dimensions.
  //@{

  /*! Check whether the arrangement is empty. */
  bool is_empty () const
  {
    return (halfedges_begin() == halfedges_end());
  }

  /*! Get the number of arrangement vertices. */
  Size number_of_vertices () const
  {
    return (dcel.size_of_vertices());
  }

  /*! Get the number of arrangement halfedges (then result is always even). */
  Size number_of_halfedges () const
  {
    return (dcel.size_of_halfedges());
  }

  /*! Get the number of arrangement edges. */
  Size number_of_edges () const
  {
    return (dcel.size_of_halfedges() / 2);
  }

  /*! Get the number of arrangement faces. */
  Size number_of_faces () const
  {
    return (dcel.size_of_faces());
  }
  //@}

  /// \name Traversal functions for the arrangement vertices.
  //@{

  /*! Get an iterator for the first vertex in the arrangement. */
  Vertex_iterator vertices_begin() 
  { 
    return (Vertex_iterator (dcel.vertices_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement vertices. */
  Vertex_iterator vertices_end()
  {
    return (Vertex_iterator(dcel.vertices_end())); 
  }

  /*! Get a const iterator for the first vertex in the arrangement. */
  Vertex_const_iterator vertices_begin() const
  { 
    return (Vertex_const_iterator (dcel.vertices_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement vertices. */
  Vertex_const_iterator vertices_end() const
  {
    return (Vertex_const_iterator (dcel.vertices_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement halfedges.
  //@{

  /*! Get an iterator for the first halfedge in the arrangement. */
  Halfedge_iterator halfedges_begin() 
  { 
    return (Halfedge_iterator (dcel.halfedges_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement halfedges. */
  Halfedge_iterator halfedges_end()
  {
    return (Halfedge_iterator(dcel.halfedges_end())); 
  }

  /*! Get a const iterator for the first halfedge in the arrangement. */
  Halfedge_const_iterator halfedges_begin() const
  { 
    return (Halfedge_const_iterator (dcel.halfedges_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement halfedges. */
  Halfedge_const_iterator halfedges_end() const
  {
    return (Halfedge_const_iterator (dcel.halfedges_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement edges.
  //@{

  /*! Get an iterator for the first edge in the arrangement. */
  Edge_iterator edges_begin() 
  { 
    return (Edge_iterator (dcel.edges_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement edges. */
  Edge_iterator edges_end()
  {
    return (Edge_iterator(dcel.edges_end())); 
  }

  /*! Get a const iterator for the first edge in the arrangement. */
  Edge_const_iterator edges_begin() const
  { 
    return (Edge_const_iterator (dcel.edges_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement edges. */
  Edge_const_iterator edges_end() const
  {
    return (Edge_const_iterator (dcel.edges_end())); 
  }
  //@}

  /// \name Traversal functions for the arrangement faces.
  //@{

  /*! Get the unbounded face (non-const version). */
  Face_handle unbounded_face ()
  {
    return (Face_handle (un_face));
  }

  /*! Get the unbounded face (const version). */
  Face_const_handle unbounded_face () const
  {
    return (Face_const_handle (un_face));
  }

  /*!
   * Get the face containing the given isolated vertex (non-const version).
   * \param The query vertex.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A handle to the face containing v.
   */
  Face_handle incident_face (Vertex_handle v)
  {
    CGAL_precondition (v.is_isolated());

    // Get the fictitious incident halfedge of the vertex.
    Halfedge  *p_he = v.p_v->halfedge();

    // If this halfedge is vald, return it incident face. Otherwise, return
    // the unbounded face.
    if (p_he != NULL)
      return (Face_handle (p_he->face()));
    else
      return (Face_handle (un_face));
  }

  /*!
   * Get the face containing the given isolated vertex (const version).
   * \param The query vertex.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A const handle to the face containing v.
   */
  Face_const_handle incident_face (Vertex_const_handle v) const
  {
    CGAL_precondition (v.is_isolated());

    // Get the fictitious incident halfedge of the vertex.
    const Halfedge  *p_he = v.p_v->halfedge();

    // If this halfedge is vald, return it incident face. Otherwise, return
    // the unbounded face.
    if (p_he != NULL)
      return (Face_const_handle (p_he->face()));
    else
      return (Face_const_handle (un_face));
  }

  /*! Get an iterator for the first face in the arrangement. */
  Face_iterator faces_begin() 
  { 
    return (Face_iterator (dcel.faces_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement faces. */
  Face_iterator faces_end()
  {
    return (Face_iterator(dcel.faces_end())); 
  }

  /*! Get a const iterator for the first face in the arrangement. */
  Face_const_iterator faces_begin() const
  { 
    return (Face_const_iterator (dcel.faces_begin())); 
  }
  
  /*! Get a past-the-end const iterator for the arrangement faces. */
  Face_const_iterator faces_end() const
  {
    return (Face_const_iterator (dcel.faces_end())); 
  }
  //@}

  /// \name Casting away constness for handle types.
  //@{
  Vertex_handle non_const_handle (Vertex_const_handle vh)
  {
    Vertex    *p_v = const_cast<Vertex*> (vh.p_v);
    return (Vertex_handle (p_v));
  }

  Halfedge_handle non_const_handle (Halfedge_const_handle hh)
  {
    Halfedge  *p_he = const_cast<Halfedge*> (hh.p_he);
    return (Halfedge_handle (p_he));
  }

  Face_handle non_const_handle (Face_const_handle fh)
  {
    Face      *p_f = const_cast<Face*> (fh.p_f);
    return (Face_handle (p_f));
  }
  //@}

  /// \name Vertex manipulation functions.
  //@{

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   * \pre p is geometrically equivalent to the current point
   *      associated with v.
   * \return A handle for a the modified vertex (same as v).
   */
  Vertex_handle modify_vertex (Vertex_handle v,
			       const Point_2& p);

  /*!
   * Insert an isolated vertex in the interior of a given face.
   * \param p The isolated point.
   * \param f The face into which we insert the new isolated vertex.
   * \return A handle for the isolated vertex that has been created.
   */
  Vertex_handle insert_isolated_vertex (const Point_2& p,
                                        Face_handle f);

  /*!
   * Remove an isolated vertex from the interior of a given face.
   * \param v The vertex to remove.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A handle for the face containing v.
   */
  Face_handle remove_isolated_vertex (Vertex_handle v);
  
  ///@}

  /// \name Specilaized insertion functions for x-monotone curves.
  //@{

  /*!
   * Insert an x-monotone curve into the arrangement as a new hole (inner
   * component) inside the given face.
   * \param cv The given x-monotone curve.
   * \param f The face into which we insert the new hole.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target equals the curve target.
   */
  Halfedge_handle insert_in_face_interior (const X_monotone_curve_2& cv, 
                                           Face_handle f);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \pre The left endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_left_vertex (const X_monotone_curve_2& cv, 
					   Vertex_handle v);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's left endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_left_vertex (const X_monotone_curve_2& cv,
					   Halfedge_handle prev);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \pre The right endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_right_vertex (const X_monotone_curve_2& cv, 
					    Vertex_handle v);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's right endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_right_vertex (const X_monotone_curve_2& cv,
					    Halfedge_handle prev);
  
  /*! 
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices.
   * \param cv The given x-monotone curve.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \pre v1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from v1 to v2.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv, 
                                      Vertex_handle v1, 
                                      Vertex_handle v2);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in one of the circular lists around a vertex.
   * \param cv The given x-monotone curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param v2 The second vertex.
   * \pre The target vertex of prev1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1 to v2.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv, 
                                      Halfedge_handle h1, 
                                      Vertex_handle v2);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in both circular lists around these two vertices.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \pre The target vertices of prev1 and prev2 are cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   */
  Halfedge_handle insert_at_vertices (const X_monotone_curve_2 & cv,
                                      Halfedge_handle prev1, 
                                      Halfedge_handle prev2);

  //@}

  /// \name Halfedge manipulation functions.
  //@{

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   * \pre cv is geometrically equivalent to the current curve
   *      associated with e.
   * \return A handle for a the modified halfedge (same as e).
   */
  Halfedge_handle modify_edge (Halfedge_handle e, 
                               const X_monotone_curve_2& cv);

  /*!
   * Split a given edge into two, and associate the given x-monotone
   * curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param cv1 The curve that should be associated with the first split edge.
   * \param cv2 The curve that should be associated with the second split edge.
   * \pre cv1's source and cv2's target equal the endpoints of the curve
   *      currently assoicated with e (respectively), and cv1's target equals
   *      cv2's target, and this is the split point (ot vice versa).
   * \return A handle for the halfedge whose source is the source of the the
   *         original halfedge e, and whose target is the split point.
   */
  Halfedge_handle split_edge (Halfedge_handle e,
                              const X_monotone_curve_2& cv1, 
                              const X_monotone_curve_2& cv2);

  /*!
   * Merge two edges to form a single edge, and associate the given x-monotone
   * curve with the merged edge.
   * \param e1 The first edge to merge (one of the pair of twin halfegdes).
   * \param e2 The second edge to merge (one of the pair of twin halfegdes).
   * \param cv The curve that should be associated with merged edge.
   * \return A handle for the merged halfedge.
   */
  Halfedge_handle merge_edge (Halfedge_handle e1, 
                              Halfedge_handle e2, 
                              const X_monotone_curve_2& cv);              

  /*!
   * Remove an edge from the arrangement.
   * \param e The edge to remove (one of the pair of twin halfegdes).
   * \return A handle for the remaining face.
   */
  Face_handle remove_edge (Halfedge_handle e);
  //@}

protected:

  /// \name Allocating and de-allocating points and curves.
  //@{

  /*! Allocate a new point. */
  Stored_point_2 *_new_point (const Point_2& p)
  {
    Stored_point_2   *p_sp = points_alloc.allocate (1);

    points_alloc.construct (p_sp, p);
    points.push_back (*p_sp);
    return (p_sp);
  }

  /*! De-allocate a point. */
  void _delete_point (Point_2& p)
  {
    Stored_point_2   *p_sp = static_cast<Stored_point_2*> (&p);

    points.erase (p_sp);
    points_alloc.destroy (p_sp);
    points_alloc.deallocate (p_sp, 1);
  }

  /*! Allocate a new curve. */
  Stored_curve_2 *_new_curve (const X_monotone_curve_2& cv)
  {
    Stored_curve_2   *p_scv = curves_alloc.allocate (1);

    curves_alloc.construct (p_scv, cv);
    curves.push_back (*p_scv);
    return (p_scv);
  }

  /*! De-allocate a curve. */
  void _delete_curve (X_monotone_curve_2& cv)
  {
    Stored_curve_2   *p_scv = static_cast<Stored_curve_2*> (&cv);

    curves.erase (p_scv);
    curves_alloc.destroy (p_scv);
    curves_alloc.deallocate (p_scv, 1);
  }
  //@}

  /// \name Converting handles to pointers (for the arrangement accessor).
  //@{

  /*! Convert a vertex handle to a pointer to a DCEL vertex. */
  Vertex* _vertex (Vertex_handle vh) const
  {
    return (vh.p_v);
  }

  /*! Convert a constant vertex handle to a pointer to a DCEL vertex. */
  const Vertex* _vertex (Vertex_const_handle vh) const
  {
    return (vh.p_v);
  }

  /*! Convert a halfedge handle to a pointer to a DCEL halfedge. */
  Halfedge* _halfedge (Halfedge_handle hh) const
  {
    return (hh.p_he);
  }

  /*! Convert a constant halfedge handle to a pointer to a DCEL halfedge. */
  const Halfedge* _halfedge (Halfedge_const_handle hh) const
  {
    return (hh.p_he);
  }

  /*! Convert a face handle to a pointer to a DCEL face. */
  Face* _face (Face_handle fh) const
  {
    return (fh.p_f);
  }

  /*! Convert a constant face handle to a pointer to a DCEL face. */
  const Face* _face (Face_const_handle fh) const
  {
    return (fh.p_f);
  }

  /*! Convert a holes iterator to a DCEL holes iterator. */
  typename Dcel::Face::Holes_iterator _holes_iterator (Holes_iterator it)
  {
    return (it.holes_it);
  }

  /*! Convert an isolated vertices iterator to a DCEL iterator. */
  typename Dcel::Face::Isolated_vertices_iterator
      _isolated_vertices_iterator (Isolated_vertices_iterator it)
  {
    return (it.iso_verts_it);
  }
  //@}

  /// \name Converting pointers to handles (for the arrangement accessor).
  //@{

  /*! Convert a pointer to a DCEL vertex to a vertex handle. */
  Vertex_handle _handle_for (Vertex *v)
  {
    return (Vertex_handle (v));
  }

  /*! Convert a pointer to a DCEL vertex to a constant vertex handle. */
  Vertex_const_handle _const_handle_for (const Vertex *v) const
  {
    return (Vertex_const_handle (v));
  }

  /*! Convert a pointer to a DCEL halfedge to a halfedge handle. */
  Halfedge_handle _handle_for (Halfedge *he)
  {
    return (Halfedge_handle (he));
  }

  /*! Convert a pointer to a DCEL halfedge to a constant halfedge handle. */
  Halfedge_const_handle _const_handle_for (const Halfedge *he) const
  {
    return (Halfedge_const_handle (he));
  }

  /*! Convert a pointer to a DCEL face to a face handle. */
  Face_handle _handle_for (Face *f)
  {
    return (Face_handle (f));
  }

  /*! Convert a pointer to a DCEL face to a constant face handle. */
  Face_const_handle _const_handle_for (const Face *f) const
  {
    return (Face_const_handle (f));
  }
  //@}

  /// \name Auxiliary (protected) functions.
  //@{

  /*!
   * Locate the place for the given curve around the given vertex.
   * \param v The given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \pre v is one of cv's endpoints.
   * \return A pointer to a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   *         A NULL return value indicates a precondition violation.
   */
  Halfedge* _locate_around_vertex (Vertex *v,
                                   const X_monotone_curve_2& cv) const;

  /*!
   * Compute the distance (in halfedges) between two halfedges.
   * \param e1 The source halfedge.
   * \param e2 The destination halfedge.
   * \return In case e1 and e2 belong to the same connected component, the 
   *         function returns number of boundary halfedges between the two 
   *         halfedges. Otherwise, it returns (-1).
   */
  int _halfedge_distance (const Halfedge *e1, const Halfedge *e2) const;

  /*!
   * Determine whether a given query halfedge lies in the interior of a new
   * face we are about to create, by connecting it with another halfedge
   * using a given x-monotone curve.
   * \param prev1 The query halfedge.
   * \param prev2 The other halfedge we are about to connect with prev1.
   * \param cv The x-monotone curve we use to connect prev1 and prev2.
   * \pre prev1 and prev2 belong to the same connected component, and by
   *      connecting them using cv we form a new face.
   * \return (true) if prev1 lies in the interior of the face we are about
   *         to create, (false) otherwise - in which case prev2 must lie
   *         inside this new face.
   */
  bool _is_inside_new_face (const Halfedge *prev1,
                            const Halfedge *prev2,
                            const X_monotone_curve_2& cv) const;

  /*!
   * Determine whether a given point lies within the region bounded by
   * a boundary of a connected component.
   * \param p The query point.
   * \param he A halfedge on the boundary of the connected component.
   * \return (true) if the point lies within region, (false) otherwise.
   */
  bool _point_is_in (const Point_2& p, 
                     const Halfedge* he) const;

  /*!
   * Move a given hole from one face to another.
   * \param from_face The face currently containing the hole.
   * \param to_face The face into which we should move the hole.
   * \param hole A DCEL holes iterator pointing at the hole.
   */
  void _move_hole (Face *from_face,
                   Face *to_face,
                   typename Dcel::Face::Holes_iterator hole);

  /*!
   * Move a given isolated vertex from one face to another.
   * \param from_face The face currently containing the isolated vertex.
   * \param to_face The face into which we should move the isolated vertex.
   * \param vit A DCEL isolated vertices iterator pointing at the vertex.
   */
  void _move_isolated_vertex
                        (Face *from_face,
                         Face *to_face,
                         typename Dcel::Face::Isolated_vertices_iterator vit);

  /*!
   * Check whether the given halfedge lies on the outer boundary of the given
   * face.
   * \param f The given face.
   * \param e The given halfedge.
   * \return A pointer to a halfedge on the outer boundary of f in case e lies
   *         on this outer boundary, or NULL if it does not.
   */
  Halfedge* _is_on_outer_boundary (Face *f, Halfedge *e) const;

  /*!
   * Check whether the given halfedge lies on the inner boundary of the given
   * face.
   * \param f The given face.
   * \param e The given halfedge.
   * \return A pointer to a halfedge on the inner boundary of f in case e lies
   *         on this outer boundary, or NULL if it does not.
   */
  Halfedge* _is_on_inner_boundary (Face *f, Halfedge *e) const;

  /*!
   * Find the hole represented by a given halfedge from the holes container
   * of a given face and earse this hole once it is found.
   * \param f The given face.
   * \param e The given halfedge.
   * \return Whether the hole was found and erased or not.
   */
  bool _find_and_erase_hole (Face *f, Halfedge* e);

  /*!
   * Find the vertex in the isolated vertices container of a given face and
   * earse this vertex once it is found.
   * \param f The given face.
   * \param v The isolated vertex.
   * \return Whether the vertex was found and erased or not.
   */
  bool _find_and_erase_isolated_vertex (Face *f, Vertex* v);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to free arrangement vertices (newly created vertices
   * or existing isolated vertices), so a new hole is formed in the face
   * that contains the two vertices.
   * \param cv The given x-monotone curve.
   * \param f The face containing the two end vertices.
   * \param v1 The free vertex that corresponds to the left endpoint of cv.
   * \param v2 The free vertex that corresponds to the right endpoint of cv.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, directed from v1 to v2.
   */
  Halfedge* _insert_in_face_interior (const X_monotone_curve_2& cv,
				      Face *f,
				      Vertex *v1, Vertex *v2);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex. The other
   * endpoint corrsponds to a free vertex (a newly created vertex or an
   * isolated vertex).
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \param v The free vertex that corresponds to the other endpoint.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, whose target is the vertex v.
   */
  Halfedge* _insert_from_vertex (const X_monotone_curve_2& cv,
				 Halfedge *prev,
				 Vertex *v);

  /*!
   * Insert an x-monotone curve into the arrangement, where the end vertices
   * are given by the target points of two given halfedges.
   * The two halfedges should be given such that in case a new face is formed,
   * it will be the incident face of the halfedge directed from the first
   * vertex to the second vertex.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \param new_face Output - whether a new face has been created.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   *         In case a new face has been created, it is given as the incident
   *         face of this halfedge.
   */
  Halfedge* _insert_at_vertices (const X_monotone_curve_2& cv,
                                 Halfedge *prev1, 
                                 Halfedge *prev2,
                                 bool& new_face);

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   */
  void _modify_vertex (Vertex *v, 
		       const Point_2& p);

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   */
  void _modify_edge (Halfedge *he, 
		     const X_monotone_curve_2& cv);

  /*!
   * Split a given edge into two at a given point, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param p The split point.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is p.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is p and its target equals e's target.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  Halfedge* _split_edge (Halfedge *e,
                         const Point_2& p,
                         const X_monotone_curve_2& cv1, 
                         const X_monotone_curve_2& cv2);

  /*!
   * Remove a pair of twin halfedges from the arrangement.
   * \param e One of the halfedges to be removed.
   * \pre In case the removal causes the creation of a new hole, e should 
   *      point at this hole.
   * \return A pointer to the remaining face.
   */
  Face *_remove_edge (Halfedge *e);

  //@}

private:

  /// \name Managing and notifying the arrangement observers.
  //@{

  /*!
   * Register a new observer (so it starts receiving notifications).
   * \param p_obs A pointer to the observer object.
   */
  void _register_observer (Observer *p_obs)
  {
    observers.push_back (p_obs);
  }

  /*!
   * Unregister a new observer (so it stops receiving notifications).
   * \param p_obs A pointer to the observer object.
   * \return Whether the observer was successfully unregistered.
   */
  bool _unregister_observer (Observer *p_obs)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
    {
      if ((*iter) == p_obs)
      {
        // Remove the p_ob pointer from the list of observers.
        observers.erase (iter);
        return (true);
      }
    }

    // If we reached here, the observer was not registered.
    return (false);
  }

  /* Notify the observers on global arrangement operations: */

  void _notify_before_assign (const Self& arr)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_assign (arr);
  }

  void _notify_after_assign ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_assign();
  }

  void _notify_before_clear ()
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_clear();
  }

  void _notify_after_clear (Face_handle u)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_clear (u);
  }

  void _notify_before_global_change ()
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_global_change();
  }

  void _notify_after_global_change ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_global_change();
  }

  /* Notify the observers on local changes in the arrangement: */
  
  void _notify_before_create_vertex (const Point_2& p)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_create_vertex (p);
  }

  void _notify_after_create_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_vertex (v);
  }

  void _notify_before_create_edge (const X_monotone_curve_2& c,
				   Vertex_handle v1, Vertex_handle v2)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_create_edge (c, v1, v2);
  }

  void _notify_after_create_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_edge (e);
  }

  void _notify_before_modify_vertex (Vertex_handle v,
				     const Point_2& p)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_vertex (v, p);
  }

  void _notify_after_modify_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_vertex (v);
  }

  void _notify_before_modify_edge (Halfedge_handle e,
                                   const X_monotone_curve_2& c)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_edge (e, c);
  }

  void _notify_after_modify_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_edge (e);
  }

  void _notify_before_split_edge (Halfedge_handle e,
                                  const X_monotone_curve_2& c1,
                                  const X_monotone_curve_2& c2)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_split_edge (e, c1, c2);
  }

  void _notify_after_split_edge (Halfedge_handle e1,
                                 Halfedge_handle e2)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_edge (e1, e2);
  }

  void _notify_before_split_face (Face_handle f,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_split_face (f, e);
  }

  void _notify_after_split_face (Face_handle f,
                                 Face_handle new_f,
                                 bool is_hole)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_face (f, new_f, is_hole);
  }

  void _notify_before_add_hole (Face_handle f,
                                Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_add_hole (f, e);
  }

  void _notify_after_add_hole (Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_hole (h);
  }

  void _notify_before_add_isolated_vertex (Face_handle f,
					   Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_add_isolated_vertex (f, v);
  }

  void _notify_after_add_isolated_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_isolated_vertex (v);
  }

  void _notify_before_merge_edge (Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  const X_monotone_curve_2& c)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_edge (e1, e2, c);
  }

  void _notify_after_merge_edge (Halfedge_handle e)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_edge (e);
  }

  void _notify_before_merge_face (Face_handle f1,
                                  Face_handle f2,
                                  Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_face (f1, f2, e);
  }

  void _notify_after_merge_face (Face_handle f)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_face (f);
  }

  void _notify_before_move_hole (Face_handle from_f,
                                 Face_handle to_f,
                                 Ccb_halfedge_circulator h)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_move_hole (from_f, to_f, h);
  }

  void _notify_after_move_hole (Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_hole (h);
  }

  void _notify_before_move_isolated_vertex (Face_handle from_f,
					    Face_handle to_f,
					    Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_move_isolated_vertex (from_f, to_f, v);
  }

  void _notify_after_move_isolated_vertex (Vertex_handle v)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_isolated_vertex (v);
  }

  void _notify_before_remove_vertex (Vertex_handle v)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_vertex (v);
  }

  void _notify_after_remove_vertex ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_vertex ();
  }

  void _notify_before_remove_edge (Halfedge_handle e)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_edge (e);
  }

  void _notify_after_remove_edge ()
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_edge ();
  }

  void _notify_before_remove_hole (Face_handle f,
				   Ccb_halfedge_circulator h)
  {
    Observers_iterator   iter;
    Observers_iterator   end = observers.end();

    for (iter = observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_hole (f, h);
  }

  void _notify_after_remove_hole (Face_handle f)
  {
    Observers_rev_iterator   iter;
    Observers_rev_iterator   end = observers.rend();

    for (iter = observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_hole (f);
  }
  //@}

};

//-----------------------------------------------------------------------------
// Declarations of the various global insertion (and deletion) functions.
//-----------------------------------------------------------------------------

/*!
 * Insert a curve into the arrangement (incremental insertion).
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param pl A point-location object associated with the arrangement.
 * \param cv The curve to be inserted.
 */
template <class Traits, class Dcel, class PointLocation>
void insert (Arrangement_2<Traits,Dcel>& arr,
	     const PointLocation& pl,
	     const typename Traits::Curve_2& c);

/*!
 * Insert a range of curves into the arrangement (aggregated insertion). 
 * The inserted curves may intersect one another and may also intersect the 
 * existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be Curve_2.
 */
template <class Traits, class Dcel, class InputIterator>
void insert (Arrangement_2<Traits,Dcel>& arr,
	     InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion).
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param pl A point-location object associated with the arrangement.
 * \param cv The x-monotone curve to be inserted.
 */
template <class Traits, class Dcel, class PointLocation>
void insert_x_monotone (Arrangement_2<Traits,Dcel>& arr,
			const PointLocation& pl,
			const typename Traits::X_monotone_curve_2& c);

/*!
 * Insert a range of x-monotone curves into the arrangement (aggregated
 * insertion). The inserted x-monotone curves may intersect one another and
 * may also intersect the existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 */
template <class Traits, class Dcel, class InputIterator>
void insert_x_monotone (Arrangement_2<Traits,Dcel>& arr,
			InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement, such that the curve
 * interior does not intersect with any existing edge or vertex in the
 * arragement (incremental insertion).
 * \param arr The arrangement.
 * \param pl A point-location object associated with the arrangement.
 * \param c The x-monotone curve to be inserted.
 * \pre The interior of c does not intersect any existing edge or vertex.
 * \return A handle for one of the new halfedges created by the insertion. 
 */
template <class Traits, class Dcel, class PointLocation>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
insert_non_intersecting (Arrangement_2<Traits,Dcel>& arr,
			 const PointLocation& pl,
			 const typename Traits::X_monotone_curve_2& c);

/*!
 * Insert a range of pairwise interior-disjoint x-monotone curves into
 * the arrangement, such that the curve interiors do not intersect with
 * any existing edge or vertex in the arragement (aggregated insertion).
 * \param arr The arrangement.
 * \param begin An iterator for the first x-monotone curve in the range.
 * \param end A past-the-end iterator for the x-monotone curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 *      The curves in the range are pairwise interior-disjoint, and their
 *      interiors do not intersect any existing edge or vertex.
 */
template <class Traits, class Dcel, class InputIterator>
void insert_non_intersecting (Arrangement_2<Traits,Dcel>& arr,
			      InputIterator begin, InputIterator end);

/*!
 * Remove an edge from the arrangement. In case it is possible to merge
 * the edges incident to the end-vertices of the removed edge after its
 * deletion, the function performs these merges as well.
 * \param arr The arrangement.
 * \param e The edge to remove (one of the pair of twin halfegdes).
 * \return A handle for the remaining face.
 */
template <class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
remove_edge (Arrangement_2<Traits,Dcel>& arr,
	     typename Arrangement_2<Traits,Dcel>::Halfedge_handle e);

CGAL_END_NAMESPACE

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_2_functions.h>
#include <CGAL/Arrangement_2/Arrangement_2_insert.h>

#endif

