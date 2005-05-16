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
  
protected:

  /*!
   * \class Representation of a point object stored in the points' container.
   */  
  class Stored_point_2 : public Point_2,
                         public In_place_list_base<Stored_point_2>
  {
  public:

    /*! Default constructor. */
    Stored_point_2 ()
    {}

    /*! Constructor from a point. */
    Stored_point_2 (const Point_2& p) :
      Point_2 (p)
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

  typedef In_place_list<Stored_point_2, false>    Points_container;
  typedef In_place_list<Stored_curve_2, false>    X_curves_container;

  typedef Arr_observer<Self>                      Observer;
  typedef std::list<Observer*>                    Observers_container;
  typedef typename Observers_container::iterator  Observers_iterater;

  // Data members:
  Dcel                dcel;     // The DCEL representing the arrangement.
  Face               *un_face;  // The unbounded face of the DCEL.
  Points_container    points;   // Storing the points that match the vertices.
  X_curves_container  curves;   // Storing the curves that match the edges.
  Observers_container observers;// Storing pointers to existing observers.
  Traits_wrapper_2   *traits;   // The traits wrapper.
  bool              own_traits; // Should we evetually free the traits object.

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

  /// \name Specilaized insertion functions.
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
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \pre v is one of cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, which is directed from left to right.
   */
  Halfedge_handle insert_from_vertex (const X_monotone_curve_2& cv, 
                                      Vertex_handle v);

  /*! 
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is one of cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_vertex (const X_monotone_curve_2& cv,
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
   * \return A handle for a new halfedge created by the split, whose target is
   *         the split point.
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
   * Find the holes represented by a given halfedge from the holes container
   * of a given face and earse this holes once it is found.
   * \param f The given face.
   * \param e The given halfedge.
   * \return Whether the hole was found and erased or not.
   */
  bool _find_and_erase_hole (Face *f, Halfedge* e);

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
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
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
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_assign (arr);
  }

  void _notify_after_assign ()
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_assign();
  }

  void _notify_before_clear ()
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_clear();
  }

  void _notify_after_clear (Face_handle u)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_clear (u);
  }

  void _notify_before_global_change ()
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_global_change();
  }

  void _notify_after_global_change ()
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_global_change();
  }

  /* Notify the observers on local changes in the arrangement: */

  void _notify_before_create_edge (const X_monotone_curve_2& c)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_create_edge (c);
  }

  void _notify_after_create_edge (Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_create_edge (e);
  }

  void _notify_before_modify_edge (Halfedge_handle e,
                                   const X_monotone_curve_2& c)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_modify_edge (e, c);
  }

  void _notify_after_modify_edge (Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_modify_edge (e);
  }

  void _notify_before_split_edge (Halfedge_handle e,
                                  const X_monotone_curve_2& c1,
                                  const X_monotone_curve_2& c2)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_split_edge (e, c1, c2);
  }

  void _notify_after_split_edge (Halfedge_handle e1,
                                 Halfedge_handle e2)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_split_edge (e1, e2);
  }

  void _notify_before_split_face (Face_handle f,
                                  Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_split_face (f, e);
  }

  void _notify_after_split_face (Face_handle f,
                                 Face_handle new_f,
                                 bool is_hole)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_split_face (f, new_f, is_hole);
  }

  void _notify_before_add_hole (Face_handle f,
                                Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_add_hole (f, e);
  }

  void _notify_after_add_hole (Ccb_halfedge_circulator h)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_add_hole (h);
  }

  void _notify_before_merge_edge (Halfedge_handle e1,
                                  Halfedge_handle e2,
                                  const X_monotone_curve_2& c)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_merge_edge (e1, e2, c);
  }

  void _notify_after_merge_edge (Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_merge_edge (e);
  }

  void _notify_before_merge_face (Face_handle f1,
                                  Face_handle f2,
                                  Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_merge_face (f1, f2, e);
  }

  void _notify_after_merge_face (Face_handle f)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_merge_face (f);
  }

  void _notify_before_move_hole (Face_handle from_f,
                                 Face_handle to_f,
                                 Ccb_halfedge_circulator h)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_move_hole (from_f, to_f, h);
  }

  void _notify_after_move_hole (Ccb_halfedge_circulator h)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->after_move_hole (h);
  }

  void _notify_before_remove_edge (Halfedge_handle e)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_remove_edge (e);
  }

  void _notify_before_remove_hole (Ccb_halfedge_circulator h)
  {
    Observers_iterater   iter;
    Observers_iterater   end = observers.end();

    for (iter = observers.begin(); iter != end; iter++)
      (*iter)->before_remove_hole (h);
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
template <class Arrangement, class PointLocation>
void arr_insert (Arrangement& arr, const PointLocation& pl,
                 const typename Arrangement::Traits_2::Curve_2& c);

/*!
 * Insert a range of curves into the arrangement (aggregated insertion). 
 * The inserted curves may intersect one another and may also intersect the 
 * existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be Curve_2.
 */
template <class Arrangement, class InputIterator>
void arr_insert (Arrangement& arr,
                 InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion).
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param pl A point-location object associated with the arrangement.
 * \param cv The x-monotone curve to be inserted.
 */
template <class Arrangement, class PointLocation>
void arr_insert_x_monotone (Arrangement& arr, const PointLocation& pl,
			    const typename Arrangement::X_monotone_curve_2& c);

/*!
 * Insert a range of x-monotone curves into the arrangement (aggregated
 * insertion). The inserted x-monotone curves may intersect one another and
 * may also intersect the existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 */
template <class Arrangement, class InputIterator>
void arr_insert_x_monotone (Arrangement& arr,
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
template <class Arrangement, class PointLocation>
typename Arrangement::Halfedge_handle
arr_insert_non_intersecting
                (Arrangement& arr, const PointLocation& pl,
                 const typename Arrangement::X_monotone_curve_2& c);

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
template <class Arrangement, class InputIterator>
void arr_insert_non_intersecting
                (Arrangement& arr,
                 InputIterator begin, InputIterator end);

/*!
 * Remove an edge from the arrangement. In case it is possible to merge
 * the edges incident to the end-vertices of the removed edge after its
 * deletion, the function performs these merges as well.
 * \param arr The arrangement.
 * \param e The edge to remove (one of the pair of twin halfegdes).
 * \return A handle for the remaining face.
 */
template <class Arrangement>
typename Arrangement::Face_handle
arr_remove_edge (Arrangement& arr,
                 typename Arrangement::Halfedge_handle e);

CGAL_END_NAMESPACE

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_2_functions.h>
#include <CGAL/Arrangement_2/Arrangement_2_insert.h>

#endif

