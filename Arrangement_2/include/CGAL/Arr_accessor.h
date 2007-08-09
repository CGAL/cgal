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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
#ifndef CGAL_ARR_ACCESSOR_H
#define CGAL_ARR_ACCESSOR_H

/*! \file
 * Definition of the Arr_accessor<Arrangement> class.
 */

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class that provides access to some of the internal arrangement operations.
 * Used mostly by the global insertion functions and by the sweep-line visitors
 * for utilizing topological and geometrical information available during the
 * algorithms they perform.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_>
class Arr_accessor
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef Arr_accessor<Arrangement_2>                   Self;

  typedef typename Arrangement_2::Size                  Size;
  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                     Ccb_halfedge_circulator;

  typedef typename Arrangement_2::_All_vertex_iterator  All_vertex_iterator;
  typedef typename Arrangement_2::_All_vertex_const_iterator
                                                   All_vertex_const_iterator;
  typedef typename Arrangement_2::_All_edge_iterator    All_edge_iterator;
  typedef typename Arrangement_2::_All_edge_const_iterator
                                                   All_edge_const_iterator;

private:

  typedef typename Arrangement_2::DVertex               DVertex;
  typedef typename Arrangement_2::DHalfedge             DHalfedge;
  typedef typename Arrangement_2::DFace                 DFace;
  typedef typename Arrangement_2::DHole                 DHole;
  typedef typename Arrangement_2::DIso_vert             DIso_vert;

private:

  Arrangement_2  *p_arr;           // The associated arrangement.

public:

  /*! Constructor with an associated arrangement. */
  Arr_accessor (Arrangement_2& arr) :
    p_arr (&arr)
  {}

  /* Get the arrangement. */
  Arrangement_2& arrangement()
  {
    return (*p_arr);
  }

  /* Get the arrangement (const version). */
  const Arrangement_2& arrangement() const
  {
    return (*p_arr);
  }


  /// \name Accessing the notification functions (for the global functions).
  //@{

  /*! Notify that a global operation is about to take place. */
  void notify_before_global_change ()
  {
    p_arr->_notify_before_global_change();
  }

  /*! Notify that a global operation was completed. */
  void notify_after_global_change ()
  {
    p_arr->_notify_after_global_change();
  }
  //@}

  /// \name Traversal functions.
  //@{

  /*!
   * Get one of the four fictitious vertices (non-const version).
   * \param inf_x, inf_y Specify the required fictitious vertex:
   *                     (-oo, -oo) is the bottom-left vertex.
   *                     (-oo, +oo) is the top-left vertex.
   *                     (+oo, -oo) is the bottom-right vertex.
   *                     (+oo, +oo) is the top-right vertex.
   * \pre inf_x and inf_y are not NO_BOUNDARY. 
   */
  Vertex_handle fictitious_vertex (Boundary_type inf_x,
                                   Boundary_type inf_y)
  {
    CGAL_precondition (inf_x != NO_BOUNDARY && inf_y != NO_BOUNDARY);

    DVertex      *v;

    if (inf_x == MINUS_INFINITY)
      v = (inf_y == MINUS_INFINITY) ? p_arr->v_bl : p_arr->v_tl;
    else
      v = (inf_y == MINUS_INFINITY) ? p_arr->v_br : p_arr->v_tr;

    return (p_arr->_handle_for (v));
  }

  /*!
   * Get one of the four fictitious vertices (const version).
   * \param inf_x, inf_y Specify the required fictitious vertex:
   *                     (-oo, -oo) is the bottom-left vertex.
   *                     (-oo, +oo) is the top-left vertex.
   *                     (+oo, -oo) is the bottom-right vertex.
   *                     (+oo, +oo) is the top-right vertex.
   * \pre inf_x and inf_y are not NO_BOUNDARY. 
   */
  Vertex_const_handle fictitious_vertex (Boundary_type inf_x,
                                         Boundary_type inf_y) const
  {
    CGAL_precondition (inf_x != NO_BOUNDARY && inf_y != NO_BOUNDARY);

    const DVertex      *v;

    if (inf_x == MINUS_INFINITY)
      v = (inf_y == MINUS_INFINITY) ? p_arr->v_bl : p_arr->v_tl;
    else
      v = (inf_y == MINUS_INFINITY) ? p_arr->v_br : p_arr->v_tr;

    return (p_arr->_const_handle_for (v));
  }

  /*!
   * Get the bottom left vertex (const version).
   */
  Vertex_const_handle bottom_left_fictitious_vertex() const
  {
    return (p_arr->_const_handle_for (p_arr->v_bl));
  }

  /*!
   * Get the top left vertex (const version).
   */
  Vertex_const_handle top_left_fictitious_vertex() const
  {
    return (p_arr->_const_handle_for (p_arr->v_tl));
  }

  /*!
   * Get the bottom right vertex (const version).
   */
  Vertex_const_handle bottom_right_fictitious_vertex() const
  {
    return (p_arr->_const_handle_for (p_arr->v_br));
  }

  /*!
   * Get the top right vertex (const version).
   */
  Vertex_const_handle top_right_fictitious_vertex() const
  {
    return (p_arr->_const_handle_for (p_arr->v_tr));
  }

  /*!
   * Get the bottom left vertex (non-const version).
   */
  Vertex_handle bottom_left_fictitious_vertex() 
  {
    return (p_arr->_handle_for (p_arr->v_bl));
  }

  /*!
   * Get the top left vertex (non-const version).
   */
  Vertex_handle top_left_fictitious_vertex() 
  {
    return (p_arr->_handle_for (p_arr->v_tl));
  }

  /*!
   * Get the bottom right vertex (non-const version).
   */
  Vertex_handle bottom_right_fictitious_vertex() 
  {
    return (p_arr->_handle_for (p_arr->v_br));
  }

  /*!
   * Get the top right vertex (non-const version).
   */
  Vertex_handle top_right_fictitious_vertex() 
  {
    return (p_arr->_handle_for (p_arr->v_tr));
  }

  /*! Get the fictitious face of the arrangement (non-const version). */
  Face_handle fictitious_face ()
  {
    return (p_arr->_handle_for (p_arr->un_face));
  }

  /*! Get the fictitious face of the arrangement (const version). */
  Face_const_handle fictitious_face () const
  {
    return (p_arr->_const_handle_for (p_arr->un_face));
  }

  /*! Get an iterator for the first vertex in the arrangement. */
  All_vertex_iterator all_vertices_begin() 
  { 
    typedef typename Arrangement_2::_Is_non_fictitious_vertex   
                                                     Is_non_fictitious_vertex;

    return (All_vertex_iterator (p_arr->dcel.vertices_begin(),
                                 p_arr->dcel.vertices_end(),
                                 Is_non_fictitious_vertex (p_arr->v_bl,
                                                           p_arr->v_tl,
                                                           p_arr->v_br,
                                                           p_arr->v_tr))); 
  }

  /*! Get a past-the-end iterator for the arrangement vertices. */
  All_vertex_iterator all_vertices_end()
  {
    typedef typename Arrangement_2::_Is_non_fictitious_vertex   
                                                     Is_non_fictitious_vertex;

    return (All_vertex_iterator (p_arr->dcel.vertices_end(),
                                 p_arr->dcel.vertices_end(),
                                 Is_non_fictitious_vertex (p_arr->v_bl,
                                                           p_arr->v_tl,
                                                           p_arr->v_br,
                                                           p_arr->v_tr))); 
  }

  /*! Get a const iterator for the first vertex in the arrangement. */
  All_vertex_const_iterator all_vertices_begin() const 
  { 
    typedef typename Arrangement_2::_Is_non_fictitious_vertex   
                                                     Is_non_fictitious_vertex;

    return (All_vertex_const_iterator
            (p_arr->dcel.vertices_begin(),
             p_arr->dcel.vertices_end(),
             Is_non_fictitious_vertex (p_arr->v_bl, p_arr->v_tl,
                                       p_arr->v_br, p_arr->v_tr))); 
  }

  /*! Get a past-the-end const iterator for the arrangement vertices. */
  All_vertex_const_iterator all_vertices_end() const
  {
    typedef typename Arrangement_2::_Is_non_fictitious_vertex   
                                                     Is_non_fictitious_vertex;

    return (All_vertex_const_iterator
            (p_arr->dcel.vertices_end(),
             p_arr->dcel.vertices_end(),
             Is_non_fictitious_vertex (p_arr->v_bl, p_arr->v_tl,
                                       p_arr->v_br, p_arr->v_tr))); 
  }

  /*! Get an iterator for the first edge in the arrangement. */
  All_edge_iterator all_edges_begin() 
  { 
    return (All_edge_iterator (p_arr->dcel.edges_begin())); 
  }

  /*! Get a past-the-end iterator for the arrangement edges. */
  All_edge_iterator all_edges_end()
  {
    return (All_edge_iterator (p_arr->dcel.edges_end())); 
  }

  /*! Get a const iterator for the first edge in the arrangement. */
  All_edge_const_iterator all_edges_begin() const
  { 
    return (All_edge_const_iterator (p_arr->dcel.edges_begin()));
  }
  
  /*! Get a past-the-end const iterator for the arrangement edges. */
  All_edge_const_iterator all_edges_end() const
  {
    return (All_edge_const_iterator (p_arr->dcel.edges_end()));
  }
  //@}

  /// \name Local operations and predicates for the arrangement.
  //@{

  /*!
   * Check whether the given infinite curve end lies on the given fictitious
   * edge.
   * \param cv The curve.
   * \param ind MIN_END if we refer to cv's minimal end;
   *            MAX_END if we refer to its maximal end.
   * \param he The fictitious edge.
   * \pre The relevent end of cv lies at infinity and he is a ficititious
   *      halfedge.
   * \return Whether the curve end lies on the fictitious edge.
   */
  bool is_on_ficitious_edge (const X_monotone_curve_2& cv, Curve_end ind,
                             Halfedge_const_handle he) const
  {
    typedef Arr_traits_basic_adaptor_2<typename Arrangement_2::Traits_2>
                                                             Traits_adaptor_2;

    // Check whether the relevant end of cv lies at infinity.
    const Traits_adaptor_2  *traits =
                   static_cast<const Traits_adaptor_2*> (p_arr->get_traits());
    const Boundary_type  inf_x = traits->boundary_in_x_2_object() (cv, ind);
    const Boundary_type  inf_y = traits->boundary_in_y_2_object() (cv, ind);
    bool                 eq_source, eq_target;

    CGAL_precondition (inf_x != NO_BOUNDARY || inf_y != NO_BOUNDARY);

    return (p_arr->_is_on_fictitious_edge (cv, ind,
                                           inf_x, inf_y,
                                           p_arr->_halfedge (he),
                                           eq_source, eq_target) &&
            ! eq_source && ! eq_target);
  }

  /*!
   * Locate an ficititious halfegde on the outer CCB of a given face which
   * contains the given curve end in its interior.
   * \param fh The face.
   * \param cv The curve.
   * \param ind MIN_END if we refer to cv's minimal end;
   *            MAX_END if we refer to its maximal end.
   * \pre The relevant end of cv must lie at infinity and f must be an
   *      unbounded face.
   * \return A handle to the halfedge along fh's outer CCB that contains the
   *         relevant end of cv in its interior.
   */
  Halfedge_handle locate_along_ccb (Face_handle fh,
                                    const X_monotone_curve_2& cv,
                                    Curve_end ind) const
  {
    typedef Arr_traits_basic_adaptor_2<typename Arrangement_2::Traits_2>
                                                             Traits_adaptor_2;

    // Check whether the relevant end of cv lies at infinity.
    const Traits_adaptor_2  *traits =
                   static_cast<const Traits_adaptor_2*> (p_arr->get_traits());
    const Boundary_type  inf_x = traits->boundary_in_x_2_object() (cv, ind);
    const Boundary_type  inf_y = traits->boundary_in_y_2_object() (cv, ind);

    CGAL_precondition (inf_x != NO_BOUNDARY || inf_y != NO_BOUNDARY);

    DHalfedge*  p_he = p_arr->_locate_along_ccb (p_arr->_face (fh),
                                                 cv, ind,
                                                 inf_x, inf_y);
    
    CGAL_assertion (p_he != NULL);
    return (p_arr->_handle_for (p_he));
  }

  /*!
   * Locate an unbounded end of a given curve along the bounding rectangle.
   * \param cv The curve.
   * \param ind MIN_END if we refer to cv's minimal end;
   *            MAX_END if we refer to its maximal end.
   * \pre The relevant end of cv lies at infinity.
   * \return An object along the bounding rectangle that contains the curve
   *         end. This object is either a Halfedge_const_hanlde (in the
   *         general case) or a Vertex_const_handle in case cv overlaps
   *         an existing edge.
   */
  CGAL::Object locate_unbounded_end (const X_monotone_curve_2& cv,
                                     Curve_end ind) const
  {
    typedef Arr_traits_basic_adaptor_2<typename Arrangement_2::Traits_2>
                                                             Traits_adaptor_2;

    // Check whether the relevant end of cv lies at infinity.
    const Traits_adaptor_2  *traits =
                   static_cast<const Traits_adaptor_2*> (p_arr->get_traits());
    const Boundary_type  inf_x = traits->boundary_in_x_2_object() (cv, ind);
    const Boundary_type  inf_y = traits->boundary_in_y_2_object() (cv, ind);

    CGAL_precondition (inf_x != NO_BOUNDARY || inf_y != NO_BOUNDARY);

    // Traverse the inner boundary of the fictitious face of the arrangement.
    const DHalfedge   *first = *(p_arr->un_face->holes_begin());
    const DHalfedge   *curr = first;
    bool               eq_source, eq_target;

    do
    {
      if (p_arr->_is_on_fictitious_edge (cv, ind,
                                         inf_x, inf_y,
                                         curr,
                                         eq_source, eq_target))
      {
        if (eq_source)
        {
          // cv's end coincides with the source vertex of the current
          // fictitious halfedge.
          Vertex_const_handle vh = p_arr->_const_handle_for 
                                                 (curr->opposite()->vertex());
          return (CGAL::make_object (vh));
        }
        else if (eq_target)
        {
          // cv's end coincides with the target vertex of the current
          // fictitious halfedge.
          Vertex_const_handle vh = p_arr->_const_handle_for (curr->vertex());
          return (CGAL::make_object (vh));
        }

        // The current ficitious edge contains cv's end in its interior.
        // Note we return curr's twin, whose incident face is a valid
        // unbounded face (whereas the incident face of curr is the fictitious
        // face).
        Halfedge_const_handle hh = p_arr->_const_handle_for (curr->opposite());
        return (CGAL::make_object (hh));
      }

      curr = curr->next();

    } while (curr != first);

    // We should never reach here.
    CGAL_assertion (false);
    return Object();
  }

  /*!
   * Locate the place for the given curve around the given vertex.
   * \param vh A handle for the arrangement vertex.
   * \param cv The given x-monotone curve.
   * \pre v is one of cv's endpoints.
   * \return A handle for a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   */
  Halfedge_handle locate_around_vertex (Vertex_handle vh,
                                        const X_monotone_curve_2& cv) const
  {
    typedef Arr_traits_basic_adaptor_2<typename Arrangement_2::Traits_2>
                                                             Traits_adaptor_2;

    const Traits_adaptor_2  *traits =
                   static_cast<const Traits_adaptor_2*> (p_arr->get_traits());
    Curve_end                ind = MIN_END;

    if (traits->boundary_in_x_2_object() (cv, MAX_END) == NO_BOUNDARY &&
        traits->boundary_in_y_2_object() (cv, MAX_END) == NO_BOUNDARY &&
        traits->equal_2_object()
        (vh->point(),
         p_arr->get_traits()->construct_max_vertex_2_object()(cv)))
    {
      ind = MAX_END;
    }

    DHalfedge*  he = p_arr->_locate_around_vertex (p_arr->_vertex (vh),
                                                   cv, ind);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Compute the distance (in halfedges) between two halfedges.
   * \param e1 A handle for the source halfedge.
   * \param e2 A handle for the destination halfedge.
   * \return In case e1 and e2 belong to the same connected component, the 
   *         function returns number of boundary halfedges between the two 
   *         halfedges. Otherwise, it returns (-1).
   */
  int halfedge_distance (Halfedge_const_handle e1,
                         Halfedge_const_handle e2) const
  {
    // If the two halfedges do not belong to the same component, return (-1).
    const DHalfedge     *he1 = p_arr->_halfedge (e1);
    const DHalfedge     *he2 = p_arr->_halfedge (e2);
    
    if (he1 == he2)
      return (0);

    const DHole         *hole1 = (he1->is_on_hole()) ? he1->hole() : NULL;
    const DHole         *hole2 = (he2->is_on_hole()) ? he2->hole() : NULL;
    const DFace         *f1 = (hole1 == NULL) ? he1->face() : hole1->face();
    const DFace         *f2 = (hole2 == NULL) ? he2->face() : hole2->face();

    if (f1 != f2 || hole1 != hole2)
      return (-1);

    // Compute the distance between the two halfedges.
    unsigned int         dist = p_arr->_halfedge_distance (he1, he2);
    return (static_cast<int> (dist));
  }

  /*!
   * Determine whether a given query halfedge lies in the interior of a new
   * face we are about to create, by connecting it with another halfedge
   * using a given x-monotone curve.
   * \param prev1 A handle for the query halfedge.
   * \param prev2 The other halfedge we are about to connect with prev1.
   * \param cv The x-monotone curve we use to connect prev1 and prev2.
   * \pre prev1 and prev2 belong to the same connected component, and by
   *      connecting them using cv we form a new face.
   * \return (true) if prev1 lies in the interior of the face we are about
   *         to create, (false) otherwise - in which case prev2 must lie
   *         inside this new face.
   */
  bool is_inside_new_face (Halfedge_handle prev1,
                           Halfedge_handle prev2,
                           const X_monotone_curve_2& cv) const
  {
    return (p_arr->_is_inside_new_face (p_arr->_halfedge (prev1),
                                        p_arr->_halfedge (prev2),
                                        cv));
  }

  /*!
   * Compare the x-coordinates of a given vertex (which may lie at infinity
   * if the traits class supports unbounded curves) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the comparison of the x-coordinates of p and v.
   */
  Comparison_result compare_x (const Point_2& p,
                               Vertex_const_handle v) const
  {
    return (p_arr->_compare_x (p, p_arr->_vertex (v)));
  }

  /*!
   * Compare the given vertex (which may lie at infinity if the traits class
   * supports unbounded curves) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the xy-lexicographic comparison of p and v.
   */
  Comparison_result compare_xy (const Point_2& p,
                                Vertex_const_handle v) const
  {
    return (p_arr->_compare_xy (p, p_arr->_vertex (v)));
  }

  /*!
   * Compare the relative y-position of the given point and the given edge
   * (which may be fictitious if the traits class supports unbounded curves).
   * \param p The point.
   * \param he The halfedge.
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  Comparison_result compare_y_at_x (const Point_2& p,
                                    Halfedge_const_handle he) const
  {
    return (p_arr->_compare_y_at_x (p, p_arr->_halfedge (he)));
  }

  /*!
   * Check if the given vertex represents one of the ends of a given curve.
   * \param v The vertex.
   * \param cv The curve.
   * \param ind MIN_END if we refer to cv's minimal end;
   *            MAX_END if we refer to its maximal end.
   * \return Whether v represents the left (or right) end of cv.
   */
  bool are_equal (Vertex_const_handle v,
                  const X_monotone_curve_2& cv, Curve_end ind) const
  {
    return (p_arr->_are_equal (p_arr->_vertex (v), cv, ind));
  }

  /*!
   * Determine whether a given point lies within the region bounded by
   * a boundary of a connected component.
   * \param p The query point.
   * \param he A handle for a halfedge on the boundary of the connected
   *           component.
   * \return (true) if the point lies within region, (false) otherwise.
   */
  bool point_is_in (const Point_2& p, 
                    Halfedge_const_handle he) const
  {
    return (p_arr->_point_is_in (p, NULL, p_arr->_halfedge (he)));
  }

  /*!
   * Check whether the given halfedge lies on the outer boundary of its
   * incident face.
   * \param he The given halfedge.
   * \return (true) in case he lies on the outer boundary of its incident face;
   *         (false) if he lies on a hole inside this face.
   */
  bool is_on_outer_boundary (Halfedge_const_handle he) const
  {
    const DHalfedge    *p_he = p_arr->_halfedge (he);

    return (! p_he->is_on_hole());
  }

  /*!
   * Check whether the given halfedge lies on the inner boundary of its
   * incident face.
   * \param he The given halfedge.
   * \return (true) in case he lies on a hole inside its incident face;
   *         (false) if he lies on the outer boundary of this face.
   */
  bool is_on_inner_boundary (Halfedge_const_handle he) const
  {
    const DHalfedge    *p_he = p_arr->_halfedge (he);

    return (p_he->is_on_hole());
  }

  /*!
   * Create a new vertex and associate it with the given point.
   * \param p The point.
   * \return A handle to the newly created vertex.
   */
  Vertex_handle create_vertex (const Point_2& p, Boundary_type, Boundary_type )
  {
    DVertex* v = p_arr->_create_vertex (p);
    
    CGAL_assertion (v != NULL);
    return (p_arr->_handle_for (v));    
  }

   Vertex_handle create_vertex (const Point_2& p )
   {
     DVertex* v = p_arr->_create_vertex (p);
     CGAL_assertion (v != NULL);
     return (p_arr->_handle_for (v));
   }
		    
  
  /*!
   * Create a new vertex at infinity.
   * \param inf_x MINUS_INFINITY if this vertex lies at x = -oo;
   *              PLUS_INFINITY if this vertex lies at x = +oo;
   *              NO_BOUNDARY if the vertex has a finite x-coordinate.
   * \param inf_y MINUS_INFINITY if this vertex lies at y = -oo;
   *              PLUS_INFINITY if this vertex lies at y = +oo;
   *              NO_BOUNDARY if the vertex has a finite y-coordinate.
   * \return A pointer to the newly created vertex.
   */
  Vertex_handle create_vertex_at_infinity (Boundary_type inf_x,
                                           Boundary_type inf_y)
  {
    DVertex   *v = p_arr->_create_vertex_at_infinity (inf_x, inf_y);

    CGAL_assertion (v != NULL);
    return (p_arr->_handle_for (v));
  }

  /*!
   * Insert an x-monotone curve into the arrangement, where the end vertices
   * are given by the target points of two given halfedges.
   * The two halfedges should be given such that in case a new face is formed,
   * it will be the incident face of the halfedge directed from the first
   * vertex to the second vertex.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \param res The comparsion result between the points associated with the
   *            target vertex of prev and the target vertex of prev2.
   * \param new_face Output - whether a new face has been created.
   * \param both_unbounded Indicates whether in case the face is split, both
   *                       resulting faces should be unbounded.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   *         In case a new face has been created, it is given as the incident
   *         face of this halfedge.
   */
  Halfedge_handle insert_at_vertices_ex (const X_monotone_curve_2& cv,
                                         Halfedge_handle prev1, 
                                         Halfedge_handle prev2,
                                         Comparison_result res,
                                         bool& new_face,
                                         bool both_unbounded = false)
  {
    DHalfedge*  he = p_arr->_insert_at_vertices (cv,
                                                 p_arr->_halfedge (prev1),
                                                 p_arr->_halfedge (prev2),
                                                 res,
                                                 new_face,
                                                 both_unbounded);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

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
   * \param res The comparsion result between the points associated with
   *            the target vertex of prev and the vertex v.
   * \return A handle to one of the halfedges corresponding to the inserted
   *         curve, whose target is the vertex v.
   */
  Halfedge_handle insert_from_vertex_ex (const X_monotone_curve_2& cv,
                                         Halfedge_handle prev,
                                         Vertex_handle v,
                                         Comparison_result res)
  {
    DHalfedge*  he = p_arr->_insert_from_vertex (cv,
                                                 p_arr->_halfedge (prev),
                                                 p_arr->_vertex (v),
                                                 res);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to free arrangement vertices (newly created vertices
   * or existing isolated vertices), so a new hole is formed in the face
   * that contains the two vertices.
   * \param cv The given x-monotone curve.
   * \param f The face containing the two end vertices.
   * \param v1 The free vertex that corresponds to the left endpoint of cv.
   * \param v2 The free vertex that corresponds to the right endpoint of cv.
   * \param res The comparsion result between the points associated with the
   *            vertices v1 and v2.
   * \return A handle to one of the halfedges corresponding to the inserted
   *         curve, directed from v1 to v2.
   */
  Halfedge_handle insert_in_face_interior_ex (const X_monotone_curve_2& cv,
                                              Face_handle f,
                                              Vertex_handle v1,
                                              Vertex_handle v2,
                                              Comparison_result res)
  {
    DHalfedge*  he = p_arr->_insert_in_face_interior (cv,
                                                      p_arr->_face (f),
                                                      p_arr->_vertex (v1),
                                                      p_arr->_vertex (v2),
                                                      res);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  
  }

  /*!
   * Insert the given vertex as an isolated vertex inside the given face.
   * \param f The face that should contain the isolated vertex.
   * \param v The isolated vertex.
   */
  void insert_isolated_vertex (Face_handle f, Vertex_handle v)
  {
    p_arr->_insert_isolated_vertex (p_arr->_face (f), p_arr->_vertex(v));
  }
  
  /*!
   * Relocate all holes and isolated vertices to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * In case insert_at_vertices_ex() was invoked and indicated that a new face
   * has been created, this function should be called with the halfedge
   * returned by insert_at_vertices_ex().
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void relocate_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_in_new_face (p_arr->_halfedge (new_he));
    return;
  }

  void relocate_isolated_vertices_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_isolated_vertices_in_new_face (p_arr->_halfedge(new_he));
    return;
  }

  void relocate_holes_in_new_face (Halfedge_handle new_he)
  {
    p_arr->_relocate_holes_in_new_face (p_arr->_halfedge(new_he));
    return;
  }

  /*!
   * Move a hole from one face to another.
   * \param from_face The source face.
   * \param to_face The destination face.
   * \param hole A CCB circulator that corresponds to the outer boundary
   *             of the hole to move.
   */
  void move_inner_ccb (Face_handle from_face, Face_handle to_face,
                  Ccb_halfedge_circulator hole)
  {
    DHalfedge        *he = p_arr->_halfedge (hole);

    p_arr->_move_hole (p_arr->_face (from_face),
                       p_arr->_face (to_face),
                       he->hole()->iterator());
  }
  
  /*!
   * Move an isolated vertex from one face to another.
   * \param from_face The source face.
   * \param to_face The destination face.
   * \param v The isolated vertex to move.
   */
  void move_isolated_vertex (Face_handle from_face, Face_handle to_face,
                             Vertex_handle v)
  {
    DVertex          *iso_v = p_arr->_vertex (v);
     
    p_arr->_move_isolated_vertex (p_arr->_face (from_face),
                                  p_arr->_face (to_face),
                                  iso_v->isolated_vertex()->iterator());
  }

  /*!
   * Remove an isolated vertex from its face.
   * \param v The isolated vertex to remove.
   */
  void remove_isolated_vertex_ex (Vertex_handle v)
  {
    CGAL_precondition (v->is_isolated());
    DVertex *iso_v = p_arr->_vertex (v);

    p_arr->_remove_isolated_vertex (iso_v);
    return;
  }

  /*!
   * Remove a vertex at infinity, causing its two incident fictitious edges
   * to merge.
   * \param v The vertex at infinity to remove.
   */
  void remove_vertex_at_infinity (Vertex_handle v)
  {
    CGAL_precondition (v->is_at_infinity());
    DVertex  *inf_v = p_arr->_vertex (v);

    p_arr->_remove_vertex_at_infinity (inf_v);
    return;
  }

  /*!
   * Modify the point associated with a given vertex. The point may be
   * geometrically different than the one currently associated with the vertex.
   * \param v The vertex to modify.
   * \param p The new point to associate with v.
   * \return A handle for the modified vertex (same as v).
   */
  Vertex_handle modify_vertex_ex (Vertex_handle v,
                                  const Point_2& p)
  {
    p_arr->_modify_vertex (p_arr->_vertex (v),
                           p);


    return (v);
  }
        
  /*!
   * Modify the x-monotone curve associated with a given edge. The curve may be
   * geometrically different than the one currently associated with the edge.
   * \param e The edge to modify.
   * \param cv The new x-monotone curve to associate with e.
   * \return A handle for the modified edge (same as e).
   */
  Halfedge_handle modify_edge_ex (Halfedge_handle e,
                                  const X_monotone_curve_2& cv)
  {
    p_arr->_modify_edge (p_arr->_halfedge (e),
                         cv);

    return (e);
  }
          
  /*!
   * Split a given edge into two at a given point, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param p The split point.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is p.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is p and its target equals e's target.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  Halfedge_handle split_edge_ex (Halfedge_handle e,
                                 const Point_2& p,
                                 const X_monotone_curve_2& cv1, 
                                 const X_monotone_curve_2& cv2)
  {
    DHalfedge*  he = p_arr->_split_edge (p_arr->_halfedge (e),
                                         p,
                                         cv1, cv2);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Split a given edge into two at the given vertex, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is v's point.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is v's point and its target equals e's target.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split vertex v.
   */
  Halfedge_handle split_edge_ex (Halfedge_handle e,
                                 Vertex_handle v,
                                 const X_monotone_curve_2& cv1,
                                 const X_monotone_curve_2& cv2)
  {
    DHalfedge*  he = p_arr->_split_edge (p_arr->_halfedge (e),
                                         p_arr->_vertex (v),
                                         cv1, cv2);

    CGAL_assertion (he != NULL);
    return (p_arr->_handle_for (he));
  }

  /*!
   * Split a given fictitious edge into two, forming a new vertex at
   * infinity.
   * \param he The edge to split (one of the pair of twin halfegdes).
   * \param inf_x Is the x-coordinate of the new vertex at infinity.
   * \param inf_y Is the y-coordinate of the new vertex at infinty.
   * \pre inf_x and inf_y cannot be both NO_BOUNDARY.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is the newly created vertex.
   */
  Halfedge_handle split_fictitious_edge (Halfedge_handle he,
                                         Boundary_type inf_x,
                                         Boundary_type inf_y)
  {
    DVertex    *p_he = p_arr->_split_fictitious_edge (p_arr->_halfedge (he),
                                                      inf_x, inf_y);
    
    return (p_arr->_handle_for (p_he));
  }

  /*!
   * Split a given fictitious edge into two at a given vertex at infinity.
   * \param he The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \return A handle for the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  Halfedge_handle split_fictitious_edge (Halfedge_handle he,
                                         Vertex_handle v)
  {
    DHalfedge*  p_he = p_arr->_split_fictitious_edge (p_arr->_halfedge (he),
                                                      p_arr->_vertex (v));

    return (p_arr->_handle_for (p_he));
  }

  /*!
   * Remove a pair of twin halfedges from the arrangement.
   * \param e A handle for one of the halfedges to be removed.
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \pre In case the removal causes the creation of a new hole, e should 
   *      point at this hole.
   * \return A handle for the remaining face.
   */
  Face_handle remove_edge_ex (Halfedge_handle e,
                              bool remove_source = true,
                              bool remove_target = true)
  {
    DFace*      f = p_arr->_remove_edge (p_arr->_halfedge (e),
                                         remove_source, remove_target);
    
    CGAL_assertion (f != NULL);
    return (p_arr->_handle_for (f));
  }


  bool are_on_same_inner_component(Halfedge_handle e1, Halfedge_handle e2)
  {
     DHalfedge  *p_prev1 = p_arr->_halfedge (e1);
     DHalfedge  *p_prev2 = p_arr->_halfedge (e2);
     DHole      *hole1 = (p_prev1->is_on_hole()) ? p_prev1->hole() : NULL;
     DHole      *hole2 = (p_prev2->is_on_hole()) ? p_prev2->hole() : NULL;

     return (hole1 == hole2 && hole1 != NULL);
  }
  //@}

  /// \name Functions used by the arrangement reader.
  //@{
  typedef DVertex                         Dcel_vertex;
  typedef DHalfedge                       Dcel_halfedge;
  typedef DFace                           Dcel_face;
  typedef DHole                           Dcel_hole;
  typedef DIso_vert                       Dcel_isolated_vertex;

  /*!
   * Clear the entire arrangment.
   */
  void clear_all ()
  {
    p_arr->clear();
    p_arr->dcel.delete_all();
    return;
  }

  /*!
   * Create a new vertex, associated with the given point.
   * \param p The point.
   * \return A pointer to the created DCEL vertex.
   */
  Dcel_vertex* new_vertex (const Point_2& p)
  {
    typename Arrangement_2::Stored_point_2  *p_pt = p_arr->_new_point (p);
    Dcel_vertex                             *new_v = p_arr->dcel.new_vertex();

    new_v->set_point (p_pt);
    return (new_v);
  }

  /*!
   * Create a new vertex at infinity.
   * \param inf_x The infinity type in x.
   * \param inf_y The infinity type in y.
   * \return A pointer to the created DCEL vertex.
   */
  Dcel_vertex* new_vertex_at_infinity (Boundary_type inf_x,
                                       Boundary_type inf_y)
  {
    Dcel_vertex                             *new_v = p_arr->dcel.new_vertex();

    new_v->set_boundary (inf_x, inf_y);
    return (new_v);
  }

  /*!
   * Create a new edge (halfedge pair), associated with the given curve.
   * \param cv The x-monotone curve.
   * \return A pointer to one of the created DCEL halfedge.
   */
  Dcel_halfedge* new_edge (const X_monotone_curve_2& cv)
  {
    typename Arrangement_2::Stored_curve_2  *p_cv = p_arr->_new_curve (cv);
    Dcel_halfedge                           *new_he = p_arr->dcel.new_edge();

    new_he->set_curve (p_cv);
    return (new_he);
  }

  /*!
   * Create a new fictitious edge (halfedge pair).
   * \return A pointer to one of the created DCEL halfedge.
   */
  Dcel_halfedge* new_fictitious_edge ()
  {
    Dcel_halfedge                           *new_he = p_arr->dcel.new_edge();

    new_he->set_curve (NULL);
    return (new_he);
  }

  /*!
   * Create a new face.
   * \return A pointer to the created DCEL face.
   */
  Dcel_face* new_face ()
  {
    return (p_arr->dcel.new_face());
  }

  /*!
   * Create a new hole.
   * \return A pointer to the created DCEL hole.
   */
  Dcel_hole* new_hole ()
  {
    return (p_arr->dcel.new_hole());
  }

  /*!
   * Create a new isolated vertex.
   * \return A pointer to the created DCEL isolated vertex.
   */
  Dcel_isolated_vertex* new_isolated_vertex ()
  {
    return (p_arr->dcel.new_isolated_vertex());
  }
  //@}
};

CGAL_END_NAMESPACE

#endif
