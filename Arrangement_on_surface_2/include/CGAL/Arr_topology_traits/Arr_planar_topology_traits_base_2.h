// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_PLANAR_TOPOLOGY_TRAITS_BASE_2_H
#define CGAL_ARR_PLANAR_TOPOLOGY_TRAITS_BASE_2_H

/*! \file
 * Definition of the Arr_planar_topology_traits_base_2<GeomTraits> class.
 */

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_basic_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_insertion_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_insertion_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_overlay_subcurve.h>
#include <CGAL/Sweep_line_2/Arr_overlay_traits_2.h>
#include <CGAL/Sweep_line_2/Arr_overlay_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_batched_pl_sl_visitor.h>
#include <CGAL/Sweep_line_2/Arr_vert_decomp_sl_visitor.h>
#include <CGAL/Arr_point_location/Arr_batched_point_location_traits_2.h>

namespace CGAL {

/*! \class Arr_planar_topology_traits_base_2
 * A base topology-traits class that encapsulates the embedding of 2D
 * arrangements of bounded or unbounded curves on the plane.
 */
template <class GeomTraits_,
          class Dcel_ = Arr_default_dcel<GeomTraits_> >
class Arr_planar_topology_traits_base_2
{
public:

  ///! \name The geometry-traits types.
  //@{
  typedef GeomTraits_                                     Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
  //@}

  ///! \name The DCEL types.
  //@{
  typedef Dcel_                                           Dcel;
  typedef typename Dcel::Size                             Size;
  typedef typename Dcel::Vertex                           Vertex;
  typedef typename Dcel::Halfedge                         Halfedge;
  typedef typename Dcel::Face                             Face;
  typedef typename Dcel::Outer_ccb                        Outer_ccb;
  typedef typename Dcel::Inner_ccb                        Inner_ccb;
  typedef typename Dcel::Isolated_vertex                  Isolated_vertex;
  //@}

  typedef Arr_planar_topology_traits_base_2<Geometry_traits_2,
                                            Dcel>         Self;

protected:

  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>    Traits_adaptor_2;

  // Data members:
  Dcel                m_dcel;       // The DCEL.

  const Traits_adaptor_2 * traits;       // The geometry-traits adaptor.
  bool                own_traits;   // Inidicate whether we should evetually
                                    // free the traits object.

  // Copy constructor and assignment operator - not supported.
  Arr_planar_topology_traits_base_2 (const Self& );
  Self& operator= (const Self& );

public:

  ///! \name Construction methods.
  //@{

  /*! Default constructor. */
  Arr_planar_topology_traits_base_2 () :
    own_traits (true)
  {
    traits = new Traits_adaptor_2;
  }

  /*! Constructor with a geometry-traits class. */
  Arr_planar_topology_traits_base_2 (const Geometry_traits_2 * geom_traits) :
    own_traits (false)
  {
    traits = static_cast<const Traits_adaptor_2*>(geom_traits);
  }

  /*! Assign the contents of another topology-traits class. */
  void assign (const Self& other);

  /*! Destructor. */
  virtual ~Arr_planar_topology_traits_base_2 ()
  {
    // Clear the DCEL.
    m_dcel.delete_all();

    if (own_traits)
      delete traits;    
  }
  //@}

  ///! \name Common topology-traits methods.
  //@{

  /*! Get the DCEL (const version). */
  const Dcel& dcel () const
  {
    return (m_dcel);
  }

  /*! Get the DCEL (non-const version). */
  Dcel& dcel ()
  {
    return (m_dcel);
  }

  /*!
   * Receive a notification on the creation of a new boundary vertex that
   * corresponds to the given curve end.
   * \param v The new boundary vertex.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param ps_x The boundary condition of the curve end in x.
   * \param ps_y The boundary condition of the curve end in y.
   */
  void notify_on_boundary_vertex_creation (Vertex *,
                                           const X_monotone_curve_2& ,
                                           Arr_curve_end,
                                           Arr_parameter_space /* ps_x */,
                                           Arr_parameter_space /* ps_y */)
  {
    // In the planar-topology traits this function should never be invoked:
    return;
  }

  /*!
   * Given two predecessor halfedges that belong to the same inner CCB of
   * a face, determine what happens when we insert an edge connecting the
   * target vertices of the two edges.
   * \param prev1 The first predecessor halfedge.
   * \param prev2 The second predecessor halfedge.
   * \param cv The curve to be inserted
   * \pre The two halfedges belong to the same inner CCB.
   * \return A pair indicating whether the insertion will cause the face
   *         to split (the first flag), and if so - whether the split face
   *         will form a hole in the original face.
   */
  std::pair<bool, bool>
  face_split_after_edge_insertion(const Halfedge *
                                    CGAL_precondition_code(prev1),
                                  const Halfedge *
                                    CGAL_precondition_code(prev2),
                                  const X_monotone_curve_2 & /* cv */) const
  {
    CGAL_precondition (prev1->is_on_inner_ccb());
    CGAL_precondition (prev2->is_on_inner_ccb());
    CGAL_precondition (prev1->inner_ccb() == prev2->inner_ccb());

    // In case of a planar topology, connecting two vertices on the same
    // inner CCB closes a new face that becomes a hole in the original face:
    return (std::make_pair (true, true));
  }

#if CGAL_NEW_FACE_SPLIT_STRATEGY
    /*!
   * Given two predecessor halfedges that belong to the same inner CCB of
   * a face, determine what happens when we insert an edge connecting the
   * target vertices of the two edges.
   * \param prev1 The first predecessor halfedge.
   * \param prev2 The second predecessor halfedge.
   * \param cv The curve to be inserted
   * \pre The two halfedges belong to the same inner CCB.
   * \return A pair indicating whether the insertion will cause the face
   *         to split (the first flag), and if so - whether the prev1 will be
   *         incident to the split face (second flag).
   */
  std::pair<bool, bool>
  face_update_upon_edge_insertion(const Halfedge *
                                    CGAL_precondition_code(prev1),
                                  const Halfedge *
                                    CGAL_precondition_code(prev2),
                                  const X_monotone_curve_2 & cv) const
  {
      // In case of a planar topology, connecting two vertices on the same
      // CCB closes a new face that becomes a hole in the original face:
      return (std::make_pair (true, true));
  }
#endif


  /*!
   * Determine whether the removal of the given edge will cause the creation
   * of a hole.
   * \param he The halfedge to be removed.
   * \pre Both he and its twin lie on an outer CCB of their incident faces.
   * \return Whether a new hole will be created.
   */
  bool hole_creation_after_edge_removal (const Halfedge *he) const
  {
    CGAL_precondition (! he->is_on_inner_ccb());
    CGAL_precondition (! he->opposite()->is_on_inner_ccb());

    // Check whether the halfedge and its twin belong to the same outer CCB
    // (and are therefore incident to the same face).
    if (he->outer_ccb() == he->opposite()->outer_ccb())
    {
      // In this case we cut an antenna, therefore we seperate a new hole
      // from the outer CCB of the face.
      return (true);
    }
    else
    {
      // The edge separates two faces. When we remove it, these two faces will
      // be merged, but no new hole will be created.
      return (false);
    }
  }

  /*!
   * Given two predecessor halfedges that will be used for inserting a
   * new halfedge pair (prev1 will be the predecessor of the halfedge he1,
   * and prev2 will be the predecessor of its twin he2), such that the
   * insertion will create a new perimetric face that forms a hole inside
   * an existing perimetric face, determine whether he1 will be incident to
   * this new face.
   * \param prev1 The first predecessor halfedge.
   * \param prev2 The second predecessor halfedge.
   * \param cv The x-monotone curve we use to connect prev1's target and
   *           prev2's target vertex.
   * \pre prev1 and prev2 belong to the same inner connected component.
   * \return true if he1 (and prev1) lies in the interior of the face we
   *         are about to create, false otherwise - in which case he2
   *         (and prev2) must be incident to this new face.
   */
  bool is_on_new_perimetric_face_boundary (const Halfedge *,
                                           const Halfedge *,
                                           const X_monotone_curve_2&) const
  {
    // We can never have perimetric faces in a planar topology:
    CGAL_error();
    return (false);
  }

  /*!
   * Determine whether the two halfedges, belonging to different outer CCBs,
   * belong to the outer boundary of the same face.
   * \param e1 The first halfedge.
   * \param e2 The second halfedge.
   * \return Whether the two halfedge belong to the outer boundary of the same
   *         face.
   */
  bool boundaries_of_same_face (const Halfedge *,
                                const Halfedge *) const
  {
    // This predicate is only used for case 3.3.2 of the insertion process,
    // therefore it should never be invoked in the planar case.
    CGAL_error();
    return (false);
  }


#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE || CGAL_NEW_FACE_SPLIT_STRATEGY
  /*!
   * Computes the sign of two halfedges approaching and leaving the
   * boundary
   * \param he1 The halfedge entering the boundary
   * \param he2 The halfedge leaving the boundary
   * \return the perimetricity of the subpath
   */
  CGAL::Sign _sign_of_subpath(const Halfedge* he1, const Halfedge* he2) 
    const {
    return CGAL::ZERO;
  }
    
  /*!
   * Computes the sign of a halfedge and a curve approaching and leaving the
   * boundary
   * \param he1 The halfedge entering the boundary
   * \param cv2 The curve leaving the boundary
   * \param end2 The end of the curve leaving the boundary
   * \return the perimetricity of the subpath
   */
  CGAL::Sign _sign_of_subpath(const Halfedge* he1, 
                              const bool target,
                              const X_monotone_curve_2& cv2,
                              const CGAL::Arr_curve_end& end2) const {
    return CGAL::ZERO;
  }
#endif

  /*!
   * Determine whether the given point lies in the interior of the given face.
   * \param f The face.
   * \param p The query point.
   * \param v The vertex associated with p (if exists).
   * \param f must not be fictitious, and v must not lie at infinity.
   * \return Whether p is contained in f's interior.
   */
  bool is_in_face (const Face *f, const Point_2& p, const Vertex *v) const;
  //@}

  /// \name Additional accessors, specialized for this topology-traits class.
  //@{

  /*! This function is used by the "walk" point-location strategy. */
  virtual const Face* initial_face () const = 0;
  //@}

  /// \name Additional predicates, specialized for this topology-traits class.
  //@{

  /*!
   * Compare the given vertex (which may lie at infinity) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the comparison of the x-coordinates of p and v.
   */
  virtual Comparison_result compare_x (const Point_2& p,
                                       const Vertex* v) const = 0;

  /*!
   * Compare the given vertex (which may lie at infinity) and the given point.
   * \param p The point.
   * \param v The vertex.
   * \return The result of the xy-lexicographic comparison of p and v.
   */
  virtual Comparison_result compare_xy (const Point_2& p,
                                        const Vertex* v) const = 0;

  /*!
   * Compare the relative y-position of the given point and the given edge
   * (which may be fictitious).
   * \param p The point.
   * \param he The edge (one of the pair of halfedges).
   * \pre p should lie in the x-range of the given edge.
   * \return The relative y-position of the point p and the edge.
   */
  virtual Comparison_result compare_y_at_x (const Point_2& p,
                                            const Halfedge* he) const = 0;
  //@}
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <class GeomTraits, class Dcel_>
void Arr_planar_topology_traits_base_2<GeomTraits, Dcel_>::assign
    (const Self& other)
{
  // Clear the current DCEL and duplicate the other DCEL.
  m_dcel.delete_all();
  m_dcel.assign (other.m_dcel);
  
  // Take care of the traits object.
  if (own_traits && traits != NULL)
    delete traits;
  
  if (other.own_traits)
    traits = new Traits_adaptor_2;
  else
    traits = other.traits;
  own_traits = other.own_traits;

  return;
}

//-----------------------------------------------------------------------------
// Determine whether the given vertex lies in the interior of the given face.
//
template <class GeomTraits, class Dcel_>
bool Arr_planar_topology_traits_base_2<GeomTraits, Dcel_>::
is_in_face(const Face *f, const Point_2& p, const Vertex *v) const
{
  CGAL_precondition (v == NULL || ! v->has_null_point());
  CGAL_precondition (v == NULL || 
                     traits->equal_2_object()(p, v->point()));

  // In case the face is unbounded and has no outer ccbs, this is the single
  // unbounded face of an arrangement of bounded curves. This face obviously
  // contains any point in its interior.
  if (f->is_unbounded() && f->number_of_outer_ccbs() == 0)
    return (true);

  // Keep a counter of the number of x-monotone curves that intersect an upward
  // vertical emanating from p (except for some degenerate cases that are
  // explained below).
  unsigned int       n_ray_intersections = 0;

  // Get a halfedge along the outer CCB of the given face, go over all curves
  // of the boundary component, and count those which are above p.
  // We begin by comparing p to the source vertex of the first halfedge.
  // Note that if p coincides with this vertex, p is obviously not in the
  // interior of the component.
  const Halfedge    *first = *(f->outer_ccbs_begin());

  
  // Some left ends of curves may not yet have the curve assigned, 
  // In case they are at TOP/BOTTOM they are not comparable to p 
  // We advance first until we find a comparable vertex to start with.
  // We always find such a vertex since on the ccb there is at least one 
  // interior vertex or a vertex at LEFT/RIGHT which are comparable. 
  while( first->vertex()->parameter_space_in_x()==ARR_INTERIOR 
      && first->has_null_curve() 
      && first->next()->has_null_curve()){ 
    first = first->next(); 
  }

  
  const Halfedge    *curr = first;
  Comparison_result  res_source;
  Comparison_result  res_target;
  Comparison_result  res_y_at_x;

  if (curr->opposite()->vertex() == v)
    return (false);

  res_source = compare_xy (p, curr->opposite()->vertex());

  do
  {
    // Compare p to the target vertex of the current halfedge.
    // If the vertex v associated with p (if v is given and is not NULL)
    // on the boundary of the component, p is obviously not in the interior
    // the component.
    if (curr->vertex() == v)
      return (false);
    
    // We jump over vertices at TOP/BOTTOM that do not yet have a curve 
    if(    curr->vertex()->parameter_space_in_x()==ARR_INTERIOR 
        && curr->has_null_curve() 
        && curr->next()->has_null_curve()){ 
      curr = curr->next(); 
      continue;
    }
    
    res_target = compare_xy (p, curr->vertex());

    // In case the current halfedge belongs to an "antenna", namely its
    // incident face is the same as its twin's, we can simply skip it
    // (in order not to count it twice).
    if (! curr->opposite()->is_on_inner_ccb() &&
        curr->outer_ccb()->face() == curr->opposite()->outer_ccb()->face())
    {
      curr = curr->next();
      res_source = res_target;
      continue;
    }

    // Check that if we shoot a "tilted" vertical ray from p upward
    // (by "tilted" we mean the angle it forms with the x-axis is
    //  PI/2 + epsilon, where epsilon is arbitrarily small), then we hit
    // the x-monotone curve associated with curr once.
    if (res_source != res_target)
    {
      res_y_at_x = compare_y_at_x (p, curr);

      if (res_y_at_x == SMALLER)
      {
        n_ray_intersections++;
      }        
      else if (res_y_at_x == EQUAL)
      {
        // In this case p lies on the current edge, so it is obviously not
        // contained in the interior of the component.
        return (false);
      }
    }

    // Proceed to the next halfedge along the component boundary.
    // Note that the source vertex of this halfedge is the current target.
    curr = curr->next();
    res_source = res_target;

  } while (curr != first);

  // The query point lies inside the connected components if and only if the
  // ray we shoot from it intersects the boundary an odd number of time.
  return ((n_ray_intersections % 2) != 0);
}

} //namespace CGAL

#endif
