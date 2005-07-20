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
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_LANDMARKS_POINT_LOCATION_H
#define CGAL_ARR_LANDMARKS_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_landmarks_point_location<Arrangement> template.
 */

//#define CGAL_LM_DEBUG
//#define LANDMARKS_CLOCK_DEBUG
//#define TRAITS_CLOCK_DEBUG
//#define RATIONAL_CIRCLES

#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/Arr_point_location/Arr_lm_vertices_generator.h>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A class that answers point-location and vertical ray-shooting queries
 * on a planar arrangement using the Landmarks algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_, 
	  class Arr_landmarks_generator_ = 
	                      Arr_landmarks_vertices_generator<Arrangement_> >
class Arr_landmarks_point_location 
{
public:

  typedef Arrangement_                         Arrangement_2;
  typedef typename Arrangement_2::Traits_2     Traits_2;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef typename Arrangement_2::Vertex_const_iterator	
	                                    Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator	
	                                    Halfedge_const_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator 
	                                    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
	                                    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Holes_const_iterator	
	                                    Holes_const_iterator;
  typedef typename Arrangement_2::Isolated_vertices_const_iterator
                                      Isolated_vertices_const_iterator;

  typedef typename Traits_2::Point_2		            Point_2;
  typedef typename Traits_2::X_monotone_curve_2	    X_monotone_curve_2;

  typedef Arr_landmarks_generator_		              Arr_landmarks_generator;

  typedef std::list<Halfedge_const_handle>          Edge_list;
  typedef typename Edge_list::iterator		          Std_edge_iterator;

protected:

  typedef Arr_traits_basic_wrapper_2<Traits_2>      Traits_wrapper_2;

  // Data members:
  const Arrangement_2     *p_arr;      // The associated arrangement.
  const Traits_wrapper_2  *traits;     // Its associated traits object.
  Arr_landmarks_generator *lm_gen;     // The associated generator of landmarks
  mutable Edge_list	   m_flipped_edges;// The list of edges that were
                                       // flipped during this point location
  mutable Halfedge_const_handle *m_start_edge; //edge the landmark is on
  bool	delete_generator;	             // Indicates whether the generator was
                                       // locally allocated.

public:

  /*! Default constructor. */
  Arr_landmarks_point_location () : 
    p_arr (NULL),
    traits (NULL),
    lm_gen(NULL), 
    m_start_edge(NULL),
    delete_generator (false)
  {}

  /*! Constructor given an arrangement only. */
  Arr_landmarks_point_location (const Arrangement_2& arr) :
    p_arr (&arr), 
    m_start_edge(NULL)
  {
    traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
    lm_gen = new Arr_landmarks_generator(arr);
    delete_generator = true;
  }

  /*! Constructor given an arrangement, and landmarks generator. */
  Arr_landmarks_point_location (const Arrangement_2& arr, 
				Arr_landmarks_generator *gen) :
    p_arr (&arr), 
    lm_gen (gen), 
    m_start_edge(NULL),
    delete_generator (false)
  {
    traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
  }

  /*! Destructor. */
  ~Arr_landmarks_point_location () 
  {
    if (delete_generator) 
      delete lm_gen;
  }
        
  /*! Attach an arrangement object. */
  void init (const Arrangement_2& arr, Arr_landmarks_generator *gen = NULL) 
  {
    p_arr = &arr;

    if (gen != NULL)
    {
      lm_gen = gen;
      delete_generator = false;
    }
    else
    {
      lm_gen = new Arr_landmarks_generator(arr);
      delete_generator = true;
    }

    traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
  }
  
  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object locate (const Point_2& p) const;

protected:

  /*!
   * Walks from the given vertex to the query point.
   * \param vh The given vertex handle.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object _walk_from_vertex(Vertex_const_handle vh, 
			   const Point_2 & p) const;

  /*!
   * Walks from a point on a given halfedge to the query point.
   * \param eh The given halfedge handle.
   * \param np The point that the walk starts from.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object _walk_from_edge(Halfedge_const_handle eh, 
			 const Point_2 & p, 
			 const Point_2 & np) const;
  
  /*!
   * Walks from a point in a face to the query point.
   * \param eh A halfedge handle that points to the face.
   * \param np The point that the walk starts from.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object _walk_from_face(Face_const_handle face, 
			 const Point_2 & p, 
			 const Point_2 & np) const;


  //IXX: this is the new find face, and not the old one.
  /*!
   * Walks from a point in a face to the query point.
   * \param eh A halfedge handle that points to the face.
   * \param np The point that the walk starts from.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object _find_face (const Point_2 & p, 
		     Vertex_const_handle vh,
		     bool & new_vertex) const;

  bool _is_point_in_face (const Point_2 & p, 
			  const Ccb_halfedge_const_circulator & face, 
			  bool & found_edge,             
			  bool & found_vertex,         
			  Halfedge_const_handle  & out_edge) const;

  bool _find_edge_to_flip (const Point_2 & p,            
			   const Point_2 & np,     
			   const Ccb_halfedge_const_circulator & face,
			   Halfedge_const_handle  & out_edge) const ;
  
  bool _check_approximate_intersection (const X_monotone_curve_2 & seg,
					const X_monotone_curve_2 & cv, 
					bool & intersect) const ;
  
};

CGAL_END_NAMESPACE

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_landmarks_pl_functions.h>

#endif
