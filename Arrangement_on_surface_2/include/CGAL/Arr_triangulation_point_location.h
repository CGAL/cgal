// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_TRIANGULATION_POINT_LOCATION_H
#define CGAL_ARR_TRIANGULATION_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_triangulation_point_location<Arrangement> template.
 */

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

namespace CGAL {

/*! \class
 * A class that answers point-location and queries
 * on a planar arrangement using the triangulation algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation.
 */
template <class Arrangement_>
class Arr_triangulation_point_location : public Arr_observer <Arrangement_>
{
public:

  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Traits_2::Kernel                     Kernel;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
	typedef typename Arrangement_2::Vertex_handle		      Vertex_handle;
	typedef typename Arrangement_2::Halfedge_handle		    Halfedge_handle;
	typedef typename Arrangement_2::Face_handle			      Face_handle;

  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Hole_const_iterator  Hole_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator  
    Halfedge_const_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator 
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator 
    Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;

  typedef typename Traits_2::Point_2                    Point_2;
  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;

  typedef std::list<Halfedge_const_handle>              Edge_list;
  typedef typename Edge_list::iterator                  Std_edge_iterator;

    //----------------------------------------------------------
  // Triangulation Types
  //----------------------------------------------------------
  typedef Triangulation_vertex_base_with_info_2<Vertex_const_handle,Kernel> 
                                                                      Vbb;
  typedef Triangulation_hierarchy_vertex_base_2<Vbb>                  Vb;
  //typedef Triangulation_face_base_with_info_2<CGAL::Color,Kernel>    Fbt;
  typedef Constrained_triangulation_face_base_2<Kernel>               Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>                       TDS;
  typedef Exact_predicates_tag                                        Itag;
  //typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>    CDT;
  typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>     CDT_t;
  typedef Triangulation_hierarchy_2<CDT_t>                            CDTH;
  typedef Constrained_triangulation_plus_2<CDTH>                      CDT;

  typedef typename CDT::Point                    CDT_Point;
  typedef typename CDT::Edge                     CDT_Edge;
  typedef typename CDT::Face_handle              CDT_Face_handle;
  typedef typename CDT::Vertex_handle            CDT_Vertex_handle;
  typedef typename CDT::Finite_faces_iterator    CDT_Finite_faces_iterator;
  typedef typename CDT::Finite_vertices_iterator CDT_Finite_vertices_iterator;
  typedef typename CDT::Finite_edges_iterator    CDT_Finite_edges_iterator;
  typedef typename CDT::Locate_type              CDT_Locate_type;

protected:

  typedef Arr_traits_basic_adaptor_2<Traits_2>  Traits_adaptor_2;

  // Data members:
  const Traits_adaptor_2  *m_traits;     // Its associated traits object.
  bool                    ignore_notifications;  
  CDT                     cdt;
  bool                    updated_cdt;

public:

  /*! Default constructor. */
  Arr_triangulation_point_location () : 
    m_traits (NULL)
  {
  }

  /*! Constructor given an arrangement. */
  Arr_triangulation_point_location (const Arrangement_2& arr) :
    Arr_observer<Arrangement_2> (const_cast<Arrangement_2 &>(arr))
  {
    build_triangulation();
  }
   
  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  Object locate (const Point_2& p) const;

   //Observer functions that are relevant to overload
   //-------------------------------------------------

    /*! Attach an arrangement object. */
    virtual void before_attach (const Arrangement_2& arr)
    {
      m_traits = static_cast<const Traits_adaptor_2*> (arr.traits());
    }

    virtual void after_attach ()
    {
      build_triangulation();
    }

    virtual void before_detach ()
    {
      clear_triangulation();
    }

   /*!
    * Notification after the arrangement has been assigned with another
    * arrangement.
    * \param u A handle to the unbounded face.
    */
    virtual void after_assign ()
    { 
      clear_triangulation();
      build_triangulation();
    }

    /*!
    * Notification after the arrangement is cleared.
    * \param u A handle to the unbounded face.
    */
    virtual void after_clear ()
    { 
      clear_triangulation();
      build_triangulation();
    }

    /*! Notification before a global operation modifies the arrangement. */
    virtual void before_global_change ()
    { 
      clear_triangulation();
      ignore_notifications = true;
    }

    /*! Notification after a global operation is completed. */
    virtual void after_global_change ()
    {
      build_triangulation();
      ignore_notifications = false;
    }

    /*!
    * Notification after the creation of a new vertex.
    * \param v A handle to the created vertex.
    */
    virtual void after_create_vertex (Vertex_handle /* v */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after the creation of a new edge.
    * \param e A handle to one of the twin halfedges that were created.
    */
    virtual void after_create_edge (Halfedge_handle /* e */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after an edge was split.
    * \param e1 A handle to one of the twin halfedges forming the first edge.
    * \param e2 A handle to one of the twin halfedges forming the second edge.
    */
    virtual void after_split_edge (Halfedge_handle /* e1 */,
      Halfedge_handle /* e2 */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after a face was split.
    * \param f A handle to the face we have just split.
    * \param new_f A handle to the new face that has been created.
    * \param is_hole Whether the new face forms a hole inside f.
    */
    virtual void after_split_face (Face_handle /* f */,
      Face_handle /* new_f */,
      bool /* is_hole */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after a hole was created inside a face.
    * \param h A circulator representing the boundary of the new hole.
    */
    virtual void after_add_hole (Ccb_halfedge_circulator /* h */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after an edge was merged.
    * \param e A handle to one of the twin halfedges forming the merged edge.
    */
    virtual void after_merge_edge (Halfedge_handle /* e */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after a face was merged.
    * \param f A handle to the merged face.
    */
    virtual void after_merge_face (Face_handle /* f */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification after a hole is moved from one face to another.
    * \param h A circulator representing the boundary of the hole.
    */
    virtual void after_move_hole (Ccb_halfedge_circulator /* h */)
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notificaion before the removal of a vertex.
    * \param v A handle to the vertex to be deleted.
    */
    virtual void after_remove_vertex ()
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification before the removal of an edge.
    * \param e A handle to one of the twin halfedges to be deleted.
    */
    virtual void after_remove_edge ()
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

    /*!
    * Notification before the removal of a hole.
    * \param h A circulator representing the boundary of the hole.
    */
    virtual void after_remove_hole ()
    {
      if (! ignore_notifications)
      {
        clear_triangulation();
        build_triangulation();
      }
    }

protected:
  void clear_triangulation ();
  void build_triangulation();  
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_triangulation_pl_functions.h>

#endif
