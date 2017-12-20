// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Naama mayer         <naamamay@post.tau.ac.il>


// This includes a function for rotating an arrangement on a sphere with
// geodesic arcs.
// It should be templated with a transformation (affine transformations) and a
// transformation traits, defines the way the faces should be rotated.


#ifndef CGAL_ARR_SGM_3_ATOS_H
#define CGAL_ARR_SGM_3_ATOS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm.h>

namespace CGAL {

template < class Arrangement, class Transformation_3, class TransformTraits>
void Arr_transform_on_sphere(Arrangement & arr,
                             const Transformation_3 & aff,
                             TransformTraits & tran_tr)
{
  typedef typename Arrangement::Geometry_traits_2         Geometry_traits_2;
  typedef typename Arrangement::Topology_traits           Topology_traits;
    
  typedef typename Geometry_traits_2::Curve_2             Curve_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
		
  typedef typename Arrangement::Vertex_handle             Vertex_handle;
  typedef typename Arrangement::Halfedge_handle           Halfedge_handle;
  typedef typename Arrangement::Face_handle               Face_handle;
  typedef typename Arrangement::Edge_iterator             Edge_iterator;
  typedef typename Arrangement::Halfedge_around_vertex_circulator
    Halfedge_around_vertex_circulator;

  const Geometry_traits_2 * geom_traits = arr.geometry_traits();
  Topology_traits * topol_traits = arr.topology_traits();

  Arr_accessor<Arrangement> m_arr_access(arr);
	
  // Preprocessing loop - merge all the edges that were splited 
  // (meaning have a common endpoint that lies on the boundary and their degree
  // is 2) on the identification curve.
  for (Vertex_handle vi1 = arr.vertices_begin() ; vi1 != arr.vertices_end() ;)
  {	
    Vertex_handle v_temp = vi1;
    ++vi1;
		
    Arr_parameter_space bx =
      geom_traits->parameter_space_in_x_2_object()(v_temp->point());
    Arr_parameter_space by =
      geom_traits->parameter_space_in_y_2_object()(v_temp->point());
		
    // use vertex->parameter_space_in_x() != interior || vertex->parameter_space_in_y() != interior)
    if ((bx != ARR_INTERIOR || by != ARR_INTERIOR) && (v_temp->degree() == 2))
    {
      Curve_2 merged_cv;
      Halfedge_around_vertex_circulator havc = v_temp->incident_halfedges(); 
      Halfedge_around_vertex_circulator havc_next = havc;
      ++havc_next;

      // Compare the normals
      bool normal_eq1 = havc->curve().normal() ==
        havc_next->twin()->curve().normal();				
      // Compare the points
      bool point_eq1 = havc->target()->point() ==
        havc_next->twin()->source()->point();
      // Compare the direction of the edges.
      bool eq_direction = havc->direction() == havc_next->twin()->direction();

      if (point_eq1 && normal_eq1 && eq_direction)
      {
        if (havc->source()->point() == havc->curve().source())
        {
          merged_cv = Curve_2(havc->source()->point(),
                              havc_next->twin()->target()->point(), 
                              havc->curve().normal());
        }
        else if (havc->source()->point() == havc->curve().target())
        {
          merged_cv = Curve_2(havc_next->twin()->target()->point(),
                              havc->source()->point(),
                              havc->curve().normal());
        }
        else 
          CGAL_error_msg
            ("One of the edge points should be equal to the surce points");
					
        // Erases the point from the list of boundary vertices (This is a
        // different data structute than the DCEL, which is handled by the
        // next call.)
        topol_traits->erase_redundant_vertex(&(*v_temp));

        // Merge the edges into a single one, and delete the vertex from the
        // DCEL. (By default, the merge_edge() funtion deletes the vertex.)
        arr.merge_edge(havc, havc_next->twin() , merged_cv);			
      }
    }
  }

  //Rotate all the vertices.
  for (Vertex_handle vi1 = arr.vertices_begin(); vi1 != arr.vertices_end() ;
       ++vi1)
  {			
    m_arr_access.modify_vertex_ex(vi1, aff.transform(vi1->point()));	
  }		
	
  unsigned int num_of_edges = arr.number_of_edges();
  Edge_iterator ei1 = arr.edges_begin();
	
  // Rotate all the halfedges.
  // The loop is over the initial edges , since new edges are created and
  // added to the edges list.
  for (unsigned int i=0 ; i < num_of_edges ; ++i)
  {
    Curve_2 new_cv;

    Halfedge_handle hei1 = ei1;
    ++ei1;
		
    // Take only the halfedge that its source is equal to the source of
    // the curve.
    bool eq1 = hei1->source()->point() == aff.transform(hei1->curve().source());
    if (!eq1)
    {
      hei1 = hei1->twin();
      eq1 = hei1->source()->point() == aff.transform(hei1->curve().source());
    }
		
    CGAL_assertion_msg
      (eq1, 
       "The new curve endpoints should be equal to the transform of the original");
		
    // Create a new curve instead of the old one.
    new_cv = Curve_2(hei1->source()->point(), hei1->target()->point(), 
                     aff.transform(hei1->curve().normal()));

    // Modify the edge with the new curve.
    m_arr_access.modify_edge_ex(hei1, new_cv);

    std::list<CGAL::Object> objects;
    // Try to split the curve into x_monotone pieces. 
    geom_traits->make_x_monotone_2_object()(new_cv , std::back_inserter(objects));
		
    // If the curve is not x-monotone - split it into 2 x_monotone parts.
    // Since the curves were x_monotone before , can assume that it will be
    // splited into 2 parts max.
    if (objects.size() == 2)
    {						
      typename std::list<CGAL::Object>::iterator it = objects.begin();
			
      // The curve that its left vertex lies on the identification curve
      const X_monotone_curve_2 * sub_cv1 =
        object_cast<X_monotone_curve_2>(&(*it));				
      ++it;
      //The curve that its rigth vertex lies on the identification curve
      const X_monotone_curve_2 * sub_cv2 =
        object_cast<X_monotone_curve_2>(&(*it));
			
      bool eq1 = (*sub_cv1).source() == hei1->source()->point();
      bool eq2 = (*sub_cv2).target() == hei1->target()->point();
			
      if (eq1 && eq2)
        m_arr_access.split_edge_ex (hei1, (*sub_cv1).target(), *sub_cv1,
                                    *sub_cv2);
      else
        CGAL_error_msg
          ("The new curve endpoints should be equal to the original ones");
    }
  }
	
  // Update all the vertices that located on the boundary after the rotation
  // with boundary conditions
  for (Vertex_handle vi = arr.vertices_begin() ; vi != arr.vertices_end() ;
       ++vi)
  {			
    Arr_parameter_space bx =
      geom_traits->parameter_space_in_x_2_object()(vi->point());
    Arr_parameter_space by =
      geom_traits->parameter_space_in_y_2_object()(vi->point());

    if (bx != ARR_INTERIOR || by != ARR_INTERIOR)
    {	
      // The target of the Halfedge_around_vertex_circulator is the relevant
      // point.
      Halfedge_around_vertex_circulator havc = vi->incident_halfedges();
			
      Arr_curve_end ind;
	  if (geom_traits->construct_min_vertex_2_object()(havc->curve()) == vi->point()) 
		ind = ARR_MIN_END;
	  else 
		ind = ARR_MAX_END;

      // Check if it was already added.
      if (topol_traits->discontinuity_vertex(havc->curve(), ind)== NULL &&
          topol_traits->south_pole() != &(*havc->target()) &&
          topol_traits->north_pole() != &(*havc->target()) )
      {
        // Update the boundary conditions for the vertex.
        topol_traits->notify_on_boundary_vertex_creation(&(*havc->target()),
                                                         havc->curve(), ind,
                                                         bx, by);
        m_arr_access.set_vertex_boundary(havc->target(), bx, by);
		
      }
    }	
  }
	
  // Transform the faces with the suitable transformation traits.
  for(Face_handle f1 = arr.faces_begin() ; f1 != arr.faces_end() ; ++f1)
    tran_tr.rotate_face(f1, aff);
}

} //namespace CGAL

#endif // CGAL_ARR_SGM_3_ATOS_H
