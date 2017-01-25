// Copyright (c) 2005,2006,2008,2009,2010,2011 Tel-Aviv University (Israel).
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

#ifndef CGAL_ARR_TRIANGULATION_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_TRIANGULATION_POINT_LOCATION_FUNCTIONS_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
* Member-function definitions for the
* Arr_triangulation_point_location<Arrangement> class.
*/

// #define CGAL_TRG_DEBUG

#ifdef CGAL_TRG_DEBUG
  #define CGAL_TRG_PRINT_DEBUG(expr)   std::cout << expr << std::endl
#else
  #define CGAL_TRG_PRINT_DEBUG(expr)
#endif

namespace CGAL {

template <typename Arrangement_2_>
typename Arr_triangulation_point_location<Arrangement_2_>::result_type
Arr_triangulation_point_location<Arrangement_2_>::
locate_in_unbounded(const Point_2& p) const
{
  CGAL_TRG_PRINT_DEBUG("unbounded face");

  //! \todo Here we assume that there is only one unbounded face.
  Face_const_handle face_found = this->arrangement()->unbounded_faces_begin();

  // Check whether the query point coincides with any of the isolated
  // vertices contained inside this face.
  typename Traits_adaptor_2::Equal_2 equal = m_traits->equal_2_object();
  Isolated_vertex_const_iterator iso_verts_it;
  for (iso_verts_it = face_found->isolated_vertices_begin();
       iso_verts_it != face_found->isolated_vertices_end(); ++iso_verts_it)
  {
    if (equal (p, iso_verts_it->point())) {
      Vertex_const_handle  vh = iso_verts_it;
      return make_result(vh);
    }
  }
  return make_result(face_found);
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <typename Arrangement_2_>
typename Arr_triangulation_point_location<Arrangement_2_>::result_type
Arr_triangulation_point_location<Arrangement_2_>::locate(const Point_2& p)
  const
{
  CGAL_TRG_PRINT_DEBUG("------ locate point " << p);

  //locate in the CDT
  CDT_Point p1 = static_cast <CDT_Point>(p);

  //locate point
  int li;
  CDT_Locate_type cdt_lt;
  CDT_Face_handle fh = m_cdt.locate(p1, cdt_lt, li);

  switch (cdt_lt) {
   case CDT::OUTSIDE_AFFINE_HULL:
   case CDT::OUTSIDE_CONVEX_HULL:
    return locate_in_unbounded(p);

   case CDT::VERTEX:
    //get the vertex from li, which is the index of the vertex
    CGAL_TRG_PRINT_DEBUG("vertex: "<< fh->vertex(li)->info()->point());
    return make_result(fh->vertex(li)->info());

   case CDT::EDGE:
    CGAL_TRG_PRINT_DEBUG("locate type = edge" << li);
    //li is the index of the vertex OPOSITE to the edge
    if (m_cdt.is_constrained(CDT_Edge(fh,li))) {
      //the edge found is an edge in the plannar map
      CGAL_TRG_PRINT_DEBUG("the edge is a constrained");
      //get the 2 vertices incident to the edge in the plannar map
      int v1_index = (li+1)%3, v2_index = (li+2)%3;
      CGAL_TRG_PRINT_DEBUG("v1 = " << v1_index << ", v2 = " << v2_index);
      Vertex_const_handle v1_of_edge = fh->vertex(v1_index)->info();
      Vertex_const_handle v2_of_edge = fh->vertex(v2_index)->info();
      //go over all halfedges incident to v1, and check if v2 is their source
      Halfedge_around_vertex_const_circulator circ1 =
        v1_of_edge->incident_halfedges();
      Halfedge_around_vertex_const_circulator circ1_done (circ1);

      Halfedge_const_handle edge_found;
      do {
        if (v2_of_edge == (*circ1).source()) {
          edge_found = circ1;
          CGAL_TRG_PRINT_DEBUG("edge_found = "
                               << edge_found->source()->point()
                               << " towards "
                               << edge_found->target()->point());
        }
      } while (++circ1 != circ1_done);
      return make_result(edge_found);
    }
    // if the edge is not a constrained - its not an edge of the
    // plannar map, which means we're inside of a pm face -
    // lets look at the face as if it was a face case.
    // no break - continue to the face caes

   case CDT::FACE:
    break;
  }

  //we're in case CDT::FACE
  CGAL_TRG_PRINT_DEBUG("FACE ");

  //get 3 pm vertices of face
  Vertex_const_handle v0 = fh->vertex(0)->info();
  Vertex_const_handle v1 = fh->vertex(1)->info();
  Vertex_const_handle v2 = fh->vertex(2)->info();

  //the vertices should not be isolated, since we do not insert the
  //isolated vertices as points in the triangulation, only edges
  // (and thus vertices inceident to this edge).
  //in the future it is possible to add isolated vertices to the
  // triangulation, and then, when found, take its incident_face
  CGAL_assertion(!v0->is_isolated());
  CGAL_assertion(!v1->is_isolated());
  CGAL_assertion(!v2->is_isolated());
  // if (v0->is_isolated()) return make_result(v0->face());
  // if (v1->is_isolated()) return make_result(v1->face());
  // if (v2->is_isolated()) return make_result(v2->face());

  // Find the face in the arrangement that contains the triangle <v0,v1,v2>
  // Loop over the incident vertices of v0, and try to find a ccb that first
  // contains v1 and then v2 in that order. If such a face does not exists,
  // the unbounded face contains the triangle.
  Halfedge_around_vertex_const_circulator havc0 = v0->incident_halfedges();
  Halfedge_around_vertex_const_circulator havc0_done(havc0);
  bool found_v2 = false;
  Halfedge_const_handle he;
  do {
    bool found_v1 = false;
    found_v2 = false;
    he = havc0->twin();
    Halfedge_const_handle he_done(he);
    do {
      if (!found_v1) {
        if (he->target() == v1) found_v1 = true;
      }
      else if (!found_v2) {
        if (he->target() == v2) {
          found_v2 = true;
          break;
        }
      }
      he = he->next();
    } while (he != he_done);
  } while ((++havc0 != havc0_done) && !found_v2);

  Face_const_handle face_found = (found_v2) ?
    he->face() : this->arrangement()->unbounded_faces_begin();

  // we still have to check whether the query point coincides with
  // any of the isolated vertices contained inside this face.
  typename Traits_adaptor_2::Equal_2 equal = m_traits->equal_2_object();
  Isolated_vertex_const_iterator iso_verts_it;
  for (iso_verts_it = face_found->isolated_vertices_begin();
      iso_verts_it != face_found->isolated_vertices_end(); ++iso_verts_it)
  {
    if (equal (p, iso_verts_it->point())) {
      Vertex_const_handle  vh = iso_verts_it;
      return make_result(vh);
    }
  }

  return make_result(face_found);
}


//----------------------------------------------------
/*! triangulate the arrangement into a cdt (Constaint Delauney Triangulation):
go over all halfedges, and insert each halfedge as a constraint to the cdt.
*/
template <typename Arrangement_2_>
void Arr_triangulation_point_location<Arrangement_2_>::clear_triangulation()
{ m_cdt.clear(); }

//----------------------------------------------------
/*! triangulate the arrangement into a cdt (Constaint Delauney Triangulation):
go over all halfedges, and insert each halfedge as a constraint to the cdt.
*/
template <typename Arrangement_2_>
void Arr_triangulation_point_location<Arrangement_2_>::build_triangulation()
{
  CGAL_TRG_PRINT_DEBUG("build_triangulation");

  //Go over the arrangement, and create a triangulation of it
  Edge_const_iterator eit = this->arrangement()->edges_begin();
  for (; eit != this->arrangement()->edges_end(); ++eit) {
    //get vertices from edge
    Vertex_const_handle pm_vh1 = (*eit).source();
    Vertex_const_handle pm_vh2 = (*eit).target();

    //get curve
    X_monotone_curve_2 cv = (*eit).curve();

    //get points from vertices
    Point_2 pm_p1 = pm_vh1->point() ;
    Point_2 pm_p2 = pm_vh2->point() ;

    //cast the points to be CDT points
    CDT_Point cdt_p1 = static_cast <CDT_Point> (pm_p1);
    CDT_Point cdt_p2 = static_cast <CDT_Point> (pm_p2);

    //check if source point is equal to destination point
    if (m_traits->equal_2_object()(pm_p1, pm_p2)) {
      std::cerr << "WARNING: source point is equal to destination point!!! "
                << pm_p1 << std::endl ;
      CDT_Vertex_handle cdt_vh1 = m_cdt.insert(cdt_p1);
      cdt_vh1->info() = pm_vh1;
      continue;
    }

    //insert vertices to the CDT
    CDT_Vertex_handle cdt_vh1 = m_cdt.insert(cdt_p1);
    CDT_Vertex_handle cdt_vh2 = m_cdt.insert(cdt_p2);

    //connect new CDT vertex with Pm vertex
    cdt_vh1->info() = pm_vh1;
    cdt_vh2->info() = pm_vh2;

    //add constraint from the two points
    m_cdt.insert_constraint(cdt_vh1, cdt_vh2);

    //print
    CGAL_TRG_PRINT_DEBUG("source = " << pm_p1 << " , target = " << pm_p2);
  }

  CGAL_assertion(m_cdt.is_valid());
  CGAL_TRG_PRINT_DEBUG("finished updating the CDT ");
}

} //namespace CGAL

#endif
