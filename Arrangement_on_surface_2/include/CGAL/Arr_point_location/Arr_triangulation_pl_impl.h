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
#ifndef CGAL_ARR_TRIANGULATION_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_TRIANGULATION_POINT_LOCATION_FUNCTIONS_H

/*! \file
* Member-function definitions for the Arr_triangulation_point_location<Arrangement>
* class.
*/

//#define CGAL_TRG_DEBUG

#ifdef CGAL_TRG_DEBUG
	#define CGAL_TRG_PRINT_DEBUG(expr)   std::cout << expr << std::endl
#else
	#define CGAL_TRG_PRINT_DEBUG(expr)
#endif

namespace CGAL {

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement_2>
Object Arr_triangulation_point_location<Arrangement_2>
::locate (const Point_2& p) const
{
  CGAL_TRG_PRINT_DEBUG("------ locate point "<< p);

  //init output
  Face_const_handle face_found = this->arrangement()->unbounded_face();

  //locate in the CDT
  CDT_Point p1 = static_cast <CDT_Point> (p);

  //locate point
  int li;
  CDT_Locate_type cdt_lt;
  CDT_Face_handle fh = cdt.locate(p1,cdt_lt,li);

  switch(cdt_lt) {
  case CDT::OUTSIDE_AFFINE_HULL:
  case CDT::OUTSIDE_CONVEX_HULL:
    {
      CGAL_TRG_PRINT_DEBUG("unbounded face" );

      // we still have to check whether the query point coincides with
      // any of the isolated vertices contained inside this face.
      Isolated_vertex_const_iterator   iso_verts_it;
      typename Traits_adaptor_2::Equal_2  equal = m_traits->equal_2_object();

      for (iso_verts_it = face_found->isolated_vertices_begin();
          iso_verts_it != face_found->isolated_vertices_end(); ++iso_verts_it)
      {
        if (equal (p, iso_verts_it->point()))
        {
          Vertex_const_handle  vh = iso_verts_it;
          return (CGAL::make_object (vh));
        }
      }

      return (CGAL::make_object(face_found));
    }
  case CDT::VERTEX:
    {
      //get the vertex from li, which is the index of the vertex
      Vertex_const_handle vertex_found = fh->vertex(li)->info();
      CGAL_TRG_PRINT_DEBUG("vertex: "<< vertex_found->point());
      return (CGAL::make_object(vertex_found));
    }
  case CDT::EDGE:
    {
      CGAL_TRG_PRINT_DEBUG("locate type = edge"<<li );
      //li is the index of the vertex OPOSITE to the edge 
      if ( cdt.is_constrained(CDT_Edge(fh,li)) )
      {  //the edge found is an edge in the plannar map
        CGAL_TRG_PRINT_DEBUG("the edge is a constrained");
        //get the 2 vertices incident to the edge in the plannar map 
        int v1_index = (li+1)%3, v2_index = (li+2)%3;
        CGAL_TRG_PRINT_DEBUG("v1 = "<<v1_index<<", v2 = "<<v2_index );
        Vertex_const_handle v1_of_edge = fh->vertex(v1_index)->info();       
        Vertex_const_handle v2_of_edge = fh->vertex(v2_index)->info();
        //go over all halfedges incident to v1, and check if v2 is their source
        Halfedge_around_vertex_const_circulator circ1 = 
            v1_of_edge->incident_halfedges(); 
        Halfedge_around_vertex_const_circulator circ1_done (circ1);
        
        Halfedge_const_handle edeg_found;

        do {
          if (v2_of_edge == (*circ1).source())
          {
            edeg_found = circ1;
            CGAL_TRG_PRINT_DEBUG("edeg_found = "<< edeg_found->source()->point()
              <<" towards "<< edeg_found->target()->point());
          }
        } while (++circ1 != circ1_done);

        return (CGAL::make_object(edeg_found)); 
      }
      //if the edge is not a constrained - its not an edge of the 
      //plannar map, which means we're inside of a pm face -
      //lets look at the face as if it was a face case.
      // no break - continue to the face caes
    }
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
  if (v0->is_isolated())
    return (CGAL::make_object(v0->face()));
  if (v1->is_isolated())
    return (CGAL::make_object(v1->face()));
  if (v2->is_isolated())
    return (CGAL::make_object(v2->face()));

  //find the face in the pm correspond to the 3 vertices
  Halfedge_around_vertex_const_circulator havc0 = v0->incident_halfedges(); 
  Halfedge_around_vertex_const_circulator havc0_done (havc0);

  Halfedge_around_vertex_const_circulator havc1 = v1->incident_halfedges(); 
  Halfedge_around_vertex_const_circulator havc1_done (havc1);

  Halfedge_around_vertex_const_circulator havc2 = v2->incident_halfedges(); 
  Halfedge_around_vertex_const_circulator havc2_done (havc2);

  //loop to find face
  bool found = false;
  bool found_unbounded = false;
  do {
    //get face from halfedge
    Face_const_handle f0 = (*havc0).face(); 
    do {
      Face_const_handle f1 = (*havc1).face(); 
      if ( f0 == f1 )
      {
        CGAL_TRG_PRINT_DEBUG("f0 == f1");
        do {
          Face_const_handle f2 = (*havc2).face();
          if ( f1 == f2 )
          {
            CGAL_TRG_PRINT_DEBUG("f1 == f2");
            if (face_found != f0) {
              face_found = f0;
              found = true;
            }
            else
              found_unbounded = true;
          }
        } while ((++havc2 != havc2_done) && !found );
      }
    } while ((++havc1 != havc1_done)&& !found );
  } while ((++havc0 != havc0_done)&& !found );

  if (face_found == this->arrangement()->unbounded_face())
  {
    if (! found_unbounded)
    {
      std::cerr<< "NOT GOOD - face not found" << std::endl;
      //debug - print some more info
      std::cout << "p = "<< p <<std::endl;
      std::cout << "v0 = "<< v0->point() 
        <<", v1 = "<< v1->point()
        <<", v2 = "<<v2->point() <<std::endl;
    }
  }

  // we still have to check whether the query point coincides with
  // any of the isolated vertices contained inside this face.
  Isolated_vertex_const_iterator   iso_verts_it;
  typename Traits_adaptor_2::Equal_2  equal = m_traits->equal_2_object();

  for (iso_verts_it = face_found->isolated_vertices_begin();
      iso_verts_it != face_found->isolated_vertices_end(); ++iso_verts_it)
  {
    if (equal (p, iso_verts_it->point()))
    {
      Vertex_const_handle  vh = iso_verts_it;
      return (CGAL::make_object (vh));
    }
  }		

  return (CGAL::make_object(face_found));
}


//----------------------------------------------------
/*! triangulate the arrangement into a cdt (Constaint Delauney Triangulation):
go over all halfedges, and insert each halfedge as a constraint to the cdt. 
*/
template <class Arrangement_2>
void Arr_triangulation_point_location<Arrangement_2>
::clear_triangulation () 
{ 
  cdt.clear();
}

//----------------------------------------------------
/*! triangulate the arrangement into a cdt (Constaint Delauney Triangulation):
go over all halfedges, and insert each halfedge as a constraint to the cdt. 
*/
template <class Arrangement_2>
void Arr_triangulation_point_location<Arrangement_2>
::build_triangulation ()
{ 
  CGAL_TRG_PRINT_DEBUG("build_triangulation");

  //Go over the arrangement, and create a triangulation of it
  Edge_const_iterator eit;

  eit = this->arrangement()->edges_begin();

  for (eit = this->arrangement()->edges_begin(); 
       eit != this->arrangement()->edges_end(); eit++)
  {
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
    if (m_traits->equal_2_object()(pm_p1, pm_p2))
    {
      std::cerr << "WARNING: source point is equal to destination point!!! " 
        << pm_p1 << std::endl ;
      CDT_Vertex_handle cdt_vh1 = cdt.insert(cdt_p1);
      cdt_vh1->info() = pm_vh1;
      continue;
    }

    //insert vertices to the CDT
    CDT_Vertex_handle cdt_vh1 = cdt.insert(cdt_p1);
    CDT_Vertex_handle cdt_vh2 = cdt.insert(cdt_p2);

    //connect new CDT vertex with Pm vertex
    cdt_vh1->info() = pm_vh1;
    cdt_vh2->info() = pm_vh2;

    //add constraint from the two points
    cdt.insert_constraint(cdt_vh1, cdt_vh2);

    //print
    CGAL_TRG_PRINT_DEBUG("source = " << pm_p1 << " , target = " << pm_p2 );
  }

  //the triangulation is now updated
  updated_cdt = true;

  CGAL_assertion(cdt.is_valid());
  CGAL_TRG_PRINT_DEBUG("finished updating the CDT " );
}


} //namespace CGAL

#endif
