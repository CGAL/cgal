// Copyright (c) 2004  Tel-Aviv University (Israel).
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
// Author(s)     : Idit Haran <haranidi@post.tau.ac.il>

#ifndef CGAL_PM_TRIANGLE_POINT_LOCATION_C
#define CGAL_PM_TRIANGLE_POINT_LOCATION_C

#include <CGAL/Pm_triangle_point_location.h>

//#define CGAL_PM_DEBUG

CGAL_BEGIN_NAMESPACE

//debug
#ifdef CGAL_PM_DEBUG
  int point_number = 0;
  int constr_number = 0;
#endif


//----------------------------------------------------
/*!
 */
//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
template <class Planar_map>
typename Pm_triangle_point_location<Planar_map>::Halfedge_const_handle
Pm_triangle_point_location<Planar_map>::
locate(const typename Planar_map::Traits::Point_2 & p, 
       Locate_type & lt) const
{
#ifdef CGAL_PM_DEBUG
  point_number++;
  std::cout << std::endl << "LOCATE NUMBER "<< point_number << std::endl;
  std::cout << "------ locate point "<< p << std::endl;
#endif

  //init output
  Face_handle f = pm->unbounded_face(), last = pm->faces_end();  
              // invariant: p always lies in f's interior or holes 
  Halfedge_handle e = pm->halfedges_end(); // closest halfedge so far
  lt = Planar_map::UNBOUNDED_FACE;

  //init CDT
  if (! updated_cdt)
  {
    std::cout << "ERROR: CDT is not updated !" << std::endl;
    return e;
  }
  //assert(cdt.is_valid());

  //locate in the CDT
  CDT_Point p1 = static_cast <CDT_Point> (p);
  //CDT_Face_handle fh1 = cdt.locate(p1);       //locate point
  //fh1->info() = CGAL::YELLOW;               //update the face's color

  //locate point
  int li;
  CDT_Locate_type cdt_lt;
  CDT_Face_handle fh = cdt.locate(p1,cdt_lt,li);
  switch(cdt_lt) {
  case CDT::OUTSIDE_AFFINE_HULL:
  case CDT::OUTSIDE_CONVEX_HULL:
    lt = Planar_map::UNBOUNDED_FACE;
#ifdef CGAL_PM_DEBUG
    std::cout << "unbounded face" << std::endl;
#endif
    return e;
  case CDT::VERTEX:
    {
#ifdef CGAL_PM_DEBUG
      std::cout << "locate type = vertex"<<li << std::endl;
#endif
      //get the vertex from li, which is the index of the vertex
      Vertex_handle vertex_found = fh->vertex(li)->info();
      e = vertex_found->incident_halfedges();
      lt = Planar_map::VERTEX;
#ifdef CGAL_PM_DEBUG
      std::cout << "vertex is "<< vertex_found->point() << std::endl;
      std::cout << "e = "<< e->source()->point()
		<<" towards "<< e->target()->point() << std::endl;
#endif
      return e;
    }
  case CDT::EDGE:
    {
#ifdef CGAL_PM_DEBUG
      std::cout << "locate type = edge"<<li << std::endl;
#endif
      //li is the index of the vertex OPOSITE to the edge 
      if ( cdt.is_constrained(CDT_Edge(fh,li)) )
      {  //the edge found is an edge in the plannar map
#ifdef CGAL_PM_DEBUG
	std::cout << "the edge is a constrained" << std::endl;
#endif
	lt = Planar_map::EDGE;
	//get the 2 vertices incident to the edge in the plannar map 
	int v1_index = (li+1)%3, v2_index = (li+2)%3;
#ifdef CGAL_PM_DEBUG
	std::cout << "v1 = "<<v1_index<<", v2 = "<<v2_index << std::endl;
#endif
	Vertex_handle v1_of_edge = fh->vertex(v1_index)->info();       
	Vertex_handle v2_of_edge = fh->vertex(v2_index)->info();
	//go over all halfedges incident to v1, and check if v2 is their source
	Halfedge_around_vertex_circulator circ1 = v1_of_edge->incident_halfedges(); 
	Halfedge_around_vertex_circulator circ1_done (circ1);

	//std::cout << "Halfedges incident to " << v1_of_edge->point() << std::endl;
	do {
	  //std::cout << "edge = "<< (*circ1).curve() << std::endl;    
	  //std::cout << "source = "<<circ1->source()->point() << std::endl; 
	  //std::cout << "target = "<<circ1->target()->point() << std::endl;
	  if (v2_of_edge == circ1->source())
	  {
	    e = circ1;
#ifdef CGAL_PM_DEBUG
	    std::cout << "e = "<< e->source()->point()
		      <<" towards "<< e->target()->point() << std::endl;
#endif
	  }
	} while (++circ1 != circ1_done);
	
	//check if e is valid
#ifdef CGAL_PM_DEBUG
	if (pm->halfedges_end() == e)
	  std::cout << "ERROR: e is invalid " << std::endl;
#endif	
	return e; // ????
      }
      //if the edge is not a constrained - its not an edge of the 
      //plannar map, which means we're inside of a pm face -
      //lets look at the face as if it was a face case.
      // no break - continue to the face caes
    }
  case CDT::FACE:
    lt = Planar_map::FACE;
    break;
  }

  //we're in case CDT::FACE
#ifdef CGAL_PM_DEBUG
  std::cout << "FACE "  << std::endl;
#endif

  //get 3 pm vertices of face
  Vertex_handle v0 = fh->vertex(0)->info();       
  Vertex_handle v1 = fh->vertex(1)->info();
  Vertex_handle v2 = fh->vertex(2)->info();

  //find the face in the pm correspond to the 3 vertices
  Halfedge_around_vertex_circulator havc0 = v0->incident_halfedges(); 
  Halfedge_around_vertex_circulator havc0_done (havc0);

  Halfedge_around_vertex_circulator havc1 = v1->incident_halfedges(); 
  Halfedge_around_vertex_circulator havc1_done (havc1);

  Halfedge_around_vertex_circulator havc2 = v2->incident_halfedges(); 
  Halfedge_around_vertex_circulator havc2_done (havc2);

#ifdef CGAL_PM_DEBUG
  //printings, debugging
  std::cout << "Halfedges incident to " << v0->point() << std::endl;
  do {
    std::cout << (*havc0).curve() << std::endl;    
  } while (++havc0 != havc0_done);

  std::cout << "Halfedges incident to " << v1->point() << std::endl;
  do {
    std::cout << (*havc1).curve() << std::endl;    
  } while (++havc1 != havc1_done);

  std::cout << "Halfedges incident to " << v2->point() << std::endl;
  do {
    std::cout << (*havc2).curve() << std::endl;    
  } while (++havc2 != havc2_done);
#endif //CGAL_PM_DEBUG

  //loop to find face
  Face_iterator found_face = pm->unbounded_face();
  bool found = false;
  bool found_unbounded = false;
  do {
    //get face from halfedge
    Face_iterator f0 = (*havc0).face(); 
    //printf ("f0 = %e \n", f0);
    //Pm_Face_iterator f0_2 = (*havc0).twin().face();
    do {
      Face_iterator f1 = (*havc1).face(); 
      //printf ("f1 = %e \n", f0);
      if ( f0 == f1 )
      {
#ifdef CGAL_PM_DEBUG
	std::cout <<"f0 == f1"<< std::endl;
#endif
	do {
	  Face_iterator f2 = (*havc2).face();
	  if ( f1 == f2 )
	  {
#ifdef CGAL_PM_DEBUG
	    std::cout << "f1 == f2"<< std::endl;
#endif
	    if (found_face != f0) {
	      found_face = f0;
	      e = havc0; //!!! this is the relevant update for the PL !!!
	      found = true;
	    }
	    else
	      found_unbounded = true;
	  }
	  //printf ("f2 = %e \n", f0);
	} while ((++havc2 != havc2_done) && !found );
      }
    } while ((++havc1 != havc1_done)&& !found );
  } while ((++havc0 != havc0_done)&& !found );

  if (found_face == pm->unbounded_face())
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
    lt = Planar_map::UNBOUNDED_FACE;
    e = pm->halfedges_end();
    return e;
  }

#ifdef CGAL_PM_DEBUG
  //print the curves of the found face
  std::cout << "found face for point "<< p1 <<" is: "<<std::endl;
  Ccb_halfedge_circulator ccb = found_face->outer_ccb();
  Ccb_halfedge_circulator ccb_end = ccb;
  do {
    std::cout << ccb->source()->point()<<" towards "
              << ccb->target()->point()<<std::endl;
  } while ((++ccb) != ccb_end);
  
    std::cout << "e = "<< e->source()->point()
	      <<" towards "<< e->target()->point() << std::endl;
#endif //CGAL_PM_DEBUG

  return e;
}

//----------------------------------------------------
/*!
 */
template <class Planar_map>
typename Pm_triangle_point_location<Planar_map>::Halfedge_handle
Pm_triangle_point_location<Planar_map>::
locate(const typename Planar_map::Traits::Point_2 & p, 
       Locate_type & lt) {

  //  if (! updated_cdt)
  //    triangulate_pm();
  //  assert(cdt.is_valid());

  Halfedge_handle h = Halfedge_handle_unconst(((cPLp)this)->locate(p,lt));
  return h;
}


//----------------------------------------------------
/*!
 */
template <class Planar_map>
typename Pm_triangle_point_location<Planar_map>::Halfedge_const_handle
Pm_triangle_point_location<Planar_map>::
vertical_ray_shoot(const typename Planar_map::Traits::Point_2 & p, 
                   Locate_type & lt, 
                   bool up) const
{
  std::cout << "ERROR: Vertical ray shoot NOT suppored in CDT point location"
	    << std::endl;
  CGAL_precondition_msg(false, 
			"Vertical ray shoot NOT suppored in CDT point location");
  assert(false);
  Halfedge_handle e = pm->halfedges_end(); // closest halfedge so far  
  lt = Planar_map::UNBOUNDED_FACE;
  return e;
}

//----------------------------------------------------
/*!
 */
template <class Planar_map>
typename Pm_triangle_point_location<Planar_map>::Halfedge_handle
Pm_triangle_point_location<Planar_map>::
vertical_ray_shoot(const typename Planar_map::Traits::Point_2 & p, 
                   Locate_type & lt, bool up)
{
  Halfedge_handle h =
    Halfedge_handle_unconst(((cPLp)this)->vertical_ray_shoot(p,lt,up));
  return h;
}

//----------------------------------------------------
/*! insert an halfedge to the plannar map - 
    insert it as a constraint to the CDT
 */
template <class Planar_map>
void Pm_triangle_point_location<Planar_map>::insert_to_cdt
(Halfedge_handle hh, const typename Planar_map::Traits::X_monotone_curve_2 &cv)
{
#ifdef CGAL_PM_DEBUG
  constr_number++; 
  std::cout << "in insert_to_cdt, num = "<< constr_number << std::endl ;
  //std::cout << *pm << std::endl;
  std::cout << "cv = "<< cv << std::endl;
  std::cout << hh->source()->point()<<" towards "
  	    << hh->target()->point()<< std::endl;
#endif //CGAL_PM_DEBUG
  
  Vertex_handle pm_vh1 = hh->source();
  Vertex_handle pm_vh2 = hh->target();

  //get points from vertices
  typename Planar_map::Traits::Point_2 pm_p1 = (*pm_vh1).point() ;
  typename Planar_map::Traits::Point_2 pm_p2 = (*pm_vh2).point() ;

  //cast the points to be CDT points
  CDT_Point cdt_p1 = static_cast <CDT_Point> (pm_p1);
  CDT_Point cdt_p2 = static_cast <CDT_Point> (pm_p2);

  //insert vertices to the CDT
  CDT_Vertex_handle cdt_vh1 = cdt.insert(cdt_p1);
  CDT_Vertex_handle cdt_vh2 = cdt.insert(cdt_p2);

  //connect new CDT vertex with Pm vertex
  cdt_vh1->info() = pm_vh1;
  cdt_vh2->info() = pm_vh2;

#ifdef CGAL_PM_DEBUG
  std::cout << "call to insert_constraint with "
	    << cdt_vh1->point() << cdt_vh2->point() << std::endl ;
#endif //CGAL_PM_DEBUG

  //add constraint from the two points
  cdt.insert_constraint(cdt_vh1, cdt_vh2);

#ifdef CGAL_PM_DEBUG
  std::cout << "finished insert_to_cdt" << std::endl ;
#endif //CGAL_PM_DEBUG
}

//----------------------------------------------------
/*! triangulate the plannar map into a cdt (Constaint
    Delauney Triangulation) :
    go over all halfedges, and insert each halfedge 
    as a constraint to the cdt. 
 */
template <class Planar_map>
void Pm_triangle_point_location<Planar_map>::triangulate_pm()
{ 
#ifdef CGAL_PM_DEBUG
  std::cout << "in triangulate_pm" << std::endl ;
  //std::cout << *pm << std::endl;
#endif
  
  cdt.clear();

  //Go over plannar map, and create a triangulation of it
  Edge_iterator eit;

#ifdef CGAL_PM_DEBUG
  std::cout << "All edges" << std::endl;
#endif

  eit=pm->edges_begin();
  //  std::cout << "0" << std::endl;

  for (eit=pm->edges_begin(); eit != pm->edges_end(); eit++)
  {
    //std::cout << "1" << std::endl;
    //get vertices from edge
    Vertex_handle pm_vh1 = (*eit).source();
    Vertex_handle pm_vh2 = (*eit).target();

    //    std::cout << "2" << std::endl;
    //get curve
    X_monotone_curve_2 cv = (*eit).curve();
    //    std::cout << "cv = "<< cv << std::endl;

    //get points from vertices
    typename Planar_map::Traits::Point_2 pm_p1 = (*pm_vh1).point() ;
    typename Planar_map::Traits::Point_2 pm_p2 = (*pm_vh2).point() ;

    //    std::cout << "3" << std::endl;
    //cast the points to be CDT points
    CDT_Point cdt_p1 = static_cast <CDT_Point> (pm_p1);
    CDT_Point cdt_p2 = static_cast <CDT_Point> (pm_p2);

    //    std::cout << "4" << std::endl;
    //insert vertices to the CDT
    CDT_Vertex_handle cdt_vh1 = cdt.insert(cdt_p1);
    CDT_Vertex_handle cdt_vh2 = cdt.insert(cdt_p2);

    //    std::cout << "5" << std::endl;
    //connect new CDT vertex with Pm vertex
    cdt_vh1->info() = pm_vh1;
    cdt_vh2->info() = pm_vh2;

    //    std::cout << "6" << std::endl;
    //    std::cout << "pm_p1 = "<< pm_p1 <<" , cdt_p1 = "<< cdt_p1 << std::endl;
    //    std::cout << "pm_p2 = "<< pm_p2 <<" , cdt_p2 = "<< cdt_p2 << std::endl;
    //std::cout << "pm_vh1 = "<< pm_vh1 <<" , cdt_vh1 = "<< cdt_vh1 << std::endl;
    //std::cout << "pm_vh2 = "<< pm_vh2 <<" , cdt_vh2 = "<< cdt_vh2 << std::endl;
   
    //add constraint from the two points
    cdt.insert_constraint(cdt_vh1, cdt_vh2);

    //print
    //    std::cout << "7" << std::endl;
#ifdef CGAL_PM_DEBUG
    std::cout << "curve = " << cv ;
    std::cout << " , source = " << pm_p1 ;
    std::cout << " , target = " << pm_p2 << std::endl;
#endif
  }
  

 //the triangulation is now updated
  updated_cdt = true;

  // --- printing and debugging of CDT
  //  std::cout << "cdt is now updated" << std::endl;
  assert(cdt.is_valid());

#ifdef CGAL_PM_DEBUG
  int count = 0;
  
  for (CDT_Finite_edges_iterator eit = cdt.finite_edges_begin();
       eit != cdt.finite_edges_end();
       ++eit)
    if (cdt.is_constrained(*eit)) ++count;
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;

  CDT_Finite_faces_iterator fc = cdt.finite_faces_begin();
  for( ; fc != cdt.finite_faces_end(); ++fc)  
  {
    fc->info() = CGAL::BLUE;
  }

  count=0;
  fc = cdt.finite_faces_begin();
  for( ; fc != cdt.finite_faces_end(); ++fc)  
  {
    count++;
    std::cout << "fc number "<<count<<" ,color = " << fc->info() << std::endl;
  }

#endif  //CGAL_PM_DEBUG

}

CGAL_END_NAMESPACE

#endif  //CGAL_PM_TRIANGLE_POINT_LOCATION_C
