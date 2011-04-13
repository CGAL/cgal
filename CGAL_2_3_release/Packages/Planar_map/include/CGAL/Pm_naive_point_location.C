// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Pm_naive_point_location.C
// package       : pm (4.08)
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
//                 Oren Nechushtan <theoren@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_NAIVE_POINT_LOCATION_C
#define CGAL_PM_NAIVE_POINT_LOCATION_C

#ifndef CGAL_PM_NAIVE_POINT_LOCATION_H
#include <CGAL/Pm_naive_point_location.h>
#endif // CGAL_PM_NAIVE_POINT_LOCATION_H

CGAL_BEGIN_NAMESPACE

//IMPLEMENTATION
//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
template <class Planar_map>
Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::locate(const Point& p, 
					    Locate_type& lt) const{
  typename Planar_map::Vertex_iterator vit=pm->vertices_begin();
  for (; vit!=pm->vertices_end(); ++vit) {
    if (traits->point_is_same(p,vit->point()) ) {
      lt = Planar_map::VERTEX; 
      Halfedge_handle h(vit->incident_halfedges());	
      return h;
    }
  }
	
  typename Planar_map::Halfedge_iterator hit=pm->halfedges_begin();
  for (; hit!=pm->halfedges_end(); ++hit) {
    if (traits->curve_get_point_status(hit->curve(),p)==Traits::ON_CURVE) {
      lt = Planar_map::EDGE; 
      return hit;
    }
  }
	
  lt=Planar_map::UNBOUNDED_FACE;
  Locate_type temp;
  Halfedge_handle h = vertical_ray_shoot(p,temp,true);
  if( temp!=Planar_map::UNBOUNDED_FACE )       
    {
      if (temp==Planar_map::VERTEX) {  
	//since h points at the vertex and is the first 
	//halfedge after the ray clockwise! then the  face
	//is to its _right_ (maybe the specs will change in the future) 
	h=h->twin();        
      }        
		
      if ( !(h->face()->is_unbounded()) ) 
	lt=Planar_map::FACE;
      return h;
    }
  else //==the vertical ray shoot returned the halfedges_end() iterator.
    {
      if (pm->unbounded_face()->holes_begin() == 
	  pm->unbounded_face()->holes_end() ) //an empty map
	return h; //return halfedges_end()
      else {
	//- returns a halfedge on an inner ccb of the unbounded face
	typename Planar_map::Holes_iterator hot=
	  pm->unbounded_face()->holes_begin();
	return (*hot);
      }
    }
}

template <class Planar_map>
Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::locate(const Point& p, Locate_type& lt){
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->locate(p,lt);
  if (!((Bounding_box*)get_bounding_box())->locate(p,lt,h))
    h=((cPLp)this)->locate(p,lt);
  return h;
}

template <class Planar_map>
Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::vertical_ray_shoot(const Point& p, 
							Locate_type& lt, 
							bool up) const{

  typename Planar_map::Halfedge_iterator it=pm->halfedges_begin(),
    eit=pm->halfedges_end(),closest_edge=eit;
  bool first = false;
  typename Traits::Curve_point_status point_above_under;
  int curve_above_under;
	
  lt=Planar_map::EDGE;
	
  // set the flags for comparison acording to the ray 
  // direction (up/down)
  if (up) 
    {
      point_above_under = Traits::UNDER_CURVE;
      curve_above_under = LARGER;
    } 
  else 
    {
      point_above_under = Traits::ABOVE_CURVE;
      curve_above_under = SMALLER;
    }
  for ( ; it != eit; ) 
    {
      if ( traits->curve_get_point_status(it->curve(), p) 
	   == point_above_under ) 
        {
	  if (!first) 
            {
	      closest_edge = it;
	      first = true;
            } 
	  else 
            {
	      if ( traits->curve_compare_at_x(closest_edge->curve(),
					      it->curve(), p) ==
		   curve_above_under) 
                {
		  closest_edge = it;
                }
            }
        }
      ++it;++it;
    }
	
  // if we didn't find any edge above p then it is the empty face
  if (!first) {
    lt=Planar_map::UNBOUNDED_FACE;
    Halfedge_handle h=pm->halfedges_end();
    return h; //==NULL
  }
  // if the closest point is a vertex then find the first clockwise 
  // edge from the vertical segment
  typename Planar_map::Vertex_handle v = pm->vertices_end();
  bool maybe_vertical=false; // BUG fix (Oren)
  if ( traits->point_is_same_x(closest_edge->target()->point(), p) ) 
    {
      v = closest_edge->target();
      maybe_vertical=true; // BUG fix (Oren)
    }
	
  if ( traits->point_is_same_x( closest_edge->source()->point(), p) ) 
    {
      if (!maybe_vertical || 
	  traits->point_is_higher(closest_edge->target()->point(),
				  closest_edge->source()->point())==up) 
                                  // BUG fix (Oren)
	v = closest_edge->source();
      /*
	special care for the vertical cases:
		  
	x             p
	|
	x     and     x
	|
	p             x
      */
    }
	
	
  //if (closest_is_vertex)
  if (v != pm->vertices_end()) 
    {
      lt=Planar_map::VERTEX;
      if (up)    
	closest_edge = find_lowest(v,traits, false);  
      else
	closest_edge = find_lowest(v,traits, true);
    }
	
  Halfedge_handle h;	
  if (lt==Planar_map::VERTEX)
    {
      h=closest_edge;
    }
  else if (up) 
    {
      // return the edge that is going from right to left
      // such that p is to the left of this edge
		
      if ( traits->point_is_right( closest_edge->source()->point(),
				   closest_edge->target()->point()) )
	{
	  h=closest_edge;  //source is right of the target
	}
      else
	h=closest_edge->twin();
    } 
  else 
    {
      if ( traits->point_is_left( closest_edge->source()->point(),
				  closest_edge->target()->point()) )
	h=closest_edge;
      else
	h=closest_edge->twin();
    }
  return h;
}

template <class Planar_map>
Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::vertical_ray_shoot(const Point& p, 
							Locate_type& lt, 
							bool up){
  /* Make sure the source point is in the bounding box on the output */
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
  /* Apply the bounding box on the output */
  if (!((Bounding_box*)get_bounding_box())->vertical_ray_shoot(p,lt,up,h))
    {
      h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
      CGAL_assertion(lt!=Planar_map::UNBOUNDED_FACE);
    }
  return h;
}


//find the first halfedge pointing at v, when going clockwise
//if highest==true - start from 12 oclock, else start from 6 oclock
template <class Planar_map>
Pm_naive_point_location<Planar_map>::Halfedge_handle 
Pm_naive_point_location<Planar_map>::
find_lowest(
	    typename Pm_naive_point_location<Planar_map>::Vertex_handle v,
	    typename Pm_naive_point_location<Planar_map>::Traits_wrap *traits, 
	    bool highest) const{
	
  Halfedge_handle lowest_left = pm->halfedges_end();
  Halfedge_handle lowest_right = pm->halfedges_end();
  Halfedge_handle vertical_up = pm->halfedges_end();
  Halfedge_handle vertical_down = pm->halfedges_end();
	
	
  typename Planar_map::Halfedge_around_vertex_circulator first = 
    v->incident_halfedges();
  typename Planar_map::Halfedge_around_vertex_circulator curr = first;
	
  do {
    if ( traits->point_is_left(curr->source()->point(), v->point())) 
      {
	if (lowest_left == pm->halfedges_end())
	  lowest_left = curr;
			
	if (traits->curve_compare_at_x_left(curr->curve(),
					    lowest_left->curve(), 
					    v->point())==SMALLER)
	  lowest_left = curr;
      }
		
    if ( traits->point_is_right(curr->source()->point(), 
				v->point()) ) 
      {
	if (lowest_right == pm->halfedges_end())
	  lowest_right = curr;
			
	if (traits->curve_compare_at_x_right(curr->curve(),
					     lowest_right->curve(), 
					     v->point())==LARGER
	    )
	  lowest_right = curr;
      }
		
		
		
    if (traits->curve_is_vertical(curr->curve())) {
      if (traits->compare_y(v->point(),
			    curr->source()->point())==LARGER)
	//debug
	//{ std::cout << "vertical up = " << curr->curve() << std::endl;
				
	vertical_up=curr; 
			
      //}//enddebug
			
      if (traits->compare_y(v->point(),
			    curr->source()->point())==SMALLER)
	//debug
	//{ std::cout << "vertical down = " << curr->curve() << std::endl;
				
	vertical_down=curr;
      //}//enddebug
			
    }        
        
  } while (++curr != first);
	
	
  /*
	vertical_down  
	|
	v   <- lowest_right      
	'v'  
	lowest_left->  ^ 
	|
	vertical_up
	*/
	
  if (!highest) {
    if (lowest_left!= pm->halfedges_end()) 
      return lowest_left;
    else 
      if (vertical_down!= pm->halfedges_end()) 
	return vertical_down;
      else
	return lowest_right;
  }
  else { //down
    if (lowest_right!=pm->halfedges_end()) 
      return lowest_right;
    else 
      if (vertical_up!=pm->halfedges_end()) 
	return vertical_up;
      else
	return lowest_left;
  }
	
}

CGAL_END_NAMESPACE

#endif // CGAL_PM_NAIVE_POINT_LOCATION_C
