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
// author(s)     : Iddo Hanniel <hanniel@post.tau.ac.il>
//                 Oren Nechushtan <theoren@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_PM_NAIVE_POINT_LOCATION_C
#define CGAL_PM_NAIVE_POINT_LOCATION_C

#include <CGAL/Pm_naive_point_location.h>

CGAL_BEGIN_NAMESPACE

//IMPLEMENTATION
//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
template <class Planar_map>
typename Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::locate(const Point & p, 
                                            Locate_type & lt) const
{
  typename Planar_map::Vertex_iterator vit=pm->vertices_begin();
  for (; vit != pm->vertices_end(); ++vit) {
    if (traits->point_is_same(p,vit->point()) ) {
      lt = Planar_map::VERTEX; 
      Halfedge_handle h(vit->incident_halfedges());
      return h;
    }
  }

  typename Planar_map::Halfedge_iterator hit=pm->halfedges_begin();
  for (; hit != pm->halfedges_end(); ++hit) {
    if (traits->curve_get_point_status(hit->curve(),p)==Traits::ON_CURVE) {
      lt = Planar_map::EDGE; 
      return hit;
    }
  }

  lt = Planar_map::UNBOUNDED_FACE;
  Locate_type temp;
  Halfedge_handle h = vertical_ray_shoot(p,temp,true);
  if (temp != Planar_map::UNBOUNDED_FACE) {
    if (temp == Planar_map::VERTEX) {  
      //since h points at the vertex and is the first 
      //halfedge after the ray clockwise! then the  face
      //is to its _right_ (maybe the specs will change in the future) 
      h = h->twin();        
    }        

    if (!(h->face()->is_unbounded())) lt = Planar_map::FACE;
    return h;
  }

  //==the vertical ray shoot returned the halfedges_end() iterator.
  if (pm->unbounded_face()->holes_begin() == 
      pm->unbounded_face()->holes_end() ) //an empty map
    return h; //return halfedges_end()

  //- returns a halfedge on an inner ccb of the unbounded face
  typename Planar_map::Holes_iterator hot =
    pm->unbounded_face()->holes_begin();
  return (*hot);
}

template <class Planar_map>
typename Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::locate(const Point & p, Locate_type & lt)
{
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->locate(p,lt);
  if (!((Bounding_box*)get_bounding_box())->locate(p,lt,h))
    h = ((cPLp)this)->locate(p,lt);
  return h;
}

template <class Planar_map>
typename Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::
vertical_ray_shoot(const Point & p, Locate_type & lt, bool up) const
{
  typename Traits::Curve_point_status point_above_under, r;
  int curve_above_under;

  lt = Planar_map::EDGE;

  // set the flags for comparison acording to the ray 
  // direction (up/down)
  if (up) {
    point_above_under = Traits::UNDER_CURVE;
    curve_above_under = LARGER;
  } else {
    point_above_under = Traits::ABOVE_CURVE;
    curve_above_under = SMALLER;
  }

  typename Planar_map::Halfedge_iterator it  = pm->halfedges_begin(),
                                         eit = pm->halfedges_end(),
                                         closest_edge = eit;
  bool first = false;
  // For each halfedge
  while (it != eit) {
    // Find if p is in the x-range of the curve and above or below it
    // according to the direction of the shoot.
    r = traits->curve_get_point_status(it->curve(), p);
    if (r == point_above_under) {
      // If the first curve in the x-range was not found yet
      if (!first) {
        closest_edge = it;
        first = true;
      } else {
        // We found another curve in the x-range and we want to remember
        // the closest
        if ( traits->curve_compare_at_x(closest_edge->curve(),
                                        it->curve(), p) == curve_above_under) 
        {
          closest_edge = it;
        }
      }
    }
    if (( r == Traits::ON_CURVE ) && (traits->curve_is_vertical(it->curve())))
    {
      // The vertical ray shoot is not including p itself,
      // thus we are interested only in vertical curves that
      // extend upwards (downwards, resp.)
      // In this case the Locate type is always EDGE
      // Remark: This treatment was originally written in the walk PL.
      //
      if (up && 
	  traits->point_is_right_top
	      (traits->curve_righttop_most(it->curve()), p) ||
	  ! up &&
	  traits->point_is_left_low
              (traits->curve_leftlow_most(it->curve()), p))
	/*
	  x       x
	  |       |
	  p=x  or  p
	  |
	  x
	*/
      {
        lt = Planar_map::EDGE;
        if (up == traits->point_is_left_low(it->target()->point(),
                                            it->source()->point()))
          return it;
        else 
          return it->twin();
      }
    }
    ++it; ++it;
  }
  
  // if we didn't find any edge above p then it is the empty face
  if ( ! first) {
    lt = Planar_map::UNBOUNDED_FACE;
    Halfedge_handle h = pm->halfedges_end();
    return h; //==NULL
  }

  // if the closest point is a vertex then find the first clockwise 
  // edge from the vertical segment
  typename Planar_map::Vertex_handle v = pm->vertices_end();
  bool maybe_vertical = false; // BUG fix (Oren)
  if (traits->point_is_same_x(closest_edge->target()->point(), p)) {
    v = closest_edge->target();
    maybe_vertical=true; // BUG fix (Oren)
  }

  if ( traits->point_is_same_x( closest_edge->source()->point(), p) ) 
    {
      if (!maybe_vertical || 
	  traits->point_is_right_top(closest_edge->target()->point(),
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
  if (v != pm->vertices_end()) {
    lt = Planar_map::VERTEX;
    closest_edge = (up) ? find_lowest(v, false) : find_lowest(v, true);
  }

  if (lt == Planar_map::VERTEX) return closest_edge;
  if (up) {
    // return the edge that is going from right to left
    // such that p is to the left of this edge

    return (traits->point_is_right(closest_edge->source()->point(),
                                   closest_edge->target()->point())) ?
      closest_edge : closest_edge->twin();
  } 
  return (traits->point_is_left(closest_edge->source()->point(),
                                closest_edge->target()->point())) ?
    closest_edge : closest_edge->twin();
}

template <class Planar_map>
typename Pm_naive_point_location<Planar_map>::Halfedge_handle
Pm_naive_point_location<Planar_map>::vertical_ray_shoot(const Point & p, 
                                                        Locate_type & lt, 
                                                        bool up)
{
  /* Make sure the source point is in the bounding box on the output */
  ((Bounding_box*)get_bounding_box())->insert(p);
  Halfedge_handle h=((cPLp)this)->vertical_ray_shoot(p,lt,up);
  /* Apply the bounding box on the output */
  if (!((Bounding_box*)get_bounding_box())->vertical_ray_shoot(p,lt,up,h)) {
    h = ((cPLp)this)->vertical_ray_shoot(p,lt,up);
    CGAL_assertion(lt!=Planar_map::UNBOUNDED_FACE);
  }
  return h;
}


//find the first halfedge pointing at v, when going clockwise
//if highest==true - start from 12 oclock, else start from 6 oclock
template <class Planar_map>
typename Pm_naive_point_location<Planar_map>::Halfedge_handle 
Pm_naive_point_location<Planar_map>::
find_lowest(typename Pm_naive_point_location<Planar_map>::Vertex_handle v,
            bool highest) const
{
  Halfedge_handle lowest_left = pm->halfedges_end();
  Halfedge_handle lowest_right = pm->halfedges_end();
  Halfedge_handle vertical_up = pm->halfedges_end();
  Halfedge_handle vertical_down = pm->halfedges_end();

  typename Planar_map::Halfedge_around_vertex_circulator first = 
    v->incident_halfedges();
  typename Planar_map::Halfedge_around_vertex_circulator curr = first;

  do {
    if (traits->point_is_left(curr->source()->point(), v->point())) {
      if (lowest_left == pm->halfedges_end())
        lowest_left = curr;
      else if (traits->curve_compare_at_x_left(curr->curve(),
                                               lowest_left->curve(), 
                                               v->point()) == SMALLER)
        lowest_left = curr;
    }

    if (traits->point_is_right(curr->source()->point(), v->point())) {
      if (lowest_right == pm->halfedges_end())
        lowest_right = curr;
      else if (traits->curve_compare_at_x_right(curr->curve(),
                                                lowest_right->curve(), 
                                                v->point()) == LARGER)
        lowest_right = curr;
    }

    if (traits->curve_is_vertical(curr->curve())) 
    {
      if (traits->compare_xy(v->point(),
			     curr->source()->point()) == LARGER)
	vertical_up=curr; 
			
      if (traits->compare_xy(v->point(),
			     curr->source()->point()) == SMALLER)
	vertical_down=curr;			
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
    if (lowest_left != pm->halfedges_end()) return lowest_left;
    if (vertical_down != pm->halfedges_end()) return vertical_down;
    return lowest_right;
  }
  if (lowest_right != pm->halfedges_end()) return lowest_right;
  if (vertical_up != pm->halfedges_end()) return vertical_up;
  return lowest_left;
}

CGAL_END_NAMESPACE

#endif
