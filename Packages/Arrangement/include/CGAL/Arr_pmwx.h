// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.2-I-13 $
// release_date  : $CGAL_Date: 2000/04/14 $
//
// file          : include/CGAL/Arr_pmwx.h
// package       : arr (1.24)
// maintainer    : Sigal Raab <raab@math.tau.ac.il>
// author(s)     : Iddo Hanniel <hanniel@math.tau.ac.il>
// coordinator   : Dan Halperin <halperin@math.tau.ac.il>
//
// ======================================================================
#ifndef CGAL_ARR_PMWX_H
#define CGAL_ARR_PMWX_H

#include <CGAL/Pm_with_intersections_misc.h>

CGAL_BEGIN_NAMESPACE

template<class Planar_map_>
class Arr_pmwx : public Planar_map_
{
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  typedef typename Planar_map::Traits_wrap Pm_traits_wrap;
  typedef typename Planar_map::Halfedge_handle Halfedge_handle;
  typedef typename Planar_map::Vertex_handle Vertex_handle;
  typedef typename Planar_map::Face_handle Face_handle;
  typedef typename Planar_map::Locate_type Locate_type;
  typedef typename Planar_map::Halfedge_around_vertex_circulator 
  Halfedge_around_vertex_circulator;
  typedef typename Planar_map::Ccb_halfedge_circulator Ccb_halfedge_circulator;
  typedef typename Planar_map::Holes_iterator Holes_iterator;
  typedef typename Traits::X_curve X_curve;
  typedef typename Traits::Point Point;
  typedef Pm_change_notification<Planar_map> Pmwx_change_notification; 

  Arr_pmwx() : Planar_map()
  { 
  }    
	
  Arr_pmwx (Pm_point_location_base<Planar_map> *pl_ptr) : Planar_map(pl_ptr)
  {
  }

  Arr_pmwx (const Traits& tr_, Pm_point_location_base<Planar_map> *pl_ptr) 
    : Planar_map(tr_, pl_ptr, NULL)
  {
  }

  Arr_pmwx(Pm_traits_wrap *tr_ptr, Pm_point_location_base<Planar_map> *pl_ptr) 
    : Planar_map(tr_ptr, pl_ptr, NULL)
  {
  }
	
  ~Arr_pmwx()
  {
  }

  Halfedge_handle insert_from_vertex(const X_curve& c, 
				     Vertex_handle v1, 
				     Pmwx_change_notification *en = NULL)
  {
    return insert(c, en);
  }

  Halfedge_handle insert(const X_curve& c, Pmwx_change_notification *en = NULL)
  {
    Halfedge_handle last_edge;
    //step 0 :ensure that cv source is left of target
    X_curve cv;
    Point last_point;
    bool left_to_right=true;
		
    //debug
    //std::cout << "source=" << traits->curve_source(c).to_point();
    //std::cout << "target=" << traits->curve_target(c).to_point() << endl;
    //std::cout << "left_low(t,s)=" << traits->point_is_left_low
    // (traits->curve_target(c),traits->curve_source(c)) << endl;
		
    if ( traits->point_is_left_low( traits->curve_target(c),
				    traits->curve_source(c)))
      {
	//flip if target is left/ if vertical target is lower
	cv=traits->curve_flip(c);
	left_to_right=false;
      }
    else
      {
	cv=c;
      }
		
    X_curve cvo=cv; //the original curve before intersections
    last_point = traits->curve_target(cvo);

    Locate_type lt;
    Halfedge_handle h=locate(traits->curve_source(cv),lt);
		
    Vertex_handle srcv,trgtv;
		
    //step 1: case 3: source(cv) is in face - find face, break
    //        case 1: source(cv) is on edge - split edge, find next face
    //        case 2: source(cv) is on vertex - find next face
    Face_handle next_face;
		
    switch (lt) 
      {
      case EDGE : {
		  
	// bug fix (Sigal):
	// use the halfedge that has the same direction as the underlying curve
	if (!traits->point_is_same(h->source()->point(),
				   traits->curve_source(h->curve()))) {
	  h = h->twin();
	}

	X_curve cv1=h->curve();
	X_curve tmp1,tmp2;

	//split curve at splitting point :	
	traits->curve_split(cv1,tmp1,tmp2, traits->curve_source(cv)); 
			
	Halfedge_handle sp = Planar_map::split_edge(h, tmp1,tmp2);
	if (en != NULL) en->split_edge(sp, sp->next_halfedge(), tmp1, tmp2);
	h= sp;

	//workaround since MSVC does not approve the automatic cast here
	lt=VERTEX; 
			
	// - we are now at the state of vertex ! (no need to break)
      }
			
      case VERTEX : {
	srcv=h->target();
	//finding the next face the curve goes through
	Halfedge_around_vertex_circulator 
	  strt,circ,circ1=srcv->incident_halfedges();
	strt=circ=circ1;
	++circ1;
			
	if (traits->curves_overlap(cv,circ->curve())) 
	  {
	    //there is an overlap, don't get into the curve_is_between_cw loop.
	    next_face=circ->face();
	    break;
	  }      
			
	if (circ!=circ1) 
	  {
	    do {
	      //check overlap
	      if (traits->curves_overlap(cv,circ1->curve())) 
		{
		  // there is an overlap, 
		  // don't get into the curve_is_between_cw func.
		  next_face=circ1->face();
		  break;        
		}
					
	      if (traits->curve_is_between_cw
		  (cv, circ->curve(), circ1->curve(),traits->curve_source(cv)))
		{
		  next_face=circ->face();
		  break;
		}
	      ++circ;
	      ++circ1;
	    } while (circ!=strt);
				
	  }
	else 
	  { //only one edge from the vertex (tip)
	    next_face=circ->face();
	  }
			
	break;
      }    
			
      case FACE : {
	next_face=h->face();
	break;
      }
		
      case UNBOUNDED_FACE: {
	next_face=unbounded_face();
	break;
      }
			
      default:
	CGAL_assertion(false);
			
      } //switch(lt)
		
    Point p1,p2;
    //maybe do it with some other ctr (Zurich decisions)
    X_curve cv1,cv2;
    Halfedge_handle hh;
		
    while (1) 
      {
	//step 2: finding the closest intersection point, p, 
	//with the faces boundary
	//if no intersection exists - insert cv (into arrangement and pm) 
	//and break
			
	//change - to prevent cascaded intersection
	//if (!(find_nearest_intersection(cv,traits->curve_source(cv),
	//next_face,p1,p2,hh)) ) 
	if (!(find_nearest_intersection(cvo,traits->curve_source(cv),
					next_face,p1,p2,hh, en)) ) 
	  {
	    if (lt==FACE || lt==UNBOUNDED_FACE) 
	      {
		Halfedge_handle e;
		if (left_to_right) 
		  {
		    e=Planar_map::insert_in_face_interior(cv,next_face);
		  }
		else 
		  {
		    e=Planar_map::insert_in_face_interior
		      (traits->curve_flip(cv),next_face);
		  }
		if (en != NULL) en->add_edge(cv,e,left_to_right);
		last_edge = e;
		if (!traits->point_is_same
		    (last_edge->target()->point(), last_point))
		  last_edge = last_edge->twin();
		return last_edge; 
	      }
	    else 
	      { //source is on a vertex and target inside next_face 
		Halfedge_handle e;
		if (left_to_right) 
		  {
		    e=Planar_map::insert_from_vertex(cv,srcv,true); 
		  }
		else 
		  {
		    e=Planar_map::insert_from_vertex
		      (traits->curve_flip(cv),srcv,false);
		  }
		if (en != NULL) en->add_edge(cv,e,left_to_right);
		last_edge = e;
		if (!traits->point_is_same(last_edge->target()->point(), 
					   last_point))
		  last_edge = last_edge->twin();
		return last_edge; 
	      }
	  } //if (!find_nearest_intersection...
		
	// bug fix (Sigal):
	// use the halfedge that has the same direction as the underlying curve
	if (! traits->point_is_same(hh->source()->point(),
				    traits->curve_source(hh->curve()))) {
	  hh = hh->twin();
	}
	
			
			
	//DEALING WITH OVERLAP:
	//1. there are 2 kinds of overlap - 
	//   a - p1 == srcv                  p1------p2
	//                                 srcv------------
			
	//   b - p1 is right of srcv         p1------p2
	//                             srcv--------------------
			
	//2. we can seperate the cases, we can make a block of steps 3,4,5 for 
	//   overlap/no overlap, or we can deal only with a, and b can be dealt
	//   by splitting at p1 (and in the step after it becomes 
	//   srcv == case a).
	//3. We take the latter way (easier to implement) , 
	//   might be changed in future!!
			
	//step 3: if p1,p2 is on edge - split the edge and find the 
	//        vertex of p.
	// if p is on a vertex find it
			
	if (traits->point_is_same(p1, hh->source()->point())) 
	  {
	    if (traits->point_is_same(p1,p2)) 
	      { //no overlap
		trgtv=hh->source();
	      }
	    else 
	      { //deal with overlap
		//we need to check first if lt==FACE because in that case 
		//srcv doesn't hold a point yet (it is null) - otherwise 
		//srcv has a value
		if ( (lt==FACE || lt==UNBOUNDED_FACE) || 
		     !traits->point_is_same(p1,srcv->point())) 
		  {
		    //deal with it as if the curve ends at p1. (the segment 
		    //(p1,p2) will be dealt with at the next iteration)
		    trgtv=hh->source();
		  }
		else 
		  { //p1 == srcv
		    if (traits->point_is_same(p2,hh->target()->point())) 
		      {
			trgtv=hh->target();
			if (en != NULL) en->add_edge(hh->curve(),hh,
						     left_to_right,true);
			last_edge = hh;
		      }
		    else 
		      { //p2 is inside hh
			X_curve crv=hh->curve();
			X_curve crv1,crv2;
			traits->curve_split(crv,crv1,crv2,p2);
			Halfedge_handle e=Planar_map::split_edge(hh,crv1,crv2);
			if (en != NULL) en->split_edge(e, e->next_halfedge(), 
						       crv1, crv2);
			trgtv=e->target();
							
			//inserting the curve (crv1/2) that is incident to p1
			if (traits->point_is_same(traits->curve_source(crv1),
						  p1) ||
			    traits->point_is_same(traits->curve_target(crv1),
						  p1)) 
			  {
			    if (en != NULL) en->add_edge(crv1,e,left_to_right,
							 true);
			  }
			else 
			  {
			    if (en != NULL) en->add_edge(crv2,e,left_to_right,
							 true);
			  }
			last_edge = e;
		      }
		  }
	      }
	  }
	else
	  {
	    if (traits->point_is_same(p1, hh->target()->point())) 
	      {
		if (traits->point_is_same(p1,p2)) 
		  { //no overlap
		    trgtv=hh->target();
		  }
		else 
		  { //deal with overlap
		    //see comment above why this condition is here
		    if ( (lt==FACE || lt==UNBOUNDED_FACE) || 
			 !traits->point_is_same(p1,srcv->point())) 
		      {
			//deal with it as if the curve ends at p1. 
			//(the segment (p1,p2) will be dealt with at the 
			//next iteration)
			trgtv=hh->target();
		      } 
		    else 
		      {
			if (traits->point_is_same(p2,hh->source()->point())) 
			  {
			    trgtv=hh->source();
			    if (en != NULL) en->add_edge(hh->curve(),hh,
							 left_to_right,true);
			    last_edge = hh;
			  }
			else 
			  { //p2 is inside hh
			    X_curve crv=hh->curve();
			    X_curve crv1,crv2;
			    traits->curve_split(crv,crv1,crv2,p2);
			    Halfedge_handle e=Planar_map::split_edge(hh,crv1,
								     crv2);
			    if (en != NULL) 
			      en->split_edge(e, e->next_halfedge(), 
					     crv1, crv2);
			    trgtv=e->target();
			    //inserting the curve (crv1/2) that is incident 
			    //to p1 
			    //bug fix (shai and iddo 17/3/00)
			    //e can be non-incident to p1 (if cv is right to 
			    //left) make it incident
			    if ( ! traits->point_is_same(e->source()->point()
							 ,p1) &&
				 ! traits->point_is_same(e->target()->point(),
							 p1))
			      e = e->next_halfedge(); // end bug fix

			    if (traits->point_is_same(traits->curve_source
						      (crv1),p1) ||
				traits->point_is_same(traits->curve_target
						      (crv1),p1)) 
			      {
				if (en != NULL) 
				  en->add_edge(crv1,e,left_to_right,true);
			      }
			    else 
			      {
				if (en != NULL) 
				  en->add_edge(crv2,e,left_to_right,true);
			      }
			    last_edge = e;
			  }
		      }
		  }
	      }
	    else 
	      { //p1 is not on a vertex 
		if (traits->point_is_same(p1,p2)) 
		  { //no overlap
		    X_curve crv=hh->curve();
		    //debug
		    //std::cout << traits->curve_source(crv).to_point() 
		    //<< " , " ;
		    //std::cout << traits->curve_target(crv).to_point() 
		    //<< std::endl ;
						
		    X_curve crv1,crv2;
		    //maybe check if it's fromleft to right in order to 
		    //assure left to right curve ???
						
		    traits->curve_split(crv,crv1,crv2,p2);
		    Halfedge_handle he = Planar_map::split_edge(hh,crv1,crv2);
		    if (en != NULL) en->split_edge(he, he->next_halfedge(), 
						   crv1, crv2);
		    trgtv=he->target();
		  }
		else 
		  { //deal with overlap
		    //will never get here since we took care in the previous 
		    //step that
		    //overlap will take place only when p1 is on a vertex 
		    //(otherwise we 
		    //create a vertex at p1 (split the edge) and in the next 
		    //iteration 
		    //p1 is on a vertex.)
		    CGAL_assertion(false);
		  }
	      }
	  }
	//step 4: if (p2==cv.target)
	//          insert cv and return. (p2 is the end of cv)
	if (traits->point_is_same(p2,traits->curve_target(cv))) 
	  {
	    if (lt==FACE || lt==UNBOUNDED_FACE) 
	      {
		// bug fix (sigal) : go on only if it is not 
		//  an overlap
		if ( traits->point_is_same(p1,p2)) 
		  { 
		    //can't get here if there is an overlap.
		    Halfedge_handle e;
		    if (left_to_right) 
		      {
			//.current_iterator() :
			e=Planar_map::insert_from_vertex(cv,trgtv,false); 
		      }
		    else 
		      {
			//.current_iterator() :
			e=Planar_map::insert_from_vertex
			  (traits->curve_flip(cv),trgtv,true);
		      }
					
		    if (en != NULL) en->add_edge(cv,e,left_to_right);
		    last_edge = e;
		    if (!traits->point_is_same(last_edge->target()->point(), 
					       last_point))
		      last_edge = last_edge->twin();
		    return last_edge; 
		  }
	      }
	    else 
	      { //cv source is on a vertex or edge (that was split at step 3
		// or at the switch statement) 
					
		if (traits->point_is_same(p1,p2)) { //no overlap
		  Halfedge_handle e;
		  if (left_to_right) 
		    {
		      //.current_iterator(), current_iterator() :
		      e=Planar_map::insert_at_vertices(cv,srcv,trgtv); 
		    }
		  else 
		    {
		      //.current_iterator(), current_iterator() :
		      e=Planar_map::insert_at_vertices(traits->curve_flip(cv),
						       srcv,trgtv);
		    }
						
		  if (en != NULL) en->add_edge(cv,e,left_to_right);
		  last_edge = e;
		  if (!traits->point_is_same(last_edge->target()->point(), 
					     last_point))
		    last_edge = last_edge->twin();
		  return last_edge; 
		}
		else 
		  {
		    //do nothing - it was dealt with
		    if (traits->point_is_same(p1,srcv->point())) 
		      {
			//wev'e alraedy pushed (p1,p2) into edge list in step 3
			if (!traits->point_is_same
			    (last_edge->target()->point(), last_point))
			  last_edge = last_edge->twin();
			return last_edge; 
		      }
		    else 
		      {
			//do nothing p2 will be dealt with at next iteration
			//we are still at p1
		      }
		  }
	      }
	  }
			
	//step 5: split cv at trgtv->point() into cv1 and cv2 and insert cv1.
	traits->curve_split(cv,cv1,cv2,trgtv->point());
			
	if (lt==FACE || lt==UNBOUNDED_FACE) 
	  {
	    Halfedge_handle e;
	    if (left_to_right) 
	      {
		//.current_iterator():
		e=Planar_map::insert_from_vertex(cv1,trgtv,false);
	      }
	    else 
	      {
		//.current_iterator() :
		e=Planar_map::insert_from_vertex(traits->curve_flip(cv1),
						 trgtv,true); 
	      }
				
	    if (en != NULL) en->add_edge(cv1,e,left_to_right);
	    last_edge = e;
	  }
	else 
	  {
	    if (traits->point_is_same(p1,p2)) 
	      { //no overlap
		Halfedge_handle e;
		if (left_to_right) 
		  {
		    //debug
		    //std::cout << "srcv=" << srcv->point().to_point();
		    //std::cout << " trgv=" << trgtv->point().to_point() << 
		    //std::endl;

		    //.current_iterator() .current_iterator():
		    e=Planar_map::insert_at_vertices(cv1,srcv,trgtv); 
		  }
		else 
		  {
		    //.current_iterator() .current_iterator():
		    e=Planar_map::insert_at_vertices(traits->curve_flip(cv1),
						     srcv,trgtv); 
		  }
					
		if (en != NULL) en->add_edge(cv1,e,left_to_right);
		last_edge = e;
	      }
	    else 
	      {//deal with overlap (don't insert into pm)
		if (!traits->point_is_same(p1,srcv->point())) 
		  {
		    Halfedge_handle e;
		    if (left_to_right) 
		      {
			//.current_iterator() .current_iterator():
			e=Planar_map::insert_at_vertices(cv1,srcv,trgtv); 
		      }
		    else 
		      {
			//.current_iterator() .current_iterator():
			e=Planar_map::insert_at_vertices
			  (traits->curve_flip(cv1),srcv,trgtv);
		      }
						
		    if (en != NULL) en->add_edge(cv1,e,left_to_right);
		    last_edge = e;
		  }            
		else 
		  { //p1 == srcv->point(), was dealt with at step 3 
		    //so do nothing
		  }
	      }
	  }
			
	//step 6: find next face and goto step 2
	srcv=trgtv;
	cv=cv2;
			
	//finding the next face the curve goes through
	Halfedge_around_vertex_circulator strt,circ,circ1 = 
	  srcv->incident_halfedges();
	strt=circ=circ1;
	++circ1;
	//check overlaps
	if (traits->curves_overlap(cv,circ->curve())) 
	  {
	    //there is an overlap, don't get into the curve_is_between_cw loop.
	    next_face=circ->face();
	  }      
	else 
	  {
				
	    if (circ!=circ1) 
	      {
		do {
		  //check overlap
		  if (traits->curves_overlap(cv,circ1->curve())) 
		    {
		      //there is an overlap, 
		      //don't get into the curve_is_between_cw func.
		      next_face=circ1->face();
		      break;        
		    }
						
		  if (traits->curve_is_between_cw(cv,circ->curve(),
						  circ1->curve(),
						  traits->curve_source(cv))) 
		    {
		      next_face=circ->face();
		      break;
		    }
		  ++circ;
		  ++circ1;
		} while (circ!=strt);
	      }
	    else
	      {
		next_face=circ->face();
	      }
				
	  }
			
	lt = VERTEX;
			
      } //while(1)
  }


  ///////////////////////////////////////////////////////////////////
  //           PRIVATE IMPLEMENTATIONS
  ////////////////////////////////////
private:

  //returns true iff there is an intersection strictly to the right of r_pt
  //with an edge of the face f (in that case the nearest intersecting 
  //halfedge will be heald at hh and the intersection points in p1,p2)

  //in the future we will move this function to the pl strategy, so that the
  //default pl can make use of its trapezoidal map to find it faster!


  bool find_closer_intersection(const X_curve& new_cv,
				const Halfedge_handle& halfedge,
				const Point& r_pt, 
				Point& p1, 
				Point& p2,
				Halfedge_handle& h, 
				Pmwx_change_notification *en) const 
  {
    Point p_cand1, p_cand2;
    bool found = false;
    const X_curve& halfedge_cv = halfedge->curve(); 
    bool have_support_cv = false;
    if (en != NULL)
      have_support_cv = en->have_support_curve();

    if (have_support_cv)
      {
	// if we are here then en must be non-NULL
	const X_curve &support_cv = en->edge_support_curve(halfedge);
		
	// this line improves the running time for leda_exact_segemnt
	// changing halfedge_cv to support_cv results in poor results but 
	//might be good in other cases???
	if (!traits->do_intersect_to_right(new_cv, halfedge_cv, r_pt))
	  return false;
		
	if (traits->nearest_intersection_to_right(new_cv, support_cv, 
						  r_pt, p_cand1, p_cand2)) 
	  {
	    // since we are checking on the parent, we 
	    // should make sure that the 
	    // intersection point is on the halfedge_cv and not only 
	    // on the parent.
	    // do not worry: we will get the same intersection point 
	    // for the correct
	    // halfedge_cv as well, and therefore we can throw it away 
	    // if it's not on halfedge_cv
	    
	    // the intersection is only one point
	    if (traits->point_is_same(p_cand1, p_cand2)) 
	      {
		if (traits->curve_get_point_status(halfedge_cv, p_cand1) == 
		    Traits::ON_CURVE) 
		  {
		    found = true;
		  }
	      }
	    else // the intersection is a segment (there is an overlap)
	      {
		Point halfedge_cv_left = traits->curve_source(halfedge_cv);
		Point halfedge_cv_right = traits->curve_target(halfedge_cv);
		if ( traits->point_is_left_low( halfedge_cv_right, 
						halfedge_cv_left)) 
		  {
		    halfedge_cv_left = traits->curve_target(halfedge_cv);
		    halfedge_cv_right = traits->curve_source(halfedge_cv);
		  }
		if ( traits->point_is_left_low( p_cand1, halfedge_cv_left)) 
		  {
		    p_cand1 = halfedge_cv_left;
		  }
		if ( traits->point_is_left_low( halfedge_cv_right, p_cand2)) 
		  {
		    p_cand2 = halfedge_cv_right;
		  }
		if (traits->point_is_left_low( p_cand1, p_cand2)) 
		  {
		    found = true;
		  }
	      }
	  }
      }
    else
      {
	if (traits->nearest_intersection_to_right(new_cv, halfedge_cv, 
						  r_pt, p_cand1, p_cand2)) 
	  {
	    if (!traits->point_is_left_low(p_cand1, p_cand2)) 
	      {
		Point ptmp = p_cand1;
		p_cand1 = p_cand2;
		p_cand2 = ptmp;
	      }
	    found = true;
	  }
      }

    if ( found && traits->point_is_left_low( p_cand1, p1)) 
      {
	p1=p_cand1;
	p2=p_cand2;
	h = halfedge;
      }
    return found;
  }

  bool find_nearest_intersection(const X_curve& cv, 
				 const Point& r_pt, 
				 Face_handle f, 
				 Point& p1, 
				 Point& p2, 
				 Halfedge_handle& h, 
				 Pmwx_change_notification *en) const
  {
    Ccb_halfedge_circulator first,ccb_circ;
    bool FindFlag=false;
    //initialize p with a far value 
    p1=traits->point_to_right(traits->curve_target(cv)); 
		
    Point p_cand1,p_cand2;
    if (f->does_outer_ccb_exist()) 
      {
	ccb_circ=first=f->outer_ccb();
	do {
	  if (find_closer_intersection(cv, ccb_circ, r_pt, p1, p2, h, en)) 
	    {
	      FindFlag = true;
	    }
	} while (++ccb_circ!=first) ;
			
      }   //outer ccb
		
    for (Holes_iterator hit=f->holes_begin(); hit!=f->holes_end(); ++hit) 
      {
	first=ccb_circ=*hit;
	do {
	  if (find_closer_intersection(cv, ccb_circ, r_pt, p1, p2, h, en)) 
	    {
	      FindFlag = true;
	    }
	  ++ccb_circ;
	} while (ccb_circ!=first);
      }
		
    return FindFlag;
  }

};


CGAL_END_NAMESPACE

#endif
