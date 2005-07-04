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
#ifndef CGAL_ARR_LANDMARKS_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_LANDMARKS_POINT_LOCATION_FUNCTIONS_H

/*! \file
 * Member-function definitions for the Arr_landmarks_point_location<Arrangement>
 * class.
 */

#ifdef LANDMARKS_CLOCK_DEBUG
	#define LM_CLOCK_DEBUG(cmd)   cmd
#else
	#define LM_CLOCK_DEBUG(cmd)
#endif

#define CGAL_LM_DEBUG

#ifdef CGAL_LM_DEBUG
	#define PRINT_DEBUG(expr)   std::cout << expr << std::endl
	#define LM_DEBUG(cmd)   cmd
#else
	#define PRINT_DEBUG(expr)
	#define LM_DEBUG(cmd) 
#endif

//#define CGAL_LM_DEBUG_RANDOM

#ifdef CGAL_LM_DEBUG_RANDOM
	#define PRINT_DEBUG_R(expr)   std::cout << expr << std::endl
	#define LM_DEBUG_R(cmd)   cmd
#else
	#define PRINT_DEBUG_R(expr)
	#define LM_DEBUG_R(cmd) 
#endif 

#define PRINT_ERROR(expr)   std::cerr << expr << std::endl


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement_2, class Arr_landmarks_generator>
Object Arr_landmarks_point_location<Arrangement_2,Arr_landmarks_generator>
::locate (const Point_2& p) const
{
	PRINT_DEBUG("------ locate point "<< p) ;

	//if this is an empty map - return the unbounded face
	if (p_arr->number_of_vertices() == 0) 
		return make_object (p_arr->unbounded_face());

	Object  landmark_location_obj; 
	Point_2 landmark_point = lm_gen->get_closest_landmark (p, 
								landmark_location_obj);
	//NN_Point_2 nearest_landmark = lm_gen->find_closest_landmark(p);
	//Point_2 landmark_point = nearest_landmarks.get_point();
	//Object  landmark_location_obj = nearest_landmarks.get_obj() ;
	PRINT_DEBUG("nearest neighbor of point "<< p << " is " << landmark_point);

	//walk from the nearest_vertex to the point p, using walk algorithm, 
	//and find the location of p.	
	LM_CLOCK_DEBUG(clock_t w_time_start = clock());

	Object	out_obj; //the output object

	//if the landmark s not found in the arangement
	Vertex_const_handle     v;
	Halfedge_const_handle   h;
	Face_const_handle       f;

	if (landmark_location_obj.is_empty())
	{
		PRINT_ERROR( "landmark_location_obj is empty" );
		CGAL_assertion (false);
		return out_obj;
	}
	else if (assign (v, landmark_location_obj))
	{
		PRINT_DEBUG( "landmark_location_obj is a vertex: "<< v->point());
		out_obj = _walk_from_vertex( v , p);
	}
	else if (assign (f, landmark_location_obj))
	{
		PRINT_DEBUG( "landmark_location_obj is a face. ");
		out_obj = _walk_from_face( f, p, landmark_point);
	}
	else if (assign (h, landmark_location_obj))
	{
		PRINT_DEBUG( "landmark_location_obj is a halfedge: "<< h->curve());
		out_obj = _walk_from_edge( h, p, landmark_point);
	}
	else 
	{
		PRINT_ERROR( "unknown object");
		CGAL_assertion (false);
		return out_obj;
	}

	LM_CLOCK_DEBUG(clock_t w_time_end = clock() );
	LM_CLOCK_DEBUG(clock_for_walk += (double) (w_time_end - w_time_start) );
	PRINT_DEBUG( "return from walk" << std::endl);

#ifdef CGAL_LM_DEBUG
	if (out_obj.is_empty())
	{
		PRINT_ERROR( "object is empty" );
		CGAL_assertion (false);
	}
	else if (assign (h, out_obj))
	{
		PRINT_DEBUG( "object is a halfedge: "<< h->curve());
	}
	else if (assign (v, out_obj))
	{
		PRINT_DEBUG( "object is a vertex: "<< v->point());
	}
	else if (assign (f, out_obj))
	{
		PRINT_DEBUG( "object is a face. ");
	}
#endif

  if (assign (f, out_obj))
	{
    // If we reached here, we did not locate the query point in any of the holes
    // inside the current face, so we conclude it is contained in this face.
    // However, we first have to check whether the query point coincides with
    // any of the isolated vertices contained inside this face.
    Isolated_vertices_const_iterator   iso_verts_it;
    typename Traits_wrapper_2::Equal_2  equal = traits->equal_2_object();

    for (iso_verts_it = f->isolated_vertices_begin();
        iso_verts_it != f->isolated_vertices_end(); ++iso_verts_it)
    {
      if (equal (p, (*iso_verts_it).point()))
        return (make_object (iso_verts_it));
    }		
	}

	return (out_obj) ;
}

//----------------------------------------------------
// walks from the vertex to the point
// \param nearest_vertex - (input) the closest vertex to point p
// \param p - (input) the point to locate.
//
template <class Arrangement_2, class Arr_landmarks_generator>
Object Arr_landmarks_point_location<Arrangement_2,Arr_landmarks_generator>
::_walk_from_vertex(Vertex_const_handle nearest_vertex,
					const Point_2 & p)   const
{ 
	PRINT_DEBUG("inside walk_from_vertex. p= "<< p	<< 
				 ", nearest_vertex = "<<nearest_vertex->point() );

	//inits
	bool new_vertex = false;
	Vertex_const_handle vh = nearest_vertex;
	Object obj;

	Vertex_const_handle     v;
	Halfedge_const_handle   h;
	Face_const_handle       f;

  if (vh->is_isolated())
  {
    f = p_arr->incident_face(vh);
    return _walk_from_face(f, p, vh->point());
  }

	//find face
	do {
		//find the edge out_edge which is the best possibly 
		//pointing to the face containing p

		flipped_edges.clear();  //clear the curves that were flipped

    new_vertex = false;
		LM_CLOCK_DEBUG(clock_t ff_time_start = clock());
		obj = _find_face (p, vh, new_vertex);
		LM_CLOCK_DEBUG(clock_t ff_time_end = clock() );
		LM_CLOCK_DEBUG(clock_ff += (double) (ff_time_end - ff_time_start) );

		if (new_vertex) {
			PRINT_DEBUG( "NEW vertex 1 " );
			//check if the new vertex is really closer 
			// I removed the check if the vertex is closer since there is no 
			// compare distance 
			//if (traits->compare_distance(p, out_vertex->point(), vh->point())  
			//	!= SMALLER) {PRINT_DEBUG("Error 2: new vertex"); return; } 
			if (assign (v, obj))
			{
				vh = v;
			}
			else
			{
				CGAL_assertion (false);	
				return Object();
			}
		}
		else if (obj.is_empty())
		{
			PRINT_ERROR( "object is empty" );
			CGAL_assertion (false);
			return obj;
		}
		else if (assign (h, obj))
		{
			PRINT_DEBUG( "_find_face found a halfedge: "<< h->curve());
			return (obj);
		}
		else if (assign (v, obj))
		{
				PRINT_DEBUG( "_find_face found a vertex: "<< v->point());
				return (obj);
		}
		else if (assign (f, obj))
		{
			PRINT_DEBUG("face that was a face ");
			return _walk_from_face (f, p, vh->point());
		}

	} while (new_vertex);	
			
	CGAL_assertion (false);
	return Object();
}

//----------------------------------------------------
// Return the edge around vh that is before cw to the segment vh->p
// 
//\param p - the input point.
//\param vh - the input vertex
//\param found_vertex_or_edge - output bool, 
//                 if the point was found on a vertex or halfedge
//\param new_vertex - output bool, if a closer vertex to p was found
//
template <class Arrangement_2, class Arr_landmarks_generator>
Object Arr_landmarks_point_location<Arrangement_2,Arr_landmarks_generator>
::_find_face (const Point_2 & p, 
			  Vertex_const_handle vh,
			  bool & new_vertex 
        ) const
{ 
	LM_CLOCK_DEBUG( entries_to_find_face++ );
	PRINT_DEBUG("inside find_face. p ="<< p <<" , vh = "<<vh->point() ); 

	new_vertex = false;

	// check if the point equals the vertex. 
	if (traits->equal_2_object()(vh->point(), p))
	{
		return make_object (vh);
	}

	//create a segment vh--p. 
	Point_2 v = vh->point();
	X_monotone_curve_2 seg(v, p);

	//get halfedges around vh
  CGAL_assertion (!vh->is_isolated());
	Halfedge_around_vertex_const_circulator circ = vh->incident_halfedges(); 
	Halfedge_around_vertex_const_circulator circ_done (circ);
	Halfedge_around_vertex_const_circulator prev = circ;

	typename Traits_wrapper_2::Compare_xy_2           compare_xy = 
                                      traits->compare_xy_2_object();

	typename Traits_wrapper_2::Compare_cw_around_point_2 compare_cw_around_point =
                                      traits->compare_cw_around_point_2_object();

	//check if cv_other_point is to the left of v,  to the right, 
	//or if the curve is vertical
	Comparison_result cv_orient = compare_xy(v, (*circ).source()->point());

	//check if p is to the left of v,  to the right, or if the segment is vertical
	Comparison_result seg_orient = compare_xy(v,p);

	//save results 
	Comparison_result res1;

	PRINT_DEBUG("seg_orient ="<< seg_orient ); 
	PRINT_DEBUG("cv_orient ="<< cv_orient << " , circ ="<< (*circ).curve());

	/////////////////////////////////// 6.12 - wrapper	
	//TODO: what if both seg and curve are verticals ? @@@@
	//what if one is vertical ?

	//curves are to different sides
	if (cv_orient != seg_orient)  
	{
		PRINT_DEBUG("seg_orient != cv_orient : "); 
		//find a curve that is on the same side as seg
		do {
			circ++;
			cv_orient = compare_xy(v,(*circ).source()->point());
		} while (seg_orient != cv_orient && circ!=circ_done);

		//if exists - go on from next "if" 
		//if not exist - find the curve that is the largest (cw) 
		//in this side, and this is the edge we're looking for
		if (seg_orient != cv_orient) 
		{
			//seg is on one side (right/left) and all curves are on the 
			//other side (left/right)
			//circ == circ_done
			do {
				prev = circ;
				circ++;
				res1 =  compare_cw_around_point((*circ).curve(), 
                                        (*prev).curve(), v);//TODO, false);
				PRINT_DEBUG("circ = " << (*circ).curve() << "  res1= " << res1 ); 
			} while (res1==LARGER && circ!=circ_done);

			//out_edge = prev;
			//found_face = true;
			PRINT_DEBUG ( "new_find_face return face = " << (*prev).curve() );
			return make_object ((*prev).face());
		}
	}

	//both curves are to the same side
	if (seg_orient == cv_orient) 
	{
		PRINT_DEBUG("seg_orient == cv_orient : "); 
		res1 = compare_cw_around_point(seg, (*circ).curve(), v);
		if (res1 == LARGER) 
		{
			//if the segment is larger than the curve cw, we will go ++ 
			//cw with the circ and find a curve that is larger than seg. 
			//then we will take the curve that was just before that. 
			PRINT_DEBUG("res1 == LARGER : "); 
			do {
				prev = circ;
				circ++;
				PRINT_DEBUG("circ++ = " << (*circ).curve() ); 
				cv_orient = compare_xy(v,(*circ).source()->point());
				if (seg_orient == cv_orient) 
				{
					res1 =  compare_cw_around_point(seg, (*circ).curve(), v);
					PRINT_DEBUG("circ = "<<(*circ).curve()<<" res1= "<<res1); 
				}
			} while (res1 == LARGER && seg_orient == cv_orient 
				     && circ!=circ_done);

			//if res1 is not larger => seg is between prev and circ,return prev
			//if seg_orient != cv_orient, then we changes side and the 
			//other side is larger than seg.
			// then also seg is between prev & circ and we have to return prev

			if (res1 == LARGER && seg_orient == cv_orient)  
			{//we 're only out the while because the circ end
				//in this case the seg is larger than ALL curves and all 
				//curves are to the same side 
				//we need to find the largest of all curves.
				PRINT_DEBUG("circ == circ_done : "); 
				do {
					prev = circ;
					circ++;
					res1 =  compare_cw_around_point((*circ).curve(), 
                                          (*prev).curve(), v);//TODO: false);
					PRINT_DEBUG("circ = "<<(*circ).curve()<<" res1= "<<res1); 
				} while (res1 == LARGER && circ!=circ_done);
				//if circ == circ_done, then prev is the largest
				//else if circ is not larger than prev, than prev is 
				//the largest. anyway, prev is the largest - return it
			}

			//out_edge = prev;
			//found_face = true;
			PRINT_DEBUG ( "new_find_face return " << (*prev).curve() );
			return make_object ((*prev).face());
		}
		else if (res1 ==SMALLER) 
		{
			//if the segment is smaller (cw) than the curve, we need to find 
			//ccw (--) the curve that seg is larger than. 
			//since we can't go --, we will go ++ few times (if necessary). 
			PRINT_DEBUG("res1 == SMALLER : "); 

			//loop 1 - until reach a curve on the other side, if exists
			do {
				prev = circ; 
				circ++;
				cv_orient = compare_xy(v,(*circ).source()->point());
			} while (circ!=circ_done  && seg_orient == cv_orient); 

			if (seg_orient != cv_orient) 
			{
				//loop 2 - until reach the same side again 
				do {
					prev = circ; 
					circ++;
					cv_orient = compare_xy(v,(*circ).source()->point());
				}	while (seg_orient != cv_orient) ;

				//prev is the last edge from the other side, 
				//and curve is the first curve in this side
				res1 =  compare_cw_around_point(seg, (*circ).curve(), v);

				//loop 3 - until find curve > seg cw.
				while (res1 == LARGER && seg_orient==cv_orient) 
				{
					prev = circ; 
					circ++;
					cv_orient = compare_xy(v,(*circ).source()->point());
					if (seg_orient == cv_orient) 
					{
            res1 =  compare_cw_around_point(seg, (*circ).curve(), v);
						PRINT_DEBUG("circ = "<<(*circ).curve()<<" res1= "<<res1);
					}
				}
                
				//now we can say that the output edge is prev
				//out_edge = prev;
				//found_face = true;
				PRINT_DEBUG ( "new_find_face return " << (*prev).curve() );
				return make_object ((*prev).face());
			}
	
			// else - (circ == circ_done)
			//there are no curves on the other side, 
			//find the smallest (cw) on this side
			do {
					prev = circ;
					circ++;
					res1 =  compare_cw_around_point((*circ).curve(), (*prev).curve(),v);//IXX, false);
					PRINT_DEBUG("circ = "<<(*circ).curve()<<" res1= "<< res1); 
			} while (res1 == LARGER && circ!=circ_done);
			//now circ < prev ==> circ is the smallest

			//if seg > smallest, smallest++,
			// otherwise, out_edge  = smallest->twin();
			res1 =  compare_cw_around_point(seg, (*circ).curve(), v);
			if (res1 == SMALLER) 
			{
				//out_edge = circ->twin();
				//found_face = true;
				PRINT_DEBUG ( "new_find_face return " << (*circ).curve() );
				return make_object ((*circ).twin()->face());
			}

			//else: seg > smallest
			circ_done = circ; 
			do  
			{
				prev = circ;
				circ++;
				res1 =  compare_cw_around_point(seg, (*circ).curve(), v);
			} while (res1 == LARGER && circ!= circ_done);

			//out_edge = prev;
			//found_face = true;
			PRINT_DEBUG ( "new_find_face return " << (*prev).curve() );
			return make_object ((*prev).face());
		}
		else //EQUAL
		{
			//TODO: specail case - new vertex or on edge ot something
			PRINT_DEBUG ( "specail case: seg is equal cw to circ " 
				<< (*circ).curve() );
			if (traits->equal_2_object()(p,(*circ).source()->point())) 
			{
				PRINT_DEBUG ( "p is on a vertex ");
				//out_edge = circ;
				//out_vertex = circ->source();
				//lt = Planar_map::VERTEX;
				//found_vertex_or_edge = true;
				return make_object ((*circ).source());
			}

			if (traits->is_in_x_range_2_object()((*circ).curve(),p) && 
				traits->compare_y_at_x_2_object()(p,(*circ).curve()) == EQUAL) 
			{
				// p lies on cv1  
				PRINT_DEBUG ( "p is on an edge ");
        Halfedge_const_handle temp_he = circ;
				//out_edge = circ;
				//lt = Planar_map::EDGE;
				//found_vertex_or_edge = true;
				return make_object (temp_he); 
			}
	
			//p does not lie on cv1 ==> 
			// the target of the equal curve is a better vertex to p 
			PRINT_ERROR("WARNING 11: found closer vertex during new_find_face");
			// out_vertex is the closer vertex
			//out_vertex = circ->source();
			new_vertex = true;
			PRINT_DEBUG( "The new vertex is: "<< (*circ).source()->point() );
			// check validity (the new vertex is vetween them on a line) @@@@
			return make_object((*circ).source());		
		}
	}

	PRINT_ERROR("ERROR 13: new_find_face did not find the face !");
	LM_DEBUG(getchar());
	return Object();
}		

//----------------------------------------------------
// walks from the edge to the point
// param eh - (input) halfedge that the closest point to p is located on
// param p - (input) the point to locate.
// param np - (input) the point on the edge to start the walk .
//
template <class Arrangement, class Arr_landmarks_generator>
Object Arr_landmarks_point_location<Arrangement, Arr_landmarks_generator>
::_walk_from_edge(Halfedge_const_handle eh,       
				  const Point_2 & p,       
				  const Point_2 & np)   const    

{ 
	PRINT_DEBUG("inside walk_from_edge. p= "<< p << ", eh = "
		<<eh->source()->point() << "-->"  <<eh->target()->point());
	X_monotone_curve_2 cv = eh->curve() ;
	Point_2 src = eh->source()->point();
	Point_2 trg = eh->target()->point();
	Comparison_result res;

	LM_DEBUG( 
		if (! traits->is_in_x_range_2_object()(cv, np)) 
			std::cout<<"WARNING 5: np is not on the edge's x_range"
			<<std::endl; 
		else if (traits->compare_y_at_x_2_object()(np,cv) != EQUAL)
			std::cout<<"WARNING 6: np is not on the edge"
			<<std::endl; 
	);
	PRINT_DEBUG("inside walk_from_edge. p= "<< p	<< 
								 ", eh = "<<eh->source()->point() << "-->"  
								 <<eh->target()->point());

	// I deleted the special case: if cv is vertical: 
	// If p ON cv - o.k
	// if p not in x range - o.k.
	// if p in the same x range and is not on (meaning below or above) 
	//   then it might have been better to take the vertex and not the face. 

	//check if p equals one of the edge's end points
	if (traits->equal_2_object()(p, src))
	{
		Vertex_const_handle vh = eh->source();
		return make_object(vh);
	}
	if (traits->equal_2_object()(p, trg))
	{
		Vertex_const_handle vh = eh->target();
		return make_object(vh);
	}

	//if p is in eh's x_range, then we need to check if it is above/below eh
	//and orient the halfedge eh accordingly, so that it will point to the face 
	//that is most likely containing p
	if (traits->is_in_x_range_2_object()(cv, p))
	{
		//check if p is above/below cv
		res = 	traits->compare_y_at_x_2_object()(p,cv);
		PRINT_DEBUG("curve compare y at x: p= "<< p << ", cv =  "<< cv 
					<<", res = "<<res);
		switch (res) { 
			case EQUAL://p is on cv - found !
				return make_object(eh);
			case LARGER:  //p is above cv
				//orient e from left to right
				if (traits->compare_x_2_object()(src,trg) == LARGER) 
					eh = eh->twin();//it is oriented from right to left
				break;
			case SMALLER: //p is below cv
				//orient e from right to left
				if (traits->compare_x_2_object()(src,trg) == SMALLER ) 
					eh = eh->twin();//it is oriented from right to left
				break;
		}	
		PRINT_DEBUG("call from walk_from_edge to walk_from face: eh= "
									<<eh->source()->point() << "-->"  
									<<eh->target()->point());
		return _walk_from_face (eh->face(), p, np);
	}

	//if p is in NOT in eh's x_range,
	//we will check if p is on the left or right to eh, 
	// and take this vertex to start with.
	else 
	{
		Vertex_const_handle vh = eh->source();
		if (traits->compare_xy_2_object()(src, trg) != 
			traits->compare_xy_2_object()(p, src)) 
		{
			vh = eh->target();
		}
		PRINT_DEBUG("call from walk_from_edge to walk_from_vertex: vh= "
					<<vh->point());
		return _walk_from_vertex(vh, p);
	}

}

//----------------------------------------------------
// walks from the face to the point
// \param eh - (input) halfedge that points to the face that the 
//                     closest point to p (np) is located in
// \param p - (input) the point to locate.
// \param np - (input) the point to start the walk .
//
// new algorithm: 
// 1. check if p is in face. 
// 	yes: 
// 		go over holes. for each hole h:
// 			is p in h?
// 				yes: face = h. goto 1.
// 	no: 
// 		call new function: 
// 			find_closest_intersection_in_face( gets face, v, p) 
//                 that  returns e and intersection point (if needed). 
// 			(the function go over all edges surrounding p. 
//           for each one- as done now - take out of the walk procedure) 
//             face = e->twin. goto1.
// 			if intersection not found --- ? error ? 
//
template <class Arrangement, class Arr_landmarks_generator>
Object Arr_landmarks_point_location<Arrangement, Arr_landmarks_generator>
::_walk_from_face(Face_const_handle face,       
				  const Point_2 & p,
				  const Point_2 & np)   const   
{ 
	PRINT_DEBUG("inside walk_from_face. p= "<< p ); 

	//inits
	//Halfedge_const_handle out_edge = eh;
	//Face_const_handle face = out_edge->face() ; //get the face that is the best potentially contains p.
	flipped_edges.clear();  //remove all elements from the flipped_edges list
	bool p_in_face = false;
	Ccb_halfedge_const_circulator h_circ;
	bool found_edge = false;
	bool found_vertex = false;
	Halfedge_const_handle out_edge;
	
	LM_CLOCK_DEBUG( bool first_hit = true; );
	LM_CLOCK_DEBUG(clock_t na_time_start = clock());//time start

	do {
		PRINT_DEBUG(std::endl << "inside loop on face ");
		p_in_face = false;
		if (face->is_unbounded())  {
			p_in_face = true;
			PRINT_DEBUG("unbounded face ");
		}
		else {		
			h_circ = face->outer_ccb();
			LM_CLOCK_DEBUG(clock_t is_point_ts = clock());
			p_in_face = _is_point_in_face(p, h_circ, found_edge, 
				found_vertex, out_edge); 
			LM_CLOCK_DEBUG( if (! p_in_face) first_hit = false );
			LM_CLOCK_DEBUG(clock_t is_point_te = clock() );
			LM_CLOCK_DEBUG(clock_is_point += 
				(double)(is_point_te - is_point_ts));
			PRINT_DEBUG("is_point_in_face returned  "<<  p_in_face );
		}
		if (found_vertex)
		{
			LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
			Vertex_const_handle v = out_edge->target(); 
			return make_object(v); //is it really the target?
		}
		else if (found_edge) 
		{
			LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
			Halfedge_const_handle h = out_edge; 
			return make_object(h);
		}
		else if (p_in_face){
			//check holes
			PRINT_DEBUG(" p in face. go over holes" );
			Holes_const_iterator hole_it  = face->holes_begin();
			Holes_const_iterator hole_end = face->holes_end();	
			bool p_in_hole;
			while (hole_it != hole_end) 
			{
				PRINT_DEBUG(" loop on holes");
				LM_CLOCK_DEBUG(clock_t is_point_ts = clock());
				p_in_hole = _is_point_in_face(p, *hole_it, 
					found_edge, found_vertex, out_edge); 
				LM_CLOCK_DEBUG(clock_t is_point_te = clock() );
				LM_CLOCK_DEBUG(clock_is_point += 
					(double)(is_point_te - is_point_ts));
				if (found_vertex)
				{
					LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
					return make_object(out_edge->target()); 
					//is it really the target?
				}
				else if (found_edge) 
				{
					LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
					return make_object(out_edge);
				}
				else if (p_in_hole) {
					LM_CLOCK_DEBUG( first_hit = false );
					h_circ = *hole_it; 
					//update the new "face " to be the hole. check its holes.			

					LM_CLOCK_DEBUG(clock_t find_edge_time_start = clock());
					if ( _find_edge_to_flip (p, np, h_circ , out_edge) ) {
						LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
						LM_CLOCK_DEBUG(clock_find_edge += 
							(double) (find_edge_time_end - find_edge_time_start));
						out_edge = out_edge->twin();
						face = out_edge->face();
					}
					else {
						LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
						LM_CLOCK_DEBUG(clock_find_edge += 
							(double) (find_edge_time_end - find_edge_time_start));
						PRINT_ERROR( "ERROR 10:  intersection not found");
						LM_DEBUG(getchar());
						return make_object (p_arr->unbounded_face());
					}

					hole_it = hole_end; //to get out of the loop
					p_in_face = false;
				} 
				else {
					++hole_it;
				}
			}
		}
		else {
			//find edge to switch face to (face is never unbounded)
			LM_CLOCK_DEBUG(clock_t find_edge_time_start = clock());
			if ( _find_edge_to_flip (p, np, face->outer_ccb() , out_edge) ) {
				LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
				LM_CLOCK_DEBUG(clock_find_edge += 
					(double) (find_edge_time_end - find_edge_time_start));
				out_edge = out_edge->twin();				
				face = out_edge->face();
				PRINT_DEBUG("after find_edge_to_flip. changed to twin @ " );
				PRINT_DEBUG("out_edge  =  "<< out_edge->source()->point() 
					        <<" -->"<< out_edge->target()->point() );
			}
			else {
				LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
				LM_CLOCK_DEBUG(clock_find_edge += 
					(double) (find_edge_time_end - find_edge_time_start));
				PRINT_ERROR("ERROR 9:  intersection not found");
				LM_DEBUG(getchar());
				//e = p_arr->halfedges_end();
				//lt = Planar_map::UNBOUNDED_FACE;
				LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
				return make_object (p_arr->unbounded_face());
			}
		}
	}while (!p_in_face); 

	LM_CLOCK_DEBUG(clock_t na_time_end = clock() );
	LM_CLOCK_DEBUG(clock_new_alg += (double) (na_time_end - na_time_start) );

	if (face == p_arr->unbounded_face()) 
	{
		PRINT_DEBUG("before return from walk from face. unbounded face.");
		LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
		return make_object (p_arr->unbounded_face());
	}

	PRINT_DEBUG("before return from walk from face. ");
	LM_CLOCK_DEBUG( if (first_hit) number_of_hits++ );
	return make_object (face);
}

//----------------------------------------------------
/*!
* Checks if p is inside the face
* 
* \param p - the input point.
* \param e - the input edge
* \param found_edge - did we found the edge p lies on
* \param new_vertex - output bool, if a closer vertex to p was found
*/

template <class Arrangement, class Arr_landmarks_generator>
bool Arr_landmarks_point_location<Arrangement, Arr_landmarks_generator>
::_is_point_in_face (const Point_2 & p,                  //input seg src 
					 const Ccb_halfedge_const_circulator & face, //input face
					 bool & found_edge,                  //out is p on e?
					 bool & found_vertex,      //out is p equals e->target?
					 Halfedge_const_handle  & out_edge) const  //output if on curve				
{
	LM_CLOCK_DEBUG( entries_is_point_in_face ++) ;
	PRINT_DEBUG("inside is_point_in_face. face = " << (*face).source()->point()
								<<"-->" << (*face).target()->point());

	found_edge = false;
	found_vertex = false;
	int number_of_edges_above_p = 0;

	//loop on all edges in this face
	Ccb_halfedge_const_circulator curr = face;
	Ccb_halfedge_const_circulator last =  curr;

	typename Traits_wrapper_2::Equal_2 equal = traits->equal_2_object();
	typename Traits_wrapper_2::Compare_xy_2     compare_xy = 
                                       traits->compare_xy_2_object();
	typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x = 
                                       traits->compare_y_at_x_2_object();


	do  {
		X_monotone_curve_2 cv = (*curr).curve();
		Point_2 p1 = (*curr).source()->point();
		Point_2 p2 = (*curr).target()->point();

		//check if p equals one of the endpoints of e
		if (equal(p, p1))   {
			found_vertex = true;
			out_edge = (*curr).twin() ;
			PRINT_DEBUG("p is "<< p );
			PRINT_DEBUG("out_edge is "<< out_edge->curve() );
			PRINT_DEBUG("target is "<< out_edge->target()->point() );
			return (true); 
		}
		if (equal(p, p2))   {
			found_vertex = true;
			out_edge = curr;
			PRINT_DEBUG("p is "<< p );
			PRINT_DEBUG("out_edge is "<< out_edge->curve() );
			PRINT_DEBUG("target is "<< out_edge->target()->point() );
			return (true); 
		}

		//check in_x_range lexicographically. 
		//This is if p is on different sides from p1 and p2
		if  (compare_xy(p, p1) != compare_xy(p, p2) ) 
		{
			//check cv to see it p is above, below or on cv
			Comparison_result compare_y_at_x_res = compare_y_at_x(p, cv);

			switch (compare_y_at_x_res) 
			{
			case EQUAL:
				found_edge = true;
				out_edge = curr;
				return (true); 
			case SMALLER: //p is below cv
				//count cv
				number_of_edges_above_p ++;
				break;
			case LARGER: //p is above cv
				//don't count cv. continue
				break;
			default: //should not be equal
				CGAL_assertion (false);
			}
		}

		++curr;
	} while (curr != last);  

	//if number_of_edges_above_p is odd, then p is inside the face, 
	//else - p is outside f. 
	//to check if number_of_edges_above_p is odd/even, 
	//simply do bitwise AND operator with 1. 
	return (number_of_edges_above_p & 1) ;
}

//----------------------------------------------------
/*!
* find the edge in the face that is intersecting with the segment v - p.  
* (out_edge)
* each edge can be chosen only once, and will get into the switch_curves list 
* returns false if there is no intersection, thus no edge to flip. 
* 
* \param p - the input point.
* \param v - the input closest point.
* \param face - the input face
* \param out_edge - the output edge
*/
template <class Arrangement, class Arr_landmarks_generator>
bool Arr_landmarks_point_location<Arrangement, Arr_landmarks_generator>
::_find_edge_to_flip (const Point_2 & p,                  //input seg src 
					  const Point_2 & np,           //input closest point
					  const Ccb_halfedge_const_circulator & face, //input face
					  Halfedge_const_handle  & out_edge) const //output edge to flip
{
	//we want to eliminate the calls to nearest. 
	//this means we will not call to find_real_intersection at all. 
	//this also means that we will take the first edge that is intersecting, 
	//and not check which is closer. 
	//what we will check is whether this edge was selected already 
	//(and flipped),  in this case we will not flip it again, 
	//but move to the next edge
	LM_CLOCK_DEBUG( entries_find_edge++ );
	PRINT_DEBUG("inside find_edge_to_flip.   " );

	Point_2 vp = np;
	X_monotone_curve_2 seg(vp, p);   //create a segment vp--p. 

	//loop on all edges in this face
	Ccb_halfedge_const_circulator curr = face;
	Ccb_halfedge_const_circulator last =  curr;
	bool p_in_x_range, v_in_x_range, p1_in_x_range, p2_in_x_range;	

	typename Traits_wrapper_2::Is_in_x_range_2         is_in_x_range = 
                                      traits->is_in_x_range_2_object();


	do  {
		X_monotone_curve_2 cv = (*curr).curve();
		Point_2 p1 = (*curr).source()->point();
		Point_2 p2 = (*curr).target()->point();

		// check if the curve was already flipped - in this case, 
		//    don't check it at all.
		Std_edge_iterator found1 = std::find (flipped_edges.begin(), 
			flipped_edges.end(), curr);
		Std_edge_iterator found2 = std::find (flipped_edges.begin(), 
			flipped_edges.end(), (*curr).twin());
		if (found1 != flipped_edges.end() || found2 != flipped_edges.end()) {
			PRINT_DEBUG("curve "<<p1<<"-->"<<p2<<" was found in the list");
		}
		else {
			PRINT_DEBUG("curr = " << p1 << "-->" << p2 );

			//check x-range
			p_in_x_range = is_in_x_range(cv, p);
			v_in_x_range = is_in_x_range(cv, vp);
			p1_in_x_range = is_in_x_range(seg, p1);
			p2_in_x_range = is_in_x_range(seg, p2);
			PRINT_DEBUG("p_in_x_range = " << p_in_x_range 
				<< " , v_in_x_range = " << v_in_x_range);
			PRINT_DEBUG("p1_in_x_range = " << p1_in_x_range 
				<< " , p2_in_x_range = " << p2_in_x_range);

			if (p_in_x_range || v_in_x_range || p1_in_x_range || p2_in_x_range)
			{		
				bool intersect;
				
				LM_CLOCK_DEBUG(clock_t check_app_ts = clock());	
				bool check_res = _check_approximate_intersection
					(seg, cv, intersect);
				LM_CLOCK_DEBUG(clock_t check_app_te = clock());
				LM_CLOCK_DEBUG(clock_check_app += 
					(double) (check_app_te - check_app_ts));

				if (check_res) 
				{
					if (intersect) {
						out_edge = curr;
						flipped_edges.push_back(out_edge);
						return (true);
					}
				}
				else {
					PRINT_ERROR("ERROR 12: check_approximate_intersection \
								did not return an answer.");
				}
			}
			//else - there is not intersection.
		}

		++curr;
	} while (curr != last);  

	return (false);
}


//----------------------------------------------------
//
// Trim the segment to the x range of the input curve
// 
// \param seg - the input seg.
// \param cv - the input curve
// \param trimmed_seg - the output trimmed segment
//
template <class Arrangement_2, class Arr_landmarks_generator>
bool Arr_landmarks_point_location<Arrangement_2,Arr_landmarks_generator>
::_check_approximate_intersection (const X_monotone_curve_2 & seg,  
								  const X_monotone_curve_2 & cv,
								  bool & intersect) const 
{
	LM_CLOCK_DEBUG( entries_to_check_app ++ );

	typename Traits_wrapper_2::Equal_2              equal = 
                                            traits->equal_2_object();
	typename Traits_wrapper_2::Compare_xy_2          compare_xy = 
                                            traits->compare_xy_2_object();
	typename Traits_wrapper_2::Compare_y_at_x_2     compare_y_at_x = 
                                            traits->compare_y_at_x_2_object();
	typename Traits_wrapper_2::Is_in_x_range_2      is_in_x_range = 
                                            traits->is_in_x_range_2_object();
	typename Traits_wrapper_2::Compare_y_at_x_right_2 compare_y_at_x_right = 
                                      traits->compare_y_at_x_right_2_object();
	typename Traits_wrapper_2::Compare_y_at_x_left_2 compare_y_at_x_left = 
                                      traits->compare_y_at_x_left_2_object();

	Point_2 seg_right = traits->construct_max_vertex_2_object()(seg);
	Point_2 seg_left = traits->construct_min_vertex_2_object()(seg);
	Point_2 cv_right = traits->construct_max_vertex_2_object()(cv);
	Point_2 cv_left = traits->construct_min_vertex_2_object()(cv);
	intersect = false;

	PRINT_DEBUG("seg_right =  " << seg_right << " , seg_left = " << seg_left);
	PRINT_DEBUG("cv_right = " << cv_right << " , cv_left = " << cv_left);
	
	//compare the 2 left end-points and the 2 right end-points
	Comparison_result comp_left_xy_res = compare_xy(seg_left, cv_left);
	Comparison_result comp_right_xy_res = compare_xy(seg_right, cv_right);

	// the left end-point of the segments is equal to the left end-point 
	//of the curve
	if (comp_left_xy_res == EQUAL) 
	{
		PRINT_DEBUG("the left end-point of the segments is equal to the \
					left end-point of the curve");
		// compare to the right of the curves
		Comparison_result curves_comp = compare_y_at_x_right(seg,cv,seg_left);
	
		if (curves_comp == EQUAL) {
			PRINT_DEBUG("overlap !!!");
			return (false); 
		}                                                                         
		if (comp_right_xy_res == SMALLER) 
		{	//the segments ends before (to the right) of the curve 
			Comparison_result curve_comp_right_res = 
				compare_y_at_x(seg_right,cv);

            if (curve_comp_right_res == EQUAL) {
				PRINT_DEBUG("2 points collide");
				return (false); 
			}
			if (curves_comp != curve_comp_right_res) 
			{ //the segment is on the other side of the segment's end-point
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
			else {
				//intersect = false;
				return (true);
			}
		}
		else if (comp_right_xy_res == LARGER) 
		{ //the curve ends before (to the right) of the segment
			Comparison_result curve_comp_right_res =
				compare_y_at_x(cv_right,seg) ;			

			if (curve_comp_right_res == EQUAL) {
				PRINT_DEBUG("2 points collide");
				return (false); 
			}
			if (curves_comp == curve_comp_right_res) 
			{//the segment is on the same side as the curve's end-point
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
			else {
				//intersect = false;
				return (true);
			}
		}
		else { //(comp_right_xy_res == EQUAL) 
				PRINT_DEBUG("2 endpoints collide");
				return (false); 
		}
	}

	// the right end-point of the segments is equal to the right end-point 
	//of the curve
	else if (comp_right_xy_res == EQUAL) 
	{
		PRINT_DEBUG("the right end-point of the segments is equal to the \
					right end-point of the curve");
		// compare to the left of the curves
		Comparison_result curves_comp = compare_y_at_x_left(seg, cv, seg_right);

		if (curves_comp == EQUAL) {
			PRINT_DEBUG("overlap !!!");
			return (false); 
		}
		if (comp_left_xy_res == SMALLER) 
		{ //the curve ends before (to the left) of the segment 
			Comparison_result curve_comp_left_res = 
				compare_y_at_x(cv_left,seg) ;

			if (curve_comp_left_res == EQUAL) {
				PRINT_DEBUG("2 points collide");
				return (false); 
			}
			if (curves_comp == curve_comp_left_res) {
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
			else {
				//intersect = false;
				return (true);
			}
		}
		else if (comp_left_xy_res == LARGER) 
		{ //the segment ends before (to the left) of the curve 
			Comparison_result curve_comp_left_res = 
				compare_y_at_x(seg_left,cv) ;

			if (curve_comp_left_res == EQUAL) {
				PRINT_DEBUG("2 points collide");
				return (false); 
			}
			if (curves_comp != curve_comp_left_res) {
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
			else {
				//intersect = false;
				return (true);
			}
		}
		else { //(comp_left_xy_res == EQUAL) 
				PRINT_DEBUG("2 endpoints collide");
				return (false); 
		}
	}

	//one curve is inside the other curve's x-range (fully) 
	else if (comp_left_xy_res != comp_right_xy_res) { 
		if (comp_left_xy_res == LARGER) 
		{	//the segment is inside the curve's x-range
			PRINT_DEBUG("the segment is inside the curve's x-range.");
			LM_DEBUG ( //checks for debugging - remove later
				if (! is_in_x_range(cv,seg_left) ) {
					PRINT_ERROR( "! is_in_x_range(cv,seg_left) "); 
					return (false);}
				if (! is_in_x_range(cv,seg_right) )  {
					PRINT_ERROR("(! is_in_x_range(cv,seg_right) "); 
					return (false);}
			)

			Comparison_result curve_comp_left_res = 
				compare_y_at_x(seg_left,cv);
			Comparison_result curve_comp_right_res =
				compare_y_at_x(seg_right,cv) ;

			if ((curve_comp_left_res == EQUAL)||
				(curve_comp_right_res == EQUAL))
			{
				//PRINT_ERROR(" WARNING 7: left or right endpoint of the 
				//segment is on the curve ");
				//this should not happen since p is on the curve, 
				//we should have find it already, and if v is on the curve, 
				//than v should have cut the curve in two
				//this can happen if we're walking from edge. 
				PRINT_DEBUG("no intersection");
				//intersect = false;
				return (true);
				//return (false);
			}
			if  (curve_comp_left_res == curve_comp_right_res) 
			{	//no intersection
				PRINT_DEBUG("no intersection");
				//intersect = false;
				return (true);
			}
			else 
			{
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
		}
		else
		{	//the curve is inside the segments's x-range
			PRINT_DEBUG("the curve is inside the segments's x-range");
			LM_DEBUG ( //checks for debugging - remove later
				if (! is_in_x_range(seg,cv_left) ) {
				 PRINT_ERROR( "! is_in_x_range(seg,cv_left)"); 
				 return (false);}
				if (! is_in_x_range(seg,cv_right) ) {
				 PRINT_ERROR("! is_in_x_range(seg,cv_right)"); 
				 return (false);}
			)

			Comparison_result curve_comp_left_res = compare_y_at_x(cv_left,seg);
			Comparison_result curve_comp_right_res = compare_y_at_x(cv_right,seg);

			if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL) 
			{
				PRINT_ERROR(" WARNING 8: left or right endpoint of the curve \
							is on the segment ");
				//this means that the endpoint of the curve is a closer 
				//vertex to p than v is. 
				return (false);
			}
			if  (curve_comp_left_res == curve_comp_right_res) 
			{	//no intersection
				PRINT_DEBUG("no intersection");
				//intersect = false;
				return (true);
			}
			else 
			{
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);
			}
		}
	}
	else 
	{//one endpoint of the first curve is inside the other curve's 
		//x-range and one endpoitn is out, and vise versa.
		if (comp_left_xy_res == SMALLER) 
		{ //if the segment right point is in and the curve's left point
			PRINT_DEBUG("the segment right point is in and the \
						curve's left point");
			LM_DEBUG ( //checks for debugging - remove later
				if (! is_in_x_range(seg,cv_left) ) {
				 PRINT_ERROR("! is_in_x_range(seg,cv_left)"); 
				 return (false);}
				if (! is_in_x_range(cv,seg_right) ) {
				 PRINT_ERROR( "! is_in_x_range(cv,seg_right)"); 
				 return (false);}
			)

			Comparison_result curve_comp_left_res = 
				compare_y_at_x(cv_left,seg);
			Comparison_result curve_comp_right_res = 
				compare_y_at_x(seg_right,cv);

			if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL) 
			{
				if (equal(seg_right, cv_left)) {
					//this case is o.k., just return that there is 
					//no intersection, because probably v is the same 
					//for both curves
					PRINT_DEBUG("segment's right endpoint = \
								curve's left endpoint ");
					PRINT_DEBUG("no intersection");
					//intersect = false;
					return (true);		
				}
				PRINT_ERROR(" WARNING 9: left or right endpoint \
							is on the curve ");
				//its either the segment's end-point is on the curve, 
				//and this should not happen since if p is on the curve, 
				//we should have find it already, and if v is on the curve, 
				//than v should have cut the curve in two.
				//another options is that the curve end-point is on 
				//the segment, and in this case this endpoint is a closer 
				//vertex to p than v      
				return (false);
			}
			if  (curve_comp_left_res == curve_comp_right_res) 
			{	//intersection
				PRINT_DEBUG("intersecting");
				intersect = true;
				return (true);				
			}
			else 
			{
				PRINT_DEBUG("no intersection");
				//intersect = false;
				return (true);			
			}
		}
		else 
		{ //if the curve's right point is in and the segment's left point
			PRINT_DEBUG("the curve's right point is in and the \
						segment's left point");
			LM_DEBUG ( //checks for debugging - remove later
				if (! is_in_x_range(cv,seg_left) ) {
					PRINT_ERROR("! is_in_x_range(cv,seg_left) ***  "); 
					return (false);	}
				if (! is_in_x_range(seg,cv_right) ) {
					PRINT_ERROR("! is_in_x_range(seg,cv_right) *** "); 
					return (false); }
			)

			Comparison_result curve_comp_left_res = 
				compare_y_at_x(seg_left,cv);
			Comparison_result curve_comp_right_res = 
				compare_y_at_x(cv_right,seg);

			PRINT_DEBUG("curve_comp_right_res="<<curve_comp_right_res);
			PRINT_DEBUG("curve_comp_left_res= "<<curve_comp_left_res);

			if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL)
			{
				if (equal(seg_left, cv_right)) {
					//this case is o.k., just return that 
					//there is no intersection, because probably 
					//v is the same for both curves
					PRINT_DEBUG("segment's left endpoint = \
								curve's right endpoint ");
					PRINT_DEBUG("no intersection");
					//intersect = false;
					return (true);		
				}
				PRINT_ERROR(" WARNING 10: left or right endpoint \
							is on the curve ");
				//same explanation as in warning 9
				return (false);
			}
			if  (curve_comp_left_res == curve_comp_right_res) 
			{	//intersection
				PRINT_DEBUG("intersecting . ");
				intersect = true;
				return (true);				
			}
			else 
			{
				PRINT_DEBUG("no intersection .");
				//intersect = false;
				return (true);			
			}
		}
	}


	PRINT_ERROR("ERROR 11: not existing option in \
				check_approximate_intersection ");
	LM_DEBUG(getchar());
	return (false);
}

CGAL_END_NAMESPACE

#endif
