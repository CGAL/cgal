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

#ifndef CGAL_PM_LANDMARKS_POINT_LOCATION_C
#define CGAL_PM_LANDMARKS_POINT_LOCATION_C

#include <CGAL/Pm_landmarks_point_location.h>

//#define CGAL_LM_DEBUG
#define LM_CLOCK_DEBUG

CGAL_BEGIN_NAMESPACE

//----------------------------------------------------
/*!
*/
//if unbounded face - returns NULL or some edge on unbounded face 
//if its a vertex returns a halfedge pointing _at_ it
template <class Planar_map, class Nearest_neighbor>
typename Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::Halfedge_const_handle
Pm_landmarks_point_location<Planar_map, Nearest_neighbor>
::locate(const Point_2 & p, Locate_type & lt) const
{

	#ifdef CGAL_LM_DEBUG
		std::cout << "------ locate point "<< p << std::endl;
	#endif //CGAL_LM_DEBUG

	//init output
	Face_handle f = pm->unbounded_face(), last = pm->faces_end();  
	// invariant: p always lies in f's interior or holes 
	Halfedge_const_handle e = pm->halfedges_end(); // closest halfedge so far
	lt = Planar_map::UNBOUNDED_FACE;

	//check that the nearest neigbor search data structure is updated
	if (pm->number_of_vertices() == 0) 
		return e;
	if (! updated_nn)
	{
		std::cout << "ERROR 1: tree is not updated ! number of vertices in pm = " << pm->number_of_vertices() << std::endl;
		return e;
	}

	#ifdef LM_CLOCK_DEBUG
		clock_t nn_time_start = clock(); 	//check time - start
	#endif

	//get nearest vertex to point p 
	double px = CGAL::to_double(p.x());
	double py = CGAL::to_double(p.y());
	NN_Point_2  nnp (px, py); 
	NN_Point_2  nearest_point = nn.find_nearest_neighbor(nnp); 
	#ifdef CGAL_LM_DEBUG 
		std::cout << "nearest neighbor of point "<< p << " is " << nearest_point.x() <<','<<nearest_point.y()  << std::endl;
	#endif //CGAL_LM_DEBUG
	Vertex_handle nearest_vertex = nearest_point.vertex();

	#ifdef LM_CLOCK_DEBUG
		clock_t nn_time_end = clock(); //check time - end	
		double nn_period = (double) (nn_time_end - nn_time_start); 
		clock_for_nn_search += nn_period;
	#endif

	//walk from the nearest_vertex to the point p, using walk algorithm, 
	//and find the location of p.
	//lt is the  location type as needed, and e is the halfedge pointing at the vertex or face 
	//that the point lies at, or if the point lies on an halfedge, its the halfedge itself.
	#ifdef CGAL_LM_DEBUG 
		std::cout << "call to walk_from_vertex_to_point" << std::endl;
	#endif //CGAL_LM_DEBUG

	#ifdef LM_CLOCK_DEBUG
		clock_t w_time_start = clock();//time start
	#endif

	// !!! the operation !!!
	walk( nearest_vertex , p, e, lt);

	#ifdef LM_CLOCK_DEBUG
		clock_t w_time_end = clock(); //time end
		double w_period = (double) (w_time_end - w_time_start);
		clock_for_walk += w_period;
	#endif

	#ifdef CGAL_LM_DEBUG
		std::cout << "return from walk_from_vertex_to_point. !!!" << std::endl;
		std::cout << "lt = "<<lt <<", e = "<<e->source()->point() << "->" 
			<< e->target()->point() << std::endl;
		std::cout  << std::endl  << std::endl  << std::endl;
		//getchar
	#endif //CGAL_LM_DEBUG

	return e;
}

//----------------------------------------------------
/*!
*/
template <class Planar_map, class Nearest_neighbor>
typename Pm_landmarks_point_location<Planar_map,Nearest_neighbor >::Halfedge_handle
Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
locate(const Point_2 & p, Locate_type & lt) 
{
	Halfedge_handle h = 
		Halfedge_handle_unconst(((const_Self_ptr)this)->locate(p,lt));
	//Halfedge_handle h=((cPLp)this)->locate(p,lt);

	return h;
}


//----------------------------------------------------
/*!
*/
template <class Planar_map, class Nearest_neighbor>
typename Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::Halfedge_const_handle
Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
vertical_ray_shoot(const Point_2 & p, 
				   Locate_type & lt, 
				   bool up) const
{
	//std::cout << "ERROR: Vertical ray shoot NOT suppored in CDT point location"
	//<< std::endl;
	CGAL_precondition_msg(false, 
		"Vertical ray shoot NOT suppored in CDT point location");
	assert(false);
	Halfedge_const_handle e = pm->halfedges_end(); // closest halfedge so far  
	lt=Planar_map::UNBOUNDED_FACE;
	return e;
}

//----------------------------------------------------
/*!
*/
template <class Planar_map, class Nearest_neighbor>
typename Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::Halfedge_handle
Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::vertical_ray_shoot(
	const Point_2& p, 
	Locate_type& lt, bool up)
{
	//Halfedge_const_handle h=((const_Self_ptr)this)->vertical_ray_shoot(p,lt,up);
	Halfedge_handle h =
		Halfedge_handle_unconst(((const_Self_ptr)this)->vertical_ray_shoot(p,lt,up));

	return h;
}

//----------------------------------------------------
/*! insert an halfedge to the plannar map - 
add its end points to the landmarks tree
*/
template <class Planar_map, class Nearest_neighbor>
void Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::insert_halfedge_to_ln_tree
(Halfedge_handle hh, const typename Planar_map::Traits::X_monotone_curve_2 &cv)
{
	#ifdef CGAL_LM_DEBUG
		std::cout << "in insert_halfedge_to_ln_tree" << std::endl ;
		std::cout << "cv = "<< cv << std::endl;
		std::cout << hh->source()->point()<<" towards "<< hh->target()->point()<< std::endl;
	#endif //CGAL_LM_DEBUG 

	create_landmarks_tree();
}

//----------------------------------------------------
/*! go over all vertices, and insert each vertex to the 
landmarks tree.
*/
template <class Planar_map, class Nearest_neighbor>
void Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::create_landmarks_tree()
{ 
	#ifdef CGAL_LM_DEBUG
		std::cout << "in create_landmarks_tree" << std::endl ;
	#endif //CGAL_LM_DEBUG

	if (pm->vertices_begin() == pm->vertices_end()) {
		updated_nn = true;
		#ifdef CGAL_LM_DEBUG
			std::cout << "empty pm. out create_landmarks_tree" << std::endl ;
		#endif //CGAL_LM_DEBUG
	}

	//Go over planar map, and create a triangulation of it
	Vertex_iterator   vit;
	NN_Point_list      plist; 

	for (vit=pm->vertices_begin(); vit != pm->vertices_end(); vit++)
	{
		//get point from vertex
		Point_2 p = vit->point() ;
		double px = CGAL::to_double(p.x());
		double py = CGAL::to_double(p.y());
		Vertex_handle vh = vit;
		NN_Point_2 np (px, py, vh); 
		//insert point into list
		plist.push_back(np); 

		//print
		#ifdef CGAL_LM_DEBUG
			//std::cout << "point is= " << p << std::endl;
		#endif //CGAL_LM_DEBUG
	} 

	#ifdef CGAL_LM_DEBUG
		std::cout << "before clean" << std::endl ;
	#endif //CGAL_LM_DEBUG
	nn.clean();

	#ifdef CGAL_LM_DEBUG
		std::cout << "before init" << std::endl ;
	#endif //CGAL_LM_DEBUG
	nn.init(plist.begin(), plist.end());

	//the triangulation is now updated
	updated_nn = true;

	#ifdef CGAL_LM_DEBUG
		std::cout << "out create_landmarks_tree" << std::endl ;
	#endif //CGAL_LM_DEBUG
}


//----------------------------------------------------
/*! walks from the vertex to the point
*
* Return the number of times that the segment s-->t intersects the curve
* in segments this number can be 0 or 1
* 
* \param cv - the input curve (intersecting curve?)
* \param s - start point of the input segment
* \param t - end of point the segment
* \param p - the output point.
* \pre s is not equal to t
* \return the number of times the curve intersects the segment.
*         0 - no intersection, 1 - intersect, 2- overlap
*
*/

template <class Planar_map, class Nearest_neighbor>
void
//typename Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::Halfedge_const_handle
Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
walk(Vertex_handle nearest_vertex,           //input
						  const Point_2 & p,          //input 
						  Halfedge_const_handle& e,   //output
						  Locate_type& lt)   const    //output
{ 

	#ifdef CGAL_LM_DEBUG 
		std::cout << "inside walk_from_vertex_to_point. p= "<< p	<< ", nearest_vertex = "<<nearest_vertex->point()  << std::endl;
	#endif //CGAL_LM_DEBUG

	bool new_vertex = false;
	bool found_face = false;
	bool found_vertex_or_edge = false;
	bool found_edge = false;
	bool found_vertex = false;
	Face_handle face;
	Vertex_handle out_vertex;
	Vertex_handle vh = nearest_vertex;
	Halfedge_handle out_edge;

	int debug_count = 0; 

	do {
		debug_count++; 
		#ifdef CGAL_LM_DEBUG 
			std::cout << "inside first do loop. debug_count = " <<debug_count<< std::endl;
		#endif //CGAL_LM_DEBUG
		//find the edge out_edge which is the best possibly 
		//pointing to the face containing p

		#ifdef LM_CLOCK_DEBUG
			clock_t ff_time_start = clock();//time start
		#endif

		find_face (p, vh, found_vertex_or_edge, new_vertex, found_face , out_vertex, out_edge, lt);

		#ifdef LM_CLOCK_DEBUG
			clock_t ff_time_end = clock();//time end
			double ff_period = (double) (ff_time_end - ff_time_start);
			clock_ff += ff_period;
		#endif

		#ifdef CGAL_LM_DEBUG 
			std::cout << "find_face returned:" <<std::endl;
			std::cout <<"found_vertex_or_edge = "<<found_vertex_or_edge << std::endl;
			std::cout	<<"new_vertex = "<< new_vertex << std::endl;
			std::cout	<< "found_face="<<found_face << std::endl;
			std::cout	<<"lt = " <<lt<< std::endl;
		#endif //CGAL_LM_DEBUG

		if (found_vertex_or_edge) {
			if (lt == Planar_map::EDGE) {
				#ifdef CGAL_LM_DEBUG 
					std::cout	<<"out_edge = " << out_edge->source()->point() <<"--> "<< out_edge->target()->point() << std::endl;
				#endif //CGAL_LM_DEBUG
				e = out_edge;
			}
			else { // (lt == Planar_map::VERTEX) 
				#ifdef CGAL_LM_DEBUG 
					std::cout	<<"out_vertex = " << out_vertex->point() << std::endl;
				#endif //CGAL_LM_DEBUG
				e = out_vertex->incident_halfedges();
			}
			return;
		}
		if (new_vertex) {
			std::cout << "NEW vertex 1 " << std::endl;
			//check if the new vertex is really closer 
			if (traits->compare_distance(p, out_vertex->point(), vh->point())  
				== SMALLER) {
					vh = out_vertex;
				}
			else {
				std::cout << "Error 2: new vertex is not closer to p than vh! ";
				std::cout << "out_vertex is:  " << out_vertex->point() << std::endl;
				getchar();
				return; 
			}
		}

	} while (new_vertex);

	//get the face that is the best potentially contains p.
	if (found_face) {
		face = out_edge->face();
	}
	else {
		std::cerr << "face not found" << std::endl;
		return;
	}

/////////////////////////////IXXXXXXXXXXXXXXXXXXX
//new algorithm should be: (after face was found) 
//1. check if p is in face. 
//	yes: 
//		go over holes. for each hole h:
//			is p in h?
//				yes: face = h. goto 1.
//	no: 
//		call new function: 
//			find_closest_intersection_in_face( gets face, v, p) that  returns e and intersection point (if needed). 
//			(the function go over all edges surrounding p. for each one- as done now - take out of the walk procedure) 
//            face = e->twin. goto1.
//			if intersection not found --- ? error ? 
//////////////////////////////IXXXXXXXXXXXXXXXXXXX

	// IXXXXXXXXXX
	bool p_in_face = false;
	Ccb_halfedge_circulator h_circ;
	
	#ifdef LM_CLOCK_DEBUG
		clock_t na_time_start = clock(); //time start
	#endif

	#ifdef CGAL_LM_DEBUG
		std::cout << "start new algo " <<std::endl;
	#endif //CGAL_LM_DEBUG

	do {
		#ifdef CGAL_LM_DEBUG
			std::cout <<std::endl << "inside loop on face " <<std::endl;
		#endif //CGAL_LM_DEBUG
		found_vertex = found_edge = p_in_face = false;
		if (face->is_unbounded())  {
			p_in_face = true;
			#ifdef CGAL_LM_DEBUG
				std::cout << "unbounded face " <<std::endl;
			#endif //CGAL_LM_DEBUG
		}
		else {		
			h_circ = face->outer_ccb();
			p_in_face = is_point_in_face(p, h_circ, found_edge, found_vertex, out_edge); 
			#ifdef CGAL_LM_DEBUG
				std::cout << "is_point_in_face returned  "<<  p_in_face <<std::endl;
			#endif //CGAL_LM_DEBUG
		}
		if (found_vertex || found_edge ) {
			lt = found_vertex ? Planar_map::VERTEX : Planar_map::EDGE;
			e = out_edge;
			return;
		}
		if (p_in_face){
			//check holes
			#ifdef CGAL_LM_DEBUG
				std::cout << " p in face. go over holes" <<std::endl;
			#endif //CGAL_LM_DEBUG
			Holes_iterator hole_it  = face->holes_begin();
			Holes_iterator hole_end = face->holes_end();	
			bool p_in_hole;
			while (hole_it != hole_end) 
			{
				#ifdef CGAL_LM_DEBUG
					std::cout << " loop on holes" <<std::endl;
				#endif //CGAL_LM_DEBUG
				p_in_hole = is_point_in_face(p, *hole_it, found_edge, found_vertex, out_edge); 
				if (found_vertex || found_edge ) {
					lt = found_vertex ? Planar_map::VERTEX : Planar_map::EDGE;
					e = out_edge;
					return;
				} 
				if (p_in_hole) {
					h_circ = *hole_it; //update the new "face " to be the hole. check its holes.
					
					#ifdef LM_CLOCK_DEBUG
						clock_t fc_time_start = clock(); //time start
					#endif

					if ( find_closest_intersection_in_face (p, vh, h_circ , out_edge) ) {
						out_edge = out_edge->twin();
						face = out_edge->face();
					}
					else {
						std::cerr << "ERROR 10:  intersection not found" << std::endl;
					}

					#ifdef LM_CLOCK_DEBUG			
						clock_t fc_time_end = clock();//time end
						double fc_period = (double) (fc_time_end - fc_time_start);
						clock_fciif += fc_period;
					#endif

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
			#ifdef LM_CLOCK_DEBUG
				clock_t fc_time_start = clock(); //time start
			#endif

			if ( find_closest_intersection_in_face (p, vh, face->outer_ccb() , out_edge) ) {
				out_edge = out_edge->twin();				
				face = out_edge->face();
				#ifdef CGAL_LM_DEBUG
					std::cout << "after find_closest_intersection_in_face . out_edge changed to twin @ " <<std::endl;
					std::cout << "out_edge  =  " << out_edge->source()->point() <<" -->" 
								     << out_edge->target()->point() <<std::endl;
				#endif //CGAL_LM_DEBUG
			}
			else {
				std::cerr << "ERROR 9:  intersection not found" << std::endl;
			}

			#ifdef LM_CLOCK_DEBUG			
				clock_t fc_time_end = clock();//time end
				double fc_period = (double) (fc_time_end - fc_time_start);
				clock_fciif += fc_period;
			#endif

		}
	}while (!p_in_face); 
	//IXXXXXXXXXXX

	#ifdef LM_CLOCK_DEBUG			
		clock_t na_time_end = clock();//time end
		double na_period = (double) (na_time_end - na_time_start);
		clock_new_alg += na_period;
	#endif


	if (face != pm->unbounded_face()) 
		lt = Planar_map::FACE;
	else 
		lt = Planar_map::UNBOUNDED_FACE;
	e = out_edge;

	return;
}


//----------------------------------------------------
/*!
* Return the edge around vh that is before cw to the segment vh->p
* 
* \param p - the input point.
* \param vh - the input vertex
* \param found_vertex_or_edge - output bool, 
*                  if the point was found on a vertex or halfedge
* \param new_vertex - output bool, if a closer vertex to p was found
* \param out_vertex - the output vertex (if closer vertex found)
* \param out_edge - the output edge, if found p on edge, 
*           or if normally return the edge around vh closer (before cw) to p
*/
template <class Planar_map, class Nearest_neighbor>
void Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_face (const Point_2 & p, 
		   Vertex_handle vh,
		   bool & found_vertex_or_edge, 
		   bool & new_vertex, 
		   bool & found_face,
		   Vertex_handle & out_vertex, 
		   Halfedge_handle & out_edge,
		   Locate_type& lt  ) const
{ 
	#ifdef CGAL_LM_DEBUG 
		std::cout << "inside find_face. p ="<< p <<" , vh = "<<vh->point() <<std::endl;
	#endif //CGAL_LM_DEBUG

	new_vertex = false;
	found_vertex_or_edge = false;
	found_face = false;	

	// check if the point equals the vertex. 
	if (traits->point_equal(vh->point(), p)) {
		lt = Planar_map::VERTEX;
		out_edge =  vh->incident_halfedges();
		out_vertex = vh;
		found_vertex_or_edge = true;
		return;
	}

	//create a segment vh--p. 
	Curve_2 seg(vh->point(), p);

	//circulate vh and find the two edges that vh-p lies between.
	bool cv_equal_cv1 = false;
	bool cv_equal_cv2 = false;
	bool is_btw = false;

	Halfedge_around_vertex_circulator circ = vh->incident_halfedges(); 
	Halfedge_around_vertex_circulator circ_done (circ);
	Halfedge_handle first = circ;
	Halfedge_handle prev = circ;
	Curve_2 cv1, cv2;
	++circ;
	int debug_count = 0;

	while (circ != circ_done) {  
		debug_count++;
		//std::cout << " first loop in find face. debug_count = " <<debug_count<<std::endl;
		cv2 = circ->curve();
		cv1 = prev->curve();
		//std::cout << "circ = "<<circ->curve() << std::endl;    
		//check if seg is between prev_cv and curr_cv
		//std::cout << " cv1 = " <<cv1 << " cv2 = " <<cv2<< std::endl;
		is_btw = traits->curve_is_between_cw(seg, cv1, cv2, vh->point(), 	cv_equal_cv1, cv_equal_cv2);
		//std::cout << " is_between returned: " << is_btw << std::endl;

		if (cv_equal_cv1 || cv_equal_cv2) {
			//std::cout << " cv_equal_cv1= " <<cv_equal_cv1 << " cv_equal_cv2= " <<cv_equal_cv2 << std::endl;

			if (cv_equal_cv1) {
				if (traits->point_equal(p,prev->source()->point())) {
					out_edge = prev;
					out_vertex = prev->source();
					lt = Planar_map::VERTEX;
					found_vertex_or_edge = true;
					return; 
				}
				if (traits->point_in_x_range(cv1,p) && 
					traits->curve_compare_y_at_x(p,cv1) == EQUAL) {
						// p lies on cv1    
						out_edge = prev;
						lt = Planar_map::EDGE;
						found_vertex_or_edge = true;
						return; 
					}
					//p does not lie on cv1 ==> 
					// the target of the equal curve is a better vertex to p 
					std::cout << " WARNING: found closer vertex during find_face" <<std::endl;
					// out_vertex =  the closer vertex
					out_vertex = prev->source();
					new_vertex = true;
					std::cout << "The new vertex is: "<< out_vertex->point() << std::endl;
					// check validity (the new vertex is vetween them on a line) @@@@
					return;
			}
			else { //cv_equal_cv2
				if (traits->point_equal(p,circ->source()->point())) {
					out_edge = circ;
					out_vertex = circ->source();
					lt = Planar_map::VERTEX;
					found_vertex_or_edge = true;
					return; 
				}
				if (traits->point_in_x_range(cv2,p) && 
					traits->curve_compare_y_at_x(p,cv2) == EQUAL) {
						// p lies on cv1    
						out_edge = circ;
						lt = Planar_map::EDGE;
						found_vertex_or_edge = true;
						return; 
					}
					//p does not lie on cv1 ==> 
					// the target of the equal curve is a better vertex to p 
					std::cout << " WARNING: found closer vertex " <<std::endl;
					// out_vertex =  the closer vertex
					out_vertex = prev->source();
					new_vertex = true;
					std::cout << "The new vertex is: "<< out_vertex->point() << std::endl;
					// check validity (the new vertex is vetween them on a line) @@@@
					return;
			}

		}

		if (is_btw) {
			#ifdef CGAL_LM_DEBUG 
				std::cout << " cv is between " << cv1 << " and " << cv2 << std::endl;
			#endif //CGAL_LM_DEBUG
			out_edge = prev;
			found_face = true;
			//std::cout << " before return from find_face. out_edge =   " << out_edge->curve() << std::endl;
			return;
		}

		prev = circ;
		//std::cout << " prev = " <<prev->curve() << std::endl;
		++circ;
	}

	//std::cout << " not found during first loop !" <<std::endl;

	//if p not found between edges so far, try prev->first
	if ( !is_btw && !cv_equal_cv1 && !cv_equal_cv2) { 
		cv1 = prev->curve();
		cv2 = first->curve();
		//std::cout << " before  curve_is_between_cw "<< std::endl;
		//std::cout << " cv1 = " <<cv1 << " , cv2 = " << cv2 << std::endl;
		//std::cout << " seg = " <<seg << " ,vh point = " << vh->point() << std::endl;
		is_btw = traits->curve_is_between_cw(seg, cv1, cv2 , vh->point(), 
			cv_equal_cv1, cv_equal_cv2);
		//std::cout << " is_btw =  " <<is_btw << std::endl;
		//std::cout << " cv_equal_cv1= " <<cv_equal_cv1 << " cv_equal_cv2= " <<cv_equal_cv2 << std::endl;

		if (cv_equal_cv1 || cv_equal_cv2) {
			//std::cout << " cv_equal_cv1= " <<cv_equal_cv1 << " cv_equal_cv2= " <<cv_equal_cv2 << std::endl;

			if (cv_equal_cv1) {
				if (traits->point_equal(p,prev->source()->point())) {
					out_edge = prev;
					out_vertex = prev->source();
					lt = Planar_map::VERTEX;
					found_vertex_or_edge = true;
					return; 
				}
				if (traits->point_in_x_range(cv1,p) && 
					traits->curve_compare_y_at_x(p,cv1) == EQUAL) {
						// p lies on cv1    
						out_edge = prev;
						lt = Planar_map::EDGE;
						found_vertex_or_edge = true;
						return; 
					}
					//p does not lie on cv1 ==> 
					// the target of the equal curve is a better vertex to p 
					std::cout << " WARNING: found closer vertex " <<std::endl;
					// out_vertex =  the closer vertex
					out_vertex = prev->source();
					new_vertex = true;
					std::cout << "The new vertex is: "<< out_vertex->point() << std::endl;
					// check validity (the new vertex is vetween them on a line) @@@@
					return;
			}
			else { //cv_equal_cv2
				if (traits->point_equal(p,first->source()->point())) {
					out_edge = first;
					out_vertex = first->source();
					lt = Planar_map::VERTEX;
					found_vertex_or_edge = true;
					return; 
				}
				if (traits->point_in_x_range(cv2,p) && 
					traits->curve_compare_y_at_x(p,cv2) == EQUAL) {
						// p lies on cv1    
						out_edge = first;
						lt = Planar_map::EDGE;
						found_vertex_or_edge = true;
						return; 
					}
					//p does not lie on cv1 ==> 
					// the target of the equal curve is a better vertex to p 
					std::cout << " WARNING: found closer vertex " <<std::endl;
					// out_vertex =  the closer vertex
					out_vertex = prev->source();
					new_vertex = true;
					std::cout << "The new vertex is: "<< out_vertex->point() << std::endl;
					// check validity (the new vertex is vetween them on a line) @@@@
					return;
			}

		}

		if (is_btw) {
			#ifdef CGAL_LM_DEBUG 
				std::cout << " cv is between " << cv1 << " and " << cv2 << std::endl;
			#endif //CGAL_LM_DEBUG
			out_edge = prev;
			found_face = true;
			//std::cout << " before return from find_face. out_edge =   " << out_edge->curve() << std::endl;
			return;
		}

	}

	std::cerr << "ERROR 4: edge above cw not found !" <<std::endl;
}	

//----------------------------------------------------
/*!
* Checks if there is an intersection between seg and e
* 
* \param p - the input point.
* \param seg - the segment
* \param e - the input edge
* \param num_of_intersections - output: 
*                               number of intersections between p-vp and e
* \param change_side - did e change side from p to vp
* \param found_edge - did we found the edge p lies on
* \param closest_interect_point - output: the closest intersection point to p
* \param new_vertex - output bool, if a closer vertex to p was found
* \param out_vertex - the output vertex (if closer vertex found)
*/
template <class Planar_map, class Nearest_neighbor>
void Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_intersection (const Point_2 & p,                  //input seg src 
				   Curve_2 &seg,                   //input seg
				   Halfedge_handle e,                  //input curve
				   int  & num_of_intersections,        //out num intersections
				   bool & change_side,                 //out did we move side?
				   bool & found_edge,                  //out is p on e?
				   Point_2 & closest_interect_point,   //out return value 
				   bool & new_vertex,                  //out if closest to p?
				   Vertex_handle & out_vertex ) const  //out the closer vertex
{
	//the idea is that the function checks whether the segment seg intersects e.
	//if it finds closer vertex to p than v, then it returns the new vertex in out_vertex, 
	//   and set the flag new_vertex to true. This can be the case if seg  intersects e in its end point. 
	//   or if there is an overlap between seg and e.
	//if p is on e, then it changes the flag found_edge to be true. 
	//if there is no intersection it returns all flags to be false, and num_of_intersections to be 0.
	//if there is a regular intersection (in a point), it returns the number of intersection (usually 1). 
	//   since this could be a specail point (like a tangent point), it also returns the flag change_side to check
	//   whether the point p is on a different side of e than seg->source() is.
	//closest_intersect_point is the intersection point between seg and e. 
	//   if there is more than one intersection, this is the closest one to p.
	entries_to_fi++;

	#ifdef CGAL_LM_DEBUG 
		std::cout << " inside find_intersection. p = " << p << ", seg = "<< seg 	<<", e = "<< e->curve() <<std::endl;
	#endif //CGAL_LM_DEBUG
	change_side = false;
	num_of_intersections = 0;
	found_edge = false;
	new_vertex = false;

	bool define_max = false;
	int max_intersections = 0;
	int intersection_counter = 0; 

	//if the define  MAX_INTERSECTIONS_BETWEEN_2_CURVES is on, then we know not to 
	//check for more intersections after the first X where found (where X is max_intersections)
#ifdef MAX_INTERSECTIONS_BETWEEN_2_CURVES
	max_intersections = MAX_INTERSECTIONS_BETWEEN_2_CURVES;
	define_max = true;
#else 
	std::cout << " MAX_INTERSECTIONS_BETWEEN_2_CURVES is not defined" <<std::endl;
#endif

	//Curve_2 seg(vh->point(), p); //seg in the reffered segment vh - p.
	Curve_2 seg1, seg2;
	Curve_2 cv = e->curve();
	//bool is_vertical = traits->curve_is_vertical(seg);  ???

	//loop on all intersections of seg and e
	do {
		CGAL::Object res_obj;
		Point_2 inter_point;
		Curve_2 overlap_seg;
		Point_2 vp = traits->curve_source(seg);
		Comparison_result comp_xy_res = traits->compare_xy(p, vp);

		#ifdef LM_CLOCK_DEBUG			
			clock_t ni_time_start = clock();//time start
		#endif

		switch (comp_xy_res) {
		case LARGER: //p is on the right of vp (or vertical and up)
			//find nearest intersection to the right of vp
			//std::cout << "check nearest intersection to the right of : "<< vp << "seg = "<<seg<< " cv = "<< cv <<std::endl;
			res_obj = traits->nearest_intersection_to_right(seg, cv, vp);
			break;
		case SMALLER: //p is on the left of vp (or vertical and down)
			//find nearest intersection to the right of p
			//std::cout << "check nearest intersection to the right of : "<< p <<" seg = "<<seg<< " cv = "<< cv <<std::endl;
			//res_obj = traits->nearest_intersection_to_right(seg, cv, p); 
			//std::cout << "check nearest intersection to the left of : "<< vp <<" seg = "<<seg<< " cv = "<< cv <<std::endl;
			res_obj = traits->nearest_intersection_to_left(seg, cv, vp); 
			break;
		default: //should not be equal
			CGAL_assertion (false);
		}

		#ifdef LM_CLOCK_DEBUG			
			clock_t ni_time_end = clock();//time end
			double ni_period = (double) (ni_time_end - ni_time_start);
			clock_ni += ni_period;
		#endif

		// Empty object is returned - no intersection.
		if (res_obj.is_empty()) {
			//if (comp_xy_res == SMALLER) {
			//	if (traits->is_on_segment (cv, p)) {
			//		closest_interect_point = p;
			//		found_edge = true;
			//		std::cout << "p is on edge" <<std::endl;
			//		return;
			//	}
			//}
			#ifdef CGAL_LM_DEBUG 
				std::cout << " no intersection" <<std::endl;
			#endif //CGAL_LM_DEBUG
			return;
		}

		// Intersection is a point
		else if (assign(inter_point, res_obj)) {
			#ifdef CGAL_LM_DEBUG 
				std::cout << " intersection is a point = "<< inter_point <<std::endl;
			#endif //CGAL_LM_DEBUG

			//if the intersection point is p, then p is on cv - found
			if (traits->point_equal(p, inter_point)) {
				closest_interect_point = inter_point;
				found_edge = true;
				//std::cout << " intersection equals p (p is on edge)" <<std::endl;
				return;
			} 
			//if the intersection point is vh->point, ignore it ! 
			if (traits->point_equal(seg.source(), inter_point)) {
				//std::cout << " intersection equals seg->source() = " << seg.source()  <<std::endl;
				return;
			}

			//check if the intersection point is an end point of e
			if (traits->point_equal(e->source()->point(), inter_point)) {
				new_vertex = true;
				out_vertex = e->source();
				//std::cout << " intersection found new vertex: e->source =  " << e->source()->point() <<std::endl;
				return;
			}
			if (traits->point_equal(e->target()->point(), inter_point)) {
				new_vertex = true;
				out_vertex = e->target();
				//std::cout << " intersection found new vertex: e->target =  " << e->target()->point()<<std::endl;
				return;
			}
			//else - regular intersection

			if ((!define_max) || (max_intersections > 1)) {
				std::cout << " warning: should not reach here if segments" <<std::endl;
				// @@@@ check if this is not a tangent point !
				// @@@@ to check tangent point need to do curves_compare_y_at_x_left and right (
				// @@@@ if tangent point - don't change side.
			}

			//if not first intersection - compare inter_point with closest_interect_point (ref p)
			//if inter_point closer - change it .
			if ((num_of_intersections == 0) || 
				(traits->compare_distance(p, inter_point, closest_interect_point) 
				== SMALLER)) {
					closest_interect_point = inter_point;

					// split seg and save the part of seg closer to p
					traits->curve_split(seg, seg1, seg2, inter_point);
					Point_2 seg1_src = traits->curve_source(seg1);
					Point_2 seg2_trg = traits->curve_target(seg2);
					if ((traits->compare_xy(seg1_src, seg2_trg)) == comp_xy_res) {
						seg = seg1; 
					}
					else { 
						seg = seg2;
					}
					//std::cout << "seg was split inside find_intersect. new seg = " <<seg << std::endl;
				}
				change_side = change_side ? false : true; 
				num_of_intersections++;	
		}

		// Intersection is a segment
		else if (assign(overlap_seg, res_obj))
		{
			#ifdef CGAL_LM_DEBUG 
				std::cout << "intersection is an overlapped segment" << std::endl;
			#endif //CGAL_LM_DEBUG

			Point_2 ov_src = traits->curve_source(overlap_seg);
			Point_2 ov_trg = traits->curve_target(overlap_seg);
			//if p equals on of the end points of the overlapped segment, then p is on e.
			if (traits->point_equal(p, ov_src) || traits->point_equal(p, ov_trg)) {
				found_edge = true;
				return;
			}
			//if the one of the end points of e equals on of the  overlapped segment's end points, 
			//then a new vertex is found (the endpoint of e)
			if (traits->point_equal(e->source()->point(), ov_src) || 
				traits->point_equal(e->source()->point(), ov_trg)) 
			{  //e->source is equal
				new_vertex = true;
				out_vertex = e->source();
				return;
			} 
			if (traits->point_equal(e->target()->point(), ov_src) ||
				traits->point_equal(e->target()->point(), ov_trg)) 
			{ //e->source is equal
				new_vertex = true;
				out_vertex = e->target();
				return;
			}

			std::cout << "WARNING: this case was not implemented !!!" <<std::endl;
			//getchar();
			// @@@@ - else, relevant only if not segments 
			//check if ov_src == vp or ov_trg == vp. if so ???@@@@
			//      - we need to compare cw vp-p and right_inter_point v.s. the rest of the curve
			//else if not equal vp - then need to compare from right and from left 
			//  of this curve (this is to find if tangent. this is only relevant if 
			//  ((!define_max) || (max_intersections > 1)). otherwise, its irrelevant
		}   

		else {// We should never reach here:
			CGAL_assertion (false);
		}

		// go on with the loop
		++intersection_counter;
	} while ((num_of_intersections > 0) && 
		((!define_max) || 
		(define_max && (intersection_counter < max_intersections)) ) );

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

template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
is_point_in_face (const Point_2 & p,                  //input seg src 
				const Ccb_halfedge_circulator & face,               //input face
				bool & found_edge,                  //out is p on e?
				bool & found_vertex,               //out is p equals e->target?
				Halfedge_handle  & out_edge) const              //output if on curve				
{
	#ifdef CGAL_LM_DEBUG
		std::cout << "inside is_point_in_face. face = " << face->source()->point()
			              <<"-->" << face->target()->point()<<std::endl;
	#endif //CGAL_LM_DEBUG

	found_edge = false;
	found_vertex = false;
	Comparison_result	  compare_y_at_x_res; 
	int number_of_edges_above_p = 0;

	//if  (face->is_unbounded())
	//	return (true);

	//loop on all edges in this face
	Ccb_halfedge_circulator curr = face;
	Ccb_halfedge_circulator last =  curr;

	do  {
		Curve_2 cv = curr->curve();
		Point_2 p1 = curr->source()->point();
		Point_2 p2 = curr->target()->point();

		//check if p equals one of the endpoints of e
		if (traits->point_equal(p, p1))   {
			found_vertex = true;
			out_edge = curr->twin() ;
			return (true); 
		}
		if (traits->point_equal(p, p2))   {
			found_vertex = true;
			out_edge = curr ;
			return (true); 
		}

		//check in_x_range lexicographically. This is if p is on different sides from p1 and p2
		if  (traits->point_is_left_low(p, p1) != traits->point_is_left_low(p, p2) ) 
		{
			//check cv to see it p is above, below or on cv
			compare_y_at_x_res =  traits->curve_compare_y_at_x(p, cv);
			switch (compare_y_at_x_res) {
		case EQUAL:
			found_edge = true;
			out_edge = curr;
			break; 
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
	//to check if number_of_edges_above_p is odd/even, simply do bitwise AND operator with 1. 
	return (number_of_edges_above_p & 1) ;
}

//----------------------------------------------------
/*!
* Finds the intersection curve between the segment v - p and the face. (out_edge)
* if there is more than one intersection - return the closest edge to p.
* returns false if there is no intersection.
* 
* \param p - the input point.
* \param v - the input vertex
* \param face - the input face
* \param out_edge - the output edge
*/
template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_closest_intersection_in_face (const Point_2 & p,                  //input seg src 
								   	Vertex_handle  v,           //input vertex
								  const Ccb_halfedge_circulator & face,                                  //input face
								  Halfedge_handle  & out_edge) const              //output if on curve				
{
	#ifdef CGAL_LM_DEBUG
	std::cout << "inside find_closest_intersection_in_face.   " << std::endl;
	#endif //CGAL_LM_DEBUG

	Point_2 vp = v->point();
	Curve_2 seg(vp, p);   //create a segment vh--p. 

	bool is_found_intersection = false;
	bool is_inter_point_updated = false;
	Point_2 inter_point;

	//loop on all edges in this face
	Ccb_halfedge_circulator curr = face;
	Ccb_halfedge_circulator last =  curr;
	bool p_in_x_range, v_in_x_range, p1_in_x_range, p2_in_x_range;	
	bool check_real_intersection; 

	do  {
		Curve_2 cv = curr->curve();
		Point_2 p1 = curr->source()->point();
		Point_2 p2 = curr->target()->point();
		#ifdef CGAL_LM_DEBUG
			std::cout << "curr = " << p1 << "-->" << p2 <<std::endl;
		#endif //CGAL_LM_DEBUG

		//@@@@ check if these operation take long time
		p_in_x_range = traits->point_in_x_range(cv, p);
		v_in_x_range = traits->point_in_x_range(cv, vp);
		p1_in_x_range = traits->point_in_x_range(seg, p1);
		p2_in_x_range = traits->point_in_x_range(seg, p2);
		//
		#ifdef CGAL_LM_DEBUG
			std::cout << "p_in_x_range = " << p_in_x_range << " , v_in_x_range = " << v_in_x_range <<std::endl;
			std::cout << "p1_in_x_range = " << p1_in_x_range << " , p2_in_x_range = " << p2_in_x_range <<std::endl;
		#endif //CGAL_LM_DEBUG

		if ( p_in_x_range ||  v_in_x_range || p1_in_x_range || p2_in_x_range)
		{
			//check_real_intersection = false;
			//if (p_in_x_range && v_in_x_range) { 
			//	//check if on different sides - there is intersection - what if v is on cv ? 
			//	Comparison_result res_p = traits->curve_compare_y_at_x(p,cv); 
			//	Comparison_result res_v = traits->curve_compare_y_at_x(vp,cv); 
			//	if (res_p == EQUAL) {
			//		std::cerr << "ERROR 5: p should not be on any edge" << std::endl;
			//	}
			//	 if (res_p != res_v && res_v != EQUAL) { 
			//		 if (! is_found_intersection) { //this is the first intersection
			//			is_found_intersection = true;
			//			out_edge = curr;
			//		 }
			//		 else { //this is not the first inetrsection - need to check who is closer (real intersection)
			//			 check_real_intersection = true;
			//		 }
			//	}
			//}
			//else { //need to check real intersection
			//	check_real_intersection = true;
			//}

			//if (check_real_intersection) 
			//{
				Point_2 temp_inter_point;
				//check intersection and update temp_inter_point
				if ( find_real_intersection(p, v, curr, temp_inter_point) ) 
				{
					#ifdef CGAL_LM_DEBUG
						std::cout << "find_real_intersection returned true "  <<std::endl;
					#endif //CGAL_LM_DEBUG

					if (! is_found_intersection) { //this is the first intersection
						is_found_intersection = true;
						is_inter_point_updated = true;
						inter_point = temp_inter_point;
						out_edge = curr;
					}
					else { //compare to previous intersection
						if (!is_inter_point_updated)  {
							//check intersection between out_edge and v-p and update inter_point
							if ( ! find_real_intersection(p, v, out_edge, inter_point) )
								std::cerr << "ERROR 6: not found intersection when exists" << std::endl;
						}
						//compare which point is closer. 
						if (traits->compare_distance(p, temp_inter_point, inter_point) 
							== SMALLER) { // the temp_inter point is closer to p
								inter_point = temp_inter_point;
								out_edge = curr;
							}
							is_inter_point_updated = true;
					}
				}
				else {
					#ifdef CGAL_LM_DEBUG
						std::cout << "find_real_intersection returned false "  <<std::endl;
					#endif //CGAL_LM_DEBUG
				}
			//}
		}

		++curr;
	} while (curr != last);  

	return (is_found_intersection);
}

//----------------------------------------------------
/*!
* Finds the intersection curve between the segment v - p and the face. (out_edge)
* if there is more than one intersection - return the closest edge to p.
* returns false if there is no intersection.
* 
* \param p - the input point.
* \param v - the input vertex
* \param e - the input curve
* \param out_point - the output point
*/
template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_real_intersection (const Point_2 & p,                  //input seg src 
						Vertex_handle  v,                        //input vertex
						Halfedge_handle e,                  //input curve
						Point_2 & out_point) const      //output intersection point						
{
	bool new_vertex, found_edge, change_side;
	Curve_2 seg(v->point(), p); //seg in the reffered segment vh - p.
	int num_intersections; 
	Vertex_handle temp_vertex; 

	#ifdef LM_CLOCK_DEBUG				
		clock_t fi_time_start = clock();//time start
	#endif

	find_intersection (p, seg, e, num_intersections, change_side, found_edge, 
		out_point, new_vertex, temp_vertex);

	#ifdef LM_CLOCK_DEBUG		
		clock_t fi_time_end = clock();//time end
		double fi_period = (double) (fi_time_end - fi_time_start);
		clock_fi += fi_period;
	#endif

	//check results
	if (new_vertex) { 
		std::cerr << "ERROR 7: new vertex in find_real_intersection " << std::endl;
		//ixx: @@@@ maybe later we will also have to deal with this case
		return (false);
	}
	if (found_edge) {
		std::cerr << "ERROR 8: new edge in find_real_intersection " << std::endl;
		//ixx: @@@@ maybe later we will also have to deal with this case
		return (false);
	}

	return (change_side) ;
}

CGAL_END_NAMESPACE

#endif  //CGAL_PM_LANDMARKS_POINT_LOCATION_C
