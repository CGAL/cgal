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
 * Member-function definitions for the 
 * Arr_landmarks_point_location<Arrangement> class.
 */

//#define CGAL_LM_DEBUG
#ifdef CGAL_LM_DEBUG
  #define PRINT_DEBUG(expr)   std::cout << expr << std::endl
  #define LM_DEBUG(cmd)   cmd
#else
  #define PRINT_DEBUG(expr)
  #define LM_DEBUG(cmd) 
#endif

//#define PRINT_ERROR(expr)   std::cerr << expr << std::endl
#define PRINT_ERROR(expr)   std::cout << expr << std::endl

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
    return (CGAL::make_object (p_arr->unbounded_face()));
  
  Object  lm_location_obj; 
  Point_2 landmark_point = lm_gen->get_closest_landmark (p, 
               lm_location_obj);

  PRINT_DEBUG("test ") ;

  //NN_Point_2 nearest_landmark = lm_gen->find_closest_landmark(p);
  //Point_2 landmark_point = nearest_landmarks.get_point();
  //Object  lm_location_obj = nearest_landmarks.get_obj() ;
  PRINT_DEBUG("nearest neighbor of point "<< p << " is " << landmark_point);
  
  //walk from the nearest_vertex to the point p, using walk algorithm, 
  //and find the location of p.   
  Object  out_obj; //the output object
  
  //if the landmark s not found in the arangement
  const Vertex_const_handle     *vh;
  const Halfedge_const_handle   *hh;
  const Face_const_handle       *fh;
  m_start_edge = NULL;
  
  if (lm_location_obj.is_empty())
  {
    PRINT_ERROR( "lm_location_obj is empty" );
    CGAL_assertion (false);
    return out_obj;
  }
  else if ((vh = object_cast<Vertex_const_handle>(&lm_location_obj)) != NULL)
  {
    PRINT_DEBUG( "lm_location_obj is a vertex: "<< (*vh)->point());
    out_obj = _walk_from_vertex (*vh, p);
  }
  else if ((fh = object_cast<Face_const_handle>(&lm_location_obj)) != NULL)
  {
    PRINT_DEBUG( "lm_location_obj is a face. ");
    out_obj = _walk_from_face (*fh, p, landmark_point);
  }
  else if ((hh = object_cast<Halfedge_const_handle>(&lm_location_obj)) != NULL)
  {
    PRINT_DEBUG( "lm_location_obj is a halfedge: "<< (*hh)->curve());
    out_obj = _walk_from_edge (*hh, p, landmark_point);
  }
  else 
  {
    PRINT_ERROR( "unknown object");
    CGAL_assertion (false);
    return out_obj;
  }
  
  PRINT_DEBUG( "return from walk" << std::endl);
  
#ifdef CGAL_LM_DEBUG
  if (out_obj.is_empty())
  {
    PRINT_ERROR( "object is empty" );
    CGAL_assertion (false);
  }
  else if ((hh = object_cast<Halfedge_const_handle>(&out_obj)) != NULL)
  {
    PRINT_DEBUG( "object is a halfedge: "<< (*hh)->curve());
  }
  else if ((vh = object_cast<Vertex_const_handle>(&out_obj)) != NULL)
  {
    PRINT_DEBUG( "object is a vertex: "<< (*vh)->point());
  }
  else if ((fh = object_cast<Face_const_handle>(&out_obj)) != NULL)
  {
    PRINT_DEBUG( "object is a face. ");
  }
#endif
  
  if ((fh = object_cast<Face_const_handle>(&out_obj)) != NULL)
  {
    // If we reached here, we did not locate the query point in any of the
    // holes inside the current face, so we conclude it is contained in this
    // face.
    // However, we first have to check whether the query point coincides with
    // any of the isolated vertices contained inside this face.
    Isolated_vertices_const_iterator   iso_verts_it;
    typename Traits_wrapper_2::Equal_2 equal = traits->equal_2_object();

    for (iso_verts_it = (*fh)->isolated_vertices_begin();
         iso_verts_it != (*fh)->isolated_vertices_end(); ++iso_verts_it)
    {
      if (equal (p, iso_verts_it->point()))
      {
        Vertex_const_handle  vh = iso_verts_it;
        return (CGAL::make_object (vh));
      }
    }    
  }

  return (out_obj);
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
  PRINT_DEBUG("inside walk_from_vertex. p= "<< p  << 
        ", nearest_vertex = "<<nearest_vertex->point() );
  
  //inits
  Vertex_const_handle       vh = nearest_vertex;
  
  if (vh->is_isolated())
  {
    Face_const_handle f = vh->face();
    return _walk_from_face(f, p, vh->point());
  }

  //find face
  bool                      new_vertex = false;
  Object                    obj;
  const Face_const_handle  *p_fh;

  do
  {
    //find the edge out_edge which is the best possibly 
    //pointing to the face containing p

    m_flipped_edges.clear();  //clear the curves that were flipped
    
    new_vertex = false;
    obj = _find_face (p, vh, new_vertex);

    if (new_vertex)
    {
      PRINT_DEBUG( "NEW vertex 1 " );
      //check if the new vertex is really closer 
      // I removed the check if the vertex is closer since there is no 
      // compare distance 
      //if (traits->compare_distance(p, out_vertex->point(), vh->point())  
      //  != SMALLER) {PRINT_DEBUG("Error 2: new vertex"); return; }
      vh = object_cast<Vertex_const_handle> (obj);
    }
    else if (obj.is_empty())
    {
      PRINT_ERROR( "object is empty" );
      CGAL_assertion (false);
      return obj;
    }
    else if (object_cast<Halfedge_const_handle>(&obj) != NULL)
    {
      PRINT_DEBUG ("_find_face found a halfedge: " << 
		   (object_cast<Halfedge_const_handle>(&obj))->curve());
      return (obj);
    }
    else if (object_cast<Vertex_const_handle>(&obj) != NULL)
    {
      PRINT_DEBUG ("_find_face found a vertex: " << 
		   (object_cast<Vertex_const_handle>(&obj))->point());
      return (obj);
    }
    else if ((p_fh = object_cast<Face_const_handle>(&obj)) != NULL)
    {
      PRINT_DEBUG ("_find_face found a face.");
      return _walk_from_face (*p_fh, p, vh->point());
    }
    
  } while (new_vertex);  
  
  // We should never reach here:
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
        bool & new_vertex) const
{ 
  PRINT_DEBUG("inside find_face. p ="<< p <<" , vh = "<<vh->point() ); 
  
  new_vertex = false;
  
  // check if the point equals the vertex. 
  if (traits->equal_2_object()(vh->point(), p))
  {
    return (CGAL::make_object (vh));
  }

  //create a segment vh--p. 
  const Point_2&     v = vh->point();
  X_monotone_curve_2 seg = traits->construct_x_monotone_curve_2_object()(v, p);
  
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
  Comparison_result cv_orient = circ->twin()->direction();
  // RWRW: compare_xy(v, (*circ).source()->point());
  
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
        cv_orient = circ->twin()->direction(); 
          // RWRW: compare_xy(v,(*circ).source()->point());
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
  return (CGAL::make_object ((*prev).face()));
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
          cv_orient = circ->twin()->direction();
            //RWRW: compare_xy(v,(*circ).source()->point());
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
        return (CGAL::make_object ((*prev).face()));
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
          cv_orient = circ->twin()->direction();
          // RWRW: compare_xy(v,(*circ).source()->point());
        } while (circ!=circ_done  && seg_orient == cv_orient); 
        
        if (seg_orient != cv_orient) 
        {
          //loop 2 - until reach the same side again 
          do {
            prev = circ; 
            circ++;
            cv_orient = circ->twin()->direction();
            // RWRW: compare_xy(v,(*circ).source()->point());
          }  while (seg_orient != cv_orient) ;
          
          //prev is the last edge from the other side, 
          //and curve is the first curve in this side
          res1 =  compare_cw_around_point(seg, (*circ).curve(), v);
          
          //loop 3 - until find curve > seg cw.
          while (res1 == LARGER && seg_orient==cv_orient) 
          {
            prev = circ; 
            circ++;
            cv_orient = circ->twin()->direction();
            // RWRW: compare_xy(v,(*circ).source()->point());
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
          return (CGAL::make_object ((*prev).face()));
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
    return (CGAL::make_object ((*circ).twin()->face()));
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
  return (CGAL::make_object ((*prev).face()));
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
    return (CGAL::make_object ((*circ).source()));
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
    return (CGAL::make_object (temp_he)); 
  }
  
  //p does not lie on cv1 ==> 
  // the target of the equal curve is a better vertex to p 
  PRINT_ERROR("WARNING 11: found closer vertex during new_find_face");
  // out_vertex is the closer vertex
  //out_vertex = circ->source();
  new_vertex = true;
  PRINT_DEBUG( "The new vertex is: "<< (*circ).source()->point() );
  // check validity (the new vertex is vetween them on a line) @@@@
  return (CGAL::make_object((*circ).source()));    
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
  const X_monotone_curve_2& cv = eh->curve() ;
  const Point_2&            src = eh->source()->point();
  const Point_2&            trg = eh->target()->point();
  Comparison_result res;

  LM_DEBUG( 
    if (! traits->is_in_x_range_2_object()(cv, np)) 
      std::cout<<"WARNING 5: np is not on the edge's x_range"
      <<std::endl; 
    else if (traits->compare_y_at_x_2_object()(np,cv) != EQUAL)
      std::cout<<"WARNING 6: np is not on the edge"
      <<std::endl; 
  );
  PRINT_DEBUG("inside walk_from_edge. p= "<< p  << 
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
    return (CGAL::make_object(vh));
  }
  if (traits->equal_2_object()(p, trg))
  {
    Vertex_const_handle vh = eh->target();
    return (CGAL::make_object(vh));
  }
  
  //save the edge we're starting from
  m_start_edge = &eh;

  //if p is in eh's x_range, then we need to check if it is above/below eh
  //and orient the halfedge eh accordingly, so that it will point to the face 
  //that is most likely containing p
  if (traits->is_in_x_range_2_object()(cv, p))
  {
    //check if p is above/below cv
    res =   traits->compare_y_at_x_2_object()(p,cv);
    PRINT_DEBUG("curve compare y at x: p= "<< p << ", cv =  "<< cv 
          <<", res = "<<res);
    switch (res) { 
      case EQUAL://p is on cv - found !
        return (CGAL::make_object(eh));
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
    if (eh->direction() != traits->compare_xy_2_object()(p, src))
      //RWRW: if (traits->compare_xy_2_object()(src, trg) != 
      //          traits->compare_xy_2_object()(p, src)) 
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
//   yes: 
//     go over holes. for each hole h:
//       is p in h?
//         yes: face = h. goto 1.
//   no: 
//     call new function: 
//       find_closest_intersection_in_face( gets face, v, p) 
//                 that  returns e and intersection point (if needed). 
//       (the function go over all edges surrounding p. 
//           for each one- as done now - take out of the walk procedure) 
//             face = e->twin. goto1.
//       if intersection not found --- ? error ? 
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
  m_flipped_edges.clear();  //remove all elements from the m_flipped_edges list
  bool p_in_face = false;
  Ccb_halfedge_const_circulator h_circ;
  bool found_edge = false;
  bool found_vertex = false;
  Halfedge_const_handle out_edge;
  
  do {
    PRINT_DEBUG(std::endl << "inside loop on face ");
    p_in_face = false;
    if (face->is_unbounded())  {
      p_in_face = true;
      PRINT_DEBUG("unbounded face ");
    }
    else {    
      h_circ = face->outer_ccb();
      p_in_face = _is_point_in_face(p, h_circ, found_edge, 
        found_vertex, out_edge); 
      PRINT_DEBUG("is_point_in_face returned  "<<  p_in_face );
    }
    if (found_vertex)
    {
      Vertex_const_handle v = out_edge->target(); 
      return (CGAL::make_object(v)); //is it really the target?
    }
    else if (found_edge) 
    {
      Halfedge_const_handle h = out_edge; 
      return (CGAL::make_object(h));
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
        p_in_hole = _is_point_in_face(p, *hole_it, 
          found_edge, found_vertex, out_edge); 
        if (found_vertex)
        {
          return (CGAL::make_object(out_edge->target())); 
          //is it really the target?
        }
        else if (found_edge) 
        {
          return (CGAL::make_object(out_edge));
        }
        else if (p_in_hole) {
          h_circ = *hole_it; 
          //update the new "face " to be the hole. check its holes.      

          if ( _find_edge_to_flip (p, np, h_circ , out_edge) ) {
            out_edge = out_edge->twin();
            face = out_edge->face();
          }
          else {
            PRINT_ERROR( "ERROR 10:  intersection not found");
            LM_DEBUG(getchar());
            return (CGAL::make_object (p_arr->unbounded_face()));
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
      if ( _find_edge_to_flip (p, np, face->outer_ccb() , out_edge) ) {
        out_edge = out_edge->twin();        
        face = out_edge->face();
        PRINT_DEBUG("after find_edge_to_flip. changed to twin @ " );
        PRINT_DEBUG("out_edge  =  "<< out_edge->source()->point() 
                  <<" -->"<< out_edge->target()->point() );
      }
      else {
        PRINT_ERROR("ERROR 9:  intersection not found");
        LM_DEBUG(getchar());
        //e = p_arr->halfedges_end();
        //lt = Planar_map::UNBOUNDED_FACE;
        return (CGAL::make_object (p_arr->unbounded_face()));
      }
    }
  }while (!p_in_face); 

  if (face == p_arr->unbounded_face()) 
  {
    PRINT_DEBUG("before return from walk from face. unbounded face.");
    return (CGAL::make_object (p_arr->unbounded_face()));
  }

  PRINT_DEBUG("before return from walk from face. ");
  return (CGAL::make_object (face));
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
    const X_monotone_curve_2& cv = (*curr).curve();
    const Point_2&            p1 = (*curr).source()->point();
    const Point_2&            p2 = (*curr).target()->point();

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
  PRINT_DEBUG("inside find_edge_to_flip.   " );

  //create a segment vp--p. 
  const Point_2&     vp = np;
  X_monotone_curve_2 seg = traits->construct_x_monotone_curve_2_object()(vp, p);

  //loop on all edges in this face
  Ccb_halfedge_const_circulator curr = face;
  Ccb_halfedge_const_circulator last =  curr;
  bool p_in_x_range, v_in_x_range, p1_in_x_range, p2_in_x_range;  

  typename Traits_wrapper_2::Is_in_x_range_2         is_in_x_range = 
                                      traits->is_in_x_range_2_object();


  do  {
    const X_monotone_curve_2& cv = (*curr).curve();
    const Point_2&            p1 = (*curr).source()->point();
    const Point_2&            p2 = (*curr).target()->point();
 
    // check if the curve was already flipped - in this case, 
    //    don't check it at all.
    Std_edge_iterator found1 = std::find (m_flipped_edges.begin(), 
      m_flipped_edges.end(), curr);
    Std_edge_iterator found2 = std::find (m_flipped_edges.begin(), 
      m_flipped_edges.end(), (*curr).twin());
    if (found1 != m_flipped_edges.end() || found2 != m_flipped_edges.end()) {
      PRINT_DEBUG("curve "<<p1<<"-->"<<p2<<" was found in the list");
    }
    else if (m_start_edge && 
             (curr == *m_start_edge || (*curr).twin() == *m_start_edge))
    { //check that curr and (*curr).twin() is not equal to m_start_edge
      PRINT_DEBUG("curve "<<p1<<"-->"<<p2<<" is the start edge");
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
        bool check_res = _check_approximate_intersection (seg, cv, intersect);

        if (check_res) 
        {
          if (intersect) {
            out_edge = curr;
            m_flipped_edges.push_back(out_edge);
            return (true);
          }
        }
        else {
          PRINT_ERROR("ERROR 12: check_approximate_intersection "\
                      <<"did not return an answer.");
        }
      }
      //else - there is not intersection.
    }  

    ++curr;
  } while (curr != last);  

  return (false);
}


//----------------------------------------------------
//check if cv intersects seg. 
// \param seg - the input seg.
// \param cv - the input curve
// \param intersect - output: if the two curves intersects odd number of times
// \returns true if the function returned a value, false if error
// if the segment intersects the curve in one of the curves endpoints: 
// if it is the curve's left endpoint => intersects. 
// the curves right endpoint => no intersection.
// if the segment intersect the curve in one of the segments end points, 
// it may because the query point is on the curve (this should not happen, 
// since we check first if the query is on an edge), or 
// beacuse the landmark is on the curve, and in this case we should not pass 
// to the edge's other side.
template <class Arrangement_2, class Arr_landmarks_generator>
bool Arr_landmarks_point_location<Arrangement_2,Arr_landmarks_generator>
::_check_approximate_intersection (const X_monotone_curve_2 & seg,  
                  const X_monotone_curve_2 & cv,
                  bool & intersect) const 
{
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

  const Point_2& seg_right = traits->construct_max_vertex_2_object()(seg);
  const Point_2& seg_left = traits->construct_min_vertex_2_object()(seg);
  const Point_2& cv_right = traits->construct_max_vertex_2_object()(cv);
  const Point_2& cv_left = traits->construct_min_vertex_2_object()(cv);
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
    PRINT_DEBUG("the left end-point of the segments is equal to the"\
                <<" left end-point of the curve");
    // compare to the right of the curves
    Comparison_result curves_comp = compare_y_at_x_right(seg,cv,seg_left);
  
    if (curves_comp == EQUAL) 
    {
      PRINT_DEBUG("overlap !!!");
      //this means: 
      // 1. the query is on the curve (should have been checked). 
      // 2. the landmark is on the curve (should have been checked).
      // 3. the curve endpoint is on the segment. if it's the left endpoint - 
      //    intersect. right - no intersection.
      if (compare_y_at_x(cv_left, seg) == EQUAL)
      {
        intersect = true;
        return (true);
      }
      else if (compare_y_at_x(cv_right, seg) == EQUAL)
      {
       //intersect = false;
        return (true);
      }
      return (false); //V
    }                                                                         
    if (comp_right_xy_res == SMALLER) 
    {  //the segments ends before (to the right) of the curve 
      Comparison_result curve_comp_right_res = 
        compare_y_at_x(seg_right,cv);

      if (curve_comp_right_res == EQUAL) {
        PRINT_DEBUG("2 points collide");
        //this means that the query point is on the curve
        return (false); //V
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
        //this means that the curve's right enpoint is on the segment
        //intersect = false;
        return (true);
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
        //this means that the query point is on the curve
        //(even on the curve's endpoint = vertex)
        return (false); //V
    }
  }

  // the right end-point of the segments is equal to the right end-point 
  //of the curve
  else if (comp_right_xy_res == EQUAL) 
  {
    PRINT_DEBUG("the right end-point of the segments is equal to the"\
                <<" right end-point of the curve");
    // compare to the left of the curves
    Comparison_result curves_comp = compare_y_at_x_left(seg, cv, seg_right);

    if (curves_comp == EQUAL) 
    {
      PRINT_DEBUG("overlap !!!");
      //this means: 
      // 1. the query is on the curve (should have been checked). 
      // 2. the landmark is on the curve (should have been checked).
      // 3. the curve endpoint is on the segment. if it's the left endpoint - 
      //    intersect. right - no intersection.
      if (compare_y_at_x(cv_left, seg) == EQUAL)
      {
        intersect = true;
        return (true);
      }
      else if (compare_y_at_x(cv_right, seg) == EQUAL)
      {
       //intersect = false;
        return (true);
      }
      return (false); //V
    }
    if (comp_left_xy_res == SMALLER) 
    { //the curve ends before (to the left) of the segment 
      Comparison_result curve_comp_left_res = 
        compare_y_at_x(cv_left,seg) ;

      if (curve_comp_left_res == EQUAL) {
        PRINT_DEBUG("2 points collide");
        //this means that the curve's left enpoint is on the segment
        intersect = true;
        return (true);
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
        //this means that the query point is the curve
        return (false); //V
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
        return (false); //V
    }
  }

  //one curve is inside the other curve's x-range (fully) 
  else if (comp_left_xy_res != comp_right_xy_res) { 
    if (comp_left_xy_res == LARGER) 
    {  //the segment is inside the curve's x-range
      PRINT_DEBUG("the segment is inside the curve's x-range.");
      LM_DEBUG ( //checks for debugging - remove later
        if (! is_in_x_range(cv,seg_left) ) {
          PRINT_ERROR( "! is_in_x_range(cv,seg_left) "); 
          return (false); } //V-debugging
        if (! is_in_x_range(cv,seg_right) )  {
          PRINT_ERROR("(! is_in_x_range(cv,seg_right) "); 
          return (false); } //V-debugging
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
      {  //no intersection
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
    {  //the curve is inside the segments's x-range
      PRINT_DEBUG("the curve is inside the segments's x-range");
      LM_DEBUG ( //checks for debugging - remove later
        if (! is_in_x_range(seg,cv_left) ) {
         PRINT_ERROR( "! is_in_x_range(seg,cv_left)"); 
         return (false);} //V-debugging
        if (! is_in_x_range(seg,cv_right) ) {
         PRINT_ERROR("! is_in_x_range(seg,cv_right)"); 
         return (false);} //V-debugging
      )

      Comparison_result curve_comp_left_res = compare_y_at_x(cv_left,seg);
      Comparison_result curve_comp_right_res = compare_y_at_x(cv_right,seg);

      if (curve_comp_right_res == EQUAL)
      {
        PRINT_DEBUG("the curve right endpoint is on the segment");
        PRINT_DEBUG("no intersection");
        //intersect = false;
        return (true);
      }
      else if (curve_comp_left_res == EQUAL) 
      {
        PRINT_DEBUG("the curve left endpoint is on the segment");
        PRINT_DEBUG("intersecting");
        intersect = true;
        return (true);
      }
      else if  (curve_comp_left_res == curve_comp_right_res) 
      {  //no intersection
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
      PRINT_DEBUG("the segment right point is in and the"\
                  <<"curve's left point");
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
          PRINT_DEBUG("segment's right endpoint = "\
                      <<"curve's left endpoint ");
          PRINT_DEBUG("no intersection");
          //intersect = false;
          return (true);    
        }
        if (curve_comp_left_res == EQUAL)
        {
          //the curve's left end-point is on the segment.
          PRINT_DEBUG("the curve's left end-point is on the segment");
          PRINT_DEBUG("intersecting");
          intersect = true;                    
          return (true);
        }
        //the segment's right endpoint is on the curve. 
        return (false); //V
      }
      if  (curve_comp_left_res == curve_comp_right_res) 
      {  //intersection
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
      PRINT_DEBUG("the curve's right point is in and the "\
                  <<"segment's left point");
      LM_DEBUG ( //checks for debugging - remove later
        if (! is_in_x_range(cv,seg_left) ) {
          PRINT_ERROR("! is_in_x_range(cv,seg_left) ***  "); 
          return (false);  } //V-debugging
        if (! is_in_x_range(seg,cv_right) ) {
          PRINT_ERROR("! is_in_x_range(seg,cv_right) *** "); 
          return (false); } //V-debugging
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
          PRINT_DEBUG("segment's left endpoint = "\
                      <<"curve's right endpoint ");
          PRINT_DEBUG("no intersection");
          //intersect = false;
          return (true);    
        }
        if (curve_comp_right_res == EQUAL)
        {
          PRINT_DEBUG("intersecting");
          intersect = true;
          return (true);        
        }
        //the segment's left endpoint is on the curve. 
        return (false); //V
      }
      if  (curve_comp_left_res == curve_comp_right_res) 
      {  //intersection
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
}

CGAL_END_NAMESPACE

#endif
