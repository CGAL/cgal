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

#ifdef LANDMARKS_CLOCK_DEBUG
  #define LM_CLOCK_DEBUG(cmd)   cmd
#else
  #define LM_CLOCK_DEBUG(cmd)
#endif

#ifdef TRAITS_CLOCK_DEBUG
  #define TR_CLOCK_DEBUG(cmd)   cmd
#else
  #define TR_CLOCK_DEBUG(cmd)
#endif

#ifdef CGAL_LM_DEBUG
  #define PRINT_DEBUG(expr)   std::cout << expr << std::endl
  #define LM_DEBUG(cmd)   cmd
#else
  #define PRINT_DEBUG(expr)
  #define LM_DEBUG(cmd) 
#endif

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
  //std::cout << "------ locate point "<< p <<std::endl;
  PRINT_DEBUG("------ locate point "<< p) ;

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
    //build the kd-tree for nearest neighbor search
    PRINT_DEBUG("landmarks tree is not updated. updating ... "<< p) ;
    LM_CLOCK_DEBUG(clock_t create_nn_time_start = clock());
    create_landmarks_tree();
    LM_CLOCK_DEBUG(clock_t create_nn_time_end = clock() );
    LM_CLOCK_DEBUG(clock_create_nn += (double) (create_nn_time_end - create_nn_time_start));
    //std::cout << "ERROR 1: tree is not updated ! number of vertices in pm = " << pm->number_of_vertices() << std::endl;
    //LM_DEBUG(getchar());
    //return e;
  }

  LM_CLOCK_DEBUG(clock_t nw_time_start = clock());

  //get nearest vertex to point p 
  LM_CLOCK_DEBUG(clock_t nn_time_start = clock());

  double px, py;
  point_to_double_coords(p, px, py);

  NN_Point_2  nnp (px, py); 
  NN_Point_2  nearest_point = nn.find_nearest_neighbor(nnp); 
  PRINT_DEBUG("nearest neighbor of point "<< p << " is " << nearest_point.x() <<','<<nearest_point.y()) ;
  Vertex_handle nearest_vertex = nearest_point.vertex();
  LM_CLOCK_DEBUG(clock_t nn_time_end = clock());
  LM_CLOCK_DEBUG(clock_for_nn_search += (double) (nn_time_end - nn_time_start) );

  //walk from the nearest_vertex to the point p, using walk algorithm, 
  //and find the location of p.
  //lt is the  location type as needed, and e is the halfedge pointing at the vertex or face 
  //that the point lies at, or if the point lies on an halfedge, its the halfedge itself.
  PRINT_DEBUG("call to walk_from_vertex_to_point");
  LM_CLOCK_DEBUG(clock_t w_time_start = clock());
  // !!! the operation !!!
  walk( nearest_vertex , p, e, lt);
  LM_CLOCK_DEBUG(clock_t w_time_end = clock() );
  LM_CLOCK_DEBUG(clock_for_walk += (double) (w_time_end - w_time_start) );

  LM_CLOCK_DEBUG(clock_t nw_time_end = clock() );
  LM_CLOCK_DEBUG(clock_nn_and_walk += (double) (nw_time_end - nw_time_start) );

  PRINT_DEBUG( "return from walk_from_vertex_to_point. !!!" );
  PRINT_DEBUG("lt = "<<lt <<", e = "<<e->source()->point() << "->"<< e->target()->point()  );
  PRINT_DEBUG(std::endl  << std::endl);
  //std::cout << "lt = "<<lt <<", e = "<<e->source()->point() << "->"<< e->target()->point() << std::endl << std::endl;
  //getchar

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
  CGAL_precondition_msg(false, 
    "Vertical ray shoot NOT suppored in CDT point location");
  CGAL_assertion (false);
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
  PRINT_DEBUG("in insert_halfedge_to_ln_tree. cv = "<< cv );
  //PRINT_DEBUG(hh->source()->point()<<" towards "<< hh->target()->point());
  
  //create_landmarks_tree();
}

//----------------------------------------------------
/*! go over all vertices, and insert each vertex to the 
landmarks tree.
*/
template <class Planar_map, class Nearest_neighbor>
void Pm_landmarks_point_location<Planar_map, Nearest_neighbor>
::create_landmarks_tree() const
{ 
  PRINT_DEBUG("in create_landmarks_tree");

  if (pm->vertices_begin() == pm->vertices_end()) {
    updated_nn = true;
    PRINT_DEBUG("empty pm. out create_landmarks_tree");    
  }

  //Go over planar map, and create a triangulation of it
  Vertex_iterator   vit;
  NN_Point_list      plist; 

  for (vit=pm->vertices_begin(); vit != pm->vertices_end(); vit++)
  {
    //get point from vertex
    Point_2 p = vit->point();
    double px, py;
    point_to_double_coords(p, px, py);
    Vertex_handle vh = vit;
    NN_Point_2 np (px, py, vh); 
    //insert point into list
    plist.push_back(np); 

    //PRINT_DEBUG("point is= " << p);
  } 

  PRINT_DEBUG("before clean");
  nn.clean();

  PRINT_DEBUG("before init");
  nn.init(plist.begin(), plist.end());

  //the triangulation is now updated
  updated_nn = true;
  PRINT_DEBUG("out create_landmarks_tree" );
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
  PRINT_DEBUG("inside walk_from_vertex_to_point. p= "<< p  << 
                 ", nearest_vertex = "<<nearest_vertex->point() );

  bool new_vertex = false;
  bool found_face = false;
  bool found_vertex_or_edge = false;
  bool found_edge = false;
  bool found_vertex = false;
  Face_handle face;
  Vertex_handle out_vertex;
  Vertex_handle vh = nearest_vertex;
  Halfedge_handle out_edge;

  flipped_edges.clear();  //remove all elements from the flipped_edges list

  int debug_count = 0; 

  do {
    debug_count++; 
    //find the edge out_edge which is the best possibly 
    //pointing to the face containing p

    LM_CLOCK_DEBUG(clock_t ff_time_start = clock());
    new_find_face (p, vh, found_vertex_or_edge, new_vertex, found_face , out_vertex, out_edge, lt);
    //find_face (p, vh, found_vertex_or_edge, new_vertex, found_face , out_vertex, out_edge, lt);
    LM_CLOCK_DEBUG(clock_t ff_time_end = clock() );
    LM_CLOCK_DEBUG(clock_ff += (double) (ff_time_end - ff_time_start) );

    if (found_vertex_or_edge) {
      if (lt == Planar_map::EDGE) {
        PRINT_DEBUG("out_edge = " << out_edge->source()->point() <<"--> "<< out_edge->target()->point() );
        e = out_edge;
      }
      else { // (lt == Planar_map::VERTEX) 
        PRINT_DEBUG("out_vertex = " << out_vertex->point());
        e = out_vertex->incident_halfedges();
      }
      return;
    }
    if (new_vertex) {
      PRINT_DEBUG( "NEW vertex 1 " );
      //check if the new vertex is really closer 
      TR_CLOCK_DEBUG(e_compare_distance ++);
      TR_CLOCK_DEBUG(clock_t compare_distance_ts = clock());
      if (traits->compare_distance(p, out_vertex->point(), vh->point())  
        == SMALLER) {
          vh = out_vertex;
          flipped_edges.clear();  //clear the curves that were flipped
        }
      else {
        PRINT_DEBUG( "Error 2: new vertex is not closer to p than vh! ");
        PRINT_DEBUG( "out_vertex is:  " << out_vertex->point());
        LM_DEBUG(getchar());
        return; 
      }
      TR_CLOCK_DEBUG(clock_t compare_distance_te = clock());
      TR_CLOCK_DEBUG(c_compare_distance += (double) (compare_distance_te - compare_distance_ts));
    }

  } while (new_vertex);

  //get the face that is the best potentially contains p.
  if (found_face) {
    face = out_edge->face();
    PRINT_DEBUG("face that was found is " << out_edge->source()->point() 
                   <<"--> "<< out_edge->target()->point() );
  }
  else {
    std::cerr << "face not found" << std::endl;
    return;
  }

/////////////////////////////IXXXXXXXXXXXXXXXXXXX
//new algorithm should be: (after face was found) 
//1. check if p is in face. 
//  yes: 
//    go over holes. for each hole h:
//      is p in h?
//        yes: face = h. goto 1.
//  no: 
//    call new function: 
//      find_closest_intersection_in_face( gets face, v, p) that  returns e and intersection point (if needed). 
//      (the function go over all edges surrounding p. for each one- as done now - take out of the walk procedure) 
//            face = e->twin. goto1.
//      if intersection not found --- ? error ? 
//////////////////////////////IXXXXXXXXXXXXXXXXXXX

  // IXXXXXXXXXX
  bool p_in_face = false;
  Ccb_halfedge_circulator h_circ;
  
  LM_CLOCK_DEBUG(clock_t na_time_start = clock());//time start
  PRINT_DEBUG("start new algo ");

  do {
    PRINT_DEBUG(std::endl << "inside loop on face ");
    found_vertex = found_edge = p_in_face = false;
    if (face->is_unbounded())  {
      p_in_face = true;
      PRINT_DEBUG("unbounded face ");
    }
    else {    
      h_circ = face->outer_ccb();
      LM_CLOCK_DEBUG(clock_t is_point_ts = clock());
      p_in_face = is_point_in_face(p, h_circ, found_edge, found_vertex, out_edge); 
      LM_CLOCK_DEBUG(clock_t is_point_te = clock() );
      LM_CLOCK_DEBUG(clock_is_point += (double) (is_point_te - is_point_ts));
      PRINT_DEBUG("is_point_in_face returned  "<<  p_in_face );
    }
    if (found_vertex || found_edge ) {
      lt = found_vertex ? Planar_map::VERTEX : Planar_map::EDGE;
      e = out_edge;
      return;
    }
    if (p_in_face){
      //check holes
      PRINT_DEBUG(" p in face. go over holes" );
      Holes_iterator hole_it  = face->holes_begin();
      Holes_iterator hole_end = face->holes_end();  
      bool p_in_hole;
      while (hole_it != hole_end) 
      {
        PRINT_DEBUG(" loop on holes");
        LM_CLOCK_DEBUG(clock_t is_point_ts = clock());
        p_in_hole = is_point_in_face(p, *hole_it, found_edge, found_vertex, out_edge); 
        LM_CLOCK_DEBUG(clock_t is_point_te = clock() );
        LM_CLOCK_DEBUG(clock_is_point += (double) (is_point_te - is_point_ts));
        if (found_vertex || found_edge ) {
          lt = found_vertex ? Planar_map::VERTEX : Planar_map::EDGE;
          e = out_edge;
          return;
        } 
        if (p_in_hole) {
          h_circ = *hole_it; //update the new "face " to be the hole. check its holes.      

          LM_CLOCK_DEBUG(clock_t find_edge_time_start = clock());
          if ( find_edge_to_flip (p, vh, h_circ , out_edge) ) {
            LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
            LM_CLOCK_DEBUG(clock_find_edge += (double) (find_edge_time_end - find_edge_time_start));
            out_edge = out_edge->twin();
            face = out_edge->face();
          }
          else {
            LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
            LM_CLOCK_DEBUG(clock_find_edge += (double) (find_edge_time_end - find_edge_time_start));
            std::cerr << "ERROR 10:  intersection not found" << std::endl;
            LM_DEBUG(getchar());
            e = pm->halfedges_end();
            lt = Planar_map::UNBOUNDED_FACE;
            return;
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
      if ( find_edge_to_flip (p, vh, face->outer_ccb() , out_edge) ) {
        LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
        LM_CLOCK_DEBUG(clock_find_edge += (double) (find_edge_time_end - find_edge_time_start));
        out_edge = out_edge->twin();        
        face = out_edge->face();
        PRINT_DEBUG("after find_edge_to_flip . out_edge changed to twin @ " );
        PRINT_DEBUG("out_edge  =  " << out_edge->source()->point() <<" -->"
                      << out_edge->target()->point() );
      }
      else {
        LM_CLOCK_DEBUG(clock_t find_edge_time_end = clock() );
        LM_CLOCK_DEBUG(clock_find_edge += (double) (find_edge_time_end - find_edge_time_start));
        std::cerr << "ERROR 9:  intersection not found" << std::endl;
        LM_DEBUG(getchar());
        e = pm->halfedges_end();
        lt = Planar_map::UNBOUNDED_FACE;
        return;
      }        

    }
  }while (!p_in_face); 
  //IXXXXXXXXXXX

  LM_CLOCK_DEBUG(clock_t na_time_end = clock() );
  LM_CLOCK_DEBUG(clock_new_alg += (double) (na_time_end - na_time_start) );

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
new_find_face (const Point_2 & p, 
       Vertex_handle vh,
       bool & found_vertex_or_edge, 
       bool & new_vertex, 
       bool & found_face,
       Vertex_handle & out_vertex, 
       Halfedge_handle & out_edge,
       Locate_type& lt  ) const
{ 
  LM_CLOCK_DEBUG( entries_to_find_face++ );
  PRINT_DEBUG("inside find_face. p ="<< p <<" , vh = "<<vh->point() ); 

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
  Point_2 v = vh->point();
  Curve_2 seg(v, p);

  //get halfedges around vh
  Halfedge_around_vertex_circulator circ = vh->incident_halfedges(); 
  Halfedge_around_vertex_circulator circ_done (circ);
  Halfedge_around_vertex_circulator prev = circ;

  //get the curve
  //Curve_2 cv = circ->curve();
  //get the point of cv that is NOT equal to vh.
  //Point_2  cv_other_point = (traits->point_equal(v,cv->source())) ? cv->target() : cv->source() ;  --- no need
  //Point_2  cv_other_point = circ->source()->point();
  //check if cv_other_point is to the left of v,  to the right, or if the curve is vertical
  Comparison_result cv_orient = traits->compare_xy(v,circ->source()->point());

  //check if p is to the left of v,  to the right, or if the segment is vertical
  Comparison_result seg_orient = traits->compare_xy(v,p);

  //save results 
  Comparison_result res1;
  // Comparison_result res2;

  PRINT_DEBUG("seg_orient ="<< seg_orient ); 
  PRINT_DEBUG("cv_orient ="<< cv_orient << ", circ ="<< circ->curve()); 

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
      cv_orient = traits->compare_xy(v,circ->source()->point());
    } while (seg_orient != cv_orient && circ!=circ_done);

    //if exists - go on from next "if" 
    //if not exist - find the curve that is the largest (cw) in this side, and this is the edge we're looking for
    if (seg_orient != cv_orient) 
    {
      //seg is on one side (right/left) and all curves are on the other side (left/right)
      //circ == circ_dome
      do {
        prev = circ;
        circ++;
        TR_CLOCK_DEBUG(e_curves_compare_cw ++);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
        res1 =  traits->curves_compare_y_at_x_cw(circ->curve(), prev->curve(), v, false);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
        TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
        PRINT_DEBUG("circ = " << circ->curve() << "  res1= " << res1 ); 
      } while (res1==LARGER && circ!=circ_done);

      out_edge = prev;
      found_face = true;
      PRINT_DEBUG ( "new_find_face return " << out_edge->curve() );
      return;
    }
  }

  //both curves are to the same side
  if (seg_orient == cv_orient) 
  {
    PRINT_DEBUG("seg_orient == cv_orient : "); 
    TR_CLOCK_DEBUG(e_curves_compare_cw ++);
    TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
    res1 = traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
    TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
    TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
    if (res1 == LARGER) 
    {
      //if the segment is larger than the curve cw, we will go ++ cw with the circ
      //and find a curve that is larger than seg. then we will take the curve that
      //was just before that. 
      PRINT_DEBUG("res1 == LARGER : "); 
      do {
        prev = circ;
        circ++;
        PRINT_DEBUG("circ++ = " << circ->curve() ); 
        cv_orient = traits->compare_xy(v,circ->source()->point());
        if (seg_orient == cv_orient) 
        {
          TR_CLOCK_DEBUG(e_curves_compare_cw ++);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
          res1 =  traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
          TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
          PRINT_DEBUG("circ = " << circ->curve() << "  res1= " << res1 ); 
        }
      } while (res1 == LARGER && seg_orient == cv_orient && circ!=circ_done);

      //if res1 is not larger ==> seg is between prev and circ, return prev
      //if seg_orient != cv_orient, then we changes side and the other side is larger than seg.
      //  then also seg is between prev and circ and we have to return prev.

      if (res1 == LARGER && seg_orient == cv_orient)  //we 're only out the while because the circ end
      {
        //in this case the seg is larger than ALL curves and all curves are to the same side 
        //we need to find the largest of all curves.
        PRINT_DEBUG("circ == circ_done : "); 
        do {
          prev = circ;
          circ++;
          TR_CLOCK_DEBUG(e_curves_compare_cw ++);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
          res1 =  traits->curves_compare_y_at_x_cw(circ->curve(), prev->curve(), v, false);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
          TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
          PRINT_DEBUG("circ = " << circ->curve() << "  res1= " << res1 ); 
        } while (res1 == LARGER && circ!=circ_done);
        //if circ == circ_done, then prev is the largest
        //else if circ is not larger than prev, than prev is the largest
        //anyway, prev is the largest - return it
      }

      out_edge = prev;
      found_face = true;
      PRINT_DEBUG ( "new_find_face return " << out_edge->curve() );
      return;
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
        cv_orient = traits->compare_xy(v,circ->source()->point());
      } while (circ!=circ_done  && seg_orient == cv_orient); 

      if (seg_orient != cv_orient) 
      {
        //loop 2 - until reach the same side again 
        do {
          prev = circ; 
          circ++;
          cv_orient = traits->compare_xy(v,circ->source()->point());
        }  while (seg_orient != cv_orient) ;

        //prev is the last edge from the other side, and curve is the first curve in this side
        TR_CLOCK_DEBUG(e_curves_compare_cw ++);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
        res1 =  traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
        TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));

        //loop 3 - until find curve > seg cw.
        while (res1 == LARGER && seg_orient==cv_orient) 
        {
          prev = circ; 
          circ++;
          cv_orient = traits->compare_xy(v,circ->source()->point());
          if (seg_orient == cv_orient) 
          {
            TR_CLOCK_DEBUG(e_curves_compare_cw ++);
            TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
                        res1 =  traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
            TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
            TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
            PRINT_DEBUG("circ = " << circ->curve() << "  res1= " << res1 ); 
          }
        }
                
        //now we can say that the output edge is prev
        out_edge = prev;
        found_face = true;
        PRINT_DEBUG ( "new_find_face return " << out_edge->curve() );
        return;
      }
  
      // else - (circ == circ_done)
      //there are no curves on the other side, find the smallest (cw) on this side
      do {
          prev = circ;
          circ++;
          TR_CLOCK_DEBUG(e_curves_compare_cw ++);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
          res1 =  traits->curves_compare_y_at_x_cw(circ->curve(), prev->curve(), v, false);
          TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
          TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
          PRINT_DEBUG("circ = " << circ->curve() << "  res1= " << res1 ); 
      } while (res1 == LARGER && circ!=circ_done);
      //now circ < prev ==> circ is the smallest

      //if seg > smallest, smallest++, otherwise, out_edge  = smallest->twin();
      TR_CLOCK_DEBUG(e_curves_compare_cw ++);
      TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
      res1 =  traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
      TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
      TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
      if (res1 == SMALLER) 
      {
        out_edge = circ->twin();
        found_face = true;
        PRINT_DEBUG ( "new_find_face return " << out_edge->curve() );
        return;
      }

      //else: seg > smallest
      circ_done = circ; 
      do  
      {
        prev = circ;
        circ++;
        TR_CLOCK_DEBUG(e_curves_compare_cw ++);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_ts = clock());
        res1 =  traits->curves_compare_y_at_x_cw(seg, circ->curve(), v);
        TR_CLOCK_DEBUG(clock_t curves_compare_cw_te = clock());
        TR_CLOCK_DEBUG(c_curves_compare_cw += (double) (curves_compare_cw_te - curves_compare_cw_ts));
      } while (res1 == LARGER && circ!= circ_done);

      out_edge = prev;
      found_face = true;
      PRINT_DEBUG ( "new_find_face return " << out_edge->curve() );
      return;
    }
    else //EQUAL
    {
      //TODO: specail case - new vertex or on edge ot something
      PRINT_DEBUG ( "specail case: seg is equal cw to circ " << circ->curve() );
      LM_DEBUG(getchar());
      if (traits->point_equal(p,circ->source()->point())) 
      {
        PRINT_DEBUG ( "p is on a vertex ");
        out_edge = circ;
        out_vertex = circ->source();
        lt = Planar_map::VERTEX;
        found_vertex_or_edge = true;
        return; 
      }

      if (traits->point_in_x_range(circ->curve(),p) && 
        traits->curve_compare_y_at_x(p,circ->curve()) == EQUAL) 
      {
        // p lies on cv1  
        PRINT_DEBUG ( "p is on an edge ");
        out_edge = circ;
        lt = Planar_map::EDGE;
        found_vertex_or_edge = true;
        return; 
      }
  
      //p does not lie on cv1 ==> 
      // the target of the equal curve is a better vertex to p 
      std::cerr << " WARNING 11: found closer vertex during new_find_face" << std::endl;
      // out_vertex =  the closer vertex
      out_vertex = circ->source();
      new_vertex = true;
      PRINT_DEBUG( "The new vertex is: "<< out_vertex->point() );
      // check validity (the new vertex is vetween them on a line) @@@@
      return;    
    }
  }

  std::cerr << "ERROR 13: new_find_face did not find the face !" <<std::endl;
  LM_DEBUG(getchar());
}    


  /////////////////////////////////// 6.12
/*
//save special cases
  bool is_equal_curve = false;
  Halfedge_handle  ident_halfedge;



  switch (seg_orient) 
  {
  case EQUAL: //vertical segment
    {
      //TODO - IV
    }
    break;
  case SMALLER: //p is to the right of v:    v -- p
    {
      if (cv_orient == LARGER)  //other_point is to the left of v. cv =  o -- v
      {
        //DID    -  I
        //find a curve to the right of v. 
        while ( (cv_orient == LARGER) && (circ != circ_done) ) 
        {
          Comparison_result cv_orient = traits->compare_x(v,circ->source()->point());
          circ++; 
        } 
                //if exist edge to the right of v, then circ = this edge, cv = circ->curve(), cv_orient = SMALLER
        if (cv_orient != LARGER) 
        {
          cv = circ->curve(); 
        }
        //else, if such edge does not exist, find the upper edge to the left and return it.
        else 
        {
          do {
            //Halfedge_handle prev = circ; 
            //circ++;
            //res1 = traits->curves_compare_y_at_x_left(prev->curve(), circ->curve(), v);
            //if (res1 != SMALLER) {
            //  std::cout << "all from left. found upper left " << prev->curve() << std::endl;
            //  out_edge = prev;
            //  found_face = true;
            //  return;
            }
          } while (circ != circ_done);
        }
      }

      if  (cv_orient == EQUAL)  //cv is vertical
      {
        //DID - II
        //need  to check if vertical up or down
        res1 = traits->compare_xy(v,circ->source()->point());
        if (res1 == SMALLER) 
        {   //if up - test if the next curve (++) is right 
          ++circ;
          cv_orient = traits->compare_x(v, circ->source()->point());
            //if so -> exist edge to the right of v, then circ = this edge, cv = circ->curve(), cv_orient = SMALLER
          if (cv_orient == SMALLER) 
          {
            cv = circ->curve(); 
          }
          //else, if such edge does not exist, return the vertical edge
          else 
          {
            std::cout << "return vertical edge up " << prev->curve() << std::endl;
            out_edge = --circ;
            found_face = true;
            return;
          }
        }
        else 
        {  //if down - check if the previous curve (--) is right
          Halfedge_handle next = circ; 
          --circ;
          cv_orient = traits->compare_x(v, circ->source()->point());
            //if so -> exist edge to the right of v, then circ = this edge, cv = circ->curve(), cv_orient = SMALLER
          if (cv_orient == SMALLER) 
          {
            cv = circ->curve(); 
          }
          //else, if such edge does not exist, return the twin of the  vertical edge down
          else 
          {
            std::cout << "return vertical tein of edge down " << next->curve() << std::endl;
            out_edge = next->twin();
            found_face = true;
            return;
          }    
        }
      }

      if (cv_orient == SMALLER)  //other_point is to the right of v. cv = v -- o
      {
        //check if seg is above / below cv
                res1 = traits->curves_compare_y_at_x_right(seg, cv, v);
        if (res1 == SMALLER)  //seg is after (cw) cv. 
        { 
          bool found_next = false;
          Comparison_result next_orient; 
          Halfedge_handle next = circ; ++next; 
          do { //++ loop until find the (next) curve
            //check orientation of next curve
            next_orient = traits->compare_x(v, next->source()->point());
            if (next_orient != SMALLER) 
            {  //next curve is not to the right - and seg is to the right and after circ =>
              //this is the next we're looking for
              found_next = true;
            }
            else 
            { //need to check if this is the next or we have to continue the loop
              res2 = traits->curves_compare_y_at_x_right(seg, next->curve(), v);
              if (res2 == LARGER) {  //great !
                found_next = true;
              }
              else if (res2 == SMALLER) { //keep going on a loop
                ++next;
              }
              else  //if (res2 == EQUAL) 
              {
                is_equal_curve = true;
                ident_halfedge = next;
              }
            }
          } while (!found_next && !is_equal_curve && next!=circ_done);
          if (found_next) {
            std::cout << "found next " << next->curve() << std::endl;
            out_edge = --next;
            found_face = true;
            return;
          }
          if (next == circ_done) 
          { //all curves are from the right of v, and seg is smaller than all of them. 
            //we need to find the lower curve to the right - this is the out edge
            circ = circ_done; 
            next = circ; ++next;
            Halfedge_handle prev = circ; --prev; 
            do {
              res2 = traits->curves_compare_y_at_x_right(circ->curve(), next->curve(), v);
              if (res2 != SMALLER) 
              { //we found the lower one, its circ. (this should never return EQUAL)
                out_edge = circ;
                found_face = true;
                return;
              }
              prev = circ;
              ++circ;
            } while (circ != circ_done) ; 
            std::cerr << "did not found lower edge to the right" <<std::endl;
          }
        }

        else if (res1 == LARGER)  //seg is before (cw) cv
        { 
          Halfedge_handle prev = circ; --prev; 
          bool found_prev = false;
          Comparison_result prev_orient; 
          do { //--loop until find the (prev) curve
            //check orientation of prev curve
            prev_orient = traits->compare_x(v, prev->source()->point());
            if (prev_orient != SMALLER) 
            {  //prev curve is not to the right - and seg is to the right and before circ =>
              //this is the prev we're looking for
              found_prev = true;
            }
            else 
            { //need to check if this is the prev or we have to continue the loop
              res2 = traits->curves_compare_y_at_x_right(seg, prev->curve(), v);
              if (res2 == SMALLER) {  //great !
                found_prev = true;
              }
              else if (res2 == LARGER) { //keep going on a loop
                --prev;
              }
              else  //if (res2 == EQUAL) 
              {
                is_equal_curve = true;
                ident_halfedge = prev;
              }
            }
          } while (!found_prev && !is_equal_curve && prev!=circ_done);
          if (found_prev) {
            std::cout << "found prev " << prev->curve() << std::endl;
            out_edge = prev;
            found_face = true;
            return;
          }
          if (prev == circ_done) 
          { //all curves are from the right of v, and seg is larger than all of them. 
            //we need to find the upper curve to the right, its twin is the out edge 
            //its twin is better than the lowest curve, because it is closer to p
            circ = circ_done; 
            prev = circ; --prev; 
            do {
              res2 = traits->curves_compare_y_at_x_right(circ->curve(), prev->curve(), v);
              if (res2 != LARGER) 
              { //we found the upper one, its prev. (this should never return EQUAL)
                out_edge = prev->twin();
                found_face = true;
                return;
              }
              circ = prev;
              --prev;
            } while (circ != circ_done) ; 
            std::cerr << "did not found upper edge to the right" <<std::endl;
          }
        }

        else if (res1 == EQUAL) { //equals curve
          is_equal_curve = true;
          ident_halfedge = circ;
        }
      }

    }
    break;

  case LARGER:  //p is to the left of v:    p -- v
    {
      //TODO - III
    }
    break;
  }
*/  
  



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
  LM_CLOCK_DEBUG( entries_to_find_face++ );
  PRINT_DEBUG("inside find_face. p ="<< p <<" , vh = "<<vh->point() ); 

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
    cv2 = circ->curve();
    cv1 = prev->curve();
    //check if seg is between prev_cv and curr_cv
    TR_CLOCK_DEBUG(e_curve_is_between_cw ++);
    TR_CLOCK_DEBUG(clock_t curve_is_between_cw_ts = clock());
    is_btw = traits->curve_is_between_cw(seg, cv1, cv2, vh->point(),   cv_equal_cv1, cv_equal_cv2);
    TR_CLOCK_DEBUG(clock_t curve_is_between_cw_te = clock());
    TR_CLOCK_DEBUG(c_curve_is_between_cw += (double) (curve_is_between_cw_te - curve_is_between_cw_ts));

    if (cv_equal_cv1 || cv_equal_cv2) {
      if (cv_equal_cv1) {
        if (traits->point_equal(p,prev->source()->point())) {
          out_edge = prev;
          out_vertex = prev->source();
          lt = Planar_map::VERTEX;
          found_vertex_or_edge = true;
          return; 
        }

        TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
        TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());
        if (traits->point_in_x_range(cv1,p) && 
          traits->curve_compare_y_at_x(p,cv1) == EQUAL) {
            TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
            TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

            // p lies on cv1    
            out_edge = prev;
            lt = Planar_map::EDGE;
            found_vertex_or_edge = true;
            return; 
          }

          TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
          TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));
    
          //p does not lie on cv1 ==> 
          // the target of the equal curve is a better vertex to p 
          std::cerr << " WARNING 1: found closer vertex during find_face" << std::endl;
          // out_vertex =  the closer vertex
          out_vertex = prev->source();
          new_vertex = true;
          PRINT_DEBUG( "The new vertex is: "<< out_vertex->point() );
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

        TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
        TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());

        if (traits->point_in_x_range(cv2,p) && 
          traits->curve_compare_y_at_x(p,cv2) == EQUAL) {
            TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
            TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));
  
            // p lies on cv1    
            out_edge = circ;
            lt = Planar_map::EDGE;
            found_vertex_or_edge = true;
            return; 
          }

          TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
          TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

          //p does not lie on cv1 ==> 
          // the target of the equal curve is a better vertex to p 
          std::cerr <<" WARNING 2: found closer vertex " << std::endl;
          // out_vertex =  the closer vertex
          out_vertex = prev->source();
          new_vertex = true;
          PRINT_DEBUG( "The new vertex is: "<< out_vertex->point() );
          // check validity (the new vertex is vetween them on a line) @@@@
          return;
      }

    }

    if (is_btw) {
      PRINT_DEBUG(" cv is between " << cv1 << " and " << cv2 ); 
      out_edge = prev;
      found_face = true;
      return;
    }

    prev = circ;
    ++circ;
  }

  //if p not found between edges so far, try prev->first
  if ( !is_btw && !cv_equal_cv1 && !cv_equal_cv2) { 
    cv1 = prev->curve();
    cv2 = first->curve();

    TR_CLOCK_DEBUG(e_curve_is_between_cw ++);
    TR_CLOCK_DEBUG(clock_t curve_is_between_cw_ts = clock());  
    is_btw = traits->curve_is_between_cw(seg, cv1, cv2 , vh->point(), cv_equal_cv1, cv_equal_cv2);
    TR_CLOCK_DEBUG(clock_t curve_is_between_cw_te = clock());
    TR_CLOCK_DEBUG(c_curve_is_between_cw +=  (double) (curve_is_between_cw_te - curve_is_between_cw_ts));
    
    if (cv_equal_cv1 || cv_equal_cv2) {

      if (cv_equal_cv1) {
        if (traits->point_equal(p,prev->source()->point())) {
          out_edge = prev;
          out_vertex = prev->source();
          lt = Planar_map::VERTEX;
          found_vertex_or_edge = true;
          return; 
        }

        TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
        TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());

        if (traits->point_in_x_range(cv1,p) && 
          traits->curve_compare_y_at_x(p,cv1) == EQUAL) {

            TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
            TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

            // p lies on cv1    
            out_edge = prev;
            lt = Planar_map::EDGE;
            found_vertex_or_edge = true;
            return; 
          }

          TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
          TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

          //p does not lie on cv1 ==> 
          // the target of the equal curve is a better vertex to p 
          std::cerr <<" WARNING 3: found closer vertex " << std::endl;
          // out_vertex =  the closer vertex
          out_vertex = prev->source();
          new_vertex = true;
          PRINT_DEBUG( "The new vertex is: "<< out_vertex->point() );
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

        TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
        TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());
        if (traits->point_in_x_range(cv2,p) && 
          traits->curve_compare_y_at_x(p,cv2) == EQUAL) {
            TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
            TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

            // p lies on cv1    
            out_edge = first;
            lt = Planar_map::EDGE;
            found_vertex_or_edge = true;
            return; 
          }
          TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
          TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

          //p does not lie on cv1 ==> 
          // the target of the equal curve is a better vertex to p 
          std::cerr <<" WARNING 4: found closer vertex " << std::endl;
          // out_vertex =  the closer vertex
          out_vertex = prev->source();
          new_vertex = true;
          PRINT_DEBUG( "The new vertex is: "<< out_vertex->point() );
          // check validity (the new vertex is vetween them on a line) @@@@
          return;
      }

    }

    if (is_btw) {
      PRINT_DEBUG(" cv is between " << cv1 << " and " << cv2 ); 
      out_edge = prev;
      found_face = true;
      return;
    }

  }

  std::cerr << "ERROR 4: edge above cw not found !" <<std::endl;
  LM_DEBUG(getchar());
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
*
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
  LM_CLOCK_DEBUG(entries_to_fi++);
  PRINT_DEBUG(" inside find_intersection. p = " << p << ", seg = "<< seg   <<", e = "<< e->curve() ); 

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
  PRINT_DEBUG(" MAX_INTERSECTIONS_BETWEEN_2_CURVES is not defined");
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
    Point_2 vp = seg.source();
    LM_CLOCK_DEBUG(clock_t ni_time_start = clock());    //time start  

    Comparison_result comp_xy_res = traits->compare_xy(p, vp);
    switch (comp_xy_res) {
    case LARGER: //p is on the right of vp (or vertical and up)
      {
        //find nearest intersection to the right of vp
        PRINT_DEBUG(" check nearest intersection to the right of : "<< vp << "seg = "<<seg<< " cv = "<< cv );

        TR_CLOCK_DEBUG(e_nearest_intersection_to_right ++);
        TR_CLOCK_DEBUG(clock_t nearest_intersection_to_right_ts = clock());

        res_obj = traits->nearest_intersection_to_right(seg, cv, vp);

        TR_CLOCK_DEBUG(clock_t nearest_intersection_to_right_te = clock());
        TR_CLOCK_DEBUG(c_nearest_intersection_to_right += (double) (nearest_intersection_to_right_te - nearest_intersection_to_right_ts) );

        PRINT_DEBUG(" out nearest intersection to the right " );
      }
      break;
    case SMALLER: //p is on the left of vp (or vertical and down)
      //find nearest intersection to the right of p
      //res_obj = traits->nearest_intersection_to_right(seg, cv, p); 
      {
        PRINT_DEBUG( "check nearest intersection to the left of : "<< vp <<" seg = "<<seg<< " cv = "<< cv );

        TR_CLOCK_DEBUG(e_nearest_intersection_to_left ++);
        TR_CLOCK_DEBUG(clock_t nearest_intersection_to_left_ts = clock());

        res_obj = traits->nearest_intersection_to_left(seg, cv, vp); 

        TR_CLOCK_DEBUG(clock_t nearest_intersection_to_left_te = clock());
        TR_CLOCK_DEBUG(c_nearest_intersection_to_left += (double) (nearest_intersection_to_left_te - nearest_intersection_to_left_ts));

        PRINT_DEBUG(" out nearest intersection to the left " );
      }
      break;
    default: //should not be equal
      CGAL_assertion (false);
    }

    LM_CLOCK_DEBUG(    
      clock_t ni_time_end = clock();//time end
      double ni_period = (double) (ni_time_end - ni_time_start);
      clock_ni += ni_period;
    );  

    // Empty object is returned - no intersection.
    if (res_obj.is_empty()) {
      //if (comp_xy_res == SMALLER) {
      //  if (traits->is_on_segment (cv, p)) {
      //    closest_interect_point = p;
      //    found_edge = true;
      //    //p is on edge;
      //    return;
      //  }
      //}
      PRINT_DEBUG( " no intersection"); 
      return;
    }

    // Intersection is a point
    else if (assign(inter_point, res_obj)) {
      PRINT_DEBUG(" intersection is a point = "<< inter_point ); 

      //if the intersection point is p, then p is on cv - found
      if (traits->point_equal(p, inter_point)) {
        closest_interect_point = inter_point;
        found_edge = true;
        return;
      } 
      //if the intersection point is vh->point, ignore it ! 
      if (traits->point_equal(seg.source(), inter_point)) {
        return;
      }

      //check if the intersection point is an end point of e
      if (traits->point_equal(e->source()->point(), inter_point)) {
        new_vertex = true;
        out_vertex = e->source();
        return;
      }
      if (traits->point_equal(e->target()->point(), inter_point)) {
        new_vertex = true;
        out_vertex = e->target();
        return;
      }
      //else - regular intersection

      if ((!define_max) || (max_intersections > 1)) {
        std::cerr <<" WARNING 5: need to check tangent point " << std::endl;
        // @@@@ check if this is not a tangent point !
        // @@@@ to check tangent point need to do curve_compare_y_at_x_left and right (
        // @@@@ if tangent point - don't change side.
      }

      //if not first intersection - compare inter_point with closest_interect_point (ref p)
      //if inter_point closer - change it .
      if ((num_of_intersections == 0) || 
        (traits->compare_distance(p, inter_point, closest_interect_point) 
        == SMALLER)) {
          closest_interect_point = inter_point;

          // split seg and save the part of seg closer to p
          TR_CLOCK_DEBUG(e_curve_split ++);
          TR_CLOCK_DEBUG(clock_t curve_split_ts = clock());
      
          traits->curve_split(seg, seg1, seg2, inter_point);
          
          TR_CLOCK_DEBUG(clock_t curve_split_te = clock());
          TR_CLOCK_DEBUG(c_curve_split += (double) (curve_split_te - curve_split_ts));
  
          Point_2 seg1_src = seg1.source();
          Point_2 seg2_trg = seg2.target();
          if ((traits->compare_xy(seg1_src, seg2_trg)) == comp_xy_res) {
            seg = seg1; 
          }
          else { 
            seg = seg2;
          }
        }
        change_side = change_side ? false : true; 
        num_of_intersections++;  
    }

    // Intersection is a segment
    else if (assign(overlap_seg, res_obj))
    {
      PRINT_DEBUG("intersection is an overlapped segment"); 
      Point_2 ov_src = overlap_seg.source();
      Point_2 ov_trg = overlap_seg.target();
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

      std::cerr <<" WARNING 6: this case was not implemented !!!" << std::endl;
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
*/

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
  LM_CLOCK_DEBUG( entries_is_point_in_face ++) ;
  PRINT_DEBUG("inside is_point_in_face. face = " << face->source()->point()
                <<"-->" << face->target()->point());

  found_edge = false;
  found_vertex = false;
  int number_of_edges_above_p = 0;

  //if  (face->is_unbounded())
  //  return (true);

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
      TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());

      Comparison_result compare_y_at_x_res =  traits->curve_compare_y_at_x(p, cv);

      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x += (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

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
*
template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_closest_intersection_in_face (const Point_2 & p,                  //input seg src 
                     Vertex_handle  v,           //input vertex
                  const Ccb_halfedge_circulator & face,                                  //input face
                  Halfedge_handle  & out_edge) const              //output if on curve        
{
  PRINT_DEBUG("inside find_closest_intersection_in_face.   " );

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
    PRINT_DEBUG("curr = " << p1 << "-->" << p2 );

    //@@@@ check if these operations take too long
    TR_CLOCK_DEBUG(e_point_in_x_range += 4);
    TR_CLOCK_DEBUG(clock_t point_in_x_range_ts = clock());  
    p_in_x_range = traits->point_in_x_range(cv, p);
    v_in_x_range = traits->point_in_x_range(cv, vp);
    p1_in_x_range = traits->point_in_x_range(seg, p1);
    p2_in_x_range = traits->point_in_x_range(seg, p2);
    TR_CLOCK_DEBUG(clock_t point_in_x_range_te = clock());
    TR_CLOCK_DEBUG(c_point_in_x_range +=  (double) (point_in_x_range_te - point_in_x_range_ts));

    //
    PRINT_DEBUG("p_in_x_range = " << p_in_x_range << " , v_in_x_range = " << v_in_x_range);
    PRINT_DEBUG("p1_in_x_range = " << p1_in_x_range << " , p2_in_x_range = " << p2_in_x_range);

    check_real_intersection = false;
    if ( p_in_x_range ||  v_in_x_range || p1_in_x_range || p2_in_x_range)
    {    
      bool intersect;
      if (check_approximate_intersection(seg, cv, intersect)) 
      {
        if (intersect) {
          if (! is_found_intersection) { //this is the first intersection
            is_found_intersection = true;
            out_edge = curr;
          }
          else { //this is not the first inetrsection - need to check who is closer (real intersection)
            check_real_intersection = true;
          }
        }
      }
      else {
        check_real_intersection = true;
      }
    }
    
    if (check_real_intersection) 
    {
      Point_2 temp_inter_point;
      //check intersection and update temp_inter_point
      if ( find_real_intersection(p, v, curr, temp_inter_point) ) 
      {
        PRINT_DEBUG("find_real_intersection returned true ");
        if (! is_found_intersection) { //this is the first intersection
          is_found_intersection = true;
          is_inter_point_updated = true;
          inter_point = temp_inter_point;
          out_edge = curr;
        }
        else { //compare to previous intersection
          if (!is_inter_point_updated)  {
            //check intersection between out_edge and v-p and update inter_point
            if ( ! find_real_intersection(p, v, out_edge, inter_point) ) {
              std::cerr << "ERROR 6: not found intersection when exists" << std::endl;
              LM_DEBUG(getchar());
            }
          }

          //compare which point is closer. 
          TR_CLOCK_DEBUG(e_compare_distance ++);
          TR_CLOCK_DEBUG(clock_t compare_distance_ts = clock());
          if (traits->compare_distance(p, temp_inter_point, inter_point) 
            == SMALLER) 
          { // the temp_inter point is closer to p
            inter_point = temp_inter_point;
            out_edge = curr;
          }
          is_inter_point_updated = true;
          TR_CLOCK_DEBUG(clock_t compare_distance_te = clock());
          TR_CLOCK_DEBUG(c_compare_distance += (double) (compare_distance_te - compare_distance_ts));
        }
      }
      else {
        PRINT_DEBUG("find_real_intersection returned false " );
      }
    }

    ++curr;
  } while (curr != last);  

  return (is_found_intersection);
}*/

//----------------------------------------------------
/*!
* find the edge in the face that is intersecting with the segment v - p.  (out_edge)
* each edge can be chosen only once, and will get into the switch_curves list 
* returns false if there is no intersection, thus no edge to flip. 
* 
* \param p - the input point.
* \param v - the input vertex
* \param face - the input face
* \param out_edge - the output edge
*/
template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
find_edge_to_flip (const Point_2 & p,                  //input seg src 
                     Vertex_handle  v,           //input vertex
                  const Ccb_halfedge_circulator & face,                                  //input face
                  Halfedge_handle  & out_edge) const              //output edge to flip      
{
  //we want to eliminate the calls to nearest. 
  //this means we will not call to find_real_intersection at all. 
  //this also means that we will take the first edge that is intersecting, 
  //and not check which is closer. 
  //what we will check is whether this edge was selected already (and flipped), 
  // in this case we will not flip it again, but move to the next edge
  LM_CLOCK_DEBUG( entries_find_edge++ );
  PRINT_DEBUG("inside find_edge_to_flip.   " );

  Point_2 vp = v->point();
  Curve_2 seg(vp, p);   //create a segment vh--p. 

  //loop on all edges in this face
  Ccb_halfedge_circulator curr = face;
  Ccb_halfedge_circulator last =  curr;
  bool p_in_x_range, v_in_x_range, p1_in_x_range, p2_in_x_range;  

  do  {
    Curve_2 cv = curr->curve();
    Point_2 p1 = curr->source()->point();
    Point_2 p2 = curr->target()->point();

    // check if the curve was already flipped - in this case, don't check it at all.
    Std_edge_iterator found1 = std::find (flipped_edges.begin(), flipped_edges.end(), curr);
    Std_edge_iterator found2 = std::find (flipped_edges.begin(), flipped_edges.end(), curr->twin());
    if (found1 != flipped_edges.end() || found2 != flipped_edges.end()) {
      PRINT_DEBUG("curve " << p1 << "-->" << p2  <<" was found in the list");      
    }
    else {
      PRINT_DEBUG("curr = " << p1 << "-->" << p2 );

      //check x-range
      TR_CLOCK_DEBUG(e_point_in_x_range += 4);
      TR_CLOCK_DEBUG(clock_t point_in_x_range_ts = clock());  
      p_in_x_range = traits->point_in_x_range(cv, p);
      v_in_x_range = traits->point_in_x_range(cv, vp);
      p1_in_x_range = traits->point_in_x_range(seg, p1);
      p2_in_x_range = traits->point_in_x_range(seg, p2);
      TR_CLOCK_DEBUG(clock_t point_in_x_range_te = clock());
      TR_CLOCK_DEBUG(c_point_in_x_range +=  (double) (point_in_x_range_te - point_in_x_range_ts));
      PRINT_DEBUG("p_in_x_range = " << p_in_x_range << " , v_in_x_range = " << v_in_x_range);
      PRINT_DEBUG("p1_in_x_range = " << p1_in_x_range << " , p2_in_x_range = " << p2_in_x_range);

      if ( p_in_x_range ||  v_in_x_range || p1_in_x_range || p2_in_x_range)
      {    
        bool intersect;
        
        LM_CLOCK_DEBUG(clock_t check_app_ts = clock());  
        bool check_res = check_approximate_intersection(seg, cv, intersect);
        LM_CLOCK_DEBUG(clock_t check_app_te = clock());
        LM_CLOCK_DEBUG(clock_check_app +=  (double) (check_app_te - check_app_ts));

        if (check_res) 
        {
          if (intersect) {
            out_edge = curr;
            flipped_edges.push_back(curr);
            return (true);
          }
        }
        else {
          std::cerr << "ERROR 12: check_approximate_intersection did not return an answer." << std::endl;
        }
      }
      //else - there is not intersection.
    }

    ++curr;
  } while (curr != last);  

  return (false);
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
*
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

  LM_CLOCK_DEBUG(clock_t fi_time_start = clock());  //time start

  find_intersection (p, seg, e, num_intersections, change_side, found_edge, 
    out_point, new_vertex, temp_vertex);

  LM_CLOCK_DEBUG(  
    clock_t fi_time_end = clock();//time end
    double fi_period = (double) (fi_time_end - fi_time_start);
    clock_fi += fi_period;
  );  

  //check results
  if (new_vertex) { 
    std::cerr << "ERROR 7: new vertex in find_real_intersection " << std::endl;
    LM_DEBUG(getchar());
    //ixx: @@@@ maybe later we will also have to deal with this case
    return (false);
  }
  if (found_edge) {
    std::cerr << "ERROR 8: new edge in find_real_intersection " << std::endl;
    LM_DEBUG(getchar());
    //ixx: @@@@ maybe later we will also have to deal with this case
    return (false);
  }

  return (change_side) ;
}

*/
//----------------------------------------------------
/*!
* Trim the segment to the x range of the input curve
* 
* \param seg - the input seg.
* \param cv - the input curve
* \param trimmed_seg - the output trimmed segment
*/
template <class Planar_map, class Nearest_neighbor>
bool Pm_landmarks_point_location<Planar_map, Nearest_neighbor>::
check_approximate_intersection (const Curve_2 & seg,  
                                const Curve_2 & cv,
                                bool & intersect) const 
{
  LM_CLOCK_DEBUG( entries_to_check_app ++ );

  Point_2 seg_right = traits->curve_righttop_most(seg);
  Point_2 seg_left = traits->curve_leftlow_most(seg);
  Point_2 cv_right = traits->curve_righttop_most(cv);
  Point_2 cv_left = traits->curve_leftlow_most(cv);
  intersect = false;

  PRINT_DEBUG("seg_right =  " << seg_right << " , seg_left = " << seg_left);
  PRINT_DEBUG("cv_right = " << cv_right << " , cv_left = " << cv_left);
  
  //compare the 2 left end-points and the 2 right end-points
  TR_CLOCK_DEBUG(e_compare_xy += 2);
  TR_CLOCK_DEBUG(clock_t compare_xy_ts = clock());
  Comparison_result comp_left_xy_res = traits->compare_xy(seg_left, cv_left);
  Comparison_result comp_right_xy_res = traits->compare_xy(seg_right, cv_right);
  TR_CLOCK_DEBUG(clock_t compare_xy_te = clock());
  TR_CLOCK_DEBUG(c_compare_xy += (double) (compare_xy_te - compare_xy_ts));


  // the left end-point of the segments is equal to the left end-point of the curve
  if (comp_left_xy_res == EQUAL) 
  {
    PRINT_DEBUG("the left end-point of the segments is equal to the left end-point of the curve");
    // compare to the right of the curves
    TR_CLOCK_DEBUG(e_curves_compare_y_at_x_right ++);
    TR_CLOCK_DEBUG(clock_t curves_compare_y_at_x_right_ts = clock());
    Comparison_result curves_comp = traits->curves_compare_y_at_x_right(seg, cv, seg_left); 
    TR_CLOCK_DEBUG(clock_t curves_compare_y_at_x_right_te = clock());
    TR_CLOCK_DEBUG(c_curves_compare_y_at_x_right += 
      (double) (curves_compare_y_at_x_right_te - curves_compare_y_at_x_right_ts));
  
    if (curves_comp == EQUAL) {
      PRINT_DEBUG("overlap !!!");
      return (false); 
    }                                                                         
    if (comp_right_xy_res == SMALLER) { //the segments ends before (to the right) of the curve 
      TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(seg_right,cv) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=   
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts));

            if (curve_comp_right_res == EQUAL) {
        PRINT_DEBUG("2 points collide");
        return (false); 
      }
      if (curves_comp != curve_comp_right_res) { //the segment is on the other side of the segment's end-point
        PRINT_DEBUG("intersecting");
        intersect = true;
        return (true);
      }
      else {
        //intersect = false;
        return (true);
      }
    }
    else if (comp_right_xy_res == LARGER) { //the curve ends before (to the right) of the segment
      TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(cv_right,seg) ;      
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

      if (curve_comp_right_res == EQUAL) {
        PRINT_DEBUG("2 points collide");
        return (false); 
      }
      if (curves_comp == curve_comp_right_res) {//the segment is on the same side as the curve's end-point
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

  // the right end-point of the segments is equal to the right end-point of the curve
  else if (comp_right_xy_res == EQUAL) 
  {
    PRINT_DEBUG("the right end-point of the segments is equal to the right end-point of the curve");
    // compare to the left of the curves
    TR_CLOCK_DEBUG(e_curves_compare_y_at_x_left ++);
    TR_CLOCK_DEBUG(clock_t curves_compare_y_at_x_left_ts = clock());  
    Comparison_result curves_comp = traits->curves_compare_y_at_x_left(seg, cv, seg_right); 
    TR_CLOCK_DEBUG(clock_t curves_compare_y_at_x_left_te = clock());
    TR_CLOCK_DEBUG(c_curves_compare_y_at_x_left +=  
      (double) (curves_compare_y_at_x_left_te - curves_compare_y_at_x_left_ts)); 

    if (curves_comp == EQUAL) {
      PRINT_DEBUG("overlap !!!");
      return (false); 
    }
    if (comp_left_xy_res == SMALLER) { //the curve ends before (to the left) of the segment 
      TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res =traits->curve_compare_y_at_x(cv_left,seg) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

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
    else if (comp_left_xy_res == LARGER) { //the segment ends before (to the left) of the curve 
      TR_CLOCK_DEBUG(e_curve_compare_y_at_x ++);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res =traits->curve_compare_y_at_x(seg_left,cv) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

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
    {  //the segment is inside the curve's x-range
      PRINT_DEBUG("the segment is inside the curve's x-range.");
      LM_DEBUG ( //checks for debugging - remove later
        if (! traits->point_in_x_range(cv,seg_left) ) {
          std::cerr << "! traits->point_in_x_range(cv,seg_left) " << std::endl; return (false);}
        if (! traits->point_in_x_range(cv,seg_right) )  {
          std::cerr << "(! traits->point_in_x_range(cv,seg_right)  " << std::endl; return (false);}
      )

      TR_CLOCK_DEBUG(e_curve_compare_y_at_x += 2);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res = traits->curve_compare_y_at_x(seg_left,cv);
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(seg_right,cv) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

      if ((curve_comp_left_res == EQUAL) || (curve_comp_right_res == EQUAL) )
      {
        std::cerr <<" WARNING 7: left or right endpoint of the segment is on the curve " << std::endl;
        //this should not happen since p is on the curve, we should have find it already, 
        //and if v is on the curve, than v should have cut the curve in two
        return (false);
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
        if (! traits->point_in_x_range(seg,cv_left) ) {
          std::cerr << "! traits->point_in_x_range(seg,cv_left)  " << std::endl; return (false);}
        if (! traits->point_in_x_range(seg,cv_right) ) {
          std::cerr << "! traits->point_in_x_range(seg,cv_right)  ---" << std::endl; return (false); }
      )

      TR_CLOCK_DEBUG(e_curve_compare_y_at_x += 2);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res = traits->curve_compare_y_at_x(cv_left,seg);
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(cv_right,seg) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

      if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL) 
      {
        std::cerr <<" WARNING 8: left or right endpoint of the curve is on the segment " << std::endl;
        //this means that the endpoint of the curve is a closer vertex to p than v is. 
        return (false);
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
  }
  else { //one endpoint of the first curve is inside the other curve's  x-range and one endpoitn is out, and vise versa.
    if (comp_left_xy_res == SMALLER) 
    { //if the segment right point is in and the curve's left point
      PRINT_DEBUG("the segment right point is in and the curve's left point");
      LM_DEBUG ( //checks for debugging - remove later
        if (! traits->point_in_x_range(seg,cv_left) ) {
          std::cerr << "! traits->point_in_x_range(seg,cv_left)  " << std::endl;  return (false);  }
        if (! traits->point_in_x_range(cv,seg_right) ) {
          std::cerr << "! traits->point_in_x_range(cv,seg_right)  " << std::endl; return (false); }
      )

      TR_CLOCK_DEBUG(e_curve_compare_y_at_x += 2);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res = traits->curve_compare_y_at_x(cv_left,seg);
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(seg_right,cv) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

      if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL) 
      {
        if (traits->point_equal(seg_right, cv_left)) {
          //this case is o.k., just return that there is no intersection, because probably v is the same for both curves
          PRINT_DEBUG("segment's right endpoint = curve's left endpoint ");
          PRINT_DEBUG("no intersection");
          //intersect = false;
          return (true);    
        }
        std::cerr <<" WARNING 9: left or right endpoint is on the curve " << std::endl;
        //its either the segment's end-point is on the curve, and this should not 
        //happen since if p is on the curve, we should have find it already, and if v 
        //is on the curve, than v should have cut the curve in two.
        //another options is that the curve end-point is on the segment, and in this 
        //case this endpoint is a closer vertex to p than v      
        return (false);
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
      PRINT_DEBUG("the curve's right point is in and the segment's left point");
      LM_DEBUG ( //checks for debugging - remove later
        if (! traits->point_in_x_range(cv,seg_left) ) {
          std::cerr << "! traits->point_in_x_range(cv,seg_left)  " << std::endl; return (false);  }
        if (! traits->point_in_x_range(seg,cv_right) ) {
          std::cerr << "! traits->point_in_x_range(seg,cv_right) *** " << std::endl; return (false); }
      )

      TR_CLOCK_DEBUG(e_curve_compare_y_at_x += 2);
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_ts = clock());  
      Comparison_result curve_comp_left_res = traits->curve_compare_y_at_x(seg_left,cv);
      Comparison_result curve_comp_right_res =traits->curve_compare_y_at_x(cv_right,seg) ;
      TR_CLOCK_DEBUG(clock_t curve_compare_y_at_x_te = clock());
      TR_CLOCK_DEBUG(c_curve_compare_y_at_x +=  
        (double) (curve_compare_y_at_x_te - curve_compare_y_at_x_ts)); 

      if (curve_comp_right_res == EQUAL || curve_comp_left_res == EQUAL)
      {
        if (traits->point_equal(seg_left, cv_right)) {
          //this case is o.k., just return that there is no intersection, because probably v is the same for both curves
          PRINT_DEBUG("segment's left endpoint = curve's right endpoint ");
          PRINT_DEBUG("no intersection");
          //intersect = false;
          return (true);    
        }
        std::cerr <<" WARNING 10: left or right endpoint is on the curve " << std::endl;
        //same explanation as in warning 9
        return (false);
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


  std::cerr << "ERROR 11: not existing option in check_approximate_intersection " << std::endl;
  LM_DEBUG(getchar());
  return (false);
}

/* */

CGAL_END_NAMESPACE

#endif  //CGAL_PM_LANDMARKS_POINT_LOCATION_C

/*
  TR_CLOCK_DEBUG(e_compare_distance ++);
  TR_CLOCK_DEBUG(clock_t compare_distance_ts = clock());

  TR_CLOCK_DEBUG(clock_t compare_distance_te = clock());
  TR_CLOCK_DEBUG(c_compare_distance += (double) (compare_distance_te - compare_distance_ts));
*/



