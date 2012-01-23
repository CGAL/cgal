// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
#ifndef CGAL_ARR_LM_MIDDLE_EDGES_GENERATOR_H
#define CGAL_ARR_LM_MIDDLE_EDGES_GENERATOR_H

/*! \file
* Definition of the Arr_middle_edges_landmarks_generator<Arrangement> template.
*/

#include <CGAL/Arr_point_location/Arr_lm_generator_base.h>

namespace CGAL {

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_lm_generator and  implements the 
* function called "void _create_point_list(Point_list &)" 
* to creates the set of landmarks, which are the middle points of the 
* arrangement edges, which must be segments !
* IMPORTANT: THIS ALGORITHM WORKS ONLY FOR SEGMENTS !!!
*/
template <class Arrangement_, 
          class Nearest_neighbor_ 
          = Arr_landmarks_nearest_neighbor <typename Arrangement_::Traits_2> >
class Arr_middle_edges_landmarks_generator 
  : public Arr_landmarks_generator_base <Arrangement_, Nearest_neighbor_>
{
public:
  typedef Arrangement_					  Arrangement_2;
  typedef Arr_middle_edges_landmarks_generator<Arrangement_2, 
                                               Nearest_neighbor_>    Self;
  typedef Arr_landmarks_generator_base<Arrangement_2, Nearest_neighbor_>  Base;

  typedef typename Arrangement_2::Traits_2		Traits_2;
  typedef typename Arrangement_2::Edge_const_iterator	Edge_const_iterator;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Vertex_handle		Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle	Halfedge_handle;
  typedef typename Arrangement_2::Face_handle		Face_handle;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
                                                       Ccb_halfedge_circulator;
  typedef typename Base::NN_Points_set                 NN_Points_set;
  typedef typename Base::NN_Point_2                    NN_Point_2;
  
  typedef typename Traits_2::Point_2		       Point_2;
  typedef std::vector<Point_2>			       Points_set;
  
private:
  
  /*! Copy constructor - not supported. */
  Arr_middle_edges_landmarks_generator (const Self& );
  
  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );
  
  
public: 
  /*! Constructor. */
  Arr_middle_edges_landmarks_generator 
  (const Arrangement_2& arr, int /* lm_num */ = -1) : 
    Base (arr)
  {
    //CGAL_PRINT_DEBUG("Arr_middle_edges_landmarks_generator constructor.");
    this->build_landmark_set();
  }
  
  //Observer functions that should be empty, because they
  // got nothing to do with middle edges
  //-------------------------------------------------
  virtual void after_create_vertex (Vertex_handle /* v */){}
  virtual void after_split_face (Face_handle /* f */,
                                 Face_handle /* new_f */, bool /* is_hole */)
  {}
  virtual void after_add_hole (Ccb_halfedge_circulator /* h */) {}
  
  virtual void after_merge_face (Face_handle /* f */){}
  virtual void after_move_hole (Ccb_halfedge_circulator /* h */){}
  virtual void after_remove_vertex () {}
  virtual void after_remove_hole (Face_handle /* f */) {}
  
  
protected:
  /*!
   * create a set of middle_edges points 
   * the number of points is equal to the number of edges in the arrangement.
   */
  virtual void _create_nn_points_set (NN_Points_set &nn_points) 
  {
    //CGAL_PRINT_DEBUG("create_middle_edges_points_list");
    Edge_const_iterator    eit;
    Halfedge_const_handle  hh;
    Arrangement_2 *arr = this->arrangement();

    if(arr->number_of_vertices() == 1)
    {
      //special treatment for arrangement with one isolated verrtex
      Vertex_const_iterator vit = arr->vertices_begin();
      CGAL::Object obj = CGAL::make_object(vit);
      Point_2 p (vit->point());
      NN_Point_2 np(p, obj);
      nn_points.push_back(np);

      return;
    }
    for (eit=arr->edges_begin(); eit != arr->edges_end(); eit++)
    {
      //get 2 endpoints of edge
      hh = eit;
      const Point_2& p1 = hh->source()->point();
      const Point_2& p2 = hh->target()->point();
      Point_2 p ((p1.x()+p2.x())/2, (p1.y()+p2.y())/2);
      
      //CGAL_PRINT_DEBUG("mid point is= " << p);
      
      CGAL::Object obj = CGAL::make_object(hh); 
      NN_Point_2 np(p, obj); 
      nn_points.push_back(np);
    } 
  }
  
  virtual void _create_points_set (Points_set & /* points */)
  {
    std::cerr << "should not reach here!"<< std::endl;
    CGAL_error();
  }
};

} //namespace CGAL

#endif
