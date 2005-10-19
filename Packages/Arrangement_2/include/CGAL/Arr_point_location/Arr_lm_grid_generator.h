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
#ifndef CGAL_ARR_LM_GRID_GENERATOR_H
#define CGAL_ARR_LM_GRID_GENERATOR_H

/*! \file
* Definition of the Arr_grid_landmarks_generator<Arrangement> template.
*/

//#include <CGAL/Arr_point_location/Arr_lm_generator.h>
#include <list>
#include <algorithm>   // for random_shuffle
#include <vector>   
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>

CGAL_BEGIN_NAMESPACE

#define CONICS

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_lm_generator and  implements the 
* function called "void _create_point_list(Point_list &)" 
* to creates the set of landmarks on the grid.
* the size of the grid is determined by the number of landmarks. 
*/
template <class Arrangement_>
class Arr_grid_landmarks_generator 
  : public Arr_observer <Arrangement_>
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef Arr_grid_landmarks_generator<Arrangement_2>   Self;
  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;

	typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
	typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
	typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
	typedef typename Arrangement_2::Vertex_handle		      Vertex_handle;
	typedef typename Arrangement_2::Halfedge_handle		    Halfedge_handle;
	typedef typename Arrangement_2::Face_handle			      Face_handle;
	typedef typename Arrangement_2::Vertex_const_iterator 
													Vertex_const_iterator;
    typedef typename Arrangement_2::Ccb_halfedge_circulator 
													Ccb_halfedge_circulator;

  typedef typename Traits_2::Approximate_number_type	  ANT;

#ifdef SEGMENTS
  typedef typename Traits_2::Kernel::FT                 FT;
#elif defined (CONICS)
  typedef CGAL::CORE_algebraic_number_traits        Nt_traits;
  typedef Nt_traits::Algebraic                             FT;
#endif

  typedef typename Traits_2::Point_2                    Point_2;

  typedef std::vector<Point_2>                          Points_set;
	typedef std::pair<Point_2,CGAL::Object>				    PL_pair;
	typedef std::vector<PL_pair>						          Pairs_set;
	typedef typename std::vector<PL_pair>::iterator		Pairs_iterator;

protected:

	typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;

  // Data members:
	const Traits_wrapper_2  *traits;    // Its associated traits object.
	bool  ignore_notifications;	
	bool  updated;

  int	        number_of_landmarks;
  Pairs_set		lm_pairs;

  //bounding box of the arrangement
  ANT x_min, x_max, y_min, y_max;
  FT step_x, step_y;
  int sqrt_n;

private:

  /*! Copy constructor - not supported. */
  Arr_grid_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

  
public: 
    /*! Constructor. */
    Arr_grid_landmarks_generator 
      (const Arrangement_2& arr, int lm_num = -1) : 
        Arr_observer<Arrangement_2> (const_cast<Arrangement_2 &>(arr)), 
		  ignore_notifications (false), 
		  updated (false),   
		  number_of_landmarks (lm_num)
    {
      PRINT_DEBUG("Arr_grid_landmarks_generator constructor. "
        <<"number_of_landmarks = "<< number_of_landmarks); 

		  traits = static_cast<const Traits_wrapper_2*> (arr.get_traits());
      build_landmarks_set();
    }

	  /*!
	  * Creates the landmarks set (choosing the landmarks) , 
	  * and saving them in the nearest neighbor search structure.
	  * This is a pure virtual function (must be implemented in 
	  * the class that derives from this one)
	  */
    virtual void build_landmarks_set ()
    {
      //Go over planar map, and insert all vertices as landmarks
      Points_set    points; 

      _create_points_set(points);

      //locate the landmarks in the arrangement using batched point location
      // global function.
      locate(*(this->arrangement()),points.begin(),points.end(),
             std::back_inserter(lm_pairs));

      Pairs_iterator pit;
      LM_DEBUG(int count =0);
      for (pit = lm_pairs.begin(); pit != lm_pairs.end(); pit++)
      {
        PRINT_DEBUG("grid point "<<count++<<" is= "<< pit->first);
      }

      updated = true;
    }

	  /*!
	  * clear the tree
	  */
	  virtual void clear_landmarks_set ()
	  {
		  PRINT_DEBUG("clear_landmarks_set.");

      //clear the database
      lm_pairs.clear();

		  updated = false;	
	  }

	  /*!
	  * get the nearest neighbor (landmark) to the given point
	  */
	  virtual Point_2 get_closest_landmark (Point_2 p, Object &obj)
	  {
		  CGAL_assertion(updated);
		  PRINT_DEBUG("step_x = "<<step_x << ", step_y = "<<step_y);
      PRINT_DEBUG("x_min = "<<x_min << ", y_min = "<<y_min);
		  PRINT_DEBUG("sqrt_n = "<<sqrt_n);

      //approximate the steps
      Point_2 step_p (step_x, step_y);
      ANT ant_step_x = traits->approximate_2_object()(step_p, 0);
      ANT ant_step_y = traits->approximate_2_object()(step_p, 1);

		  PRINT_DEBUG("ant_step_x = "<<ant_step_x << ", ant_step_y = "<<ant_step_y);

      //calculate the index of the point 
      ANT x = traits->approximate_2_object()(p, 0);
      ANT y = traits->approximate_2_object()(p, 1);

		  PRINT_DEBUG("x = "<<x << ", y = "<<y);

      int i = static_cast<int>(((x-x_min)/ant_step_x)+0.5);
      int j = static_cast<int>(((y-y_min)/ant_step_y)+0.5);

      if (x > x_max || x < x_min || y > y_max || y< y_min)
      {
        PRINT_DEBUG(" query out of range ");
        if (x > x_max)
          i = sqrt_n-1;
        if (y > y_max)
          j = sqrt_n-1;
        if (x < x_min)
          i = 0;
        if (y < y_min)
          j = 0;
      }


		  PRINT_DEBUG("i = "<<i << ", j = "<<j);

      int index = sqrt_n * i + j;

		  PRINT_DEBUG("index = "<<index);

      obj = lm_pairs[index].second;

      return lm_pairs[index].first;
	  }

protected:

   /*!
   * create a set of landmark points on a grid. 
   * the number of points is given as a parametr to this class' constructor.
   * We first calculate the Arrangement's bounding rectangle. 
   * This is actually the bounding rectangle of the Arrangement's vertices.
   * then, we divide the size of each rectangle edge (corresponding to x and y
   * axis) with the number of landmarks, to get the step in x and in y.
   */
  virtual void _create_points_set (Points_set & points)
  {
    PRINT_DEBUG("create_grid_points_list");
    
    //init min/max
    Arrangement_2 *arr = this->arrangement();
    Vertex_const_iterator vit = arr->vertices_begin();
    x_min = x_max = traits->approximate_2_object()(vit->point(), 0);
    y_min = y_max = traits->approximate_2_object()(vit->point(), 1);

    //find bounding box
    ANT x, y;
    Point_2 left, right, top, bottom;
    left = right = top = bottom = vit->point();

    for (vit=arr->vertices_begin(); vit != arr->vertices_end(); vit++)
    {
      x = traits->approximate_2_object()(vit->point(), 0);
      y = traits->approximate_2_object()(vit->point(), 1);
      if (x < x_min) { x_min = x; left = vit->point();}
      if (x > x_max) { x_max = x; right = vit->point();}
      if (y < y_min) { y_min = y; bottom = vit->point();}
      if (y > y_max) { y_max = y; top = vit->point();}
    }

    PRINT_DEBUG( "x_max= " << x_max <<" x_min = "<< x_min );
    PRINT_DEBUG( "y_max= " << y_max <<" y_min = "<< y_min );

    // n is the number of grid points.
    //if n is not given in the constructor then this number
    //is set to be the number of vertices in the arrangement.
    int n; 
    if (number_of_landmarks > 0)
      n = number_of_landmarks;
    else
      n= arr->number_of_vertices();

    //calculate the step size
    sqrt_n =  static_cast<int> (::sqrt(n) + 0.99999);
    FT delta_x = right.x() - left.x();
    FT delta_y = top.y() - bottom.y();
    step_x = delta_x / (sqrt_n-1);
    step_y = delta_y / (sqrt_n-1);

    PRINT_DEBUG( "n= " << n <<" sqrt_n = "<< sqrt_n );
    PRINT_DEBUG( "step_x= " << step_x <<" step_y = "<< step_y );
    PRINT_DEBUG( "left= " << left <<" right = "<< right );
    PRINT_DEBUG( "top= " << top <<" bottom = "<< bottom );

    int i, j; //i : x-indedx, j : y-index

    for (i=0; i< sqrt_n; i++)
    {
      for (j=0; j< sqrt_n; j++)
      {
        Point_2 p(left.x() + i*step_x, bottom.y() + j*step_y);

        //put in a list 
        points.push_back(p); 

        PRINT_DEBUG("grid point ("<<i<<','<<j<<") is= " << p);
      }
    }

    //double px,py;
    //int count; count = 0;

    // create the grid points.
    //for (px = x_min; px < x_max; px += step_x) 
    //{
    //  for (py = y_min; py < y_max; py += step_y) 
    //  {
    //    Point_2 p(px, py);

    //    //put in a list 
    //    points.push_back(p); 

    //    PRINT_DEBUG("grid point "<< count++ <<" is= " << p);        
    //  }
    //}

    PRINT_DEBUG("end create_grid_points_list");
  }


public:
  	  //Observer functions that are relevant to overload
	  //-------------------------------------------------

    /*! 
    * Notification before the arrangement is assigned with another
    * arrangement.
    * \param arr The arrangement to be copied.
    */
    virtual void before_assign (const Arrangement_2& arr)
    {
      clear_landmarks_set();
      traits = static_cast<const Traits_wrapper_2*> (arr.get_traits());
		  ignore_notifications = true;   
    }
    /*!
    * Notification after the arrangement has been assigned with another
    * arrangement.
    * \param u A handle to the unbounded face.
    */
    virtual void after_assign ()
    { 
      build_landmarks_set();
		  ignore_notifications = false;
    }

    /*! 
    * Notification before the observer is attached to an arrangement.
    * \param arr The arrangement we are about to attach the observer to.
    */
    virtual void before_attach (const Arrangement_2& arr)
    {
		  clear_landmarks_set();
		  traits = static_cast<const Traits_wrapper_2*> (arr.get_traits());
		  ignore_notifications = true;
    }

    /*!
    * Notification after the observer has been attached to an arrangement.
    */
    virtual void after_attach ()
    {
		  build_landmarks_set();
		  ignore_notifications = false;
    }

    /*! 
    * Notification before the observer is detached from the arrangement.
    */
    virtual void before_detach ()
    {
		  clear_landmarks_set();
    }

    /*!
    * Notification after the arrangement is cleared.
    * \param u A handle to the unbounded face.
    */
    virtual void after_clear (Face_handle /* u */)
    { 
      clear_landmarks_set();
      build_landmarks_set();
    }

    /*! Notification before a global operation modifies the arrangement. */
    virtual void before_global_change ()
    { 
      clear_landmarks_set();
      ignore_notifications = true;
    }

    /*! Notification after a global operation is completed. */
    virtual void after_global_change ()
    {
      build_landmarks_set();
      ignore_notifications = false;
    }

    /*!
    * Notification after the creation of a new vertex.
    * \param v A handle to the created vertex.
    */
    virtual void after_create_vertex (Vertex_handle /* v */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after the creation of a new edge.
    * \param e A handle to one of the twin halfedges that were created.
    */
    virtual void after_create_edge (Halfedge_handle /* e */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after an edge was split.
    * \param e1 A handle to one of the twin halfedges forming the first edge.
    * \param e2 A handle to one of the twin halfedges forming the second edge.
    */
    virtual void after_split_edge (Halfedge_handle /* e1 */,
      Halfedge_handle /* e2 */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after a face was split.
    * \param f A handle to the face we have just split.
    * \param new_f A handle to the new face that has been created.
    * \param is_hole Whether the new face forms a hole inside f.
    */
    virtual void after_split_face (Face_handle /* f */,
      Face_handle /* new_f */,
      bool /* is_hole */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after a hole was created inside a face.
    * \param h A circulator representing the boundary of the new hole.
    */
    virtual void after_add_hole (Ccb_halfedge_circulator /* h */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after an edge was merged.
    * \param e A handle to one of the twin halfedges forming the merged edge.
    */
    virtual void after_merge_edge (Halfedge_handle /* e */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after a face was merged.
    * \param f A handle to the merged face.
    */
    virtual void after_merge_face (Face_handle /* f */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification after a hole is moved from one face to another.
    * \param h A circulator representing the boundary of the hole.
    */
    virtual void after_move_hole (Ccb_halfedge_circulator /* h */)
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notificaion before the removal of a vertex.
    * \param v A handle to the vertex to be deleted.
    */
    virtual void after_remove_vertex ()
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification before the removal of an edge.
    * \param e A handle to one of the twin halfedges to be deleted.
    */
    virtual void after_remove_edge ()
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

    /*!
    * Notification before the removal of a hole.
    * \param h A circulator representing the boundary of the hole.
    */
    virtual void after_remove_hole ()
    {
      if (! ignore_notifications)
      {
        clear_landmarks_set();
        build_landmarks_set();
      }
    }

};

CGAL_END_NAMESPACE


#endif
