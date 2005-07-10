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
#ifndef CGAL_ARR_LANDMARKS_GENERATOR_H
#define CGAL_ARR_LANDMARKS_GENERATOR_H

/*! \file
* Definition of the Arr_landmarks_generator<Arrangement> template.
*/

#include <list>
#include <algorithm>   // for random_shuffle
#include <vector>   
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2/Arr_traits_wrapper_2.h>
#include <CGAL/Arr_point_location/Arr_lm_nearest_neighbor.h>
#include <CGAL/Arr_batched_point_location.h>


//#define CGAL_LM_DEBUG

//#ifdef CGAL_LM_DEBUG
//	#define PRINT_DEBUG(expr)   std::cout << expr << std::endl
//	#define LM_DEBUG(cmd)   cmd
//#else
//	#define PRINT_DEBUG(expr)
//	#define LM_DEBUG(cmd) 
//#endif


CGAL_BEGIN_NAMESPACE

/*! \class
* This class is related to the Landmarks point location, and given as 
* a parameter (or template parameter) to it. 
* It inherites from Arr_observer and updates this list whenever the 
* arrangement changes.
* This is a pure virtual class that handles the changes in a general manner 
* i.e. it builds the search structure for each small change in the arrangement. 
* The class tghat inherites from it must implement at least one virtual 
* function called "void _create_point_list(Point_list &)" which 
* actually creates the list of landmarks.
* I uses the batched point location to locate these landmarks inside 
*/
template <class Arrangement_, 
		  class Nearest_neighbor_ 
			  = Arr_landmarks_nearest_neighbor <typename Arrangement_::Traits_2> >
class Arr_landmarks_generator 
	: public Arr_observer <Arrangement_>
{
public:

	typedef Arrangement_							Arrangement_2;
	typedef typename Arrangement_2::Traits_2		Traits_2;

	typedef Arr_landmarks_generator<Arrangement_2, Nearest_neighbor_> 	Self;

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

	typedef typename Traits_2::Point_2					      Point_2;
	typedef typename Traits_2::X_monotone_curve_2		  X_monotone_curve_2;	

	typedef Nearest_neighbor_							            Nearest_neighbor;
	typedef typename Nearest_neighbor_::NN_Point_2		NN_Point_2;
	typedef std::list<NN_Point_2>						          NN_Points_set;

	typedef std::vector<Point_2>						          Points_set;
	typedef std::pair<Point_2,CGAL::Object>				    PL_pair;
	typedef std::vector<PL_pair>						          Pairs_set;
	typedef typename std::vector<PL_pair>::iterator		Pairs_iterator;

protected:

	typedef Arr_traits_basic_wrapper_2<Traits_2>  Traits_wrapper_2;

	// Data members:
	const Arrangement_2     *p_arr;     // The associated arrangement.
	const Traits_wrapper_2  *traits;    // Its associated traits object.
	Nearest_neighbor		nn;			// The associated nearest neighbor obj
	bool  ignore_notifications;	
	bool  updated;
	int	  num_small_not_updated_changes;
	int	  num_landmarks;

private:

  /*! Copy constructor - not supported. */
  Arr_landmarks_generator (const Self& );

  /*! Assignment operator - not supported. */
  Self& operator= (const Self& );

	
public: 
	  /*! Constructor. */
	  Arr_landmarks_generator (const Arrangement_2& arr) : 
	      Arr_observer<Arrangement_2> (const_cast<Arrangement_2 &>(arr)), 
		  p_arr(&arr),
		  ignore_notifications (false), 
		  updated (false), 
		  num_small_not_updated_changes(0), 
		  num_landmarks(0)
	  {
		  PRINT_DEBUG("Arr_landmarks_generator constructor"); 
		  traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
		  // need to call this in the inherited class constructor
		  // build_landmarks_set();
	  }

	  /*!
	  * Creates the landmarks set (choosing the landmarks) , 
	  * and saving them in the nearest neighbor search structure.
	  * This is a pure virtual function (must be implemented in 
	  * the class that derives from this one)
	  */
	  virtual void build_landmarks_set ()
	  {
		 PRINT_DEBUG("build_landmarks_set."); 
		 NN_Points_set     nn_points; 

		 //Go over planar map, and insert all vertices as landmarks
		 _create_nn_points_set(nn_points);

		 //the search structure is now updated
		 PRINT_DEBUG("call to initialize the nearest neighbor search."); 
		 nn.clean();
		 nn.init(nn_points.begin(), nn_points.end());

		 // num_landmarks = ?
		 num_small_not_updated_changes = 0;
		 updated = true;
	  }

	  /*!
	  * clear the tree
	  */
	  virtual void clear_landmarks_set ()
	  {
		  PRINT_DEBUG("clear_landmarks_set.");

		  nn.clean();

		  num_landmarks = 0;
		  num_small_not_updated_changes = 0;
		  updated = false;		  
	  }

	  /*!
	  * get the nearest neighbor (landmark) to the given point
	  */
	  virtual Point_2 & get_closest_landmark (Point_2 p, Object &obj)
	  {
		  CGAL_assertion(updated);
		  return nn.find_nearest_neighbor(p, obj);
	  }

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
      p_arr = &arr;
      traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
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
		  p_arr = &arr; 
		  traits = static_cast<const Traits_wrapper_2*> (p_arr->get_traits());
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

protected:
  /*!
  * This function creates the list of landmarks with their location.
  * This is a pure virtual function, and the class that inherites from 
  * this generator must implement it.
  */
  virtual void _create_points_set (Points_set &) = 0;

  virtual void _create_nn_points_set (NN_Points_set &nn_points) 
  {
    Points_set		points;
    Pairs_set			pairs;

    //call the function that creates the landmarks 
    _create_points_set(points);

		PRINT_DEBUG("before batched point location."); 

    //locate the landmarks in the arrangement using batched point location
    // global function.
    locate(*p_arr,points.begin(),points.end(),std::back_inserter(pairs));

    //random shuffle of the points since the batched p.l. sorts them
    std::random_shuffle ( pairs.begin (), pairs.end ());

		PRINT_DEBUG("after batched point location + shuffle."); 

    //create the nn set 
    Pairs_iterator itr;
    for(itr = pairs.begin(); itr != pairs.end(); ++itr)
    {
      NN_Point_2 np(itr->first, itr->second); 
      nn_points.push_back(np);
    }
  }


};

CGAL_END_NAMESPACE


#endif
